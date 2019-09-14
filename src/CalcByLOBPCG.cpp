/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/**@file
@brief Functions to perform calculations with the
localy optimal block (preconditioned) conjugate gradient method.
*/
#include "CalcByLOBPCG.hpp"
#include "xsetmem.hpp"
#include "mltply.hpp"
#include "FileIO.hpp"
#include "wrapperMPI.hpp"
#include "expec_cisajs.hpp"
#include "expec_cisajscktaltdc.hpp"
#include "expec_totalspin.hpp"
#include "expec_energy_flct.hpp"
#include "phys.hpp"
#include "mltplyCommon.hpp"
#include "./common/setmemory.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include "dSFMT.hpp"
#include "log.hpp"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

int     initial_mode;/**< mode to get initial state (0: use same random generator for MPI, 1: use each random generator for MPI)*/

extern "C" {
  extern void zheevd_(char *jobz, char *uplo, int *n, std::complex<double> *a, int *lda, double *w, std::complex<double> *work, int *lwork, double *rwork, int * lrwork, int *iwork, int *liwork, int *info);
  extern void zgemm_(char *transa, char *transb, int *m, int *n, int *k, std::complex<double> *alpha, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb, std::complex<double> *beta, std::complex<double> *c, int *ldc);
}
/**
@brief Solve the generalized eigenvalue problem
@f[
{\hat H} |\phi\rangle = \varepsilon {\hat O} |\phi\rangle
@f]
with the Lowdin's orthogonalization
@return the truncated dimension, nsub2
*/
int CalcByLOBPCG::diag_ovrp(
  int nsub,//!<[in] Original dimension of subspace
  std::complex<double> *hsub,//!<[inout] (nsub*nsub) subspace hamiltonian -> eigenvector
  std::complex<double> *ovlp,//!<[inout] (nsub*nsub) Overrap matrix -> @f${\hat O}^{1/2}@f$
  double *eig//!<[out] (nsub) Eigenvalue
)
{
  int *iwork, info, isub, jsub, nsub2;
  char jobz = 'V', uplo = 'U', transa = 'N', transb = 'N';
  double *rwork;
  std::complex<double> *work, *mat;
  int liwork, lwork, lrwork;
  std::complex<double> one = 1.0, zero = 0.0;

  liwork = 5 * nsub + 3;
  lwork = nsub*nsub + 2 * nsub;
  lrwork = 3 * nsub*nsub + (4 + (int)log2(nsub) + 1) * nsub + 1;

  iwork = (int*)malloc(liwork * sizeof(int));
  rwork = (double*)malloc(lrwork * sizeof(double));
  work = (std::complex<double>*)malloc(lwork * sizeof(std::complex<double>));
  mat = (std::complex<double>*)malloc(nsub*nsub * sizeof(std::complex<double>));
  for (isub = 0; isub < nsub*nsub; isub++)mat[isub] = 0.0;
  /**@brief
  (1) Compute @f${\hat O}^{-1/2}@f$ with diagonalizing overrap matrix
  */
  zheevd_(&jobz, &uplo, &nsub, ovlp, &nsub, eig, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  /**@brief
  @f[
   {\hat O}^{-1/2} = \left(\frac{|O_1\rangle}{\sqrt{o_1}}, \frac{|O_2\rangle}{\sqrt{o_2}},
   ...\right)
   \\
   {\hat O} |O_i\rangle = o_i |O_i\rangle
  @f]
  if @f$o_i@f$ is very small that dimension is ignored. Therefore @f${\hat O}^{-1/2}@f$
  is nsub*nsub2 matrix.
  */
  nsub2 = 0;
  for (isub = 0; isub < nsub; isub++) {
    if (eig[isub] > 1.0e-14) {
      for (jsub = 0; jsub < nsub; jsub++)
        ovlp[jsub + nsub*nsub2] = ovlp[jsub + nsub*isub] / sqrt(eig[isub]);
      nsub2 += 1;
    }
  }
  for (isub = nsub2; isub < nsub; isub++) 
    for (jsub = 0; jsub < nsub; jsub++)
      ovlp[jsub + nsub*isub] = 0.0;
  /**
  (2) Transform @f${\hat H}'\equiv {\hat O}^{-1/2 \dagger}{\hat H}{\hat O}^{-1/2}@f$.
  @f${\hat H}'@f$ is nsub2*nsub2 matrix.
  */
  transa = 'N';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, hsub, &nsub, ovlp, &nsub, &zero, mat, &nsub);
  transa = 'C';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, ovlp, &nsub, mat, &nsub, &zero, hsub, &nsub);
  /**
  (3) Diagonalize @f${\hat H}'@f$. It is the standard eigenvalue problem.
  @f[
  {\hat H}' |\phi'_i\rangle = \varepsilon_i |\phi'_i\rangle
  @f]
  */
  zheevd_(&jobz, &uplo, &nsub2, hsub, &nsub, eig, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  /**
  (4) Transform eigenvector into the original nsub space as
  @f[
  |\phi_i\rangle = {\hat O}^{-1/2} |\phi'_i\rangle
  @f]
  */
  transa = 'N';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, ovlp, &nsub, hsub, &nsub, &zero, mat, &nsub);
 // printf("%d %d %15.5f %15.5f %15.5f\n", info, nsub2, eig[0], eig[1], eig[2]);
  for (isub = 0; isub < nsub*nsub; isub++)hsub[isub] = mat[isub];

  free(mat);
  free(work);
  free(rwork);
  free(iwork);

  return(nsub2);
}/*void diag_ovrp*/
/**@brief
Compute adaptively shifted preconditionar written in
S. Yamada, et al., Transactions of JSCES, Paper No. 20060027 (2006).
@return adaptive shift for preconditioning
*/
static double calc_preshift(
  double eig,//!<[in] Eigenvalue in this step
  double res,//!<[in] Residual 2-norm in this step
  double eps_LOBPCG//!<[in] Convergence threshold
)
{
  double k, i;
  double preshift;

  if (fabs(eig) > 10.0) k = trunc(log10(fabs(eig)));
  else k = 1.0;

  if (res < 1.0) {
    if (eps_LOBPCG > res) i = ceil(log10(eps_LOBPCG));
    else i = ceil(log10(res));

    preshift = trunc(eig / pow(10.0, k + i))*pow(10.0, k + i);
  }
  else preshift = 0.0;

  return(preshift);
}/*void calc_preshift*/
/*
@brief Compute initial guess for LOBPCG.
If this is resuterting run, read from files.
*/
static void Initialize_wave(
  std::complex<double> **wave//!<[out] [CheckList::idim_max][exct] initial eigenvector
) 
{
  FILE *fp;
  char sdt[D_FileNameMax];
  size_t byte_size;
  std::complex<double> *vin;
  int ie;
  int iproc, ierr;
  long int iv;
  long int i_max_tmp, sum_i_max, idim, i_max;
  int mythread;
  double *dnorm;
  /*
  For DSFMT
  */
  long int u_long_i;
  dsfmt_t dsfmt;
  /**@brief
  (A) For restart: Read saved eigenvector files (as binary files) from each processor
  */
  if (Def::iReStart == DC::RESTART_INOUT || Def::iReStart == DC::RESTART_IN) {
    //StartTimer(3600);
    //TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector starts: %s\n", "a", rand_i, step_i);
    fprintf(MP::STDOUT, "%s", "  Start:  Input vector.\n");

    ierr = 0;
    vin = cd_1d_allocate(Check::idim_max);
    for (ie = 0; ie < Def::k_exct; ie++) {

      sprintf(sdt, "tmpvec_set%d_rank_%d.dat", ie, MP::myrank);
      childfopenALL(sdt, "rb", &fp);
      if (fp == NULL) {
        fprintf(stdout, "Restart file is not found.\n");
        fprintf(stdout, "Start from scratch.\n");
        ierr = 1;
        break;
      }
      else {
        byte_size = fread(&iproc, sizeof(int), 1, fp);
        byte_size = fread(&i_max, sizeof(long int), 1, fp);
        //fprintf(MP::STDOUT, "Debug: i_max=%ld, step_i=%d\n", i_max, step_i);
        if (i_max != Check::idim_max) {
          fprintf(stderr, "Error: Invalid restart file.\n");
          wrapperMPI::Exit(-1);
        }
        byte_size = fread(vin, sizeof(std::complex<double>), Check::idim_max, fp);
        for (idim = 0; idim < i_max; idim++) wave[idim][ie] = vin[idim];
        fclose(fp);
      }
    }/*for (ie = 0; ie < Def::k_exct; ie++)*/
    free_cd_1d_allocate(vin);

    if (ierr == 0) {
      //TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector finishes: %s\n", "a", rand_i, step_i);
      fprintf(MP::STDOUT, "%s", "  End  :  Input vector.\n");
      //StopTimer(3600);
      if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
      return;
    }/*if (ierr == 0)*/

  }/*Def::iReStart == DC::RESTART_INOUT || Def::iReStart == DC::RESTART_IN*/

  /**@brief
  (B) For scratch (including the case that restart files are not found):
  initialize eigenvectors in the same way as TPQ and Lanczos.
  */
  i_max = Check::idim_max;

  if (initial_mode == 0) {

    for (ie = 0; ie < Def::k_exct; ie++) {

      sum_i_max = wrapperMPI::Sum_li(Check::idim_max);
      Large::iv = (sum_i_max / 2 + Def::initial_iv + ie) % sum_i_max + 1;
      iv = Large::iv;
      fprintf(MP::STDOUT, "  initial_mode=%d normal: iv = %ld i_max=%ld k_exct =%d\n\n", 
        initial_mode, iv, i_max, Def::k_exct);
#pragma omp parallel for default(none) private(idim) shared(wave,i_max,ie)
      for (idim = 0; idim < i_max; idim++) wave[idim][ie] = 0.0;
      
      sum_i_max = 0;
      for (iproc = 0; iproc < MP::nproc; iproc++) {

        i_max_tmp = wrapperMPI::Bcast_li(iproc, i_max);
        if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp) {

          if (MP::myrank == iproc) {
            wave[iv - sum_i_max][ie] = 1.0;
            if (Def::iInitialVecType == 0) {
              wave[iv - sum_i_max][ie] += std::complex < double>(0.0, 1.0);
              wave[iv - sum_i_max][ie] /= sqrt(2.0);
            }
          }/*if (MP::myrank == iproc)*/
        }/*if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp)*/

        sum_i_max += i_max_tmp;

      }/*for (iproc = 0; iproc < MP::nproc; iproc++)*/
    }/*for (ie = 0; ie < Def::k_exct; ie++)*/
  }/*if(initial_mode == 0)*/
  else if (initial_mode == 1) {
    iv = Def::initial_iv;
    fprintf(MP::STDOUT, "  initial_mode=%d (random): iv = %ld i_max=%ld k_exct =%d\n\n",
      initial_mode, iv, i_max, Def::k_exct);
#pragma omp parallel default(none) private(idim, u_long_i, mythread, dsfmt, ie) \
shared(wave, iv, MP::nthreads, MP::myrank, i_max,Def::k_exct,Def::iInitialVecType)
    {
      /*
       Initialise MT
       */
#ifdef _OPENMP
      mythread = omp_get_thread_num();
#else
      mythread = 0;
#endif
      u_long_i = 123432 + labs(iv) + mythread + MP::nthreads * MP::myrank;
      dsfmt_init_gen_rand(&dsfmt, u_long_i);

      for (ie = 0; ie < Def::k_exct; ie++) {
        if (Def::iInitialVecType == 0) {
#pragma omp for
          for (idim = 0; idim < i_max; idim++)
            wave[idim][ie] = std::complex < double>(
              2.0 * (dsfmt_genrand_close_open(&dsfmt) - 0.5), 
              2.0 * (dsfmt_genrand_close_open(&dsfmt) - 0.5));
        }
        else {
#pragma omp for
          for (idim = 0; idim < i_max; idim++)
            wave[idim][ie] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
        }
      }/*for (ie = 0; ie < Def::k_exct; ie++)*/

    }/*#pragma omp parallel*/

    dnorm = d_1d_allocate(Def::k_exct);
    wrapperMPI::Norm_dv(i_max, Def::k_exct, wave, dnorm);
#pragma omp parallel for default(none) private(idim,ie) \
shared(i_max,wave,dnorm,Def::k_exct)
    for (idim = 0; idim < i_max; idim++) 
      for (ie = 0; ie < Def::k_exct; ie++) wave[idim][ie] /= dnorm[ie];
    free_d_1d_allocate(dnorm);
  }/*else if(initial_mode==1)*/
}/*static void Initialize_wave*/
/**
@brief Output eigenvectors for restart LOBPCG method
*/
void CalcByLOBPCG::Output_restart(
  std::complex<double> **wave//!<[in] [exct][CheckList::idim_max] initial eigenvector
)
{
  FILE *fp;
  size_t byte_size;
  char sdt[D_FileNameMax];
  int ie;
  long int idim;
  std::complex<double> *vout;

  //TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector starts: %s\n", "a", rand_i, step_i);
  fprintf(MP::STDOUT, "%s", "  Start:  Output vector.\n");
  
  vout = cd_1d_allocate(Check::idim_max);
  for (ie = 0; ie < Def::k_exct; ie++) {
    sprintf(sdt, "tmpvec_set%d_rank_%d.dat", ie, MP::myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) wrapperMPI::Exit(-1);
    byte_size = fwrite(&Large::itr, sizeof(Large::itr), 1, fp);
    byte_size = fwrite(&Check::idim_max, sizeof(Check::idim_max), 1, fp);
    for (idim = 0; idim < Check::idim_max; idim++) vout[idim] = wave[idim][ie];
    byte_size = fwrite(vout, sizeof(std::complex<double>), Check::idim_max, fp);
    fclose(fp);
  }/*for (ie = 0; ie < Def::k_exct; ie++)*/
  free_cd_1d_allocate(vout);

  //TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector finishes: %s\n", "a", rand_i, step_i);
  fprintf(MP::STDOUT, "%s", "  End  :  Output vector.\n");
  if(byte_size == 0) printf("byte_size : %d\n", (int)byte_size);
}/*static void Output_restart*/
/**@brief
Core routine for the LOBPCG method
This method is introduced in 
-# S. Yamada, et al., Transactions of JSCES, Paper No. 20060027 (2006).
   https://www.jstage.jst.go.jp/article/jsces/2006/0/2006_0_20060027/_pdf
-# A. V. Knyazev, SIAM J. Sci.  Compute. 23, 517 (2001).
   http://epubs.siam.org/doi/pdf/10.1137/S1064827500366124
*/
int CalcByLOBPCG::LOBPCG_Main()
{
  char sdt[D_FileNameMax], sdt_2[D_FileNameMax];
  FILE *fp;
  int iconv = -1, i4_max;
  long int idim, i_max;
  int ie, stp;
  int ii, jj, nsub, nsub_cut, nstate;
  std::complex<double> ***wxp/*[0] w, [1] x, [2] p of Ref.1*/,
    ***hwxp/*[0] h*w, [1] h*x, [2] h*p of Ref.1*/,
    ****hsub, ****ovlp; /*Subspace Hamiltonian and Overlap*/
  double *eig, *dnorm, eps_LOBPCG, eigabs_max, preshift, precon, dnormmax, *eigsub;
  int do_precon = 0;//If = 1, use preconditioning (experimental)
  char tN = 'N', tC = 'C';
  std::complex<double> one = 1.0, zero = 0.0;

  nsub = 3 * Def::k_exct;
  nstate = Def::k_exct;

  eig = d_1d_allocate(Def::k_exct);
  dnorm = d_1d_allocate(Def::k_exct);
  eigsub = d_1d_allocate(nsub);
  hsub = cd_4d_allocate(3, Def::k_exct, 3, Def::k_exct);
  ovlp = cd_4d_allocate(3, Def::k_exct, 3, Def::k_exct);

  i_max = Check::idim_max;
  i4_max = (int)i_max;

  free_cd_2d_allocate(Wave::v0);
  free_cd_2d_allocate(Wave::v1);
  wxp = cd_3d_allocate(3, Check::idim_max, Def::k_exct);
  hwxp = cd_3d_allocate(3, Check::idim_max, Def::k_exct);
  /**@brief
  <ul>
  <li>Set initial guess of wavefunction: 
  @f${\bf x}=@f$initial guess</li>
  */
  Initialize_wave(wxp[1]);

  TimeKeeper("%s_TimeKeeper.dat", "Lanczos Eigen Value start:    %s", "a");

  zclear(i_max*Def::k_exct, &hwxp[1][0][0]);
  mltply::main (Def::k_exct, hwxp[1], wxp[1]);
  stp = 1;
  TimeKeeperWithStep("%s_TimeKeeper.dat", "%3d th Lanczos step: %s", "a", 0);

  zclear(i_max*Def::k_exct, &wxp[2][0][0]);
  zclear(i_max*Def::k_exct, &hwxp[2][0][0]);
  for (ie = 0; ie < Def::k_exct; ie++) eig[ie] = 0.0;
  for (idim = 0; idim < i_max; idim++) {
    for (ie = 0; ie < Def::k_exct; ie++) {
      wxp[2][idim][ie] = 0.0;
      hwxp[2][idim][ie] = 0.0;
      eig[ie] += real(conj(wxp[1][idim][ie]) * hwxp[1][idim][ie]);
    }
  }
  wrapperMPI::Sum_dv(Def::k_exct, eig);

  sprintf(sdt_2, "%s_Lanczos_Step.dat", Def::CDataFileHead);
  childfopenMPI(sdt_2, "w", &fp);
  fprintf(MP::STDOUT, "    Step   Residual-2-norm     Threshold      Energy\n");
  fprintf(fp, "    Step   Residual-2-norm     Threshold      Energy\n");
  fclose(fp);

  nsub_cut = nsub;
  /**@brief
  <li>@b DO LOBPCG loop</li>
  <ul>
  */
  for (stp = 1; stp <= Def::Lanczos_max; stp++) {
    /**@brief
    <li>Scale convergence threshold with the absolute value of eigenvalue
    for numerical stability</li>
    */
    eigabs_max = 0.0;
    for (ie = 0; ie < Def::k_exct; ie++)
      if (fabs(eig[ie]) > eigabs_max) eigabs_max = fabs(eig[ie]);
    eps_LOBPCG = pow(10, -0.5 *Def::LanczosEps);
    if (eigabs_max > 1.0) eps_LOBPCG *= eigabs_max;
    /**@brief
    <li>@b DO each eigenvector</li>
    <ul>
    */
    /**@brief
     <li>Compute residual vectors: @f${\bf w}={\bf X}-\mu {\bf x}@f$</li>
    */
#pragma omp parallel for default(none) private(idim,ie) \
shared(i_max,wxp,hwxp,eig,Def::k_exct)
    for (idim = 0; idim < i_max; idim++) {
      for (ie = 0; ie < Def::k_exct; ie++) {
        wxp[0][idim][ie] = hwxp[1][idim][ie] - eig[ie] * wxp[1][idim][ie];
      }
    }        
    wrapperMPI::Norm_dv(i_max, Def::k_exct, wxp[0], dnorm);

    dnormmax = 0.0;
    for (ie = 0; ie < Def::k_exct; ie++) 
      if (dnorm[ie] > dnormmax) dnormmax = dnorm[ie];
    /**@brief
    <li>Preconditioning (Point Jacobi): @f${\bf w}={\hat T}^{-1} {\bf w}@f$</li>
    */
    if (stp /= 1) {
      if (do_precon == 1) {
        for (ie = 0; ie < Def::k_exct; ie++) 
          preshift = calc_preshift(eig[ie], dnorm[ie], eps_LOBPCG);
#pragma omp parallel for default(none) private(idim,precon,ie) \
shared(wxp,List::Diagonal,preshift,i_max,eps_LOBPCG,Def::k_exct)
        for (idim = 0; idim < i_max; idim++) {
          for (ie = 0; ie < Def::k_exct; ie++){
            precon = List::Diagonal[idim] - preshift;
            if (fabs(precon) > eps_LOBPCG) wxp[0][idim][ie] /= precon;
          }
        }
      }/*if(do_precon == 1)*/
      /**@brief
        <li>Normalize residual vector: @f${\bf w}={\bf w}/|w|@f$
      */
      wrapperMPI::Norm_dv(i_max, Def::k_exct, wxp[0], dnorm);
#pragma omp parallel for default(none) private(idim,ie) \
shared(i_max,wxp,dnorm,Def::k_exct)
      for (idim = 0; idim < i_max; idim++)
        for (ie = 0; ie < Def::k_exct; ie++)
          wxp[0][idim][ie] /= dnorm[ie];
    }/*if (stp /= 1)*/
    /**@brief
    </ul>
    <li>@b END @b DO each eigenvector</li>
    <li>Convergence check</li>
    */
    childfopenMPI(sdt_2, "a", &fp);
    fprintf(MP::STDOUT, "%9d %15.5e %15.5e      ", stp, dnormmax, eps_LOBPCG);
    fprintf(fp, "%9d %15.5e %15.5e      ", stp, dnormmax, eps_LOBPCG);
    for (ie = 0; ie < Def::k_exct; ie++) {
      fprintf(MP::STDOUT, " %15.5e", eig[ie]);
      fprintf(fp, " %15.5e", eig[ie]);
    }
    if(nsub_cut == 0) printf("nsub_cut : %d", nsub_cut);
    fprintf(MP::STDOUT, "\n");
    fprintf(fp, "\n");
    fclose(fp);

    if (dnormmax < eps_LOBPCG) {
      iconv = 0;
      break;
    }
    /**@brief
    <li>@f${\bf W}={\hat H}{\bf w}@f$</li>
    */
    zclear(i_max*Def::k_exct, &hwxp[0][0][0]);
    mltply::main(Def::k_exct, hwxp[0], wxp[0]);

    TimeKeeperWithStep("%s_TimeKeeper.dat", "%3d th Lanczos step: %s", "a", stp);
    /**@brief
    <li>Compute subspace Hamiltonian and overrap matrix:
    @f${\hat H}_{\rm sub}=\{{\bf w},{\bf x},{\bf p}\}^\dagger \{{\bf W},{\bf X},{\bf P}\}@f$, 
    @f${\hat O}=\{{\bf w},{\bf x},{\bf p}\}^\dagger \{{\bf w},{\bf x},{\bf p}\}@f$,
    </li>
    */
    for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {
        zgemm_(&tN, &tC, &nstate, &nstate, &i4_max, &one,
          &wxp[ii][0][0], &nstate, &wxp[jj][0][0], &nstate, &zero, &ovlp[jj][0][ii][0], &nsub);
        zgemm_(&tN, &tC, &nstate, &nstate, &i4_max, &one,
          &wxp[ii][0][0], &nstate, &hwxp[jj][0][0], &nstate, &zero, &hsub[jj][0][ii][0], &nsub);
      }
    }
    wrapperMPI::Sum_cv(nsub*nsub, &ovlp[0][0][0][0]);
    wrapperMPI::Sum_cv(nsub*nsub, &hsub[0][0][0][0]);

    for (ie = 0; ie < Def::k_exct; ie++)
      eig[ie] = real(hsub[1][ie][1][ie]);
    /**@brief
    <li>Subspace diagonalization with the Lowdin's orthogonalization for
        generalized eigenvalue problem: @f${\hat H}_{\rm sub}{\bf v}={\hat O}\mu_{\rm sub}{\bf v}@f$,
        @f${\bf v}=(\alpha, \beta, \gamma)@f$</li>
    */
    nsub_cut = CalcByLOBPCG::diag_ovrp(nsub, &hsub[0][0][0][0], &ovlp[0][0][0][0], eigsub);
    /**@brief
    <li>Update @f$\mu=(\mu+\mu_{\rm sub})/2@f$</li>
    */
    for (ie = 0; ie < Def::k_exct; ie++)
      eig[ie] = 0.5 * (eig[ie] + eigsub[ie]);
    /**@brief
      <li>@f${\bf x}=\alpha {\bf w}+\beta {\bf x}+\gamma {\bf p}@f$,
      Normalize @f${\bf x}@f$</li>
    */
    zclear(i_max*Def::k_exct, &Wave::v1buf[0][0]);
    for (ii = 0; ii < 3; ii++) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &wxp[ii][0][0], &nstate, &one, &Wave::v1buf[0][0], &nstate);
    }
    for (idim = 0; idim < i_max; idim++) for (ie = 0; ie < Def::k_exct; ie++)
      wxp[1][idim][ie] = Wave::v1buf[idim][ie];
    /**@brief
    <li>@f${\bf X}=\alpha {\bf W}+\beta {\bf X}+\gamma {\bf P}@f$,
    Normalize @f${\bf X}@f$</li>
    */
    zclear(i_max*Def::k_exct, &Wave::v1buf[0][0]);
    for (ii = 0; ii < 3; ii++) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &hwxp[ii][0][0], &nstate, &one, &Wave::v1buf[0][0], &nstate);
    }
    for (idim = 0; idim < i_max; idim++) for (ie = 0; ie < Def::k_exct; ie++)
      hwxp[1][idim][ie] = Wave::v1buf[idim][ie];
    /**@brief
    <li>@f${\bf p}=\alpha {\bf w}+\gamma {\bf p}@f$,
    Normalize @f${\bf p}@f$</li>
    */
    zclear(i_max*Def::k_exct, &Wave::v1buf[0][0]);
    for (ii = 0; ii < 3; ii += 2) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &wxp[ii][0][0], &nstate, &one, &Wave::v1buf[0][0], &nstate);
    }
    for (idim = 0; idim < i_max; idim++) for (ie = 0; ie < Def::k_exct; ie++)
      wxp[2][idim][ie] = Wave::v1buf[idim][ie];
    /**@brief
    <li>@f${\bf P}=\alpha {\bf W}+\gamma {\bf P}@f$,
    Normalize @f${\bf P}@f$</li>
    */
    zclear(i_max*Def::k_exct, &Wave::v1buf[0][0]);
    for (ii = 0; ii < 3; ii += 2) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &hwxp[ii][0][0], &nstate, &one, &Wave::v1buf[0][0], &nstate);
    }
    for (idim = 0; idim < i_max; idim++) for (ie = 0; ie < Def::k_exct; ie++)
      hwxp[2][idim][ie] = Wave::v1buf[idim][ie];
    /**@brief
    <li>Normalize @f${\bf w}@f$ and @f${\bf W}@f$</li>
    */
    for (ii = 1; ii < 3; ii++) {
      wrapperMPI::Norm_dv(i_max, Def::k_exct, wxp[ii], dnorm);
#pragma omp parallel for default(none) private(idim,ie) \
shared(i_max,wxp,hwxp,dnorm,ii, Def::k_exct)
      for (idim = 0; idim < i_max; idim++) {
        for (ie = 0; ie < Def::k_exct; ie++) {
          wxp[ii][idim][ie] /= dnorm[ie];
          hwxp[ii][idim][ie] /= dnorm[ie];
        }/* for (ie = 0; ie < Def::k_exct; ie++)*/
      }
    }/*for (ii = 1; ii < 3; ii++)*/

  }/*for (stp = 1; stp <= Def::Lanczos_max; stp++)*/
  /**@brief
  </ul>
  <li>@b END @b DO LOBPCG iteration
  */
  //fclose(fp);

  Large::itr = stp;
  sprintf(sdt, "%s_TimeKeeper.dat", Def::CDataFileHead);

  TimeKeeper("%s_TimeKeeper.dat", "Lanczos Eigenvalue finishes:  %s", "a");
  fprintf(MP::STDOUT, "%s", "\n######  End  : Calculate Lanczos EigenValue.  ######\n\n");

  free_d_1d_allocate(eig);
  free_d_1d_allocate(dnorm);
  free_d_1d_allocate(eigsub);
  free_cd_4d_allocate(hsub);
  free_cd_4d_allocate(ovlp);
  free_cd_3d_allocate(hwxp);
  /**@brief
  <li>Output resulting vectors for restart</li>
  */
  if (Def::iReStart == DC::RESTART_OUT || Def::iReStart == DC::RESTART_INOUT){
      Output_restart(wxp[1]);
      if(iconv != 0) {
          sprintf(sdt, "%s", "Lanczos Eigenvalue is not converged in this process.");
          return 1;
      }
  }
  /**@brief
  <li>Just Move wxp[1] into ::Wave::v1. The latter must be start from 0-index (the same as FullDiag)</li>
  </ul>
  */
  Wave::v0 = cd_2d_allocate(Check::idim_max, Def::k_exct);
#pragma omp parallel for default(none) shared(i_max,wxp,Wave::v0,Def::k_exct) private(idim,ie)
  for (idim = 0; idim < i_max; idim++)
    for (ie = 0; ie < Def::k_exct; ie++) 
      Wave::v0[idim][ie] = wxp[1][idim][ie];
  free_cd_3d_allocate(wxp);
  Wave::v1 = cd_2d_allocate(Check::idim_max, Def::k_exct);

  if (iconv != 0) {
    sprintf(sdt, "%s", "Lanczos Eigenvalue is not converged in this process.");
    return -1;
  }
  else {
    return 0;
  }
}/*int LOBPCG_Main*/
/**
@brief Driver routine for LOB(P)CG method.
*/
int CalcByLOBPCG::main(
)
{
  char sdt[D_FileNameMax];
  size_t byte_size;
  long int ie;
  long int i_max = 0;
  long int idim;
  FILE *fp;
  std::complex<double> *vin;

  fprintf(MP::STDOUT, "######  Eigenvalue with LOBPCG  #######\n\n");

  if (Def::iInputEigenVec == FALSE) {

    // this part will be modified
    switch (Def::iCalcModel) {
    case DC::HubbardGC:
    case DC::SpinGC:
    case DC::KondoGC:
    case DC::SpinlessFermionGC:
      initial_mode = 1; // 1 -> random initial vector
      break;
    case DC::Hubbard:
    case DC::Kondo:
    case DC::Spin:
    case DC::SpinlessFermion:

      if (Def::iFlgGeneralSpin == TRUE) {
        initial_mode = 1;
      }
      else {
        if (Def::initial_iv>0) {
          initial_mode = 0; // 0 -> only v[iv] = 1
        }
        else {
          initial_mode = 1; // 1 -> random initial vector
        }
      }
      break;
    default:
      //fclose(fp);
      wrapperMPI::Exit(-1);
    }

    int iret = LOBPCG_Main();
    if (iret != 0) {
      if(iret ==1) return (TRUE);
      else{
          fprintf(MP::STDOUT, "  LOBPCG is not converged in this process.\n");
          return(FALSE);
      }
    }
  }/*if (Def::iInputEigenVec == FALSE)*/
  else {// Def::iInputEigenVec=true :input Wave::v1:
    /**@brief
    If this run is for spectrum calculation, eigenvectors are not computed
    and read from files.
    */
    fprintf(MP::STDOUT, "An Eigenvector is inputted.\n");
    vin = cd_1d_allocate(Check::idim_max);
    for (ie = 0; ie < Def::k_exct; ie++) {
      TimeKeeper("%s_TimeKeeper.dat", "Read Eigenvector starts:          %s", "a");
      sprintf(sdt, "%s_eigenvec_%ld_rank_%d.dat", Def::CDataFileHead, ie, MP::myrank);
      childfopenALL(sdt, "rb", &fp);
      if (fp == NULL) {
        fprintf(stderr, "Error: Inputvector file is not found.\n");
        wrapperMPI::Exit(-1);
      }
      byte_size = fread(&Step::step_i, sizeof(int), 1, fp);
      byte_size = fread(&i_max, sizeof(long int), 1, fp);
      if (i_max != Check::idim_max) {
        fprintf(stderr, "Error: Invalid Inputvector file.\n");
        wrapperMPI::Exit(-1);
      }
      byte_size = fread(vin, sizeof(std::complex<double>), Check::idim_max, fp);
#pragma omp parallel for default(none) shared(Wave::v1,vin, i_max, ie), private(idim)
      for (idim = 0; idim < i_max; idim++) {
        Wave::v1[ie][idim] = vin[idim];
      }
    }/*for (ie = 0; ie < Def::k_exct; ie++)*/
    fclose(fp);
    free_cd_1d_allocate(vin);
    TimeKeeper("%s_TimeKeeper.dat", "Read Eigenvector finishes:        %s", "a");

    if(byte_size == 0) printf("byte_size : %d\n", (int)byte_size);
  }/*Def::iInputEigenVec == TRUE*/

  fprintf(MP::STDOUT, "%s", "\n######  End  : Calculate Lanczos EigenVec.  ######\n\n");
  /**@brief
    Compute & Output physical variables to a file
    the same function as FullDiag [phys()] is used.
  */
  phys(Def::k_exct);

  Def::St=1;
  if (Def::St == 0) {
    sprintf(sdt, "%s_energy.dat", Def::CDataFileHead);
  }
  else if (Def::St == 1) {
    sprintf(sdt, "%s_energy.dat", Def::CDataFileHead);
  }

  if (childfopenMPI(sdt, "w", &fp) != 0) {
    wrapperMPI::Exit(-1);
  }
  for (ie = 0; ie < Def::k_exct; ie++) {
    //phys(ie);
    fprintf(fp, "State %ld\n", ie);
    fprintf(fp, "  Energy  %.16lf \n", Phys::energy[ie]);
    fprintf(fp, "  Doublon  %.16lf \n", Phys::doublon[ie]);
    fprintf(fp, "  Sz  %.16lf \n", Phys::Sz[ie]);
    //fprintf(fp, "  S^2  %.16lf \n", Phys::s2[ie]);
    //fprintf(fp, "  N_up  %.16lf \n", Phys::num_up[ie]);
    //fprintf(fp, "  N_down  %.16lf \n", Phys::num_down[ie]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  /*
   Output Eigenvector to a file
  */
  if (Def::iOutputEigenVec == TRUE) {
    TimeKeeper("%s_TimeKeeper.dat", "Output Eigenvector starts:          %s", "a");

    vin = cd_1d_allocate(Check::idim_max);
    for (ie = 0; ie < Def::k_exct; ie++) {

#pragma omp parallel for default(none) shared(Wave::v1,ie,vin,Check::idim_max) private(idim)
      for (idim = 0; idim < Check::idim_max; idim++)
        vin[idim] = Wave::v1[idim][ie];
      
      sprintf(sdt, "%s_eigenvec_%ld_rank_%d.dat", Def::CDataFileHead, ie, MP::myrank);
      if (childfopenALL(sdt, "wb", &fp) != 0) wrapperMPI::Exit(-1);
      byte_size = fwrite(&Large::itr, sizeof(Large::itr), 1, fp);
      byte_size = fwrite(&Check::idim_max, sizeof(Check::idim_max), 1, fp);
      byte_size = fwrite(vin, sizeof(std::complex<double>), Check::idim_max, fp);
      fclose(fp);
    }/*for (ie = 0; ie < Def::k_exct; ie++)*/
    free_cd_1d_allocate(vin);

    TimeKeeper("%s_TimeKeeper.dat", "Output Eigenvector starts:          %s", "a");
  }/*if (Def::iOutputEigenVec == TRUE)*/

  return TRUE;

}/*int CalcByLOBPCG*/
