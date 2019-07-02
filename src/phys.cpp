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
#include "phys.hpp"
#include "expec_energy_flct.hpp"
#include "expec_totalspin.hpp"
#include "expec_cisajs.hpp"
#include "expec_cisajscktaltdc.hpp"
#include "wrapperMPI.hpp"
#ifdef _SCALAPACK
#include "matrixscalapack.hpp"
#endif

/**
 * @file
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for giving a parent function to calculate physical quantities  by full diagonalization method 
 * 
 * 
 */

/** 
 * 
 * @brief A main function to calculate physical quantities by full diagonalization method.
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @param neig number of eigenvalues
 * @version 0.2
 * @details add output process of calculation results for general spin
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void phys( //!<[inout]
          long int neig //!<[in]
) {
  long int i;
  double tmp_N;
#ifdef _SCALAPACK
  std::complex<double> *vec_tmp;
  int ictxt, ierr, rank;
  long int j, i_max;

  i_max = Check::idim_max;

  if(use_scalapack){
  fprintf(stdoutMPI, "In scalapack fulldiag, total spin is not calculated !\n");
  vec_tmp = malloc(i_max*sizeof(std::complex<double>));
  }
  for (i = 0; i < neig; i++) {
    for (j = 0; j < i_max; j++) {
      v0[j + 1] = 0.0;
    }
    if (use_scalapack) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);
      if (rank == 0) {
        for (j = 0; j < i_max; j++) {
          v0[j + 1] = vec_tmp[j];
        }
      }
      else {
        for (j = 0; j < i_max; j++) {
          v0[j + 1] = 0.0;
        }
      }
    }
    else {
      if (Def::iCalcType == FullDiag) {
        if (myrank == 0) {
          for (j = 0; j < i_max; j++) {
            v0[j + 1] = v1[i][j];
          }
        }
      }
      else {
        for (j = 0; j < i_max; j++) {
          v0[j + 1] = v1[i][j];
        }
      }
    }
  }/*for (i = 0; i < neig; i++)*/
#endif

  if (expec_energy_flct(neig, v0, v1) != 0) {
    fprintf(stderr, "Error: calc expec_energy.\n");
    exitMPI(-1);
  }
  if (expec_cisajs(neig, v0, v1) != 0) {
    fprintf(stderr, "Error: calc OneBodyG.\n");
    exitMPI(-1);
  }
  if (expec_cisajscktaltdc(neig, v0, v1) != 0) {
    fprintf(stderr, "Error: calc TwoBodyG.\n");
    exitMPI(-1);
  }
    
#ifdef _SCALAPACK
  if (use_scalapack) {
    if (Def::iCalcType == FullDiag) {
      Phys::s2 = 0.0;
      Phys::Sz = 0.0;
    }
  }
  else {
    if (Def::iCalcType == FullDiag) {
      if (expec_totalspin(v1) != 0) {
        fprintf(stderr, "Error: calc TotalSpin.\n");
        exitMPI(-1);
      }
    }
  }
#else
  if (Def::iCalcType == FullDiag) {
    if (expec_totalspin(neig, v1) != 0) {
      fprintf(stderr, "Error: calc TotalSpin.\n");
      exitMPI(-1);
    }
  }
#endif

  for (i = 0; i < neig; i++) {
    if (Def::iCalcModel == Spin || Def::iCalcModel == SpinGC) {
      tmp_N = Def::NsiteMPI;
    }
    else {
      tmp_N = Phys::num_up[i] + Phys::num_down[i];
    }
    if (Def::iCalcType == FullDiag) {
#ifdef _SCALAPACK
      if (use_scalapack) {
        fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf Doublon=%10lf \n", i, Phys::energy, tmp_N,
          Phys::Sz, Phys::doublon);
      }
      else {
        fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n", i, Phys::energy, tmp_N,
          Phys::Sz, Phys::s2, Phys::doublon);
      }
#else
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n",
        i, Phys::energy[i], tmp_N, Phys::Sz[i], Phys::s2[i], Phys::doublon[i]);
#endif      
    }
    else if (Def::iCalcType == CG)
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf Doublon=%10lf \n",
        i, Phys::energy[i], tmp_N, Phys::Sz[i], Phys::doublon[i]);
  }
#ifdef _SCALAPACK
  if(use_scalapack) free(vec_tmp);
#endif  
}
