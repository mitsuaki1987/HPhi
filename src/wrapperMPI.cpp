/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**@file
@brief MPI wrapper for init, finalize, bcast, etc.
*/
#ifdef __MPI
#include <mpi.h>
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "wrapperMPI.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <complex>
#include "splash.hpp"
#include "global.hpp"
#include "common/setmemory.hpp"

/**
@brief MPI initialization wrapper
Process ID (::MP::myrank), Number of processes (::MP::nproc), 
Number of threads (::MP::nthreads), and pointer to the standard output
(::MP::STDOUT) are specified here.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void wrapperMPI::Initialize(int argc, char *argv[]){
#ifdef __MPI
  int ierr;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &MP::nproc);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MP::myrank);
  if(ierr != 0) wrapperMPI::Exit(ierr);
#else
  MP::nproc = 1;
  MP::myrank = 0;
#endif
  if (MP::myrank == 0) MP::STDOUT = stdout;
  else MP::STDOUT = fopen("/dev/null", "w");
  splash();

#pragma omp parallel default(none) shared(MP::nthreads)
#pragma omp master
#ifdef _OPENMP
  MP::nthreads = omp_get_num_threads();
#else
  MP::nthreads=1;
#endif
  fprintf(MP::STDOUT, "\n\n#####  Parallelization Info.  #####\n\n");
  fprintf(MP::STDOUT, "  OpenMP threads : %d\n", MP::nthreads);
  fprintf(MP::STDOUT, "  MPI PEs : %d \n\n", MP::nproc);
}/*void wrapperMPI::Initialize(int argc, char *argv[])*/
/**
@brief MPI Finitialization wrapper
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void wrapperMPI::Finalize(){
#ifdef __MPI
  int ierr;
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  if (MP::myrank != 0) fclose(MP::STDOUT);
}
/**
@brief MPI Abortation wrapper
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void wrapperMPI::Exit(
  int errorcode//!<[in] Error-code to be returned as that of this program
)
{
  fflush(stdout);
#ifdef __MPI
  fprintf(stdout,"\n\n #######  [HPhi] You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
  int ierr;
  ierr = MPI_Abort(MPI_COMM_WORLD, errorcode);
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  exit(errorcode);
}/*void wrapperMPI::Exit*/
/**
@brief MPI file I/O (open) wrapper.
Only the root node (::MP::myrank = 0) should be open/read/write (small) parameter files.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
FILE* wrapperMPI::Fopen(
  const char* FileName,//!<[in] Input/output file
  const char* mode//!<[in] "w", "r", etc.
){
  FILE* fp;

  if (MP::myrank == 0) fp = fopen(FileName, mode);
  else fp = fopen("/dev/null", "w");

  return fp;
}/*FILE* wrapperMPI::Fopen*/
/**
@brief MPI file I/O (get a line, fgets) wrapper.
Only the root node (::MP::myrank = 0) reads and broadcast string.
@return The same as that of fgets
@author Mitsuaki Kawamura (The University of Tokyo)
*/
char* wrapperMPI::Fgets(
  char* InputString,//!<[out] read line.
  int maxcount,//!<[in] Length of string
  FILE* fp//!<[in] file pointer
){
  int inull;
  char *ctmp;

  ctmp = InputString;
  inull = 0;
  if (MP::myrank == 0) {
    ctmp = fgets(InputString, maxcount, fp);
    if (ctmp == NULL){
      inull = 1;
    }
    
    while(*InputString == '\n' || strncmp(InputString, "#", 1)==0){
      ctmp = fgets(InputString, maxcount, fp);
      if (ctmp == NULL){
        inull=1;
        break;
      }
    }
  }
#ifdef __MPI
  MPI_Bcast(InputString, maxcount, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&inull, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (MP::myrank != 0 && inull == 1) {
    ctmp = NULL;
  }

  return ctmp;
}/*char* wrapperMPI::Fgets*/
/**
@brief MPI barrier wrapper.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void wrapperMPI::Barrier(){
#ifdef __MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}/*void wrapperMPI::Barrier()*/
/**
@brief MPI wrapper function to obtain maximum unsigned
long integer across processes.
@return Maximum value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
long int wrapperMPI::Max_li(
  long int idim//!<[in] Value to be maximized
){
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  if(ierr != 0) wrapperMPI::Exit(-1);
#endif
  return(idim);
}/*long int wrapperMPI::Max_li*/
/**
@brief MPI wrapper function to obtain maximum Double
across processes.
@return Maximum value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double wrapperMPI::Max_d(
  double dvalue//!<[in] Value to be maximized
){
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &dvalue, 1,
    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(ierr != 0) wrapperMPI::Exit(-1);
#endif
  return(dvalue);
}/*double wrapperMPI::Max_d*/
/**
@brief MPI wrapper function to obtain sum of Double
complex across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
std::complex<double> wrapperMPI::Sum_dc(
  std::complex<double> norm//!<[in] Value to be summed
){
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) wrapperMPI::Exit(-1);
#endif
  return(norm);
}/*std::complex<double> wrapperMPI::Sum_dc*/
/**
@brief MPI wrapper function to obtain sum of Double
across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double wrapperMPI::Sum_d(
  double norm//!<[in] Value to be summed
){
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) wrapperMPI::Exit(-1);
#endif
  return(norm);
}/*double wrapperMPI::Sum_d*/
/**
@brief MPI wrapper function to obtain sum of Double array
across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void wrapperMPI::Sum_dv(
  int nnorm,
  double *norm//!<[in] Value to be summed
) {
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, norm, nnorm,
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
  if (ierr != 0) wrapperMPI::Exit(-1);
#endif
}/*void wrapperMPI::Sum_dv*/
/**
@brief MPI wrapper function to obtain sum of Double array
across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void wrapperMPI::Sum_cv(
  int nnorm,
  std::complex<double> *norm//!<[in] Value to be summed
) {
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, norm, nnorm,
    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  if (ierr != 0) wrapperMPI::Exit(-1);
#endif
}/*void wrapperMPI::Sum_cv*/
/**
@brief MPI wrapper function to obtain sum of unsigned
long integer across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
long int wrapperMPI::Sum_li(
  long int idim//!<[in] Value to be summed
){
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) wrapperMPI::Exit(-1);
#endif
  return(idim);
}/*long int wrapperMPI::Sum_li*/
/**
@brief MPI wrapper function to obtain sum of
integer across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int wrapperMPI::Sum_i(
  int idim//!<[in] Value to be summed
) {
#ifdef __MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
                       MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) wrapperMPI::Exit(-1);
#endif
  return(idim);
}/*int wrapperMPI::Sum_i*/
/**
@brief MPI wrapper function to broadcast long
integer across processes.
@return Broadcasted value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
long int wrapperMPI::Bcast_li(
  int root,//!<[in] The source process of the broadcast
  long int idim//!<[in] Value to be broadcasted
) {
  long int idim0;
  idim0 = idim;
#ifdef __MPI
    MPI_Bcast(&idim0, 1, MPI_LONG, root, MPI_COMM_WORLD);
#endif
  return(idim0);
}/*long int wrapperMPI::Bcast_li*/
/**
@brief Compute norm of process-distributed vector
@f$|{\bf v}_1|^2@f$
@return Norm @f$|{\bf v}_1|^2@f$
*/
double wrapperMPI::Norm_dc(
  long int idim,//!<[in] Local dimension of vector
  std::complex<double> *_v1//!<[in] [idim] vector to be producted
){
  double dnorm =0;
  long int i;

  dnorm = 0.0;
#pragma omp parallel for default(none) private(i) \
shared(_v1, idim) reduction(+:dnorm)
    for (i = 1; i <= idim; i++) 
      dnorm += real(conj(_v1[i])*_v1[i]);
 
#ifdef __MPI
  dnorm = wrapperMPI::Sum_d(dnorm);
#endif
  return dnorm;
}/*double wrapperMPI::Norm_dc*/
/**
@brief Compute norm of process-distributed vector
@f$|{\bf v}_1|^2@f$
@return Norm @f$|{\bf v}_1|^2@f$
*/
void wrapperMPI::Norm_dv(
  long int ndim,//!<[in] Local dimension of vector
  int nstate,
  std::complex<double> **_v1,//!<[in] [idim] vector to be producted
  double *dnorm
) {
  long int idim;
  int istate;

  for (istate = 0; istate < nstate; istate++) dnorm[istate] = 0.0;
  for (idim = 1; idim <= ndim; idim++) {
    for (istate = 0; istate < nstate; istate++) {
      dnorm[istate] += real(conj(_v1[idim][istate])*_v1[idim][istate]);
    }
  }
  wrapperMPI::Sum_dv(nstate, dnorm);
  for (istate = 0; istate < nstate; istate++) dnorm[istate] = sqrt(dnorm[istate]);
}/*double NormMPI_cv*/
/**
@brief Compute conjugate scaler product of process-distributed vector
@f${\bf v}_1^* \cdot {\bf v}_2@f$
@return Conjugate scaler product @f${\bf v}_1^* \cdot {\bf v}_2@f$
*/
std::complex<double> wrapperMPI::VecProd(
  long int ndim,//!<[in] Local dimension of vector
  std::complex<double> *v1,//!<[in] [ndim] vector to be producted
  std::complex<double> *v2//!<[in] [ndim] vector to be producted
){
  long int idim;
  std::complex<double> prod, *prod_thr;
  int mythread;

  prod_thr = cd_1d_allocate(MP::nthreads);
#pragma omp parallel default(none) shared(v1,v2,ndim,prod,prod_thr) private(idim,mythread)
  {
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
#pragma omp for
    for (idim = 1; idim <= ndim; idim++) 
      prod_thr[mythread] += conj(v1[idim]) * v2[idim];
  }
  prod = 0.0;
  for (mythread = 0; mythread < MP::nthreads; mythread++)
    prod += prod_thr[mythread];
  free_cd_1d_allocate(prod_thr);

  prod = wrapperMPI::Sum_dc(prod);

  return(prod);
}/*std::complex<double> wrapperMPI::VecProd*/
/**
@brief Compute conjugate scaler product of process-distributed vector
@f${\bf v}_1^* \cdot {\bf v}_2@f$
*/
void wrapperMPI::MultiVecProd(
  long int ndim,//!<[in] Local dimension of vector
  int nstate,
  std::complex<double> **v1,//!<[in] [ndim] vector to be producted
  std::complex<double> **v2,//!<[in] [ndim] vector to be producted
  std::complex<double> *prod
) {
  long int idim;
  int istate;

  for (istate = 0; istate < nstate; istate++) prod[istate] = 0.0;
  for (idim = 1; idim <= ndim; idim++) {
    for (istate = 0; istate < nstate; istate++) {
      prod[istate] += conj(v1[idim][istate])*v2[idim][istate];
    }
  }
  wrapperMPI::Sum_cv(nstate, prod);
}/*void wrapperMPI::MultiVecProd*/
/**
@brief Wrapper of MPI_Sendrecv for std::complex<double> number.
When we pass a message longer than 2^31-1 
(max of int: 2147483647), we need to divide it.
*/
void wrapperMPI::SendRecv_cv(
  int origin,
  long int nMsgS,
  long int nMsgR,
  std::complex<double> *vecs,
  std::complex<double> *vecr
) {
#ifdef __MPI
  int ierr, two31m1 = 2147483647, modMsg, nMsgS2, nMsgR2;
  long int nMsg, nnMsg, iMsg, sMsgR, sMsgS;
  MPI_Status statusMPI;

  if (nMsgS > nMsgR) nMsg = nMsgS;
  else nMsg = nMsgR;
  nnMsg = nMsg / two31m1;
  modMsg = nMsg % two31m1;
  if (modMsg != 0) nnMsg += 1;

  sMsgS = 0;
  sMsgR = 0;
  for (iMsg = 0; iMsg < nnMsg; iMsg++) {
    nMsgS2 = nMsgS / nnMsg;
    nMsgR2 = nMsgR / nnMsg;
    if (iMsg < nMsgS % nnMsg) nMsgS2 += 1;
    if (iMsg < nMsgR % nnMsg) nMsgR2 += 1;

    ierr = MPI_Sendrecv(&vecs[sMsgS], nMsgS2, MPI_DOUBLE_COMPLEX, origin, 0,
                        &vecr[sMsgR], nMsgR2, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) wrapperMPI::Exit(-1);

    sMsgS += nMsgS2;
    sMsgR += nMsgR2;
  }
#endif
}/*void wrapperMPI::SendRecv_cv*/
/**
@brief Wrapper of MPI_Sendrecv for long integer number.
When we pass a message longer than 2^31-1
(max of int: 2147483647), we need to divide it.
*/
void wrapperMPI::SendRecv_iv(
  int origin,
  long int nMsgS,
  long int nMsgR,
  long int *vecs,
  long int *vecr
) {
#ifdef __MPI
  int ierr, two31m1 = 2147483647, modMsg, nMsgS2, nMsgR2;
  long int nMsg, nnMsg, iMsg, sMsgR, sMsgS;
  MPI_Status statusMPI;

  if (nMsgS > nMsgR) nMsg = nMsgS;
  else nMsg = nMsgR;
  nnMsg = nMsg / two31m1;
  modMsg = nMsg % two31m1;
  if (modMsg != 0) nnMsg += 1;

  sMsgS = 0;
  sMsgR = 0;
  for (iMsg = 0; iMsg < nnMsg; iMsg++) {
    nMsgS2 = nMsgS / nnMsg;
    nMsgR2 = nMsgR / nnMsg;
    if (iMsg < nMsgS % nnMsg) nMsgS2 += 1;
    if (iMsg < nMsgR % nnMsg) nMsgR2 += 1;

    ierr = MPI_Sendrecv(&vecs[sMsgS], nMsgS2, MPI_LONG, origin, 0,
                        &vecr[sMsgR], nMsgR2, MPI_LONG, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) wrapperMPI::Exit(-1);

    sMsgS += nMsgS2;
    sMsgR += nMsgR2;
  }
#endif
}/*void wrapperMPI::SendRecv_iv*/
/**
@brief Wrapper of MPI_Sendrecv for long integer number.
*/
long int wrapperMPI::SendRecv_i(
  int origin,
  long int isend
) {
#ifdef __MPI
  int ierr;
  MPI_Status statusMPI;
  long int ircv;
  ierr = MPI_Sendrecv(&isend, 1, MPI_LONG, origin, 0,
                      &ircv,  1, MPI_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) wrapperMPI::Exit(ierr);
  return ircv;
#else
  return isend;
#endif
}/*void wrapperMPI::SendRecv_i*/
