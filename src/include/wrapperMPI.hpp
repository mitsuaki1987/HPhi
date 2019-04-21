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
@brief Global variables related to MPI and OpenMP
*/
#ifndef HPHI_WRAPPER_H
#define HPHI_WRAPPER_H
#include <complex>

void InitializeMPI(int argc, char *argv[]);
void FinalizeMPI();
void exitMPI(int errorcode);
FILE* fopenMPI(const char* FileName, const char* mode);
char* fgetsMPI(char* InputString, int maxcount,FILE* fp);
void BarrierMPI();
long int MaxMPI_li(long int idim);
double MaxMPI_d(double dvalue);
std::complex<double> SumMPI_dc(std::complex<double> norm);
double SumMPI_d(double norm);
void SumMPI_dv(int nnorm, double *norm);
void SumMPI_cv(int nnorm, std::complex<double> *norm);
long int SumMPI_li(long int idim);
int SumMPI_i(int idim);
long int BcastMPI_li(int root, long int idim);
double NormMPI_dc(long int idim, std::complex<double> *_v1);
void NormMPI_dv(long int ndim, int nstate, std::complex<double> **_v1, double *dnorm);
std::complex<double> VecProdMPI(long int ndim, std::complex<double> *v1, std::complex<double> *v2);
void MultiVecProdMPI(long int ndim, int nstate, std::complex<double> **v1, std::complex<double> **v2, std::complex<double> *prod);
void SendRecv_cv(int origin, long int nMsgS, long int nMsgR,
  std::complex<double> *vecs, std::complex<double> *vecr);
void SendRecv_iv(int origin, long int nMsgS, long int nMsgR,
  long int *vecs, long int *vecr);
long int SendRecv_i(int origin, long int isend);
#endif
