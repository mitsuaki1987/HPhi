/*
HPhi-mVMC-StdFace - Common input generator
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
//
// Created by Kazuyoshi Yoshimi on 2019-01-09.
//

#include "setmemory.hpp"

///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int *i_1d_allocate(const long int N) {
  int *A;
  A = (int*)calloc((N), sizeof(int));
  return A;
}
///
/// \brief Function to free 1d array (int)
/// \param A Pointer of 1d array A
void free_i_1d_allocate(int *A) {
  free(A);
}
///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
long int *li_1d_allocate(const long int N) {
  long int *A;
  A = (long int*)calloc((N), sizeof(long int));
  return A;
}
///
/// \brief Function to free 1d array (int)
/// \param A Pointer of 1d array A
void free_li_1d_allocate(long int *A){
    free(A);
}
///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
long int **li_2d_allocate(const long int N, const long int M) {
  long int **A;
  long int int_i;
  A = (long int **)calloc((N), sizeof(long int *));
  A[0] = (long int *)calloc((M * N), sizeof(long int));
  for (int_i = 0; int_i < N; int_i++) {
    A[int_i] = A[0] + int_i * M;
  }
  return A;
}
///
/// \brief Function to free 2d array (int)
/// \param A Pointer of 2d array A
void free_li_2d_allocate(long int **A){
    free(A[0]);
    free(A);
}
///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int **i_2d_allocate(const long int N, const long int M) {
    int **A;
    long int int_i;
    A = (int **) calloc((N) , sizeof(int *));
    A[0] = (int *) calloc((M * N) , sizeof(int));
    for (int_i = 0; int_i < N; int_i++) {
        A[int_i] = A[0] + int_i * M;
    }
    //memset(A[0], 0, sizeof(int)*M*N);
    return A;
}
///
/// \brief Function to free 2d array (int)
/// \param A Pointer of 2d array A
void free_i_2d_allocate(int **A) {
  free(A[0]);
  free(A);
}
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int***i_3d_allocate(const long int N, const long int M, const long int L){
    long int int_i, int_j;
    int*** A;
    A     = (int***)calloc((N),sizeof(int**));
    A[0]  = (int**)calloc((M*N),sizeof(int*));
    A[0][0] = (int*)calloc((L*M*N),sizeof(int));
    for(int_i=0;int_i<N; int_i++) {
        A[int_i] = A[0] + int_i*M;
        for(int_j = 0; int_j<M; int_j++){
            A[int_i][int_j]= A[0][0] + int_i*M*L + int_j*L;
        }
    }
    return A;
}
///
/// \brief Function to free 3d array (int)
/// \param A A pointer of 3d array A
void free_i_3d_allocate(int ***A) {
  free(A[0][0]);
  free(A[0]);
  free(A);
}
///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double *d_1d_allocate(const long int N){
    double *A;
    A     = (double*)calloc((N),sizeof(double));
    return A;
}
///
/// \brief Function to free 1d array (double)
/// \param A Pointer of 1d array A
void free_d_1d_allocate(double *A) {
  free(A);
}
///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double **d_2d_allocate(const long int N, const long int M) {
  long int int_i;
  double **A;
  A = (double**)calloc((N), sizeof(double*));
  A[0] = (double*)calloc((M*N), sizeof(double));
  for (int_i = 0; int_i < N; int_i++) {
    A[int_i] = A[0] + int_i * M;
  }
  return A;
}
///
/// \brief Function to free 2d array (double)
/// \param A Pointer of 2d array A
void free_d_2d_allocate(double **A) {
  free(A[0]);
  free(A);
}
///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
std::complex<double> *cd_1d_allocate(const long int N){
    std::complex<double>*A;
    A     = (std::complex<double>*)calloc((N),sizeof(std::complex<double>));
    return A;
}
///
/// \brief Function to free 1d array (std::complex<double>)
/// \param A Pointer of 1d array A
void free_cd_1d_allocate(std::complex<double> *A) {
  free(A);
}
///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
std::complex<double> **cd_2d_allocate(const long int N, const long int M) {
  long int int_i;
  std::complex<double> **A;
  A = (std::complex<double>**)calloc((N), sizeof(std::complex<double>));
  A[0] = (std::complex<double>*)calloc((M*N), sizeof(std::complex<double>));
  for (int_i = 0; int_i < N; int_i++) {
    A[int_i] = A[0] + int_i * M;
  }
  return A;
}
///
/// \brief Function to free 2d array (std::complex<double>)
/// \param A Pointer of 2d array A
void free_cd_2d_allocate(std::complex<double>**A) {
  free(A[0]);
  free(A);
}
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
std::complex<double>***cd_3d_allocate(const long int N, const long int M, const long int L) {
  long int int_i, int_j;
  std::complex<double>***A;
  A = (std::complex<double>***)calloc((N), sizeof(std::complex<double>**));
  A[0] = (std::complex<double>**)calloc((M*N), sizeof(std::complex<double>*));
  A[0][0] = (std::complex<double>*)calloc((L*M*N), sizeof(std::complex<double>));
  for (int_i = 0; int_i < N; int_i++) {
    A[int_i] = A[0] + int_i * M;
    for (int_j = 0; int_j < M; int_j++) {
      A[int_i][int_j] = A[0][0] + int_i * M*L + int_j * L;
    }
  }
  return A;
}
///
/// \brief Function to free 3d array (std::complex<double>)
/// \param A A pointer of 3d array A
void free_cd_3d_allocate(std::complex<double>***A) {
  free(A[0][0]);
  free(A[0]);
  free(A);
}
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
std::complex<double>****cd_4d_allocate(const long int N, const long int M, const long int L, const long int K) {
  long int int_i, int_j, int_k;
  std::complex<double>****A;
  A = (std::complex<double>****)calloc((N), sizeof(std::complex<double>***));
  A[0] = (std::complex<double>***)calloc((M*N), sizeof(std::complex<double>**));
  A[0][0] = (std::complex<double>**)calloc((L*M*N), sizeof(std::complex<double>*));
  A[0][0][0] = (std::complex<double>*)calloc((K*L*M*N), sizeof(std::complex<double>));
  for (int_i = 0; int_i < N; int_i++) {
    A[int_i] = A[0] + int_i * M;
    for (int_j = 0; int_j < M; int_j++) {
      A[int_i][int_j] = A[0][0] + int_i * M*L + int_j * L;
      for (int_k = 0; int_k < L; int_k++) {
        A[int_i][int_j][int_k] = A[0][0][0] + int_i * M*L*K + int_j * L*K + int_k * K;
      }
    }
  }
  return A;
}
///
/// \brief Function to free 3d array (std::complex<double>)
/// \param A A pointer of 3d array A
void free_cd_4d_allocate(std::complex<double>****A) {
  free(A[0][0][0]);
  free(A[0][0]);
  free(A[0]);
  free(A);
}
