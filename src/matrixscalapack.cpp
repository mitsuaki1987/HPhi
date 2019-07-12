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
#include "matrixscalapack.hpp"
/**
 * @file
 * @version 3.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 * 
 * @brief File for diagonalization using scalapack with functions for getting indices of global and local array
 * 
 *
 */

#ifdef _SCALAPACK
int use_scalapack = 0;

/**
 * @brief compute block size for scalapack
 * @param[in] Msize size of matrix (Msize x Msize)
 * @param[in] MP::nproc number of processes
 * @return block size for scalapack
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetBlockSize(long int Msize, long int nproc) {
  long int block_size = 16;
  if (Msize * Msize / nproc > block_size * block_size)
    return block_size;
  return 1;
}
/**
 * @brief get processor array index from global array index
 * @param[in] i global array index
 * @param[in] np processor array dimension
 * @param[in] nb block size
 * @return processor array index
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetPArrayIndex(long int i, long int np, long int nb) {
  return (i/nb)%np;
}

/**
 * @brief get local array index from global array index
 * @param[in] i global array index
 * @param[in] np processor array dimension
 * @param[in] nb block size
 * @return local array index
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetLocalIndex(long int i, long int np, long int nb) {
  return (i/(np*nb))*nb + i%nb;
}

/**
 * @brief get global array index from local array index and processor array index
 * @param[in] il local array index
 * @param[in] p processor array index
 * @param[in] np processor array dimension
 * @param[in] nb block size
 * @return global array index
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetGlobalIndex(long int il, long int p, long int np, long int nb){
  return ((il/nb)*np+p)*nb + il%nb;
}

/**
 * @brief get rank of processor from indices of global matrix
 * @param[in] i index of global matrix
 * @param[in] j index of global matrix
 * @param[in] nprow processor array dimension for row
 * @param[in] npcol processor array dimension for column
 * @param[in] nb block size
 * @return rank of processor
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int MatToRank(long int i, long int j, long int nprow, long int npcol, long int nb){
  long int iproc, jproc;
  iproc = GetPArrayIndex(i, nprow, nb);
  jproc = GetPArrayIndex(j, npcol, nb);
  return iproc+jproc*nprow;
}

/**
 * @brief get column index of global matrix from given rank of processor and column index of local matrix
 * @param[in] lj column index of local matrix
 * @param[in] rank rank of processor
 * @param[in] npcol processor array dimension for column
 * @param[in] nb block size
 * @return column index of global matrix
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetMatRawInRank(long int lj, long int rank, long int npcol, long int nb){
  long int pcol = rank/npcol;
  return GetGlobalIndex(lj, pcol, npcol, nb);
}

/**
 * @brief get indices of local matrix from indices 
 (need to free memory after this function used)
 * @return indices array of local matrix
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int* GetMatElementInRank(
  long int i, /**<[in] index of local matrix*/
  long int j, /**<[in] index of local matrix*/
  long int nprow, /**<[in] processor array dimension for row*/
  long int npcol, /**<[in] processor array dimension for column*/
  long int nb/**<[in] block size*/
) {
  long int* ij;
  ij = (long int*)malloc(2 * sizeof(long int));
  ij[0] = GetLocalIndex(i, nprow, nb);
  ij[1] = GetLocalIndex(j, npcol, nb);
  return ij;
}
/**
 * @brief divide matrix
 * @param[in] m index of column of matrix
 * @param[in] n index of row of matrix
 * @param[in] Aorgmn value of matrix in (i, j)
 * @param[in, out] A distribution matrix
 * @param[in] desca descriptor of distribution matrix A
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
void DivMat(
  long int m,
  long int n,
  std::complex<double> Aorgmn,
  std::complex<double>* A,
  int* desca
) {
  long int mp = m + 1, np = n + 1;
  pzelset_(A, &mp, &np, desca, &Aorgmn);
}
/**
 * @brief get eigenvector from distributed matrix
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
void GetEigenVector(
  long int i,/**<[in] index of eigenvector*/
  long int m,/**<[in] size of eigenvector*/
  std::complex<double>* Z,/**<[in] distribution matrix of eigenvector*/
  int* descZ,/**<[in] descriptor for Z*/
  std::complex<double>* vec/**<[in, out] eigenvector*/
) {
  std::complex<double> alpha;
  long int j, ip, jp;
  int rank;
  char scope = 'A', top = ' ';

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ip = i + 1;
  for (j = 0; j < m; j++) {
    jp = j + 1;
    pzelget_(&scope, &top, &alpha, Z, &jp, &ip, descZ);
    if (rank == 0) {
      vec[j] = alpha;
    }
  }
}
/**
 * @brief diagonalization using scalapack
 * @return this returns 0 when it finished normally
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
int diag_scalapack_cmp(
  long int xNsize, /**<[in] size of matrix*/
  std::complex<double>** A,/**<[in] input matrix*/
  double* r, /**<[in, out] eigenvalue*/
  std::complex<double>* Z,/**<[in, out] distribution matrix of eigenvector*/
  int* descZ/**<[in, out] descriptor for Z*/
) {
  const int i_one = 1, i_zero = 0;
  const long int i_negone = -1;
  const double zero = 0.0, one = 1.0;
  long int m, n, mb, nb;
  int nprow, npcol;
  int myrow, mycol, info, lld;
  long int mp, nq;
  int ictxt;
  std::complex<double>* A_distr, * work, * rwork;
  double* W;
  int descA_distr[9];
  int rank, size, iam, nprocs;
  long int lwork, lrwork;
  int dims[2] = { 0,0 };
  long int i, j, ip, jp;
  char jobz = 'V', uplo = 'U', order = 'R';

  m = n = xNsize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, dims);
  nprow = dims[0]; npcol = dims[1];

  blacs_pinfo_(&iam, &nprocs);
  blacs_get_((int*)& i_negone, &i_zero, &ictxt);
  blacs_gridinit_(&ictxt, &order, &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

  mb = GetBlockSize(m, size);
  nb = GetBlockSize(n, size);

  mp = numroc_(&m, &mb, &myrow, &i_zero, &nprow);
  nq = numroc_(&n, &nb, &mycol, &i_zero, &npcol);
  W = (double*)malloc(n * sizeof(double));
  A_distr = (std::complex<double>*)malloc(mp * nq * sizeof(std::complex<double>));

  lld = (mp > 0) ? mp : 1;
  descinit_(descA_distr, &m, &n, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);
  descinit_(descZ, &m, &n, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      DivMat(i, j, A[i][j], A_distr, descA_distr);
    }
  }

  std::complex<double> wkopt, rwkopt;
  pzheev_(&jobz, &uplo, &n, A_distr, &i_one, &i_one, descA_distr, W, Z,
    &i_one, &i_one, descZ, &wkopt, &i_negone, &rwkopt, &i_negone, &info);

  lwork = (long int)wkopt.real();
  lrwork = (long int)rwkopt.real();
  work = (std::complex<double>*)malloc(lwork * sizeof(std::complex<double>));
  rwork = (std::complex<double>*)malloc(lrwork * sizeof(std::complex<double>));

  pzheev_(&jobz, &uplo, &n, A_distr, &i_one, &i_one, descA_distr, W, Z,
    &i_one, &i_one, descZ, work, &lwork, rwork, &lrwork, &info);

  if (rank == 0) {
    for (i = 0; i < n; i++) {
      r[i] = W[i];
    }
  }

  free(A_distr);
  free(work);
  free(rwork);
  free(W);

  use_scalapack = 1;

  return 0;
}

#endif
