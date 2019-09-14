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
#include "lapack_diag.hpp"
#include "matrixlapack.hpp"
#include "FileIO.hpp"
#ifdef _MAGMA
#include "matrixlapack_magma.hpp"
#endif
#ifdef _SCALAPACK
#include "matrixscalapack.hpp"
#endif
#include "global.hpp"
#include <cstring>

/** 
 * @brief performing full diagonalization using lapack
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int lapack_diag()
{
  FILE* fp;
  char sdt[D_FileNameMax] = "";
  long int i, j, i_max, xMsize;
#ifdef _SCALAPACK
  int rank, size, nprocs, nprow, npcol, myrow, mycol, ictxt;
  int i_negone = -1, i_zero = 0, iam;
  long int mb, nb, mp, nq;
  int dims[2] = { 0,0 };
  char order = 'R';
#endif

  i_max = Check::idim_max;

  xMsize = i_max;
  if (Def::iNGPU == 0) {
#ifdef _SCALAPACK
    if (MP::nproc > 1) {
      fprintf(MP::STDOUT, "Using SCALAPACK\n\n");
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Dims_create(size, 2, dims);
      nprow = dims[0]; npcol = dims[1];

      blacs_pinfo_(&iam, &nprocs);
      blacs_get_(&i_negone, &i_zero, &ictxt);
      blacs_gridinit_(&ictxt, &order, &nprow, &npcol);
      blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

      mb = GetBlockSize(xMsize, size);

      mp = numroc_(&xMsize, &mb, &myrow, &i_zero, &nprow);
      nq = numroc_(&xMsize, &mb, &mycol, &i_zero, &npcol);
      Z_vec = (std::complex<double>*)malloc(mp * nq * sizeof(std::complex<double>));
  
      diag_scalapack_cmp(xMsize, Wave::v0, Phys::energy, Z_vec, descZ_vec);
    }
    else {
      ZHEEVall(xMsize, Wave::v0, Phys::energy, Wave::v1);
    }
#else
    ZHEEVall(xMsize, Wave::v0, Phys::energy, Wave::v1);
#endif
  }
  else {
#ifdef _MAGMA
    if (MP::myrank == 0) {
      if (diag_magma_cmp(xMsize, Wave::v0, v0, v1, Def::iNGPU) != 0) {
        return -1;
      }
    }
#else
    fprintf(MP::STDOUT, "Warning: MAGMA is not used in this calculation.");
    ZHEEVall(xMsize, Wave::v0, Phys::energy, Wave::v1);
#endif
  }
  for (i = 0; i < i_max; i++) {
    for (j = 0; j < i_max; j++) {
      Wave::v0[i][j] = Wave::v1[i][j];
    }
  }

  strcpy(sdt, "Eigenvalue.dat");
  if (childfopenMPI(sdt, "w", &fp) != 0) {
    return -1;
  }
  for (i = 0; i < i_max; i++) {
    fprintf(fp, " %ld %.10lf \n", i, Phys::energy[i]);
  }
  fclose(fp);
  return 0;
}
