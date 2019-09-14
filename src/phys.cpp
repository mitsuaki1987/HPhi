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
#include "phys.hpp"
#include "expec_energy_flct.hpp"
#include "expec_totalspin.hpp"
#include "expec_cisajs.hpp"
#include "expec_cisajscktaltdc.hpp"
#include "wrapperMPI.hpp"
#ifdef _SCALAPACK
#include "matrixscalapack.hpp"
#endif
#include "global.hpp"
#include "DefCommon.hpp"
/** 
 * 
 * @brief A main function to calculate physical quantities by full diagonalization method.
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @version 0.2
 * @details add output process of calculation results for general spin
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void phys( //!<[inout]
  long int neig //!<[in] number of eigenvalues
) {
  long int i;
  double tmp_N;
#ifdef _SCALAPACK
  std::complex<double>* vec_tmp;
  int ictxt, ierr, rank;
  long int j, i_max;

  i_max = Check::idim_max;

  if (use_scalapack) {
    fprintf(MP::STDOUT, "In scalapack fulldiag, total spin is not calculated !\n");
    vec_tmp = (std::complex<double>*)malloc(i_max * sizeof(std::complex<double>));

    for (i = 0; i < neig; i++) {
      for (j = 0; j < i_max; j++) {
        Wave::v0[j][i] = 0.0;
      }
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);
      if (rank == 0) {
        for (j = 0; j < i_max; j++) {
          Wave::v0[j][i] = vec_tmp[j];
        }
      }
      else {
        for (j = 0; j < i_max; j++) {
          Wave::v0[j][i] = 0.0;
        }
      }
    }/*for (i = 0; i < neig; i++)*/
  }/*if (use_scalapack)*/
#endif

  if (expec::energy_flct::main(neig, Wave::v0, Wave::v1) != 0) {
    fprintf(stderr, "Error: calc expec_energy.\n");
    wrapperMPI::Exit(-1);
  }
  if (expec::cisajs::main(neig, Wave::v0, Wave::v1) != 0) {
    fprintf(stderr, "Error: calc OneBodyG.\n");
    wrapperMPI::Exit(-1);
  }
  if (expec::cisajscktalt::main(neig, Wave::v0, Wave::v1) != 0) {
    fprintf(stderr, "Error: calc TwoBodyG.\n");
    wrapperMPI::Exit(-1);
  }

  if (Def::iCalcType == DC::FullDiag) {
    if (expec_totalspin(neig, Wave::v1) != 0) {
      fprintf(stderr, "Error: calc TotalSpin.\n");
      wrapperMPI::Exit(-1);
    }
  }

  for (i = 0; i < neig; i++) {
    if (Def::iCalcModel == DC::Spin || Def::iCalcModel == DC::SpinGC) {
      tmp_N = Def::NsiteMPI;
    }
    else {
      tmp_N = Phys::num_up[i] + Phys::num_down[i];
    }
    if (Def::iCalcType == DC::FullDiag) {
#ifdef _SCALAPACK
      if (use_scalapack) {
        fprintf(MP::STDOUT, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf Doublon=%10lf \n", i, Phys::energy[i], tmp_N,
          Phys::Sz[i], Phys::doublon[i]);
      }
      else {
        fprintf(MP::STDOUT, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n", i, Phys::energy[i], tmp_N,
          Phys::Sz[i], Phys::s2[i], Phys::doublon[i]);
      }
#else
      fprintf(MP::STDOUT, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n",
        i, Phys::energy[i], tmp_N, Phys::Sz[i], Phys::s2[i], Phys::doublon[i]);
#endif      
    }
    else if (Def::iCalcType == DC::CG)
      fprintf(MP::STDOUT, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf Doublon=%10lf \n",
        i, Phys::energy[i], tmp_N, Phys::Sz[i], Phys::doublon[i]);
  }
#ifdef _SCALAPACK
  if (use_scalapack) free(vec_tmp);
#endif  
}
