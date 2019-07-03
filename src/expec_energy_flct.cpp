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

#include "bitcalc.hpp"
#include "mltplyCommon.hpp"
#include "mltply.hpp"
#include "expec_energy_flct.hpp"
#include "wrapperMPI.hpp"
#include "CalcTime.hpp"
#include "common/setmemory.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include "log.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
///
/// \brief Calculate expected values of energies and physical quantities for Hubbard GC model
/// \param X [in, out] X Struct to get information about file header names, dimension of hirbert space, calc type and output physical quantities.
/// \retval 0 normally finished.
/// \retval -1 abnormally finished.
int expec_energy_flct_HubbardGC(
  int nstate,
  std::complex<double> **tmp_v0
) {
  long int j;
  long int isite1;
  long int is1_up_a, is1_up_b;
  long int is1_down_a, is1_down_b;
  int bit_up, bit_down, bit_D, istate, mythread;
  long int ibit_up, ibit_down, ibit_D;
  double D, N, Sz;
  double *tmp_v02;
  long int i_max;
  int l_ibit1, u_ibit1, i_32;
  double **doublon_t, **doublon2_t, **num_t, **num2_t, **Sz_t, **Sz2_t;

  i_max = Check::idim_max;
  doublon_t = d_2d_allocate(nthreads, nstate);
  doublon2_t = d_2d_allocate(nthreads, nstate);
  num_t = d_2d_allocate(nthreads, nstate);
  num2_t = d_2d_allocate(nthreads, nstate);
  Sz_t = d_2d_allocate(nthreads, nstate);
  Sz2_t = d_2d_allocate(nthreads, nstate);
  i_32 = 0xFFFFFFFF; //2^32 - 1

  //[s] for bit count
  is1_up_a = 0;
  is1_up_b = 0;
  is1_down_a = 0;
  is1_down_b = 0;
  for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
    if (isite1 > Def::Nsite) {
      is1_up_a += Def::Tpow[2 * isite1 - 2];
      is1_down_a += Def::Tpow[2 * isite1 - 1];
    }
    else {
      is1_up_b += Def::Tpow[2 * isite1 - 2];
      is1_down_b += Def::Tpow[2 * isite1 - 1];
    }
  }
  //[e]
#pragma omp parallel default(none) \
private(j,tmp_v02,D,N,Sz,isite1,bit_up,bit_down,bit_D, \
        u_ibit1,l_ibit1, ibit_up,ibit_down,ibit_D,mythread,istate) \
shared(tmp_v0,list_1,doublon_t,doublon2_t,num_t,num2_t,Sz_t,Sz2_t,nstate, \
       i_max,myrank,is1_up_a,is1_down_a,is1_up_b,is1_down_b,i_32) \

  {
    tmp_v02 = d_1d_allocate(nstate);
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      for (istate = 0; istate < nstate; istate++)
        tmp_v02[istate] = real(conj(tmp_v0[j][istate]) * tmp_v0[j][istate]);
      bit_up = 0;
      bit_down = 0;
      bit_D = 0;
      // isite1 > Def::Nsite
      ibit_up = (long int) myrank & is1_up_a;
      u_ibit1 = ibit_up >> 32;
      l_ibit1 = ibit_up & i_32;
      bit_up += pop(u_ibit1);
      bit_up += pop(l_ibit1);

      ibit_down = (long int) myrank & is1_down_a;
      u_ibit1 = ibit_down >> 32;
      l_ibit1 = ibit_down & i_32;
      bit_down += pop(u_ibit1);
      bit_down += pop(l_ibit1);

      ibit_D = (ibit_up) & (ibit_down >> 1);
      u_ibit1 = ibit_D >> 32;
      l_ibit1 = ibit_D & i_32;
      bit_D += pop(u_ibit1);
      bit_D += pop(l_ibit1);

      // isite1 <= Def::Nsite
      ibit_up = (long int) (j - 1) & is1_up_b;
      u_ibit1 = ibit_up >> 32;
      l_ibit1 = ibit_up & i_32;
      bit_up += pop(u_ibit1);
      bit_up += pop(l_ibit1);

      ibit_down = (long int) (j - 1) & is1_down_b;
      u_ibit1 = ibit_down >> 32;
      l_ibit1 = ibit_down & i_32;
      bit_down += pop(u_ibit1);
      bit_down += pop(l_ibit1);

      ibit_D = (ibit_up) & (ibit_down >> 1);
      u_ibit1 = ibit_D >> 32;
      l_ibit1 = ibit_D & i_32;
      bit_D += pop(u_ibit1);
      bit_D += pop(l_ibit1);

      D = bit_D;
      N = bit_up + bit_down;
      Sz = bit_up - bit_down;

      for (istate = 0; istate < nstate; istate++) {
        doublon_t[mythread][istate] += tmp_v02[istate] * D;
        doublon2_t[mythread][istate] += tmp_v02[istate] * D * D;
        num_t[mythread][istate] += tmp_v02[istate] * N;
        num2_t[mythread][istate] += tmp_v02[istate] * N * N;
        Sz_t[mythread][istate] += tmp_v02[istate] * Sz;
        Sz2_t[mythread][istate] += tmp_v02[istate] * Sz * Sz;
      }
    }/*for (j = 1; j <= i_max; j++)*/
    free_d_1d_allocate(tmp_v02);
  }/*end of parallel region*/

  for (istate = 0; istate < nstate; istate++) {
    Phys::doublon[istate] = 0.0;
    Phys::doublon2[istate] = 0.0;
    Phys::num[istate] = 0.0;
    Phys::num2[istate] = 0.0;
    Phys::Sz[istate] = 0.0;
    Phys::Sz2[istate] = 0.0;
    for (mythread = 0; mythread < nthreads; mythread++) {
      Phys::doublon[istate] += doublon_t[mythread][istate];
      Phys::doublon2[istate] += doublon2_t[mythread][istate];
      Phys::num[istate] += num_t[mythread][istate];
      Phys::num2[istate] += num2_t[mythread][istate];
      Phys::Sz[istate] += Sz_t[mythread][istate];
      Phys::Sz2[istate] += Sz2_t[mythread][istate];
    }
  }
  SumMPI_dv(nstate, Phys::doublon);
  SumMPI_dv(nstate, Phys::doublon2);
  SumMPI_dv(nstate, Phys::num);
  SumMPI_dv(nstate, Phys::num2);
  SumMPI_dv(nstate, Phys::Sz);
  SumMPI_dv(nstate, Phys::Sz2);

  for (istate = 0; istate < nstate; istate++) {
    Phys::Sz[istate] *= 0.5;
    Phys::Sz2[istate] *= 0.25;
    Phys::num_up[istate] = 0.5*(Phys::num[istate] + Phys::Sz[istate]);
    Phys::num_down[istate] = 0.5*(Phys::num[istate] - Phys::Sz[istate]);
  }

  free_d_2d_allocate(doublon_t);
  free_d_2d_allocate(doublon2_t);
  free_d_2d_allocate(num_t);
  free_d_2d_allocate(num2_t);
  free_d_2d_allocate(Sz_t);
  free_d_2d_allocate(Sz2_t);
  return 0;
}
///
/// \brief Calculate expected values of energies and physical quantities for Hubbard model
/// \param X [in, out] X Struct to get information about file header names, dimension of hirbert space, calc type and output physical quantities.
/// \retval 0 normally finished.
/// \retval -1 abnormally finished.
int expec_energy_flct_Hubbard(
  
  int nstate,
  std::complex<double> **tmp_v0
) {
  long int j;
  long int isite1;
  long int is1_up_a, is1_up_b;
  long int is1_down_a, is1_down_b;
  int bit_up, bit_down, bit_D, istate, mythread;
  long int ibit_up, ibit_down, ibit_D;
  double **doublon_t, **doublon2_t, **num_t, **num2_t, **Sz_t, **Sz2_t;
  double D, N, Sz;
  double *tmp_v02;
  long int i_max, tmp_list_1;
  int l_ibit1, u_ibit1, i_32;

  i_max = Check::idim_max;

  doublon_t = d_2d_allocate(nthreads, nstate);
  doublon2_t = d_2d_allocate(nthreads, nstate);
  num_t = d_2d_allocate(nthreads, nstate);
  num2_t = d_2d_allocate(nthreads, nstate);
  Sz_t = d_2d_allocate(nthreads, nstate);
  Sz2_t = d_2d_allocate(nthreads, nstate);
  i_32 = 0xFFFFFFFF;

  //[s] for bit count
  is1_up_a = 0;
  is1_up_b = 0;
  is1_down_a = 0;
  is1_down_b = 0;
  for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
    if (isite1 > Def::Nsite) {
      is1_up_a += Def::Tpow[2 * isite1 - 2];
      is1_down_a += Def::Tpow[2 * isite1 - 1];
    }
    else {
      is1_up_b += Def::Tpow[2 * isite1 - 2];
      is1_down_b += Def::Tpow[2 * isite1 - 1];
    }
  }
  //[e]
#pragma omp parallel default(none) \
shared(tmp_v0,list_1,doublon_t,doublon2_t,num_t,num2_t,Sz_t,Sz2_t,nstate, \
       i_max, myrank,is1_up_a,is1_down_a,is1_up_b,is1_down_b,i_32) \
private(j,tmp_v02,D,N,Sz,isite1,tmp_list_1,bit_up,bit_down,bit_D,u_ibit1, \
        l_ibit1,ibit_up,ibit_down,ibit_D,mythread,istate)
  {
    tmp_v02 = d_1d_allocate(nstate);
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      for (istate = 0; istate < nstate; istate++)
        tmp_v02[istate] = real(conj(tmp_v0[j][istate]) * tmp_v0[j][istate]);
      bit_up = 0;
      bit_down = 0;
      bit_D = 0;
      tmp_list_1 = list_1[j];
      // isite1 > Def::Nsite
      ibit_up = (long int) myrank & is1_up_a;
      u_ibit1 = ibit_up >> 32;
      l_ibit1 = ibit_up & i_32;
      bit_up += pop(u_ibit1);
      bit_up += pop(l_ibit1);

      ibit_down = (long int) myrank & is1_down_a;
      u_ibit1 = ibit_down >> 32;
      l_ibit1 = ibit_down & i_32;
      bit_down += pop(u_ibit1);
      bit_down += pop(l_ibit1);

      ibit_D = (ibit_up) & (ibit_down >> 1);
      u_ibit1 = ibit_D >> 32;
      l_ibit1 = ibit_D & i_32;
      bit_D += pop(u_ibit1);
      bit_D += pop(l_ibit1);

      // isite1 <= Def::Nsite
      ibit_up = (long int) tmp_list_1 & is1_up_b;
      u_ibit1 = ibit_up >> 32;
      l_ibit1 = ibit_up & i_32;
      bit_up += pop(u_ibit1);
      bit_up += pop(l_ibit1);

      ibit_down = (long int) tmp_list_1 & is1_down_b;
      u_ibit1 = ibit_down >> 32;
      l_ibit1 = ibit_down & i_32;
      bit_down += pop(u_ibit1);
      bit_down += pop(l_ibit1);

      ibit_D = (ibit_up) & (ibit_down >> 1);
      u_ibit1 = ibit_D >> 32;
      l_ibit1 = ibit_D & i_32;
      bit_D += pop(u_ibit1);
      bit_D += pop(l_ibit1);

      D = bit_D;
      N = bit_up + bit_down;
      Sz = bit_up - bit_down;

      for (istate = 0; istate < nstate; istate++) {
        doublon_t[mythread][istate] += tmp_v02[istate] * D;
        doublon2_t[mythread][istate] += tmp_v02[istate] * D * D;
        num_t[mythread][istate] += tmp_v02[istate] * N;
        num2_t[mythread][istate] += tmp_v02[istate] * N * N;
        Sz_t[mythread][istate] += tmp_v02[istate] * Sz;
        Sz2_t[mythread][istate] += tmp_v02[istate] * Sz * Sz;
      }
    }/*for (j = 1; j <= i_max; j++)*/
    free_d_1d_allocate(tmp_v02);
  }/*end of parallel region*/

  for (istate = 0; istate < nstate; istate++) {
    Phys::doublon[istate] = 0.0;
    Phys::doublon2[istate] = 0.0;
    Phys::num[istate] = 0.0;
    Phys::num2[istate] = 0.0;
    Phys::Sz[istate] = 0.0;
    Phys::Sz2[istate] = 0.0;
    for (mythread = 0; mythread < nthreads; mythread++) {
      Phys::doublon[istate] += doublon_t[mythread][istate];
      Phys::doublon2[istate] += doublon2_t[mythread][istate];
      Phys::num[istate] += num_t[mythread][istate];
      Phys::num2[istate] += num2_t[mythread][istate];
      Phys::Sz[istate] += Sz_t[mythread][istate];
      Phys::Sz2[istate] += Sz2_t[mythread][istate];
    }
  }
  SumMPI_dv(nstate, Phys::doublon);
  SumMPI_dv(nstate, Phys::doublon2);
  SumMPI_dv(nstate, Phys::num);
  SumMPI_dv(nstate, Phys::num2);
  SumMPI_dv(nstate, Phys::Sz);
  SumMPI_dv(nstate, Phys::Sz2);

  for (istate = 0; istate < nstate; istate++) {
    Phys::Sz[istate] *= 0.5;
    Phys::Sz2[istate] *= 0.25;
    Phys::num_up[istate] = 0.5*(Phys::num[istate] + Phys::Sz[istate]);
    Phys::num_down[istate] = 0.5*(Phys::num[istate] - Phys::Sz[istate]);
  }

  free_d_2d_allocate(doublon_t);
  free_d_2d_allocate(doublon2_t);
  free_d_2d_allocate(num_t);
  free_d_2d_allocate(num2_t);
  free_d_2d_allocate(Sz_t);
  free_d_2d_allocate(Sz2_t);
  return 0;
}
///
/// \brief Calculate expected values of energies and physical quantities for Half-SpinGC model
/// \param X [in, out] X Struct to get information about file header names, dimension of hirbert space, calc type and output physical quantities.
/// \retval 0 normally finished.
/// \retval -1 abnormally finished.
int expec_energy_flct_HalfSpinGC(
  
  int nstate,
  std::complex<double> **tmp_v0
) {
  long int j;
  long int isite1;
  long int is1_up_a, is1_up_b;

  long int ibit1;
  double Sz;
  double *tmp_v02;
  long int i_max;
  int l_ibit1, u_ibit1, i_32;
  int istate, mythread;
  double **Sz_t, **Sz2_t;

  i_max = Check::idim_max;

  Sz_t = d_2d_allocate(nthreads, nstate);
  Sz2_t = d_2d_allocate(nthreads, nstate);
  i_32 = 0xFFFFFFFF; //2^32 - 1

  //[s] for bit count
  is1_up_a = 0;
  is1_up_b = 0;
  for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
    if (isite1 > Def::Nsite) {
      is1_up_a += Def::Tpow[isite1 - 1];
    }
    else {
      is1_up_b += Def::Tpow[isite1 - 1];
    }
  }
  //[e]
#pragma omp parallel default(none) \
private(j,Sz,ibit1,isite1,tmp_v02,u_ibit1,l_ibit1,mythread,istate) \
shared(tmp_v0,Sz_t,Sz2_t,nstate, i_max,myrank,i_32,is1_up_a,is1_up_b,Def::NsiteMPI)

  {
    tmp_v02 = d_1d_allocate(nstate);
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      for (istate = 0; istate < nstate; istate++)
        tmp_v02[istate] = real(conj(tmp_v0[j][istate]) * tmp_v0[j][istate]);
      Sz = 0.0;

      // isite1 > Def::Nsite
      ibit1 = (long int) myrank & is1_up_a;
      u_ibit1 = ibit1 >> 32;
      l_ibit1 = ibit1 & i_32;
      Sz += pop(u_ibit1);
      Sz += pop(l_ibit1);
      // isite1 <= Def::Nsite
      ibit1 = (long int) (j - 1)&is1_up_b;
      u_ibit1 = ibit1 >> 32;
      l_ibit1 = ibit1 & i_32;
      Sz += pop(u_ibit1);
      Sz += pop(l_ibit1);
      Sz = 2 * Sz - Def::NsiteMPI;

      for (istate = 0; istate < nstate; istate++) {
        Sz_t[mythread][istate] += tmp_v02[istate] * Sz;
        Sz2_t[mythread][istate] += tmp_v02[istate] * Sz * Sz;
      }
    }/*for (j = 1; j <= i_max; j++)*/
    free_d_1d_allocate(tmp_v02);
  }/*End of parallel region*/
  for (istate = 0; istate < nstate; istate++) {
    Phys::Sz[istate] = 0.0;
    Phys::Sz2[istate] = 0.0;
    for (mythread = 0; mythread < nthreads; mythread++) {
      Phys::Sz[istate] += Sz_t[mythread][istate];
      Phys::Sz2[istate] += Sz2_t[mythread][istate];
    }
  }
  SumMPI_dv(nstate, Phys::Sz);
  SumMPI_dv(nstate, Phys::Sz2);

  for (istate = 0; istate < nstate; istate++) {
    Phys::doublon[istate] = 0.0;
    Phys::doublon2[istate] = 0.0;
    Phys::num[istate] = Def::NsiteMPI;
    Phys::num2[istate] = Def::NsiteMPI*Def::NsiteMPI;
    Phys::Sz[istate] *= 0.5;
    Phys::Sz2[istate] *= 0.25;
    Phys::num_up[istate] = 0.5*(Phys::num[istate] + Phys::Sz[istate]);
    Phys::num_down[istate] = 0.5*(Phys::num[istate] - Phys::Sz[istate]);
  }

  free_d_2d_allocate(Sz_t);
  free_d_2d_allocate(Sz2_t);
  return 0;
}
///
/// \brief Calculate expected values of energies and physical quantities for General-SpinGC model
/// \param X [in, out] X Struct to get information about file header names, dimension of hirbert space, calc type and output physical quantities.
/// \retval 0 normally finished.
/// \retval -1 abnormally finished.
int expec_energy_flct_GeneralSpinGC(
  
  int nstate,
  std::complex<double> **tmp_v0
) {
  long int j;
  long int isite1;
  int istate, mythread;
  double Sz;
  double *tmp_v02;
  long int i_max;
  double **Sz_t, **Sz2_t;

  Sz_t = d_2d_allocate(nthreads, nstate);
  Sz2_t = d_2d_allocate(nthreads, nstate);
  i_max = Check::idim_max;

#pragma omp parallel default(none) \
private(j,istate,tmp_v02,Sz,isite1,mythread)\
shared(i_max,nstate,myrank,Sz_t,Sz2_t,tmp_v0, \
Def::SiteToBit, Def::Tpow,Def::Nsite,Def::NsiteMPI)
  {
    tmp_v02 = d_1d_allocate(nstate);
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      for (istate = 0; istate < nstate; istate++) \
        tmp_v02[istate] = real(conj(tmp_v0[j][istate]) * tmp_v0[j][istate]);
      Sz = 0.0;
      for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
        //prefactor 0.5 is added later.
        if (isite1 > Def::Nsite) {
          Sz += GetLocal2Sz(isite1, myrank, Def::SiteToBit, Def::Tpow);
        }
        else {
          Sz += GetLocal2Sz(isite1, j - 1, Def::SiteToBit, Def::Tpow);
        }
      }
      for (istate = 0; istate < nstate; istate++) {
        Sz_t[mythread][istate] += tmp_v02[istate] * Sz;
        Sz2_t[mythread][istate] += tmp_v02[istate] * Sz * Sz;
      }
    }/*for (j = 1; j <= i_max; j++)*/
    free_d_1d_allocate(tmp_v02);
  }/*End of parallel region*/
  for (istate = 0; istate < nstate; istate++) {
    Phys::Sz[istate] = 0.0;
    Phys::Sz2[istate] = 0.0;
    for (mythread = 0; mythread < nthreads; mythread++) {
      Phys::Sz[istate] += Sz_t[mythread][istate];
      Phys::Sz2[istate] += Sz2_t[mythread][istate];
    }
  }
  SumMPI_dv(nstate, Phys::Sz);
  SumMPI_dv(nstate, Phys::Sz2);

  for (istate = 0; istate < nstate; istate++) {
    Phys::doublon[istate] = 0.0;
    Phys::doublon2[istate] = 0.0;
    Phys::num[istate] = Def::NsiteMPI;
    Phys::num2[istate] = Def::NsiteMPI*Def::NsiteMPI;
    Phys::Sz[istate] *= 0.5;
    Phys::Sz2[istate] *= 0.25;
    Phys::num_up[istate] = 0.5*(Phys::num[istate] + Phys::Sz[istate]);
    Phys::num_down[istate] = 0.5*(Phys::num[istate] - Phys::Sz[istate]);
  }

  free_d_2d_allocate(Sz_t);
  free_d_2d_allocate(Sz2_t);
  return 0;
}
///
/// \brief Calculate expected values of energies and physical quantities for Half-Spin model
/// \param X [in, out] X Struct to get information about file header names, dimension of hirbert space, calc type and output physical quantities.
/// \retval 0 normally finished.
/// \retval -1 abnormally finished.
int expec_energy_flct_HalfSpin(
  
  int nstate,
  std::complex<double> **tmp_v0
) {
  long int j;
  long int isite1;
  long int is1_up_a, is1_up_b;

  long int ibit1;
  double Sz;
  double *tmp_v02;
  long int i_max, tmp_list_1;
  int l_ibit1, u_ibit1, i_32;
  int istate, mythread;
  double **Sz_t, **Sz2_t;

  i_max = Check::idim_max;
  Sz_t = d_2d_allocate(nthreads, nstate);
  Sz2_t = d_2d_allocate(nthreads, nstate);
  i_32 = 0xFFFFFFFF; //2^32 - 1

  //[s] for bit count
  is1_up_a = 0;
  is1_up_b = 0;
  for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
    if (isite1 > Def::Nsite) {
      is1_up_a += Def::Tpow[isite1 - 1];
    }
    else {
      is1_up_b += Def::Tpow[isite1 - 1];
    }
  }
  //[e]
#pragma omp parallel default(none) \
private(j,Sz,ibit1,isite1,tmp_v02,u_ibit1,l_ibit1, tmp_list_1,mythread,istate) \
shared(tmp_v0, list_1,Sz_t,Sz2_t,nstate,i_max,myrank,i_32,is1_up_a,is1_up_b,Def::NsiteMPI) 

  {
    tmp_v02 = d_1d_allocate(nstate);
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      for (istate = 0; istate < nstate; istate++) \
        tmp_v02[istate] = real(conj(tmp_v0[j][istate]) * tmp_v0[j][istate]);
      Sz = 0.0;
      tmp_list_1 = list_1[j];

      // isite1 > Def::Nsite
      ibit1 = (long int) myrank & is1_up_a;
      u_ibit1 = ibit1 >> 32;
      l_ibit1 = ibit1 & i_32;
      Sz += pop(u_ibit1);
      Sz += pop(l_ibit1);
      // isite1 <= Def::Nsite
      ibit1 = (long int) tmp_list_1 &is1_up_b;
      u_ibit1 = ibit1 >> 32;
      l_ibit1 = ibit1 & i_32;
      Sz += pop(u_ibit1);
      Sz += pop(l_ibit1);
      Sz = 2 * Sz - Def::NsiteMPI;

      for (istate = 0; istate < nstate; istate++) {
        Sz_t[mythread][istate] += tmp_v02[istate] * Sz;
        Sz2_t[mythread][istate] += tmp_v02[istate] * Sz * Sz;
      }
    }/*for (j = 1; j <= i_max; j++)*/
    free_d_1d_allocate(tmp_v02);
  }/*End of parallel region*/
  for (istate = 0; istate < nstate; istate++) {
    Phys::Sz[istate] = 0.0;
    Phys::Sz2[istate] = 0.0;
    for (mythread = 0; mythread < nthreads; mythread++) {
      Phys::Sz[istate] += Sz_t[mythread][istate];
      Phys::Sz2[istate] += Sz2_t[mythread][istate];
    }
  }
  SumMPI_dv(nstate, Phys::Sz);
  SumMPI_dv(nstate, Phys::Sz2);

  for (istate = 0; istate < nstate; istate++) {
    Phys::doublon[istate] = 0.0;
    Phys::doublon2[istate] = 0.0;
    Phys::num[istate] = Def::NsiteMPI;
    Phys::num2[istate] = Def::NsiteMPI*Def::NsiteMPI;
    Phys::Sz[istate] *= 0.5;
    Phys::Sz2[istate] *= 0.25;
    Phys::num_up[istate] = 0.5*(Phys::num[istate] + Phys::Sz[istate]);
    Phys::num_down[istate] = 0.5*(Phys::num[istate] - Phys::Sz[istate]);
  }

  free_d_2d_allocate(Sz_t);
  free_d_2d_allocate(Sz2_t); 
  return 0;
}
///
/// \brief Calculate expected values of energies and physical quantities for General-Spin model
/// \param X [in, out] X Struct to get information about file header names, dimension of hirbert space, calc type and output physical quantities.
/// \retval 0 normally finished.
/// \retval -1 abnormally finished.
int expec_energy_flct_GeneralSpin(
  
  int nstate,
  std::complex<double> **tmp_v0
) {
  long int j;
  long int isite1;
  int istate, mythread;
  double Sz;
  double *tmp_v02;
  long int i_max, tmp_list1;
  double **Sz_t, **Sz2_t;

  Sz_t = d_2d_allocate(nthreads, nstate);
  Sz2_t = d_2d_allocate(nthreads, nstate);
  i_max = Check::idim_max;

#pragma omp parallel default(none) \
private(j,Sz,isite1,tmp_v02, tmp_list1,istate,mythread) \
shared(tmp_v0, list_1,Sz_t,Sz2_t,nstate,i_max,myrank, \
Def::NsiteMPI, Def::SiteToBit, Def::Tpow,Def::Nsite)
  {
    tmp_v02 = d_1d_allocate(nstate);
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      for (istate = 0; istate < nstate; istate++)
        tmp_v02[istate] = real(conj(tmp_v0[j][istate]) * tmp_v0[j][istate]);
      Sz = 0.0;
      tmp_list1 = list_1[j];
      for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
        //prefactor 0.5 is added later.
        if (isite1 > Def::Nsite) {
          Sz += GetLocal2Sz(isite1, myrank, Def::SiteToBit, Def::Tpow);
        }
        else {
          Sz += GetLocal2Sz(isite1, tmp_list1, Def::SiteToBit, Def::Tpow);
        }
      }
      for (istate = 0; istate < nstate; istate++) {
        Sz_t[mythread][istate] += tmp_v02[istate] * Sz;
        Sz2_t[mythread][istate] += tmp_v02[istate] * Sz * Sz;
      }
    }/*for (j = 1; j <= i_max; j++)*/
    free_d_1d_allocate(tmp_v02);
  }/*End of parallel region*/
  for (istate = 0; istate < nstate; istate++) {
    Phys::Sz[istate] = 0.0;
    Phys::Sz2[istate] = 0.0;
    for (mythread = 0; mythread < nthreads; mythread++) {
      Phys::Sz[istate] += Sz_t[mythread][istate];
      Phys::Sz2[istate] += Sz2_t[mythread][istate];
    }
  }
  SumMPI_dv(nstate, Phys::Sz);
  SumMPI_dv(nstate, Phys::Sz2);

  for (istate = 0; istate < nstate; istate++) {
    Phys::doublon[istate] = 0.0;
    Phys::doublon2[istate] = 0.0;
    Phys::num[istate] = Def::NsiteMPI;
    Phys::num2[istate] = Def::NsiteMPI*Def::NsiteMPI;
    Phys::Sz[istate] *= 0.5;
    Phys::Sz2[istate] *= 0.25;
    Phys::num_up[istate] = 0.5*(Phys::num[istate] + Phys::Sz[istate]);
    Phys::num_down[istate] = 0.5*(Phys::num[istate] - Phys::Sz[istate]);
  }

  free_d_2d_allocate(Sz_t);
  free_d_2d_allocate(Sz2_t);
  return 0;
}
/**
 * @brief Parent function to calculate expected values of energy and physical quantities.
 *
 * @param X [in,out] X Struct to get information about file header names, dimension of hirbert space, calc type, physical quantities.
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * \retval 0 normally finished.
 * \retval -1 abnormally finished.
 */
int expec_energy_flct(
  
  int nstate,
  std::complex<double> **tmp_v0,
  std::complex<double> **tmp_v1
) {

  long int i, j;
  long int irght, ilft, ihfbit;
  long int i_max;
  int istate;

  switch (Def::iCalcType) {
  case TPQCalc:
  case TimeEvolution:
#ifdef _DEBUG
    fprintf(stdoutMPI, "%s", "  Start: Calculate Energy.\n");
    TimeKeeperWithStep("%s_TimeKeeper.dat", "step %d: Calculate energy begins:      %s", "a", step_i);
#endif
    break;
  case FullDiag:
  case CG:
    break;
  default:
    return -1;
  }

  i_max = Check::idim_max;
  if (GetSplitBitByModel(Def::Nsite, Def::iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }

  Large::i_max = i_max;
  Large::irght = irght;
  Large::ilft = ilft;
  Large::ihfbit = ihfbit;
  Large::mode = M_ENERGY;
  for (istate = 0; istate < nstate; istate++) Phys::energy[istate] = 0.0;

  int nCalcFlct;
  if (Def::iCalcType == TPQCalc) {
    nCalcFlct = 3201;
  }
  else {//For FullDiag
    nCalcFlct = 5301;
  }
  StartTimer(nCalcFlct);

  switch (Def::iCalcModel) {
  case HubbardGC:
    expec_energy_flct_HubbardGC(nstate, tmp_v0);
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    expec_energy_flct_Hubbard(nstate, tmp_v0);
    break;

  case SpinGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      expec_energy_flct_HalfSpinGC(nstate, tmp_v0);
    }
    else {//for generalspin
      expec_energy_flct_GeneralSpinGC(nstate, tmp_v0);
    }
    break;/*case SpinGC*/
    /* SpinGCBoost */
  case Spin:
    /*
    if(Def::iFlgGeneralSpin == FALSE){
      expec_energy_flct_HalfSpin(X);
    }
    else{
      expec_energy_flct_GeneralSpin(X);
    }
     */
    for (istate = 0; istate < nstate; istate++) {
      Phys::doublon[istate] = 0.0;
      Phys::doublon2[istate] = 0.0;
      Phys::num[istate] = Def::NsiteMPI;
      Phys::num2[istate] = Def::NsiteMPI*Def::NsiteMPI;
      Phys::Sz[istate] = 0.5 * (double)Def::Total2SzMPI;
      Phys::Sz2[istate] = Phys::Sz[istate] * Phys::Sz[istate];
    }
    break;
  default:
    return -1;
  }

  StopTimer(nCalcFlct);

#pragma omp parallel for default(none) private(i,istate) \
shared(tmp_v1,tmp_v0,nstate,i_max)
  for (i = 1; i <= i_max; i++) {
    for (istate = 0; istate < nstate; istate++) {
      tmp_v1[i][istate] = tmp_v0[i][istate];
      tmp_v0[i][istate] = 0.0;
    }
  }

  int nCalcExpec;
  if (Def::iCalcType == TPQCalc) {
    nCalcExpec = 3202;
  }
  else {//For FullDiag
    nCalcExpec = 5302;
  }
  StartTimer(nCalcExpec);
  mltply(nstate, tmp_v0, tmp_v1); // v0+=H*v1
  StopTimer(nCalcExpec);
  /* switch -> SpinGCBoost */

  for (istate = 0; istate < nstate; istate++) {
    Phys::energy[istate] = 0.0;
    Phys::var[istate] = 0.0;
  }
  for (j = 1; j <= i_max; j++) {
    for (istate = 0; istate < nstate; istate++) {
      Phys::energy[istate] += real(conj(tmp_v1[j][istate])*tmp_v0[j][istate]); // E   = <v1|H|v1>=<v1|v0>
      Phys::var[istate] += real(conj(tmp_v0[j][istate])*tmp_v0[j][istate]); // E^2 = <v1|H*H|v1>=<v0|v0>
    }
  }
  SumMPI_dv(nstate, Phys::energy);
  SumMPI_dv(nstate, Phys::var);

  switch (Def::iCalcType) {
  case TPQCalc:
  case TimeEvolution:
#ifdef _DEBUG
    fprintf(stdoutMPI, "%s", "  End  : Calculate Energy.\n");
    TimeKeeperWithStep("%s_TimeKeeper.dat", "step %d: Calculate energy finishes:    %s", "a", step_i);
#endif
    break;
  default:
    break;
  }
  return 0;
}
