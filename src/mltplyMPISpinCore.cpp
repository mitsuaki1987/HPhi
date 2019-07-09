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
@brief Functions for spin Hamiltonian + MPI (Core)

General two body term:
<table>
  <tr>
    <td></td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td></td>
    <td>MPI single</td>
    <td>MPI double</td>
    <td>MPI single</td>
    <td>MPI double</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{ju}@f$</td>
    <td>::GC_child_CisAisCjuAju_spin_MPIsingle, ::X_GC_child_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_child_CisAisCjuAju_spin_MPIdouble, ::X_GC_child_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::X_child_CisAisCjuAju_GeneralSpin_MPIsingle, ::X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
    <td>::X_child_CisAisCjuAju_GeneralSpin_MPIsingle, ::X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{jv}@f$</td>
    <td>::GC_child_CisAisCjuAjv_spin_MPIsingle, ::X_GC_child_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_child_CisAisCjuAjv_spin_MPIdouble, ::X_GC_child_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::X_child_CisAisCjuAjv_GeneralSpin_MPIsingle, ::X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
    <td>::X_child_CisAisCjuAjv_GeneralSpin_MPIsingle, ::X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{ju}@f$</td>
    <td>::GC_child_CisAitCjuAju_spin_MPIsingle, ::X_GC_child_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_child_CisAitCjuAju_spin_MPIdouble, ::X_GC_child_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::X_child_CisAitCjuAju_GeneralSpin_MPIsingle, ::X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle</td>
    <td>::X_child_CisAitCjuAju_GeneralSpin_MPIsingle, ::X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{jv}@f$</td>
    <td>::GC_child_CisAitCjuAjv_spin_MPIsingle, ::X_GC_child_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_child_CisAitCjuAjv_spin_MPIdouble, ::X_GC_child_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::X_child_CisAitCjuAjv_GeneralSpin_MPIsingle, ::X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle</td>
    <td>::X_child_CisAitCjuAjv_GeneralSpin_MPIsingle, ::X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle</td>
  </tr>
</table>
*/
#include "mltplyCommon.hpp"
#include "mltplySpinCore.hpp"
#include "mltplyMPISpinCore.hpp"
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "global.hpp"
/**
@brief Exchange and Pairlifting term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_child_CisAitCiuAiv_spin_MPIdouble(
  long int i_int /**< [in] Interaction ID*/,
  int nstate, std::complex<double> **tmp_v0 /**< [out] Result v0 = H v1*/,
  std::complex<double> **tmp_v1 /**< [in] v0 = H v1*/)
{
   X_GC_child_CisAitCiuAiv_spin_MPIdouble(
    Def::InterAll_OffDiagonal[i_int][0],  Def::InterAll_OffDiagonal[i_int][1], 
    Def::InterAll_OffDiagonal[i_int][3],  Def::InterAll_OffDiagonal[i_int][4], 
    Def::InterAll_OffDiagonal[i_int][5],  Def::InterAll_OffDiagonal[i_int][7],
    Def::ParaInterAll_OffDiagonal[i_int],nstate, tmp_v0, tmp_v1);
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief @f$c_{is}^\dagger c_{it} c_{iu}^\dagger c_{iv}@f$ term
       in Spin model + GC.
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void X_GC_child_CisAitCiuAiv_spin_MPIdouble(
  int org_isite1,//!<[in] site i
  int org_ispin1,//!<[in] spin s
  int org_ispin2,//!<[in] spin t
  int org_isite3,//!<[in] site i?
  int org_ispin3,//!<[in] spin u
  int org_ispin4,//!<[in] spin v
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  int mask1, mask2, state1, state2, origin;
  long int idim_max_buf;
  std::complex<double> Jint;

  mask1 = (int)Def::Tpow[org_isite1];
  mask2 = (int)Def::Tpow[org_isite3];
  if (org_isite1 != org_isite3) {
    origin = MP::myrank ^ (mask1 + mask2);
  }
  else {
    if (org_ispin1 == org_ispin4 && org_ispin2 == org_ispin3) { //CisAitCitAis=CisAis
      X_GC_child_CisAis_spin_MPIdouble(org_isite1, org_ispin1, tmp_J, nstate, tmp_v0, tmp_v1);
      return;
    }
    else { //CisAitCisAit=0
      return;
    }
  }

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == org_ispin2 && state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (state1 == org_ispin1 && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else {
    return;
  }

  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(nstate * idim_max_buf, Jint, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief Wrapper for calculating CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAisCjuAjv_spin_MPIdouble(
  long int i_int /**< [in] Interaction ID*/,
  int nstate, std::complex<double> **tmp_v0 /**< [out] Result v0 = H v1*/,
  std::complex<double> **tmp_v1 /**< [in] v0 = H v1*/
){
  X_GC_child_CisAisCjuAjv_spin_MPIdouble(
    Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
    Def::InterAll_OffDiagonal[i_int][4], Def::InterAll_OffDiagonal[i_int][5],
    Def::InterAll_OffDiagonal[i_int][7], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAisCjuAjv_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  int mask1, mask2, state2;
  long int origin, num1;
  long int idim_max_buf;
  std::complex<double> Jint;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4) {//CisAisCitAis
      return;
  }

  mask1 = (int)Def::Tpow[org_isite1];
  mask2 = (int)Def::Tpow[org_isite3];
  origin = MP::myrank ^ mask2;
  state2 = (origin & mask2) / mask2;
  num1 = X_SpinGC_CisAis((long int) MP::myrank + 1, mask1, org_ispin1);
  if (num1 != 0 && state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (X_SpinGC_CisAis(origin + 1, mask1, org_ispin1) == TRUE && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) Jint = 0;
  }
  else {
    return;
  }

  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(Check::idim_max*nstate, Jint, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAisCjuAjv_spin_MPIdouble*/
/**
@brief Wrapper for calculating CisAitCjuAju term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAitCjuAju_spin_MPIdouble(
  long int i_int,//!<[in] Interaction ID
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
)
{
  X_GC_child_CisAitCjuAju_spin_MPIdouble(
    Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
    Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4], 
    Def::InterAll_OffDiagonal[i_int][5], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAitCjuAju_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  int mask1, mask2, state1, num1;
  long int origin;
  long int idim_max_buf;
  std::complex<double> Jint;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {//cisaitcisais
    return;
  }

  mask1 = (int)Def::Tpow[org_isite1];
  origin = MP::myrank ^ mask1;
  state1 = (origin & mask1) / mask1;
  mask2 = (int)Def::Tpow[org_isite3];
  num1 = X_SpinGC_CisAis(origin + 1, mask2, org_ispin3);
  if (state1 == org_ispin2) {
    if (num1 != 0) {
      Jint = tmp_J;
    }
    else {
      return;
    }
  }/*if (state1 == org_ispin2)*/
  else {//state1 = org_ispin1
    num1 = X_SpinGC_CisAis((long int) MP::myrank + 1, mask2, org_ispin3);
    if (num1 != 0) {
      Jint = conj(tmp_J);
      if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
        Jint = 0;
      }
    }
    else {
      return;
    }
  }

  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(Check::idim_max*nstate, Jint, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAisCjuAjv_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAisCjuAju_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
){
  long int mask1, mask2, num1,num2;
  long int  j;
  std::complex<double> dmv;
  int one = 1;
  mask1 = (int)Def::Tpow[org_isite1];
  mask2 = (int)Def::Tpow[org_isite3];
  num1 = X_SpinGC_CisAis((long int)MP::myrank + 1, mask1, org_ispin1);
  num2 = X_SpinGC_CisAis((long int)MP::myrank + 1, mask2, org_ispin3);
  
#pragma omp parallel for default(none) private(j, dmv) \
shared(tmp_J, num1, num2,tmp_v1, tmp_v0,nstate,one,Check::idim_max)
  for (j = 1; j <= Check::idim_max; j++) {
    dmv = (std::complex<double>)(num1 * num2) * tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }/*for (j = 1; j <= Check::idim_max; j++) */
}/*std::complex<double> X_GC_child_CisAisCjuAju_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAisCjuAju_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  long int mask1, mask2, num1, num2;
  long int j;
  std::complex<double> Jint, dmv;
  int one = 1;
  Jint = tmp_J;
  mask1 = (int)Def::Tpow[org_isite1];
  mask2 = (int)Def::Tpow[org_isite3];
  num2 = X_SpinGC_CisAis((long int) MP::myrank + 1, mask2, org_ispin3);

#pragma omp parallel for default(none) private(j, dmv, num1) \
shared(Jint, num2, mask1, org_ispin1, tmp_v1, tmp_v0,nstate,one,Check::idim_max)
  for (j = 1; j <= Check::idim_max; j++) {
    num1 = X_SpinGC_CisAis(j, mask1, org_ispin1);
    dmv = Jint * (std::complex<double>)(num1 * num2);
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }/*for (j = 1; j <= Check::idim_max; j++)*/
}/*std::complex<double> X_GC_child_CisAisCjuAju_spin_MPIdouble*/
/**
@brief Exchange and Pairlifting term in Spin model + GC
       When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_child_CisAitCiuAiv_spin_MPIsingle(
  long int i_int,//!<[in] Interaction ID
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  X_GC_child_CisAitCiuAiv_spin_MPIsingle(
    Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
    Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
    Def::InterAll_OffDiagonal[i_int][5], Def::InterAll_OffDiagonal[i_int][7],
    Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/
/**
@brief Exchange and Pairlifting term in Spin model + GC
       When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void X_GC_child_CisAitCiuAiv_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  int mask2, state2, origin;
  long int mask1, idim_max_buf, j, ioff, state1, state1check;
  std::complex<double> Jint;
  int one = 1;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)Def::Tpow[org_isite3];
  origin = MP::myrank ^ mask2;
  state2 = (origin & mask2) / mask2;

  if (state2 == org_ispin4) {
    state1check = (long int) org_ispin2;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (long int) org_ispin1;
    Jint = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return;

  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);
  /*
  Index in the intra PE
  */
  mask1 = Def::Tpow[org_isite1];

#pragma omp parallel for default(none) private(j, state1, ioff) \
shared(idim_max_buf, Jint, state1check, mask1,Wave::v1buf, tmp_v1, tmp_v0,nstate,one)
  for (j = 0; j < idim_max_buf; j++) {
    state1 = X_SpinGC_CisAit(j + 1, mask1, state1check, &ioff);
    if (state1 != 0) {
      zaxpy_(&nstate, &Jint, &Wave::v1buf[j + 1][0], &one, &tmp_v0[ioff + 1][0], &one);
    }/*if (state1 != 0)*/
  }/*for (j = 0; j < idim_max_buf; j++)*/
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/
/**
@brief Wrapper for CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAisCjuAjv_spin_MPIsingle(
  long int i_int,//!<[in] Interaction ID
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  X_GC_child_CisAisCjuAjv_spin_MPIsingle(
    Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
    Def::InterAll_OffDiagonal[i_int][4], Def::InterAll_OffDiagonal[i_int][5],
    Def::InterAll_OffDiagonal[i_int][7], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
}/*void GC_child_CisAisCjuAjv_spin_MPIsingle*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAisCjuAjv_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 1
  int org_ispin3,//!<[in] Spin 2
  int org_ispin4,//!<[in] Spin 2
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  int mask2, state2, origin;
  long int mask1, idim_max_buf, j, state1, state1check;
  std::complex<double> Jint;
  int one = 1;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)Def::Tpow[org_isite3];
  origin = MP::myrank ^ mask2;
  state2 = (origin & mask2) / mask2;
  if (state2 == org_ispin4) {
    state1check = (long int) org_ispin1;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (long int) org_ispin1;
    Jint = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return;

  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);
  /*
  Index in the intra PE
  */
  mask1 = Def::Tpow[org_isite1];

#pragma omp parallel for default(none) private(j, state1) \
shared(idim_max_buf, Jint, state1check, mask1, Wave::v1buf, tmp_v1, tmp_v0,nstate,one)
  for (j = 0; j < idim_max_buf; j++) {
    state1 = (j & mask1) / mask1;
    if (state1 == state1check) {
      zaxpy_(&nstate, &Jint, &Wave::v1buf[j + 1][0], &one, &tmp_v0[j + 1][0], &one);
    }/*if (state1 == state1check)*/
  }/*for (j = 0; j < idim_max_buf; j++)*/
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/
/**
@brief Wrapper for CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAitCjuAju_spin_MPIsingle(
  long int i_int,//!<[in] Interaction ID
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  X_GC_child_CisAitCjuAju_spin_MPIsingle(
    Def::InterAll_OffDiagonal[i_int][0], 
    Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
    Def::InterAll_OffDiagonal[i_int][5], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
}/*void GC_child_CisAisCjuAjv_spin_MPIsingle*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAitCjuAju_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  int mask2, state2;
  long int mask1, j, ioff, state1, state1check;
  std::complex<double> Jint, dmv;
  int one = 1;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)Def::Tpow[org_isite3];
  state2 = (MP::myrank & mask2) / mask2;

  if (state2 == org_ispin3) {
    state1check = org_ispin2;
    Jint = tmp_J;
  }
  else {
    return;
  }

  mask1 = (int)Def::Tpow[org_isite1];

#pragma omp parallel for default(none) private(j, dmv, state1, ioff) \
shared(Jint, state1check, mask1,tmp_v1, tmp_v0,nstate,one,Check::idim_max)
  for (j = 0; j < Check::idim_max; j++) {

    state1 = (j & mask1) / mask1;
    ioff = j ^ mask1;
    if (state1 == state1check) {
      dmv = Jint;
    }
    else {
      dmv = conj(Jint);
    }
    zaxpy_(&nstate, &dmv, &tmp_v1[j + 1][0], &one, &tmp_v0[ioff + 1][0], &one);
  }/*for (j = 0; j < Check::idim_max; j++)*/
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/
/**
@brief @f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{jv}@f$ term in Spin model.
 When both site1 and site3 are in the inter process region.
*/
void X_GC_child_CisAisCjuAjv_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  long int off;
  int origin;
  std::complex<double> tmp_V;
  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4) {//cisaisciuais=0 && cisaiucisais=0
    return;
  }

  if (BitCheckGeneral(MP::myrank, org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow) == TRUE
    && GetOffCompGeneralSpin((long int) MP::myrank, org_isite3 + 1, org_ispin3, org_ispin4,
      &off, Def::SiteToBit, Def::Tpow) == TRUE)
    tmp_V = tmp_J;
  else {
    if (GetOffCompGeneralSpin((long int) MP::myrank, org_isite3 + 1, org_ispin4, org_ispin3,
      &off, Def::SiteToBit, Def::Tpow) == TRUE)
    {
      if (BitCheckGeneral((long int)off, org_isite1 + 1, org_ispin1, Def::SiteToBit,
        Def::Tpow) == TRUE)
      {
        tmp_V = conj(tmp_J);
        if(Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0.0;
      }/*BitCheckGeneral(off, org_isite1 + 1, org_ispin1)*/
      else return;
    }/*GetOffCompGeneralSpin(MP::myrank, org_isite3 + 1, org_ispin4, org_ispin3, &off)*/
    else return;
  }
  origin = (int)off;
  SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(Check::idim_max*nstate, tmp_V, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAisCjuAjv_GeneralSpin_MPIdouble*/
/**
@brief @f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{ju}@f$ term in Spin model.
When both site1 and site3 are in the inter process region.
*/
void X_GC_child_CisAitCjuAju_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1//!<[in] Vector to be producted
) {
  long int off;
  int origin;
  std::complex<double> tmp_V;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {//cisaitcisais=0 && cisaiscitais=0
    return;
  }

  if (BitCheckGeneral(MP::myrank, org_isite3 + 1, org_ispin3, Def::SiteToBit, Def::Tpow) == TRUE
    && GetOffCompGeneralSpin((long int) MP::myrank, org_isite1 + 1, org_ispin2, org_ispin1, &off,
      Def::SiteToBit, Def::Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0.0;
  }
  else if (GetOffCompGeneralSpin((long int) MP::myrank, org_isite1 + 1, org_ispin1, org_ispin2,
    &off, Def::SiteToBit, Def::Tpow) == TRUE)
  {
    if (BitCheckGeneral((long int)off, org_isite3 + 1, org_ispin3,
      Def::SiteToBit, Def::Tpow) == TRUE) {
      tmp_V = tmp_J;
    }
    else return;
  }
  else return;

  origin = (int)off;

  SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(Check::idim_max*nstate, tmp_V, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAitCjuAju_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{jv}@f$ term in the
grandcanonical general spin system when both site is in the inter process region
*/
void X_GC_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
) {
  long int tmp_off, off;
  int origin, ihermite;
  std::complex<double> tmp_V;

  ihermite = TRUE;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4 &&
    org_ispin2 == org_ispin3) { //cisaitcitais=cisais && cisaitcitais =cisais
    X_GC_child_CisAis_GeneralSpin_MPIdouble(org_isite1, org_ispin1, tmp_J, nstate, tmp_v0, tmp_v1);
    return;
  }
  //cisaitcisait
  if (GetOffCompGeneralSpin((long int) MP::myrank, org_isite1 + 1, org_ispin1, org_ispin2,
    &tmp_off, Def::SiteToBit, Def::Tpow) == TRUE) {

    if (GetOffCompGeneralSpin(tmp_off, org_isite3 + 1, org_ispin3, org_ispin4,
      &off, Def::SiteToBit, Def::Tpow) == TRUE) {

      tmp_V = tmp_J;
    }
    else ihermite = FALSE;
  }
  else {
    ihermite = FALSE;
  }

  if (ihermite == FALSE) {
    if (GetOffCompGeneralSpin((long int) MP::myrank, org_isite3 + 1, org_ispin4, org_ispin3, &tmp_off,
      Def::SiteToBit, Def::Tpow) == TRUE) {

      if (GetOffCompGeneralSpin(tmp_off, org_isite1 + 1, org_ispin2, org_ispin1, &off, Def::SiteToBit,
                                      Def::Tpow) == TRUE) {
        tmp_V = conj(tmp_J);
        if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0.0;
      }
      else return;
    }
    else return;
  }

  origin = (int)off;

  SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(Check::idim_max*nstate, tmp_V, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAitCjuAjv_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{ju}@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 */
void X_GC_child_CisAisCjuAju_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
) {
  long int num1;
  std::complex<double> tmp_V;

  num1 = BitCheckGeneral((long int) MP::myrank, org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow);

  if (num1 == TRUE) {
    num1 = BitCheckGeneral((long int) MP::myrank, org_isite3 + 1, org_ispin3, Def::SiteToBit, Def::Tpow);
    if (num1 == TRUE) {
      tmp_V = tmp_J;
    }
    else return;
  }
  else return;

  zaxpy_long(Check::idim_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAisCjuAju_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{it}@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 */
void X_GC_child_CisAit_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
) {
  long int off;
  int origin;
  std::complex<double> tmp_V;

  if (GetOffCompGeneralSpin((long int) MP::myrank, org_isite1 + 1, org_ispin1, org_ispin2,
    &off, Def::SiteToBit, Def::Tpow) == TRUE) {
    tmp_V = tmp_trans;
  }
  else if (GetOffCompGeneralSpin((long int) MP::myrank,
    org_isite1 + 1, org_ispin2, org_ispin1, &off,
    Def::SiteToBit, Def::Tpow) == TRUE) {
    tmp_V = conj(tmp_trans);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0.0;
  }
  else return;

  origin = (int)off;

  SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(Check::idim_max*nstate, tmp_V, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAit_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{is}@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 */
void X_GC_child_CisAis_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
) {
  long int num1;
  std::complex<double> tmp_V;

  num1 = BitCheckGeneral((long int) MP::myrank,
    org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow);
  if (num1 != 0) {
    tmp_V = tmp_trans;
  }
  else return;

  zaxpy_long(Check::idim_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAis_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is} c_{is}^\dagger@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 */
void X_GC_child_AisCis_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
) {
  long int num1;
  std::complex<double> tmp_V;

  num1 = BitCheckGeneral((long int) MP::myrank,
    org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow);
  if (num1 == 0) {
    tmp_V = tmp_trans;
  }
  else return;

  zaxpy_long(Check::idim_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_AisCis_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{it}@f$ term in the
canonical general spin system when both site is in the inter process region
*/
void X_child_CisAit_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  //!<[inout]
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Input wavefunction
  long int idim_max//!<[in] Similar to CheckList::idim_max
)
{
  long int off, j, tmp_off,idim_max_buf;
  int origin, one = 1;
  std::complex<double> tmp_V;
  
  if (GetOffCompGeneralSpin((long int) MP::myrank, org_isite1 + 1, org_ispin1, org_ispin2,
                            &off, Def::SiteToBit, Def::Tpow) == TRUE) {
    tmp_V = tmp_trans;
  }
  else if (GetOffCompGeneralSpin((long int) MP::myrank,
                                 org_isite1 + 1, org_ispin2, org_ispin1, &off,
                                 Def::SiteToBit, Def::Tpow) == TRUE) {
    tmp_V = conj(tmp_trans);
    if (Large::mode == M_CORR || Large::mode ==M_CALCSPEC) tmp_V = 0.0;
  }
  else return;
  
  origin = (int) off;

  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_iv(origin, idim_max + 1, idim_max_buf + 1, List::c1_org, List::c1buf_org);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

#pragma omp parallel for default(none) private(j, tmp_off) \
shared(tmp_V, idim_max_buf, List::c1buf_org,tmp_v0, tmp_v1, Wave::v1buf,nstate,one, Large::ihfbit)
  for (j = 1; j <= idim_max_buf; j++) {
    ConvertToList1GeneralSpin(List::c1buf_org[j], Large::ihfbit, &tmp_off);
    zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[tmp_off][0], &one);
  }/*for (j = 1; j <= idim_max_buf; j++)*/
}/*std::complex<double> X_child_CisAit_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{jv}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
*/
void X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
){
  long int off, j, num1;
  int origin, isite, IniSpin, one = 1;
  std::complex<double> tmp_V;

  if (GetOffCompGeneralSpin((long int)MP::myrank,
    org_isite3 + 1, org_ispin3, org_ispin4, &off,
    Def::SiteToBit, Def::Tpow) == TRUE)
  {
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
  }
  else if (GetOffCompGeneralSpin((long int)MP::myrank,
    org_isite3 + 1, org_ispin4, org_ispin3, &off,
    Def::SiteToBit, Def::Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0.0;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
  }
  else return;
  
  origin = (int)off;
  
  SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

#pragma omp parallel default(none) private(j, num1) \
shared(tmp_v0, tmp_v1, Wave::v1buf,nstate,one,tmp_V, isite, IniSpin, Def::SiteToBit, Def::Tpow,Check::idim_max)
  {
#pragma omp for
    for (j = 1; j <= Check::idim_max; j++) {
      num1 = BitCheckGeneral(j - 1, isite, IniSpin, Def::SiteToBit, Def::Tpow);
      if (num1 != 0) zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[j][0], &one);
    }/*for (j = 1; j <= Check::idim_max; j++)*/
  }/*End of parallel region*/
}/*std::complex<double> X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{it}c_{ju}^\dagger c_{ju}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
*/
void X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
){
  long int num1, j, off;
  int isite, IniSpin, FinSpin, one = 1;
  std::complex<double> tmp_V, dmv;

  num1 = BitCheckGeneral((long int)MP::myrank, 
    org_isite3+1, org_ispin3, Def::SiteToBit, Def::Tpow);
  if(num1 != 0){
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else return;

#pragma omp parallel for default(none) private(j, dmv, num1, off) \
shared(tmp_V, isite, IniSpin, FinSpin,tmp_v0, tmp_v1, Wave::v1buf,nstate,one, \
Check::idim_max,Def::SiteToBit, Def::Tpow)
  for (j = 1; j <= Check::idim_max; j++) {
    if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
      Def::SiteToBit, Def::Tpow) == TRUE)
    {
      dmv = tmp_V;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[off + 1][0], &one);
    }
    else if (GetOffCompGeneralSpin(j - 1, isite, FinSpin, IniSpin, &off,
      Def::SiteToBit, Def::Tpow) == TRUE)
    {
      dmv = conj(tmp_V);
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[off + 1][0], &one);
    }
  }/*for (j = 1; j <= Check::idim_max; j++)*/
}/*std::complex<double> X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{jv}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
*/
void X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
){
  long int off, j;
  int origin, isite, IniSpin, FinSpin, one = 1;
  std::complex<double> tmp_V;

  if (GetOffCompGeneralSpin((long int)MP::myrank,
    org_isite3 + 1, org_ispin3, org_ispin4, &off,
    Def::SiteToBit, Def::Tpow) == TRUE)
  {
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else if (GetOffCompGeneralSpin((long int)MP::myrank,
    org_isite3 + 1, org_ispin4, org_ispin3, &off,
    Def::SiteToBit, Def::Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0.0;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
    FinSpin = org_ispin2;
  }
  else return;

  origin = (int)off;

  SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

#pragma omp parallel for default(none) private(j, off) \
shared(tmp_V, isite, IniSpin, FinSpin, tmp_v0, tmp_v1, Wave::v1buf,nstate,one, \
Check::idim_max,Def::SiteToBit, Def::Tpow)
  for (j = 1; j <= Check::idim_max; j++) {
    if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
      Def::SiteToBit, Def::Tpow) == TRUE)
    {
      zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[off + 1][0], &one);
    }
  }/*for (j = 1; j <= Check::idim_max; j++)*/
}/*std::complex<double> X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{ju}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
*/
void X_GC_child_CisAisCjuAju_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
){
  long int j, num1;
  std::complex<double> tmp_V, dmv;
  int one = 1;

  num1 = BitCheckGeneral((long int)MP::myrank, org_isite3+1, org_ispin3, Def::SiteToBit, Def::Tpow);
  if (num1 != FALSE) {
    tmp_V = tmp_J;
  }
  else return;
  
#pragma omp parallel for default(none) private(j, dmv, num1) \
shared(tmp_V, org_isite1, org_ispin1,tmp_v0, tmp_v1,nstate,one, \
Check::idim_max, Def::SiteToBit, Def::Tpow)
  for (j = 1; j <= Check::idim_max; j++) {
    num1 = BitCheckGeneral(j - 1, org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow);

    dmv = tmp_V * (std::complex<double>)num1;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }/*for (j = 1; j <= Check::idim_max; j++)*/
}/*std::complex<double> X_GC_child_CisAisCjuAju_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{it}c_{ju}^\dagger c_{jv}@f$ term in the
canonical general spin system when both sites are in the inter process region
*/
void X_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
){
  long int tmp_off, off, j, idim_max_buf;
  int origin, one = 1;
  std::complex<double> tmp_V;
  int ihermite=TRUE;

  if (GetOffCompGeneralSpin((long int)MP::myrank, org_isite1 + 1, org_ispin1, org_ispin2, &tmp_off, Def::SiteToBit, Def::Tpow) == TRUE)
  {
    if (GetOffCompGeneralSpin(tmp_off, org_isite3 + 1, org_ispin3, org_ispin4, &off, Def::SiteToBit, Def::Tpow) == TRUE)
    {
      tmp_V = tmp_J;
    }
    else{
      ihermite =FALSE;
    }
  }
  else{
    ihermite=FALSE;
  }
  
  if (ihermite == FALSE) {
    if (GetOffCompGeneralSpin((long int)MP::myrank, org_isite3 + 1, org_ispin4, org_ispin3, &tmp_off, Def::SiteToBit, Def::Tpow) == TRUE)
    {
      if (GetOffCompGeneralSpin(tmp_off, org_isite1 + 1, org_ispin2, org_ispin1, &off, Def::SiteToBit, Def::Tpow) == TRUE)
      {
        tmp_V = conj(tmp_J);
        if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
          tmp_V = 0.0;
        }
      }
      else return;
    }
    else return;
  }
  
  origin = (int)off;

  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_iv(origin, Check::idim_max + 1, idim_max_buf + 1, List::c1, List::c1buf);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

#pragma omp parallel for default(none) private(j, off) \
shared(tmp_v0, tmp_v1, List::c1buf, Wave::v1buf,nstate,one,tmp_V, idim_max_buf, Check::sdim)
  for (j = 1; j <= idim_max_buf; j++) {
    ConvertToList1GeneralSpin(List::c1buf[j], Check::sdim, &off);
    zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[off][0], &one);
  }/*for (j = 1; j <= idim_max_buf; j++)*/
}/*std::complex<double> X_child_CisAitCjuAjv_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{ju}@f$ term in the
canonical general spin system when both sites are in the inter process region
*/
void X_child_CisAisCjuAju_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
) {
  long int num1;
  std::complex<double> tmp_V;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {
    num1 = BitCheckGeneral((long int) MP::myrank, org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow);
    if (num1 != FALSE) {
      tmp_V = tmp_J;
    }
    else {
      return;
    }
  }
  else {
    num1 = BitCheckGeneral((long int) MP::myrank, org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow);
    if (num1 != FALSE) {
      num1 = BitCheckGeneral((long int) MP::myrank, org_isite3 + 1, org_ispin3, Def::SiteToBit,
        Def::Tpow);
      if (num1 != FALSE) {
        tmp_V = tmp_J;
      }
      else {
        return;
      }
    }
    else {
      return;
    }
  }

  zaxpy_long(Check::idim_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_child_CisAisCjuAju_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{ju}@f$ term in the
canonical general spin system when one of these sites is in the inter process region
*/
void X_child_CisAisCjuAju_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
)
{
  long int j, num1;
  std::complex<double> tmp_V, dmv;
  int one = 1;

  num1 = BitCheckGeneral((long int) MP::myrank, org_isite3 + 1, org_ispin3, Def::SiteToBit, Def::Tpow);
  if (num1 != FALSE) {
    tmp_V = tmp_J;
  }
  else return;

#pragma omp parallel for default(none) private(j, dmv, num1) \
shared(tmp_V, org_isite1, org_ispin1,tmp_v0, tmp_v1, List::c1,nstate, \
one,Check::idim_max, Def::SiteToBit, Def::Tpow)
  for (j = 1; j <= Check::idim_max; j++) {
    num1 = BitCheckGeneral(List::c1[j], org_isite1 + 1, org_ispin1, Def::SiteToBit, Def::Tpow);

    dmv = tmp_V * (std::complex<double>)num1;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }/*for (j = 1; j <= Check::idim_max; j++)*/
}/*std::complex<double> X_child_CisAisCjuAju_GeneralSpin_MPIsingle*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{it}c_{ju}^\dagger c_{jv}@f$ term in the
 canonical general spin system when one of these sites is in the inter process region
 */
void X_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Input wavefunction
){
  long int tmp_off, off, j, idim_max_buf;
  int origin, isite, IniSpin, FinSpin, one = 1;
  std::complex<double> tmp_V;
  
  if (GetOffCompGeneralSpin((long int)MP::myrank,
    org_isite3 + 1, org_ispin3, org_ispin4, &off,
    Def::SiteToBit, Def::Tpow) == TRUE)
  {
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else if (GetOffCompGeneralSpin((long int)MP::myrank,
    org_isite3 + 1, org_ispin4, org_ispin3, &off, Def::SiteToBit, Def::Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0.0;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
    FinSpin = org_ispin2;
  }
  else return;

  origin = (int)off;
  
  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_iv(origin, Check::idim_max + 1, idim_max_buf + 1, List::c1, List::c1buf);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

#pragma omp parallel for default(none) private(j, off, tmp_off) \
shared(tmp_V, idim_max_buf, IniSpin, FinSpin, isite, tmp_v0, tmp_v1, List::c1buf, \
Wave::v1buf,nstate,one,Def::SiteToBit, Def::Tpow, Check::sdim)
  for (j = 1; j <= idim_max_buf; j++) {
    if (GetOffCompGeneralSpin(List::c1buf[j], isite, IniSpin, FinSpin, &tmp_off,
      Def::SiteToBit, Def::Tpow) == TRUE)
    {
      ConvertToList1GeneralSpin(tmp_off, Check::sdim, &off);
      zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[off][0], &one);
    }
  }/*for (j = 1; j <= idim_max_buf; j++)*/
}/*std::complex<double> X_child_CisAitCjuAjv_GeneralSpin_MPIsingle*/
/**
@brief Hopping term in Spin + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAit_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0 /**< [out] Result v0 = H v1*/,
  std::complex<double> **tmp_v1 /**< [in] v0 = H v1*/)
{
  int mask1, state1, origin;
  long int idim_max_buf;
  std::complex<double> trans;
  
  mask1 = (int)Def::Tpow[org_isite1];
  origin = MP::myrank ^ mask1;
  state1 = (origin & mask1)/mask1;

  if(state1 ==  org_ispin2){
    trans = tmp_trans;
  }
  else if(state1 == org_ispin1) {
    trans = conj(tmp_trans);
    if(Large::mode == M_CORR|| Large::mode ==M_CALCSPEC){
      trans = 0.0;
    }
  }
  else{
    return;
  }

  idim_max_buf = SendRecv_i(origin, Check::idim_max);
  SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);

  zaxpy_long(Check::idim_max*nstate, trans, &Wave::v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double>  X_GC_child_CisAit_spin_MPIdouble*/
/**
@brief Hopping term in Spin + Canonical for CalcSpectrum
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_child_CisAit_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0 /**< [out] Result v0 = H v1*/,
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  long int idim_max//!<[in] Similar to CheckList::idim_max
){
  int mask1, state1, origin, one = 1;
  long int idim_max_buf, j;
  long int tmp_off;
  std::complex<double> trans;
  
  mask1 = (int)Def::Tpow[org_isite1];
  origin = MP::myrank ^ mask1;
  state1 = (origin & mask1)/mask1;

  if(state1 ==  org_ispin2){
    trans = tmp_trans;
  }
  else{
    trans =0.0;
  }

  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_iv(origin, idim_max + 1, idim_max_buf + 1, List::c1_org, List::c1buf_org);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &Wave::v1buf[1][0]);
    
#pragma omp parallel for default(none) private(j, tmp_off) \
shared(idim_max_buf, trans, List::c1buf_org, List::c2_1, List::c2_2,Wave::v1buf, \
tmp_v0,nstate,one, Large::irght, Large::ilft, Large::ihfbit)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(List::c2_1, List::c2_2, List::c1buf_org[j], Large::irght, Large::ilft, Large::ihfbit, &tmp_off);
    zaxpy_(&nstate, &trans, &Wave::v1buf[j][0], &one, &tmp_v0[tmp_off][0], &one);
  }
}/*std::complex<double>  X_child_CisAit_spin_MPIdouble*/
/**
@brief Hopping term in Spin + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_CisAis_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0 /**< [out] Result v0 = H v1*/,
 std::complex<double> **tmp_v1 /**< [in] v0 = H v1*/
){
  int mask1, ibit1;
  mask1 = (int)Def::Tpow[org_isite1];
  ibit1 = (((long int)MP::myrank& mask1)/mask1)^(1-org_ispin1);
  if (ibit1 != 0) 
    zaxpy_long(Check::idim_max*nstate, tmp_trans, &tmp_v1[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_child_CisAis_spin_MPIdouble*/
/**
@brief Hopping term in Spin + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void X_GC_child_AisCis_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0 /**< [out] Result v0 = H v1*/,
  std::complex<double> **tmp_v1 /**< [in] v0 = H v1*/
){
  int mask1;
  int ibit1;
  mask1 = (int)Def::Tpow[org_isite1];
  ibit1 = (((long int)MP::myrank& mask1) / mask1) ^ (1 - org_ispin1);

  if (ibit1 == 0) {
    zaxpy_long(Check::idim_max*nstate, tmp_trans, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (ibit1 == 0)*/
}/*std::complex<double> X_GC_child_AisCis_spin_MPIdouble*/
