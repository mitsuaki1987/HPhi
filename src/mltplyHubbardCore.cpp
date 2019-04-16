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
@brief Functions for Hubbard hamiltonian (Core)
*/
#include <bitcalc.hpp>
#include "xsetmem.hpp"
#include "wrapperMPI.hpp"
#include "mltplyCommon.hpp"
#include "mltplyHubbardCore.hpp"

/******************************************************************************/
//[s] GetInfo functions
/******************************************************************************/

/**
@brief Compute mask for bit operation of hopping term.
@return Error-code, always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_general_hopp_GetInfo(
  struct BindStruct *X,//!<[inout]
  unsigned long int isite1,//!<[in] Site index
  unsigned long int isite2,//!<[in] Site index
  unsigned long int sigma1,//!<[in] Spin index
  unsigned long int sigma2//!<[in] Spin index
) {
  /**
  Compute mask for checking occupations of @f$(i_1,\sigma_1)@f$ (LargeList::is1_spin)
  and @f$(i_2,\sigma_2)@f$ (LargeList::is2_spin)
  */
  X->Large.is1_spin = X->Def.Tpow[2 * isite1 - 2 + sigma1];
  X->Large.is2_spin = X->Def.Tpow[2 * isite2 - 2 + sigma2];
  /**
  Compute mask for Fermion sign (LargeList::A_spin)
  */
  if (isite1 > isite2) {
    X->Large.A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
  }
  else if (isite1 < isite2) {
    X->Large.A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
  }
  else {
    if (sigma1 > sigma2) {
      X->Large.A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
    }
    else {
      X->Large.A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
    }
  }
  /**
  Compute mask for hopping (LargeList::isA_spin)
  */
  X->Large.isA_spin = X->Large.is1_spin + X->Large.is2_spin;
  return 0;
}/*int child_general_hopp_GetInfo*/
/**
@brief Compute mask for bit operation of general interaction term.
@return Error-code, always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_general_int_GetInfo(
  int iInterAll,//!<[in] It is not used
  struct BindStruct *X,//!<[inout]
  long unsigned int isite1,//!<[in] Site index
  long unsigned int isite2,//!<[in] Site index
  long unsigned int isite3,//!<[in] Site index
  long unsigned int isite4,//!<[in] Site index
  long unsigned int sigma1,//!<[in] Spin index
  long unsigned int sigma2,//!<[in] Spin index
  long unsigned int sigma3,//!<[in] Spin index
  long unsigned int sigma4,//!<[in] Spin index
  std::complex<double> tmp_V//!<[in] Coupling constant
) {
  long unsigned int is1_spin, is2_spin, is3_spin, is4_spin;
  long unsigned int A_spin, B_spin;
  long unsigned int isA_spin, isB_spin;
  /**
  Compute mask for checking occupations of @f$(i_1,\sigma_1)@f$ (LargeList::is1_spin)
  and @f$(i_2,\sigma_2)@f$ (LargeList::is2_spin)
  */
  is1_spin = X->Def.Tpow[2 * isite1 - 2 + sigma1];
  is2_spin = X->Def.Tpow[2 * isite2 - 2 + sigma2];
  /**
  Compute mask for Fermion sign (LargeList::A_spin)
  */
  if (isite1 > isite2) {
    A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
  }
  else if (isite2 > isite1) {
    A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
  }
  else {//isite1=isite2
    if (sigma1 > sigma2) {
      A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
    }
    else {
      A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
    }
  }
  /**
  Compute mask for checking occupations of @f$(i_3,\sigma_3)@f$ (LargeList::is3_spin)
  and @f$(i_4,\sigma_4)@f$ (LargeList::is4_spin)
  */
  is3_spin = X->Def.Tpow[2 * isite3 - 2 + sigma3];
  is4_spin = X->Def.Tpow[2 * isite4 - 2 + sigma4];
  /**
  Compute mask for Fermion sign (LargeList::B_spin)
  */
  if (isite3 > isite4) {
    B_spin = (X->Def.Tpow[2 * isite3 - 2 + sigma3] - X->Def.Tpow[2 * isite4 - 1 + sigma4]);
  }
  else if (isite3 < isite4) {
    B_spin = (X->Def.Tpow[2 * isite4 - 2 + sigma4] - X->Def.Tpow[2 * isite3 - 1 + sigma3]);
  }
  else {//isite3=isite4
    if (sigma3 > sigma4) {
      B_spin = (X->Def.Tpow[2 * isite3 - 2 + sigma3] - X->Def.Tpow[2 * isite4 - 1 + sigma4]);
    }
    else {
      B_spin = (X->Def.Tpow[2 * isite4 - 2 + sigma4] - X->Def.Tpow[2 * isite3 - 1 + sigma3]);
    }
  }
  /**
  Compute mask for hopping (LargeList::isA_spin, LargeList::isB_spin)
  */
  isA_spin = is1_spin + is2_spin;
  isB_spin = is3_spin + is4_spin;

  X->Large.is1_spin = is1_spin;
  X->Large.is2_spin = is2_spin;
  X->Large.is3_spin = is3_spin;
  X->Large.is4_spin = is4_spin;
  X->Large.isA_spin = isA_spin;
  X->Large.isB_spin = isB_spin;
  X->Large.A_spin = A_spin;
  X->Large.B_spin = B_spin;
  /**
  Copy coupling constant (LargeList::tmp_V)
  */
  X->Large.tmp_V = tmp_V;
  X->Large.isite1 = isite1;
  X->Large.isite2 = isite2;
  X->Large.isite3 = isite3;
  X->Large.isite4 = isite4;

  return 0;
}/*int child_general_int_GetInfo*/
/**
@brief Compute mask for bit operation of pairhop term.
@return Error-code, always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_pairhopp_GetInfo(
  int iPairHopp,//!<[in] Index of pairhopp interaction
  struct BindStruct *X//!<[inout]
) {
  int isite1 = X->Def.PairHopping[iPairHopp][0] + 1;
  int isite2 = X->Def.PairHopping[iPairHopp][1] + 1;
  /**
  Copy coupling constant (LargeList::tmp_J)
  */
  X->Large.tmp_J = X->Def.ParaPairHopping[iPairHopp];
  /**
  Compute mask for checking occupations of 
  @f$(i_1,\uparrow)@f$ (LargeList::is1_up), @f$(i_1,\downarrow)@f$ (LargeList::is1_down)
  @f$(i_2,\uparrow)@f$ (LargeList::is2_up), @f$(i_2,\downarrow)@f$ (LargeList::is2_down)
  */
  X->Large.is1_up = X->Def.Tpow[2 * isite1 - 2];
  X->Large.is1_down = X->Def.Tpow[2 * isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[2 * isite2 - 2];
  X->Large.is2_down = X->Def.Tpow[2 * isite2 - 1];

  return 0;
}/*int child_pairhopp_GetInfo*/
/**
@brief Compute mask for bit operation of exchange term.
@return Error-code, always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_exchange_GetInfo(
  int iExchange,//!<[in] Index of exchange interaction
  struct BindStruct *X//!<[inout]
) {
  int isite1 = X->Def.ExchangeCoupling[iExchange][0] + 1;
  int isite2 = X->Def.ExchangeCoupling[iExchange][1] + 1;
  /**
  Copy coupling constant (LargeList::tmp_J)
  */
  X->Large.tmp_J = -X->Def.ParaExchangeCoupling[iExchange];
  /**
  Compute mask for checking occupations of
  @f$(i_1,\uparrow)@f$ (LargeList::is1_up), @f$(i_1,\downarrow)@f$ (LargeList::is1_down)
  @f$(i_2,\uparrow)@f$ (LargeList::is2_up), @f$(i_2,\downarrow)@f$ (LargeList::is2_down)
  */
  X->Large.is1_up = X->Def.Tpow[2 * isite1 - 2];
  X->Large.is1_down = X->Def.Tpow[2 * isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[2 * isite2 - 2];
  X->Large.is2_down = X->Def.Tpow[2 * isite2 - 1];

  return 0;
}/*int child_exchange_GetInfo*/

/******************************************************************************/
//[e] GetInfo functions
/******************************************************************************/

/******************************************************************************/
//[s] core routines
/******************************************************************************/
/**
@brief Operation of @f$t c_{i\sigma}^\dagger c_{i\sigma}@f$ (Grandcanonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void GC_CisAis(
  long unsigned int j,//!<[in] Index of element of wavefunction
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of @f$(i \sigma)@f$
  std::complex<double> tmp_trans//!<[in] Transfer integral
) {
  long unsigned int A_ibit_tmp;
  long unsigned int list_1_j;
  std::complex<double> dmv;
  int one = 1;

  list_1_j = j - 1;
  A_ibit_tmp = (list_1_j & is1_spin) / is1_spin;
  dmv = tmp_trans * (std::complex<double>)A_ibit_tmp;
  zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
}/*std::complex<double> GC_CisAis*/
/**
@brief Operation of @f$t c_{i\sigma} c_{i\sigma}^\dagger@f$ (Grandcanonical)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_AisCis(
  long unsigned int j,//!<[in] Index of element of wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of @f$(i \sigma)@f$
  std::complex<double> tmp_trans//!<[in] Transfer integral
) {
  long unsigned int A_ibit_tmp;
  long unsigned int list_1_j;
  std::complex<double> dmv;
  int one = 1;

  list_1_j = j - 1;
  A_ibit_tmp = (list_1_j & is1_spin) / is1_spin;
  dmv = tmp_trans * (std::complex<double>)(1 - A_ibit_tmp);
  zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
}/*std::complex<double> GC_AisCis*/
/**
@brief @f$c_{is}\\dagger c_{is}@f$ term in Hubbard (canonical) 
@return
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int X_CisAis(
  long unsigned int list_1_j,
  struct BindStruct *X,
  long unsigned int is1_spin
) {
  int A_ibit_tmp;

  // off = j
  A_ibit_tmp = (list_1_j & is1_spin) / is1_spin;
  return A_ibit_tmp;
}
/**
@brief @f$c_{is}^\dagger c_{jt}@f$ term for canonical Hubbard
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void CisAjt(
  long unsigned int j,//!<[in] Index of wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f$v_0 = H v_1@f$
  std::complex<double> **tmp_v1,//!<[in] Vector to be producted
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  std::complex<double> tmp_V//!<[in] Hopping integral
) {
  long unsigned int ibit_tmp_1, ibit_tmp_2;
  long unsigned int bit, iexchg, off;
  int sgn;
  std::complex<double> dmv;
  int one = 1;

  ibit_tmp_1 = (list_1[j] & is1_spin);
  ibit_tmp_2 = (list_1[j] & is2_spin);
  if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
    bit = list_1[j] & diff_spin;
    SgnBit(bit, &sgn); // Fermion sign
    iexchg = list_1[j] ^ sum_spin;

    if(GetOffComp(list_2_1, list_2_2, iexchg, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &off)==FALSE){
      return;
    }
    dmv = (std::complex<double>)sgn * tmp_V;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[off][0], &one);
  }
  else {
    return;
  }
}
/**
@brief @f$c_{is}^\dagger c_{jt}@f$ term for grandcanonical Hubbard
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_CisAjt(
  long unsigned int j,//!<[in] Index of wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[in] @f$v_0 = H v_1@f$
  std::complex<double> **tmp_v1,//!<[in]Vector to be producted
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Hopping
  long unsigned int *tmp_off//!<[in] Index of wavefunction of final state
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1, ibit_tmp_2;
  long unsigned int bit;
  int sgn;
  std::complex<double> dmv;

  list_1_j = j - 1;
  ibit_tmp_1 = (list_1_j & is1_spin);
  ibit_tmp_2 = (list_1_j & is2_spin);
  *tmp_off = 0;
  int one = 1;

  if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
    bit = list_1_j & diff_spin;
    SgnBit(bit, &sgn); // Fermion sign
    list_1_off = list_1_j ^ sum_spin;
    *tmp_off = list_1_off;
    dmv = (std::complex<double>)sgn * tmp_V;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[list_1_off + 1][0], &one);
  }
  else {
    return;
  }
}/*std::complex<double> GC_CisAjt*/
/**
@brief Compute index of wavefunction of final state
@return Fermion sign
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int X_CisAjt(
  long unsigned int list_1_j,//!<[in] Similer to ::list_1 ?
  struct BindStruct *X,//!<[in]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  long unsigned int *tmp_off//!<[in] Index of wavefunction of final state
) {
  long unsigned int off;
  int sgn = 1;

  sgn = X_GC_CisAjt(list_1_j, X, is1_spin, is2_spin, sum_spin, diff_spin, tmp_off);
  if (sgn != 0) {
    if(GetOffComp(list_2_1, list_2_2, *tmp_off, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &off)!=TRUE){
      *tmp_off = 0;
      return 0;
    }
    *tmp_off = off;
    return sgn;
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int X_CisAjt*/
/**
@brief Compute index of wavefunction of final state
@return Fermion sign
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int X_GC_CisAjt(
  long unsigned int list_1_j,//!<[in] ::list_1 ?
  struct BindStruct *X,//!<[in]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  long unsigned int *tmp_off//!<[out] Index of wavefunction of final state
) {
  long unsigned int ibit_tmp_1, ibit_tmp_2;
  long unsigned int bit, off;
  int sgn = 1;

  ibit_tmp_1 = (list_1_j & is1_spin);
  ibit_tmp_2 = (list_1_j & is2_spin);

  if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
    bit = list_1_j & diff_spin;
    SgnBit(bit, &sgn); // Fermion sign
    off = list_1_j ^ sum_spin;
    *tmp_off = off;
    return sgn; // pm 1
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}
/******************************************************************************/
//[e] core routines
/******************************************************************************/

/******************************************************************************/
//[s] child element functions
/******************************************************************************/
/**
@brief Compute exchange term of canonical-Hubbard
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_exchange_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f$v_0 = H v_1@f$
  std::complex<double> **tmp_v1,//!<[in] Vector to be producted
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[off] Index of wavefunction of final state
) {
  long unsigned int off;
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  std::complex<double> tmp_J = X->Large.tmp_J;
  int one = 1;

  ibit1_up = list_1[j] & is1_up;
  ibit2_up = list_1[j] & is2_up;
  ibit1_down = list_1[j] & is1_down;
  ibit2_down = list_1[j] & is2_down;

  if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {
    iexchg = list_1[j] - (is1_down + is2_up);
    iexchg += (is1_up + is2_down);
    if(GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
      return;
    }
    *tmp_off = off;
    dmv = tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[off][0], &one);
  }
  else if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
    iexchg = list_1[j] - (is1_up + is2_down);
    iexchg += (is1_down + is2_up);
    if(GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
      return;
    }
    *tmp_off = off;
    dmv = tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[off][0], &one);
  }
}/*std::complex<double> child_exchange_element*/
/**
@brief Compute pairhopp term of canonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_pairhopp_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int off;
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  std::complex<double> tmp_J = X->Large.tmp_J;
  int one = 1;

  ibit1_up = list_1[j] & is1_up;
  ibit2_up = list_1[j] & is2_up;
  ibit1_down = list_1[j] & is1_down;
  ibit2_down = list_1[j] & is2_down;

  if (ibit1_up == 0 && ibit1_down == 0 && ibit2_up != 0 && ibit2_down != 0) {
    iexchg = list_1[j] - (is2_up + is2_down);
    iexchg += (is1_up + is1_down);

    if(GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
      return;
    }
    *tmp_off = off;
    dmv = tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[off][0], &one);
  }
}/*std::complex<double> child_pairhopp_element*/
/**
@brief Compute exchange term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_exchange_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int list_1_j, list_1_off;
  std::complex<double> tmp_J = X->Large.tmp_J;
  int one = 1;

  list_1_j = j - 1;
  ibit1_up = list_1_j & is1_up;
  ibit2_up = list_1_j & is2_up;
  ibit1_down = list_1_j & is1_down;
  ibit2_down = list_1_j & is2_down;

  if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {

    iexchg = list_1_j - (is1_down + is2_up);
    iexchg += (is1_up + is2_down);
    list_1_off = iexchg;
    *tmp_off = list_1_off;

    dmv = tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[list_1_off + 1][0], &one);
  }
  else if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
    iexchg = list_1_j - (is1_up + is2_down);
    iexchg += (is1_down + is2_up);
    list_1_off = iexchg;
    *tmp_off = list_1_off;

    dmv = tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[list_1_off + 1][0], &one);
  }
}/*std::complex<double> GC_child_exchange_element*/
/**
@brief Compute pairhopp term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_pairhopp_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int list_1_j, list_1_off;
  std::complex<double> tmp_J = X->Large.tmp_J;
  int one = 1;

  list_1_j = j - 1;

  ibit1_up = list_1_j & is1_up;

  ibit2_up = list_1_j & is2_up;

  ibit1_down = list_1_j & is1_down;

  ibit2_down = list_1_j & is2_down;

  if (ibit1_up == 0 && ibit1_down == 0 && ibit2_up != 0 && ibit2_down != 0) {
    iexchg = list_1_j - (is2_up + is2_down);
    iexchg += (is1_up + is1_down);
    list_1_off = iexchg;
    *tmp_off = list_1_off;
    dmv = tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[list_1_off + 1][0], &one);
  }
}
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$
term of canonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_CisAisCisAis_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAis(list_1[j], X, isite3);
  tmp_sgn *= X_CisAis(list_1[j], X, isite1);
  dmv = tmp_V * (std::complex<double>)tmp_sgn;
  zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
}/*std::complex<double> child_CisAisCisAis_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of canonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_CisAisCjtAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAjt(list_1[j], X, isite3, isite4, Bsum, Bdiff, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_CisAis(list_1[*tmp_off], X, isite1);
    if (tmp_sgn != 0) {
      dmv = tmp_V * (std::complex<double>)tmp_sgn;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off][0], &one);
    }
  }
}/*std::complex<double> child_CisAisCjtAku_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of canonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_CisAjtCkuAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAis(list_1[j], X, isite3);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_CisAjt(list_1[j], X, isite1, isite2, Asum, Adiff, tmp_off);
    if (tmp_sgn != 0) {
      dmv = tmp_V * (std::complex<double>)tmp_sgn;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off][0], &one);
    }
  }
}/*std::complex<double> child_CisAjtCkuAku_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of canonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_CisAjtCkuAlv_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long unsigned int tmp_off_1;
  int one = 1;

  std::complex<double> dmv;
  tmp_sgn = X_GC_CisAjt(list_1[j], X, isite3, isite4, Bsum, Bdiff, &tmp_off_1);

  if (tmp_sgn != 0) {
    tmp_sgn *= X_CisAjt(tmp_off_1, X, isite1, isite2, Asum, Adiff, tmp_off_2);
    if (tmp_sgn != 0) {
      dmv = tmp_V * (std::complex<double>)tmp_sgn;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off_2][0], &one);
    }
  }
}/*std::complex<double> child_CisAjtCkuAlv_element*/
//[s] Grand Canonical
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$
term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAisCisAis_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAis(j - 1, X, isite3);
  tmp_sgn *= X_CisAis(j - 1, X, isite1);
  if (tmp_sgn != 0) {
    dmv = tmp_V * (std::complex<double>)tmp_sgn;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }
}/*std::complex<double> GC_child_CisAisCisAis_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAisCjtAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_GC_CisAjt((j - 1), X, isite3, isite4, Bsum, Bdiff, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_CisAis(*tmp_off, X, isite1);
    if (tmp_sgn != 0) {
      dmv = tmp_V * (std::complex<double>)tmp_sgn;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off + 1][0], &one);
    }
  }
}/*std::complex<double> GC_child_CisAisCjtAku_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAjtCkuAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAis((j - 1), X, isite3);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_GC_CisAjt((j - 1), X, isite1, isite2, Asum, Adiff, tmp_off);
    if (tmp_sgn != 0) {
      dmv = tmp_V * (std::complex<double>)tmp_sgn;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off + 1][0], &one);
    }/*if (tmp_sgn != 0)*/
  }/*if (tmp_sgn != 0)*/
}/*std::complex<double> GC_child_CisAjtCkuAku_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_CisAjtCkuAlv_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long unsigned int tmp_off_1;
  std::complex<double> dmv;
  int one = 1;

  tmp_sgn = X_GC_CisAjt((j - 1), X, isite3, isite4, Bsum, Bdiff, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_GC_CisAjt(tmp_off_1, X, isite1, isite2, Asum, Adiff, tmp_off_2);
    if (tmp_sgn != 0) {
      dmv = tmp_V * (std::complex<double>)tmp_sgn;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off_2 + 1][0], &one);
    }
  }
}/*std::complex<double> GC_child_CisAjtCkuAlv_element*/
//[e] Grand Canonical
/**
@brief Compute @f$c_{is}^\dagger@f$
term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void GC_Cis(
  long unsigned int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long unsigned int is1_spin,//!<[in] Bit mask 
  std::complex<double> tmp_V,//!<[in] Coupling constant
  long unsigned int *tmp_off//!<[in] Index of final wavefunction
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;
  std::complex<double> dmv;
  int one = 1;

  list_1_j = j - 1;

  ibit_tmp_1 = (list_1_j & is1_spin);
  // is1_spin >= 1
  // is1_spin = Tpow[2*isite + ispin]

  *tmp_off = 0;

  if (ibit_tmp_1 == 0) {
    // able to create an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef __MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j | is1_spin; // OR
    *tmp_off = list_1_off;
    dmv = (std::complex<double>)(ipsgn * sgn) * tmp_V;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[list_1_off + 1][0], &one);
  }
  else {
    return;
  }
}/*std::complex<double> GC_Cis*/
/**
@brief Compute @f$c_{jt}@f$
term of grandcanonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void GC_Ajt(
  long unsigned int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long unsigned int is1_spin,//!<[in] Bit mask
  std::complex<double> tmp_V,//!<[in] Coupling constant
  long unsigned int *tmp_off//!<[in] Index of final wavefunction
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;
  std::complex<double> dmv;
  int one = 1;

  list_1_j = j - 1;

  ibit_tmp_1 = (list_1_j & is1_spin);
  // is1_spin >= 1

  *tmp_off = 0;

  if (ibit_tmp_1 == is1_spin) {
    // able to create an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef __MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j ^ is1_spin;
    *tmp_off = list_1_off;
    dmv = (std::complex<double>)(ipsgn * sgn) * tmp_V;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[list_1_off + 1][0], &one);
  }
  else {
    return;
  }
}/*std::complex<double> GC_Ajt*/
/**
@brief Compute index of final wavefunction associatesd to @f$c_{is}^\dagger@f$
term of canonical Hubbard system
@return 1 if electron (is) absent, 0 if not.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
int X_Cis(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int is1_spin,//!<[in] Bit mask
  long unsigned int *tmp_off,//!<[out] Index of final wavefunction
  long unsigned int *list_1_org,//!<[in] Similar to ::list_1
  long unsigned int *list_2_1_target,//!<[in] Similar to ::list_2_1
  long unsigned int *list_2_2_target,//!<[in] Similar to ::list_2_2
  long unsigned int _irght,//!<[in] Similar to LargeList::irght
  long unsigned int _ilft,//!<[in] Similar to LargeList::ilft
  long unsigned int _ihfbit//!<[in] Similar to LargeList::ihfbit
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;

  list_1_j = list_1_org[j];

  ibit_tmp_1 = (list_1_j & is1_spin);
  // is1_spin >= 1
  // is1_spin = Tpow[2*isite + ispin]

  *tmp_off = 0;

  if (ibit_tmp_1 == 0) {
    // able to create an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef __MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j | is1_spin; // OR

    if(GetOffComp(list_2_1_target, list_2_2_target, list_1_off, _irght, _ilft, _ihfbit, tmp_off)!=TRUE){
      *tmp_off=0;
      return 0;
    }
    sgn *= ipsgn;
    return (sgn);
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int X_Cis*/
/**
@brief Compute index of final wavefunction associatesd to @f$c_{jt}@f$
term of canonical Hubbard system
@return 1 if electron (jt) exist, 0 if not.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
int X_Ajt(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int is1_spin,//!<[in] Bit mask
  long unsigned int *tmp_off,//!<[out] Index of final wavefunction
  long unsigned int *list_1_org,//!<[in] Similar to ::list_1
  long unsigned int *list_2_1_target,//!<[in] Similar to ::list_2_1
  long unsigned int *list_2_2_target,//!<[in] Similar to ::list_2_2
  long unsigned int _irght,//!<[in] Similar to LargeList::irght
  long unsigned int _ilft,//!<[in] Similar to LargeList::ilft
  long unsigned int _ihfbit//!<[in] Similar to LargeList::ihfbit
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;

  list_1_j = list_1_org[j];

  ibit_tmp_1 = (list_1_j & is1_spin);
// is1_spin >= 1
// is1_spin = Tpow[2*isite + ispin]

  *tmp_off = 0;
  if (ibit_tmp_1 != 0) {
    // able to delete an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef __MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j ^ is1_spin;
    if(GetOffComp(list_2_1_target, list_2_2_target, list_1_off, _irght, _ilft, _ihfbit, tmp_off)!=TRUE){
      *tmp_off=0;
      return 0;
    }
    sgn *= ipsgn;
    return(sgn);
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*std::complex<double> X_Ajt*/

/******************************************************************************/
//[e] child element functions
/******************************************************************************/
