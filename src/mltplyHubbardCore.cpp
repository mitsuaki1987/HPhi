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
#include "global.hpp"

/******************************************************************************/
//[s] GetInfo functions
/******************************************************************************/

/**
@brief Compute mask for bit operation of hopping term.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_general_hopp_GetInfo(
  long int isite1,//!<[in] Site index
  long int isite2,//!<[in] Site index
  long int sigma1,//!<[in] Spin index
  long int sigma2//!<[in] Spin index
) {
  /**
  Compute mask for checking occupations of @f$(i_1,\sigma_1)@f$ (LargeList::is1_spin)
  and @f$(i_2,\sigma_2)@f$ (LargeList::is2_spin)
  */
  Large::is1_spin = Def::Tpow[2 * isite1 - 2 + sigma1];
  Large::is2_spin = Def::Tpow[2 * isite2 - 2 + sigma2];
  /**
  Compute mask for Fermion sign (LargeList::A_spin)
  */
  if (isite1 > isite2) {
    Large::A_spin = (Def::Tpow[2 * isite1 - 2 + sigma1] - Def::Tpow[2 * isite2 - 1 + sigma2]);
  }
  else if (isite1 < isite2) {
    Large::A_spin = (Def::Tpow[2 * isite2 - 2 + sigma2] - Def::Tpow[2 * isite1 - 1 + sigma1]);
  }
  else {
    if (sigma1 > sigma2) {
      Large::A_spin = (Def::Tpow[2 * isite1 - 2 + sigma1] - Def::Tpow[2 * isite2 - 1 + sigma2]);
    }
    else {
      Large::A_spin = (Def::Tpow[2 * isite2 - 2 + sigma2] - Def::Tpow[2 * isite1 - 1 + sigma1]);
    }
  }
  /**
  Compute mask for hopping (LargeList::isA_spin)
  */
  Large::isA_spin = Large::is1_spin + Large::is2_spin;
}/*int child_general_hopp_GetInfo*/
/**
@brief Compute mask for bit operation of general interaction term.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_general_int_GetInfo(
  long int isite1,//!<[in] Site index
  long int isite2,//!<[in] Site index
  long int isite3,//!<[in] Site index
  long int isite4,//!<[in] Site index
  long int sigma1,//!<[in] Spin index
  long int sigma2,//!<[in] Spin index
  long int sigma3,//!<[in] Spin index
  long int sigma4,//!<[in] Spin index
  std::complex<double> tmp_V//!<[in] Coupling constant
) {
  long int is1_spin, is2_spin, is3_spin, is4_spin;
  long int A_spin, B_spin;
  long int isA_spin, isB_spin;
  /**
  Compute mask for checking occupations of @f$(i_1,\sigma_1)@f$ (LargeList::is1_spin)
  and @f$(i_2,\sigma_2)@f$ (LargeList::is2_spin)
  */
  is1_spin = Def::Tpow[2 * isite1 - 2 + sigma1];
  is2_spin = Def::Tpow[2 * isite2 - 2 + sigma2];
  /**
  Compute mask for Fermion sign (LargeList::A_spin)
  */
  if (isite1 > isite2) {
    A_spin = (Def::Tpow[2 * isite1 - 2 + sigma1] - Def::Tpow[2 * isite2 - 1 + sigma2]);
  }
  else if (isite2 > isite1) {
    A_spin = (Def::Tpow[2 * isite2 - 2 + sigma2] - Def::Tpow[2 * isite1 - 1 + sigma1]);
  }
  else {//isite1=isite2
    if (sigma1 > sigma2) {
      A_spin = (Def::Tpow[2 * isite1 - 2 + sigma1] - Def::Tpow[2 * isite2 - 1 + sigma2]);
    }
    else {
      A_spin = (Def::Tpow[2 * isite2 - 2 + sigma2] - Def::Tpow[2 * isite1 - 1 + sigma1]);
    }
  }
  /**
  Compute mask for checking occupations of @f$(i_3,\sigma_3)@f$ (LargeList::is3_spin)
  and @f$(i_4,\sigma_4)@f$ (LargeList::is4_spin)
  */
  is3_spin = Def::Tpow[2 * isite3 - 2 + sigma3];
  is4_spin = Def::Tpow[2 * isite4 - 2 + sigma4];
  /**
  Compute mask for Fermion sign (LargeList::B_spin)
  */
  if (isite3 > isite4) {
    B_spin = (Def::Tpow[2 * isite3 - 2 + sigma3] - Def::Tpow[2 * isite4 - 1 + sigma4]);
  }
  else if (isite3 < isite4) {
    B_spin = (Def::Tpow[2 * isite4 - 2 + sigma4] - Def::Tpow[2 * isite3 - 1 + sigma3]);
  }
  else {//isite3=isite4
    if (sigma3 > sigma4) {
      B_spin = (Def::Tpow[2 * isite3 - 2 + sigma3] - Def::Tpow[2 * isite4 - 1 + sigma4]);
    }
    else {
      B_spin = (Def::Tpow[2 * isite4 - 2 + sigma4] - Def::Tpow[2 * isite3 - 1 + sigma3]);
    }
  }
  /**
  Compute mask for hopping (LargeList::isA_spin, LargeList::isB_spin)
  */
  isA_spin = is1_spin + is2_spin;
  isB_spin = is3_spin + is4_spin;

  Large::is1_spin = is1_spin;
  Large::is2_spin = is2_spin;
  Large::is3_spin = is3_spin;
  Large::is4_spin = is4_spin;
  Large::isA_spin = isA_spin;
  Large::isB_spin = isB_spin;
  Large::A_spin = A_spin;
  Large::B_spin = B_spin;
  /**
  Copy coupling constant (LargeList::tmp_V)
  */
  Large::tmp_V = tmp_V;
  Large::isite1 = isite1;
  Large::isite2 = isite2;
  Large::isite3 = isite3;
  Large::isite4 = isite4;
}/*int child_general_int_GetInfo*/
/**
@brief Compute mask for bit operation of pairhop term.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_pairhopp_GetInfo(
  int iPairHopp//!<[in] Index of pairhopp interaction
) {
  int isite1 = Def::PairHopping[iPairHopp][0] + 1;
  int isite2 = Def::PairHopping[iPairHopp][1] + 1;
  /**
  Copy coupling constant (LargeList::tmp_J)
  */
  Large::tmp_J = Def::ParaPairHopping[iPairHopp];
  /**
  Compute mask for checking occupations of 
  @f$(i_1,\uparrow)@f$ (LargeList::is1_up), @f$(i_1,\downarrow)@f$ (LargeList::is1_down)
  @f$(i_2,\uparrow)@f$ (LargeList::is2_up), @f$(i_2,\downarrow)@f$ (LargeList::is2_down)
  */
  Large::is1_up = Def::Tpow[2 * isite1 - 2];
  Large::is1_down = Def::Tpow[2 * isite1 - 1];
  Large::is2_up = Def::Tpow[2 * isite2 - 2];
  Large::is2_down = Def::Tpow[2 * isite2 - 1];
}/*int child_pairhopp_GetInfo*/
/**
@brief Compute mask for bit operation of exchange term.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_exchange_GetInfo(
  int iExchange//!<[in] Index of exchange interaction
) {
  int isite1 = Def::ExchangeCoupling[iExchange][0] + 1;
  int isite2 = Def::ExchangeCoupling[iExchange][1] + 1;
  /**
  Copy coupling constant (LargeList::tmp_J)
  */
  Large::tmp_J = -Def::ParaExchangeCoupling[iExchange];
  /**
  Compute mask for checking occupations of
  @f$(i_1,\uparrow)@f$ (LargeList::is1_up), @f$(i_1,\downarrow)@f$ (LargeList::is1_down)
  @f$(i_2,\uparrow)@f$ (LargeList::is2_up), @f$(i_2,\downarrow)@f$ (LargeList::is2_down)
  */
  Large::is1_up = Def::Tpow[2 * isite1 - 2];
  Large::is1_down = Def::Tpow[2 * isite1 - 1];
  Large::is2_up = Def::Tpow[2 * isite2 - 2];
  Large::is2_down = Def::Tpow[2 * isite2 - 1];
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
void mltply::HubbardGC::CisAis(
  long int j,//!<[in] Index of element of wavefunction
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1,//!<[in] Input producted vector
  long int is1_spin,//!<[in] Mask for occupation of @f$(i \sigma)@f$
  std::complex<double> tmp_trans//!<[in] Transfer integral
) {
  long int A_ibit_tmp;
  long int list_1_j;
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
void mltply::HubbardGC::AisCis(
  long int j,//!<[in] Index of element of wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1,//!<[in] Input producted vector
  long int is1_spin,//!<[in] Mask for occupation of @f$(i \sigma)@f$
  std::complex<double> tmp_trans//!<[in] Transfer integral
) {
  long int A_ibit_tmp;
  long int list_1_j;
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
int mltply::Hubbard::X_CisAis(
  long int list_1_j,
  long int is1_spin
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
void mltply::Hubbard::CisAjt(
  long int j,//!<[in] Index of wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f$v_0 = H v_1@f$
  std::complex<double> **tmp_v1,//!<[in] Vector to be producted
  //!<[inout]
  long int is1_spin,//!<[in] Mask for occupation of (is)
  long int is2_spin,//!<[in] Mask for occupation of (jt)
  long int sum_spin,//!<[in] Mask for hopping
  long int diff_spin,//!<[in] Mask for Fermion sign
  std::complex<double> tmp_V//!<[in] Hopping integral
) {
  long int ibit_tmp_1, ibit_tmp_2;
  long int bit, iexchg, off;
  int sgn;
  std::complex<double> dmv;
  int one = 1;

  ibit_tmp_1 = (List::c1[j] & is1_spin);
  ibit_tmp_2 = (List::c1[j] & is2_spin);
  if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
    bit = List::c1[j] & diff_spin;
    SgnBit(bit, &sgn); // Fermion sign
    iexchg = List::c1[j] ^ sum_spin;

    if(GetOffComp(List::c2_1, List::c2_2, iexchg, Large::irght, Large::ilft, Large::ihfbit, &off)==FALSE){
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
void mltply::HubbardGC::CisAjt(
  long int j,//!<[in] Index of wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[in] @f$v_0 = H v_1@f$
  std::complex<double> **tmp_v1,//!<[in]Vector to be producted
  long int is1_spin,//!<[in] Mask for occupation of (is)
  long int is2_spin,//!<[in] Mask for occupation of (jt)
  long int sum_spin,//!<[in] Mask for hopping
  long int diff_spin,//!<[in] Mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Hopping
  long int *tmp_off//!<[in] Index of wavefunction of final state
) {
  long int list_1_j, list_1_off;
  long int ibit_tmp_1, ibit_tmp_2;
  long int bit;
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
int mltply::Hubbard::X_CisAjt(
  long int list_1_j,//!<[in] Similer to ::List::c1 ?
  //!<[in]
  long int is1_spin,//!<[in] Mask for occupation of (is)
  long int is2_spin,//!<[in] Mask for occupation of (jt)
  long int sum_spin,//!<[in] Mask for hopping
  long int diff_spin,//!<[in] Mask for Fermion sign
  long int *tmp_off//!<[in] Index of wavefunction of final state
) {
  long int off;
  int sgn = 1;

  sgn = mltply::HubbardGC::X_CisAjt(list_1_j, is1_spin, is2_spin, sum_spin, diff_spin, tmp_off);
  if (sgn != 0) {
    if(GetOffComp(List::c2_1, List::c2_2, *tmp_off, Large::irght, Large::ilft, Large::ihfbit, &off)!=TRUE){
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
int mltply::HubbardGC::X_CisAjt(
  long int list_1_j,//!<[in] ::List::c1 ?
  long int is1_spin,//!<[in] Mask for occupation of (is)
  long int is2_spin,//!<[in] Mask for occupation of (jt)
  long int sum_spin,//!<[in] Mask for hopping
  long int diff_spin,//!<[in] Mask for Fermion sign
  long int *tmp_off//!<[out] Index of wavefunction of final state
) {
  long int ibit_tmp_1, ibit_tmp_2;
  long int bit, off;
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
void mltply::Hubbard::exchange_element(
  long int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] @f$v_0 = H v_1@f$
  std::complex<double> **tmp_v1,//!<[in] Vector to be producted
  //!<[inout]
  long int *tmp_off//!<[off] Index of wavefunction of final state
) {
  long int off;
  long int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long int iexchg;
  long int is1_up = Large::is1_up;
  long int is2_up = Large::is2_up;
  long int is1_down = Large::is1_down;
  long int is2_down = Large::is2_down;
  long int irght = Large::irght;
  long int ilft = Large::ilft;
  long int ihfbit = Large::ihfbit;
  std::complex<double> tmp_J = Large::tmp_J;
  int one = 1;

  ibit1_up = List::c1[j] & is1_up;
  ibit2_up = List::c1[j] & is2_up;
  ibit1_down = List::c1[j] & is1_down;
  ibit2_down = List::c1[j] & is2_down;

  if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {
    iexchg = List::c1[j] - (is1_down + is2_up);
    iexchg += (is1_up + is2_down);
    if(GetOffComp(List::c2_1, List::c2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
      return;
    }
    *tmp_off = off;
    dmv = tmp_J;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[off][0], &one);
  }
  else if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
    iexchg = List::c1[j] - (is1_up + is2_down);
    iexchg += (is1_down + is2_up);
    if(GetOffComp(List::c2_1, List::c2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
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
void mltply::Hubbard::pairhopp_element(
  long int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  long int off;
  long int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long int iexchg;
  long int is1_up = Large::is1_up;
  long int is2_up = Large::is2_up;
  long int is1_down = Large::is1_down;
  long int is2_down = Large::is2_down;
  long int irght = Large::irght;
  long int ilft = Large::ilft;
  long int ihfbit = Large::ihfbit;
  std::complex<double> tmp_J = Large::tmp_J;
  int one = 1;

  ibit1_up = List::c1[j] & is1_up;
  ibit2_up = List::c1[j] & is2_up;
  ibit1_down = List::c1[j] & is1_down;
  ibit2_down = List::c1[j] & is2_down;

  if (ibit1_up == 0 && ibit1_down == 0 && ibit2_up != 0 && ibit2_down != 0) {
    iexchg = List::c1[j] - (is2_up + is2_down);
    iexchg += (is1_up + is1_down);

    if(GetOffComp(List::c2_1, List::c2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
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
void mltply::HubbardGC::exchange_element(
  long int j,//!<[in] Index of initial wavefunction
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  long int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long int iexchg;
  long int is1_up = Large::is1_up;
  long int is2_up = Large::is2_up;
  long int is1_down = Large::is1_down;
  long int is2_down = Large::is2_down;
  long int list_1_j, list_1_off;
  std::complex<double> tmp_J = Large::tmp_J;
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
void mltply::HubbardGC::pairhopp_element(
  long int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  //!<[inout]
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  long int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> dmv;
  long int iexchg;
  long int is1_up = Large::is1_up;
  long int is2_up = Large::is2_up;
  long int is1_down = Large::is1_down;
  long int is2_down = Large::is2_down;
  long int list_1_j, list_1_off;
  std::complex<double> tmp_J = Large::tmp_J;
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
void mltply::Hubbard::CisAisCisAis_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite3,//!<[in] Site 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Wavefunction to be multiplied
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAis(List::c1[j], isite3);
  tmp_sgn *= X_CisAis(List::c1[j], isite1);
  dmv = tmp_V * (std::complex<double>)tmp_sgn;
  zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
}/*std::complex<double> child_CisAisCisAis_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of canonical Hubbard system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::CisAisCjtAku_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite3,//!<[in] Site 3
  long int isite4,//!<[in] Site 4
  long int Bsum,//!<[in] Bit mask for hopping
  long int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  //!<[inout]
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAjt(List::c1[j], isite3, isite4, Bsum, Bdiff, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_CisAis(List::c1[*tmp_off], isite1);
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
void mltply::Hubbard::CisAjtCkuAku_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite2,//!<[in] Site 2
  long int isite3,//!<[in] Site 3
  long int Asum,//!<[in] Bit mask for hopping
  long int Adiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  //!<[inout]
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = X_CisAis(List::c1[j], isite3);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_CisAjt(List::c1[j], isite1, isite2, Asum, Adiff, tmp_off);
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
void mltply::Hubbard::CisAjtCkuAlv_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite2,//!<[in] Site 2
  long int isite3,//!<[in] Site 3
  long int isite4,//!<[in] Site 4
  long int Asum,//!<[in] Bit mask for hopping
  long int Adiff,//!<[in] Bit mask for Fermion sign
  long int Bsum,//!<[in] Bit mask for hopping
  long int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  //!<[inout]
  long int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int tmp_off_1;
  int one = 1;

  std::complex<double> dmv;
  tmp_sgn = mltply::HubbardGC::X_CisAjt(List::c1[j], isite3, isite4, Bsum, Bdiff, &tmp_off_1);

  if (tmp_sgn != 0) {
    tmp_sgn *= X_CisAjt(tmp_off_1, isite1, isite2, Asum, Adiff, tmp_off_2);
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
void mltply::HubbardGC::CisAisCisAis_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite3,//!<[in] Site 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Wavefunction to be multiplied
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = mltply::Hubbard::X_CisAis(j - 1, isite3);
  tmp_sgn *= mltply::Hubbard::X_CisAis(j - 1, isite1);
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
void mltply::HubbardGC::CisAisCjtAku_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite3,//!<[in] Site 3
  long int isite4,//!<[in] Site 4
  long int Bsum,//!<[in] Bit mask for hopping
  long int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = mltply::HubbardGC::X_CisAjt((j - 1), isite3, isite4, Bsum, Bdiff, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= mltply::Hubbard::X_CisAis(*tmp_off, isite1);
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
void mltply::HubbardGC::CisAjtCkuAku_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite2,//!<[in] Site 2
  long int isite3,//!<[in] Site 3
  long int Asum,//!<[in] Bit mask for hopping
  long int Adiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = mltply::Hubbard::X_CisAis(j - 1, isite3);
  if (tmp_sgn != 0) {
    tmp_sgn *= mltply::HubbardGC::X_CisAjt(j - 1, isite1, isite2, Asum, Adiff, tmp_off);
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
void mltply::HubbardGC::CisAjtCkuAlv_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isite1,//!<[in] Site 1
  long int isite2,//!<[in] Site 2
  long int isite3,//!<[in] Site 3
  long int isite4,//!<[in] Site 4
  long int Asum,//!<[in] Bit mask for hopping
  long int Adiff,//!<[in] Bit mask for Fermion sign
  long int Bsum,//!<[in] Bit mask for hopping
  long int Bdiff,//!<[in] Bit mask for Fermion sign
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int tmp_off_1;
  std::complex<double> dmv;
  int one = 1;

  tmp_sgn = mltply::HubbardGC::X_CisAjt(j - 1, isite3, isite4, Bsum, Bdiff, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= mltply::HubbardGC::X_CisAjt(tmp_off_1, isite1, isite2, Asum, Adiff, tmp_off_2);
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
void mltply::HubbardGC::Cis(
  long int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int is1_spin,//!<[in] Bit mask 
  std::complex<double> tmp_V,//!<[in] Coupling constant
  long int *tmp_off//!<[in] Index of final wavefunction
) {
  long int list_1_j, list_1_off;
  long int ibit_tmp_1;
  long int bit;
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
    SgnBit(MP::myrank, &ipsgn); // Fermion sign
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
void mltply::HubbardGC::Ajt(
  long int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int is1_spin,//!<[in] Bit mask
  std::complex<double> tmp_V,//!<[in] Coupling constant
  long int *tmp_off//!<[in] Index of final wavefunction
) {
  long int list_1_j, list_1_off;
  long int ibit_tmp_1;
  long int bit;
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
    SgnBit(MP::myrank, &ipsgn); // Fermion sign
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
int mltply::Hubbard::X_Cis(
  long int j,//!<[in] Index of initial wavefunction
  long int is1_spin,//!<[in] Bit mask
  long int *tmp_off,//!<[out] Index of final wavefunction
  long int *list_1_org,//!<[in] Similar to ::List::c1
  long int *list_2_1_target,//!<[in] Similar to ::List::c2_1
  long int *list_2_2_target,//!<[in] Similar to ::List::c2_2
  long int _irght,//!<[in] Similar to LargeList::irght
  long int _ilft,//!<[in] Similar to LargeList::ilft
  long int _ihfbit//!<[in] Similar to LargeList::ihfbit
) {
  long int list_1_j, list_1_off;
  long int ibit_tmp_1;
  long int bit;
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
    SgnBit(MP::myrank, &ipsgn); // Fermion sign
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
int mltply::Hubbard::X_Ajt(
  long int j,//!<[in] Index of initial wavefunction
  long int is1_spin,//!<[in] Bit mask
  long int *tmp_off,//!<[out] Index of final wavefunction
  long int *list_1_org,//!<[in] Similar to ::List::c1
  long int *list_2_1_target,//!<[in] Similar to ::List::c2_1
  long int *list_2_2_target,//!<[in] Similar to ::List::c2_2
  long int _irght,//!<[in] Similar to LargeList::irght
  long int _ilft,//!<[in] Similar to LargeList::ilft
  long int _ihfbit//!<[in] Similar to LargeList::ihfbit
) {
  long int list_1_j, list_1_off;
  long int ibit_tmp_1;
  long int bit;
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
    SgnBit(MP::myrank, &ipsgn); // Fermion sign
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
