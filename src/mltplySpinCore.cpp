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
@brief Functions for spin Hamiltonian (Core)
*/

#include <bitcalc.hpp>
#include "xsetmem.hpp"
#include "wrapperMPI.hpp"
#include "mltplyCommon.hpp"
#include "mltplySpinCore.hpp"
#include "global.hpp"

/******************************************************************************/
//[s] GetInfo functions
/******************************************************************************/

/**
@brief Set parameters for the bit operation of spin-exchange term
@return Always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::exchange_GetInfo(
  int iExchange//!<[in] Counter of exchange interaction
) {
  int isite1 = Def::ExchangeCoupling[iExchange][0];
  int isite2 = Def::ExchangeCoupling[iExchange][1];
  /**
   Set the exchange coupling constant (LargeList::tmp_J)
  */
  Large::tmp_J = Def::ParaExchangeCoupling[iExchange];
  /**
  Set the bit mask for computing spin state of both site
  (LargeList::is1_up, LargeList::is2_up)
  */
  Large::is1_up = Def::Tpow[isite1];
  Large::is2_up = Def::Tpow[isite2];
  /**
  Set the bit mask for exchange 2 spins (LargeList::isA_spin)
  */
  Large::isA_spin = Large::is1_up + Large::is2_up;
  return 0;
}/*int child_exchange_spin_GetInfo*/
/**
@brief Set parameters for the bit operation of spin-pairlift term
@return Always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::pairlift_GetInfo(
  int iPairLift
) {
  int isite1 = Def::PairLiftCoupling[iPairLift][0];
  int isite2 = Def::PairLiftCoupling[iPairLift][1];
  /**
  Set the pairlift coupling constant (LargeList::tmp_J)
  */
  Large::tmp_J = Def::ParaPairLiftCoupling[iPairLift];
  /**
  Set the bit mask for computing spin state of both site
  (LargeList::is1_up, LargeList::is2_up)
  */
  Large::is1_up = Def::Tpow[isite1];
  Large::is2_up = Def::Tpow[isite2];
  /**
  Set the bit mask for exchange 2 spins (LargeList::isA_spin)
  */
  Large::isA_spin = Large::is1_up + Large::is2_up;
  return 0;
}/*int child_pairlift_spin_GetInfo*/
/**
@brief Set parameters for the bit operation of spin-general interaction term
@return Always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::general_int_GetInfo(
  long int isite1,//!<[in] Site 1
  long int isite2,//!<[in] Site 2
  long int sigma1,//!<[in] Sigma 1, final state of site 1
  long int sigma2,//!<[in] Sigma 3, initial state of site 1
  long int sigma3,//!<[in] Sigma 3, final state of site 2
  long int sigma4,//!<[in] Sigma 4, initial state of site 2
  std::complex<double> tmp_V//!<[in] General interaction constatnt
) {
  /**
  Set the pairlift coupling constant (LargeList::tmp_J)
  */
  Large::tmp_V = tmp_V;
  Large::isite1 = isite1;
  Large::isite2 = isite2;
  /**
  Set the bit mask for computing spin state of both site
  (LargeList::is1_up, LargeList::is2_up)
  */
  Large::is1_up = Def::Tpow[isite1];
  Large::is2_up = Def::Tpow[isite2];
  /**
  Set the bit mask for general interaction 
  (LargeList::is1_spin, LargeList::is2_spin, LargeList::is3_spin, LargeList::is4_spin)
  */
  Large::is1_spin = sigma1;
  Large::is2_spin = sigma2;
  Large::is3_spin = sigma3;
  Large::is4_spin = sigma4;
  return 0;
}/*int child_general_int_spin_GetInfo*/

/******************************************************************************/
//[e] GetInfo functions
/******************************************************************************/

/******************************************************************************/
//[s] core routines
/******************************************************************************/

/**
@brief Compute index of final wavefunction by @f$c_{is}^\dagger c_{it}@f$ term.
@return 1 if bit-mask is1_spin is sigma2, 0 for the other
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::C::Half::X_CisAit(
  long int j,//!<[in] Index of initial wavefunction
  long int is1_spin,//!<[in] Bit mask for computing spin state
  long int sigma2,//!<[in] Spin state at site 2
  long int *tmp_off,//!<[out] Index of final wavefunction
  long int *list_1, 
  long int *list_2_1,
  long int *list_2_2
) {
  long int list_1_j;
  long int off;
  list_1_j = list_1[j];
  if (mltply::Spin::GC::Half::X_CisAit(list_1_j, is1_spin, sigma2, &off) != 0) {
    GetOffComp(list_2_1, list_2_2, off, Large::irght, Large::ilft, Large::ihfbit, tmp_off);
    return 1;
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int X_Spin_CisAit*/
/**
@brief Compute the spin state with bit mask is1_spin
@return 1 if the spin state with bit mask is1_spin is sigma1, 0 for the other.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::C::Half::X_CisAis(
  long int j,//!<[in] Index of wavefunction
  long int is1_spin,//!<[in] Bit mask
  long int sigma1,//!<[in] Target spin state
  long int *list_1
) {
  int A_ibit_tmp;
  // off = j
  A_ibit_tmp = ((list_1[j] & is1_spin) / is1_spin) ^ (1 - sigma1);
  return A_ibit_tmp;
}/*int X_Spin_CisAis*/
/**
@brief Compute the grandcanonical spin state with bit mask is1_spin
@return 1 if the spin state with bit mask is1_spin is sigma1, 0 for the other.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::GC::Half::X_CisAis(
  long int j,//!<[in] Index of wavefunction
  long int is1_spin,//!<[in] Bit mask
  long int sigma1//!<[in] Target spin state
) {
  int A_ibit_tmp;
  long int list_1_j;
  // off = j
  list_1_j = j;
  A_ibit_tmp = ((list_1_j & is1_spin) / is1_spin) ^ (1 - sigma1);
  return A_ibit_tmp;
}/*int X_SpinGC_CisAis*/
/**
@brief Compute index of final wavefunction by @f$c_{is}^\dagger c_{it}@f$ term
(grandcanonical).
@return 1 if bit-mask is1_spin is sigma2, 0 for the other
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::GC::Half::X_CisAit(
  long int j,//!<[in] Index of wavefunction
  long int is1_spin,//!<[in] Bit mask for computing spin state
  long int sigma2,//!<[in] Spin state at site 2
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  long int list_1_j, ibit_tmp_1;

  list_1_j = j;

  ibit_tmp_1 = list_1_j & is1_spin;
  if (ibit_tmp_1 == 0 && sigma2 == 0) {    // down -> up
    *tmp_off = list_1_j + is1_spin;
    return 1;
  }
  else if (ibit_tmp_1 != 0 && sigma2 == 1) { // up -> down
    *tmp_off = list_1_j - is1_spin;
    return 1;
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int X_SpinGC_CisAit*/

/******************************************************************************/
//[e] core routines
/******************************************************************************/

/**
@brief Compute index of final wavefunction associated to spin-exchange term
@return 1 if spin of site 1 is sigmaA and spin of site 2 is sigmaB. 0 if not. 
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Spin::C::Half::X_exchange_element(
  long int j,//!<[in] Index of initial wavefunction
  //!<[inout]
  long int isA_up,//!<[in] Bit mask for spin 1
  long int isB_up,//!<[in] Bit mask for spin 2
  long int sigmaA,//!<[in] Target of spin 1
  long int sigmaB,//!<[in] Target of spin 2
  long int *tmp_off,//!<[out] Index of final wavefunction
  long int* list_1, 
  long int* list_2_1,
  long int* list_2_2
) {
  long int iexchg, off;
  long int irght = Large::irght;
  long int ilft = Large::ilft;
  long int ihfbit = Large::ihfbit;
  long int ibit_tmp_A, ibit_tmp_B;

  ibit_tmp_A = ((list_1[j] & isA_up) / isA_up);
  ibit_tmp_B = ((list_1[j] & isB_up) / isB_up);
  if (ibit_tmp_A == sigmaA && ibit_tmp_B == sigmaB) {
    iexchg = list_1[j] ^ (isA_up + isB_up);
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
    *tmp_off = off;
    return 1;
  }
  else {
    *tmp_off = 0; // just tentative
    return 0;
  }
}/*int X_child_exchange_spin_element*/
/**
@brief Multiply Hamiltonian of exchange term of canonical spin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::C::Half::exchange_element(
  long int j,//!<[in] Index of initial wavefunction
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off,//!<[out] Index of final wavefunction
  long int* list_1,
  long int* list_2_1, 
  long int* list_2_2
) {
  long int off;
  long int iexchg;
  long int is_up = Large::isA_spin;
  long int irght = Large::irght;
  long int ilft = Large::ilft;
  long int ihfbit = Large::ihfbit;
  std::complex<double> tmp_J = Large::tmp_J;
  long int ibit_tmp;
  int one = 1;

  ibit_tmp = (list_1[j] & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    return;
  }
  else {
    iexchg = list_1[j] ^ is_up;
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
    *tmp_off = off;
    zaxpy_(&nstate, &tmp_J, &tmp_v1[j][0], &one, &tmp_v0[off][0], &one);
  }
}/*std::complex<double> child_exchange_spin_element*/
/**
@brief Multiply Hamiltonian of exchange term of grandcanonical spin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::GC::Half::exchange_element(
  long int j,//!<[in] Index of initial wavefunction
  int nstate, std::complex<double> **tmp_v0,//!<[out] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  //!<[inout]
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  long int is_up = Large::isA_spin;
  std::complex<double> tmp_J = Large::tmp_J;
  long int list_1_j, list_1_off;
  int one = 1;

  list_1_j = j;

  long int ibit_tmp;
  ibit_tmp = (list_1_j & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    return;
  }
  else {
    list_1_off = list_1_j ^ is_up;
    *tmp_off = list_1_off;
    zaxpy_(&nstate, &tmp_J, &tmp_v1[j][0], &one, &tmp_v0[list_1_off][0], &one);
  }
}/*std::complex<double> GC_child_exchange_spin_element*/
/**
@brief Multiply Hamiltonian of pairlift term of grandcanonical spin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::GC::Half::pairlift_element(
  long int j,//!<[in] Index of initial wavefunction
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  long int is_up = Large::isA_spin;
  std::complex<double> tmp_J = Large::tmp_J;
  int one = 1;
  long int list_1_off;
  long int list_1_j = j;
  long int ibit_tmp;
  //ibit_tmp = ((list_1_j & is1_up) / is1_up) ^ ((list_1_j & is2_up) / is2_up);
  ibit_tmp = (list_1_j & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    list_1_off = list_1_j ^ is_up; //Change: ++ -> -- or -- -> ++
    *tmp_off = list_1_off;
    zaxpy_(&nstate, &tmp_J, &tmp_v1[j][0], &one, &tmp_v0[list_1_off][0], &one);
  }
}/*std::complex<double> GC_child_pairlift_spin_element*/
//[s]Spin
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$ term of 
canonical spsin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::C::Half::CisAisCisAis_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isA_up,//!<[in] Bit mask for spin 1
  long int isB_up,//!<[in] Bit mask for spin 2
  long int org_sigma2,//!<[in] Target for spin 1
  long int org_sigma4,//!<[in] Target for spin 2
  std::complex<double> tmp_V,//!<[in] Coupling constatnt
  int nstate, 
  std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *list_1
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;

  tmp_sgn = mltply::Spin::C::Half::X_CisAis(j, isB_up, org_sigma4, list_1);
  tmp_sgn *= mltply::Spin::C::Half::X_CisAis(j, isA_up, org_sigma2, list_1);
  dmv = (std::complex<double>)tmp_sgn * tmp_V;
  zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
}/*std::complex<double> child_CisAisCisAis_spin_element*/

//[e]Spin

//[s]GC Spin
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$ term of 
grandcanonical spsin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::GC::Half::CisAisCisAis_element(
  long int j,//!<[in] Index of initial wavefunction
  long int isA_up,//!<[in] Bit mask for spin 1
  long int isB_up,//!<[in] Bit mask for spin 2
  long int org_sigma2,//!<[in] Target for spin 1
  long int org_sigma4,//!<[in] Target for spin 2
  std::complex<double> tmp_V,//!<[in] Coupling constatnt
  int nstate, std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[in] Wavefunction to be multiplied
) {
  int tmp_sgn;
  std::complex<double> dmv = 0;
  int one = 1;

  tmp_sgn = mltply::Spin::GC::Half::X_CisAis(j, isB_up, org_sigma4);
  tmp_sgn *= mltply::Spin::GC::Half::X_CisAis(j, isA_up, org_sigma2);
  if (tmp_sgn != 0) {
    dmv = (std::complex<double>)tmp_sgn * tmp_V;
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }
}/*std::complex<double> GC_child_CisAisCisAis_spin_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{it}^\dagger c_{iu}@f$ term of 
grandcanonical spsin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::GC::Half::CisAisCitAiu_element(
  long int j,//!<[in] Index of initial wavefunction
  long int org_sigma2,//!<[in] Target for spin 1
  long int org_sigma4,//!<[in] Target for spin 2
  long int isA_up,//!<[in] Bit mask for spin 1
  long int isB_up,//!<[in] Bit mask for spin 2
  std::complex<double> tmp_V,//!<[in] Coupling constatnt
  int nstate, 
  std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = mltply::Spin::GC::Half::X_CisAit(j, isB_up, org_sigma4, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= mltply::Spin::GC::Half::X_CisAis((*tmp_off), isA_up, org_sigma2);
    if (tmp_sgn != 0) {
      dmv = (std::complex<double>)tmp_sgn * tmp_V;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off][0], &one);
    }/*if (tmp_sgn != 0)*/
  }/*if (tmp_sgn != 0)*/
}/*std::complex<double> GC_child_CisAisCitAiu_spin_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{it} c_{iu}^\dagger c_{iu}@f$ term of 
grandcanonical spsin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::GC::Half::CisAitCiuAiu_element(
  long int j,//!<[in] Index of initial wavefunction
  long int org_sigma2,//!<[in] Target for spin 1
  long int org_sigma4,//!<[in] Target for spin 2
  long int isA_up,//!<[in] Bit mask for spin 1
  long int isB_up,//!<[in] Bit mask for spin 2
  std::complex<double> tmp_V,//!<[in] Coupling constatnt
  int nstate,
  std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = mltply::Spin::GC::Half::X_CisAis(j, isB_up, org_sigma4);
  if (tmp_sgn != 0) {
    tmp_sgn *= mltply::Spin::GC::Half::X_CisAit(j, isA_up, org_sigma2, tmp_off);
    if (tmp_sgn != 0) {
      dmv = (std::complex<double>)tmp_sgn * tmp_V;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off][0], &one);
    }/*if (tmp_sgn != 0)*/
  }/*if (tmp_sgn != 0)*/
}/*std::complex<double> GC_child_CisAitCiuAiu_spin_element*/
/**
@brief Compute @f$c_{is}^\dagger c_{it} c_{iu}^\dagger c_{iv}@f$ term of
grandcanonical spsin system
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Spin::GC::Half::CisAitCiuAiv_element(
  long int j,//!<[in] Index of initial wavefunction
  long int org_sigma2,//!<[in] Target for spin 1
  long int org_sigma4,//!<[in] Target for spin 2
  long int isA_up,//!<[in] Bit mask for spin 1
  long int isB_up,//!<[in] Bit mask for spin 2
  std::complex<double> tmp_V,//!<[in] Coupling constatnt
  int nstate,
  std::complex<double> **tmp_v0,//!<[in] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[in] Wavefunction to be multiplied
  long int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int tmp_off_1;
  std::complex<double> dmv;
  int one = 1;
  tmp_sgn = mltply::Spin::GC::Half::X_CisAit(j, isB_up, org_sigma4, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= mltply::Spin::GC::Half::X_CisAit(tmp_off_1, isA_up, org_sigma2, tmp_off_2);
    if (tmp_sgn != 0) {
      dmv = (std::complex<double>)tmp_sgn * tmp_V;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[*tmp_off_2][0], &one);
    }/*if (tmp_sgn != 0)*/
  }/*if (tmp_sgn != 0)*/
}/*std::complex<double> GC_child_CisAitCiuAiv_spin_element*/
//[e]GC Spin
