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
@brief Functions for spin Hamiltonian

- mltplySpin() : Main routine of spin Hamiltonian (canonical)
  - mltplyHalfSpin() : 1/2 spin
  - mltplyGeneralSpin() : general spin
- mltplySpinGC() : Main routine of spin Hamiltonian (grandcanonical)
  - mltplyHalfSpinGC() : 1/2 spin
  - mltplyGeneralSpinGC() : general spin

Hub routines
<table>
  <tr>
    <td></td>
    <td>Get info</td>
    <td>Canonical</td>
    <td>Grandcanonical</td>
  </tr>
  <tr>
    <td>Exchange</td>
    <td>::child_exchange_spin_GetInfo</td>
    <td>::child_exchange_spin, ::child_exchange_spin_element</td>
    <td>::GC_child_exchange_spin, ::GC_child_exchange_spin_element</td>
  </tr>
  <tr>
    <td>Pair lift</td>
    <td>::child_pairlift_spin_GetInfo</td>
    <td></td>
    <td>::GC_child_pairlift_spin, ::GC_child_pairlift_spin_element</td>
  </tr>
  <tr>
    <td>General int.</td>
    <td>::child_general_int_spin_GetInfo</td>
    <td>::child_general_int_spin, ::child_general_int_spin_MPIsingle
    ::X_child_general_int_spin_MPIsingle, ::child_general_int_spin_MPIdouble,
    ::X_child_general_int_spin_MPIdouble</td>
    <td>::GC_child_general_int_spin, ::GC_child_general_int_spin_MPIsingle,
    ::GC_child_general_int_spin_MPIdouble</td>
  </tr>
  <tr>
    <td>General int for 1/2 spin</td>
    <td>::child_general_int_spin_GetInfo</td>
    <td>::child_general_int_spin, ::child_general_int_spin_MPIsingle
    ::X_child_general_int_spin_MPIsingle, ::child_general_int_spin_MPIdouble,
    ::X_child_general_int_spin_MPIdouble</td>
    <td>::GC_child_general_int_spin, ::GC_child_general_int_spin_MPIsingle,
    ::GC_child_general_int_spin_MPIdouble</td>
  </tr>
  <tr>
    <td>General int for general spin</td>
    <td></td>
    <td>::child_general_int_GeneralSpin_MPIsingle,
    ::child_general_int_GeneralSpin_MPIdouble</td>
    <td>::GC_child_general_int_GeneralSpin_MPIsingle,
    ::GC_child_general_int_GeneralSpin_MPIdouble</td>
  </tr>
</table>

General on-site term
<table>
  <tr>
    <td></td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td>General spin</td>
    <td>General spin</td>
  </tr>
  <tr>
    <td></td>
    <td>Canonical</td>
    <td>Canonical</td>
    <td>Grand canonical</td>
    <td>Grand canonical</td>
    <td>Canonical</td>
    <td>Grand canonical</td>
  </tr>
  <tr>
    <td></td>
    <td>In process</td>
    <td>Across process</td>
    <td>In process</td>
    <td>Across process</td>
    <td>Across process</td>
    <td>Across process</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s}@f$</td>
    <td>::X_Spin_CisAis</td>
    <td>::X_child_CisAis_spin_MPIdouble</td>
    <td>::X_SpinGC_CisAis</td>
    <td>::X_GC_child_CisAis_spin_MPIdouble</td>
    <td>::X_child_CisAis_GeneralSpin_MPIdouble</td>
    <td>::X_GC_child_CisAis_GeneralSpin_MPIdouble</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i t}@f$</td>
    <td>::X_Spin_CisAit</td>
    <td>::X_child_CisAit_spin_MPIdouble</td>
    <td>::X_SpinGC_CisAit</td>
    <td>::X_GC_child_CisAit_spin_MPIdouble</td>
    <td>::X_child_CisAit_GeneralSpin_MPIdouble</td>
    <td>::X_GC_child_CisAit_GeneralSpin_MPIdouble</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s} c_{i s}^\dagger c_{i s}@f$</td>
    <td>::child_CisAisCisAis_spin_element</td>
    <td></td>
    <td>::GC_child_CisAisCisAis_spin_element</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s} c_{i t}^\dagger c_{i u}@f$</td>
    <td></td>
    <td></td>
    <td>::GC_child_CisAisCitAiu_spin_element</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i t} c_{i u}^\dagger c_{i u}@f$</td>
    <td></td>
    <td></td>
    <td>::GC_child_CisAitCiuAiu_spin_element</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i t} c_{i u}^\dagger c_{i v}@f$</td>
    <td></td>
    <td></td>
    <td>::GC_child_CisAitCiuAiv_spin_element</td>
    <td>::GC_child_CisAitCiuAiv_spin_MPIsingle, ::X_GC_child_CisAitCiuAiv_spin_MPIsingle,
    ::GC_child_CisAitCiuAiv_spin_MPIdouble, ::X_GC_child_CisAitCiuAiv_spin_MPIdouble</td>
    <td></td>
    <td></td>
  </tr>
</table>
*/
#include "bitcalc.hpp"
#include "common/setmemory.hpp"
#include "mltplyCommon.hpp"
#include "mltplySpin.hpp"
#include "CalcTime.hpp"
#include "mltplySpinCore.hpp"
#include "mltplyHubbardCore.hpp"
#include "mltplyMPISpin.hpp"
#include "mltplyMPISpinCore.hpp"
#include "global.hpp"
/**
@brief Driver function for Spin hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplySpin(
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  int iret=0;
  if (Def::iFlgGeneralSpin == FALSE)
    iret = mltplyHalfSpin(nstate, tmp_v0, tmp_v1);
  else
    iret = mltplyGeneralSpin(nstate, tmp_v0, tmp_v1);
  return iret;
}/*int mltplySpin*/
/**
@brief Driver function for Spin 1/2 hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyHalfSpin(
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int i;
  long int isite1, isite2, sigma1, sigma2;
  long int sigma3, sigma4;

  /*[s] For InterAll */
  std::complex<double> tmp_V;
  /*[e] For InterAll */
  int ihermite=0;
  int idx=0;

  StartTimer(400);
  /**
  Transfer absorbed in Diagonal term.
  InterAll
  */
  StartTimer(410);
  for (i = 0; i < Def::NInterAll_OffDiagonal; i+=2) {
    if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite &&
        Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(411);
      child_general_int_spin_MPIdouble(i, nstate, tmp_v0, tmp_v1);
      StopTimer(411);
    }
    else if (Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(412);
      child_general_int_spin_MPIsingle(i, nstate, tmp_v0, tmp_v1);
      StopTimer(412);
    }
    else if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite) {
      StartTimer(413);
      child_general_int_spin_MPIsingle(i + 1, nstate, tmp_v0, tmp_v1);
      StopTimer(413);
    }
    else {
      StartTimer(414);
      for (ihermite = 0; ihermite<2; ihermite++) {
        idx = i + ihermite;
        isite1 = Def::InterAll_OffDiagonal[idx][0] + 1;
        isite2 = Def::InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = Def::InterAll_OffDiagonal[idx][1];
        sigma2 = Def::InterAll_OffDiagonal[idx][3];
        sigma3 = Def::InterAll_OffDiagonal[idx][5];
        sigma4 = Def::InterAll_OffDiagonal[idx][7];
        tmp_V = Def::ParaInterAll_OffDiagonal[idx];
        child_general_int_spin_GetInfo(isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
        child_general_int_spin(nstate, tmp_v0, tmp_v1);
      }/*for (ihermite = 0; ihermite<2; ihermite++)*/
      StopTimer(414);
    }
  }/*for (i = 0; i < Def::NInterAll_OffDiagonal; i+=2)*/
  StopTimer(410);
  /**
  Exchange 
  */
  StartTimer(420);   
  for (i = 0; i < Def::NExchangeCoupling; i++) {
    sigma1=0; sigma2=1;
    if (Def::ExchangeCoupling[i][0] + 1 > Def::Nsite &&
        Def::ExchangeCoupling[i][1] + 1 > Def::Nsite) {
      StartTimer(421);
      X_child_general_int_spin_MPIdouble(
        Def::ExchangeCoupling[i][0], sigma1, sigma2, 
        Def::ExchangeCoupling[i][1], sigma2, sigma1, 
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(421);
    }
    else if (Def::ExchangeCoupling[i][1] + 1 > Def::Nsite) {
      StartTimer(422);
      X_child_general_int_spin_MPIsingle(
        Def::ExchangeCoupling[i][0], sigma1, sigma2, 
        Def::ExchangeCoupling[i][1], sigma2, sigma1,
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(422);
    }
    else if (Def::ExchangeCoupling[i][0] + 1 > Def::Nsite) {
      StartTimer(423);
      X_child_general_int_spin_MPIsingle(
        Def::ExchangeCoupling[i][1], sigma2, sigma1, 
        Def::ExchangeCoupling[i][0], sigma1, sigma2, 
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(423);
    }
    else {
      StartTimer(424);
      child_exchange_spin_GetInfo(i);
      child_exchange_spin(nstate, tmp_v0, tmp_v1);
      StopTimer(424);
    }
  }/*for (i = 0; i < Def::NExchangeCoupling; i += 2)*/
  StopTimer(420);

  StopTimer(400);
  return 0;
}/*int mltplyHalfSpin*/
/**
@brief Driver function for General Spin hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyGeneralSpin(
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
){
  long int j;
  long int i;
  long int off = 0;
  long int tmp_off = 0;
  long int tmp_off2 = 0;
  long int ihfbit=0;
  long int isite1, isite2, sigma1, sigma2;
  long int sigma3, sigma4;

  long int tmp_sgn;
  /*[s] For InterAll */
  std::complex<double> tmp_V;
  int one = 1;
  /*[e] For InterAll */

  long int i_max;
  i_max = Check::idim_max;
  int ihermite=0;
  int idx=0;

  StartTimer(400);
  /**
  Transfer absorbed in Diagonal term.
  InterAll
  */
  StartTimer(410);
  ihfbit =Check::sdim;
  for (i = 0; i < Def::NInterAll_OffDiagonal; i += 2) {
    if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite &&
        Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(411);
      child_general_int_GeneralSpin_MPIdouble(i, nstate, tmp_v0, tmp_v1);
      StopTimer(411);
    }
    else if (Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(412);
      child_general_int_GeneralSpin_MPIsingle(i, nstate, tmp_v0, tmp_v1);
      StopTimer(412);
    }
    else if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite) {
      StartTimer(412);
      child_general_int_GeneralSpin_MPIsingle(i + 1, nstate, tmp_v0, tmp_v1);
      StopTimer(412);
    }
    else {
      StartTimer(413);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        isite1 = Def::InterAll_OffDiagonal[idx][0] + 1;
        isite2 = Def::InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = Def::InterAll_OffDiagonal[idx][1];
        sigma2 = Def::InterAll_OffDiagonal[idx][3];
        sigma3 = Def::InterAll_OffDiagonal[idx][5];
        sigma4 = Def::InterAll_OffDiagonal[idx][7];
        tmp_V = Def::ParaInterAll_OffDiagonal[idx];
#pragma omp parallel for default(none) private(j,tmp_sgn,off,tmp_off,tmp_off2) \
shared(i_max,isite1,isite2,sigma1,sigma2,sigma3,sigma4,tmp_V,ihfbit, tmp_v0,tmp_v1,list_1,list_2_1,list_2_2, \
one,nstate, Def::SiteToBit, Def::Tpow)
        for (j = 1; j <= i_max; j++) {
          tmp_sgn = GetOffCompGeneralSpin(list_1[j], isite2, sigma4, sigma3, &tmp_off, Def::SiteToBit, Def::Tpow);
          if (tmp_sgn == TRUE) {
            tmp_sgn = GetOffCompGeneralSpin(tmp_off, isite1, sigma2, sigma1, &tmp_off2, Def::SiteToBit, Def::Tpow);
            if (tmp_sgn == TRUE) {
              ConvertToList1GeneralSpin(tmp_off2, ihfbit, &off);
              zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[off][0], &one);
            }
          }/*if (tmp_sgn == TRUE)*/
        }/*for (j = 1; j <= i_max; j++)*/
      }/*for (ihermite = 0; ihermite < 2; ihermite++)*/
      StopTimer(413);
    }
  }/*for (i = 0; i < Def::NInterAll_OffDiagonal; i += 2)*/
  StopTimer(410);
  StopTimer(400);
  return 0;  
}/*int mltplyGeneralSpin*/
/**
@brief Driver function for Spin hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplySpinGC(
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  int iret=0;
  if (Def::iFlgGeneralSpin == FALSE) 
    iret = mltplyHalfSpinGC(nstate, tmp_v0, tmp_v1);
  else 
    iret = mltplyGeneralSpinGC(nstate, tmp_v0, tmp_v1);

  if(iret != 0) return iret;
  
  return iret;
}/*int mltplySpinGC*/
/**
@brief Driver function for Spin 1/2 Hamiltonian (grandcanonical)
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyHalfSpinGC(
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i;
  long int off = 0;
  long int is1_spin = 0;
  long int isite1, isite2, sigma1, sigma2;
  long int sigma3, sigma4;
  std::complex<double> tmp_trans;
  long int tmp_sgn;
  /*[s] For InterAll */
  std::complex<double> tmp_V;
  int one = 1;
  /*[e] For InterAll */

  long int i_max;
  i_max = Check::idim_max;

  int ihermite=0;
  int idx=0;

  StartTimer(500);

  StartTimer(510);
  for (i = 0; i < Def::EDNTransfer; i+=2 ) {
    if(Def::EDGeneralTransfer[i][0]+1 > Def::Nsite){
      if(Def::EDGeneralTransfer[i][1]==Def::EDGeneralTransfer[i][3]){
        fprintf(stderr, "Transverse_OffDiagonal component is illegal.\n");
      }
      else{
        StartTimer(511);
        X_GC_child_CisAit_spin_MPIdouble(
          Def::EDGeneralTransfer[i][0], Def::EDGeneralTransfer[i][1], 
          Def::EDGeneralTransfer[i][3], -Def::EDParaGeneralTransfer[i], 
          nstate, tmp_v0, tmp_v1);
        StopTimer(511);
      }
    }/*if(Def::EDGeneralTransfer[i][0]+1 > Def::Nsite)*/
    else{
      StartTimer(512);
      for(ihermite=0; ihermite<2; ihermite++){
        idx=i+ihermite;
        isite1 = Def::EDGeneralTransfer[idx][0] + 1;
        isite2 = Def::EDGeneralTransfer[idx][2] + 1;
        sigma1 = Def::EDGeneralTransfer[idx][1];
        sigma2 = Def::EDGeneralTransfer[idx][3];
        tmp_trans = -Def::EDParaGeneralTransfer[idx];
        if (child_general_hopp_GetInfo(isite1, isite2, sigma1, sigma2) != 0) {
          return -1;
        }
       
        if(sigma1==sigma2){
          fprintf(stderr, "Transverse_OffDiagonal component is illegal.\n");
        }
        else{
          // longitudinal magnetic field (considerd in diagonalcalc.cpp)
          // transverse magnetic field
          is1_spin = Def::Tpow[isite1 - 1];
#pragma omp parallel for default(none) private(j, tmp_sgn) \
shared(i_max, is1_spin, sigma2, off, tmp_trans,tmp_v0, tmp_v1,one,nstate)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = X_SpinGC_CisAit(j, is1_spin, sigma2, &off);
            if(tmp_sgn !=0){
              zaxpy_(&nstate, &tmp_trans, &tmp_v1[j][0], &one, &tmp_v0[off + 1][0], &one);
            }/*if(tmp_sgn !=0)*/
          }/*for (j = 1; j <= i_max; j++)*/
        }//sigma1 != sigma2
      }/*for(ihermite=0; ihermite<2; ihermite++)*/
      StopTimer(512);
    }
  }/*for (i = 0; i < Def::EDNTransfer; i+=2 )*/
  StopTimer(510);
  /**
  InterAll 
  */
  StartTimer(520);
  for (i = 0; i < Def::NInterAll_OffDiagonal; i += 2) {
    if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite &&
        Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(521);
      GC_child_general_int_spin_MPIdouble(i, nstate, tmp_v0, tmp_v1);
      StopTimer(521);
    }
    else if (Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(522);
      GC_child_general_int_spin_MPIsingle(i, nstate, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite) {
      StartTimer(522);
      GC_child_general_int_spin_MPIsingle(i + 1, nstate, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else {
      StartTimer(523);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        isite1 = Def::InterAll_OffDiagonal[idx][0] + 1;
        isite2 = Def::InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = Def::InterAll_OffDiagonal[idx][1];
        sigma2 = Def::InterAll_OffDiagonal[idx][3];
        sigma3 = Def::InterAll_OffDiagonal[idx][5];
        sigma4 = Def::InterAll_OffDiagonal[idx][7];
        tmp_V = Def::ParaInterAll_OffDiagonal[idx];
        child_general_int_spin_GetInfo(isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
        GC_child_general_int_spin(nstate, tmp_v0, tmp_v1);
      }
      StopTimer(523);
    }
  }/*for (i = 0; i < Def::NInterAll_OffDiagonal; i += 2)*/
  StopTimer(520);
  /**
  Exchange
  */
  StartTimer(530);
  for (i = 0; i < Def::NExchangeCoupling; i++) {
    sigma1=0; sigma2=1;
    if (Def::ExchangeCoupling[i][0] + 1 > Def::Nsite &&
        Def::ExchangeCoupling[i][1] + 1 > Def::Nsite){
      StartTimer(531);
      X_GC_child_CisAitCiuAiv_spin_MPIdouble(
        Def::ExchangeCoupling[i][0], sigma1, sigma2, 
        Def::ExchangeCoupling[i][1], sigma2, sigma1, 
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(531);
    }
    else if (Def::ExchangeCoupling[i][1] + 1 > Def::Nsite) {
      StartTimer(532);
      X_GC_child_CisAitCiuAiv_spin_MPIsingle(
        Def::ExchangeCoupling[i][0], sigma1, sigma2,
        Def::ExchangeCoupling[i][1], sigma2, sigma1,
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(532);
    }
    else if (Def::ExchangeCoupling[i][0] + 1 > Def::Nsite) {
      StartTimer(532);
      X_GC_child_CisAitCiuAiv_spin_MPIsingle(
        Def::ExchangeCoupling[i][1], sigma2, sigma1,
        Def::ExchangeCoupling[i][0], sigma1, sigma2,
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(532);
    }
    else {
      StartTimer(533);
      child_exchange_spin_GetInfo(i);
      GC_child_exchange_spin(nstate, tmp_v0, tmp_v1);
      StopTimer(533);
    }
  }/* for (i = 0; i < Def::NExchangeCoupling; i ++) */
  StopTimer(530);
  /**
  PairLift
  */
  StartTimer(540);
  for (i = 0; i < Def::NPairLiftCoupling; i++) {
    sigma1 =0; sigma2=1;
    if (Def::PairLiftCoupling[i][0] + 1 > Def::Nsite &&
        Def::PairLiftCoupling[i][1] + 1 > Def::Nsite) {
      StartTimer(541);
      X_GC_child_CisAitCiuAiv_spin_MPIdouble(
        Def::PairLiftCoupling[i][0], sigma1, sigma2, 
        Def::PairLiftCoupling[i][1], sigma1, sigma2,
        Def::ParaPairLiftCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(541);
    }
    else if (Def::PairLiftCoupling[i][1] + 1 > Def::Nsite) {
      StartTimer(542);
      X_GC_child_CisAitCiuAiv_spin_MPIsingle(
        Def::PairLiftCoupling[i][0], sigma1, sigma2, 
        Def::PairLiftCoupling[i][1], sigma1, sigma2, 
        Def::ParaPairLiftCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(542);
    }
    else if (Def::PairLiftCoupling[i][0] + 1 > Def::Nsite) {
      StartTimer(542);
      X_GC_child_CisAitCiuAiv_spin_MPIsingle(
        Def::PairLiftCoupling[i][1], sigma1, sigma2,
        Def::PairLiftCoupling[i][0], sigma1, sigma2,
        Def::ParaPairLiftCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(542);
    }
    else {
      StartTimer(543);
      child_pairlift_spin_GetInfo(i);
      GC_child_pairlift_spin(nstate, tmp_v0, tmp_v1);
      StopTimer(543);
    }
  }/*for (i = 0; i < Def::NPairLiftCoupling; i += 2)*/
  StopTimer(540);

  StopTimer(500);
  return 0;
}/*int mltplyHalfSpinGC*/
/**
@brief Driver function for General Spin hamiltonian (grandcanonical)
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyGeneralSpinGC(
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i;
  long int off = 0;
  long int tmp_off = 0;
  long int isite1, isite2, sigma1, sigma2;
  long int sigma3, sigma4;
  std::complex<double> tmp_trans;
  long int tmp_sgn;
  double num1 = 0;
  /*[s] For InterAll */
  std::complex<double> tmp_V;
  int one = 1;
  /*[e] For InterAll */

  long int i_max;
  i_max = Check::idim_max;

  int ihermite=0;
  int idx=0;

  StartTimer(500);

  StartTimer(510);
  for (i = 0; i < Def::EDNTransfer; i += 2) {
    isite1 = Def::EDGeneralTransfer[i][0] + 1;
    isite2 = Def::EDGeneralTransfer[i][2] + 1;
    sigma1 = Def::EDGeneralTransfer[i][1];
    sigma2 = Def::EDGeneralTransfer[i][3];
    tmp_trans = -Def::EDParaGeneralTransfer[idx];
    if (isite1 == isite2) {
      if (sigma1 != sigma2) {
        if (isite1 > Def::Nsite) {
          X_GC_child_CisAit_GeneralSpin_MPIdouble(
            isite1 - 1, sigma1, sigma2, tmp_trans, nstate, tmp_v0, tmp_v1);
        }/*if (isite1 > Def::Nsite)*/
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            isite1 = Def::EDGeneralTransfer[idx][0] + 1;
            isite2 = Def::EDGeneralTransfer[idx][2] + 1;
            sigma1 = Def::EDGeneralTransfer[idx][1];
            sigma2 = Def::EDGeneralTransfer[idx][3];
            tmp_trans = -Def::EDParaGeneralTransfer[idx];
    
            // transverse magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn, num1) \
shared(i_max, isite1, sigma1, sigma2, off, tmp_trans,tmp_v0, \
tmp_v1,one,nstate, Def::SiteToBit, Def::Tpow)
            for (j = 1; j <= i_max; j++) {
              num1 = GetOffCompGeneralSpin(
                j - 1, isite1, sigma2, sigma1, &off, Def::SiteToBit, Def::Tpow);
              if (num1 != 0) { // for multply
                zaxpy_(&nstate, &tmp_trans, &tmp_v1[j][0], &one, &tmp_v0[off + 1][0], &one);
              }/*if (num1 != 0)*/
            }/*for (j = 1; j <= i_max; j++)*/
          }/*for (ihermite = 0; ihermite<2; ihermite++)*/
        }
      }// sigma1 != sigma2          
      else{ // sigma1 = sigma2
        fprintf(stderr, "Error: Transverse_Diagonal component must be absorbed !");
      }
    }//isite1 = isite2
    else { // isite1 != isite2
      // hopping is not allowed in localized spin system
      return -1;
    }
  }/*for (i = 0; i < Def::EDNTransfer; i += 2)*/
  StopTimer(510);
  /**
  InterAll
  */
  StartTimer(520);
  for (i = 0; i< Def::NInterAll_OffDiagonal; i += 2) {
    if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite &&
        Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(521);
      GC_child_general_int_GeneralSpin_MPIdouble(i, nstate, tmp_v0, tmp_v1);
      StopTimer(521);
    }
    else if (Def::InterAll_OffDiagonal[i][4] + 1 > Def::Nsite) {
      StartTimer(522);
      GC_child_general_int_GeneralSpin_MPIsingle(i, nstate, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else if (Def::InterAll_OffDiagonal[i][0] + 1 > Def::Nsite) {
      StartTimer(522);
      GC_child_general_int_GeneralSpin_MPIsingle(i + 1, nstate, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else {
      StartTimer(523);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        isite1 = Def::InterAll_OffDiagonal[idx][0] + 1;
        isite2 = Def::InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = Def::InterAll_OffDiagonal[idx][1];
        sigma2 = Def::InterAll_OffDiagonal[idx][3];
        sigma3 = Def::InterAll_OffDiagonal[idx][5];
        sigma4 = Def::InterAll_OffDiagonal[idx][7];
        tmp_V = Def::ParaInterAll_OffDiagonal[idx];

        if (sigma1 == sigma2) {
          if (sigma3 == sigma4) {
            fprintf(stderr, "InterAll_OffDiagonal component is illegal.\n");
            return -1;
          }/*if (sigma3 == sigma4)*/
          else {
            //sigma3=sigma4 term is considerd as a diagonal term.
#pragma omp parallel for default(none) private(j, tmp_sgn, off) \
shared(i_max, isite1, isite2, sigma1, sigma3, sigma4, tmp_V,tmp_v0, tmp_v1, \
one,nstate, Def::SiteToBit, Def::Tpow)
            for (j = 1; j <= i_max; j++) {
              tmp_sgn = GetOffCompGeneralSpin(
                j - 1, isite2, sigma4, sigma3, &off, Def::SiteToBit, Def::Tpow);
              if (tmp_sgn == TRUE) {
                tmp_sgn = BitCheckGeneral(off, isite1, sigma1, Def::SiteToBit, Def::Tpow);
                if (tmp_sgn == TRUE) {
                  zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[off + 1][0], &one);
               }/*if (tmp_sgn == TRUE)*/
              }/*if (tmp_sgn == TRUE)*/
            }/*for (j = 1; j <= i_max; j++)*/
          }
        }/*if (sigma1 == sigma2)*/
        else if (sigma3 == sigma4) {
          //sigma1=sigma2 term is considerd as a diagonal term.
#pragma omp parallel for default(none) private(j, tmp_sgn, off, tmp_off)                                \
shared(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V, \
tmp_v0, tmp_v1,one,nstate, Def::SiteToBit, Def::Tpow)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = BitCheckGeneral(j - 1, isite2, sigma3, Def::SiteToBit, Def::Tpow);
            if (tmp_sgn == TRUE) {
              tmp_sgn = GetOffCompGeneralSpin(
                j - 1, isite1, sigma2, sigma1, &off, Def::SiteToBit, Def::Tpow);
              if (tmp_sgn == TRUE) {
                zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[off + 1][0], &one);
              }/*if (tmp_sgn == TRUE)*/
            }/*if (tmp_sgn == TRUE)*/
          }/*for (j = 1; j <= i_max; j++)*/
        }/*else if (sigma3 == sigma4)*/
        else {
#pragma omp parallel for default(none) private(j, tmp_sgn, off, tmp_off) \
shared(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V, \
tmp_v0, tmp_v1,one,nstate, Def::SiteToBit, Def::Tpow)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = GetOffCompGeneralSpin(
              j - 1, isite2, sigma4, sigma3, &tmp_off, Def::SiteToBit, Def::Tpow);
            if (tmp_sgn == TRUE) {
              tmp_sgn = GetOffCompGeneralSpin(
                tmp_off, isite1, sigma2, sigma1, &off, Def::SiteToBit, Def::Tpow);
              if (tmp_sgn == TRUE) {
                zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[off + 1][0], &one);
              }/*if (tmp_sgn == TRUE)*/
            }/*if (tmp_sgn == TRUE)*/
          }/*for (j = 1; j <= i_max; j++)*/
        }
      }
      StopTimer(523);
    }
  }/*for (i = 0; i< Def::NInterAll_OffDiagonal; i += 2)*/
  StopTimer(520);

  StopTimer(500);
  return 0;
}/*int mltplyGeneralSpinGC*/

/******************************************************************************/
//[s] child functions
/******************************************************************************/

/**
@brief Compute exchange term of spin Hamiltonian (canonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_exchange_spin(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i_max = Large::i_max;
  long int off = 0;

#pragma omp parallel for default(none) private(j) \
shared(i_max, off, tmp_v0, tmp_v1,nstate)
  for (j = 1; j <= i_max; j++) 
    child_exchange_spin_element(j, nstate, tmp_v0, tmp_v1, &off);
}/*std::complex<double> child_exchange_spin*/
/**
@brief Compute exchange term of spin Hamiltonian (grandcanonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_exchange_spin(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i_max = Large::i_max;
  long int off = 0;

#pragma omp parallel for default(none) private(j) \
shared(i_max, off,tmp_v0, tmp_v1,nstate)
  for (j = 1; j <= i_max; j++)
    GC_child_exchange_spin_element(j, nstate, tmp_v0, tmp_v1, &off);
}/*std::complex<double> GC_child_exchange_spin*/
/**
@brief Compute pair-lift term of spin Hamiltonian (grandcanonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_pairlift_spin(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i_max = Large::i_max;
  long int off = 0;

#pragma omp parallel for default(none) private(j) \
shared(i_max, off, tmp_v0, tmp_v1,nstate)
  for (j = 1; j <= i_max; j++) 
    GC_child_pairlift_spin_element(j, nstate, tmp_v0, tmp_v1, &off);
}/*std::complex<double> GC_child_pairlift_spin*/
/**
@brief Compute Inter-All term of spin Hamiltonian (canonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void child_general_int_spin(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  std::complex<double> tmp_V, dmv;
  long int j, i_max;
  long int org_sigma2, org_sigma4;
  long int isA_up, isB_up;
  long int tmp_off = 0;
  int tmp_sgn;
  int one = 1;

  i_max = Large::i_max;
  org_sigma2 = Large::is2_spin;
  org_sigma4 = Large::is4_spin;
  tmp_V = Large::tmp_V;
  isA_up = Large::is1_up;
  isB_up = Large::is2_up;

#pragma omp parallel for default(none) private(j, tmp_sgn, dmv) \
shared(i_max,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,tmp_v1, tmp_v0,one,nstate)
  for (j = 1; j <= i_max; j++) {
    tmp_sgn = X_child_exchange_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, &tmp_off);
    if (tmp_sgn != 0) {
      dmv = (std::complex<double>)tmp_sgn * tmp_V;
      zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[tmp_off][0], &one);
    }/*if (tmp_sgn != 0)*/
  }/*for (j = 1; j <= i_max; j++)*/
}/*std::complex<double> child_general_int_spin*/
/**
@brief Compute Inter-All term of spin Hamiltonian (grandcanonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_child_general_int_spin(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  std::complex<double> tmp_V;
  long int j, i_max;
  long int org_isite1, org_isite2;
  long int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long int isA_up, isB_up;
  long int tmp_off = 0;

  i_max = Large::i_max;
  org_isite1 = Large::isite1;
  org_isite2 = Large::isite2;
  org_sigma1 = Large::is1_spin;
  org_sigma2 = Large::is2_spin;
  org_sigma3 = Large::is3_spin;
  org_sigma4 = Large::is4_spin;
  tmp_V = Large::tmp_V;
  isA_up = Def::Tpow[org_isite1 - 1];
  isB_up = Def::Tpow[org_isite2 - 1];

#pragma omp parallel default(none) private(j) \
shared(tmp_v0,tmp_v1,nstate,i_max,isA_up,isB_up,org_sigma1,org_sigma2,org_sigma3,org_sigma4,tmp_off,tmp_V) 
  {
    if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp for
      for (j = 1; j <= i_max; j++)
        GC_child_CisAisCisAis_spin_element(
          j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, nstate, tmp_v0, tmp_v1);
    }
    else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        GC_child_CisAisCitAiu_spin_element(
          j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off);
    }
    else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        GC_child_CisAitCiuAiu_spin_element(
          j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off);
    }
    else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        GC_child_CisAitCiuAiv_spin_element(
          j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off);
    }
  }/*End of parallel region*/
}/*std::complex<double> GC_child_general_int_spin*/
/******************************************************************************/
//[e] child functions
/******************************************************************************/
