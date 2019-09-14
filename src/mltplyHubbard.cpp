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
@brief Function for Hubbard Hamitonian
<ul>
<li>mltplyHubbard(): Main routine of hubbard hamiltonian (canonical).</li>
<li>mltplyHubbardGC(): Main routine of hubbard hamiltonian (grandcanonical).</li>
</ul>
<table>
  <tr>
  <td></td>
    <td>Canonical</td>
    <td>Canonical</td>
    <td>Grand Canonical</td>
    <td>Grand Canonical</td>
  </tr>
  <tr>
    <td></td>
    <td>Intra Process</td>
    <td>Inter Process</td>
    <td>Intra Process</td>
    <td>Inter Process</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger@f$</td>
    <td>::Cis, ::X_Cis</td>
    <td>::Cis_MPI, ::mltply::Hubbard::C::X_Cis_MPI</td>
    <td>::GC_Cis, ::X_GC_Cis</td>
    <td>::GC_Cis_MPI, ::mltply::Hubbard::GC::X_Cis_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{j t}@f$</td>
    <td>::Ajt, ::X_Ajt</td>
    <td>::Ajt_MPI, ::mltply::Hubbard::C::X_Ajt_MPI</td>
    <td>::GC_Ajt, ::X_GC_Ajt</td>
    <td>::GC_Ajt_MPI, ::mltply::Hubbard::GC::X_Ajt_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s}@f$</td>
    <td>::CisAis, ::X_CisAis</td>
    <td>::child_CisAis_Hubbard_MPI, ::mltply::Hubbard::C::X_CisAis_MPI</td>
    <td>::GC_CisAis, ::X_CisAis</td>
    <td>::GC_child_CisAis_Hubbard_MPI, ::mltply::Hubbard::GC::X_CisAis_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{j t}@f$</td>
    <td>::CisAjt, ::X_CisAjt</td>
    <td>::child_CisAjt_MPIsingle, ::child_CisAjt_MPIdouble, 
    ::X_child_CisAjt_MPIsingle, ::mltply::Hubbard::C::X_CisAjt_MPIdouble</td>
    <td>::GC_CisAjt, ::X_GC_CisAjt</td>
    <td>::GC_child_CisAjt_Hubbard_MPI, ::mltply::Hubbard::GC::X_CisAjt_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s} c_{i s}^\dagger c_{i s}@f$</td>
    <td>::child_CisAisCisAis_element</td>
    <td>::X_child_CisAisCisAis_Hubbard_MPI</td>
    <td>::GC_child_CisAisCisAis_element</td>
    <td>::X_GC_child_CisAisCisAis_Hubbard_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s} c_{j t}^\dagger c_{j t}@f$</td>
    <td>::child_CisAisCjtAjt_element</td>
    <td>::mltply::Hubbard::C::X_CisAisCjtAjt_MPI</td>
    <td>::GC_child_CisAisCjtAjt_element</td>
    <td>::mltply::Hubbard::GC::X_CisAisCjtAjt_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s} c_{j t}^\dagger c_{k u}@f$</td>
    <td>::child_CisAisCjtAku_element</td>
    <td>::mltply::Hubbard::C::X_CisAisCjtAku_MPI</td>
    <td>::GC_child_CisAisCjtAku_element</td>
    <td>::mltply::Hubbard::GC::X_CisAisCjtAku_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{j t} c_{k u}^\dagger c_{k u}@f$</td>
    <td>::child_CisAjtCkuAku_element</td>
    <td>::mltply::Hubbard::C::X_CisAjtCkuAku_MPI</td>
    <td>::GC_child_CisAjtCkuAku_element</td>
    <td>::mltply::Hubbard::GC::X_CisAjtCkuAku_MPI</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{j t} c_{k u}^\dagger c_{l v}@f$</td>
    <td>::child_CisAjtCkuAlv_element</td>
    <td>::mltply::Hubbard::C::X_CisAjtCkuAlv_MPI</td>
    <td>::GC_child_CisAjtCkuAlv_element</td>
    <td>::mltply::Hubbard::GC::X_CisAjtCkuAlv_MPI</td>
  </tr>
</table>
Other
<table>
  <tr>
    <td></td>
    <td>Get info</td>
    <td>Canonical</td>
    <td>Grandcanonical</td>
  </tr>
  <tr>
    <td>Exchange</td>
    <td>::child_exchange_GetInfo</td>
    <td>::child_exchange, ::child_exchange_element</td>
    <td>::GC_child_exchange, ::GC_child_exchange_element</td>
  </tr>
  <tr>
    <td>Pair hop</td>
    <td>::child_pairhopp_GetInfo</td>
    <td>::child_pairhopp, ::child_pairhopp_element</td>
    <td>::GC_child_pairhopp, ::GC_child_pairhopp_element</td>
  </tr>
  <tr>
    <td>General int.</td>
    <td>::child_general_int_GetInfo</td>
    <td>::child_general_int</td>
    <td>::GC_child_general_int</td>
  </tr>
  <tr>
    <td>General hop</td>
    <td>::child_general_hopp_GetInfo</td>
    <td>::child_general_hopp, ::mltply::Hubbard::C::general_hopp_MPIsingle, 
    ::mltply::Hubbard::C::general_hopp_MPIdouble</td>
    <td>::GC_child_general_hopp, ::mltply::Hubbard::GC::general_hopp_MPIsingle, 
    ::mltply::Hubbard::GC::general_hopp_MPIdouble</td>
  </tr>
</table>
*/
#include <bitcalc.hpp>
#include "mltplyCommon.hpp"
#include "mltplyHubbard.hpp"
#include "mltplyMPIHubbard.hpp"
#include "CalcTime.hpp"
#include "mltplyHubbardCore.hpp"
#include "mltplyMPIHubbardCore.hpp"
#include "global.hpp"
/**
@brief perform Hamiltonian vector product for (extended) Hubbard type model.
@f${\bf v}_0 = {\hat H}{\bf v}_1@f$
@return errorcode. 0 for normal, other error
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltply::Hubbard::C::main(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
){
  long int i;
  long int isite1, isite2, sigma1, sigma2;
  long int isite3, isite4, sigma3, sigma4;
  long int ibitsite1, ibitsite2, ibitsite3, ibitsite4;

  std::complex<double> tmp_trans;
  /*[s] For InterAll */
  std::complex<double> tmp_V;
  /*[e] For InterAll */

  int ihermite=0;
  int idx=0;

  StartTimer(300);
  /**
  Transfer
  */
  StartTimer(310);
  for (i = 0; i < Def::EDNTransfer; i+=2) {
    if (Def::EDGeneralTransfer[i][0] >= Def::Nsite &&
        Def::EDGeneralTransfer[i][2] >= Def::Nsite) {
      StartTimer(311);
      mltply::Hubbard::C::general_hopp_MPIdouble(i, nstate, tmp_v0, tmp_v1);
      StopTimer(311);
    }
    else if (Def::EDGeneralTransfer[i][2] >= Def::Nsite) {
      StartTimer(312);
      mltply::Hubbard::C::general_hopp_MPIsingle(i, nstate, tmp_v0, tmp_v1);
      StopTimer(312);
    }
    else if (Def::EDGeneralTransfer[i][0] >= Def::Nsite) {
      StartTimer(312);
      mltply::Hubbard::C::general_hopp_MPIsingle(i + 1, nstate, tmp_v0, tmp_v1);
      StopTimer(312);
    }
    else {
      StartTimer(313);
      for (ihermite = 0; ihermite<2; ihermite++) {
        idx = i + ihermite;
        isite1 = Def::EDGeneralTransfer[idx][0];
        isite2 = Def::EDGeneralTransfer[idx][2];
        sigma1 = Def::EDGeneralTransfer[idx][1];
        sigma2 = Def::EDGeneralTransfer[idx][3];
        mltply::Hubbard::general_hopp_GetInfo(isite1, isite2, sigma1, sigma2);
        tmp_trans = -Def::EDParaGeneralTransfer[idx];
        Large::tmp_trans = tmp_trans;
        mltply::Hubbard::C::general_hopp(nstate, tmp_v0, tmp_v1, tmp_trans);
      }
      StopTimer(313);
    }
  }/*for (i = 0; i < Def::EDNTransfer; i+=2)*/
  StopTimer(310);
  /**
  InterAll
  */
  StartTimer(320);
  for (i = 0; i < Def::NInterAll_OffDiagonal; i += 2) {
        
    isite1 = Def::InterAll_OffDiagonal[i][0];
    isite2 = Def::InterAll_OffDiagonal[i][2];
    isite3 = Def::InterAll_OffDiagonal[i][4];
    isite4 = Def::InterAll_OffDiagonal[i][6];
    sigma1 = Def::InterAll_OffDiagonal[i][1];
    sigma2 = Def::InterAll_OffDiagonal[i][3];
    sigma3 = Def::InterAll_OffDiagonal[i][5];
    sigma4 = Def::InterAll_OffDiagonal[i][7];
    tmp_V = Def::ParaInterAll_OffDiagonal[i];

    if (mltply::Hubbard::CheckPE(isite1) == TRUE || mltply::Hubbard::CheckPE(isite2) == TRUE ||
        mltply::Hubbard::CheckPE(isite3) == TRUE || mltply::Hubbard::CheckPE(isite4) == TRUE) {
      StartTimer(321);
      ibitsite1 = Def::OrgTpow[2 * isite1 + sigma1];
      ibitsite2 = Def::OrgTpow[2 * isite2 + sigma2];
      ibitsite3 = Def::OrgTpow[2 * isite3 + sigma3];
      ibitsite4 = Def::OrgTpow[2 * isite4 + sigma4];
      if (ibitsite1 == ibitsite2 && ibitsite3 == ibitsite4) {
        mltply::Hubbard::C::X_CisAisCjtAjt_MPI(isite1, sigma1,
          isite3, sigma3,
          tmp_V, nstate, tmp_v0, tmp_v1);
      }
      else if (ibitsite1 == ibitsite2 && ibitsite3 != ibitsite4) {
        mltply::Hubbard::C::X_CisAisCjtAku_MPI(isite1, sigma1,
          isite3, sigma3, isite4, sigma4,
          tmp_V, nstate, tmp_v0, tmp_v1);
      }
      else if (ibitsite1 != ibitsite2 && ibitsite3 == ibitsite4) {
        mltply::Hubbard::C::X_CisAjtCkuAku_MPI(isite1, sigma1, isite2, sigma2,
          isite3, sigma3, tmp_V, nstate, tmp_v0, tmp_v1);
      }
      else if (ibitsite1 != ibitsite2 && ibitsite3 != ibitsite4) {
        mltply::Hubbard::C::X_CisAjtCkuAlv_MPI(isite1, sigma1, isite2, sigma2,
          isite3, sigma3, isite4, sigma4, tmp_V, nstate, tmp_v0, tmp_v1);
      }
      StopTimer(321);
    }
    else {
      StartTimer(322);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        isite1 = Def::InterAll_OffDiagonal[idx][0];
        isite2 = Def::InterAll_OffDiagonal[idx][2];
        isite3 = Def::InterAll_OffDiagonal[idx][4];
        isite4 = Def::InterAll_OffDiagonal[idx][6];
        sigma1 = Def::InterAll_OffDiagonal[idx][1];
        sigma2 = Def::InterAll_OffDiagonal[idx][3];
        sigma3 = Def::InterAll_OffDiagonal[idx][5];
        sigma4 = Def::InterAll_OffDiagonal[idx][7];
        tmp_V = Def::ParaInterAll_OffDiagonal[idx];

        mltply::Hubbard::general_int_GetInfo(isite1, isite2, isite3, isite4,
          sigma1, sigma2, sigma3, sigma4, tmp_V);

        mltply::Hubbard::C::general_int(nstate, tmp_v0, tmp_v1);
      }/*for (ihermite = 0; ihermite < 2; ihermite++)*/
      StopTimer(322);
    }
  }/*for (i = 0; i < Def::NInterAll_OffDiagonal; i+=2)*/
  StopTimer(320);
  /**
  Pair hopping
  */
  StartTimer(330);
  for (i = 0; i < Def::NPairHopping; i +=2) {
    sigma1=0;
    sigma2=1;
        
    if ( Def::PairHopping[i][0] >= Def::Nsite 
      || Def::PairHopping[i][1] >= Def::Nsite)
    {
      StartTimer(331);
      mltply::Hubbard::C::X_CisAjtCkuAlv_MPI(
        Def::PairHopping[i][0], sigma1, Def::PairHopping[i][1], sigma1, 
        Def::PairHopping[i][0], sigma2, Def::PairHopping[i][1], sigma2, 
        Def::ParaPairHopping[i], nstate, tmp_v0, tmp_v1);
        StopTimer(331);
    }
    else {
      StartTimer(332);
      for (ihermite = 0; ihermite<2; ihermite++) {
        idx = i + ihermite;
        mltply::Hubbard::pairhopp_GetInfo(idx);
        mltply::Hubbard::C::pairhopp(nstate, tmp_v0, tmp_v1);
      }/*for (ihermite = 0; ihermite<2; ihermite++)*/
      StopTimer(332);
    }
  }/*for (i = 0; i < Def::NPairHopping; i += 2)*/
  StopTimer(330);  
  /**
  Exchange
  */
  StartTimer(340);
  for (i = 0; i < Def::NExchangeCoupling; i ++) {
    sigma1 = 0;
    sigma2 = 1;
    if (Def::ExchangeCoupling[i][0] >= Def::Nsite ||
        Def::ExchangeCoupling[i][1] >= Def::Nsite) 
    {
      StartTimer(341);
      mltply::Hubbard::C::X_CisAjtCkuAlv_MPI(
        Def::ExchangeCoupling[i][0], sigma1, Def::ExchangeCoupling[i][1], sigma1,
        Def::ExchangeCoupling[i][1], sigma2, Def::ExchangeCoupling[i][0], sigma2,
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(341);
    }
    else {
      StartTimer(342);
      mltply::Hubbard::exchange_GetInfo(i);
      mltply::Hubbard::C::exchange(nstate, tmp_v0, tmp_v1);
      StopTimer(342);
    }
  }/*for (i = 0; i < Def::NExchangeCoupling; i ++)*/
  StopTimer(340);

  StopTimer(300);
  return 0;
}/*int mltplyHubbard*/
/**
@brief perform Hamiltonian vector product for (extended) Hubbard type model
(Grandcanonical).
@f${\bf v}_0 = {\hat H}{\bf v}_1@f$
@return errorcode. 0 for normal, other error
*/
int mltply::Hubbard::GC::main(
  int nstate,/**<[in] Number of states*/
  std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
){
  long int i;
  long int isite1, isite2, sigma1, sigma2;
  long int isite3, isite4, sigma3, sigma4;
  long int ibitsite1, ibitsite2, ibitsite3, ibitsite4;

  std::complex<double> tmp_trans;
  /*[s] For InterAll */
  std::complex<double> tmp_V;
  /*[e] For InterAll */

  int ihermite=0;
  int idx=0;

  StartTimer(200);
  /**
  Transfer
  */
  StartTimer(210);
  for (i = 0; i < Def::EDNTransfer; i += 2) {
    if (Def::EDGeneralTransfer[i][0] >= Def::Nsite &&
        Def::EDGeneralTransfer[i][2] >= Def::Nsite) {
      StartTimer(211);
      mltply::Hubbard::GC::general_hopp_MPIdouble(i, nstate, tmp_v0, tmp_v1);
      StopTimer(211);
    }
    else if (Def::EDGeneralTransfer[i][2] >= Def::Nsite){
      StartTimer(212);
      mltply::Hubbard::GC::general_hopp_MPIsingle(i, nstate, tmp_v0, tmp_v1);
      StopTimer(212);
    }
    else if (Def::EDGeneralTransfer[i][0] >= Def::Nsite) {
      StartTimer(212);
      mltply::Hubbard::GC::general_hopp_MPIsingle(i+1, nstate, tmp_v0, tmp_v1);
      StopTimer(212);
    }
    else {
      StartTimer(213);
      for (ihermite = 0; ihermite<2; ihermite++) {
        idx = i + ihermite;
        isite1 = Def::EDGeneralTransfer[idx][0];
        isite2 = Def::EDGeneralTransfer[idx][2];
        sigma1 = Def::EDGeneralTransfer[idx][1];
        sigma2 = Def::EDGeneralTransfer[idx][3];
        mltply::Hubbard::general_hopp_GetInfo(isite1, isite2, sigma1, sigma2);
        tmp_trans = -Def::EDParaGeneralTransfer[idx];
        mltply::Hubbard::GC::general_hopp(nstate, tmp_v0, tmp_v1, tmp_trans);
      }
      StopTimer(213);
    }
  }/*for (i = 0; i < Def::EDNTransfer; i += 2)*/
  StopTimer(210);
  /**
  Inter All
  */
  StartTimer(220);
  for (i = 0; i < Def::NInterAll_OffDiagonal; i+=2) {
    isite1 = Def::InterAll_OffDiagonal[i][0];
    isite2 = Def::InterAll_OffDiagonal[i][2];
    isite3 = Def::InterAll_OffDiagonal[i][4];
    isite4 = Def::InterAll_OffDiagonal[i][6];
    sigma1 = Def::InterAll_OffDiagonal[i][1];
    sigma2 = Def::InterAll_OffDiagonal[i][3];
    sigma3 = Def::InterAll_OffDiagonal[i][5];
    sigma4 = Def::InterAll_OffDiagonal[i][7];
    tmp_V = Def::ParaInterAll_OffDiagonal[i];

    if ( mltply::Hubbard::CheckPE(isite1) == TRUE || mltply::Hubbard::CheckPE(isite2) == TRUE
      || mltply::Hubbard::CheckPE(isite3) == TRUE || mltply::Hubbard::CheckPE(isite4) == TRUE) 
    {
      StartTimer(221);
      ibitsite1 = Def::OrgTpow[2 * isite1 + sigma1];
      ibitsite2 = Def::OrgTpow[2 * isite2 + sigma2];
      ibitsite3 = Def::OrgTpow[2 * isite3 + sigma3];
      ibitsite4 = Def::OrgTpow[2 * isite4 + sigma4];
      if (ibitsite1 == ibitsite2 && ibitsite3 == ibitsite4) 
        mltply::Hubbard::GC::X_CisAisCjtAjt_MPI(
          isite1, sigma1, isite3, sigma3, tmp_V, nstate, tmp_v0, tmp_v1);
      else if (ibitsite1 == ibitsite2 && ibitsite3 != ibitsite4) 
        mltply::Hubbard::GC::X_CisAisCjtAku_MPI(
          isite1, sigma1, isite3, sigma3, isite4, sigma4, tmp_V, nstate, tmp_v0, tmp_v1);
      else if (ibitsite1 != ibitsite2 && ibitsite3 == ibitsite4) 
        mltply::Hubbard::GC::X_CisAjtCkuAku_MPI(
          isite1, sigma1, isite2, sigma2, isite3, sigma3, tmp_V, nstate, tmp_v0, tmp_v1);
      else if (ibitsite1 != ibitsite2 && ibitsite3 != ibitsite4) 
        mltply::Hubbard::GC::X_CisAjtCkuAlv_MPI(
          isite1, sigma1, isite2, sigma2, isite3, sigma3, isite4, sigma4, tmp_V, nstate, tmp_v0, tmp_v1);
      StopTimer(221);
    }//InterPE
    else{
      StartTimer(222);
      for(ihermite=0; ihermite<2; ihermite++){
        idx=i+ihermite;
        isite1 = Def::InterAll_OffDiagonal[idx][0];
        isite2 = Def::InterAll_OffDiagonal[idx][2];
        isite3 = Def::InterAll_OffDiagonal[idx][4];
        isite4 = Def::InterAll_OffDiagonal[idx][6];
        sigma1 = Def::InterAll_OffDiagonal[idx][1];
        sigma2 = Def::InterAll_OffDiagonal[idx][3];
        sigma3 = Def::InterAll_OffDiagonal[idx][5];
        sigma4 = Def::InterAll_OffDiagonal[idx][7];
        tmp_V = Def::ParaInterAll_OffDiagonal[idx];
          
        mltply::Hubbard::general_int_GetInfo(isite1, isite2, isite3, isite4,
                                        sigma1, sigma2, sigma3, sigma4, tmp_V); 
        mltply::Hubbard::GC::general_int(nstate, tmp_v0, tmp_v1);
      }/*for(ihermite=0; ihermite<2; ihermite++)*/
      StopTimer(222);
    }
  }/*for (i = 0; i < Def::NInterAll_OffDiagonal; i+=2)*/
  StopTimer(220);
  /**
  Pair hopping
  */
  StartTimer(230);
  for (i = 0; i < Def::NPairHopping; i +=2) {
    sigma1 = 0;
    sigma2 = 1;
    if ( Def::PairHopping[i][0] >= Def::Nsite
      || Def::PairHopping[i][1] >= Def::Nsite) 
    {
      StartTimer(231);
      mltply::Hubbard::GC::X_CisAjtCkuAlv_MPI(
        Def::PairHopping[i][0], sigma1, Def::PairHopping[i][1], sigma1,
        Def::PairHopping[i][0], sigma2, Def::PairHopping[i][1], sigma2,
        Def::ParaPairHopping[i], nstate, tmp_v0, tmp_v1);
      StopTimer(231);
    }
    else {
      StartTimer(232);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        mltply::Hubbard::pairhopp_GetInfo(idx);
        mltply::Hubbard::GC::pairhopp(nstate, tmp_v0, tmp_v1);
      }/*for (ihermite = 0; ihermite<2; ihermite++)*/
      StopTimer(232);
    }
  }/*for (i = 0; i < Def::NPairHopping; i += 2)*/
  StopTimer(230);
  /**
  Exchange
  */
  StartTimer(240);
  for (i = 0; i < Def::NExchangeCoupling; i++) {
    sigma1=0;
    sigma2=1;
    if ( Def::ExchangeCoupling[i][0] >= Def::Nsite
      || Def::ExchangeCoupling[i][1] >= Def::Nsite) 
    {
      StartTimer(241);
      mltply::Hubbard::GC::X_CisAjtCkuAlv_MPI(
        Def::ExchangeCoupling[i][0], sigma1, Def::ExchangeCoupling[i][1], sigma1,
        Def::ExchangeCoupling[i][1], sigma2, Def::ExchangeCoupling[i][0], sigma2,
        Def::ParaExchangeCoupling[i], nstate, tmp_v0, tmp_v1);
      StopTimer(241);
    }
    else {
      StartTimer(242);
      mltply::Hubbard::exchange_GetInfo(i);
      mltply::Hubbard::GC::exchange(nstate, tmp_v0, tmp_v1);
      StopTimer(242);
    }
  }/*for (i = 0; i < Def::NExchangeCoupling; i++)*/
  StopTimer(240);

  StopTimer(200);
  return 0;
}/*int mltplyHubbardGC*/

/******************************************************************************/
//[s] child functions
/******************************************************************************/

/**
@brief Compute pairhopp term (canonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::C::pairhopp(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i_max = Large::i_max;
  long int off = 0;

#pragma omp parallel for default(none) private(j, off) \
shared(tmp_v0, tmp_v1,nstate,i_max)
  for (j = 0; j < i_max; j++) 
    mltply::Hubbard::C::pairhopp_element(j, nstate, tmp_v0, tmp_v1, &off);
  return;
}/*std::complex<double> child_pairhopp*/
/**
@brief Compute Exchange term (canonical) in single process
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::C::exchange(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i_max = Large::i_max;
  long int off = 0;

#pragma omp parallel for default(none) private(j, off) \
shared(tmp_v0, tmp_v1,nstate,i_max)
  for (j = 0; j < i_max; j++) 
    mltply::Hubbard::C::exchange_element(j, nstate, tmp_v0, tmp_v1, &off);
  return;
}/*std::complex<double> child_exchange*/
/**
@brief Compute hopping (canonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::C::general_hopp(
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1,//!<[in] Input producted vector
  std::complex<double> trans//!<[in] Hopping integral
) {
  long int j, isite1, isite2, Asum, Adiff;
  long int i_max = Large::i_max;

  isite1 = Large::is1_spin;
  isite2 = Large::is2_spin;
  Asum = Large::isA_spin;
  Adiff = Large::A_spin;
#pragma omp parallel for default(none) private(j) \
shared(tmp_v0, tmp_v1,nstate,i_max,Asum,Adiff,isite1,isite2,trans) 
  for (j = 0; j < i_max; j++)
    CisAjt(j, nstate, tmp_v0, tmp_v1, isite1, isite2, Asum, Adiff, trans);
  return;
}/*std::complex<double> child_general_hopp*/
/**
@brief Commpute hopping term (grandcanonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::GC::general_hopp(
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1,//!<[in] Input producted vector
  std::complex<double> trans//!<[in] Hopping integral
) {
  long int j, isite1, isite2, Asum, Adiff;
  long int tmp_off;
  long int i_max = Large::i_max;

  isite1 = Large::is1_spin;
  isite2 = Large::is2_spin;
  Asum = Large::isA_spin;
  Adiff = Large::A_spin;

  if (isite1 == isite2) {
#pragma omp parallel for default(none) private(j) \
shared(tmp_v0, tmp_v1,nstate,i_max,isite1, trans)
    for (j = 0; j < i_max; j++)
      mltply::Hubbard::GC::CisAis(j, nstate, tmp_v0, tmp_v1, isite1, trans);
  }/*if (isite1 == isite2)*/
  else {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(tmp_v0,tmp_v1,nstate,i_max,Asum,Adiff,isite1,isite2,trans)
    for (j = 0; j < i_max; j++) 
      mltply::Hubbard::GC::CisAjt(j, nstate, tmp_v0, tmp_v1, isite1, isite2, Asum, Adiff, trans, &tmp_off);
  }
  return;
}/*std::complex<double> GC_child_general_hopp*/
/**
@brief Compute inter-all term (canonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::C::general_int(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  std::complex<double> tmp_V;
  long int j, i_max;
  long int isite1, isite2, isite3, isite4;
  long int Asum, Bsum, Adiff, Bdiff;
  long int tmp_off = 0;
  long int tmp_off_2 = 0;

  //note: this site is labeled for not only site but site with spin.
  i_max = Large::i_max;
  isite1 = Large::is1_spin;
  isite2 = Large::is2_spin;
  Asum = Large::isA_spin;
  Adiff = Large::A_spin;

  isite3 = Large::is3_spin;
  isite4 = Large::is4_spin;
  Bsum = Large::isB_spin;
  Bdiff = Large::B_spin;

  tmp_V = Large::tmp_V;

#pragma omp parallel default(none) private(j, tmp_off, tmp_off_2) \
shared(i_max, isite1, isite2, isite3, isite4, Asum, Bsum, Adiff, Bdiff, tmp_V) \
  shared(tmp_v0, tmp_v1,nstate)
  {
    if (isite1 == isite2 && isite3 == isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::C::CisAisCisAis_element(j, isite1, isite3, tmp_V, nstate, tmp_v0, tmp_v1);
    }/*if (isite1 == isite2 && isite3 == isite4)*/
    else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::C::CisAisCjtAku_element(
          j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off);
    }/*if (isite1 == isite2 && isite3 != isite4)*/
    else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::C::CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off);
    }/*if (isite1 != isite2 && isite3 == isite4)*/
    else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::C::CisAjtCkuAlv_element(
          j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off_2);
    }/*if (isite1 != isite2 && isite3 != isite4)*/
  }/*End of parallel region*/
  return;
}/*std::complex<double> child_general_int*/
/**
@brief Compute inter-all term (canonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::GC::general_int(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  std::complex<double> tmp_V;
  long int j, i_max;
  long int isite1, isite2, isite3, isite4;
  long int Asum, Bsum, Adiff, Bdiff;
  long int tmp_off = 0;

  i_max = Large::i_max;
  isite1 = Large::is1_spin;
  isite2 = Large::is2_spin;
  Asum = Large::isA_spin;
  Adiff = Large::A_spin;

  isite3 = Large::is3_spin;
  isite4 = Large::is4_spin;
  Bsum = Large::isB_spin;
  Bdiff = Large::B_spin;

  tmp_V = Large::tmp_V;

#pragma omp parallel default(none)  private(j, tmp_off) \
shared(i_max, isite1, isite2, isite4, isite3, Asum, Bsum, Adiff, Bdiff, tmp_V) \
shared(tmp_v0, tmp_v1,nstate)
  {
    if (isite1 == isite2 && isite3 == isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::GC::CisAisCisAis_element(j, isite1, isite3, tmp_V, 
          nstate, tmp_v0, tmp_v1);
    }/*if (isite1 == isite2 && isite3 == isite4)*/
    else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::GC::CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V,
          nstate, tmp_v0, tmp_v1, &tmp_off);
    }/*if (isite1 == isite2 && isite3 != isite4)*/
    else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::GC::CisAjtCkuAku_element(
          j, isite1, isite2, isite3, Asum, Adiff, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off);
    }/*if (isite1 != isite2 && isite3 == isite4)*/
    else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::GC::CisAjtCkuAlv_element(
          j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, nstate, tmp_v0, tmp_v1, &tmp_off);
    }/*if (isite1 != isite2 && isite3 != isite4)*/
  }/*End of parallel region*/
  return;
}/*std::complex<double> GC_child_general_int*/
/**
@brief Compute pairhopp term (grandcanonical)
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::GC::pairhopp(
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Result vector
  std::complex<double> **tmp_v1//!<[in] Input producted vector
) {
  long int j;
  long int i_max = Large::i_max;
  long int off = 0;

#pragma omp parallel for default(none) private(j,off) \
shared(tmp_v0, tmp_v1,nstate,i_max)
  for (j = 0; j < i_max; j++) 
    mltply::Hubbard::GC::pairhopp_element(j, nstate, tmp_v0, tmp_v1, &off);

  return;
}/*std::complex<double> GC_child_pairhopp*/
/**
@brief Compute Exchange term (grandcanonical) in single process
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::GC::exchange(
  int nstate, std::complex<double> **tmp_v0,
  std::complex<double> **tmp_v1
) {
  long int j;
  long int i_max = Large::i_max;
  long int off = 0;

#pragma omp parallel for default(none) private(j, off) \
shared(tmp_v0, tmp_v1,nstate,i_max)
  for (j = 0; j < i_max; j++) 
    mltply::Hubbard::GC::exchange_element(j, nstate, tmp_v0, tmp_v1, &off);
  return;
}/*std::complex<double> GC_child_exchange*/
/******************************************************************************/
//[e] child functions
/******************************************************************************/
