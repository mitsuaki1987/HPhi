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
@brief Functions for Hubbar + MPI (Core)
*/
#include "mltplyCommon.hpp"
#include "mltplyHubbardCore.hpp"
#include "mltplyMPIHubbard.hpp"
#include "mltplyMPIHubbardCore.hpp"
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "global.hpp"
/**
@brief Check whether this site is in the inter process region or not
@return 1 if it is inter-process region, 0 if not.
*/
int mltply::Hubbard::CheckPE(
  int org_isite//!<[in] Site index
){
  if (org_isite >= Def::Nsite) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}/*int mltply::Hubbard::CheckPE*/
/**
@brief Check the occupation of @f$(i,s)@f$ state,
and compute the index of final wavefunction associated to 
@f$c^\dagger_{is}@f$
@return 1 if unoccupied, 0 if occupied
*/
int mltply::Hubbard::CheckBit_Cis(
  long int is1_spin,//!<[in] Index of site+spin
  long int orgbit,//!<[in] Index of initial wavefunction
  long int *offbit//!<[out] Index of final wavefunction
) {
  long int ibit_tmp;
  ibit_tmp = orgbit & is1_spin;
  if (ibit_tmp == 0) {
    *offbit = orgbit + is1_spin;
    return TRUE;
  }
  *offbit = 0;
  return FALSE;
}/*int mltply::Hubbard::CheckBit_Cis*/
/**
@brief Check the occupation of @f$(i,s)@f$ state,
and compute the index of final wavefunction associated to
@f$c_{jt}@f$
@return 1 if occupied, 0 if unoccupied
*/
int mltply::Hubbard::CheckBit_Ajt(
  long int is1_spin,//!<[in] Index of site+spin
  long int orgbit,//!<[in] Index of initial wavefunction
  long int *offbit//!<[out] Index of final wavefunction
) {
  long int ibit_tmp;
  ibit_tmp = orgbit & is1_spin;
  if (ibit_tmp != 0) {
    *offbit = orgbit - is1_spin;
    return TRUE;
  }
  *offbit = 0;
  return FALSE;
}/*int mltply::Hubbard::CheckBit_Ajt*/
/**
@brief Compute the index of final wavefunction associated to
@f$c_{4}^\dagger c_{3}c_{2}^\dagger c_{1}@f$, and
check whether this operator is relevant or not
@return 1 if relevant, 0 if irrelevant
*/
int mltply::Hubbard::CheckBit_InterAllPE(
  int org_isite1,//!<[in] Site 1
  int org_isigma1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_isigma2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_isigma3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_isigma4,//!<[in] Spin 4
  long int orgbit,//!<[in] Index of initial wavefunction
  long int *offbit//!<[out] Index of final wavefunction
){
  long int tmp_ispin;
  long int tmp_org, tmp_off;
  int iflgBitExist = TRUE;
  tmp_org=orgbit;
  tmp_off=0;

  if (mltply::Hubbard::CheckPE(org_isite1) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite1 + org_isigma1];
    if (mltply::Hubbard::CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (mltply::Hubbard::CheckPE(org_isite2) == TRUE ) {
    tmp_ispin = Def::Tpow[2 * org_isite2 + org_isigma2];
    if (mltply::Hubbard::CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (mltply::Hubbard::CheckPE(org_isite3) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite3 + org_isigma3];
    if (mltply::Hubbard::CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (mltply::Hubbard::CheckPE(org_isite4) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite4 + org_isigma4];
    if (mltply::Hubbard::CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if(iflgBitExist != TRUE){
    *offbit=0;
    return FALSE;
  }
  
  *offbit=tmp_org;
  return TRUE;
}/*int mltply::Hubbard::CheckBit_InterAllPE*/
/**
@brief Check the occupation of both site 1 and site 3
@return 1 if both sites are occupied, 0 if not
*/
int mltply::Hubbard::CheckBit_PairPE(
  int org_isite1,//!<[in] Site 1
  int org_isigma1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_isigma3,//!<[in] Spin 4
  long int orgbit//!<[in] Index pf intial wavefunction
){
  long int tmp_ispin;
  long int tmp_org, tmp_off;
  int iflgBitExist = TRUE;
  tmp_org=orgbit;
  
  if(mltply::Hubbard::CheckPE(org_isite1)==TRUE){
    tmp_ispin = Def::Tpow[2 * org_isite1 + org_isigma1];
    if (mltply::Hubbard::CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist=FALSE;
    }
  }
  
  if (mltply::Hubbard::CheckPE(org_isite3) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite3 + org_isigma3];
    if (mltply::Hubbard::CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
  }
  
  if(iflgBitExist != TRUE){
    return FALSE;
  }

  return TRUE;
}/*int mltply::Hubbard::CheckBit_PairPE*/
/**
@brief Compute the index of final wavefunction associated to
@f$c_{4}^\dagger c_{3}c_{2}^\dagger c_{1}@f$, and
Fermion sign
@return 1 if relevant, 0 if irrelevant
*/
int mltply::Hubbard::GetSgnInterAll(
  long int isite1,//!<[in] Site 1
  long int isite2,//!<[in] Site 2
  long int isite3,//!<[in] Site 3
  long int isite4,//!<[in] Site 4
  int *Fsgn,//!<[out] Fermion sign
  long int orgbit,//!<[in] Index of the initial state
  long int *offbit//!<[out] Index of the final state
){
  long int diffA;
  long int tmp_off;
  long int tmp_ispin1, tmp_ispin2;
  int tmp_sgn=0;

  tmp_ispin1=isite2;
  tmp_ispin2=isite1;

  if (tmp_ispin1 == tmp_ispin2) {
    if ((orgbit & tmp_ispin1) == 0) {
      *offbit = 0;
      *Fsgn = tmp_sgn;
      return FALSE;
    }
    tmp_sgn=1;
    tmp_off = orgbit;
  }
  else {
    if (tmp_ispin2 > tmp_ispin1) diffA = tmp_ispin2 - tmp_ispin1 * 2;
    else diffA = tmp_ispin1-tmp_ispin2*2;

    tmp_sgn= mltply::Hubbard::GC::X_CisAjt(orgbit, tmp_ispin1, tmp_ispin2,
      tmp_ispin1+tmp_ispin2, diffA, &tmp_off);
    if(tmp_sgn ==0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
  }

  tmp_ispin1 = isite4;
  tmp_ispin2 = isite3;

  if(tmp_ispin1 == tmp_ispin2){
    if( (tmp_off & tmp_ispin1) == 0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
    *offbit=tmp_off;
  }
  else{
    if(tmp_ispin2 > tmp_ispin1) diffA = tmp_ispin2 - tmp_ispin1*2;
    else diffA = tmp_ispin1-tmp_ispin2*2;
    tmp_sgn *= mltply::Hubbard::GC::X_CisAjt(tmp_off, tmp_ispin1, tmp_ispin2, 
      tmp_ispin1+tmp_ispin2, diffA, offbit);
    
    if(tmp_sgn ==0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
  }
  
  *Fsgn =tmp_sgn;
  *offbit = *offbit%Def::OrgTpow[2*Def::Nsite];
  return TRUE;
}/*int mltply::Hubbard::GetSgnInterAll*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{jt}@f$
term of grandcanonical Hubbard system
*/
void mltply::Hubbard::GC::X_CisAisCjtAjt_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  int iCheck;
  long int tmp_ispin1;
  long int i_max = Check::idim_max;
  long int tmp_off, j;
  int one = 1;

  iCheck=mltply::Hubbard::CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, (long int) MP::myrank);
  if(iCheck != TRUE){
    return;
  }

  if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite) {
    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[0][0], &tmp_v0[0][0]);
  }/*if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite)*/
  else if (org_isite1 >= Def::Nsite || org_isite3 >= Def::Nsite) {
    if (org_isite1 > org_isite3) tmp_ispin1 = Def::Tpow[2 * org_isite3 + org_ispin3];
    else                         tmp_ispin1 = Def::Tpow[2 * org_isite1 + org_ispin1];

#pragma omp parallel for default(none) private(j,tmp_off) \
shared(org_isite1,org_ispin1,org_isite3,org_ispin3,nstate,one,tmp_v0,tmp_v1,tmp_ispin1, i_max,tmp_V)
    for (j = 0; j < i_max; j++) {
      if (mltply::Hubbard::CheckBit_Ajt(tmp_ispin1, j, &tmp_off) == TRUE) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }
    }/*for (j = 0; j < i_max; j++)*/
  }
}/*std::complex<double> mltply::Hubbard::GC::X_CisAisCjtAjt_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
*/
void mltply::Hubbard::GC::X_CisAjtCkuAku_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  int iCheck, Fsgn;
  long int isite1, isite2, isite3;
  long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  long int j, Asum, Adiff;
  std::complex<double> dmv;
  long int origin, tmp_off;
  long int org_rankbit;
  int one = 1;

  iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3, org_ispin3, (long int) MP::myrank, &origin);
  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  isite2 = Def::Tpow[2 * org_isite2 + org_ispin2];
  isite3 = Def::Tpow[2 * org_isite3 + org_ispin3];

  if (iCheck == TRUE) {
    tmp_isite1 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
    Asum = tmp_isite1 + tmp_isite2;
    if (tmp_isite2 > tmp_isite1) Adiff = tmp_isite2 - tmp_isite1 * 2;
    else Adiff = tmp_isite1 - tmp_isite2 * 2;
  }
  else {
    iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1, org_ispin1, (long int) MP::myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
      Asum = tmp_isite3 + tmp_isite4;
      if (tmp_isite4 > tmp_isite3) Adiff = tmp_isite4 - tmp_isite3 * 2;
      else Adiff = tmp_isite3 - tmp_isite4 * 2;
      if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
        tmp_V = 0;
      }
    }
    else {
      return;
    }
  }

  if (MP::myrank == origin) {// only k is in PE

    if (mltply::Hubbard::CheckBit_Ajt(isite3, MP::myrank, &tmp_off) == FALSE) return;

#pragma omp parallel default(none) private(j,tmp_off) \
shared(i_max,Asum,Adiff,isite1,isite2, tmp_V,tmp_v0, tmp_v1,nstate,Large::mode)
    {
#pragma omp for
      for (j = 0; j < i_max; j++) 
        mltply::Hubbard::GC::CisAjt(j, nstate, tmp_v0, tmp_v1, isite2, isite1, Asum, Adiff, 
          tmp_V, &tmp_off);

      if (Large::mode != M_CORR) {
#pragma omp for
        for (j = 0; j < i_max; j++) 
          mltply::Hubbard::GC::CisAjt(j, nstate, tmp_v0, tmp_v1, isite1, isite2, Asum, Adiff, 
            tmp_V, &tmp_off);
      }/*if (Large::mode != M_CORR)*/
    }/*End of paralle region*/
    return;
  }//MP::myrank =origin
  else {
    wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

#pragma omp parallel default(none) private(j,dmv,tmp_off,Fsgn,org_rankbit,Adiff) \
shared(Wave::v1buf,tmp_v1,nstate,one,tmp_v0,MP::myrank,origin,isite3,org_isite3,isite1,isite2, Def::OrgTpow, \
org_isite2,org_isite1,Check::idim_max,tmp_V,tmp_isite1,tmp_isite2,tmp_isite3,tmp_isite4, Def::Nsite)
    {
      if (org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite) {
        if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
        else Adiff = isite1 - isite2 * 2;
        SgnBit(((long int) MP::myrank & Adiff), &Fsgn);
        tmp_V *= Fsgn;

        if (org_isite3 >= Def::Nsite) {
#pragma omp for
          for (j = 0; j < Check::idim_max; j++) {
            zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[j][0], &one);
          }/*for (j = 0; j < idim_max_buf; j++)*/
        }
        else { //org_isite3 <= Def::Nsite
#pragma omp for
          for (j = 0; j < Check::idim_max; j++) {
            if (mltply::Hubbard::CheckBit_Ajt(isite3, j, &tmp_off) == TRUE) {
              zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[j][0], &one);
            }
          }/*for (j = 0; j < idim_max_buf; j++)*/
        }
      }/*if (org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite)*/
      else {
        org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;
#pragma omp for
        for (j = 0; j < Check::idim_max; j++) {
          if (mltply::Hubbard::GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, j + org_rankbit, &tmp_off) == TRUE) {
            dmv = tmp_V * (std::complex<double>)Fsgn;
            zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[tmp_off][0], &one);
          }
        }/*for (j = 0; j < idim_max_buf; j++)*/
      }
    }/*End of parallel region*/
  }/*MP::myrank != origin*/
}/*std::complex<double> mltply::Hubbard::GC::X_CisAjtCkuAku_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
*/
void mltply::Hubbard::GC::X_CisAisCjtAku_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  mltply::Hubbard::GC::X_CisAjtCkuAku_MPI(
    org_isite4, org_ispin4, org_isite3, org_ispin3,
    org_isite1, org_ispin1, conj(tmp_V), nstate, tmp_v0, tmp_v1);
}/*std::complex<double> mltply::Hubbard::GC::X_CisAisCjtAku_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of grandcanonical Hubbard system
*/
void mltply::Hubbard::GC::X_CisAjtCkuAlv_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  int iCheck, Fsgn;
  long int isite1, isite2, isite3, isite4;
  long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  long int j, Adiff, Bdiff;
  std::complex<double> dmv;
  long int origin, tmp_off, tmp_off2;
  long int org_rankbit;
  int iFlgHermite = FALSE;
  int one = 1;

  iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                               org_isite3, org_ispin3, org_isite4, org_ispin4,
                               (long int) MP::myrank, &origin);
  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  isite2 = Def::Tpow[2 * org_isite2 + org_ispin2];
  isite3 = Def::Tpow[2 * org_isite3 + org_ispin3];
  isite4 = Def::Tpow[2 * org_isite4 + org_ispin4];

  if (iCheck == TRUE) {
    tmp_isite1 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = Def::OrgTpow[2 * org_isite4 + org_ispin4];
  }
  else {
    iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                 org_isite2, org_ispin2, org_isite1, org_ispin1,
                                 (long int) MP::myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = Def::OrgTpow[2 * org_isite4 + org_ispin4];
      iFlgHermite = TRUE;
      if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
        tmp_V = 0;
      }
    }
    else {
      return;
    }
  }

  if (MP::myrank == origin) {
    if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
      mltply::Hubbard::GC::X_CisAis_MPI(org_isite1, org_ispin1, tmp_V, nstate, tmp_v0, tmp_v1);
      //calc -nisniv
      mltply::Hubbard::GC::X_CisAisCjtAjt_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);
    }/*if (isite1 == isite4 && isite2 == isite3)*/
    else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
            //calc CisAku
      if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
      else Adiff = isite1 - isite4 * 2;

#pragma omp parallel for default(none) private(j, tmp_off) \
shared(i_max, tmp_V, isite1, isite4, Adiff, tmp_v1, tmp_v0,nstate)
      for (j = 0; j < i_max; j++) 
        mltply::Hubbard::GC::CisAjt(j, nstate, tmp_v0, tmp_v1, 
          isite1, isite4, (isite1 + isite4), Adiff, tmp_V, &tmp_off);
      
      //calc -CisAku njv
      mltply::Hubbard::GC::X_CisAjtCkuAku_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4, 
                                          org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);
      if (Large::mode != M_CORR) { //for hermite
#pragma omp parallel for default(none) private(j, tmp_off) \
shared(i_max, tmp_V, isite1, isite4, Adiff, tmp_v1, tmp_v0,nstate)
        for (j = 0; j < i_max; j++) 
          mltply::Hubbard::GC::CisAjt(j, nstate, tmp_v0, tmp_v1, 
            isite4, isite1, (isite1 + isite4), Adiff, tmp_V, &tmp_off);
        
        //calc -njvCkuAis
        mltply::Hubbard::GC::X_CisAisCjtAku_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4,
                                            org_isite1, org_ispin1, -tmp_V, nstate, tmp_v0, tmp_v1);
      }/*if (Large::mode != M_CORR)*/
    }/*if (isite2 == isite3)*/
    else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE
      mltply::Hubbard::GC::X_CisAisCjtAku_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3,
                                          org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);
      if (Large::mode != M_CORR) { //for hermite
        mltply::Hubbard::GC::X_CisAisCjtAku_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                            org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);
      }/*if (Large::mode != M_CORR)*/
    }/*if (isite2 != isite3)*/
    return;
  }//MP::myrank =origin
  else {
    wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

    if (org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite
     && org_isite3 >= Def::Nsite && org_isite4 >= Def::Nsite) {

      if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
      else Adiff = isite1 - isite2 * 2;
      if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
      else Bdiff = isite3 - isite4 * 2;

      if (iFlgHermite == FALSE) {
        Fsgn = mltply::Hubbard::GC::X_CisAjt((long int) MP::myrank, 
          isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
        Fsgn *= mltply::Hubbard::GC::X_CisAjt(tmp_off2, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == FALSE)*/
      else {
        Fsgn = mltply::Hubbard::GC::X_CisAjt((long int) MP::myrank,
          isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
        Fsgn *= mltply::Hubbard::GC::X_CisAjt(tmp_off2, isite1, isite2, 
          (isite1 + isite2), Adiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == TRUE)*/

      zaxpy_long(i_max*nstate, tmp_V, &Wave::v1buf[0][0], &tmp_v0[0][0]);
    }
    else {
      org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;
#pragma omp parallel for default(none) private(j,dmv,tmp_off,Fsgn) \
shared(Check::idim_max,tmp_V,tmp_isite1,tmp_isite2,tmp_isite3,tmp_isite4,org_rankbit,Wave::v1buf,tmp_v1,tmp_v0,nstate,one)
      for (j = 0; j < Check::idim_max; j++) {
        if (mltply::Hubbard::GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, j + org_rankbit, &tmp_off) == TRUE) {
          dmv = tmp_V * (std::complex<double>)Fsgn;
          zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[tmp_off][0], &one);
        }
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }
  }/*MP::myrank != origin*/
}/*std::complex<double> mltply::Hubbard::GC::X_CisAjtCkuAlv_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}@f$
term of grandcanonical Hubbard system
*/
void mltply::Hubbard::GC::X_CisAis_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  long int j, isite1, tmp_off;
  int one = 1;

  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  if (org_isite1 >= Def::Nsite) {
    if (mltply::Hubbard::CheckBit_Ajt(isite1, (long int) MP::myrank, &tmp_off) == FALSE) return;

    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[0][0], &tmp_v0[0][0]);
  }/*if (org_isite1 >= Def::Nsite)*/
  else {
#pragma omp parallel for default(none) private(j, tmp_off) \
shared(tmp_v0, tmp_v1,nstate,one, i_max, tmp_V, isite1)
    for (j = 0; j < i_max; j++) {
      if (mltply::Hubbard::CheckBit_Ajt(isite1, j, &tmp_off) == TRUE) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }/*if (mltply::Hubbard::CheckBit_Ajt(isite1, j, &tmp_off) == TRUE)*/
    }/*for (j = 0; j < i_max; j++)*/
  }/*if (org_isite1 <= Def::Nsite)*/
}/*std::complex<double> mltply::Hubbard::GC::X_CisAis_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt}@f$
term of grandcanonical Hubbard system
*/
void mltply::Hubbard::GC::X_CisAjt_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {

  if (org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite) {
    mltply::Hubbard::GC::X_general_hopp_MPIdouble(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, nstate, tmp_v0, tmp_v1);
  }
  else if (org_isite1 >= Def::Nsite || org_isite2 >= Def::Nsite) {
    mltply::Hubbard::GC::X_general_hopp_MPIsingle(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, nstate, tmp_v0, tmp_v1);
  }
  else {
    //error message will be added.
    wrapperMPI::Exit(-1);
  }
}/*std::complex<double> mltply::Hubbard::GC::X_CisAjt_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{jt}@f$
term of canonical Hubbard system
*/
void mltply::Hubbard::C::X_CisAisCjtAjt_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  int iCheck;
  long int tmp_ispin1;
  long int i_max = Check::idim_max;
  long int tmp_off, j;
  int one = 1;

  iCheck = mltply::Hubbard::CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, (long int) MP::myrank);
  if (iCheck != TRUE) return;
  
  if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite) {
    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[0][0], &tmp_v0[0][0]);
  }/*if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite)*/
  else if (org_isite1 >= Def::Nsite || org_isite3 >= Def::Nsite) {
    if (org_isite1 > org_isite3) tmp_ispin1 = Def::Tpow[2 * org_isite3 + org_ispin3];
    else                         tmp_ispin1 = Def::Tpow[2 * org_isite1 + org_ispin1];

#pragma omp parallel for default(none) private(j,tmp_off) \
shared(tmp_v0,tmp_v1,List::c1,org_isite1,org_ispin1,org_isite3,org_ispin3,nstate,one, i_max,tmp_V,tmp_ispin1)
    for (j = 0; j < i_max; j++) {
      if (mltply::Hubbard::CheckBit_Ajt(tmp_ispin1, List::c1[j], &tmp_off) == TRUE) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }
    }/*for (j = 0; j < i_max; j++)*/
  }/*if (org_isite1 >= Def::Nsite || org_isite3 >= Def::Nsite)*/
}/*std::complex<double> mltply::Hubbard::C::X_CisAisCjtAjt_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of canonical Hubbard system
*/
void mltply::Hubbard::C::X_CisAjtCkuAlv_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  long int idim_max_buf;
  int iCheck, Fsgn;
  long int isite1, isite2, isite3, isite4;
  long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  long int j, Adiff, Bdiff;
  std::complex<double> dmv;
  long int origin, tmp_off, tmp_off2;
  long int org_rankbit, ioff;
  int iFlgHermite = FALSE;
  int one = 1;

  iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                               org_isite3, org_ispin3, org_isite4, org_ispin4,
                               (long int) MP::myrank, &origin);
  //printf("iCheck=%d, MP::myrank=%d, origin=%d\n", iCheck, MP::myrank, origin);
  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  isite2 = Def::Tpow[2 * org_isite2 + org_ispin2];
  isite3 = Def::Tpow[2 * org_isite3 + org_ispin3];
  isite4 = Def::Tpow[2 * org_isite4 + org_ispin4];

  if (iCheck == TRUE) {
    tmp_isite1 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = Def::OrgTpow[2 * org_isite4 + org_ispin4];
  }/*if (iCheck == TRUE)*/
  else {
    iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                 org_isite2, org_ispin2, org_isite1, org_ispin1,
                                 (long int) MP::myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = Def::OrgTpow[2 * org_isite4 + org_ispin4];
      iFlgHermite = TRUE;
      if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0;     
    }/*if (iCheck == TRUE)*/
    else return;
  }/*if (iCheck == FALSE)*/

  if (MP::myrank == origin) {
    if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
      mltply::Hubbard::C::X_CisAis_MPI(org_isite1, org_ispin1, tmp_V, nstate, tmp_v0, tmp_v1);
      //calc -nisniv
      mltply::Hubbard::C::X_CisAisCjtAjt_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);
    }/* if (isite1 == isite4 && isite2 == isite3)*/
    else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
      if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
      else Adiff = isite1 - isite4 * 2;

      //calc CisAku
#pragma omp parallel for default(none) private(j, tmp_off) \
shared(i_max, tmp_V, isite1, isite4, Adiff, tmp_v1, nstate, tmp_v0, List::c1)
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::C::CisAjt(j, nstate, tmp_v0, tmp_v1,
          isite1, isite4, (isite1 + isite4), Adiff, tmp_V);
      
      //calc -CisAku njv
      mltply::Hubbard::C::X_CisAjtCkuAku_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4,
                                       org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);

      if (Large::mode != M_CORR) {  //for hermite
#pragma omp parallel for default(none) private(j, tmp_off) \
shared(i_max, tmp_V, isite1, isite4, Adiff,tmp_v1, tmp_v0,nstate)
        for (j = 0; j < i_max; j++) 
          mltply::Hubbard::C::CisAjt(j, nstate, tmp_v0, tmp_v1,
            isite4, isite1, (isite1 + isite4), Adiff, tmp_V);
                
        //calc -njvCkuAis
        mltply::Hubbard::C::X_CisAisCjtAku_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4, 
                                         org_isite1, org_ispin1, -tmp_V, nstate, tmp_v0, tmp_v1);
      }/*if (Large::mode != M_CORR)*/
    }/*if (isite2 == isite3)*/
    else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE
      mltply::Hubbard::C::X_CisAisCjtAku_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, 
                                       org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);

      if (Large::mode != M_CORR) //for hermite: CisAkuCjtAis=-CisAisCjtAku
        mltply::Hubbard::C::X_CisAisCjtAku_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                         org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);    
    }/*if (isite2 != isite3)*/
    return;
  }//MP::myrank =origin
  else {
    idim_max_buf = wrapperMPI::SendRecv_i(origin, Check::idim_max);
    wrapperMPI::SendRecv_iv(origin, Check::idim_max, idim_max_buf, List::c1, List::c1buf);

    wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);
    if (org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite
     && org_isite3 >= Def::Nsite && org_isite4 >= Def::Nsite)
    {
      if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
      else Adiff = isite1 - isite2 * 2;
      if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
      else Bdiff = isite3 - isite4 * 2;

      if (iFlgHermite == FALSE) {
        Fsgn = mltply::Hubbard::GC::X_CisAjt((long int) MP::myrank, isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
        Fsgn *= mltply::Hubbard::GC::X_CisAjt(tmp_off2, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == FALSE)*/
      else {
        Fsgn = mltply::Hubbard::GC::X_CisAjt((long int) MP::myrank, isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
        Fsgn *= mltply::Hubbard::GC::X_CisAjt(tmp_off2, isite1, isite2, (isite1 + isite2), Adiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == TRUE)*/
#pragma omp parallel default(none) private(j,ioff) \
shared(idim_max_buf,tmp_V,Large::irght, Large::ilft, Large::ihfbit, \
Wave::v1buf,tmp_v1,nstate,tmp_v0,List::c2_1,List::c2_2,List::c1buf,one)
      {
#pragma omp for
        for (j = 0; j < idim_max_buf; j++) {
          if (GetOffComp(List::c2_1, List::c2_2, List::c1buf[j],
            Large::irght, Large::ilft, Large::ihfbit, &ioff) == TRUE)
          {
            zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }
        }/*for (j = 0; j < idim_max_buf; j++)*/
      }/*End of parallel region*/
    }//org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite
            // && org_isite3 >= Def::Nsite && org_isite4 >= Def::Nsite
    else {
      org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;

#pragma omp parallel default(none) private(j, dmv, tmp_off, Fsgn, ioff) \
shared(MP::myrank, idim_max_buf, tmp_V, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, \
org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite4, org_ispin4, \
Wave::v1buf, tmp_v1, nstate,one, tmp_v0, List::c1buf, List::c2_1, List::c2_2, Large::irght, Large::ilft, Large::ihfbit)
      {
#pragma omp for
        for (j = 0; j < idim_max_buf; j++) {
          if (mltply::Hubbard::GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, 
            List::c1buf[j] + org_rankbit, &tmp_off) == TRUE)
          {
            if (GetOffComp(List::c2_1, List::c2_2, tmp_off, Large::irght, Large::ilft, Large::ihfbit, &ioff) == TRUE)
            {
              dmv = tmp_V * (std::complex<double>)Fsgn;
              zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
            }
          }
        }/*for (j = 0; j < idim_max_buf; j++)*/
      }/*End of parallel region*/
    }
  }/*if (MP::myrank != origin)*/
}/*std::complex<double> mltply::Hubbard::C::X_CisAjtCkuAlv_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of canonical Hubbard system
*/
void mltply::Hubbard::C::X_CisAjtCkuAku_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  long int idim_max_buf, ioff;
  int iCheck, Fsgn;
  long int isite1, isite2, isite3;
  long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  long int j, Asum, Adiff;
  std::complex<double> dmv;
  long int origin, tmp_off;
  long int org_rankbit;
  int one = 1;
  //printf("Deubg0-0: org_isite1=%d, org_ispin1=%d, org_isite2=%d, org_ispin2=%d, org_isite3=%d, org_ispin3=%d\n", org_isite1, org_ispin1,org_isite2, org_ispin2,org_isite3, org_ispin3);
  iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3, org_ispin3, (long int) MP::myrank, &origin);
  //printf("iCheck=%d, MP::myrank=%d, origin=%d\n", iCheck, MP::myrank, origin);

  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  isite2 = Def::Tpow[2 * org_isite2 + org_ispin2];
  isite3 = Def::Tpow[2 * org_isite3 + org_ispin3];

  if (iCheck == TRUE) {
    tmp_isite1 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
    Asum = tmp_isite1 + tmp_isite2;
    if (tmp_isite2 > tmp_isite1) Adiff = tmp_isite2 - tmp_isite1 * 2;
    else Adiff = tmp_isite1 - tmp_isite2 * 2;
  }/*if (iCheck == TRUE)*/
  else {
    iCheck = mltply::Hubbard::CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1, org_ispin1, (long int) MP::myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = Def::OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = Def::OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = Def::OrgTpow[2 * org_isite3 + org_ispin3];
      Asum = tmp_isite3 + tmp_isite4;
      if (tmp_isite4 > tmp_isite3) Adiff = tmp_isite4 - tmp_isite3 * 2;
      else Adiff = tmp_isite3 - tmp_isite4 * 2;
      if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) tmp_V = 0;
      //printf("tmp_isite1=%ld, tmp_isite2=%ld, Adiff=%ld\n", tmp_isite1, tmp_isite2, Adiff);
    }/*if (iCheck == TRUE)*/
    else return;   
  }/*if (iCheck == FALSE)*/

  if (MP::myrank == origin) {// only k is in PE
    //for hermite
#pragma omp parallel default(none) private(j) \
shared(i_max, Asum, Adiff, isite1, isite2, tmp_V, Large::mode, tmp_v0, tmp_v1,nstate)
    {
#pragma omp for
      for (j = 0; j < i_max; j++)
        mltply::Hubbard::C::CisAjt(j, nstate, tmp_v0, tmp_v1, 
          isite1, isite2, Asum, Adiff, tmp_V);

      if (Large::mode != M_CORR) {
#pragma omp for
        for (j = 0; j < i_max; j++)
          mltply::Hubbard::C::CisAjt(j, nstate, tmp_v0, tmp_v1, 
            isite2, isite1, Asum, Adiff, tmp_V);
      }/*if (Large::mode != M_CORR)*/
    }/*End of parallel region*/
    return;
  }//MP::myrank =origin
  else {
    idim_max_buf = wrapperMPI::SendRecv_i(origin, Check::idim_max);
    wrapperMPI::SendRecv_iv(origin, Check::idim_max, idim_max_buf, List::c1, List::c1buf);
    wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

#pragma omp parallel default(none) private(j,dmv,ioff,tmp_off,Fsgn,Adiff,org_rankbit) \
shared(idim_max_buf,tmp_V,tmp_isite1,tmp_isite2,tmp_isite3,tmp_isite4,isite3,Wave::v1buf,tmp_v1, \
nstate,one,tmp_v0,List::c1buf,List::c2_1,List::c2_2,origin,org_isite3,MP::myrank,isite1,isite2, \
org_isite1,org_isite2,Def::Nsite,Large::irght, Large::ilft, Large::ihfbit,Def::OrgTpow)
    {

      if (org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite) {
        if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
        else Adiff = isite1 - isite2 * 2;
        SgnBit(((long int) MP::myrank & Adiff), &Fsgn);
        tmp_V *= Fsgn;

        if (org_isite3 >= Def::Nsite) {
#pragma omp for
          for (j = 0; j < idim_max_buf; j++) {
            GetOffComp(List::c2_1, List::c2_2, List::c1buf[j],
              Large::irght, Large::ilft, Large::ihfbit, &ioff);
            zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }/*for (j = 0; j < idim_max_buf; j++)*/
        }/*if (org_isite3 >= Def::Nsite)*/
        else { //org_isite3 <= Def::Nsite
#pragma omp for
          for (j = 0; j < idim_max_buf; j++) {
            if (mltply::Hubbard::CheckBit_Ajt(isite3, List::c1buf[j], &tmp_off) == TRUE) {
              GetOffComp(List::c2_1, List::c2_2, List::c1buf[j],
                Large::irght, Large::ilft, Large::ihfbit, &ioff);
              zaxpy_(&nstate, &tmp_V, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
            }
          }/*for (j = 0; j < idim_max_buf; j++)*/
        }/*if (org_isite3 <= Def::Nsite)*/
      }/*if (org_isite1 >= Def::Nsite && org_isite2 >= Def::Nsite)*/
      else {
        org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;
#pragma omp for
        for (j = 0; j < idim_max_buf; j++) {
          if (mltply::Hubbard::GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, 
            List::c1buf[j] + org_rankbit, &tmp_off) == TRUE) {
            dmv = tmp_V * (std::complex<double>)Fsgn;
            GetOffComp(List::c2_1, List::c2_2, tmp_off,
              Large::irght, Large::ilft, Large::ihfbit, &ioff);
            zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }
        }/*for (j = 0; j < idim_max_buf; j++)*/
      }
    }/*End of parallel region*/
  }/*if (MP::myrank != origin)*/
}/*std::complex<double> mltply::Hubbard::C::X_CisAjtCkuAku_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of canonical Hubbard system
*/
void mltply::Hubbard::C::X_CisAisCjtAku_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  mltply::Hubbard::C::X_CisAjtCkuAku_MPI(
    org_isite4, org_ispin4, org_isite3, org_ispin3,
    org_isite1, org_ispin1, conj(tmp_V), nstate, tmp_v0, tmp_v1);
}/*std::complex<double> mltply::Hubbard::C::X_CisAisCjtAku_MPI*/

void mltply::Hubbard::C::X_CisAis_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_V,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  long int j, isite1, tmp_off;
  int one = 1;

  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  if (org_isite1 >= Def::Nsite) {
    if (mltply::Hubbard::CheckBit_Ajt(isite1, (long int) MP::myrank, &tmp_off) == FALSE)
      return;

    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[0][0], &tmp_v0[0][0]);
  }/*if (org_isite1 >= Def::Nsite)*/
  else {
#pragma omp parallel for default(none) private(j, tmp_off) \
shared(tmp_v0, tmp_v1, List::c1,nstate,one, i_max, tmp_V, isite1)
    for (j = 0; j < i_max; j++) {
      if (mltply::Hubbard::C::X_CisAis(List::c1[j], isite1) != 0) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }/*if (X_CisAis(List::c1[j], isite1) != 0)*/
    }/*for (j = 0; j < i_max; j++)*/
  }/*if (org_isite1 <= Def::Nsite)*/
}/*std::complex<double> mltply::Hubbard::C::X_CisAis_MPI*/
/**
@brief Single creation/annihilation operator
 in the inter process region for HubbardGC.
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void mltply::Hubbard::GC::X_Cis_MPI(
  int org_isite,//!<[in] Site i
  int org_ispin,//!<[in] Spin s
  std::complex<double> tmp_trans,//!<[in] Coupling constant//!<[in]
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 += H v1*/,
  std::complex<double> **tmp_v1,//!<[in] v0 += H v1*/,
  long int idim_max//!<[in] Similar to CheckList::idim_max
) {
  int mask2, state2, origin, bit2diff, Fsgn;
  std::complex<double> trans;

  // org_isite >= Nsite
  mask2 = (int)Def::Tpow[2 * org_isite + org_ispin];

  origin = MP::myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in MP::myrank
  //origin: if the state (org_isite, org_ispin) is occupied in MP::myrank, the state is not occupied in origin.

  bit2diff = MP::myrank - ((2 * mask2 - 1) & MP::myrank);

  //SgnBit((long int) (origin & bit2diff), &Fsgn); // Fermion sign
  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign

  wrapperMPI::SendRecv_cv(origin, idim_max*nstate, Check::idim_max*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

  if (state2 == mask2) {
    trans = 0;
  }
  else if (state2 == 0) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

  zaxpy_long(Check::idim_max*nstate, trans, &Wave::v1buf[0][0], &tmp_v0[0][0]);
}/*std::complex<double> mltply::Hubbard::GC::X_Cis_MPI*/
/**
@brief Single creation/annihilation operator
  in the inter process region for HubbardGC.
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void mltply::Hubbard::GC::X_Ajt_MPI(
  int org_isite,//!<[in] Site j
  int org_ispin,//!<[in] Spin t
  std::complex<double> tmp_trans,//!<[in] Coupling constant//!<[in]
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 += H v1*/,
  std::complex<double> **tmp_v1,//!<[in] v0 += H v1*/,
  long int idim_max//!<[in] Similar to CheckList::idim_max
) {
  int mask2, state2, origin, bit2diff, Fsgn;
  std::complex<double> trans;

  // org_isite >= Nsite
  mask2 = (int)Def::Tpow[2 * org_isite + org_ispin];

  origin = MP::myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in MP::myrank
  //origin: if the state (org_isite, org_ispin) is occupied in MP::myrank, the state is not occupied in origin.

  bit2diff = MP::myrank - ((2 * mask2 - 1) & MP::myrank);

  //SgnBit((long int) (origin & bit2diff), &Fsgn); // Fermion sign
  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign

  wrapperMPI::SendRecv_cv(origin, idim_max*nstate, idim_max *nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

  if (     state2 == 0    ) trans = 0;
  else if (state2 == mask2) trans = (double)Fsgn * tmp_trans;
  else return;

  zaxpy_long(idim_max *nstate, trans, &Wave::v1buf[0][0], &tmp_v0[0][0]);
}/*std::complex<double> mltply::Hubbard::GC::X_Ajt_MPI*/
/**
@brief Compute @f$c_{is}^\dagger@f$
term of canonical Hubbard system
*/
void mltply::Hubbard::C::X_Cis_MPI(
  int org_isite,//!<[in] Site i
  int org_ispin,//!<[in] Spin s
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[inout] Initial wavefunction
  long int idim_max//!<[in] Similar to CheckList::idim_max
) {
  int mask2, state2, origin, bit2diff, Fsgn;
  long int idim_max_buf, j, ioff;
  std::complex<double> trans;
  int one = 1;

  // org_isite >= Nsite
  mask2 = (int)Def::Tpow[2 * org_isite + org_ispin];

  origin = MP::myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in MP::myrank
  //origin: if the state (org_isite, org_ispin) is occupied in MP::myrank, the state is not occupied in origin.

  bit2diff = MP::myrank - ((2 * mask2 - 1) & MP::myrank);

  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = wrapperMPI::SendRecv_i(origin, idim_max);
  wrapperMPI::SendRecv_iv(origin, idim_max, idim_max_buf, List::c1_org, List::c1buf_org);
  wrapperMPI::SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

  if (state2 == mask2) {
    trans = 0;
  }
  else if (state2 == 0) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

#pragma omp parallel for default(none) private(j) \
shared(idim_max_buf, trans, ioff, Large::irght, Large::ilft, Large::ihfbit, List::c2_1, List::c2_2, \
Wave::v1buf, tmp_v1, nstate,one, tmp_v0, List::c1buf_org)
  for (j = 0; j < idim_max_buf; j++) {//idim_max_buf -> original
    GetOffComp(List::c2_1, List::c2_2, List::c1buf_org[j],
      Large::irght, Large::ilft, Large::ihfbit, &ioff);
    zaxpy_(&nstate, &trans, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 0; j < idim_max_buf; j++)*/
}/*std::complex<double> mltply::Hubbard::GC::X_Cis_MPI*/
/**
@brief Compute @f$c_{jt}@f$
term of canonical Hubbard system
*/
void mltply::Hubbard::C::X_Ajt_MPI(
  int org_isite,//!<[in] Site j
  int org_ispin,//!<[in] Spin t
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[inout] Initial wavefunction
  long int idim_max//!<[in] Similar to CheckList::idim_max
){
  int mask2, state2, origin, bit2diff, Fsgn;
  long int idim_max_buf, j, ioff;
  std::complex<double> trans;
  int one = 1;

  // org_isite >= Nsite
  mask2 = (int)Def::Tpow[2 * org_isite + org_ispin];

  origin = MP::myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in MP::myrank
  //origin: if the state (org_isite, org_ispin) is occupied in MP::myrank, the state is not occupied in origin.

  bit2diff = MP::myrank - ((2 * mask2 - 1) & MP::myrank);

  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign
  idim_max_buf = wrapperMPI::SendRecv_i(origin, idim_max);
  wrapperMPI::SendRecv_iv(origin, idim_max, idim_max_buf, List::c1_org, List::c1buf_org);
  wrapperMPI::SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

  if (state2 == 0) {
    trans = 0;
  }
  else if (state2 == mask2) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

#pragma omp parallel for default(none) private(j, ioff) \
shared(idim_max_buf, trans, Large::irght, Large::ilft, Large::ihfbit, List::c2_1, List::c2_2, \
Wave::v1buf, tmp_v1, nstate,one, tmp_v0, List::c1buf_org)
  for (j = 0; j < idim_max_buf; j++) {
    GetOffComp(List::c2_1, List::c2_2, List::c1buf_org[j],
      Large::irght, Large::ilft, Large::ihfbit, &ioff);
    zaxpy_(&nstate, &trans, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }
}/*std::complex<double> mltply::Hubbard::C::X_Ajt_MPI*/
