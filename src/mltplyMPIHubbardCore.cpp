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
int CheckPE(
  int org_isite//!<[in] Site index
){
  if (org_isite + 1 > Def::Nsite) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}/*int CheckPE*/
/**
@brief Check the occupation of @f$(i,s)@f$ state,
and compute the index of final wavefunction associated to 
@f$c^\dagger_{is}@f$
@return 1 if unoccupied, 0 if occupied
*/
int CheckBit_Cis(
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
}/*int CheckBit_Cis*/
/**
@brief Check the occupation of @f$(i,s)@f$ state,
and compute the index of final wavefunction associated to
@f$c_{jt}@f$
@return 1 if occupied, 0 if unoccupied
*/
int CheckBit_Ajt(
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
}/*int CheckBit_Ajt*/
/**
@brief Compute the index of final wavefunction associated to
@f$c_{4}^\dagger c_{3}c_{2}^\dagger c_{1}@f$, and
check whether this operator is relevant or not
@return 1 if relevant, 0 if irrelevant
*/
int CheckBit_InterAllPE(
  int org_isite1,//!<[in] Site 1
  int org_isigma1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_isigma2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_isigma3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_isigma4,//!<[in] Spin 4
  //!<[inout]
  long int orgbit,//!<[in] Index of initial wavefunction
  long int *offbit//!<[out] Index of final wavefunction
){
  long int tmp_ispin;
  long int tmp_org, tmp_off;
  int iflgBitExist = TRUE;
  tmp_org=orgbit;
  tmp_off=0;

  if (CheckPE(org_isite1) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite1 + org_isigma1];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (CheckPE(org_isite2) == TRUE ) {
    tmp_ispin = Def::Tpow[2 * org_isite2 + org_isigma2];
    if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (CheckPE(org_isite3) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite3 + org_isigma3];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (CheckPE(org_isite4) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite4 + org_isigma4];
    if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
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
}/*int CheckBit_InterAllPE*/
/**
@brief Check the occupation of both site 1 and site 3
@return 1 if both sites are occupied, 0 if not
*/
int CheckBit_PairPE(
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
  
  if(CheckPE(org_isite1)==TRUE){
    tmp_ispin = Def::Tpow[2 * org_isite1 + org_isigma1];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist=FALSE;
    }
  }
  
  if (CheckPE(org_isite3) == TRUE) {
    tmp_ispin = Def::Tpow[2 * org_isite3 + org_isigma3];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
  }
  
  if(iflgBitExist != TRUE){
    return FALSE;
  }

  return TRUE;
}/*int CheckBit_PairPE*/
/**
@brief Compute the index of final wavefunction associated to
@f$c_{4}^\dagger c_{3}c_{2}^\dagger c_{1}@f$, and
Fermion sign
@return 1 if relevant, 0 if irrelevant
*/
int GetSgnInterAll(
  long int isite1,//!<[in] Site 1
  long int isite2,//!<[in] Site 2
  long int isite3,//!<[in] Site 3
  long int isite4,//!<[in] Site 4
  int *Fsgn,//!<[out] Fermion sign
  //!<[inout]
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

    tmp_sgn=X_GC_CisAjt(orgbit, tmp_ispin1, tmp_ispin2, tmp_ispin1+tmp_ispin2, diffA, &tmp_off);
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
    tmp_sgn *=X_GC_CisAjt(tmp_off, tmp_ispin1, tmp_ispin2, tmp_ispin1+tmp_ispin2, diffA, offbit);
    
    if(tmp_sgn ==0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
  }
  
  *Fsgn =tmp_sgn;
  *offbit = *offbit%Def::OrgTpow[2*Def::Nsite];
  return TRUE;
}/*int GetSgnInterAll*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{jt}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAisCjtAjt_Hubbard_MPI(
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

  iCheck=CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, (long int) myrank);
  if(iCheck != TRUE){
    return;
  }

  if (org_isite1 + 1 > Def::Nsite && org_isite3 + 1 > Def::Nsite) {
    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > Def::Nsite && org_isite3 + 1 > Def::Nsite)*/
  else if (org_isite1 + 1 > Def::Nsite || org_isite3 + 1 > Def::Nsite) {
    if (org_isite1 > org_isite3) tmp_ispin1 = Def::Tpow[2 * org_isite3 + org_ispin3];
    else                         tmp_ispin1 = Def::Tpow[2 * org_isite1 + org_ispin1];

#pragma omp parallel  default(none) \
  shared(org_isite1,org_ispin1,org_isite3,org_ispin3,nstate,one,tmp_v0,tmp_v1,tmp_ispin1) \
  firstprivate(i_max,tmp_V,X) private(j,tmp_off)
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      if (CheckBit_Ajt(tmp_ispin1, j - 1, &tmp_off) == TRUE) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }
    }/*for (j = 1; j <= i_max; j++)*/
  }
}/*std::complex<double> X_GC_child_CisAisCjtAjt_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAjtCkuAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  long int idim_max_buf;
  int iCheck, Fsgn;
  long int isite1, isite2, isite3;
  long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  long int j, Asum, Adiff;
  std::complex<double> dmv;
  long int origin, tmp_off;
  long int org_rankbit;
  int one = 1;

  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3, org_ispin3, (long int) myrank, &origin);
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
    iCheck = CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1, org_ispin1, (long int) myrank, &origin);
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

  if (myrank == origin) {// only k is in PE

    if (CheckBit_Ajt(isite3, myrank, &tmp_off) == FALSE) return;

#pragma omp parallel default(none)  \
firstprivate(i_max,Asum,Adiff,isite1,isite2, tmp_V) \
  private(j,tmp_off) shared(tmp_v0, tmp_v1,nstate)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++) 
        GC_CisAjt(j, nstate, tmp_v0, tmp_v1, isite2, isite1, Asum, Adiff, tmp_V, &tmp_off);

      if (Large::mode != M_CORR) {
#pragma omp for
        for (j = 1; j <= i_max; j++) 
          GC_CisAjt(j, nstate, tmp_v0, tmp_v1, isite1, isite2, Asum, Adiff, tmp_V, &tmp_off);
      }/*if (Large::mode != M_CORR)*/
    }/*End of paralle region*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, Check::idim_max);
    SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

#pragma omp parallel default(none) private(j,dmv,tmp_off,Fsgn,org_rankbit,Adiff) \
  shared(v1buf,tmp_v1,nstate,one,tmp_v0,myrank,origin,isite3,org_isite3,isite1,isite2,org_isite2,org_isite1) \
firstprivate(idim_max_buf,tmp_V,tmp_isite1,tmp_isite2,tmp_isite3,tmp_isite4)
    {
      if (org_isite1 + 1 > Def::Nsite && org_isite2 + 1 > Def::Nsite) {
        if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
        else Adiff = isite1 - isite2 * 2;
        SgnBit(((long int) myrank & Adiff), &Fsgn);
        tmp_V *= Fsgn;

        if (org_isite3 + 1 > Def::Nsite) {
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[j][0], &one);
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }
        else { //org_isite3 <= Def::Nsite
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            if (CheckBit_Ajt(isite3, j - 1, &tmp_off) == TRUE) {
              zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[j][0], &one);
            }
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }
      }/*if (org_isite1 + 1 > Def::Nsite && org_isite2 + 1 > Def::Nsite)*/
      else {
        org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, (j - 1) + org_rankbit, &tmp_off) == TRUE) {
            dmv = tmp_V * (std::complex<double>)Fsgn;
            zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[tmp_off + 1][0], &one);
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }
    }/*End of parallel region*/
  }/*myrank != origin*/
}/*std::complex<double> X_GC_child_CisAjtCkuAku_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAisCjtAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  X_GC_child_CisAjtCkuAku_Hubbard_MPI(
    org_isite4, org_ispin4, org_isite3, org_ispin3,
    org_isite1, org_ispin1, conj(tmp_V), nstate, tmp_v0, tmp_v1);
}/*std::complex<double> X_GC_child_CisAisCjtAku_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAjtCkuAlv_Hubbard_MPI(
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
  long int org_rankbit;
  int iFlgHermite = FALSE;
  int one = 1;

  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                               org_isite3, org_ispin3, org_isite4, org_ispin4,
                               (long int) myrank, &origin);
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
    iCheck = CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                 org_isite2, org_ispin2, org_isite1, org_ispin1,
                                 (long int) myrank, &origin);
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

  if (myrank == origin) {
    if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
      X_GC_child_CisAis_Hubbard_MPI(org_isite1, org_ispin1, tmp_V, nstate, tmp_v0, tmp_v1);
      //calc -nisniv
      X_GC_child_CisAisCjtAjt_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);
    }/*if (isite1 == isite4 && isite2 == isite3)*/
    else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
            //calc CisAku
      if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
      else Adiff = isite1 - isite4 * 2;

#pragma omp parallel for default(none)  private(j, tmp_off) \
  firstprivate(i_max, tmp_V, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0,nstate)
      for (j = 1; j <= i_max; j++) 
        GC_CisAjt(j - 1, nstate, tmp_v0, tmp_v1, isite1, isite4, (isite1 + isite4), Adiff, tmp_V, &tmp_off);
      
      //calc -CisAku njv
      X_GC_child_CisAjtCkuAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4, 
                                          org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);
      if (Large::mode != M_CORR) { //for hermite
#pragma omp parallel for default(none)  private(j, tmp_off) \
  firstprivate(i_max, tmp_V, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0,nstate)
        for (j = 1; j <= i_max; j++) 
          GC_CisAjt(j - 1, nstate, tmp_v0, tmp_v1, isite4, isite1, (isite1 + isite4), Adiff, tmp_V, &tmp_off);
        
        //calc -njvCkuAis
        X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4,
                                            org_isite1, org_ispin1, -tmp_V, nstate, tmp_v0, tmp_v1);
      }/*if (Large::mode != M_CORR)*/
    }/*if (isite2 == isite3)*/
    else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE
      X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3,
                                          org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);
      if (Large::mode != M_CORR) { //for hermite
        X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                            org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);
      }/*if (Large::mode != M_CORR)*/
    }/*if (isite2 != isite3)*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, Check::idim_max);
    SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

    if (org_isite1 + 1 > Def::Nsite && org_isite2 + 1 > Def::Nsite
     && org_isite3 + 1 > Def::Nsite && org_isite4 + 1 > Def::Nsite) {

      if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
      else Adiff = isite1 - isite2 * 2;
      if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
      else Bdiff = isite3 - isite4 * 2;

      if (iFlgHermite == FALSE) {
        Fsgn = X_GC_CisAjt((long int) myrank, isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == FALSE)*/
      else {
        Fsgn = X_GC_CisAjt((long int) myrank, isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, isite1, isite2, (isite1 + isite2), Adiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == TRUE)*/

      zaxpy_long(i_max*nstate, tmp_V, &v1buf[1][0], &tmp_v0[1][0]);
    }
    else {
      org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;
#pragma omp parallel for default(none)  private(j,dmv,tmp_off,Fsgn) \
  firstprivate(idim_max_buf,tmp_V,tmp_isite1,tmp_isite2,tmp_isite3,tmp_isite4,org_rankbit) \
  shared(v1buf,tmp_v1,tmp_v0,nstate,one)
      for (j = 1; j <= idim_max_buf; j++) {
        if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, (j - 1) + org_rankbit, &tmp_off) == TRUE) {
          dmv = tmp_V * (std::complex<double>)Fsgn;
          zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[tmp_off + 1][0], &one);
        }
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
  }/*myrank != origin*/
}/*std::complex<double> X_GC_child_CisAjtCkuAlv_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAis_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  long int j, isite1, tmp_off;
  int one = 1;

  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  if (org_isite1 + 1 > Def::Nsite) {
    if (CheckBit_Ajt(isite1, (long int) myrank, &tmp_off) == FALSE) return;

    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > Def::Nsite)*/
  else {
#pragma omp parallel  default(none) shared(tmp_v0, tmp_v1,nstate,one) \
  firstprivate(i_max, tmp_V, isite1) private(j, tmp_off)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++) {
        if (CheckBit_Ajt(isite1, j - 1, &tmp_off) == TRUE) {
          zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
        }/*if (CheckBit_Ajt(isite1, j - 1, &tmp_off) == TRUE)*/
      }/*for (j = 1; j <= i_max; j++)*/
    }/*End of parallel region*/
  }/*if (org_isite1 + 1 <= Def::Nsite)*/
}/*std::complex<double> X_GC_child_CisAis_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAjt_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {

  if (org_isite1 + 1 > Def::Nsite && org_isite2 + 1 > Def::Nsite) {
    X_GC_child_general_hopp_MPIdouble(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, nstate, tmp_v0, tmp_v1);
  }
  else if (org_isite1 + 1 > Def::Nsite || org_isite2 + 1 > Def::Nsite) {
    X_GC_child_general_hopp_MPIsingle(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, nstate, tmp_v0, tmp_v1);
  }
  else {
    //error message will be added.
    exitMPI(-1);
  }
}/*std::complex<double> X_GC_child_CisAjt_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{jt}@f$
term of canonical Hubbard system
*/
void X_child_CisAisCjtAjt_Hubbard_MPI(
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

  iCheck = CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, (long int) myrank);
  if (iCheck != TRUE) return;
  
  if (org_isite1 + 1 > Def::Nsite && org_isite3 + 1 > Def::Nsite) {
    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > Def::Nsite && org_isite3 + 1 > Def::Nsite)*/
  else if (org_isite1 + 1 > Def::Nsite || org_isite3 + 1 > Def::Nsite) {
    if (org_isite1 > org_isite3) tmp_ispin1 = Def::Tpow[2 * org_isite3 + org_ispin3];
    else                         tmp_ispin1 = Def::Tpow[2 * org_isite1 + org_ispin1];

#pragma omp parallel for default(none) \
shared(tmp_v0,tmp_v1,list_1,org_isite1,org_ispin1,org_isite3,org_ispin3,nstate,one) \
  firstprivate(i_max,tmp_V,tmp_ispin1) private(j,tmp_off)
    for (j = 1; j <= i_max; j++) {
      if (CheckBit_Ajt(tmp_ispin1, list_1[j], &tmp_off) == TRUE) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }
    }/*for (j = 1; j <= i_max; j++)*/
  }/*if (org_isite1 + 1 > Def::Nsite || org_isite3 + 1 > Def::Nsite)*/
}/*std::complex<double> X_child_CisAisCjtAjt_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of canonical Hubbard system
*/
void X_child_CisAjtCkuAlv_Hubbard_MPI(
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

  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                               org_isite3, org_ispin3, org_isite4, org_ispin4,
                               (long int) myrank, &origin);
  //printf("iCheck=%d, myrank=%d, origin=%d\n", iCheck, myrank, origin);
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
    iCheck = CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                 org_isite2, org_ispin2, org_isite1, org_ispin1,
                                 (long int) myrank, &origin);
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

  if (myrank == origin) {
    if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
      X_child_CisAis_Hubbard_MPI(org_isite1, org_ispin1, tmp_V, nstate, tmp_v0, tmp_v1);
      //calc -nisniv
      X_child_CisAisCjtAjt_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);
    }/* if (isite1 == isite4 && isite2 == isite3)*/
    else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
      if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
      else Adiff = isite1 - isite4 * 2;

      //calc CisAku
#pragma omp parallel for default(none)  private(j, tmp_off) \
firstprivate(i_max, tmp_V, isite1, isite4, Adiff) shared(tmp_v1, nstate, tmp_v0, list_1)
      for (j = 1; j <= i_max; j++)
        CisAjt(j, nstate, tmp_v0, tmp_v1, isite1, isite4, (isite1 + isite4), Adiff, tmp_V);
      
      //calc -CisAku njv
      X_child_CisAjtCkuAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4,
                                       org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);

      if (Large::mode != M_CORR) {  //for hermite
#pragma omp parallel for default(none)  private(j, tmp_off) \
  firstprivate(i_max, tmp_V, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0,nstate)
        for (j = 1; j <= i_max; j++) 
          CisAjt(j, nstate, tmp_v0, tmp_v1, isite4, isite1, (isite1 + isite4), Adiff, tmp_V);
                
        //calc -njvCkuAis
        X_child_CisAisCjtAku_Hubbard_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4, 
                                         org_isite1, org_ispin1, -tmp_V, nstate, tmp_v0, tmp_v1);
      }/*if (Large::mode != M_CORR)*/
    }/*if (isite2 == isite3)*/
    else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE
      X_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, 
                                       org_isite2, org_ispin2, -tmp_V, nstate, tmp_v0, tmp_v1);

      if (Large::mode != M_CORR) //for hermite: CisAkuCjtAis=-CisAisCjtAku
        X_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                         org_isite3, org_ispin3, -tmp_V, nstate, tmp_v0, tmp_v1);    
    }/*if (isite2 != isite3)*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, Check::idim_max);
    SendRecv_iv(origin, Check::idim_max + 1, idim_max_buf + 1, list_1, list_1buf);

    SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);
    if (org_isite1 + 1 > Def::Nsite && org_isite2 + 1 > Def::Nsite
     && org_isite3 + 1 > Def::Nsite && org_isite4 + 1 > Def::Nsite)
    {
      if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
      else Adiff = isite1 - isite2 * 2;
      if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
      else Bdiff = isite3 - isite4 * 2;

      if (iFlgHermite == FALSE) {
        Fsgn = X_GC_CisAjt((long int) myrank, isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == FALSE)*/
      else {
        Fsgn = X_GC_CisAjt((long int) myrank, isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, isite1, isite2, (isite1 + isite2), Adiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == TRUE)*/
#pragma omp parallel default(none) private(j,ioff) \
firstprivate(idim_max_buf,tmp_V,X) \
  shared(v1buf,tmp_v1,nstate,tmp_v0,list_2_1,list_2_2,list_1buf,one)
      {
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetOffComp(list_2_1, list_2_2, list_1buf[j],
            Large::irght, Large::ilft, Large::ihfbit, &ioff) == TRUE)
          {
            zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }/*End of parallel region*/
    }//org_isite1+1 > Def::Nsite && org_isite2+1 > Def::Nsite
            // && org_isite3+1 > Def::Nsite && org_isite4+1 > Def::Nsite
    else {
      org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;

#pragma omp parallel default(none)  private(j, dmv, tmp_off, Fsgn, ioff) \
firstprivate(myrank, idim_max_buf, tmp_V, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, \
org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite4, org_ispin4) \
  shared(v1buf, tmp_v1, nstate,one, tmp_v0, list_1buf, list_2_1, list_2_2)
      {
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, 
            list_1buf[j] + org_rankbit, &tmp_off) == TRUE)
          {
            if (GetOffComp(list_2_1, list_2_2, tmp_off, Large::irght, Large::ilft, Large::ihfbit, &ioff) == TRUE)
            {
              dmv = tmp_V * (std::complex<double>)Fsgn;
              zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
            }
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }/*End of parallel region*/
    }
  }/*if (myrank != origin)*/
}/*std::complex<double> X_child_CisAjtCkuAlv_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of canonical Hubbard system
*/
void X_child_CisAjtCkuAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
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
  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3, org_ispin3, (long int) myrank, &origin);
  //printf("iCheck=%d, myrank=%d, origin=%d\n", iCheck, myrank, origin);

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
    iCheck = CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1, org_ispin1, (long int) myrank, &origin);
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

  if (myrank == origin) {// only k is in PE
    //for hermite
#pragma omp parallel default(none)  \
firstprivate(i_max, Asum, Adiff, isite1, isite2, tmp_V, X) \
  private(j) shared(tmp_v0, tmp_v1,nstate)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        CisAjt(j, nstate, tmp_v0, tmp_v1, isite1, isite2, Asum, Adiff, tmp_V);

      if (Large::mode != M_CORR) {
#pragma omp for
        for (j = 1; j <= i_max; j++)
          CisAjt(j, nstate, tmp_v0, tmp_v1, isite2, isite1, Asum, Adiff, tmp_V);
      }/*if (Large::mode != M_CORR)*/
    }/*End of parallel region*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, Check::idim_max);
    SendRecv_iv(origin, Check::idim_max + 1, idim_max_buf + 1, list_1, list_1buf);
    SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

#pragma omp parallel default(none)  private(j,dmv,ioff,tmp_off,Fsgn,Adiff,org_rankbit) \
firstprivate(idim_max_buf,tmp_V,tmp_isite1,tmp_isite2,tmp_isite3,tmp_isite4,isite3) \
  shared(v1buf,tmp_v1,nstate,one,tmp_v0,list_1buf,list_2_1,list_2_2,origin,org_isite3,myrank,isite1,isite2,org_isite1,org_isite2)
    {

      if (org_isite1 + 1 > Def::Nsite && org_isite2 + 1 > Def::Nsite) {
        if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
        else Adiff = isite1 - isite2 * 2;
        SgnBit(((long int) myrank & Adiff), &Fsgn);
        tmp_V *= Fsgn;

        if (org_isite3 + 1 > Def::Nsite) {
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            GetOffComp(list_2_1, list_2_2, list_1buf[j],
              Large::irght, Large::ilft, Large::ihfbit, &ioff);
            zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }/*if (org_isite3 + 1 > Def::Nsite)*/
        else { //org_isite3 <= Def::Nsite
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            if (CheckBit_Ajt(isite3, list_1buf[j], &tmp_off) == TRUE) {
              GetOffComp(list_2_1, list_2_2, list_1buf[j],
                Large::irght, Large::ilft, Large::ihfbit, &ioff);
              zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
            }
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }/*if (org_isite3 + 1 <= Def::Nsite)*/
      }/*if (org_isite1 + 1 > Def::Nsite && org_isite2 + 1 > Def::Nsite)*/
      else {
        org_rankbit = Def::OrgTpow[2 * Def::Nsite] * origin;
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, 
            list_1buf[j] + org_rankbit, &tmp_off) == TRUE) {
            dmv = tmp_V * (std::complex<double>)Fsgn;
            GetOffComp(list_2_1, list_2_2, tmp_off,
              Large::irght, Large::ilft, Large::ihfbit, &ioff);
            zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }
    }/*End of parallel region*/
  }/*if (myrank != origin)*/
}/*std::complex<double> X_child_CisAjtCkuAku_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of canonical Hubbard system
*/
void X_child_CisAisCjtAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  X_child_CisAjtCkuAku_Hubbard_MPI(
    org_isite4, org_ispin4, org_isite3, org_ispin3,
    org_isite1, org_ispin1, conj(tmp_V), nstate, tmp_v0, tmp_v1);
}/*std::complex<double> X_child_CisAisCjtAku_Hubbard_MPI*/

void X_child_CisAis_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  std::complex<double> tmp_V,//!<[in] Coupling constant
  //!<[inout]
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1//!<[inout] Initial wavefunction
) {
  long int i_max = Check::idim_max;
  long int j, isite1, tmp_off;
  int one = 1;

  isite1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  if (org_isite1 + 1 > Def::Nsite) {
    if (CheckBit_Ajt(isite1, (long int) myrank, &tmp_off) == FALSE)
      return;

    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > Def::Nsite)*/
  else {
#pragma omp parallel  default(none) shared(tmp_v0, tmp_v1, list_1,nstate,one) \
  firstprivate(i_max, tmp_V, isite1) private(j, tmp_off)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++) {
        if (X_CisAis(list_1[j], isite1) != 0) {
          zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
        }/*if (X_CisAis(list_1[j], isite1) != 0)*/
      }/*for (j = 1; j <= i_max; j++)*/
    }/*End of parallel region*/
  }/*if (org_isite1 + 1 <= Def::Nsite)*/
}/*std::complex<double> X_child_CisAis_Hubbard_MPI*/
/**
@brief Single creation/annihilation operator
 in the inter process region for HubbardGC.
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void X_GC_Cis_MPI(
  int org_isite,//!<[in] Site i
  int org_ispin,//!<[in] Spin s
  std::complex<double> tmp_trans,//!<[in] Coupling constant//!<[in]
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 += H v1*/,
  std::complex<double> **tmp_v1,//!<[in] v0 += H v1*/,
  long int idim_max,//!<[in] Similar to CheckList::idim_max
  long int *Tpow//!<[in] Similar to DefineList::Tpow
) {
  int mask2, state2, origin, bit2diff, Fsgn;
  long int idim_max_buf;
  std::complex<double> trans;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  //SgnBit((long int) (origin & bit2diff), &Fsgn); // Fermion sign
  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

  if (state2 == mask2) {
    trans = 0;
  }
  else if (state2 == 0) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

  zaxpy_long(idim_max_buf*nstate, trans, &v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_Cis_MPI*/
/**
@brief Single creation/annihilation operator
  in the inter process region for HubbardGC.
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void X_GC_Ajt_MPI(
  int org_isite,//!<[in] Site j
  int org_ispin,//!<[in] Spin t
  std::complex<double> tmp_trans,//!<[in] Coupling constant//!<[in]
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 += H v1*/,
  std::complex<double> **tmp_v1,//!<[in] v0 += H v1*/,
  long int idim_max,//!<[in] Similar to CheckList::idim_max
  long int *Tpow//!<[in] Similar to DefineList::Tpow
) {
  int mask2, state2, origin, bit2diff, Fsgn;
  long int idim_max_buf;
  std::complex<double> trans;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  //SgnBit((long int) (origin & bit2diff), &Fsgn); // Fermion sign
  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

  if (     state2 == 0    ) trans = 0;
  else if (state2 == mask2) trans = (double)Fsgn * tmp_trans;
  else return;

  zaxpy_long(idim_max_buf*nstate, trans, &v1buf[1][0], &tmp_v0[1][0]);
}/*std::complex<double> X_GC_Ajt_MPI*/
/**
@brief Compute @f$c_{is}^\dagger@f$
term of canonical Hubbard system
*/
void X_Cis_MPI(
  int org_isite,//!<[in] Site i
  int org_ispin,//!<[in] Spin s
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[inout] Initial wavefunction
  long int idim_max,//!<[in] Similar to CheckList::idim_max
  long int *Tpow,//!<[in] Similar to DefineList::Tpow
  long int _irght,//!<[in] Similer to LargeList::irght
  long int _ilft,//!<[in] Similer to LargeList::ilft
  long int _ihfbit//!<[in] Similer to LargeList::ihfbit
) {
  int mask2, state2, origin, bit2diff, Fsgn;
  long int idim_max_buf, j, ioff;
  std::complex<double> trans;
  int one = 1;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_iv(origin, idim_max + 1, idim_max_buf + 1, list_1_org, list_1buf_org);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

  if (state2 == mask2) {
    trans = 0;
  }
  else if (state2 == 0) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

#pragma omp parallel for default(none) private(j) \
firstprivate(idim_max_buf, trans, ioff, _irght, _ilft, _ihfbit, list_2_1, list_2_2) \
  shared(v1buf, tmp_v1, nstate,one, tmp_v0, list_1buf_org)
  for (j = 1; j <= idim_max_buf; j++) {//idim_max_buf -> original
    GetOffComp(list_2_1, list_2_2, list_1buf_org[j],
      _irght, _ilft, _ihfbit, &ioff);
    zaxpy_(&nstate, &trans, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 1; j <= idim_max_buf; j++)*/
}/*std::complex<double> X_GC_Cis_MPI*/
/**
@brief Compute @f$c_{jt}@f$
term of canonical Hubbard system
*/
void X_Ajt_MPI(
  int org_isite,//!<[in] Site j
  int org_ispin,//!<[in] Spin t
  std::complex<double> tmp_trans,//!<[in] Coupling constant
  int nstate, std::complex<double> **tmp_v0,//!<[inout] Resulting wavefunction
  std::complex<double> **tmp_v1,//!<[inout] Initial wavefunction
  long int idim_max,//!<[in] Similar to CheckList::idim_max
  long int *Tpow,//!<[in] Similar to DefineList::Tpow
  long int _irght,//!<[in] Similer to LargeList::irght
  long int _ilft,//!<[in] Similer to LargeList::ilft
  long int _ihfbit//!<[in] Similer to LargeList::ihfbit
){
  int mask2, state2, origin, bit2diff, Fsgn;
  long int idim_max_buf, j, ioff;
  std::complex<double> trans;
  int one = 1;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  SgnBit((long int) (bit2diff), &Fsgn); // Fermion sign
  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_iv(origin, idim_max + 1, idim_max_buf + 1, list_1_org, list_1buf_org);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

  if (state2 == 0) {
    trans = 0;
  }
  else if (state2 == mask2) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

#pragma omp parallel for default(none) private(j) \
firstprivate(idim_max_buf, trans, ioff, _irght, _ilft, _ihfbit, list_2_1, list_2_2) \
  shared(v1buf, tmp_v1, nstate,one, tmp_v0, list_1buf_org)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, list_1buf_org[j],
      _irght, _ilft, _ihfbit, &ioff);
    zaxpy_(&nstate, &trans, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }
}/*std::complex<double> X_Ajt_MPI*/
