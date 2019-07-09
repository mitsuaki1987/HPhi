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
#include "SingleExHubbard.hpp"
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "common/setmemory.hpp"
#include "mltplyHubbardCore.hpp"
#include "mltplyMPIHubbardCore.hpp"
#include "mltplyCommon.hpp"
#include "global.hpp"
#ifdef __MPI
#include "mpi.h"
#endif
/**@file
@brief Functions to compute singly excited state in Hubbard model
*/
/**
@brief Calculation of Single excited state for Hubbard canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Single::Hubbard(
  //!<define list to get and put information of calculation
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1,//!<[in] v0 = H v1
  int iEx
) {
  long int idim_max;
  long int i, j;
  long int org_isite, ispin, itype;
  long int is1_spin;
  int isgn = 1, one = 1;
  std::complex<double> tmpphi, dmv;
  long int tmp_off = 0;
  //tmp_v0
  if (Def::NSingleExcitationOperator[iEx] == 0) {
    return TRUE;
  }

  idim_max = Check::idim_maxOrg;
  for (i = 0; i < Def::NSingleExcitationOperator[iEx]; i++) {
    org_isite = Def::SingleExcitationOperator[iEx][i][0];
    ispin = Def::SingleExcitationOperator[iEx][i][1];
    itype = Def::SingleExcitationOperator[iEx][i][2];
    tmpphi = Def::ParaSingleExcitationOperator[iEx][i];
    is1_spin = Def::Tpow[2 * org_isite + ispin];
    if (itype == 1) {
      if (org_isite >= Def::Nsite) {
        X_Cis_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, idim_max, 
          Def::Tpow, 
          Large::irght, Large::ilft, Large::ihfbit);
      }
      else {
#pragma omp parallel for default(none) private(j,  isgn,tmp_off,dmv) \
shared(Large::irght, Large::ilft, Large::ihfbit,nstate,tmp_v0, tmp_v1, List::c1_org,one, \
idim_max, tmpphi, org_isite, ispin, List::c2_1, List::c2_2, is1_spin)
        for (j = 1; j <= idim_max; j++) {//idim_max -> original dimension
          isgn = X_Cis(j, is1_spin, &tmp_off, List::c1_org, List::c2_1, List::c2_2,
            Large::irght, Large::ilft, Large::ihfbit);
          dmv = (std::complex<double>)isgn * tmpphi;
          zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
        }
      }
    }
    else if (itype == 0) {
      if (org_isite >= Def::Nsite) {
        X_Ajt_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, 
          idim_max, Def::Tpow, Large::irght, Large::ilft, Large::ihfbit);
      }
      else {
#pragma omp parallel for default(none) private(j, isgn, tmp_off,dmv) \
shared(tmp_v0,tmp_v1,List::c1_org,List::c1,one,nstate,idim_max,tmpphi,org_isite,ispin, \
List::c2_1,List::c2_2,is1_spin,MP::myrank, Large::irght, Large::ilft, Large::ihfbit)
        for (j = 1; j <= idim_max; j++) {//idim_max -> original dimension
          isgn = X_Ajt(j, is1_spin, &tmp_off, List::c1_org, List::c2_1, List::c2_2,
            Large::irght, Large::ilft, Large::ihfbit);
          dmv = (std::complex<double>)isgn * tmpphi;
          zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
        }
      }
    }
  }
  return TRUE;
}/*int GetSingleExcitedStateHubbard*/
/**
@brief Calculation of Single excited state for Hubbard Grand canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Single::HubbardGC(
  //!<define list to get and put information of calculation
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1,//!<[in] v0 = H v1
  int iEx
) {
  long int idim_max;
  long int i, j;
  long int org_isite, ispin, itype;
  long int is1_spin;
  std::complex<double> tmpphi;
  long int tmp_off = 0;
  //idim_max = Check::idim_max;
  idim_max = Check::idim_maxOrg;
  //tmp_v0
  if (Def::NSingleExcitationOperator[iEx] == 0) {
    return TRUE;
  }

  // SingleEx
  for (i = 0; i < Def::NSingleExcitationOperator[iEx]; i++) {
    org_isite = Def::SingleExcitationOperator[iEx][i][0];
    ispin = Def::SingleExcitationOperator[iEx][i][1];
    itype = Def::SingleExcitationOperator[iEx][i][2];
    tmpphi = Def::ParaSingleExcitationOperator[iEx][i];
    if (itype == 1) {
      if (org_isite >= Def::Nsite) {
        X_GC_Cis_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, 
          idim_max, Def::Tpow);
      }
      else {
#pragma omp parallel for default(none) private(j, is1_spin, tmp_off) \
shared(tmp_v0,tmp_v1,nstate,idim_max, tmpphi, org_isite, ispin,Def::Tpow)
        for (j = 1; j <= idim_max; j++) {
          is1_spin = Def::Tpow[2 * org_isite + ispin];
          GC_Cis(j, nstate, tmp_v0, tmp_v1, is1_spin, tmpphi, &tmp_off);
        }/*for (j = 1; j <= idim_max; j++)*/
      }
    }
    else if (itype == 0) {
      if (org_isite >= Def::Nsite) {
        X_GC_Ajt_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, 
          idim_max, Def::Tpow);
      }
      else {
#pragma omp parallel for default(none) private(j, is1_spin, tmp_off) \
shared(tmp_v0,tmp_v1,nstate, idim_max, tmpphi, org_isite, ispin,Def::Tpow)
        for (j = 1; j <= idim_max; j++) {
          is1_spin = Def::Tpow[2 * org_isite + ispin];
          GC_Ajt(j, nstate, tmp_v0, tmp_v1, is1_spin, tmpphi, &tmp_off);
        }/*for (j = 1; j <= idim_max; j++)*/
      }
    }
  }
  return TRUE;
}/*int GetSingleExcitedStateHubbardGC*/
