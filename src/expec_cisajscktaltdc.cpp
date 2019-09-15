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

#include "expec_cisajscktaltdc.hpp"
#include "mltply.hpp"
#include "mltplyCommon.hpp"
#include "FileIO.hpp"
#include "bitcalc.hpp"
#include "expec_cisajscktaltdc.hpp"
#include "mltplySpinCore.hpp"
#include "mltplyHubbardCore.hpp"
#include "wrapperMPI.hpp"
#include "mltplyMPISpin.hpp"
#include "mltplyMPISpinCore.hpp"
#include "mltplyMPIHubbardCore.hpp"
#include "common/setmemory.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include "log.hpp"
/**
 * @file
 * 
 * @brief  File for calculating two-body green's functions
 * 
 * @version 0.2
 * @details add function to treat the case of general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
/**
@brief Rearray interactions
@return 0 normally finished
@return -1 unnormally finished
*/
int Rearray_Interactions(
  int i,/**<[in]*/
  long int* org_isite1,/**<[in] site number on the site 1*/
  long int* org_isite2,/**<[in] site number on the site 2*/
  long int* org_isite3,/**<[in] site number on the site 3*/
  long int* org_isite4,/**<[in] site number on the site 4*/
  long int* org_sigma1,/**<[in] spin index on the site 1*/
  long int* org_sigma2,/**<[in] spin index on the site 2*/
  long int* org_sigma3,/**<[in] spin index on the site 3*/
  long int* org_sigma4,/**<[in] spin index on the site 4*/
  std::complex<double>* tmp_V/**<[in] value of interaction*/
)
{
  long int tmp_org_isite1, tmp_org_isite2, tmp_org_isite3, tmp_org_isite4;
  long int tmp_org_sigma1, tmp_org_sigma2, tmp_org_sigma3, tmp_org_sigma4;

  tmp_org_isite1 = Def::CisAjtCkuAlvDC[i][0];
  tmp_org_sigma1 = Def::CisAjtCkuAlvDC[i][1];
  tmp_org_isite2 = Def::CisAjtCkuAlvDC[i][2];
  tmp_org_sigma2 = Def::CisAjtCkuAlvDC[i][3];
  tmp_org_isite3 = Def::CisAjtCkuAlvDC[i][4];
  tmp_org_sigma3 = Def::CisAjtCkuAlvDC[i][5];
  tmp_org_isite4 = Def::CisAjtCkuAlvDC[i][6];
  tmp_org_sigma4 = Def::CisAjtCkuAlvDC[i][7];

  if (tmp_org_isite1 == tmp_org_isite2 && tmp_org_isite3 == tmp_org_isite4) {
    if (tmp_org_isite1 > tmp_org_isite3) {
      *org_isite1 = tmp_org_isite3;
      *org_sigma1 = tmp_org_sigma3;
      *org_isite2 = tmp_org_isite4;
      *org_sigma2 = tmp_org_sigma4;
      *org_isite3 = tmp_org_isite1;
      *org_sigma3 = tmp_org_sigma1;
      *org_isite4 = tmp_org_isite2;
      *org_sigma4 = tmp_org_sigma2;
    }
    else {
      *org_isite1 = tmp_org_isite1;
      *org_sigma1 = tmp_org_sigma1;
      *org_isite2 = tmp_org_isite2;
      *org_sigma2 = tmp_org_sigma2;
      *org_isite3 = tmp_org_isite3;
      *org_sigma3 = tmp_org_sigma3;
      *org_isite4 = tmp_org_isite4;
      *org_sigma4 = tmp_org_sigma4;
    }
    *tmp_V = 1.0;

  }
  else if (tmp_org_isite1 == tmp_org_isite4 && tmp_org_isite3 == tmp_org_isite2) {
    if (tmp_org_isite1 > tmp_org_isite3) {
      *org_isite1 = tmp_org_isite3;
      *org_sigma1 = tmp_org_sigma3;
      *org_isite2 = tmp_org_isite2;
      *org_sigma2 = tmp_org_sigma2;
      *org_isite3 = tmp_org_isite1;
      *org_sigma3 = tmp_org_sigma1;
      *org_isite4 = tmp_org_isite4;
      *org_sigma4 = tmp_org_sigma4;
    }
    else {
      *org_isite1 = tmp_org_isite1;
      *org_sigma1 = tmp_org_sigma1;
      *org_isite2 = tmp_org_isite4;
      *org_sigma2 = tmp_org_sigma4;
      *org_isite3 = tmp_org_isite3;
      *org_sigma3 = tmp_org_sigma3;
      *org_isite4 = tmp_org_isite2;
      *org_sigma4 = tmp_org_sigma2;
    }
    *tmp_V = -1.0;
  }
  else {
    return -1;
  }
  return 0;
}
/**
 * @brief Child function to calculate two-body green's functions for Hubbard GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
void expec::cisajscktalt::Hubbard::GC(
  int nstate, 
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
  std::complex<double> **prod
) {
  long int i, j;
  long int isite1, isite2, isite3, isite4;
  long int org_isite1, org_isite2, org_isite3, org_isite4;
  long int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long int Asum, Bsum, Adiff, Bdiff;
  long int tmp_off = 0;
  std::complex<double> tmp_V = 1.0;
  long int i_max;

  for (i = 0; i < Def::NCisAjtCkuAlvDC; i++) {
    zclear(Large::i_max*nstate, &Xvec[0][0]);
    org_isite1 = Def::CisAjtCkuAlvDC[i][0];
    org_sigma1 = Def::CisAjtCkuAlvDC[i][1];
    org_isite2 = Def::CisAjtCkuAlvDC[i][2];
    org_sigma2 = Def::CisAjtCkuAlvDC[i][3];
    org_isite3 = Def::CisAjtCkuAlvDC[i][4];
    org_sigma3 = Def::CisAjtCkuAlvDC[i][5];
    org_isite4 = Def::CisAjtCkuAlvDC[i][6];
    org_sigma4 = Def::CisAjtCkuAlvDC[i][7];

    if (mltply::Hubbard::CheckPE(org_isite1) == TRUE || mltply::Hubbard::CheckPE(org_isite2) == TRUE ||
      mltply::Hubbard::CheckPE(org_isite3) == TRUE || mltply::Hubbard::CheckPE(org_isite4) == TRUE) {
      isite1 = Def::OrgTpow[2 * org_isite1 + org_sigma1];
      isite2 = Def::OrgTpow[2 * org_isite2 + org_sigma2];
      isite3 = Def::OrgTpow[2 * org_isite3 + org_sigma3];
      isite4 = Def::OrgTpow[2 * org_isite4 + org_sigma4];
      if (isite1 == isite2 && isite3 == isite4) {
        mltply::Hubbard::GC::X_CisAisCjtAjt_MPI(org_isite1, org_sigma1, org_isite3, org_sigma3,
          1.0, nstate, Xvec, vec);
      }
      else if (isite1 == isite2 && isite3 != isite4) {
        mltply::Hubbard::GC::X_CisAisCjtAku_MPI(
          org_isite1, org_sigma1, org_isite3, org_sigma3, org_isite4, org_sigma4,
          1.0, nstate, Xvec, vec);
      }
      else if (isite1 != isite2 && isite3 == isite4) {
        mltply::Hubbard::GC::X_CisAjtCkuAku_MPI(org_isite1, org_sigma1, org_isite2, org_sigma2,
          org_isite3, org_sigma3, 1.0, nstate, Xvec, vec);
      }
      else if (isite1 != isite2 && isite3 != isite4) {
        mltply::Hubbard::GC::X_CisAjtCkuAlv_MPI(org_isite1, org_sigma1, org_isite2, org_sigma2,
          org_isite3, org_sigma3, org_isite4, org_sigma4, 1.0, nstate, Xvec, vec);
      }
    }//InterPE
    else {
      mltply::Hubbard::general_int_GetInfo(org_isite1, org_isite2, org_isite3, org_isite4,
        org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_V);

      i_max = Large::i_max;
      isite1 = Large::is1_spin;
      isite2 = Large::is2_spin;
      Asum = Large::isA_spin;
      Adiff = Large::A_spin;

      isite3 = Large::is3_spin;
      isite4 = Large::is4_spin;
      Bsum = Large::isB_spin;
      Bdiff = Large::B_spin;

      if (isite1 == isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate, \
i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_V)
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::GC::CisAisCisAis_element(j, isite1, isite3, tmp_V, nstate, Xvec, vec);
        }
      }
      else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,Xvec,nstate,i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_V)
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::GC::CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff,
            tmp_V, nstate, Xvec, vec, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,Xvec,nstate,i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_V) 
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::GC::CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff,
            tmp_V, nstate, Xvec, vec, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,Xvec,nstate,i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_V) 
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::GC::CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff,
            Bsum, Bdiff, tmp_V, nstate, Xvec, vec, &tmp_off);
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }//Intra PE
}
/**
 * @brief Child function to calculate two-body green's functions for Hubbard model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
void expec::cisajscktalt::Hubbard::C(
  int nstate, 
  std::complex<double> **Xvec, 
  std::complex<double> **vec,
  std::complex<double> **prod
){
  long int i, j;
  long int isite1, isite2, isite3, isite4;
  long int org_isite1, org_isite2, org_isite3, org_isite4;
  long int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long int Asum, Bsum, Adiff, Bdiff;
  long int tmp_off = 0;
  std::complex<double> tmp_V;
  long int i_max;

  for (i = 0; i < Def::NCisAjtCkuAlvDC; i++) {
    zclear(Large::i_max*nstate, &Xvec[0][0]);
    org_isite1 = Def::CisAjtCkuAlvDC[i][0];
    org_sigma1 = Def::CisAjtCkuAlvDC[i][1];
    org_isite2 = Def::CisAjtCkuAlvDC[i][2];
    org_sigma2 = Def::CisAjtCkuAlvDC[i][3];
    org_isite3 = Def::CisAjtCkuAlvDC[i][4];
    org_sigma3 = Def::CisAjtCkuAlvDC[i][5];
    org_isite4 = Def::CisAjtCkuAlvDC[i][6];
    org_sigma4 = Def::CisAjtCkuAlvDC[i][7];
    tmp_V = 1.0;

    if (Def::iFlgSzConserved == TRUE) {
      if (org_sigma1 + org_sigma3 != org_sigma2 + org_sigma4) {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (mltply::Hubbard::CheckPE(org_isite1) == TRUE || mltply::Hubbard::CheckPE(org_isite2) == TRUE ||
      mltply::Hubbard::CheckPE(org_isite3) == TRUE || mltply::Hubbard::CheckPE(org_isite4) == TRUE) {
      isite1 = Def::OrgTpow[2 * org_isite1 + org_sigma1];
      isite2 = Def::OrgTpow[2 * org_isite2 + org_sigma2];
      isite3 = Def::OrgTpow[2 * org_isite3 + org_sigma3];
      isite4 = Def::OrgTpow[2 * org_isite4 + org_sigma4];
      if (isite1 == isite2 && isite3 == isite4) {
        mltply::Hubbard::C::X_CisAisCjtAjt_MPI(org_isite1, org_sigma1,
          org_isite3, org_sigma3, 1.0, nstate, Xvec, vec);
      }
      else if (isite1 == isite2 && isite3 != isite4) {
        //printf("org_isite1=%d, org_isite2=%d, org_isite3=%d, org_isite4=%d\n", org_isite1, org_isite2, org_isite3, org_isite4);
        mltply::Hubbard::C::X_CisAisCjtAku_MPI(org_isite1, org_sigma1,
          org_isite3, org_sigma3, org_isite4, org_sigma4, 1.0, nstate, Xvec, vec);
      }
      else if (isite1 != isite2 && isite3 == isite4) {
        mltply::Hubbard::C::X_CisAjtCkuAku_MPI(org_isite1, org_sigma1, org_isite2, org_sigma2,
          org_isite3, org_sigma3, 1.0, nstate, Xvec, vec);

      }
      else if (isite1 != isite2 && isite3 != isite4) {
        mltply::Hubbard::C::X_CisAjtCkuAlv_MPI(org_isite1, org_sigma1, org_isite2, org_sigma2,
          org_isite3, org_sigma3, org_isite4, org_sigma4, 1.0, nstate, Xvec, vec);
      }
    }//InterPE
    else {
      mltply::Hubbard::general_int_GetInfo(
        org_isite1, org_isite2, org_isite3, org_isite4,
        org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_V
      );

      i_max = Large::i_max;
      isite1 = Large::is1_spin;
      isite2 = Large::is2_spin;
      Asum = Large::isA_spin;
      Adiff = Large::A_spin;

      isite3 = Large::is3_spin;
      isite4 = Large::is4_spin;
      Bsum = Large::isB_spin;
      Bdiff = Large::B_spin;

      tmp_V = 1.0;
      if (isite1 == isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V,Xvec,nstate, \
i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff)
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::C::CisAisCisAis_element(j, isite1, isite3, tmp_V, nstate, Xvec, vec);
        }
      }
      else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,tmp_V,Xvec,nstate,i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff)
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::C::CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff,
            tmp_V, nstate, Xvec, vec, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,tmp_V,Xvec,nstate,i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff)
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::C::CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff,
            tmp_V, nstate, Xvec, vec, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,tmp_V,Xvec,nstate,i_max,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff)
        for (j = 0; j < i_max; j++) {
          mltply::Hubbard::C::CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff,
            tmp_V, nstate, Xvec, vec, &tmp_off);
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief Child function to calculate two-body green's functions for 1/2 Spin model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
void expec::cisajscktalt::Spin::C::Half(
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec,
  std::complex<double> **prod
){
  long int i, j;
  long int org_isite1, org_isite2, org_isite3, org_isite4;
  long int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long int isA_up, isB_up;
  long int is1_up, is2_up;
  long int tmp_off = 0;
  int tmp_sgn, num1, num2, one = 1;
  std::complex<double> tmp_V;
  long int i_max;
  std::complex<double> dmv;

  i_max = Check::idim_max;
  Large::mode = M_CORR;

  for (i = 0; i < Def::NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[0][0]);
    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4,
      &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        is1_up = Def::Tpow[org_isite1];
        is2_up = Def::Tpow[org_isite3];
        num1 = mltply::Spin::GC::Half::X_CisAis((long int)MP::myrank, is1_up, org_sigma1);
        num2 = mltply::Spin::GC::Half::X_CisAis((long int)MP::myrank, is2_up, org_sigma3);
        zaxpy_long(i_max*nstate, tmp_V * (std::complex<double>)(num1*num2),
          &vec[0][0], &Xvec[0][0]);
      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {
        is1_up = Def::Tpow[org_isite1];
        num1 = mltply::Spin::GC::Half::X_CisAis((long int)MP::myrank, is1_up, org_sigma1);
        zaxpy_long(i_max*nstate, tmp_V * (std::complex<double>)num1, &vec[0][0], &Xvec[0][0]);
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {//exchange
        mltply::Spin::C::Half::X_general_int_MPIdouble(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4, 
          tmp_V, nstate, Xvec, vec);
      }
      else {  // other process is not allowed
                // error message will be added
      }
    }
    else if (org_isite1 >= Def::Nsite || org_isite3 >= Def::Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        is1_up = Def::Tpow[org_isite1];
        is2_up = Def::Tpow[org_isite3];
        num2 = mltply::Spin::GC::Half::X_CisAis((long int)MP::myrank, is2_up, org_sigma3);
#pragma omp parallel for default(none)shared(vec,Xvec,nstate,one, \
i_max, tmp_V, is1_up, org_sigma1, num2) private(j, num1,dmv)
        for (j = 0; j < i_max; j++) {
          num1 = mltply::Spin::C::Half::X_CisAis(j, is1_up, org_sigma1);
          dmv = tmp_V * (std::complex<double>)(num1*num2);
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
        }
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {//exchange
        mltply::Spin::C::Half::X_general_int_MPIsingle(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4, 
          tmp_V, nstate, Xvec, vec);
      }
      else {  // other process is not allowed
                // error message will be added
      }
    }
    else {
      isA_up = Def::Tpow[org_isite1];
      isB_up = Def::Tpow[org_isite3];
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,Xvec,nstate,i_max,isA_up,isB_up,org_sigma2,org_sigma4, tmp_V)
        for (j = 0; j < i_max; j++) {
          mltply::Spin::C::Half::CisAisCisAis_element(j, isA_up, isB_up, org_sigma2, org_sigma4,
            tmp_V, nstate, Xvec, vec);
        }
      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma3 == org_sigma2) {
#pragma omp parallel for default(none) private(j, dmv) \
shared(i_max,isA_up,org_sigma1, tmp_V,vec,Xvec,nstate,one)
        for (j = 0; j < i_max; j++) {
          dmv = tmp_V * (std::complex<double>)mltply::Spin::C::Half::X_CisAis(j, 
            isA_up, org_sigma1);
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
        }
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) { // exchange
#pragma omp parallel for default(none) private(j, tmp_sgn, dmv,tmp_off) \
shared(vec,Xvec,nstate,one, i_max,isA_up,isB_up,org_sigma2,org_sigma4,tmp_V)
        for (j = 0; j < i_max; j++) {
          tmp_sgn = mltply::Spin::C::Half::X_exchange_element(j, isA_up, isB_up, org_sigma2, org_sigma4, &tmp_off);
          dmv = tmp_sgn;
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off][0], &one);
        }
      }
      else {  // other process is not allowed
                // error message will be added
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief Child function to calculate two-body green's functions for General Spin model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
void expec::cisajscktalt::Spin::C::General(
  int nstate, 
  std::complex<double> **Xvec, 
  std::complex<double> **vec, 
  std::complex<double> **prod
){
  long int i, j;
  long int org_isite1, org_isite2, org_isite3, org_isite4;
  long int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long int tmp_off = 0;
  long int tmp_off_2 = 0;
  long int list1_off = 0;
  int num1, one = 1;
  std::complex<double> tmp_V;
  long int i_max;
  int tmp_Sz;
  long int tmp_org = 0;
  i_max = Check::idim_max;
  Large::mode = M_CORR;

  for (i = 0; i < Def::NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[0][0]);

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, 
      &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V) != 0) {
      zclear(nstate, prod[i]);
      continue;
    }
    tmp_Sz = 0;

    for (j = 0; j < 2; j++) {
      tmp_org = Def::CisAjtCkuAlvDC[i][4 * j + 1] * Def::Tpow[Def::CisAjtCkuAlvDC[i][4 * j]];
      tmp_Sz += GetLocal2Sz(Def::CisAjtCkuAlvDC[i][4 * j], tmp_org, Def::SiteToBit, Def::Tpow);
      tmp_org = Def::CisAjtCkuAlvDC[i][4 * j + 3] * Def::Tpow[Def::CisAjtCkuAlvDC[i][4 * j + 2]];
      tmp_Sz -= GetLocal2Sz(Def::CisAjtCkuAlvDC[i][4 * j + 2], tmp_org, Def::SiteToBit, Def::Tpow);
    }
    if (tmp_Sz != 0) { // not Sz conserved
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        mltply::Spin::C::General::X_CisAisCjuAju_MPIdouble(
          org_isite1, org_sigma1, org_isite3, org_sigma3, 
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::C::General::X_CisAitCjuAjv_MPIdouble(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4,
          tmp_V, nstate, Xvec, vec);
      }
      else {
      }
    }
    else if (org_isite3 >= Def::Nsite || org_isite1 >= Def::Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        mltply::Spin::C::General::X_CisAisCjuAju_MPIsingle(
          org_isite1, org_sigma1, org_isite3, org_sigma3, tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::C::General::X_CisAitCjuAjv_MPIsingle(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4,
          tmp_V, nstate, Xvec, vec);
      }
      else {
      }
    }
    else {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j, num1) \
shared(vec,List::c1,Xvec,nstate,one, i_max,org_isite1, org_sigma1,org_isite3, org_sigma3, \
tmp_V, Def::SiteToBit, Def::Tpow)
        for (j = 0; j < i_max; j++) {
          num1 = BitCheckGeneral(List::c1[j], org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
          if (num1 != FALSE) {
            num1 = BitCheckGeneral(List::c1[j], org_isite3, org_sigma3, Def::SiteToBit, Def::Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[j][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(list1_off,j,num1,tmp_off,tmp_off_2) \
shared(i_max,org_isite1,org_isite3,org_sigma1,org_sigma2,org_sigma3,org_sigma4, \
MP::myrank,tmp_V,vec,List::c1,Xvec,nstate,one,Def::SiteToBit, Def::Tpow, Check::sdim)
        for (j = 0; j < i_max; j++) {
          num1 = GetOffCompGeneralSpin(List::c1[j], org_isite3, org_sigma4, org_sigma3, &tmp_off, Def::SiteToBit, Def::Tpow);
          if (num1 != FALSE) {
            num1 = GetOffCompGeneralSpin(tmp_off, org_isite1, org_sigma2, org_sigma1, &tmp_off_2,
              Def::SiteToBit, Def::Tpow);
            if (num1 != FALSE) {
              ConvertToList1GeneralSpin(tmp_off_2, Check::sdim, &list1_off);
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[list1_off][0], &one);
            }
          }
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief Child function to calculate two-body green's functions for 1/2 Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
void expec::cisajscktalt::Spin::GC::Half(
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
  std::complex<double> **prod
){
  long int i, j;
  long int org_isite1, org_isite2, org_isite3, org_isite4;
  long int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long int isA_up, isB_up;
  long int tmp_off = 0;
  std::complex<double> tmp_V;
  long int i_max;
  i_max = Check::idim_max;

  for (i = 0; i < Def::NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[0][0]);

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, 
      &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite) { //org_isite3 >= org_isite1 > Nsite

      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        mltply::Spin::GC::Half::X_CisAisCjuAju_MPIdouble(
          org_isite1, org_sigma1, org_isite3, org_sigma3, tmp_V, nstate, Xvec, vec);

      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) { //diagonal (for spin: cuadcdau=cuau)
        mltply::Spin::GC::Half::X_CisAis_MPIdouble(
          org_isite1, org_sigma1, tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::Half::X_CisAisCjuAjv_MPIdouble(
          org_isite1, org_sigma1, org_isite3, org_sigma3, org_sigma4, 
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        mltply::Spin::GC::Half::X_CisAitCjuAju_MPIdouble(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3,
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::Half::X_CisAitCiuAiv_MPIdouble(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4, 
          tmp_V, nstate, Xvec, vec);
      }
    }
    else if (org_isite3 >= Def::Nsite || org_isite1 >= Def::Nsite) { //org_isite3 > Nsite >= org_isite1
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        mltply::Spin::GC::Half::X_CisAisCjuAju_MPIsingle(
          org_isite1, org_sigma1, org_isite3, org_sigma3, tmp_V, nstate, Xvec, vec);

      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::Half::X_CisAisCjuAjv_MPIsingle(
          org_isite1, org_sigma1, org_isite3, org_sigma3, org_sigma4,
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        mltply::Spin::GC::Half::X_CisAitCjuAju_MPIsingle(
          org_isite1, org_sigma2, org_isite3, org_sigma3, 
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::Half::X_CisAitCiuAiv_MPIsingle(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4,
          tmp_V, nstate, Xvec, vec);
      }
    }
    else {
      if (org_isite1 == org_isite2 && org_isite3 == org_isite4) {
        isA_up = Def::Tpow[org_isite2];
        isB_up = Def::Tpow[org_isite4];
        if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j) \
shared(vec,Xvec,nstate, i_max,isA_up,isB_up,org_sigma2,org_sigma4,tmp_V)
          for (j = 0; j < i_max; j++) {
            mltply::Spin::GC::Half::CisAisCisAis_element(j, isA_up, isB_up, org_sigma2, org_sigma4,
              tmp_V, nstate, Xvec, vec);
          }
        }
        else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,Xvec,nstate,i_max,isA_up,isB_up,org_sigma2,org_sigma4,tmp_V)
          for (j = 0; j < i_max; j++) {
            mltply::Spin::GC::Half::CisAisCitAiu_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
              tmp_V, nstate, Xvec, vec, &tmp_off);
          }
        }
        else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,Xvec,nstate,i_max,isA_up,isB_up,org_sigma2,org_sigma4,tmp_V)
          for (j = 0; j < i_max; j++) {
            mltply::Spin::GC::Half::CisAitCiuAiu_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
              tmp_V, nstate, Xvec, vec, &tmp_off);
          }
        }
        else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j,tmp_off) \
shared(vec,Xvec,nstate,i_max,isA_up,isB_up,org_sigma2,org_sigma4,tmp_V)
          for (j = 0; j < i_max; j++) {
            mltply::Spin::GC::Half::CisAitCiuAiv_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
              tmp_V, nstate, Xvec, vec, &tmp_off);
          }
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief Child function to calculate two-body green's functions for General Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
void expec::cisajscktalt::Spin::GC::General(
  int nstate, 
  std::complex<double> **Xvec, 
  std::complex<double> **vec, 
  std::complex<double> **prod
){
  long int i, j;
  long int org_isite1, org_isite2, org_isite3, org_isite4;
  long int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long int tmp_off = 0;
  long int tmp_off_2 = 0;
  int  num1, one = 1;
  std::complex<double> tmp_V;
  long int i_max;
  i_max = Check::idim_max;
  Large::mode = M_CORR;

  for(i=0;i<Def::NCisAjtCkuAlvDC;i++){
    zclear(i_max*nstate, &Xvec[0][0]);

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4,
      &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 >= Def::Nsite && org_isite3 >= Def::Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        mltply::Spin::GC::General::X_CisAisCjuAju_MPIdouble(
          org_isite1, org_sigma1, org_isite3, org_sigma3, tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::General::X_CisAisCjuAjv_MPIdouble(
          org_isite1, org_sigma1, org_isite3, org_sigma3, org_sigma4, 
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        mltply::Spin::GC::General::X_CisAitCjuAju_MPIdouble(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, 
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::General::X_CisAitCjuAjv_MPIdouble(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4,
          tmp_V, nstate, Xvec, vec);
      }
    }
    else if (org_isite3 >= Def::Nsite || org_isite1 >= Def::Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        mltply::Spin::GC::General::X_CisAisCjuAju_MPIsingle(
          org_isite1, org_sigma1, org_isite3, org_sigma3, tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::General::X_CisAisCjuAjv_MPIsingle(
          org_isite1, org_sigma1, org_isite3, org_sigma3, org_sigma4, 
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        mltply::Spin::GC::General::X_CisAitCjuAju_MPIsingle(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, 
          tmp_V, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        mltply::Spin::GC::General::X_CisAitCjuAjv_MPIsingle(
          org_isite1, org_sigma1, org_sigma2, org_isite3, org_sigma3, org_sigma4, 
          tmp_V, nstate, Xvec, vec);
      }
    }
    else {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j, num1) \
shared(vec,Xvec,nstate,one,i_max,org_isite1, org_sigma1,org_isite3, org_sigma3, \
tmp_V, Def::SiteToBit, Def::Tpow)
        for (j = 0; j < i_max; j++) {
          num1 = BitCheckGeneral(j, org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
          if (num1 != FALSE) {
            num1 = BitCheckGeneral(j, org_isite3, org_sigma3, Def::SiteToBit, Def::Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[j][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j, tmp_off, num1) \
shared(vec,Xvec,nstate,one,i_max,org_isite1, org_isite3, org_sigma1,org_sigma3,org_sigma4, \
tmp_V, Def::SiteToBit, Def::Tpow)
        for (j = 0; j < i_max; j++) {
          num1 = GetOffCompGeneralSpin(j, org_isite3, org_sigma4, org_sigma3,
            &tmp_off, Def::SiteToBit, Def::Tpow);
          if (num1 != FALSE) {
            num1 = BitCheckGeneral(tmp_off, org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[tmp_off][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp parallel for default(none) private(j, num1, tmp_off) \
shared(vec,Xvec,nstate,one,i_max,org_isite1, org_isite3, org_sigma1,org_sigma2, \
org_sigma3, tmp_V, Def::SiteToBit, Def::Tpow)
        for (j = 0; j < i_max; j++) {
          num1 = BitCheckGeneral(j, org_isite3, org_sigma3, Def::SiteToBit, Def::Tpow);
          if (num1 != FALSE) {
            num1 = GetOffCompGeneralSpin(j, org_isite1, org_sigma2, org_sigma1,
              &tmp_off, Def::SiteToBit, Def::Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[tmp_off][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(tmp_off,tmp_off_2,j, num1) \
shared(i_max,org_isite1,org_isite3,org_sigma1,org_sigma2,org_sigma3,org_sigma4, \
tmp_V,vec,Xvec,nstate,one, Def::SiteToBit, Def::Tpow)
        for (j = 0; j < i_max; j++) {
          num1 = GetOffCompGeneralSpin(j, org_isite3, org_sigma4, org_sigma3,
            &tmp_off, Def::SiteToBit, Def::Tpow);
          if (num1 != FALSE) {
            num1 = GetOffCompGeneralSpin(tmp_off, org_isite1, org_sigma2, org_sigma1, 
              &tmp_off_2, Def::SiteToBit, Def::Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[tmp_off_2][0], &one);
            }
          }
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief Parent function to calculate two-body green's functions
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 *
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 * @note The origin of function's name cisajscktalt comes from c=creation, i=ith site, s=spin, a=annihiration, j=jth site and so on.
 *
 * @version 0.2
 * @details add function to treat the case of general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int expec::cisajscktalt::main(
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec
) {
  FILE *fp;
  char sdt[D_FileNameMax];
  long int icaca;
  std::complex<double> **prod;
  //For TPQ
  int step = 0, rand_i = 0, istate;

  if (Def::NCisAjtCkuAlvDC < 1) return 0;
  Large::mode = M_CORR;

  //Make File Name for output
  prod = cd_2d_allocate(Def::NCisAjtCkuAlvDC, nstate);
  switch (Def::iCalcType) {
  case DC::TPQCalc:
    step = Def::istep;
    TimeKeeperWithRandAndStep("%s_TimeKeeper.dat", "set %d step %d:expec_cisajscktaltdc finishes: %s", "a", 0, step);
    break;
  case DC::TimeEvolution:
    step = Def::istep;
    TimeKeeperWithStep("%s_TimeKeeper.dat", "step %d:expec_cisajscktaltdc finishes: %s", "a", step);
    break;
  case DC::FullDiag:
  case DC::CG:
    break;
  }

  switch (Def::iCalcModel) {
  case DC::HubbardGC:
    expec::cisajscktalt::Hubbard::GC(nstate, Xvec, vec, prod);
    break;

  case DC::KondoGC:
  case DC::Hubbard:
  case DC::Kondo:
    expec::cisajscktalt::Hubbard::C(nstate, Xvec, vec, prod);
    break;

  case DC::Spin:
    if (Def::iFlgGeneralSpin == FALSE) {
      expec::cisajscktalt::Spin::C::Half(nstate, Xvec, vec, prod);
    }
    else {
      expec::cisajscktalt::Spin::C::General(nstate, Xvec, vec, prod);
    }
    break;

  case DC::SpinGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      expec::cisajscktalt::Spin::GC::Half(nstate, Xvec, vec, prod);
    }
    else {
      expec::cisajscktalt::Spin::GC::General(nstate, Xvec, vec, prod);
    }
    break;

  default:
    break;
  }

  for (istate = 0; istate < nstate; istate++) {
    switch (Def::iCalcType) {
    case DC::TPQCalc:
      step = Def::istep;
      sprintf(sdt, "%s_cisajscktalt_set%dstep%d.dat", Def::CDataFileHead, istate, step);
      break;
    case DC::TimeEvolution:
      step = Def::istep;
      sprintf(sdt, "%s_cisajscktalt_step%d.dat", Def::CDataFileHead, step);
      break;
    case DC::FullDiag:
    case DC::CG:
      sprintf(sdt, "%s_cisajscktalt_eigen%d.dat", Def::CDataFileHead, istate);
      break;
    }
    if (childfopenMPI(sdt, "w", &fp) == 0) {
      for (icaca = 0; icaca < Def::NCisAjtCkuAlvDC; icaca++) {
        fprintf(fp, " %4d %4d %4d %4d %4d %4d %4d %4d %.10lf %.10lf\n",
          Def::CisAjtCkuAlvDC[icaca][0], Def::CisAjtCkuAlvDC[icaca][1],
          Def::CisAjtCkuAlvDC[icaca][2], Def::CisAjtCkuAlvDC[icaca][3],
          Def::CisAjtCkuAlvDC[icaca][4], Def::CisAjtCkuAlvDC[icaca][5],
          Def::CisAjtCkuAlvDC[icaca][6], Def::CisAjtCkuAlvDC[icaca][7],
          real(prod[icaca][istate]), imag(prod[icaca][istate]));
      }
      fclose(fp);
    }
    else return -1;
  }/*for (istate = 0; istate < nstate; istate++)*/

  if (Def::iCalcType == DC::TPQCalc) {
    TimeKeeperWithRandAndStep("%s_TimeKeeper.dat", "set %d step %d:expec_cisajscktaltdc finishes: %s", "a", rand_i, step);
  }
  else if (Def::iCalcType == DC::TimeEvolution) {
    TimeKeeperWithStep("%s_TimeKeeper.dat", "step %d:expec_cisajscktaltdc finishes: %s", "a", step);
  }
  //[s] this part will be added
  /* For FullDiag, it is convinient to calculate the total spin for each vector.
     Such functions will be added
     if(Def::iCalcType==FullDiag){
     if(Def::iCalcModel==Spin){
     expec_cisajscktaltdc_alldiag_spin(vec);
     }else if(Def::iCalcModel==Hubbard || Def::iCalcModel==Kondo){
     expec_cisajscktaltdc_alldiag(vec);
     }else{//
     Phys::s2=0.0;
     }
     }
  */
  //[e]
  free_cd_2d_allocate(prod);
  return 0;
}
