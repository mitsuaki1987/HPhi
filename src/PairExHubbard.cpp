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
/*-------------------------------------------------------------*/
#include "PairExHubbard.hpp"
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "mltplyCommon.hpp"
#include "mltplyHubbard.hpp"
#include "mltplyHubbardCore.hpp"
#include "mltplyMPIHubbard.hpp"
#include "mltplyMPIHubbardCore.hpp"
#include "global.hpp"
#ifdef __MPI
#include "common/setmemory.hpp"
#endif

/**
@breif Calculation of pair excited state for Hubbard Grand canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Pair::Hubbard::GC(
  int nstate, 
  std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {
  long int i, j;
  long int isite1;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;

  std::complex<double> tmp_trans = 0;
  long int ibit;
  long int is;

  for (i = 0; i < Def::NPairExcitationOperator[iEx]; i++) {
    org_isite1 = Def::PairExcitationOperator[iEx][i][0];
    org_isite2 = Def::PairExcitationOperator[iEx][i][2];
    org_sigma1 = Def::PairExcitationOperator[iEx][i][1];
    org_sigma2 = Def::PairExcitationOperator[iEx][i][3];
    tmp_trans = Def::ParaPairExcitationOperator[iEx][i];

    if (org_isite1 >= Def::Nsite &&
      org_isite2 >= Def::Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        if (Def::PairExcitationOperator[iEx][i][4] == 0) {
          if (org_sigma1 == 0) {
            is = Def::Tpow[2 * org_isite1];
          }
          else {
            is = Def::Tpow[2 * org_isite1 + 1];
          }
          ibit = (long int) MP::myrank & is;
          if (ibit != is) {
            //minus sign comes from negative tmp_trans due to readdef
            zaxpy_long(Check::idim_max*nstate, -tmp_trans, &tmp_v1[0][0], &tmp_v0[0][0]);
          }
        }
        else {//Def::PairExcitationOperator[iEx][i][4]==1
          if (org_sigma1 == 0) {
            is = Def::Tpow[2 * org_isite1];
          }
          else {
            is = Def::Tpow[2 * org_isite1 + 1];
          }
          ibit = (long int) MP::myrank & is;
          if (ibit == is) {
            zaxpy_long(Check::idim_max *nstate, tmp_trans, &tmp_v1[0][0], &tmp_v0[0][0]);
          }
        }
      }
      else {
        mltply::Hubbard::GC::X_general_hopp_MPIdouble(org_isite1, org_sigma1, org_isite2, org_sigma2,
          -tmp_trans, nstate, tmp_v0, tmp_v1);
      }
    }
    else if (org_isite2 >= Def::Nsite || org_isite1 >= Def::Nsite) {
      if (org_isite1 < org_isite2) {
        mltply::Hubbard::GC::X_general_hopp_MPIsingle(org_isite1, org_sigma1, org_isite2, org_sigma2,
          -tmp_trans, nstate, tmp_v0, tmp_v1);
      }
      else {
        mltply::Hubbard::GC::X_general_hopp_MPIsingle(org_isite2, org_sigma2, org_isite1, org_sigma1,
          -conj(tmp_trans), nstate, tmp_v0, tmp_v1);
      }
    }
    else {

      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2 && Def::PairExcitationOperator[iEx][i][4] == 0) {
        isite1 = Def::Tpow[2 * org_isite1 + org_sigma1];
#pragma omp parallel for default(none) private(j) \
shared(Check::idim_max,isite1, tmp_trans,tmp_v0,tmp_v1,nstate)
        for (j = 0; j < Check::idim_max; j++) {
          mltply::Hubbard::GC::AisCis(j, nstate, tmp_v0, tmp_v1, isite1, -tmp_trans);
        }
      }
      else {
        mltply::Hubbard::general_hopp_GetInfo(org_isite1, org_isite2, org_sigma1, org_sigma2);
        mltply::Hubbard::GC::general_hopp(nstate, tmp_v0, tmp_v1, tmp_trans);
      }
    }
  }
  return TRUE;
}
/**
@breif Calculation of pair excited state for Hubbard canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Pair::Hubbard::C(
  int nstate, 
  std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {
  long int i, j;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int tmp_off = 0;

  std::complex<double> tmp_trans = 0, dmv;
  int tmp_sgn, num1, one = 1;
  long int ibit;
  long int is, Asum, Adiff;
  long int ibitsite1, ibitsite2;

  //  i_max = Check::idim_max;
  Large::mode = M_CALCSPEC;

  for (i = 0; i < Def::NPairExcitationOperator[iEx]; i++) {
    org_isite1 = Def::PairExcitationOperator[iEx][i][0];
    org_isite2 = Def::PairExcitationOperator[iEx][i][2];
    org_sigma1 = Def::PairExcitationOperator[iEx][i][1];
    org_sigma2 = Def::PairExcitationOperator[iEx][i][3];
    tmp_trans = Def::ParaPairExcitationOperator[iEx][i];
    ibitsite1 = Def::OrgTpow[2 * org_isite1 + org_sigma1];
    ibitsite2 = Def::OrgTpow[2 * org_isite2 + org_sigma2];
    mltply::Hubbard::general_hopp_GetInfo(org_isite1, org_isite2, org_sigma1, org_sigma2);
    Asum = Large::isA_spin;
    Adiff = Large::A_spin;

    if (Def::iFlagListModified == TRUE // Not to adopt HubbrdNConserved
      && org_sigma1 != org_sigma2) {
      if (org_isite1 >= Def::Nsite &&
        org_isite2 >= Def::Nsite)
      {
        mltply::Hubbard::C::X_CisAjt_MPIdouble(org_isite1, org_sigma1, org_isite2, org_sigma2,
          -tmp_trans, nstate, tmp_v0, tmp_v1, Check::idim_max, List::a1, List::b2_1, List::b2_2);
      }
      else if (org_isite2 >= Def::Nsite
        || org_isite1 >= Def::Nsite)
      {
        if (org_isite1 < org_isite2) {
          mltply::Hubbard::C::X_CisAjt_MPIsingle(org_isite1, org_sigma1, org_isite2, org_sigma2,
            -tmp_trans, nstate, tmp_v0, tmp_v1, Check::idim_max, List::a1, List::b2_1, List::b2_2);
        }
        else {
          mltply::Hubbard::C::X_CisAjt_MPIsingle(org_isite2, org_sigma2, org_isite1, org_sigma1,
            -conj(tmp_trans), nstate, tmp_v0, tmp_v1, Check::idim_max, List::a1, List::b2_1, List::b2_2);
        }
      }
      else {
#pragma omp parallel for default(none) private(j,tmp_sgn,tmp_off,dmv) \
shared(tmp_v0,tmp_v1,one,nstate,Check::idim_max,tmp_trans,Asum,Adiff, \
ibitsite1,ibitsite2,List::a1,List::b2_1,List::b2_2,MP::myrank)
        for (j = 0; j < Check::idim_max; j++) {
          tmp_sgn = mltply::Hubbard::C::X_CisAjt(List::a1[j], 
            ibitsite1, ibitsite2, Asum, Adiff, &tmp_off, List::b2_1, List::b2_2);
          dmv = tmp_trans * (std::complex<double>)tmp_sgn;
          zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
        }
      }
    }
    else {
      if (org_isite1 >= Def::Nsite &&
        org_isite2 >= Def::Nsite) {
        if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {//diagonal
          is = Def::Tpow[2 * org_isite1 + org_sigma1];
          ibit = (long int) MP::myrank & is;
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
            if (ibit != is) {
              dmv = -tmp_trans;
#pragma omp parallel for default(none) private(j) \
shared(tmp_v0, tmp_v1,one,dmv,nstate,Check::idim_max, tmp_trans)
              for (j = 0; j < Check::idim_max; j++) {
                zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
              }
            }
          }
          else {
            if (ibit == is) {
              zaxpy_long(Check::idim_max *nstate, tmp_trans, &tmp_v1[0][0], &tmp_v0[0][0]);
            }
          }
        }
        else {
          mltply::Hubbard::C::X_general_hopp_MPIdouble(org_isite1, org_sigma1, org_isite2, org_sigma2,
            -tmp_trans, nstate, tmp_v0, tmp_v1, Check::idim_max, List::a1, List::b2_1, List::b2_2);
        }
      }
      else if (org_isite2 >= Def::Nsite || org_isite1 >= Def::Nsite) {
        if (org_isite1 < org_isite2) {
          mltply::Hubbard::C::X_general_hopp_MPIsingle(org_isite1, org_sigma1, org_isite2, org_sigma2,
            -tmp_trans, nstate, tmp_v0, tmp_v1, Check::idim_max, List::a1, List::b2_1, List::b2_2);
        }
        else {
          mltply::Hubbard::C::X_general_hopp_MPIsingle(org_isite2, org_sigma2, org_isite1, org_sigma1,
            -conj(tmp_trans), nstate, tmp_v0, tmp_v1, Check::idim_max, List::a1, List::b2_1, List::b2_2);
        }
      }
      else {
        mltply::Hubbard::general_hopp_GetInfo(org_isite1, org_isite2, org_sigma1, org_sigma2);
        if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
          is = Def::Tpow[2 * org_isite1 + org_sigma1];
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
#pragma omp parallel for default(none) private(num1,ibit,dmv) \
shared(List::a1,nstate,tmp_v0,tmp_v1,one,Check::idim_max,is,tmp_trans)
            for (j = 0; j < Check::idim_max; j++) {
              ibit = List::a1[j] & is;
              num1 = (1 - ibit / is);
              dmv = -tmp_trans * (std::complex<double>)num1;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
          else {
#pragma omp parallel for default(none) private(num1,ibit,dmv) \
shared(List::a1,nstate,tmp_v0,tmp_v1,one,Check::idim_max,is,tmp_trans)
            for (j = 0; j < Check::idim_max; j++) {
              ibit = List::a1[j] & is;
              num1 = ibit / is;
              dmv = tmp_trans * (std::complex<double>)num1;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
        }
        else {
          mltply::Hubbard::C::general_hopp(nstate, tmp_v0, tmp_v1, tmp_trans, 
            Check::idim_max, List::a1, List::b2_1, List::b2_2);
        }
      }
    }
  }
  return TRUE;
}
