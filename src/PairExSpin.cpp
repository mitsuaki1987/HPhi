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
#include "PairExSpin.hpp"
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "mltplyMPISpinCore.hpp"
#include "mltplySpinCore.hpp"
#include "mltplyCommon.hpp"
#include "global.hpp"
#ifdef __MPI
#include "common/setmemory.hpp"
#endif
#include <complex>

/**
@breif Calculation of pair excited state for Half Spin Grand canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Pair::Spin::GC::Half(
  int nstate,
  std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {
  long int i, j;
  long int isite1;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int tmp_off = 0;

  std::complex<double> tmp_trans = 0, dmv;
  long int i_max;
  int tmp_sgn, one = 1;

  i_max = Check::idim_maxOrg;

  for (i = 0; i < Def::NPairExcitationOperator[iEx]; i++) {
    org_isite1 = Def::PairExcitationOperator[iEx][i][0] + 1;
    org_isite2 = Def::PairExcitationOperator[iEx][i][2] + 1;
    org_sigma1 = Def::PairExcitationOperator[iEx][i][1];
    org_sigma2 = Def::PairExcitationOperator[iEx][i][3];
    tmp_trans = Def::ParaPairExcitationOperator[iEx][i];
    if (org_isite1 == org_isite2) {
      if (org_isite1 > Def::Nsite) {
        if (org_sigma1 == org_sigma2) {  // longitudinal magnetic field
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
            mltply::Spin::GC::Half::X_AisCis_MPIdouble(org_isite1 - 1, org_sigma1, -tmp_trans, nstate, tmp_v0, tmp_v1);
          }
          else {
            mltply::Spin::GC::Half::X_CisAis_MPIdouble(org_isite1 - 1, org_sigma1, tmp_trans, nstate, tmp_v0, tmp_v1);
          }
        }
        else {  // transverse magnetic field
          mltply::Spin::GC::Half::X_CisAit_MPIdouble(org_isite1 - 1, org_sigma1, org_sigma2, tmp_trans, nstate, tmp_v0, tmp_v1);
        }
      }
      else {
        isite1 = Def::Tpow[org_isite1 - 1];
        if (org_sigma1 == org_sigma2) {
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
            // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn,dmv) \
shared(one,nstate,tmp_v0, tmp_v1,i_max, isite1, org_sigma1, tmp_trans)
            for (j = 1; j <= i_max; j++) {
              dmv = (1.0 - mltply::Spin::GC::Half::X_CisAis(j, isite1, org_sigma1))
                * (-tmp_trans);
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
          else {
            // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn,dmv)             \
shared(tmp_v0, tmp_v1,one,nstate,i_max, isite1, org_sigma1, tmp_trans)
            for (j = 1; j <= i_max; j++) {
              dmv = (std::complex<double>)mltply::Spin::GC::Half::X_CisAis(j, isite1, org_sigma1)
                * tmp_trans;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
        }
        else {
          // transverse magnetic field
          // fprintf(MP::STDOUT, "Debug: isite1=%d, org_sigma2=%d\n", isite1, org_sigma2);
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off,dmv)    \
shared(tmp_v0, tmp_v1,one,nstate,i_max, isite1, org_sigma2, tmp_trans)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = mltply::Spin::GC::Half::X_CisAit(j, isite1, org_sigma2, &tmp_off);
            if (tmp_sgn != 0) {
              dmv = (std::complex<double>)tmp_sgn * tmp_trans;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off + 1], &one);
            }
          }
        }
      }
    }
    else {
      fprintf(MP::STDOUT, "ERROR: hopping is not allowed in localized spin system\n");
      return FALSE;
    }
  }
  return TRUE;
}
/**
@breif Calculation of pair excited state for general Spin Grand canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Pair::Spin::GC::General(
  int nstate,
  std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {
  long int i, j;
  int num1;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int tmp_off = 0;
  int one = 1;
  std::complex<double> tmp_trans = 0, dmv;
  long int i_max;
  i_max = Check::idim_maxOrg;

  for (i = 0; i < Def::NPairExcitationOperator[iEx]; i++) {
    org_isite1 = Def::PairExcitationOperator[iEx][i][0] + 1;
    org_isite2 = Def::PairExcitationOperator[iEx][i][2] + 1;
    org_sigma1 = Def::PairExcitationOperator[iEx][i][1];
    org_sigma2 = Def::PairExcitationOperator[iEx][i][3];
    tmp_trans = Def::ParaPairExcitationOperator[iEx][i];
    if (org_isite1 == org_isite2) {
      if (org_isite1 > Def::Nsite) {
        if (org_sigma1 == org_sigma2) {
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
            // longitudinal magnetic field
            mltply::Spin::GC::General::X_AisCis_MPIdouble(org_isite1 - 1, org_sigma1, -tmp_trans, nstate, tmp_v0, tmp_v1);
          }
          else {
            mltply::Spin::GC::General::X_CisAis_MPIdouble(org_isite1 - 1, org_sigma1, tmp_trans, nstate, tmp_v0, tmp_v1);
          }
        }
        else {
          // transverse magnetic field
          mltply::Spin::GC::General::X_CisAit_MPIdouble(org_isite1 - 1, org_sigma1, org_sigma2, tmp_trans, nstate, tmp_v0, tmp_v1);
        }
      }
      else {//org_isite1 <= Def::Nsite
        if (org_sigma1 == org_sigma2) {
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
            // longitudinal magnetic field
#pragma omp parallel for default(none) private(j,num1,dmv) \
shared(tmp_v0,tmp_v1,one,nstate,i_max,org_isite1,org_sigma1,tmp_trans,Def::SiteToBit, Def::Tpow)
            for (j = 1; j <= i_max; j++) {
              num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
              dmv = -tmp_trans * (1.0 - num1);
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
          else {
            // longitudinal magnetic field
#pragma omp parallel for default(none) private(j,num1,dmv) \
shared(tmp_v0,tmp_v1,one,nstate,i_max,org_isite1,org_sigma1,tmp_trans,Def::SiteToBit,Def::Tpow)
            for (j = 1; j <= i_max; j++) {
              num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
              dmv = tmp_trans * (std::complex<double>)num1;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
        }
        else {
          // transverse magnetic field
#pragma omp parallel for default(none) private(j,num1,dmv,tmp_off) \
shared(tmp_v0,tmp_v1,one,nstate,i_max,org_isite1,org_sigma1,org_sigma2,tmp_trans, \
Def::SiteToBit, Def::Tpow)
          for (j = 1; j <= i_max; j++) {
            num1 = GetOffCompGeneralSpin(j - 1, org_isite1, org_sigma2, org_sigma1, &tmp_off, Def::SiteToBit, Def::Tpow);
            if (num1 != 0) {
              dmv = tmp_trans * (std::complex<double>)num1;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off + 1], &one);
            }
          }
        }
      }
    }
    else {
      fprintf(MP::STDOUT, "ERROR: hopping is not allowed in localized spin system\n");
      return FALSE;
    }
  }
  return TRUE;
}
/// \brief Calculation of pair excited state for Spin Grand canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetExcitedState::Pair::Spin::GC::main(
  /**< [in,out] define list to get and put information of calculation*/
  int nstate, std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {

  int iret = 0;
  if (Def::iFlgGeneralSpin == FALSE) {
    iret = GetExcitedState::Pair::Spin::GC::Half(nstate, tmp_v0, tmp_v1, iEx);
  }
  else {
    iret = GetExcitedState::Pair::Spin::GC::General(nstate, tmp_v0, tmp_v1, iEx);
  }
  return iret;
}
/**
@breif Calculation of pair excited state for Half Spin canonical system
returns TRUE: Normally finished
returns FALSE: Abnormally finished
author Kazuyoshi Yoshimi
version 1.2
*/
int GetExcitedState::Pair::Spin::C::Half(
  int nstate, 
  std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
)
{
  long int i, j;
  long int isite1;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int tmp_off = 0;
  std::complex<double> tmp_trans = 0, dmv;
  long int i_max;
  int num1, one = 1;
  long int ibit1;
  long int is1_up;

  i_max = Check::idim_maxOrg;

  for (i = 0; i < Def::NPairExcitationOperator[iEx]; i++) {
    org_isite1 = Def::PairExcitationOperator[iEx][i][0] + 1;
    org_isite2 = Def::PairExcitationOperator[iEx][i][2] + 1;
    org_sigma1 = Def::PairExcitationOperator[iEx][i][1];
    org_sigma2 = Def::PairExcitationOperator[iEx][i][3];
    tmp_trans = Def::ParaPairExcitationOperator[iEx][i];
    if (org_sigma1 == org_sigma2) {
      if (org_isite1 == org_isite2) {
        if (org_isite1 > Def::Nsite) {
          is1_up = Def::Tpow[org_isite1 - 1];
          ibit1 = mltply::Spin::GC::Half::X_CisAis((long int) MP::myrank + 1, is1_up, org_sigma1);
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
            if (ibit1 == 0) {
              dmv = -tmp_trans;
#pragma omp parallel for default(none) private(j) \
shared(tmp_v0, tmp_v1,one,nstate,dmv,i_max, tmp_trans)
              for (j = 1; j <= i_max; j++) {
                zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
              }
            }
          }
          else {
            if (ibit1 != 0) {
#pragma omp parallel for default(none) private(j) \
shared(tmp_v0, tmp_v1,one,nstate,i_max, tmp_trans)
              for (j = 1; j <= i_max; j++)
                zaxpy_(&nstate, &tmp_trans, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
        }// org_isite1 > Def::Nsite
        else {
          isite1 = Def::Tpow[org_isite1 - 1];
          if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2 &&
            Def::PairExcitationOperator[iEx][i][4] == 0) {
#pragma omp parallel for default(none) private(j,dmv) \
shared(tmp_v0,tmp_v1,one,nstate,i_max,isite1,org_sigma1,tmp_trans)
            for (j = 1; j <= i_max; j++) {
              dmv = (1.0 - mltply::Spin::C::Half::X_CisAis(j, isite1, org_sigma1)) * (-tmp_trans);
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
          else {
#pragma omp parallel for default(none) private(j,dmv) \
shared(tmp_v0,tmp_v1,one,nstate,i_max,isite1,org_sigma1,tmp_trans)
            for (j = 1; j <= i_max; j++) {
              dmv = (std::complex<double>)mltply::Spin::C::Half::X_CisAis(j, isite1, org_sigma1) 
                * tmp_trans;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
        }
      }
      else {
        fprintf(MP::STDOUT, "Error: isite1 must be equal to isite2 for Spin system. \n");
        return FALSE;
      }
    }
    else { //org_sigma1 != org_sigma2             // for the canonical case
      if (org_isite1 > Def::Nsite) {//For MPI
        mltply::Spin::C::Half::X_CisAit_MPIdouble(org_isite1 - 1, org_sigma2, tmp_trans, 
          nstate, tmp_v0, tmp_v1, i_max);
      }
      else {
        isite1 = Def::Tpow[org_isite1 - 1];
#pragma omp parallel for default(none) private(j,tmp_off,num1,dmv) \
shared(tmp_v0,tmp_v1,one,nstate,i_max,isite1,org_sigma2,tmp_trans)
        for (j = 1; j <= i_max; j++) {
          num1 = mltply::Spin::C::Half::X_CisAit(j, isite1, org_sigma2, &tmp_off);
          if (num1 != 0) {
            dmv = tmp_trans*(double)num1;
            zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
          }
        }
      }
    }
  }
  return TRUE;
}
/**
@breif Calculation of pair excited state for general Spin canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Pair::Spin::C::General(
  int nstate,
  std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
)
{
  long int i, j;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int tmp_off = 0;
  long int off = 0;
  std::complex<double> tmp_trans = 0, dmv;
  long int i_max;
  int tmp_sgn, num1, one = 1;
  i_max = Check::idim_maxOrg;

  for (i = 0; i < Def::NPairExcitationOperator[iEx]; i++) {
    org_isite1 = Def::PairExcitationOperator[iEx][i][0] + 1;
    org_isite2 = Def::PairExcitationOperator[iEx][i][2] + 1;
    org_sigma1 = Def::PairExcitationOperator[iEx][i][1];
    org_sigma2 = Def::PairExcitationOperator[iEx][i][3];
    tmp_trans = Def::ParaPairExcitationOperator[iEx][i];
    if (org_isite1 == org_isite2) {
      if (org_isite1 > Def::Nsite) {
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          num1 = BitCheckGeneral((long int) MP::myrank,
            org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
            if (num1 == 0) {
#pragma omp parallel for default(none) private(j,dmv)  \
shared(tmp_v0, tmp_v1,one,nstate,i_max, tmp_trans)
              for (j = 1; j <= i_max; j++) {
                dmv = -tmp_trans;
                zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
              }
            }
          }
          else {
            if (num1 != 0) {
#pragma omp parallel for default(none) private(j) \
shared(tmp_v0, tmp_v1,one,nstate,i_max, tmp_trans)
              for (j = 1; j <= i_max; j++) {
                zaxpy_(&nstate, &tmp_trans, tmp_v1[j], &one, tmp_v0[j], &one);
              }
            }
          }
        }//org_sigma1=org_sigma2
        else {//org_sigma1 != org_sigma2
          mltply::Spin::C::General::X_CisAit_MPIdouble(org_isite1 - 1, org_sigma1, org_sigma2, 
            tmp_trans, nstate, tmp_v0, tmp_v1, i_max);
        }
      }
      else {//org_isite1 <= Def::Nsite
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          if (Def::PairExcitationOperator[iEx][i][4] == 0) {
#pragma omp parallel for default(none) private(j, num1,dmv) \
shared(tmp_v0, tmp_v1, List::c1,one,nstate,i_max, org_isite1, org_sigma1, tmp_trans, \
Def::SiteToBit, Def::Tpow)
            for (j = 1; j <= i_max; j++) {
              num1 = BitCheckGeneral(List::c1[j], org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
              dmv = -tmp_trans * (1.0 - num1);
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
          else {
#pragma omp parallel for default(none) private(j, num1,dmv) \
shared(tmp_v0, tmp_v1, List::c1,one,nstate,i_max, org_isite1, \
org_sigma1, tmp_trans, Def::SiteToBit, Def::Tpow)
            for (j = 1; j <= i_max; j++) {
              num1 = BitCheckGeneral(List::c1[j], org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
              dmv = tmp_trans * (std::complex<double>)num1;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
        }//org_sigma1=org_sigma2
        else {//org_sigma1 != org_sigma2
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off, off) \
shared(tmp_v0, tmp_v1, List::c1_org,one,nstate, i_max, org_isite1, Large::ihfbit, \
org_sigma1, org_sigma2, tmp_trans, MP::myrank, Def::SiteToBit, Def::Tpow)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = GetOffCompGeneralSpin(List::c1_org[j], org_isite1, org_sigma2, org_sigma1, &off,
              Def::SiteToBit, Def::Tpow);
            if (tmp_sgn != FALSE) {
              ConvertToList1GeneralSpin(off, Large::ihfbit, &tmp_off);
              zaxpy_(&nstate, &tmp_trans, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
            }
          }
        }
      }
    }
    else {
      fprintf(MP::STDOUT, "ERROR: hopping is not allowed in localized spin system\n");
      return FALSE;
    }//org_isite1 != org_isite2
  }

  return TRUE;
}
/// Calculation of pair excited state for Spin canonical system
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetExcitedState::Pair::Spin::C::main(
  int nstate,
  std::complex<double> **tmp_v0, /**< [out] Result v0 = H v1*/
  std::complex<double> **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {
  int iret = 0;
  if (Def::iFlgGeneralSpin == FALSE) {
    iret = GetExcitedState::Pair::Spin::C::Half(nstate, tmp_v0, tmp_v1, iEx);
  }
  else {
    iret = GetExcitedState::Pair::Spin::C::General(nstate, tmp_v0, tmp_v1, iEx);
  }
  return iret;
}
