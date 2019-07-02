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

/**
 * @file
 * @version 2.1
 * @details add functions to calculate diagonal components for Time evolution.
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @version 0.2
 * @details modify functions to calculate diagonal components for general spin.
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  Calculate diagonal components, i.e. @f$ H_d |\phi_0> = E_d |\phi_0> @f$.
 * 
 * 
 */
#include "bitcalc.hpp"
#include "FileIO.hpp"
#include "diagonalcalc.hpp"
#include "mltplySpinCore.hpp"
#include "wrapperMPI.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include "log.hpp"
#include <iostream>
/**
 * @brief Update the vector by the general two-body diagonal interaction, \f$ H_{i\sigma_1 j\sigma_2} n_ {i\sigma_1}n_{j\sigma_2}\f$.\n
 * (Using in Time Evolution mode).
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param isigma1 [in] a spin index at \f$i \f$ site.
 * @param isigma2 [in] a spin index at \f$j \f$ site.
 * @param dtmp_V [in] A value of general two-body diagonal interaction \f$ H_{i\sigma_1 j\sigma_2} \f$
 * @param X [in] Define list to get the operator information.
 * @param tmp_v0 [in,out] Result vector
 * @param tmp_v1 [in] Input produced vector
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 2.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalTEInterAll(
  long int isite1,
  long int isite2,
  long int isigma1,
  long int isigma2,
  double dtmp_V,
  
  std::complex<double> *tmp_v0,
  std::complex<double> *tmp_v1
) {
  long int is1_spin;
  long int is2_spin;
  long int is1_up;
  long int is2_up;

  long int ibit1_spin;
  long int ibit2_spin;

  long int num1;
  long int num2;

  long int j;
  long int i_max = Check::idim_max;
  /*
  Forse isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
    j = isigma2;
    isigma2 = isigma1;
    isigma1 = j;
  }
  /*
  When isite1 & site2 are in the inter process regino
  */
  if (isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:
      is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];
      num1 = 0;
      ibit1_spin = (long int)myrank&is1_spin;
      num1 += ibit1_spin / is1_spin;
      num2 = 0;
      ibit2_spin = (long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;
      break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:
      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        is2_up = Def::Tpow[isite2 - 1];
        num1 = X_SpinGC_CisAis((long int) myrank + 1, is1_up, isigma1);
        num2 = X_SpinGC_CisAis((long int) myrank + 1, is2_up, isigma2);
      }/*if (Def::iFlgGeneralSpin == FALSE)*/
      else {//start:generalspin
        num1 = BitCheckGeneral((long int) myrank, isite1, isigma1,
          Def::SiteToBit, Def::Tpow);
        num2 = BitCheckGeneral((long int) myrank, isite2, isigma2,
          Def::SiteToBit, Def::Tpow);
      }
      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
    }/*if (isite1 > Def::Nsite)*/

    if (num1 * num2 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V) private(j)
      for (j = 1; j <= i_max; j++) {
        tmp_v0[j] += dtmp_V * tmp_v1[j];
      }
    }
    return 0;
  }/*if (isite1 > Def::Nsite)*/
  else if (isite2 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:

      is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];

      num2 = 0;
      ibit2_spin = (long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;
      if (num2 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)\
                     firstprivate(i_max, dtmp_V, is1_spin) private(num1, ibit1_spin, j)
        for (j = 1; j <= i_max; j++) {
          num1 = 0;
          ibit1_spin = (j - 1) & is1_spin;
          num1 += ibit1_spin / is1_spin;
          tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
        }
      }
      break;/*case HubbardGC:*/

    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];

      num2 = 0;
      ibit2_spin = (long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;
      if (num2 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, list_1)\
                     firstprivate(i_max, dtmp_V, is1_spin) private(num1, ibit1_spin, j)
        for (j = 1; j <= i_max; j++) {
          num1 = 0;
          ibit1_spin = list_1[j] & is1_spin;
          num1 += ibit1_spin / is1_spin;
          tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
        }
      }
      break;/*case KondoGC, Hubbard, Kondo:*/

    case SpinGC:

      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        is2_up = Def::Tpow[isite2 - 1];
        num2 = X_SpinGC_CisAis((long int)myrank + 1, is2_up, isigma2);

        if (num2 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)\
                     firstprivate(i_max, dtmp_V, is1_up, isigma1, X) private(num1, j)
          for (j = 1; j <= i_max; j++) {
            num1 = X_SpinGC_CisAis(j, is1_up, isigma1);
            tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
          }
        }
      }/* if (Def::iFlgGeneralSpin == FALSE)*/
      else {//start:generalspin
        num2 = BitCheckGeneral((long int)myrank, isite2, isigma2,
          Def::SiteToBit, Def::Tpow);
        if (num2 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(j - 1, isite1, isigma1, Def::SiteToBit, Def::Tpow);
            tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
          }
        }
      }/* if (Def::iFlgGeneralSpin == TRUE)*/

      break;/*case SpinGC:*/

    case Spin:

      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        is2_up = Def::Tpow[isite2 - 1];
        num2 = X_SpinGC_CisAis((long int)myrank + 1, is2_up, isigma2);

        if (num2 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, num2) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = X_Spin_CisAis(j, is1_up, isigma1);
            tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
          }
        }
      }/* if (Def::iFlgGeneralSpin == FALSE)*/
      else /* if (Def::iFlgGeneralSpin == TRUE)*/ {
        num2 = BitCheckGeneral((long int)myrank, isite2, isigma2, \
          Def::SiteToBit, Def::Tpow);
        if (num2 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, list_1)\
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(list_1[j], isite1, isigma1, Def::SiteToBit, Def::Tpow);
            tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
          }
        }
      } /* if (Def::iFlgGeneralSpin == TRUE)*/

      break;/*case Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    }/*switch (Def::iCalcModel)*/
    return 0;
  }/*else if (isite2 > Def::Nsite)*/

  switch (Def::iCalcModel) {
  case HubbardGC: //list_1[j] -> j-1
    is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
    is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];
#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) \
private(num1, ibit1_spin, num2, ibit2_spin)
    for (j = 1; j <= i_max; j++) {
      num1 = 0;
      num2 = 0;
      ibit1_spin = (j - 1)&is1_spin;
      num1 += ibit1_spin / is1_spin;
      ibit2_spin = (j - 1)&is2_spin;
      num2 += ibit2_spin / is2_spin;
      tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
    }
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
    is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];

#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1, list_1) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) \
private(num1, ibit1_spin, num2, ibit2_spin)
    for (j = 1; j <= i_max; j++) {
      num1 = 0;
      num2 = 0;
      ibit1_spin = list_1[j] & is1_spin;
      num1 += ibit1_spin / is1_spin;

      ibit2_spin = list_1[j] & is2_spin;
      num2 += ibit2_spin / is2_spin;
      tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
    }
    break;

  case Spin:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1)  firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) \
private(j, num1, num2)
      for (j = 1; j <= i_max; j++) {
        num1 = X_Spin_CisAis(j, is1_up, isigma1);
        num2 = X_Spin_CisAis(j, is2_up, isigma2);
        tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
      }
    }
    else {
#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1, list_1) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) \
private(j, num1)
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(list_1[j], isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          num1 = BitCheckGeneral(list_1[j], isite2, isigma2, Def::SiteToBit, Def::Tpow);
          tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
        }
      }

    }
    break;

  case SpinGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) \
private(j, num1, num2)
      for (j = 1; j <= i_max; j++) {
        num1 = X_SpinGC_CisAis(j, is1_up, isigma1);
        num2 = X_SpinGC_CisAis(j, is2_up, isigma2);
        tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
      }
    }
    else {//start:generalspin
#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) \
private(j, num1)
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(j - 1, isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          num1 = BitCheckGeneral(j - 1, isite2, isigma2, Def::SiteToBit, Def::Tpow);
          tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
        }
      }
    }
    break;

  default:
    fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
    return -1;
  }
  return 0;
}
/**
 * @brief Update the vector by the chemical potential \f$ \mu_{i \sigma_1} n_ {i \sigma_1} \f$ \n
 * generated by the commutation relation in terms of the general two-body interaction, \n
 * \f$ c_ {i \sigma_1} a_{j\sigma_2}c_ {j \sigma_2}a_ {i \sigma_1} = c_ {i \sigma_1}a_ {i \sigma_1}-c_ {i \sigma_1} a_ {i \sigma_1} c_ {j \sigma_2}a_{j\sigma_2}\f$ .
 * (Using in Time Evolution mode).
 * @param isite1 [in] a site number
 * @param spin [in] a spin number
 * @param dtmp_V [in] A value of coulombintra interaction \f$ \mu_{i \sigma_1} \f$
 * @param X [in] Define list to get dimension number
 * @param tmp_v0 [in,out] Result vector
 * @param tmp_v1 [in] Input produced vector
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 2.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalTEChemi(
  long int isite1,
  long int spin,
  double dtmp_V,
  
  std::complex<double> *tmp_v0,
  std::complex<double> *tmp_v1
) {
  long int is1_up;
  long int num1;
  long int isigma1 = spin;
  long int is1, ibit1;

  long int j;
  long int i_max = Check::idim_max;

  /*
    When isite1 is in the inter process region
  */
  if (isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      if (spin == 0) {
        is1 = Def::Tpow[2 * isite1 - 2];
      }
      else {
        is1 = Def::Tpow[2 * isite1 - 1];
      }
      ibit1 = (long int)myrank & is1;
      num1 = ibit1 / is1;
      break;/*case HubbardGC, case KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:

      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        num1 = (((long int)myrank& is1_up) / is1_up) ^ (1 - spin);
      } /*if (Def::iFlgGeneralSpin == FALSE)*/
      else /*if (Def::iFlgGeneralSpin == TRUE)*/ {
        num1 = BitCheckGeneral((long int)myrank,
          isite1, isigma1, Def::SiteToBit, Def::Tpow);
      }/*if (Def::iFlgGeneralSpin == TRUE)*/
      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    } /*switch (Def::iCalcModel)*/
    if (num1 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)  \
firstprivate(i_max, dtmp_V) private(j)
      for (j = 1; j <= i_max; j++) {
        tmp_v0[j] += dtmp_V * tmp_v1[j];
      }
    }/*if (num1 != 0)*/
    return 0;

  }/*if (isite1 >= Def::Nsite*/

  switch (Def::iCalcModel) {
  case HubbardGC:
    if (spin == 0) {
      is1 = Def::Tpow[2 * isite1 - 2];
    }
    else {
      is1 = Def::Tpow[2 * isite1 - 1];
    }

#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
    for (j = 1; j <= i_max; j++) {
      ibit1 = (j - 1)&is1;
      num1 = ibit1 / is1;
      tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
    }
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    if (spin == 0) {
      is1 = Def::Tpow[2 * isite1 - 2];
    }
    else {
      is1 = Def::Tpow[2 * isite1 - 1];
    }

#pragma omp parallel for default(none) \
shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
    for (j = 1; j <= i_max; j++) {
      ibit1 = list_1[j] & is1;
      num1 = ibit1 / is1;
      tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
    }
    break;

  case SpinGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
#pragma omp parallel for default(none) \
shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1)
      for (j = 1; j <= i_max; j++) {
        num1 = (((j - 1)& is1_up) / is1_up) ^ (1 - spin);
        tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
      }
    }
    else {
#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(j - 1, isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          tmp_v0[j] += dtmp_V * tmp_v1[j];
        }
      }
    }
    break;

  case Spin:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
#pragma omp parallel for default(none) \
shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1)
      for (j = 1; j <= i_max; j++) {
        num1 = ((list_1[j] & is1_up) / is1_up) ^ (1 - spin);
        tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
      }
    }
    else {
#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1, list_1) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(list_1[j], isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          tmp_v0[j] += dtmp_V * tmp_v1[j];

        }
      }
    }

    break;
  default:
    fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
    return -1;
  }
  return 0;
}
/**
 * @brief Update the vector by the general one-body diagonal interaction, \f$ \mu_{i\sigma_1} n_ {i\sigma_1}\f$.\n
 * (Using in Time Evolution mode).
 * @param isite1 [in] a site number \f$i \f$
 * @param dtmp_V [in] A value of general one-body diagonal interaction \f$ \mu_{i\sigma_1} \f$
 * @param spin [in] a spin index at \f$i \f$ site.
 * @param X [in] Define list to get the operator information.
 * @param tmp_v0 [in,out] Result vector
 * @param tmp_v1 [in] Input produced vector
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 2.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalTETransfer
(
  long int isite1,
  double dtmp_V,
  long int spin,
  
  std::complex<double> *tmp_v0,
  std::complex<double> *tmp_v1
) {
  long int is1_up;
  long int ibit1_up;
  long int num1;
  long int isigma1 = spin;
  long int is1, ibit1;

  long int j;
  long int i_max = Check::idim_max;

  /*
    When isite1 is in the inter process region
  */
  if (isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:
      if (spin == 0) {
        is1 = Def::Tpow[2 * isite1 - 2];
      }
      else {
        is1 = Def::Tpow[2 * isite1 - 1];
      }
      ibit1 = (long int)myrank & is1;
      num1 = ibit1 / is1;
      break;/*case HubbardGC, case KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:
      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        num1 = (((long int)myrank& is1_up) / is1_up) ^ (1 - spin);
      } /*if (Def::iFlgGeneralSpin == FALSE)*/
      else /*if (Def::iFlgGeneralSpin == TRUE)*/ {
        num1 = BitCheckGeneral((long int)myrank,
          isite1, isigma1, Def::SiteToBit, Def::Tpow);
      }/*if (Def::iFlgGeneralSpin == TRUE)*/
      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    } /*switch (Def::iCalcModel)*/

    if (num1 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)\
                     firstprivate(i_max, dtmp_V) private(j)
      for (j = 1; j <= i_max; j++) {
        tmp_v0[j] += dtmp_V * tmp_v1[j];
      }
    }
  }/*if (isite1 >= Def::Nsite*/
  else {//(isite1 < Def::Nsite)
    switch (Def::iCalcModel) {
    case HubbardGC:
      if (spin == 0) {
        is1 = Def::Tpow[2 * isite1 - 2];
      }
      else {
        is1 = Def::Tpow[2 * isite1 - 1];
      }
#pragma omp parallel for default(none) shared(list_1, tmp_v0, tmp_v1) \
        firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
      for (j = 1; j <= i_max; j++) {
        ibit1 = (j - 1) & is1;
        num1 = ibit1 / is1;
        tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
      }
      break;

    case KondoGC:
    case Hubbard:
    case Kondo:
      if (spin == 0) {
        is1 = Def::Tpow[2 * isite1 - 2];
      }
      else {
        is1 = Def::Tpow[2 * isite1 - 1];
      }
#pragma omp parallel for default(none) shared(list_1, tmp_v0, tmp_v1) \
        firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
      for (j = 1; j <= i_max; j++) {
        ibit1 = list_1[j] & is1;
        num1 = ibit1 / is1;
        tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
      }
      break;

    case SpinGC:
      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
#pragma omp parallel for default(none) \
shared(list_1, tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
        for (j = 1; j <= i_max; j++) {
          ibit1_up = (((j - 1) & is1_up) / is1_up) ^ (1 - spin);
          tmp_v0[j] += dtmp_V * ibit1_up*tmp_v1[j];
        }
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
        for (j = 1; j <= i_max; j++) {
          num1 = BitCheckGeneral(j - 1, isite1, isigma1, Def::SiteToBit, Def::Tpow);
          if (num1 != 0) {
            tmp_v0[j] += dtmp_V * tmp_v1[j];
          }
        }
      }
      break;

    case Spin:
      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
#pragma omp parallel for default(none) shared(list_1, tmp_v0, tmp_v1)\
           firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
        for (j = 1; j <= i_max; j++) {
          ibit1_up = ((list_1[j] & is1_up) / is1_up) ^ (1 - spin);
          tmp_v0[j] += dtmp_V * ibit1_up * tmp_v1[j];
        }
      }
      else {
#pragma omp parallel for default(none) shared(list_1, tmp_v0, tmp_v1)\
          firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
        for (j = 1; j <= i_max; j++) {
          num1 = BitCheckGeneral(list_1[j], isite1, isigma1, Def::SiteToBit, Def::Tpow);
          tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
        }
      }
      break;

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
    }
  }
  return 0;
}
///  @fn diagonalcalcForTE() Update the vector for diagonal operators ( using in Time Evolution mode).
/// \param X [in] Struct to get the information of the diagonal operators.
/// \param tmp_v0 [in,out] Result vector
/// \param tmp_v1 [in] Input produced vector
/// \retval -1 fail to update the vector.
/// \retval  0 succeed to update the vector.
/// \version 2.1
int diagonalcalcForTE
(
  const int _istep,
  
  std::complex<double> *tmp_v0,
  std::complex<double> *tmp_v1
) {

  long int i;
  long int isite1, isite2;
  long int A_spin, B_spin;
  double tmp_V;

  if (Def::NTETransferDiagonal[_istep] > 0) {
    for (i = 0; i < Def::NTETransferDiagonal[_istep]; i++) {
      isite1 = Def::TETransferDiagonal[_istep][i][0] + 1;
      A_spin = Def::TETransferDiagonal[_istep][i][1];
      tmp_V = -Def::ParaTETransferDiagonal[_istep][i];
      SetDiagonalTETransfer(isite1, tmp_V, A_spin, tmp_v0, tmp_v1);
    }
  }
  else if (Def::NTEInterAllDiagonal[_istep] > 0) {
    for (i = 0; i < Def::NTEInterAllDiagonal[_istep]; i++) {
      //Assume n_{1\sigma_1} n_{2\sigma_2}
      isite1 = Def::TEInterAllDiagonal[_istep][i][0] + 1;
      A_spin = Def::TEInterAllDiagonal[_istep][i][1];
      isite2 = Def::TEInterAllDiagonal[_istep][i][2] + 1;
      B_spin = Def::TEInterAllDiagonal[_istep][i][3];
      tmp_V = Def::ParaTEInterAllDiagonal[_istep][i];
      if (SetDiagonalTEInterAll(isite1, isite2, A_spin, B_spin, tmp_V, tmp_v0, tmp_v1) != 0) {
        return -1;
      }
    }

    if (Def::NTEChemi[_istep] > 0) {
      for (i = 0; i < Def::NTEChemi[_istep]; i++) {
        isite1 = Def::TEChemi[_istep][i] + 1;
        A_spin = Def::SpinTEChemi[_istep][i];
        tmp_V = -Def::ParaTEChemi[_istep][i];
        if (SetDiagonalTEChemi(isite1, A_spin, tmp_V, tmp_v0, tmp_v1) != 0) {
          return -1;
        }
      }
    }
  }
  return 0;
}
/**
 * @brief Calculate the components for Coulombintra interaction, \f$ U_i n_ {i \uparrow}n_{i \downarrow} \f$
 * @param isite1 [in] a site number
 * @param dtmp_V [in] A value of coulombintra interaction \f$ U_i \f$
 * @param X [in] Define list to get dimension number
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalCoulombIntra
(
  long int isite1,
  double dtmp_V
) {
  long int is;
  long int ibit;
  long int is1_up, is1_down;

  long int j;
  long int i_max = Check::idim_max;

  /*
   When isite1 is in the inter process region
  */
  if (isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is = is1_up + is1_down;
      ibit = (long int)myrank & is;
      if (ibit == is) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
                                       firstprivate(i_max, dtmp_V) private(j) 
        for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V;
      }

      break; /*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case Spin:
    case SpinGC:
      /*
       They do not have the Coulomb term
      */
      break;

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
      //break;

    }/*switch (Def::iCalcModel)*/

    return 0;

  }/*if (isite1 >= Def::Nsite*/
  else {
    switch (Def::iCalcModel) {
    case HubbardGC:
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is = is1_up + is1_down;
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, is, dtmp_V) private(ibit) 
      for (j = 1; j <= i_max; j++) {
        ibit = (j - 1)&is;
        if (ibit == is) {
          list_Diagonal[j] += dtmp_V;
        }
      }

      break;
    case KondoGC:
    case Hubbard:
    case Kondo:
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is = is1_up + is1_down;
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, is, dtmp_V) private(ibit) 
      for (j = 1; j <= i_max; j++) {
        ibit = list_1[j] & is;
        if (ibit == is) {
          list_Diagonal[j] += dtmp_V;
        }
      }
      break;

    case Spin:
    case SpinGC:
      break;

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
      //break;
    }
  }
  return 0;
}
/**
 * @brief Calculate the components for the chemical potential \f$ \mu_{i \sigma_1} n_ {i \sigma_1} \f$
 * @param isite1 [in] a site number
 * @param dtmp_V [in] A value of coulombintra interaction \f$ \mu_{i \sigma_1} \f$
 * @param spin [in] Spin index for the chemical potential
 * @param X [in] Define list to get dimension number
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalChemi
(
  long int isite1,
  double dtmp_V,
  long int spin
) {
  long int is1_up;
  long int ibit1_up;
  long int num1;
  long int isigma1 = spin;
  long int is1, ibit1;

  long int j;
  long int i_max = Check::idim_max;

  /*
    When isite1 is in the inter process region
  */
  if (isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      if (spin == 0) {
        is1 = Def::Tpow[2 * isite1 - 2];
      }
      else {
        is1 = Def::Tpow[2 * isite1 - 1];
      }
      ibit1 = (long int)myrank & is1;
      num1 = ibit1 / is1;
#pragma omp parallel for default(none) shared(list_Diagonal) \
                     firstprivate(i_max, dtmp_V, num1) private(j)
      for (j = 1; j <= i_max; j++) list_Diagonal[j] += num1 * dtmp_V;

      break;/*case HubbardGC, case KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:

      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        ibit1_up = (((long int)myrank& is1_up) / is1_up) ^ (1 - spin);
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, ibit1_up) private(j)
        for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V * ibit1_up;
      } /*if (Def::iFlgGeneralSpin == FALSE)*/
      else /*if (Def::iFlgGeneralSpin == TRUE)*/ {
        num1 = BitCheckGeneral((long int)myrank,
          isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V) private(j)
          for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V;
        }/*if (num1 != 0)*/
      }/*if (Def::iFlgGeneralSpin == TRUE)*/
      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    } /*switch (Def::iCalcModel)*/

    return 0;

  }/*if (isite1 >= Def::Nsite*/

  switch (Def::iCalcModel) {
  case HubbardGC:
    if (spin == 0) {
      is1 = Def::Tpow[2 * isite1 - 2];
    }
    else {
      is1 = Def::Tpow[2 * isite1 - 1];
    }
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
    for (j = 1; j <= i_max; j++) {

      ibit1 = (j - 1)&is1;
      num1 = ibit1 / is1;
      //fprintf(stdoutMPI, "DEBUG: spin=%ld  is1=%ld: isite1=%ld j=%ld num1=%ld \n",spin,is1,isite1,j,num1);

      list_Diagonal[j] += num1 * dtmp_V;
    }
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    if (spin == 0) {
      is1 = Def::Tpow[2 * isite1 - 2];
    }
    else {
      is1 = Def::Tpow[2 * isite1 - 1];
    }

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
    for (j = 1; j <= i_max; j++) {

      ibit1 = list_1[j] & is1;
      num1 = ibit1 / is1;
      list_Diagonal[j] += num1 * dtmp_V;
    }
    break;

  case SpinGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
      for (j = 1; j <= i_max; j++) {
        ibit1_up = (((j - 1)& is1_up) / is1_up) ^ (1 - spin);
        list_Diagonal[j] += dtmp_V * ibit1_up;
      }
    }
    else {
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(j - 1, isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          list_Diagonal[j] += dtmp_V;
        }
      }
    }
    break;

  case Spin:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
      for (j = 1; j <= i_max; j++) {
        ibit1_up = ((list_1[j] & is1_up) / is1_up) ^ (1 - spin);
        list_Diagonal[j] += dtmp_V * ibit1_up;
      }
    }
    else {
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1) 
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(list_1[j], isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          list_Diagonal[j] += dtmp_V;
        }
      }
    }

    break;
  default:
    fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
    return -1;
  }

  return 0;
}
/**
 * @brief Calculate the components for Coulombinter interaction, \f$ V_{ij} n_ {i}n_{j} \f$
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param dtmp_V [in] A value of coulombinter interaction \f$ V_{ij} \f$
 * @param X [in] Define list to get the operator information.
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalCoulombInter
(
  long int isite1,
  long int isite2,
  double dtmp_V
) {
  long int is1_up, is1_down;
  long int ibit1_up, ibit1_down;
  long int num1;
  long int is2_up, is2_down;
  long int ibit2_up, ibit2_down;
  long int num2;

  long int j;
  long int i_max = Check::idim_max;

  /*
   Force isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
  }/*if (isite2 < isite1)*/
  /*
    When isite1 & site2 are in the inter process region
  */
  if (/*isite2 => */ isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

      num1 = 0;
      num2 = 0;

      ibit1_up = (long int)myrank&is1_up;
      num1 += ibit1_up / is1_up;
      ibit1_down = (long int)myrank&is1_down;
      num1 += ibit1_down / is1_down;

      ibit2_up = (long int)myrank&is2_up;
      num2 += ibit2_up / is2_up;
      ibit2_down = (long int)myrank&is2_down;
      num2 += ibit2_down / is2_down;

#pragma omp parallel for default(none) shared(list_Diagonal) \
      firstprivate(i_max, dtmp_V, num1, num2) private(j)
      for (j = 1; j <= i_max; j++) list_Diagonal[j] += num1 * num2*dtmp_V;

      break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case Spin:
    case SpinGC:
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V)
      for (j = 1; j <= i_max; j++) {
        list_Diagonal[j] += dtmp_V;
      }
      break;/*case Spin, SpinGC*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    }/*switch (Def::iCalcModel)*/

    return 0;

  }/*if (isite1 > Def::Nsite)*/
  else if (isite2 > Def::Nsite /* => isite1 */) {

    switch (Def::iCalcModel) {
    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];
      num2 = 0;
      ibit2_up = (long int)myrank&is2_up;
      num2 += ibit2_up / is2_up;
      ibit2_down = (long int)myrank&is2_down;
      num2 += ibit2_down / is2_down;
      break;

    case Spin:
    case SpinGC:
      break;

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
    }

    switch (Def::iCalcModel) {

    case HubbardGC:

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, num2, is1_up, is1_down) \
private(num1, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_up = (j - 1)&is1_up;
        num1 += ibit1_up / is1_up;
        ibit1_down = (j - 1)&is1_down;
        num1 += ibit1_down / is1_down;

        list_Diagonal[j] += num1 * num2*dtmp_V;
      }

      break;/*case HubbardGC*/

    case KondoGC:
    case Hubbard:
    case Kondo:

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, is1_down, num2) \
private(num1, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_up = list_1[j] & is1_up;
        num1 += ibit1_up / is1_up;
        ibit1_down = list_1[j] & is1_down;
        num1 += ibit1_down / is1_down;

        list_Diagonal[j] += num1 * num2*dtmp_V;
      }
      break;/*case KondoGC, Hubbard, Kondo:*/

    case Spin:
    case SpinGC:
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V)
      for (j = 1; j <= i_max; j++) {
        list_Diagonal[j] += dtmp_V;
      }
      break;/* case Spin, SpinGC:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    }/*switch (Def::iCalcModel)*/

    return 0;

  }/*else if (isite2 > Def::Nsite)*/
  else {
    switch (Def::iCalcModel) {
    case HubbardGC: //list_1[j] -> j-1
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];
#pragma omp parallel for default(none) shared( list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1, ibit1_up, ibit1_down, num2, ibit2_up, ibit2_down)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        num2 = 0;
        ibit1_up = (j - 1)&is1_up;
        num1 += ibit1_up / is1_up;
        ibit1_down = (j - 1)&is1_down;
        num1 += ibit1_down / is1_down;

        ibit2_up = (j - 1)&is2_up;
        num2 += ibit2_up / is2_up;
        ibit2_down = (j - 1)&is2_down;
        num2 += ibit2_down / is2_down;

        list_Diagonal[j] += num1 * num2*dtmp_V;
      }
      break;
    case KondoGC:
    case Hubbard:
    case Kondo:
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1, ibit1_up, ibit1_down, num2, ibit2_up, ibit2_down)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        num2 = 0;
        ibit1_up = list_1[j] & is1_up;
        num1 += ibit1_up / is1_up;
        ibit1_down = list_1[j] & is1_down;
        num1 += ibit1_down / is1_down;

        ibit2_up = list_1[j] & is2_up;
        num2 += ibit2_up / is2_up;
        ibit2_down = list_1[j] & is2_down;
        num2 += ibit2_down / is2_down;

        list_Diagonal[j] += num1 * num2*dtmp_V;
      }
      break;

    case Spin:
    case SpinGC:
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V)
      for (j = 1; j <= i_max; j++) {
        list_Diagonal[j] += dtmp_V;
      }
      break;
    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
    }
  }

  return 0;
}
/**
 * @brief Calculate the components for Hund interaction, \f$ H_{ij}(n_ {i\uparrow}n_{j\uparrow}+ n_ {i\downarrow}n_{j\downarrow})\f$
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param dtmp_V [in] A value of Hund interaction \f$ H_{ij} \f$
 * @param X [in] Define list to get the operator information.
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalHund
(
  long int isite1,
  long int isite2,
  double dtmp_V
) {

  long int is1_up, is1_down;
  long int ibit1_up, ibit1_down;
  long int num1_up, num1_down;
  long int is2_up, is2_down;
  long int ibit2_up, ibit2_down;
  long int num2_up, num2_down;

  long int is_up;
  long int ibit;
  long int j;
  long int i_max = Check::idim_max;
  /*
  Force isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
  }
  /*
  When isite1 & site2 are in the inter process region
  */
  if (/*isite2 >= */ isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

      num1_up = 0;
      num1_down = 0;
      num2_up = 0;
      num2_down = 0;

      ibit1_up = (long int)myrank &is1_up;
      num1_up = ibit1_up / is1_up;
      ibit1_down = (long int)myrank &is1_down;
      num1_down = ibit1_down / is1_down;

      ibit2_up = (long int)myrank &is2_up;
      num2_up = ibit2_up / is2_up;
      ibit2_down = (long int)myrank &is2_down;
      num2_down = ibit2_down / is2_down;

#pragma omp parallel for default(none) shared(list_Diagonal) \
  firstprivate(i_max, dtmp_V, num1_up, num1_down, num2_up, num2_down) private(j)
      for (j = 1; j <= i_max; j++)
        list_Diagonal[j] += dtmp_V * (num1_up*num2_up + num1_down * num2_down);

      break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:

      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
      is_up = is1_up + is2_up;
      ibit = (long int)myrank & is_up;
      if (ibit == 0 || ibit == is_up) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V) private(j) 
        for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V;
      }
      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
    }

    return 0;

  }/*if (isite1 > Def::Nsite)*/
  else if (isite2 > Def::Nsite /* >= isite1 */) {

    switch (Def::iCalcModel) {

    case HubbardGC:

      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

      num2_up = 0;
      num2_down = 0;

      ibit2_up = (long int)myrank &is2_up;
      num2_up = ibit2_up / is2_up;
      ibit2_down = (long int)myrank &is2_down;
      num2_down = ibit2_down / is2_down;

#pragma omp parallel for default(none) shared( list_Diagonal) \
firstprivate(i_max, dtmp_V, num2_up, num2_down, is1_up, is1_down) \
private(num1_up, num1_down, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1_up = 0;
        num1_down = 0;

        ibit1_up = (j - 1)&is1_up;
        num1_up = ibit1_up / is1_up;
        ibit1_down = (j - 1)&is1_down;
        num1_down = ibit1_down / is1_down;

        list_Diagonal[j] += dtmp_V * (num1_up*num2_up + num1_down * num2_down);
      }
      break;/*case HubbardGC:*/

    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

      num2_up = 0;
      num2_down = 0;

      ibit2_up = (long int)myrank&is2_up;
      num2_up = ibit2_up / is2_up;
      ibit2_down = (long int)myrank&is2_down;
      num2_down = ibit2_down / is2_down;

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, num2_up, num2_down, is1_up, is1_down) \
private(num1_up, num1_down, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1_up = 0;
        num1_down = 0;

        ibit1_up = list_1[j] & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit1_down = list_1[j] & is1_down;
        num1_down = ibit1_down / is1_down;

        list_Diagonal[j] += dtmp_V * (num1_up*num2_up + num1_down * num2_down);
      }
      break;/*case KondoGC, Hubbard, Kondo:*/

    case SpinGC:
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
      ibit2_up = (long int)myrank & is2_up;

      if (ibit2_up == is2_up) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = (j - 1) & is1_up;
          if (ibit1_up == is1_up) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      else if (ibit2_up == 0) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = (j - 1) & is1_up;
          if (ibit1_up == 0) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      break;/*case SpinGC:*/

    case Spin:
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
      ibit2_up = (long int)myrank & is2_up;

      if (ibit2_up == is2_up) {
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = list_1[j] & is1_up;
          if (ibit1_up == is1_up) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      else if (ibit2_up == 0) {
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = list_1[j] & is1_up;
          if (ibit1_up == 0) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      break;/*case Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    }/*switch (Def::iCalcModel)*/

    return 0;

  }/*else if (isite2 > Def::Nsite)*/
  else {
    switch (Def::iCalcModel) {
    case HubbardGC: // list_1[j] -> j-1
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

#pragma omp parallel for default(none) shared( list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1_up, num1_down, num2_up, num2_down, ibit1_up, ibit1_down, ibit2_up, ibit2_down)
      for (j = 1; j <= i_max; j++) {
        num1_up = 0;
        num1_down = 0;
        num2_up = 0;
        num2_down = 0;

        ibit1_up = (j - 1)&is1_up;
        num1_up = ibit1_up / is1_up;
        ibit1_down = (j - 1)&is1_down;
        num1_down = ibit1_down / is1_down;

        ibit2_up = (j - 1)&is2_up;
        num2_up = ibit2_up / is2_up;
        ibit2_down = (j - 1)&is2_down;
        num2_down = ibit2_down / is2_down;

        list_Diagonal[j] += dtmp_V * (num1_up*num2_up + num1_down * num2_down);
      }
      break;
    case KondoGC:
    case Hubbard:
    case Kondo:
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1_up, num1_down, num2_up, num2_down, ibit1_up, ibit1_down, ibit2_up, ibit2_down)
      for (j = 1; j <= i_max; j++) {
        num1_up = 0;
        num1_down = 0;
        num2_up = 0;
        num2_down = 0;

        ibit1_up = list_1[j] & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit1_down = list_1[j] & is1_down;
        num1_down = ibit1_down / is1_down;

        ibit2_up = list_1[j] & is2_up;
        num2_up = ibit2_up / is2_up;
        ibit2_down = list_1[j] & is2_down;
        num2_down = ibit2_down / is2_down;

        list_Diagonal[j] += dtmp_V * (num1_up*num2_up + num1_down * num2_down);
      }
      break;

    case SpinGC:
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
      is_up = is1_up + is2_up;
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, is_up) private(j, ibit) 
      for (j = 1; j <= i_max; j++) {
        ibit = (j - 1) & is_up;
        if (ibit == 0 || ibit == is_up) {
          list_Diagonal[j] += dtmp_V;
        }
      }
      break;

    case Spin:
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
      is_up = is1_up + is2_up;
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, is_up) private(j, ibit) 
      for (j = 1; j <= i_max; j++) {
        ibit = list_1[j] & is_up;
        if (ibit == 0 || ibit == is_up) {
          list_Diagonal[j] += dtmp_V;
        }
      }
      break;
    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;
    }
  }
  return 0;
}
/**
 * @brief Calculate the components for general two-body diagonal interaction, \f$ H_{i\sigma_1 j\sigma_2} n_ {i\sigma_1}n_{j\sigma_2}\f$
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param isigma1 [in] a spin index at \f$i \f$ site.
 * @param isigma2 [in] a spin index at \f$j \f$ site.
 * @param dtmp_V [in] A value of general two-body diagonal interaction \f$ H_{i\sigma_1 j\sigma_2} \f$
 * @param X [in] Define list to get the operator information.
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalInterAll
(
  long int isite1,
  long int isite2,
  long int isigma1,
  long int isigma2,
  double dtmp_V
)
{
  long int is1_spin;
  long int is2_spin;
  long int is1_up;
  long int is2_up;

  long int ibit1_spin;
  long int ibit2_spin;

  long int num1;
  long int num2;

  long int j;
  long int i_max = Check::idim_max;

  /*
  Forse isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
    j = isigma2;
    isigma2 = isigma1;
    isigma1 = j;
  }
  /*
  When isite1 & site2 are in the inter process regino
  */
  if (isite1 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];

      num1 = 0;
      ibit1_spin = (long int)myrank&is1_spin;
      num1 += ibit1_spin / is1_spin;

      num2 = 0;
      ibit2_spin = (long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, num2, num1) private(ibit1_spin, j)
      for (j = 1; j <= i_max; j++) list_Diagonal[j] += num1 * num2*dtmp_V;

      break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:

      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        is2_up = Def::Tpow[isite2 - 1];
        num1 = X_SpinGC_CisAis((long int)myrank + 1, is1_up, isigma1);
        num2 = X_SpinGC_CisAis((long int)myrank + 1, is2_up, isigma2);

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, num1, num2) private(j)
        for (j = 1; j <= i_max; j++) {
          list_Diagonal[j] += num1 * num2*dtmp_V;
        }
      }/*if (Def::iFlgGeneralSpin == FALSE)*/
      else {//start:generalspin
        num1 = BitCheckGeneral((long int)myrank, isite1, isigma1,
          Def::SiteToBit, Def::Tpow);
        num2 = BitCheckGeneral((long int)myrank, isite2, isigma2,
          Def::SiteToBit, Def::Tpow);
        if (num1 != 0 && num2 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, num1, X) private(j)
          for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V * num1;
        }
      }/*if (Def::iFlgGeneralSpin == TRUE)*/

      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    }/*if (isite1 > Def::Nsite)*/

    return 0;

  }/*if (isite1 > Def::Nsite)*/
  else if (isite2 > Def::Nsite) {

    switch (Def::iCalcModel) {

    case HubbardGC:

      is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];

      num2 = 0;
      ibit2_spin = (long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_spin, num2) private(num1, ibit1_spin, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_spin = (j - 1)&is1_spin;
        num1 += ibit1_spin / is1_spin;
        list_Diagonal[j] += num1 * num2*dtmp_V;
      }
      break;/*case HubbardGC:*/

    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];

      num2 = 0;
      ibit2_spin = (long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;

#pragma omp parallel for default(none) shared(list_Diagonal, list_1) \
firstprivate(i_max, dtmp_V, is1_spin, num2) private(num1, ibit1_spin, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_spin = list_1[j] & is1_spin;
        num1 += ibit1_spin / is1_spin;
        list_Diagonal[j] += num1 * num2*dtmp_V;
      }
      break;/*case KondoGC, Hubbard, Kondo:*/

    case SpinGC:

      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        is2_up = Def::Tpow[isite2 - 1];
        num2 = X_SpinGC_CisAis((long int)myrank + 1, is2_up, isigma2);

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, num2) private(j, num1)
        for (j = 1; j <= i_max; j++) {
          num1 = X_SpinGC_CisAis(j, is1_up, isigma1);
          list_Diagonal[j] += num1 * num2*dtmp_V;
        }
      }/* if (Def::iFlgGeneralSpin == FALSE)*/
      else {//start:generalspin
        num2 = BitCheckGeneral((long int)myrank, isite2, isigma2,
          Def::SiteToBit, Def::Tpow);
        if (num2 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(j - 1, isite1, isigma1, Def::SiteToBit, Def::Tpow);
            list_Diagonal[j] += dtmp_V * num1;
          }
        }
      }/* if (Def::iFlgGeneralSpin == TRUE)*/

      break;/*case SpinGC:*/

    case Spin:

      if (Def::iFlgGeneralSpin == FALSE) {
        is1_up = Def::Tpow[isite1 - 1];
        is2_up = Def::Tpow[isite2 - 1];
        num2 = X_SpinGC_CisAis((long int)myrank + 1, is2_up, isigma2);

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, num2) private(j, num1)
        for (j = 1; j <= i_max; j++) {
          num1 = X_Spin_CisAis(j, is1_up, isigma1);
          list_Diagonal[j] += num1 * num2*dtmp_V;
        }
      }/* if (Def::iFlgGeneralSpin == FALSE)*/
      else /* if (Def::iFlgGeneralSpin == TRUE)*/ {
        num2 = BitCheckGeneral((long int)myrank, isite2, isigma2, \
          Def::SiteToBit, Def::Tpow);
        if (num2 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) \
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(list_1[j], isite1, isigma1, Def::SiteToBit, Def::Tpow);
            list_Diagonal[j] += dtmp_V * num1;
          }
        }
      } /* if (Def::iFlgGeneralSpin == TRUE)*/

      break;/*case Spin:*/

    default:
      fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
      return -1;

    }/*switch (Def::iCalcModel)*/

    return 0;

  }/*else if (isite2 > Def::Nsite)*/

  switch (Def::iCalcModel) {
  case HubbardGC: //list_1[j] -> j-1
    is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
    is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) private(num1, ibit1_spin, num2, ibit2_spin)
    for (j = 1; j <= i_max; j++) {
      num1 = 0;
      num2 = 0;
      ibit1_spin = (j - 1)&is1_spin;
      num1 += ibit1_spin / is1_spin;
      ibit2_spin = (j - 1)&is2_spin;
      num2 += ibit2_spin / is2_spin;
      list_Diagonal[j] += num1 * num2*dtmp_V;
    }
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    is1_spin = Def::Tpow[2 * isite1 - 2 + isigma1];
    is2_spin = Def::Tpow[2 * isite2 - 2 + isigma2];

#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) private(num1, ibit1_spin, num2, ibit2_spin)
    for (j = 1; j <= i_max; j++) {
      num1 = 0;
      num2 = 0;
      ibit1_spin = list_1[j] & is1_spin;
      num1 += ibit1_spin / is1_spin;

      ibit2_spin = list_1[j] & is2_spin;
      num2 += ibit2_spin / is2_spin;
      list_Diagonal[j] += num1 * num2*dtmp_V;
    }
    break;

  case Spin:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
      for (j = 1; j <= i_max; j++) {
        num1 = X_Spin_CisAis(j, is1_up, isigma1);
        num2 = X_Spin_CisAis(j, is2_up, isigma2);
        list_Diagonal[j] += num1 * num2*dtmp_V;
      }
    }
    else {
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) private(j, num1)
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(list_1[j], isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          num1 = BitCheckGeneral(list_1[j], isite2, isigma2, Def::SiteToBit, Def::Tpow);
          list_Diagonal[j] += dtmp_V * num1;
        }
      }

    }
    break;

  case SpinGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      is1_up = Def::Tpow[isite1 - 1];
      is2_up = Def::Tpow[isite2 - 1];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
      for (j = 1; j <= i_max; j++) {
        num1 = X_SpinGC_CisAis(j, is1_up, isigma1);
        num2 = X_SpinGC_CisAis(j, is2_up, isigma2);
        list_Diagonal[j] += num1 * num2*dtmp_V;
      }
    }
    else {//start:generalspin
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) private(j, num1)
      for (j = 1; j <= i_max; j++) {
        num1 = BitCheckGeneral(j - 1, isite1, isigma1, Def::SiteToBit, Def::Tpow);
        if (num1 != 0) {
          num1 = BitCheckGeneral(j - 1, isite2, isigma2, Def::SiteToBit, Def::Tpow);
          list_Diagonal[j] += dtmp_V * num1;
        }
      }
    }
    break;

  default:
    fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
    return -1;
  }
  return 0;
}
/**
 * @brief Calculate diagonal components and obtain the list, list_diagonal.
 *
 * @param X [in] Struct to get the information of the diagonal operators.
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @retval -1 fail to calculate diagonal components.
 * @retval 0 succeed to calculate diagonal components.
 */
int diagonalcalc
(
  
) {

  FILE *fp;
  long int i, j;
  long int isite1, isite2;
  long int spin;
  double tmp_V;

  /*[s] For InterAll*/
  long int A_spin, B_spin;
  /*[e] For InterAll*/
  long int i_max = Check::idim_max;

  fprintf(stdoutMPI, "%s", "  Start: Calculate diagaonal components of Hamiltonian. \n");
  TimeKeeper("%s_TimeKeeper.dat", "diagonal calculation starts:   %s", "a");

#pragma omp parallel for default(none) private(j) shared(list_Diagonal) firstprivate(i_max)
  for (j = 1; j <= i_max; j++) {
    list_Diagonal[j] = 0.0;
  }

  if (Def::NCoulombIntra > 0) {
    if (childfopenMPI("CHECK_CoulombIntra.dat", "w", &fp) != 0) {
      return -1;
    }
    for (i = 0; i < Def::NCoulombIntra; i++) {
      isite1 = Def::CoulombIntra[i][0] + 1;
      tmp_V = Def::ParaCoulombIntra[i];
      fprintf(fp, "i=%ld isite1=%ld tmp_V=%lf \n", i, isite1, tmp_V);
      SetDiagonalCoulombIntra(isite1, tmp_V);
    }
    fclose(fp);
  }

  if (Def::EDNChemi > 0) {
    if (childfopenMPI("CHECK_Chemi.dat", "w", &fp) != 0) {
      return -1;
    }
    for (i = 0; i < Def::EDNChemi; i++) {
      isite1 = Def::EDChemi[i] + 1;
      spin = Def::EDSpinChemi[i];
      tmp_V = -Def::EDParaChemi[i];
      fprintf(fp, "i=%ld spin=%ld isite1=%ld tmp_V=%lf \n", i, spin, isite1, tmp_V);
      if (SetDiagonalChemi(isite1, tmp_V, spin) != 0) {
        return -1;
      }
    }
    fclose(fp);
  }

  if (Def::NCoulombInter > 0) {
    if (childfopenMPI("CHECK_INTER_U.dat", "w", &fp) != 0) {
      return -1;
    }
    for (i = 0; i < Def::NCoulombInter; i++) {
      isite1 = Def::CoulombInter[i][0] + 1;
      isite2 = Def::CoulombInter[i][1] + 1;
      tmp_V = Def::ParaCoulombInter[i];
      fprintf(fp, "i=%ld isite1=%ld isite2=%ld tmp_V=%lf \n", i, isite1, isite2, tmp_V);
      if (SetDiagonalCoulombInter(isite1, isite2, tmp_V) != 0) {
        return -1;
      }
    }
    fclose(fp);
  }
  if (Def::NHundCoupling > 0) {
    if (childfopenMPI("CHECK_Hund.dat", "w", &fp) != 0) {
      return -1;
    }
    for (i = 0; i < Def::NHundCoupling; i++) {
      isite1 = Def::HundCoupling[i][0] + 1;
      isite2 = Def::HundCoupling[i][1] + 1;
      tmp_V = -Def::ParaHundCoupling[i];
      if (SetDiagonalHund(isite1, isite2, tmp_V) != 0) {
        return -1;
      }
      fprintf(fp, "i=%ld isite1=%ld isite2=%ld tmp_V=%lf \n", i, isite1, isite2, tmp_V);
    }
    fclose(fp);
  }

  if (Def::NInterAll_Diagonal > 0) {
    if (childfopenMPI("CHECK_InterAll.dat", "w", &fp) != 0) {
      return -1;
    }
    for (i = 0; i < Def::NInterAll_Diagonal; i++) {
      isite1 = Def::InterAll_Diagonal[i][0] + 1;
      A_spin = Def::InterAll_Diagonal[i][1];
      isite2 = Def::InterAll_Diagonal[i][2] + 1;
      B_spin = Def::InterAll_Diagonal[i][3];
      tmp_V = Def::ParaInterAll_Diagonal[i];
      fprintf(fp, "i=%ld isite1=%ld A_spin=%ld isite2=%ld B_spin=%ld tmp_V=%lf \n", i, isite1, A_spin, isite2, B_spin, tmp_V);
      SetDiagonalInterAll(isite1, isite2, A_spin, B_spin, tmp_V);
    }
    fclose(fp);
  }

  TimeKeeper("%s_TimeKeeper.dat", "diagonal calculation finishes: %s", "a");
  fprintf(stdoutMPI, "%s", "  End  : Calculate diagaonal components of Hamiltonian. \n\n");
  return 0;
}
