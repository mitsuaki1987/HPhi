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

#include "mltply.hpp"
#include "FileIO.hpp"
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "mltplyHubbard.hpp"
#include "mltplyHubbardCore.hpp"
#include "mltplySpinCore.hpp"
#include "mltplyMPIHubbard.hpp"
#include "mltplyMPISpinCore.hpp"
#include "common/setmemory.hpp"
#include "mltplyCommon.hpp"

/**
 * @file
 * 
 * @brief  File for calculation of one body green's function
 *
 * @version 0.1, 0.2
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
/**
 * @brief function of calculation for one body green's function for Hubbard GC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_HubbardGC(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
  std::complex<double> **prod
){
  long int i;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int i_max;
  long int ibit;
  long int is;
  std::complex<double> tmp_OneGreen = 1.0;
  int complex_conj, istate;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    complex_conj = 0;
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];
    if (org_isite1 > X->Def.Nsite &&
      org_isite2 > X->Def.Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        if (org_sigma1 == 0) {
          is = X->Def.Tpow[2 * org_isite1 - 2];
        }
        else {
          is = X->Def.Tpow[2 * org_isite1 - 1];
        }
        ibit = (long int)myrank & is;
        if (ibit == is) {
          zaxpy_long(i_max*nstate, tmp_OneGreen, &vec[1][0], &Xvec[1][0]);
        }
      }
      else {
        X_GC_child_general_hopp_MPIdouble(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite2 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_isite1 < org_isite2) {
        X_GC_child_general_hopp_MPIsingle(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
      else {
        X_GC_child_general_hopp_MPIsingle(org_isite2 - 1, org_sigma2, org_isite1 - 1, org_sigma1, 
          -tmp_OneGreen, X, nstate, Xvec, vec);
        complex_conj = 1;
      }
    }
    else {
      if (child_general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2) != 0) {
        return -1;
      }
      GC_child_general_hopp(nstate, Xvec, vec, X, tmp_OneGreen);
    }

    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
    if (complex_conj == 1) 
      for (istate = 0; istate < nstate; istate++) prod[i][istate] = conj(prod[i][istate]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for Hubbard model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_Hubbard(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
  std::complex<double> **prod
) {
  long int i, j;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int i_max;
  int num1, one = 1, complex_conj, istate;
  long int ibit;
  long int is;
  std::complex<double> tmp_OneGreen = 1.0, dmv;

  i_max = X->Check.idim_max;
  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    complex_conj = 0;
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];

    if (X->Def.iFlgSzConserved == TRUE) {
      if (org_sigma1 != org_sigma2) {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (X->Def.iCalcModel == Kondo || X->Def.iCalcModel == KondoGC) {
      if ((X->Def.LocSpn[org_isite1 - 1] == 1 && X->Def.LocSpn[org_isite2 - 1] == 0) ||
          (X->Def.LocSpn[org_isite1 - 1] == 0 && X->Def.LocSpn[org_isite2 - 1] == 1)
        )
      {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (org_isite1 > X->Def.Nsite &&
      org_isite2 > X->Def.Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {//diagonal
        is = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];
        ibit = (long int)myrank & is;
        if (ibit == is) {
          zaxpy_long(i_max*nstate, tmp_OneGreen, &vec[1][0], &Xvec[1][0]);
        }
      }
      else {
        X_child_general_hopp_MPIdouble(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2, 
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite2 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_isite1 < org_isite2) {
        X_child_general_hopp_MPIsingle(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
      else {
        X_child_general_hopp_MPIsingle(org_isite2 - 1, org_sigma2, org_isite1 - 1, org_sigma1, 
          -tmp_OneGreen, X, nstate, Xvec, vec);
        complex_conj = 1;
      }
    }
    else {
      if (child_general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2) != 0) {
        return -1;
      }
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        is = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];

#pragma omp parallel for default(none) shared(list_1, vec,Xvec,nstate,one,tmp_OneGreen) \
firstprivate(i_max, is) private(num1, ibit, dmv)
        for (j = 1; j <= i_max; j++) {
          ibit = list_1[j] & is;
          num1 = ibit / is;
          dmv = (std::complex<double>)num1;
          zaxpy_(&nstate, &dmv, vec[j], &one, Xvec[j], &one);
        }
      }
      else {
        child_general_hopp(nstate, Xvec, vec, X, tmp_OneGreen);
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
    if (complex_conj == 1)
      for (istate = 0; istate < nstate; istate++) prod[i][istate] = conj(prod[i][istate]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for Half-Spin model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinHalf(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
  std::complex<double> **prod
) {
  long int i, j;
  long int isite1;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  std::complex<double> dmv;
  long int i_max;
  long int ibit1;
  long int is1_up;
  int one = 1;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];

    if (org_sigma1 == org_sigma2) {
      if (org_isite1 == org_isite2) {
        if (org_isite1 > X->Def.Nsite) {
          is1_up = X->Def.Tpow[org_isite1 - 1];
          ibit1 = X_SpinGC_CisAis((long int)myrank + 1, is1_up, org_sigma1);
          if (ibit1 != 0) {
            zaxpy_long(i_max*nstate, 1.0, &vec[1][0], &Xvec[1][0]);
          }
        }// org_isite1 > X->Def.Nsite
        else {
          isite1 = X->Def.Tpow[org_isite1 - 1];
#pragma omp parallel for default(none) private(j,dmv)  \
  firstprivate(i_max, isite1, org_sigma1, X) shared(vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            dmv = X_Spin_CisAis(j, isite1, org_sigma1);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for General-Spin model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGeneral(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
  std::complex<double> **prod
) {
  long int i, j;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  std::complex<double> dmv;
  long int i_max;
  int num1, one = 1;
  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];

    if (org_isite1 == org_isite2) {
      if (org_isite1 > X->Def.Nsite) {
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          num1 = BitCheckGeneral((long int)myrank,
                                           org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
          if (num1 != 0) {
            zaxpy_long(i_max*nstate, 1.0, &vec[1][0], &Xvec[1][0]);
          }
        }
      }
      else {//org_isite1 <= X->Def.Nsite
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1,dmv) \
  firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec,Xvec, list_1,nstate,one)
          for (j = 1; j <= i_max; j++) {
            dmv = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for Half-SpinGC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGCHalf(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
    std::complex<double> **prod
) {
  long int i, j;
  long int isite1;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  std::complex<double> dmv;
  long int i_max;
  int tmp_sgn, one = 1;
  long int tmp_off = 0;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];

    if (org_isite1 == org_isite2) {
      if (org_isite1 > X->Def.Nsite) {
        if (org_sigma1 == org_sigma2) {  // longitudinal magnetic field
          X_GC_child_CisAis_spin_MPIdouble(org_isite1 - 1, org_sigma1, 1.0, X, nstate, Xvec, vec);
        }
        else {  // transverse magnetic field
          X_GC_child_CisAit_spin_MPIdouble(org_isite1 - 1, org_sigma1, org_sigma2, 1.0, X, nstate, Xvec, vec);
        }
      }
      else {
        isite1 = X->Def.Tpow[org_isite1 - 1];

        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn,dmv) \
  firstprivate(i_max, isite1, org_sigma1, X) shared(vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            dmv = X_SpinGC_CisAis(j, isite1, org_sigma1);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
        else {
          // transverse magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off,dmv) \
  firstprivate(i_max, isite1, org_sigma2, X) shared(vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = X_SpinGC_CisAit(j, isite1, org_sigma2, &tmp_off);
            if (tmp_sgn != 0) {
              dmv = (std::complex<double>)tmp_sgn;
              zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for General SpinGC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGCGeneral(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec,
  std::complex<double> **prod
) {
  long int i, j;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  std::complex<double> dmv;
  long int i_max;
  long int tmp_off = 0;
  int num1, one = 1;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];
    if (org_isite1 == org_isite2) {
      if (org_isite1 > X->Def.Nsite) {
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          X_GC_child_CisAis_GeneralSpin_MPIdouble(org_isite1 - 1, org_sigma1,
            1.0, X, nstate, Xvec, vec);
        }
        else {
          // transverse magnetic field
          X_GC_child_CisAit_GeneralSpin_MPIdouble(
            org_isite1 - 1, org_sigma1, org_sigma2, 1.0, X, nstate, Xvec, vec);
        }
      }
      else {//org_isite1 <= X->Def.Nsite
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1,dmv) \
  firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
            dmv = (std::complex<double>)num1;
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
        else {
          // transverse magnetic field
#pragma omp parallel for default(none) private(j, num1,dmv) \
  firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X,tmp_off) shared(vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            num1 = GetOffCompGeneralSpin(
              j - 1, org_isite1, org_sigma2, org_sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != 0) {
              dmv = (std::complex<double>)num1;
              zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for Spin model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_Spin(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec,
  std::complex<double> **prod
) {
  int info = 0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    info = expec_cisajs_SpinHalf(X, nstate, Xvec, vec, prod);
  }
  else {
    info = expec_cisajs_SpinGeneral(X, nstate, Xvec, vec, prod);
  }
  return info;
}

/**
 * @brief function of calculation for one body green's function for SpinGC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGC(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec,
  std::complex<double> **prod
) {
  int info = 0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    info = expec_cisajs_SpinGCHalf(X, nstate, Xvec, vec, prod);
  }
  else {
    info = expec_cisajs_SpinGCGeneral(X, nstate, Xvec, vec, prod);
  }
  return info;
}
/**
 * @brief function of calculation for one body green's function
 *
 * @param X [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvectors.
 *
 * @version 0.2
 * @details add calculation one body green's functions for general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs(
  struct BindStruct *X,
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec
) {
  FILE *fp;
  char sdt[D_FileNameMax];
  std::complex<double> **prod;
  long int irght, ilft, ihfbit, ica;
  long int i_max;
  //For TPQ
  int step = 0, rand_i = 0, istate;

  if (X->Def.NCisAjt < 1) return 0;

  i_max = X->Check.idim_max;
  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }
  X->Large.i_max = i_max;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_CORR;

  switch (X->Def.iCalcType) {
  case TPQCalc:
    step = X->Def.istep;
    TimeKeeperWithRandAndStep(X, "%s_TimeKeeper.dat", "set %d step %d:expec_cisajs begins: %s", "a", 0, step);
    break;
  case TimeEvolution:
    step = X->Def.istep;
    TimeKeeperWithStep(X, "%s_TimeKeeper.dat", "step %d:expec_cisajs begins: %s", "a", step);
    break;
  case FullDiag:
  case CG:
    break;
  }

  prod = cd_2d_allocate(X->Def.NCisAjt, nstate);
  switch (X->Def.iCalcModel) {
  case HubbardGC:
    if (expec_cisajs_HubbardGC(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  case KondoGC:
  case Hubbard:
  case Kondo:
    if (expec_cisajs_Hubbard(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  case Spin: // for the Sz-conserved spin system
    if (expec_cisajs_Spin(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  case SpinGC:
    if (expec_cisajs_SpinGC(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  default:
    return -1;
  }

  for (istate = 0; istate < nstate; istate++) {
    switch (X->Def.iCalcType) {
    case TPQCalc:
      step = X->Def.istep;
      sprintf(sdt, "%s_cisajs_set%dstep%d.dat", X->Def.CDataFileHead, istate, step);
      break;
    case TimeEvolution:
      step = X->Def.istep;
      sprintf(sdt, "%s_cisajs_step%d.dat", X->Def.CDataFileHead, step);
      break;
    case FullDiag:
    case CG:
      sprintf(sdt, "%s_cisajs_eigen%d.dat", X->Def.CDataFileHead, istate);
      break;
    }
    if (childfopenMPI(sdt, "w", &fp) == 0) {
      for (ica = 0; ica < X->Def.NCisAjt; ica++) {
        fprintf(fp, " %4d %4d %4d %4d %.10lf %.10lf\n",
          X->Def.CisAjt[ica][0], X->Def.CisAjt[ica][1], X->Def.CisAjt[ica][2], X->Def.CisAjt[ica][3],
          real(prod[ica][istate]), imag(prod[ica][istate]));
      }
      fclose(fp);
    }
    else return -1;
  }/*for (istate = 0; istate < nstate; istate++)*/

  if (X->Def.St == 0) {
    if (X->Def.iCalcType == TPQCalc) {
      TimeKeeperWithRandAndStep(X, "%s_TimeKeeper.dat", "set %d step %d:expec_cisajs finishes: %s", "a", rand_i, step);
    }
    else if (X->Def.iCalcType == TimeEvolution) {
      TimeKeeperWithStep(X, "%s_TimeKeeper.dat", "step %d:expec_cisajs finishes: %s", "a", step);
    }
  }
  else if (X->Def.St == 1) {
    TimeKeeper(X, "%s_TimeKeeper.dat", "CG expec_cisajs finishes:     %s", "a");
    fprintf(stdoutMPI, "%s", "  End  : Calculate one body Green functions.\n\n");
  }
  free_cd_2d_allocate(prod);
  return 0;
}
