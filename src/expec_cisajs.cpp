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

#include "expec_cisajs.hpp"
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
#include "DefCommon.hpp"
#include "global.hpp"
#include "log.hpp"
#include <complex>

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
 */
void expec::cisajs::Hubbard::GC(
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

  i_max = Check::idim_max;

  for (i = 0; i < Def::NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    complex_conj = 0;
    org_isite1 = Def::CisAjt[i][0];
    org_isite2 = Def::CisAjt[i][2];
    org_sigma1 = Def::CisAjt[i][1];
    org_sigma2 = Def::CisAjt[i][3];
    if (org_isite1 >= Def::Nsite &&
      org_isite2 >= Def::Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        if (org_sigma1 == 0) {
          is = Def::Tpow[2 * org_isite1];
        }
        else {
          is = Def::Tpow[2 * org_isite1];
        }
        ibit = (long int)MP::myrank & is;
        if (ibit == is) {
          zaxpy_long(i_max*nstate, tmp_OneGreen, &vec[1][0], &Xvec[1][0]);
        }
      }
      else {
        mltply::Hubbard::GC::X_general_hopp_MPIdouble(org_isite1, org_sigma1, org_isite2, org_sigma2,
          -tmp_OneGreen, nstate, Xvec, vec);
      }
    }
    else if (org_isite2 >= Def::Nsite || org_isite1 >= Def::Nsite) {
      if (org_isite1 < org_isite2) {
        mltply::Hubbard::GC::X_general_hopp_MPIsingle(org_isite1, org_sigma1, org_isite2, org_sigma2,
          -tmp_OneGreen, nstate, Xvec, vec);
      }
      else {
        mltply::Hubbard::GC::X_general_hopp_MPIsingle(org_isite2, org_sigma2, org_isite1, org_sigma1, 
          -tmp_OneGreen, nstate, Xvec, vec);
        complex_conj = 1;
      }
    }
    else {
      mltply::Hubbard::general_hopp_GetInfo(org_isite1, org_isite2, org_sigma1, org_sigma2);
      mltply::Hubbard::GC::general_hopp(nstate, Xvec, vec, tmp_OneGreen);
    }

    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
    if (complex_conj == 1) 
      for (istate = 0; istate < nstate; istate++) prod[i][istate] = conj(prod[i][istate]);
  }
}
/**
 * @brief function of calculation for one body green's function for Hubbard model.
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
void expec::cisajs::Hubbard::C(
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

  i_max = Check::idim_max;
  for (i = 0; i < Def::NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    complex_conj = 0;
    org_isite1 = Def::CisAjt[i][0];
    org_isite2 = Def::CisAjt[i][2];
    org_sigma1 = Def::CisAjt[i][1];
    org_sigma2 = Def::CisAjt[i][3];

    if (Def::iFlgSzConserved == TRUE) {
      if (org_sigma1 != org_sigma2) {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (Def::iCalcModel == DC::Kondo || Def::iCalcModel == DC::KondoGC) {
      if ((Def::LocSpn[org_isite1] == 1 && Def::LocSpn[org_isite2] == 0) ||
          (Def::LocSpn[org_isite1] == 0 && Def::LocSpn[org_isite2] == 1)
        )
      {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (org_isite1 >= Def::Nsite &&
      org_isite2 >= Def::Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {//diagonal
        is = Def::Tpow[2 * org_isite1 + org_sigma1];
        ibit = (long int)MP::myrank & is;
        if (ibit == is) {
          zaxpy_long(i_max*nstate, tmp_OneGreen, &vec[1][0], &Xvec[1][0]);
        }
      }
      else {
        mltply::Hubbard::C::X_general_hopp_MPIdouble(org_isite1, org_sigma1, org_isite2, org_sigma2, 
          -tmp_OneGreen, nstate, Xvec, vec);
      }
    }
    else if (org_isite2 >= Def::Nsite || org_isite1 >= Def::Nsite) {
      if (org_isite1 < org_isite2) {
        mltply::Hubbard::C::X_general_hopp_MPIsingle(org_isite1, org_sigma1, org_isite2, org_sigma2,
          -tmp_OneGreen, nstate, Xvec, vec);
      }
      else {
        mltply::Hubbard::C::X_general_hopp_MPIsingle(org_isite2, org_sigma2, org_isite1, org_sigma1, 
          -tmp_OneGreen, nstate, Xvec, vec);
        complex_conj = 1;
      }
    }
    else {
      mltply::Hubbard::general_hopp_GetInfo(org_isite1, org_isite2, org_sigma1, org_sigma2);
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        is = Def::Tpow[2 * org_isite1 + org_sigma1];

#pragma omp parallel for default(none) private(num1, ibit, dmv) \
shared(List::c1, vec,Xvec,nstate,one,tmp_OneGreen, i_max, is)
        for (j = 1; j <= i_max; j++) {
          ibit = List::c1[j] & is;
          num1 = ibit / is;
          dmv = (std::complex<double>)num1;
          zaxpy_(&nstate, &dmv, vec[j], &one, Xvec[j], &one);
        }
      }
      else {
        mltply::Hubbard::C::general_hopp(nstate, Xvec, vec, tmp_OneGreen);
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
    if (complex_conj == 1)
      for (istate = 0; istate < nstate; istate++) prod[i][istate] = conj(prod[i][istate]);
  }
}
/**
 * @brief function of calculation for one body green's function for Half-Spin model.
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
void expec::cisajs::Spin::C::Half(
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

  i_max = Check::idim_max;

  for (i = 0; i < Def::NCisAjt; i++) {
    zclear(i_max * nstate, &Xvec[1][0]);
    org_isite1 = Def::CisAjt[i][0];
    org_isite2 = Def::CisAjt[i][2];
    org_sigma1 = Def::CisAjt[i][1];
    org_sigma2 = Def::CisAjt[i][3];

    if (org_sigma1 == org_sigma2) {
      if (org_isite1 == org_isite2) {
        if (org_isite1 >= Def::Nsite) {
          is1_up = Def::Tpow[org_isite1];
          ibit1 = mltply::Spin::GC::Half::X_CisAis((long int)MP::myrank + 1, is1_up, org_sigma1);
          if (ibit1 != 0) {
            zaxpy_long(i_max*nstate, 1.0, &vec[1][0], &Xvec[1][0]);
          }
        }// org_isite1 >= Def::Nsite
        else {
          isite1 = Def::Tpow[org_isite1];
#pragma omp parallel for default(none) private(j,dmv) \
shared(i_max, isite1, org_sigma1, vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            dmv = mltply::Spin::C::Half::X_CisAis(j, isite1, org_sigma1);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief function of calculation for one body green's function for General-Spin model.
 */
void expec::cisajs::Spin::C::General(
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
  i_max = Check::idim_max;

  for (i = 0; i < Def::NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = Def::CisAjt[i][0];
    org_isite2 = Def::CisAjt[i][2];
    org_sigma1 = Def::CisAjt[i][1];
    org_sigma2 = Def::CisAjt[i][3];

    if (org_isite1 == org_isite2) {
      if (org_isite1 >= Def::Nsite) {
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          num1 = BitCheckGeneral((long int)MP::myrank,
                                           org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
          if (num1 != 0) {
            zaxpy_long(i_max*nstate, 1.0, &vec[1][0], &Xvec[1][0]);
          }
        }
      }
      else {//org_isite1 <= Def::Nsite
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1,dmv) \
shared(i_max, org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow, vec,Xvec, List::c1,nstate,one)
          for (j = 1; j <= i_max; j++) {
            dmv = BitCheckGeneral(List::c1[j], org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief function of calculation for one body green's function for Half-SpinGC model.
 */
void expec::cisajs::Spin::GC::Half(
  int nstate,
  std::complex<double> **Xvec,
  std::complex<double> **vec, 
  std::complex<double>** prod
) {
  long int i, j;
  long int isite1;
  long int org_isite1, org_isite2, org_sigma1, org_sigma2;
  std::complex<double> dmv;
  long int i_max;
  int tmp_sgn, one = 1;
  long int tmp_off = 0;

  i_max = Check::idim_max;

  for (i = 0; i < Def::NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = Def::CisAjt[i][0];
    org_isite2 = Def::CisAjt[i][2];
    org_sigma1 = Def::CisAjt[i][1];
    org_sigma2 = Def::CisAjt[i][3];

    if (org_isite1 == org_isite2) {
      if (org_isite1 >= Def::Nsite) {
        if (org_sigma1 == org_sigma2) {  // longitudinal magnetic field
          mltply::Spin::GC::Half::X_CisAis_MPIdouble(org_isite1, org_sigma1, 1.0, nstate, Xvec, vec);
        }
        else {  // transverse magnetic field
          mltply::Spin::GC::Half::X_CisAit_MPIdouble(org_isite1, org_sigma1, org_sigma2, 1.0, nstate, Xvec, vec);
        }
      }
      else {
        isite1 = Def::Tpow[org_isite1];

        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn,dmv) \
shared(i_max, isite1, org_sigma1,vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            dmv = mltply::Spin::GC::Half::X_CisAis(j, isite1, org_sigma1);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
        else {
          // transverse magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off,dmv) \
shared(i_max, isite1, org_sigma2, vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = mltply::Spin::GC::Half::X_CisAit(j, isite1, org_sigma2, &tmp_off);
            if (tmp_sgn != 0) {
              dmv = (std::complex<double>)tmp_sgn;
              zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
}
/**
 * @brief function of calculation for one body green's function for General SpinGC model.
 */
void expec::cisajs::Spin::GC::General(
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

  i_max = Check::idim_max;

  for (i = 0; i < Def::NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = Def::CisAjt[i][0];
    org_isite2 = Def::CisAjt[i][2];
    org_sigma1 = Def::CisAjt[i][1];
    org_sigma2 = Def::CisAjt[i][3];
    if (org_isite1 == org_isite2) {
      if (org_isite1 >= Def::Nsite) {
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          mltply::Spin::GC::General::X_CisAis_MPIdouble(org_isite1, org_sigma1,
            1.0, nstate, Xvec, vec);
        }
        else {
          // transverse magnetic field
          mltply::Spin::GC::General::X_CisAit_MPIdouble(
            org_isite1, org_sigma1, org_sigma2, 1.0, nstate, Xvec, vec);
        }
      }
      else {//org_isite1 <= Def::Nsite
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1,dmv) \
shared(i_max, org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow, vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, Def::SiteToBit, Def::Tpow);
            dmv = (std::complex<double>)num1;
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
        else {
          // transverse magnetic field
#pragma omp parallel for default(none) private(j, num1,dmv) \
shared(i_max, org_isite1, org_sigma1, org_sigma2, tmp_off,Def::SiteToBit, Def::Tpow,vec,Xvec,nstate,one)
          for (j = 1; j <= i_max; j++) {
            num1 = GetOffCompGeneralSpin(
              j - 1, org_isite1, org_sigma2, org_sigma1, &tmp_off, Def::SiteToBit, Def::Tpow);
            if (num1 != 0) {
              dmv = (std::complex<double>)num1;
              zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
    }
    wrapperMPI::MultiVecProd(i_max, nstate, vec, Xvec, prod[i]);
  }
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
int expec::cisajs::main(
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

  if (Def::NCisAjt < 1) return 0;

  i_max = Check::idim_max;
  if (GetSplitBitByModel(Def::Nsite, Def::iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }
  Large::i_max = i_max;
  Large::irght = irght;
  Large::ilft = ilft;
  Large::ihfbit = ihfbit;
  Large::mode = M_CORR;

  switch (Def::iCalcType) {
  case DC::TPQCalc:
    step = Def::istep;
    TimeKeeperWithRandAndStep("%s_TimeKeeper.dat", "set %d step %d:expec_cisajs begins: %s", "a", 0, step);
    break;
  case DC::TimeEvolution:
    step = Def::istep;
    TimeKeeperWithStep("%s_TimeKeeper.dat", "step %d:expec_cisajs begins: %s", "a", step);
    break;
  case DC::FullDiag:
  case DC::CG:
    break;
  }

  prod = cd_2d_allocate(Def::NCisAjt, nstate);
  switch (Def::iCalcModel) {
  case DC::HubbardGC:
    expec::cisajs::Hubbard::GC(nstate, Xvec, vec, prod);
    break;

  case DC::KondoGC:
  case DC::Hubbard:
  case DC::Kondo:
    expec::cisajs::Hubbard::C(nstate, Xvec, vec, prod);
    break;

  case DC::Spin: // for the Sz-conserved spin system
    if (Def::iFlgGeneralSpin == FALSE) {
      expec::cisajs::Spin::C::Half(nstate, Xvec, vec, prod);
    }
    else {
      expec::cisajs::Spin::C::General(nstate, Xvec, vec, prod);
    }
    break;

  case DC::SpinGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      expec::cisajs::Spin::GC::Half(nstate, Xvec, vec, prod);
    }
    else {
      expec::cisajs::Spin::GC::General(nstate, Xvec, vec, prod);
    }
    break;

  default:
    break;
  }

  for (istate = 0; istate < nstate; istate++) {
    switch (Def::iCalcType) {
    case DC::TPQCalc:
      step = Def::istep;
      sprintf(sdt, "%s_cisajs_set%dstep%d.dat", Def::CDataFileHead, istate, step);
      break;
    case DC::TimeEvolution:
      step = Def::istep;
      sprintf(sdt, "%s_cisajs_step%d.dat", Def::CDataFileHead, step);
      break;
    case DC::FullDiag:
    case DC::CG:
      sprintf(sdt, "%s_cisajs_eigen%d.dat", Def::CDataFileHead, istate);
      break;
    }
    if (childfopenMPI(sdt, "w", &fp) == 0) {
      for (ica = 0; ica < Def::NCisAjt; ica++) {
        fprintf(fp, " %4d %4d %4d %4d %.10lf %.10lf\n",
          Def::CisAjt[ica][0], Def::CisAjt[ica][1], Def::CisAjt[ica][2], Def::CisAjt[ica][3],
          real(prod[ica][istate]), imag(prod[ica][istate]));
      }
      fclose(fp);
    }
    else return -1;
  }/*for (istate = 0; istate < nstate; istate++)*/

  if (Def::St == 0) {
    if (Def::iCalcType == DC::TPQCalc) {
      TimeKeeperWithRandAndStep("%s_TimeKeeper.dat", "set %d step %d:expec_cisajs finishes: %s", "a", rand_i, step);
    }
    else if (Def::iCalcType == DC::TimeEvolution) {
      TimeKeeperWithStep("%s_TimeKeeper.dat", "step %d:expec_cisajs finishes: %s", "a", step);
    }
  }
  else if (Def::St == 1) {
    TimeKeeper("%s_TimeKeeper.dat", "CG expec_cisajs finishes:     %s", "a");
    fprintf(MP::STDOUT, "%s", "  End  : Calculate one body Green functions.\n\n");
  }
  free_cd_2d_allocate(prod);
  return 0;
}
