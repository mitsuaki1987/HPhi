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
#include <bitcalc.hpp>
#include "mltplyCommon.hpp"
#include "mltplySpinCore.hpp"
#include "wrapperMPI.hpp"
#include "mltplyMPISpin.hpp"
#include "mltplyMPISpinCore.hpp"
#include "expec_totalspin.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
/**
 * @file
 *
 * @brief  File for calculating total spin
 *
 * @version 0.2
 * @details modify to treat the case of general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 */
/**
 * @brief function of calculating totalspin for Hubbard model
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_Hubbard(
  
  int nstate,
  std::complex<double> **vec
) {
  long int j;
  long int irght, ilft, ihfbit;
  long int isite1, isite2;
  long int is1_up, is2_up, is1_down, is2_down;
  long int iexchg, off;
  int num1_up, num2_up, istate;
  int num1_down, num2_down;
  long int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  std::complex<double> tmp_spn_z;
  long i_max;
  i_max = Check::idim_max;

  GetSplitBitByModel(Def::Nsite, Def::iCalcModel, &irght, &ilft, &ihfbit);
  for (istate = 0; istate < nstate; istate++) {
    Phys::s2[istate] = 0.0;
    Phys::Sz[istate] = 0.0;
  }
  for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
    is1_up = Def::Tpow[2 * isite1 - 2];
    is1_down = Def::Tpow[2 * isite1 - 1];

    for (isite2 = 1; isite2 <= Def::NsiteMPI; isite2++) {
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

      for (j = 1; j <= i_max; j++) {

        ibit1_up = List::c1[j] & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit2_up = List::c1[j] & is2_up;
        num2_up = ibit2_up / is2_up;

        ibit2_down = List::c1[j] & is2_down;
        num2_down = ibit2_down / is2_down;
        ibit1_down = List::c1[j] & is1_down;
        num1_down = ibit1_down / is1_down;

        tmp_spn_z = (num1_up - num1_down) * (num2_up - num2_down);
        for (istate = 0; istate < nstate; istate++)
          Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate] * tmp_spn_z / 4.0);
        if (isite1 == isite2) {
          for (istate = 0; istate < nstate; istate++) {
            Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate]) * (num1_up + num1_down - 2 * num1_up * num1_down) / 2.0;
            Phys::Sz[istate] += real(conj(vec[j][istate]) * vec[j][istate]) * (num1_up - num1_down) / 2.0;
          }
        }
        else {
          if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
            iexchg = List::c1[j] - (is1_up + is2_down);
            iexchg += (is2_up + is1_down);
            GetOffComp(List::c2_1, List::c2_2, iexchg, irght, ilft, ihfbit, &off);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[off][istate]) / 2.0;
          }
          else if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {
            iexchg = List::c1[j] - (is1_down + is2_up);
            iexchg += (is2_down + is1_up);
            GetOffComp(List::c2_1, List::c2_2, iexchg, irght, ilft, ihfbit, &off);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[off][istate]) / 2.0;
          }
        }
      }
    }
  }
  SumMPI_dv(nstate, Phys::s2);
  SumMPI_dv(nstate, Phys::Sz);
}
/**
 * @brief function of calculating totalspin for Hubbard model in grand canonical ensemble
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_HubbardGC(
  
  int nstate,
  std::complex<double> **vec
) {
  long int j;
  long int isite1, isite2;
  long int is1_up, is2_up, is1_down, is2_down;
  long int iexchg, off;
  int num1_up, num2_up, istate;
  int num1_down, num2_down;
  long int ibit1_up, ibit2_up, ibit1_down, ibit2_down, list_1_j;
  std::complex<double> tmp_spn_z;
  long int i_max;

  i_max = Check::idim_max;

  for (istate = 0; istate < nstate; istate++) {
    Phys::s2[istate] = 0.0;
    Phys::Sz[istate] = 0.0;
  }
  for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
    for (isite2 = 1; isite2 <= Def::NsiteMPI; isite2++) {
      is1_up = Def::Tpow[2 * isite1 - 2];
      is1_down = Def::Tpow[2 * isite1 - 1];
      is2_up = Def::Tpow[2 * isite2 - 2];
      is2_down = Def::Tpow[2 * isite2 - 1];

      for (j = 1; j <= i_max; j++) {
        list_1_j = j - 1;
        ibit1_up = list_1_j & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit2_up = list_1_j & is2_up;
        num2_up = ibit2_up / is2_up;

        ibit1_down = list_1_j & is1_down;
        num1_down = ibit1_down / is1_down;
        ibit2_down = list_1_j & is2_down;
        num2_down = ibit2_down / is2_down;

        tmp_spn_z = (num1_up - num1_down) * (num2_up - num2_down);
        for (istate = 0; istate < nstate; istate++)
          Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate] * tmp_spn_z / 4.0);
        if (isite1 == isite2) {
          Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate]) * (num1_up + num1_down - 2 * num1_up * num1_down) / 2.0;
          Phys::Sz[istate] += real(conj(vec[j][istate]) * vec[j][istate]) * (num1_up - num1_down) / 2.0;
        }
        else {
          if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
            iexchg = list_1_j - (is1_up + is2_down);
            iexchg += (is2_up + is1_down);
            off = iexchg + 1;
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[off][istate] / 2.0);
          }
          else if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {
            iexchg = list_1_j - (is1_down + is2_up);
            iexchg += (is2_down + is1_up);
            off = iexchg + 1;
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[off][istate] / 2.0);
          }
        }
      }
    }
  }
  SumMPI_dv(nstate, Phys::s2);
  SumMPI_dv(nstate, Phys::Sz);
}
/**
 * @brief function of calculating totalspin for spin model
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.2
 * @details modify for hybrid parallel
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_Spin(
   
  int nstate,
  std::complex<double> **vec
) {
  long int j;
  long int irght, ilft, ihfbit;
  long int isite1, isite2;
  long int tmp_isite1, tmp_isite2;

  long int is1_up, is2_up;
  long int iexchg, off, off_2;

  int num1_up, num2_up;
  int num1_down, num2_down;
  int sigma_1, sigma_2, istate;
  long int ibit1_up, ibit2_up, ibit_tmp, is_up;
  std::complex<double> spn_z1, spn_z2;
  long int i_max;
  double spn_z;

  i_max = Check::idim_max;
  for (istate = 0; istate < nstate; istate++) {
    Phys::s2[istate] = 0.0;
    Phys::Sz[istate] = 0.0;
  }
  if (Def::iFlgGeneralSpin == FALSE) {
    GetSplitBitByModel(Def::Nsite, Def::iCalcModel, &irght, &ilft, &ihfbit);
    for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
      for (isite2 = 1; isite2 <= Def::NsiteMPI; isite2++) {

        if (isite1 > Def::Nsite && isite2 > Def::Nsite) {
#ifdef __MPI
          is1_up = Def::Tpow[isite1 - 1];
          is2_up = Def::Tpow[isite2 - 1];
          is_up = is1_up + is2_up;
          num1_up = X_SpinGC_CisAis((long int) MP::myrank + 1, is1_up, 1);
          num1_down = 1 - num1_up;
          num2_up = X_SpinGC_CisAis((long int) MP::myrank + 1, is2_up, 1);
          num2_down = 1 - num2_up;
          spn_z = (num1_up - num1_down) * (num2_up - num2_down);

          for (j = 1; j <= i_max; j++) {
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate] * spn_z) / 4.0;
          }
          if (isite1 == isite2) {
            for (j = 1; j <= i_max; j++) {
              for (istate = 0; istate < nstate; istate++)
                Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate]) / 2.0;
            }
          }
          else {//off diagonal
            //debug spn += X_child_general_int_spin_TotalS_MPIdouble(isite1 - 1, isite2 - 1, nstate, vec, vec);
          }
#endif
        }
        else if (isite1 > Def::Nsite || isite2 > Def::Nsite) {
#ifdef __MPI
          if (isite1 < isite2) {
            tmp_isite1 = isite1;
            tmp_isite2 = isite2;
          }
          else {
            tmp_isite1 = isite2;
            tmp_isite2 = isite1;
          }

          is1_up = Def::Tpow[tmp_isite1 - 1];
          is2_up = Def::Tpow[tmp_isite2 - 1];
          num2_up = X_SpinGC_CisAis((long int) MP::myrank + 1, is2_up, 1);
          num2_down = 1 - num2_up;

          //diagonal
          for (j = 1; j <= i_max; j++) {
            ibit1_up = List::c1[j] & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            spn_z = (num1_up - num1_down) * (num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate] * spn_z) / 4.0;
          }
          if (isite1 < isite2) {
            //debug spn += X_child_general_int_spin_MPIsingle(isite1 - 1, 0, 1, isite2 - 1, 1, 0, 1.0, nstate, vec, vec);
          }
          else {
            //debug spn += conj(X_child_general_int_spin_MPIsingle(isite2 - 1, 1, 0, isite1 - 1, 0, 1, 1.0, nstate, vec, vec));
          }
#endif
        }//isite1 > Nsite || isite2 > Nsite
        else {
          is1_up = Def::Tpow[isite1 - 1];
          is2_up = Def::Tpow[isite2 - 1];
          is_up = is1_up + is2_up;

          for (j = 1; j <= i_max; j++) {
            ibit1_up = List::c1[j] & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            ibit2_up = List::c1[j] & is2_up;
            num2_up = ibit2_up / is2_up;
            num2_down = 1 - num2_up;

            spn_z = (num1_up - num1_down) * (num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate] * spn_z) / 4.0;

            if (isite1 == isite2) {
              for (istate = 0; istate < nstate; istate++)
                Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate]) / 2.0;
            }
            else {
              ibit_tmp = (num1_up) ^ (num2_up);
              if (ibit_tmp != 0) {
                iexchg = List::c1[j] ^ (is_up);
                GetOffComp(List::c2_1, List::c2_2, iexchg, irght, ilft, ihfbit, &off);
                for (istate = 0; istate < nstate; istate++)
                  Phys::s2[istate] += real(conj(vec[j][istate]) * vec[off][istate]) / 2.0;
              }
            }
          }// j
        }
      }//isite2
    }//isite1
  }//generalspin=FALSE
  else {
    double S1 = 0;
    double S2 = 0;
    for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
      for (isite2 = 1; isite2 <= Def::NsiteMPI; isite2++) {
        S1 = 0.5 * (Def::SiteToBit[isite1 - 1] - 1);
        S2 = 0.5 * (Def::SiteToBit[isite2 - 1] - 1);
        if (isite1 == isite2) {
          for (j = 1; j <= i_max; j++) {
            spn_z1 = 0.5 * GetLocal2Sz(isite1, List::c1[j], Def::SiteToBit, Def::Tpow);
            for (istate = 0; istate < nstate; istate++) {
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate]) * S1 * (S1 + 1.0);
              Phys::Sz[istate] += real(conj(vec[j][istate]) * vec[j][istate] * spn_z1);
            }
          }
        }
        else {
          for (j = 1; j <= i_max; j++) {
            spn_z1 = 0.5 * GetLocal2Sz(isite1, List::c1[j], Def::SiteToBit, Def::Tpow);
            spn_z2 = 0.5 * GetLocal2Sz(isite2, List::c1[j], Def::SiteToBit, Def::Tpow);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate]) * vec[j][istate] * spn_z1 * spn_z2);

            sigma_1 = GetBitGeneral(isite1, List::c1[j], Def::SiteToBit, Def::Tpow);
            sigma_2 = GetBitGeneral(isite2, List::c1[j], Def::SiteToBit, Def::Tpow);

            ibit_tmp = GetOffCompGeneralSpin(List::c1[j], isite2, sigma_2, sigma_2 + 1, &off, Def::SiteToBit,
              Def::Tpow);
            if (ibit_tmp == TRUE) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 - 1, &off_2, Def::SiteToBit,
                Def::Tpow);
              if (ibit_tmp == TRUE) {
                ConvertToList1GeneralSpin(off_2, Check::sdim, &off);
                for (istate = 0; istate < nstate; istate++)
                  Phys::s2[istate] += real(conj(vec[j][istate]) * vec[off][istate])
                  * sqrt(S2 * (S2 + 1) - real(spn_z2) * (real(spn_z2) + 1)) 
                  * sqrt(S1 * (S1 + 1) - real(spn_z1) * (real(spn_z1) - 1)) / 2.0;
              }
            }

            ibit_tmp = GetOffCompGeneralSpin(List::c1[j], isite2, sigma_2, sigma_2 - 1, &off, Def::SiteToBit,
              Def::Tpow);
            if (ibit_tmp == TRUE) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 + 1, &off_2, Def::SiteToBit,
                Def::Tpow);
              if (ibit_tmp == TRUE) {
                ConvertToList1GeneralSpin(off_2, Check::sdim, &off);
                for (istate = 0; istate < nstate; istate++)
                  Phys::s2[istate] += real(conj(vec[j][istate]) * vec[off][istate])
                  * sqrt(S2 * (S2 + 1) - real(spn_z2) * (real(spn_z2) - 1.0)) 
                  * sqrt(S1 * (S1 + 1) - real(spn_z1) * (real(spn_z1) + 1)) / 2.0;
              }
            }
          }
        }
      }
    }
  }
  SumMPI_dv(nstate, Phys::s2);
  SumMPI_dv(nstate, Phys::Sz);
}
/**
 * @brief function of calculating totalspin for spin model in grand canonical ensemble
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.2
 * @details add function to treat a calculation of total spin for general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_SpinGC(
  
  int nstate,
  std::complex<double> **vec
) {
  long int j;
  long int isite1, isite2, tmp_isite1, tmp_isite2;
  long int is1_up, is2_up;
  long int iexchg, off, off_2;
  int num1_up, num2_up, istate;
  int num1_down, num2_down;
  int sigma_1, sigma_2;
  long int ibit1_up, ibit2_up, ibit_tmp, is_up;
  std::complex<double> spn_z1, spn_z2;
  long int list_1_j;
  long int i_max;

  i_max = Check::idim_max;
  Large::mode = M_TOTALS;
  for (istate = 0; istate < nstate; istate++) {
    Phys::s2[istate] = 0.0;
    Phys::Sz[istate] = 0.0;
  }
  if (Def::iFlgGeneralSpin == FALSE) {
    for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
      if (isite1 > Def::Nsite) {
        is1_up = Def::Tpow[isite1 - 1];
        ibit1_up = MP::myrank & is1_up;
        num1_up = ibit1_up / is1_up;
        num1_down = 1 - num1_up;
        for (j = 1; j <= i_max; j++) {
          for (istate = 0; istate < nstate; istate++)
            Phys::Sz[istate] += real(conj(vec[j][istate])*vec[j][istate]) * (num1_up - num1_down) / 2.0;
        }
      }
      else {
        is1_up = Def::Tpow[isite1 - 1];
        for (j = 1; j <= i_max; j++) {
          list_1_j = j - 1;
          ibit1_up = list_1_j & is1_up;
          num1_up = ibit1_up / is1_up;
          num1_down = 1 - num1_up;
          for (istate = 0; istate < nstate; istate++)
            Phys::Sz[istate] += real(conj(vec[j][istate])*vec[j][istate]) * (num1_up - num1_down) / 2.0;
        }
      }
      for (isite2 = 1; isite2 <= Def::NsiteMPI; isite2++) {

        if (isite1 > Def::Nsite && isite2 > Def::Nsite) {
          is1_up = Def::Tpow[isite1 - 1];
          is2_up = Def::Tpow[isite2 - 1];
          num1_up = X_SpinGC_CisAis((long int)MP::myrank + 1, is1_up, 1);
          num1_down = 1 - num1_up;
          num2_up = X_SpinGC_CisAis((long int)MP::myrank + 1, is2_up, 1);
          num2_down = 1 - num2_up;
          spn_z2 = (num1_up - num1_down)*(num2_up - num2_down) / 4.0;
          for (j = 1; j <= i_max; j++) {
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate] * spn_z2);
          }
          if (isite1 == isite2) {
            for (j = 1; j <= i_max; j++) {
              for (istate = 0; istate < nstate; istate++)
                Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate]) / 2.0;
            }
          }//isite1 = isite2
          else {//off diagonal
            //debug spn += X_GC_child_CisAitCiuAiv_spin_MPIdouble(
            //debug   isite1 - 1, 0, 1, isite2 - 1, 1, 0, 1.0, nstate, vec, vec) / 2.0;
          }
        }
        else if (isite1 > Def::Nsite || isite2 > Def::Nsite) {
          if (isite1 < isite2) {
            tmp_isite1 = isite1;
            tmp_isite2 = isite2;
          }
          else {
            tmp_isite1 = isite2;
            tmp_isite2 = isite1;
          }
          is1_up = Def::Tpow[tmp_isite1 - 1];
          is2_up = Def::Tpow[tmp_isite2 - 1];
          num2_up = X_SpinGC_CisAis((long int)MP::myrank + 1, is2_up, 1);
          num2_down = 1 - num2_up;
          //diagonal
          for (j = 1; j <= i_max; j++) {
            list_1_j = j - 1;
            ibit1_up = list_1_j & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            spn_z2 = (num1_up - num1_down)*(num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate] * spn_z2) / 4.0;
          }
          if (isite1 < isite2) {
            //debug spn += X_GC_child_CisAitCiuAiv_spin_MPIsingle(isite1 - 1, 0, 1, isite2 - 1, 1, 0, 1.0, nstate, vec, vec) / 2.0;
          }
          else {
            //debug spn += conj(X_GC_child_CisAitCiuAiv_spin_MPIsingle(isite2 - 1, 1, 0, isite1 - 1, 0, 1, 1.0, nstate, vec, vec)) / 2.0;
          }
        }
        else {
          is2_up = Def::Tpow[isite2 - 1];
          is_up = is1_up + is2_up;
          for (j = 1; j <= i_max; j++) {
            list_1_j = j - 1;
            ibit1_up = list_1_j & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            ibit2_up = list_1_j & is2_up;
            num2_up = ibit2_up / is2_up;
            num2_down = 1 - num2_up;

            spn_z2 = (num1_up - num1_down)*(num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate] * spn_z2) / 4.0;

            if (isite1 == isite2) {
              for (istate = 0; istate < nstate; istate++)
                Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate]) / 2.0;
            }
            else {
              ibit_tmp = (num1_up) ^ (num2_up);
              if (ibit_tmp != 0) {
                iexchg = list_1_j ^ (is_up);
                off = iexchg + 1;
                for (istate = 0; istate < nstate; istate++)
                  Phys::s2[istate] += real(conj(vec[j][istate])*vec[off][istate]) / 2.0;
              }
            }
          }//j  
        }//else
      }
    }
  }
  else {//general spin
    double S1 = 0;
    double S2 = 0;
    for (isite1 = 1; isite1 <= Def::NsiteMPI; isite1++) {
      S1 = 0.5*(Def::SiteToBit[isite1 - 1] - 1);
      if (isite1 > Def::Nsite) {
        spn_z1 = 0.5*GetLocal2Sz(isite1, (long int) MP::myrank, Def::SiteToBit, Def::Tpow);
        for (j = 1; j <= i_max; j++) {
          for (istate = 0; istate < nstate; istate++) {
            Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate]) * S1*(S1 + 1.0);
            Phys::Sz[istate] += real(conj(vec[j][istate])*vec[j][istate] * spn_z1);
          }
        }
      }
      else {
        for (j = 1; j <= i_max; j++) {
          spn_z1 = 0.5*GetLocal2Sz(isite1, j - 1, Def::SiteToBit, Def::Tpow);
          for (istate = 0; istate < nstate; istate++) {
            Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate]) * S1*(S1 + 1.0);
            Phys::Sz[istate] += real(conj(vec[j][istate])*vec[j][istate] * spn_z1);
          }
        }
      }
      for (isite2 = 1; isite2 <= Def::NsiteMPI; isite2++) {
        if (isite1 == isite2) continue;
        S2 = 0.5*(Def::SiteToBit[isite2 - 1] - 1);
        if (isite1 > Def::Nsite && isite2 > Def::Nsite) {
        }
        else if (isite1 > Def::Nsite || isite2 > Def::Nsite) {
        }
        else { //inner-process
          for (j = 1; j <= i_max; j++) {
            spn_z1 = 0.5*GetLocal2Sz(isite1, j - 1, Def::SiteToBit, Def::Tpow);
            spn_z2 = 0.5*GetLocal2Sz(isite2, j - 1, Def::SiteToBit, Def::Tpow);
            for (istate = 0; istate < nstate; istate++)
              Phys::s2[istate] += real(conj(vec[j][istate])*vec[j][istate] * spn_z1*spn_z2);

            sigma_1 = GetBitGeneral(isite1, j - 1, Def::SiteToBit, Def::Tpow);
            sigma_2 = GetBitGeneral(isite2, j - 1, Def::SiteToBit, Def::Tpow);

            ibit_tmp = GetOffCompGeneralSpin(j - 1, isite2, sigma_2, sigma_2 + 1, &off, Def::SiteToBit, Def::Tpow);
            if (ibit_tmp != 0) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 - 1, &off_2, Def::SiteToBit, Def::Tpow);
              if (ibit_tmp != 0) {
                for (istate = 0; istate < nstate; istate++)
                  Phys::s2[istate] += real(conj(vec[j][istate])*vec[off_2 + 1][istate])
                  * sqrt(S2*(S2 + 1) - real(spn_z2) * (real(spn_z2) + 1))
                  * sqrt(S1*(S1 + 1) - real(spn_z1) * (real(spn_z1) - 1)) / 2.0;
              }
            }
            ibit_tmp = GetOffCompGeneralSpin(j - 1, isite2, sigma_2, sigma_2 - 1, &off, Def::SiteToBit, Def::Tpow);
            if (ibit_tmp != 0) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 + 1, &off_2, Def::SiteToBit, Def::Tpow);
              if (ibit_tmp != 0) {
                for (istate = 0; istate < nstate; istate++)
                  Phys::s2[istate] += real(conj(vec[j][istate])*vec[off_2 + 1][istate])
                  * sqrt(S2*(S2 + 1) - real(spn_z2) * (real(spn_z2) - 1.0))
                  * sqrt(S1*(S1 + 1) - real(spn_z1) * (real(spn_z1) + 1)) / 2.0;
              }
            }
          }//j  
        }//inner-process
      }//isite2
    }//isite1
  }
  SumMPI_dv(nstate, Phys::s2);
  SumMPI_dv(nstate, Phys::Sz);
}
/**
 * @brief Parent function of calculation of total spin
 *
 * @param[in,out] X data list of calculation parameters
 * @param[in] vec eigenvector
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @retval 0 calculation is normally finished
 */
int expec_totalspin
(
  
  int nstate,
  std::complex<double> **vec
)
{
  int istate;

  Large::mode = M_TOTALS;
  switch (Def::iCalcModel) {
  case Spin:
    totalspin_Spin(nstate, vec);
    for (istate = 0; istate < nstate; istate++)
      Phys::Sz[istate] = Def::Total2SzMPI / 2.;
    break;
  case SpinGC:
    totalspin_SpinGC(nstate, vec);
    break;
  case Hubbard:
  case Kondo:
    totalspin_Hubbard(nstate, vec);
    break;
  case HubbardGC:
  case KondoGC:
    totalspin_HubbardGC(nstate, vec);
    break;
  default:
    for (istate = 0; istate < nstate; istate++) {
      Phys::s2[istate] = 0.0;
      Phys::Sz[istate] = 0.0;
    }
  }
  return 0;
}
