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

#ifndef HPHI_MLTPLYSPINCORE_H
#define HPHI_MLTPLYSPINCORE_H

namespace mltply {
  namespace Spin {
    int general_int_GetInfo(
      long int isite1, long int isite2, long int sigma1, long int sigma2, long int sigma3,
      long int sigma4, std::complex<double> tmp_V);
    int exchange_GetInfo(int iExchange);
    int pairlift_GetInfo(int iPairLift);
    namespace C {
      namespace Half {
        void exchange_element(long int j, int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1, long int* tmp_off,
          long int* list_1, long int* list_2_1, long int* list_2_2);
        int X_exchange_element(long int j, long int isA_up, long int isB_up, long int sigmaA, long int sigmaB, long int* tmp_off,
          long int* list_1, long int* list_2_1, long int* list_2_2);
        void CisAisCisAis_element(long int j, long int isA_up, long int isB_up, long int org_sigma2, long int org_sigma4,
          std::complex<double> tmp_V, int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1, long int *list_1);
        int X_CisAit(long int j, long int is1_spin, long int sigma2, long int* tmp_off, long int *list_1, long int *list_2_1, long int *list_2_2);
        int X_CisAis(long int j, long int is1_spin, long int sigma1, long int *list_1);
      }
    }
    namespace GC {
      namespace Half {
        void pairlift_element(long int j, int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1, long int* tmp_off);
        void exchange_element(long int j, int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1, long int* tmp_off);
        void CisAisCisAis_element(
          long int j, long int isA_up, long int isB_up, long int org_sigma2, long int org_sigma4,
          std::complex<double> tmp_V, int nstate, std::complex<double>** tmp_v0,
          std::complex<double>** tmp_v1);
        void CisAisCitAiu_element(
          long int j, long int org_sigma2, long int org_sigma4, long int isA_up, long int isB_up,
          std::complex<double> tmp_V, int nstate, std::complex<double>** tmp_v0,
          std::complex<double>** tmp_v1, long int* tmp_off);
        void CisAitCiuAiu_element(
          long int j, long int org_sigma2, long int org_sigma4, long int isA_up, long int isB_up,
          std::complex<double> tmp_V, int nstate, std::complex<double>** tmp_v0,
          std::complex<double>** tmp_v1, long int* tmp_off);
        void CisAitCiuAiv_element(
          long int j, long int org_sigma2, long int org_sigma4, long int isA_up, long int isB_up,
          std::complex<double> tmp_V, int nstate, std::complex<double>** tmp_v0,
          std::complex<double>** tmp_v1, long int* tmp_off_2);
        int X_CisAit(long int j, long int is1_spin, long int sigma2, long int* tmp_off);
        int X_CisAis(long int j, long int is1_spin, long int sigma1);
      }
      namespace General {

      }
    }
  }
}

#endif
