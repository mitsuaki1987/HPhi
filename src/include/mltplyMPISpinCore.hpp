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
#include <complex>

namespace mltply {
  namespace Spin {
    namespace C {
      namespace Half {
        void X_CisAit_MPIdouble(int org_isite1, int org_ispin2,
          std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
          long int idim_max, long int* list_1, long int *list_2_1, long int *list_2_2);
      }
      namespace General {
        void X_CisAit_MPIdouble(int org_isite1, int org_ispin1, int org_ispin2,
          std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
          long int idim_max, long int *list_1, long int *list_2_1, long int *list_2_2);
        void X_CisAisCjuAju_MPIdouble(
          int org_isite1, int org_ispin1, int org_isite3, int org_ispin3,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAjv_MPIdouble(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
          long int i_max, long int *list_1, long int *list_2_1, long int *list_2_2);
        void X_CisAisCjuAju_MPIsingle(
          int org_isite1, int org_ispin1, int org_isite3, int org_ispin3,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1, long int *list_1);
        void X_CisAitCjuAjv_MPIsingle(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1, 
          long int i_max, long int *list_1, long int *list_2_1, long int *list_2_2);
      }
    }
    namespace GC {
      namespace Half {
        void X_CisAitCiuAiv_MPIdouble(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAjv_MPIdouble(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAju_MPIdouble(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAju_MPIdouble(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCiuAiv_MPIsingle(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAjv_MPIsingle(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAju_MPIsingle(int org_isite1, int org_ispin2,
          int org_isite3, int org_ispin3, std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAju_MPIsingle(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAju_MPIsingle(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAit_MPIdouble(int org_isite1, int org_ispin1, int org_ispin2,
          std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAis_MPIdouble(int org_isite1, int org_ispin1,
          std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_AisCis_MPIdouble(int org_isite1, int org_ispin1,
          std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void CisAisCjuAjv_MPIdouble(long int i_int, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void CisAitCjuAju_MPIdouble(long int i_int, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void CisAitCiuAiv_MPIdouble(long int i_int, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAjv_MPIsingle(long int i_int, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAju_MPIsingle(long int i_int, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCiuAiv_MPIsingle(long int i_int, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
      }
      namespace General {
        void X_CisAisCjuAjv_MPIdouble(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAju_MPIdouble(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAjv_MPIdouble(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAjv_MPIsingle(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAju_MPIsingle(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAitCjuAjv_MPIsingle(int org_isite1, int org_ispin1,
          int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAit_MPIdouble(int org_isite1, int org_ispin1, int org_ispin2,
          std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAis_MPIdouble(int org_isite1, int org_ispin1,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_AisCis_MPIdouble(int org_isite1, int org_ispin1,
          std::complex<double> tmp_J, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAju_MPIdouble(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void X_CisAisCjuAju_MPIsingle(int org_isite1, int org_ispin1,
          int org_isite3, int org_ispin3, std::complex<double> tmp_trans, int nstate,
          std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
      }
    }
  }
}
