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
  namespace Hubbard {
    namespace C {
      void general_hopp_MPIdouble(long int itrans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
        long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      void X_general_hopp_MPIdouble(int org_isite1, int org_ispin1,
        int org_isite2, int org_ispin2, std::complex<double> tmp_trans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
        long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      void general_hopp_MPIsingle(long int itrans,
        int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
        long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      void X_general_hopp_MPIsingle(int org_isite1, int org_ispin1,
        int org_isite2, int org_ispin2, std::complex<double> tmp_trans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
        long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      void X_CisAjt_MPIsingle(int org_isite1, int org_ispin1,
        int org_isite2, int org_ispin2, std::complex<double> tmp_trans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
        long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      void X_CisAjt_MPIdouble(int org_isite1, int org_ispin1,
        int org_isite2, int org_ispin2, std::complex<double> tmp_trans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
        long int i_max, long int *list_1, long int *list_2_1, long int *list_2_2);
    }
    namespace GC {
      void general_hopp_MPIdouble(long int itrans,
        int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
      void X_general_hopp_MPIdouble(
        int org_isite1, int org_ispin1, int org_isite2, int org_ispin2,
        std::complex<double> tmp_trans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
      void general_hopp_MPIsingle(long int itrans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
      void X_general_hopp_MPIsingle(
        int org_isite1, int org_ispin1, int org_isite2, int org_ispin2,
        std::complex<double> tmp_trans, int nstate,
        std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    }
  }
}
