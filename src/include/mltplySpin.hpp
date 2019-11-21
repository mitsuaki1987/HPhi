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

#ifndef HPHI_MLTPLYSPIN_H
#define HPHI_MLTPLYSPIN_H

namespace mltply {
  namespace Spin{
    namespace C{
      int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
        long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      namespace Half {
        int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
          long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
        void general_int(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
          long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
        void exchange(
          int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
          long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      }
      namespace General {
        int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
          long int i_max, long int* list_1, long int* list_2_1, long int* list_2_2);
      }
    }
    namespace GC {
      int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
      namespace Half {
        int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
        void general_int(
          int nstate, std::complex<double>** tmp_v0,
          std::complex<double>** tmp_v1);

        void exchange(
          int nstate, std::complex<double>** tmp_v0,
          std::complex<double>** tmp_v1);

        void pairlift(int nstate, std::complex<double>** tmp_v0,
          std::complex<double>** tmp_v1);
      }
      namespace General {
        int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
      }
    }
    namespace GCBoost {
      int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    }

  }
}
#endif
