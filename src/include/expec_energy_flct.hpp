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

namespace expec {
  namespace energy_flct {
    namespace Hubbard{
      void GC(int nstate, std::complex<double>** tmp_v0);
      void C(int nstate, std::complex<double>** tmp_v0);
    }
    namespace Spin {
      namespace GC {
        void Half(int nstate, std::complex<double>** tmp_v0);
        void General(int nstate, std::complex<double>** tmp_v0);
      }
    }
    int main(int nstate, std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
  }
}
