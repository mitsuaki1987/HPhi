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
  namespace cisajscktalt {
    int main(int nstate, std::complex<double>** Xvec, std::complex<double>** vec);
    namespace Hubbard {
      void GC(int nstate, std::complex<double>** Xvec, 
        std::complex<double>** vec, std::complex<double>** prod);
      void C(int nstate, std::complex<double>** Xvec, 
        std::complex<double>** vec, std::complex<double>** prod);
    }
    namespace Spin {
      namespace C{
        void Half(int nstate, std::complex<double>** Xvec,
          std::complex<double>** vec, std::complex<double>** prod);
        void General(int nstate, std::complex<double>** Xvec,
          std::complex<double>** vec, std::complex<double>** prod);
      }
      namespace GC {
        void Half(int nstate, std::complex<double>** Xvec, 
          std::complex<double>** vec, std::complex<double>** prod);
        void General(int nstate, std::complex<double>** Xvec,
          std::complex<double>** vec, std::complex<double>** prod);
      }
    }
  }
}
//TODO void expec_cisajscktaltdc_alldiag_spin(std::complex<double>* vec);
