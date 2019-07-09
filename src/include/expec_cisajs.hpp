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
  namespace cisajs {
    int main(int nstate, std::complex<double>** Xvec, std::complex<double>** vec);
    void HubbardGC(int nstate, std::complex<double>** Xvec, std::complex<double>** vec, std::complex<double>** prod);
    void Hubbard(int nstate, std::complex<double>** Xvec, std::complex<double>** vec, std::complex<double>** prod);
    void SpinHalf(int nstate, std::complex<double>** Xvec, std::complex<double>** vec, std::complex<double>** prod);
    void SpinGeneral(int nstate, std::complex<double>** Xvec, std::complex<double>** vec, std::complex<double>** prod);
    void SpinGCHalf(int nstate, std::complex<double>** Xvec, std::complex<double>** vec, std::complex<double>** prod);
    void SpinGCGeneral(int nstate, std::complex<double>** Xvec, std::complex<double>** vec, std::complex<double>** prod);
  }
}