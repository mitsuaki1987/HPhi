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

std::complex<double> child_Cor_1(long int j,
  long int is_A, long int is1, long int is2, long int is3,
  long int b_sigma_A, long int irght, long int ilft, long int ihfbit, 
  std::complex<double>* vec);

std::complex<double> child_Cor_2(long int j,
  long int is_B, long int is1,
  long int is3, long int is4,
  long int b_sigma_B, long int irght, long int ilft,
  long int ihfbit, std::complex<double>* vec);

std::complex<double> child_Cor_3(long int j,
  long int is_A, long int is_B, long int is1,
  long int is2, long int is3, long int is4,
  long int b_sigma_A, long int b_sigma_B,
  long int irght, long int ilft, long int ihfbit, std::complex<double>* vec);
int expec_cisajscktaltdc(int nstate, std::complex<double>** Xvec, std::complex<double>** vec);
void expec_cisajscktaltdc_alldiag_spin(std::complex<double>* vec);
