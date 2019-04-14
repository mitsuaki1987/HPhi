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
/*-------------------------------------------------------------*/
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <cstring>

int OutputRealEigenValue(int xNsize, double *ene, char *filename);
int OutputCmpEigenValue(int xNsize, std::complex<double> *ene, char *filename);
int OutputRealEigenVec(int xNsize, const int nene, double **vec, const int nproc, char *filename);
int OutputCmpEigenVec(int xNsize, const int nene, std::complex<double> **vec, const int nproc, char *filename);
int InputRealEigenValue(int xNsize, double *ene, char *filename);
int InputCmpEigenValue(int xNsize, std::complex<double> *ene, char *filename);
int InputRealEigenVec(int xNsize, const int nene, double **vec, const int nproc, char *filename);
int InputCmpEigenVec(int xNsize, const int nene, std::complex<double> **vec, const int nproc, char *filename);
