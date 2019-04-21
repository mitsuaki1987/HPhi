/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//
// Created by Kazuyoshi Yoshimi on 2019-01-09.
//

#ifndef MVMC_SETMEMORY_H
#define MVMC_SETMEMORY_H

#include <cstdlib>
#include <cstring>
#include <complex>


int *i_1d_allocate(const long int N);
void free_i_1d_allocate(int *A);
long int *li_1d_allocate(const long int N);
void free_li_1d_allocate(long int *A);
long int **li_2d_allocate(const long int N, const long int M);
void free_li_2d_allocate(long int **A);
void free_li_1d_allocate(long int *A);
int **i_2d_allocate(const long int N, const long int M);
void free_i_2d_allocate(int **A);
int ***i_3d_allocate(const long int N, const long int M, const long int L);
void free_i_3d_allocate(int ***A);
double *d_1d_allocate(const long int N);
void free_d_1d_allocate(double *A);
double **d_2d_allocate(const long int N, const long int M);
void free_d_2d_allocate(double **A);
std::complex<double> *cd_1d_allocate(const long int N);
void free_cd_1d_allocate(std::complex<double>*A);
std::complex<double> **cd_2d_allocate(const long int N, const long int M);
void free_cd_2d_allocate(std::complex<double>**A);
std::complex<double>***cd_3d_allocate(const long int N, const long int M, const long int L);
void free_cd_3d_allocate(std::complex<double>***A);
std::complex<double>****cd_4d_allocate(const long int N, const long int M, const long int L, const long int K);
void free_cd_4d_allocate(std::complex<double>****A);

#endif //MVMC_SETMEMORY_H
