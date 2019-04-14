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

#ifndef HPHI_MLTPLYHUBBARD_H
#define HPHI_MLTPLYHUBBARD_H

#include "Common.hpp"

int mltplyHubbard(struct BindStruct *X, int nstate, std::complex<double> **tmp_v0,std::complex<double> **tmp_v1);

int mltplyHubbardGC(struct BindStruct *X, int nstate, std::complex<double> **tmp_v0,std::complex<double> **tmp_v1);

void GC_child_general_hopp
(
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 std::complex<double> trans
 );


void GC_child_general_int(
                         int nstate, std::complex<double> **tmp_v0,
                         std::complex<double> **tmp_v1,
                         struct BindStruct *X
                         );


void child_general_int
(
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X
 );

void child_general_hopp
(int nstate, std::complex<double> **tmp_v0, std::complex<double> **tmp_v1, struct BindStruct *X, std::complex<double> trans);

void child_exchange
(
 int nstate, 
  std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X
 );

void child_pairhopp
(
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X
 );

void GC_child_exchange
(
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X
 );

void GC_child_pairlift
(
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X
 );

void GC_child_pairhopp
(
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X
 );

#endif
