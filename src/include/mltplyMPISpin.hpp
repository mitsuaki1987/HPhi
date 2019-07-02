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

//Define Mode for mltply
// complex version

#pragma once
#include <complex>

void child_general_int_spin_MPIdouble
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_child_general_int_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 std::complex<double> tmp_J,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );


void X_child_general_int_spin_TotalS_MPIdouble
(
 int org_isite1,
 int org_isite3,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void child_general_int_spin_MPIsingle
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_child_general_int_spin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 std::complex<double> tmp_J,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void GC_child_general_int_spin_MPIdouble
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void GC_child_general_int_spin_MPIsingle
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void GC_child_general_int_GeneralSpin_MPIdouble
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void GC_child_general_int_GeneralSpin_MPIsingle
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void child_general_int_GeneralSpin_MPIdouble
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void child_general_int_GeneralSpin_MPIsingle
(
 long int i_int,
 
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );
