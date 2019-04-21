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

#pragma once
#include <complex>
#include "struct.hpp"

int CheckPE
(
 int isite,
 struct BindStruct *X
 );

int CheckBit_Cis
(
 long int is1_spin,
 long int orgbit,
 long int *offbit
 );

int CheckBit_Ajt
(
 long int is1_spin,
 long int orgbit,
 long int *offbit
 );

int CheckBit_InterAllPE
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 struct BindStruct *X,
 long int orgbit,
 long int *offbit
 );

int CheckBit_PairPE
(
 int isite1,
 int isigma1,
 int isite3,
 int isigma3,
 struct BindStruct *X,
 long int orgbit
 );

int GetSgnInterAll
(
 long int isite1,
 long int isite2,
 long int isite3,
 long int isite4,
 int *Fsgn,
 struct BindStruct *X,
 long int orgbit,
 long int *offbit
 );

void X_GC_child_CisAisCjtAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_GC_child_CisAjtCkuAlv_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_GC_child_CisAjtCkuAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_GC_child_CisAisCjtAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_GC_child_CisAis_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
);

void X_GC_child_CisAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
);

void X_child_CisAisCjtAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_child_CisAjtCkuAlv_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_child_CisAjtCkuAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_child_CisAisCjtAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );

void X_child_CisAis_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 std::complex<double> tmp_V,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
);

void X_child_CisAjt_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 std::complex<double> tmp_trans,
 struct BindStruct *X,
 int nstate, 
  std::complex<double> **tmp_v0, 
 std::complex<double> **tmp_v1
 );

void X_child_CisAjt_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 std::complex<double> tmp_trans,
 struct BindStruct *X,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1
 );


void X_GC_Cis_MPI
(
 int org_isite,
 int org_ispin,
 std::complex<double> tmp_trans,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 long int idim_max,
 long int *Tpow
 );

void X_GC_Ajt_MPI
(
 int org_isite,
 int org_ispin,
 std::complex<double> tmp_trans,
 int nstate, 
  std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 long int idim_max,
 long int *Tpow
 );

void X_Cis_MPI
(
 int org_isite,
 int org_ispin,
 std::complex<double> tmp_trans,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 long int idim_max,
 long int *Tpow,
 long int _irght,
 long int _ilft,
 long int _ihfbit
 );

void X_Ajt_MPI
(
 int org_isite,
 int org_ispin,
 std::complex<double> tmp_trans,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 long int idim_max,
 long int *Tpow,
 long int _irght,
 long int _ilft,
 long int _ihfbit
 );
