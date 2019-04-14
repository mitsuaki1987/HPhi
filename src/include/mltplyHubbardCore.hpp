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

#ifndef HPHI_MLTPLYHUBBARDCORE_H
#define HPHI_MLTPLYHUBBARDCORE_H

#include "Common.hpp"

void child_pairhopp_element
(
 long unsigned int j,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

void GC_child_exchange_element
(
 long unsigned int j,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

void GC_child_pairhopp_element
(
 long unsigned int j,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

void child_exchange_element
(
 long unsigned int j,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

void child_CisAisCisAis_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite3,
 std::complex<double> tmp_V,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

void child_CisAisCjtAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                std::complex<double> tmp_V,
                int nstate, std::complex<double> **tmp_v0,
                std::complex<double> **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void child_CisAjtCkuAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int Asum,
                long unsigned int Adiff,
                std::complex<double> tmp_V,
                int nstate, std::complex<double> **tmp_v0,
                std::complex<double> **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void child_CisAjtCkuAlv_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Asum,
                long unsigned int Adiff,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                std::complex<double> tmp_V,
                int nstate, std::complex<double> **tmp_v0,
                std::complex<double> **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off_2
        );
//[s]Grand canonical
void GC_child_CisAisCisAis_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite3,
                std::complex<double> tmp_V,
                int nstate, std::complex<double> **tmp_v0,
                std::complex<double> **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void GC_child_CisAisCjtAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                std::complex<double> tmp_V,
                int nstate, std::complex<double> **tmp_v0,
                std::complex<double> **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void GC_child_CisAjtCkuAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int Asum,
                long unsigned int Adiff,
                std::complex<double> tmp_V,
                int nstate, std::complex<double> **tmp_v0,
                std::complex<double> **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void GC_child_CisAjtCkuAlv_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Asum,
                long unsigned int Adiff,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                std::complex<double> tmp_V,
                int nstate, std::complex<double> **tmp_v0,
                std::complex<double> **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off_2
        );
//[e]Grand canonical

void GC_CisAis
(
 long unsigned int j,
 int nstate,
 std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 std::complex<double> tmp_trans
);

void GC_AisCis(
        long unsigned int j,
        int nstate, std::complex<double> **tmp_v0,
        std::complex<double> **tmp_v1,
        struct BindStruct *X,
        long unsigned int is1_spin,
        std::complex<double> tmp_trans
);

int X_CisAis
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin
 );

int X_CisAjt
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 long unsigned int *tmp_off
 );


int X_GC_CisAjt
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 long unsigned int *tmp_off
 );


void CisAjt
(
 long unsigned int j,
 int nstate, 
 std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 std::complex<double> tmp_V
 );


void GC_CisAjt
(
 long unsigned int j,
 int nstate, 
 std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 std::complex<double> tmp_V,
 long unsigned int *tmp_off
);


int child_general_hopp_GetInfo
(
 struct BindStruct *X,
 unsigned long int isite1,
 unsigned long int isite2,
 unsigned long int sigma1,
 unsigned long int sigma2
 );

int child_general_int_GetInfo
(
 int iInterAll,
 struct BindStruct *X,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isite3,
 long unsigned int isite4,
 long unsigned int sigma1,
 long unsigned int sigma2,
 long unsigned int sigma3,
 long unsigned int sigma4,
 std::complex<double> tmp_V
 );


int child_pairhopp_GetInfo
(
 int iPairHopp,
 struct BindStruct *X 
  );

int child_exchange_GetInfo
(
 int iExchange,
 struct BindStruct *X 
 );


void GC_Ajt
(
 long unsigned int j,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 long unsigned int is1_spin,
 std::complex<double> tmp_V,
 long unsigned int *tmp_off
 );

void GC_Cis
(
 long unsigned int j,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 long unsigned int is1_spin,
 std::complex<double> tmp_V,
 long unsigned int *tmp_off
 );



void GC_Ajt
(
 long unsigned int j,
 int nstate, std::complex<double> **tmp_v0,
 std::complex<double> **tmp_v1,
 long unsigned int is1_spin,
 std::complex<double> tmp_V,
 long unsigned int *tmp_off
 );

int X_Cis
(
 long unsigned int j,
 long unsigned int is1_spin,
 long unsigned int *tmp_off,
 long unsigned int *list_1_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target,
 long unsigned int _irght,
 long unsigned int _ilft,
 long unsigned int _ihfbit
 );



int X_Ajt
(
 long unsigned int j,
 long unsigned int is1_spin,
 long unsigned int *tmp_off,
 long unsigned int *list_1_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target,
 long unsigned int _irght,
 long unsigned int _ilft,
 long unsigned int _ihfbit
 );

#endif
