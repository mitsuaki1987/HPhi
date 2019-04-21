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

#ifndef HPHI_MLTPLYSPINCORE_H
#define HPHI_MLTPLYSPINCORE_H

#include "Common.hpp"

void child_exchange_spin_element(long int j, int nstate, std::complex<double> **tmp_v0, std::complex<double> **tmp_v1, struct BindStruct *X, long int *tmp_off);

void GC_child_pairlift_spin_element(long int j, int nstate, std::complex<double> **tmp_v0, std::complex<double> **tmp_v1, struct BindStruct *X, long int *tmp_off);
void GC_child_exchange_spin_element(long int j, int nstate, std::complex<double> **tmp_v0, std::complex<double> **tmp_v1, struct BindStruct *X, long int *tmp_off);
int X_child_exchange_spin_element(long int j, struct BindStruct *X, long int isA_up, long int isB_up, long int sigmaA, long int sigmaB, long int *tmp_off);
//[s]Spin
void child_CisAisCisAis_spin_element(long int j, long int isA_up, long int isB_up, long int org_sigma2, long int org_sigma4, std::complex<double> tmp_V, int nstate, std::complex<double> **tmp_v0, std::complex<double> **tmp_v1);
//[e]Spin

//[s]GC Spin
void GC_child_CisAisCisAis_spin_element(
  long int j, long int isA_up, long int isB_up, long int org_sigma2, long int org_sigma4, 
  std::complex<double> tmp_V, int nstate, std::complex<double> **tmp_v0,
  std::complex<double> **tmp_v1);
void GC_child_CisAisCitAiu_spin_element(
  long int j, long int org_sigma2, long int org_sigma4, long int isA_up, long int isB_up, 
  std::complex<double> tmp_V, int nstate, std::complex<double> **tmp_v0, 
  std::complex<double> **tmp_v1, long int *tmp_off);
void GC_child_CisAitCiuAiu_spin_element(
  long int j, long int org_sigma2, long int org_sigma4, long int isA_up, long int isB_up, 
  std::complex<double> tmp_V, int nstate, std::complex<double> **tmp_v0, 
  std::complex<double> **tmp_v1, long int *tmp_off);
void GC_child_CisAitCiuAiv_spin_element(
  long int j, long int org_sigma2, long int org_sigma4, long int isA_up, long int isB_up, 
  std::complex<double> tmp_V, int nstate, std::complex<double> **tmp_v0,
  std::complex<double> **tmp_v1, long int *tmp_off_2);
//[e]GC Spin

int child_general_int_spin_GetInfo(struct BindStruct *X, 
  long int isite1, long int isite2, long int sigma1, long int sigma2, long int sigma3, 
  long int sigma4, std::complex<double> tmp_V);
int child_exchange_spin_GetInfo(int iExchange, struct BindStruct *X);
int child_pairlift_spin_GetInfo(int iPairLift, struct BindStruct *X);
int X_SpinGC_CisAit(long int j,long int is1_spin,long int sigma2,long int *tmp_off);
int X_Spin_CisAit(long int j, struct BindStruct *X, long int is1_spin, long int sigma2, long int *tmp_off);
int X_Spin_CisAis(long int j,long int is1_spin,long int sigma1);
int X_SpinGC_CisAis(long int j,long int is1_spin,long int sigma1);

#endif
