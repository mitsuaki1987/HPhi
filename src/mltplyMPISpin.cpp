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
/**@file
@brief Functions for spin Hamiltonian + MPI
*/

#include "mltplyCommon.hpp"
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "mltplyMPISpin.hpp"
#include "mltplyMPISpinCore.hpp"
#include "global.hpp"
/**
@brief Exchange term in Spin model
  When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::C::Half::general_int_MPIdouble(
  long int i_int,//!<[in] Interaction ID
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1,//!<[in] v0 = H v1
  long int i_max,
  long int* list_1, 
  long int* list_2_1,
  long int* list_2_2
){
  mltply::Spin::C::Half::X_general_int_MPIdouble(
    (int)Def::InterAll_OffDiagonal[i_int][0], (int)Def::InterAll_OffDiagonal[i_int][1],
    (int)Def::InterAll_OffDiagonal[i_int][3], (int)Def::InterAll_OffDiagonal[i_int][4],
    (int)Def::InterAll_OffDiagonal[i_int][5], (int)Def::InterAll_OffDiagonal[i_int][7],
    Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1, i_max, list_1, list_2_1, list_2_2);
}/*void child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::C::Half::X_general_int_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1,//!<[in] Vector to be producted
  long int i_max, 
  long int* list_1,
  long int* list_2_1,
  long int* list_2_2
) {
  int mask1, mask2, state1, state2, origin;
  long int idim_max_buf, j, ioff;
  std::complex<double> Jint;
  int one = 1;

  mask1 = (int)Def::Tpow[org_isite1];
  mask2 = (int)Def::Tpow[org_isite3];
  origin = MP::myrank ^ (mask1 + mask2);

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == org_ispin2 && state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (state1 == org_ispin1 && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return;

  idim_max_buf = wrapperMPI::SendRecv_i(origin, i_max);
  wrapperMPI::SendRecv_iv(origin, i_max, idim_max_buf, list_1, List::buf);
  wrapperMPI::SendRecv_cv(origin, i_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

#pragma omp parallel for default(none) private(j, ioff) \
shared(idim_max_buf,Jint,Large::irght, Large::ilft, Large::ihfbit, \
list_2_1,list_2_2,List::buf,Wave::v1buf,tmp_v1,tmp_v0,nstate,one)
  for (j = 0; j < idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, List::buf[j],
        Large::irght, Large::ilft, Large::ihfbit, &ioff);
    zaxpy_(&nstate, &Jint, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 0; j < idim_max_buf; j++)*/
}/*std::complex<double> X_child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
  When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::C::Half::X_general_int_TotalS_MPIdouble(
  int org_isite1,//!<[in] site 1
  int org_isite3,//!<[in] site 3
  int nstate, 
  std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1,//!<[in] Vector to be producted
  long int i_max, 
  long int* list_1, 
  long int* list_2_1, 
  long int* list_2_2
){
  int mask1, mask2, num1_up, num2_up, origin, one = 1;
  long int idim_max_buf, j, ioff, ibit_tmp;
  std::complex<double> dmv;

  mask1 = (int)Def::Tpow[org_isite1];
  mask2 = (int)Def::Tpow[org_isite3];
  if (mask1 == mask2) origin = MP::myrank ^ mask1;
  else origin = MP::myrank ^ (mask1 + mask2);
  num1_up = (origin & mask1) / mask1;
  num2_up = (origin & mask2) / mask2;

  ibit_tmp = (num1_up) ^ (num2_up);
  if (ibit_tmp == 0) return;

  idim_max_buf = wrapperMPI::SendRecv_i(origin, i_max);
  wrapperMPI::SendRecv_iv(origin, i_max, idim_max_buf, list_1, List::buf);
  wrapperMPI::SendRecv_cv(origin, i_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

#pragma omp parallel for default(none)  private(j, dmv, ioff) \
shared(idim_max_buf, Large::irght, Large::ilft, Large::ihfbit, \
list_2_1, list_2_2, List::buf, Wave::v1buf, tmp_v1, tmp_v0,nstate,one)
  for (j = 0; j < idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, List::buf[j],
      Large::irght, Large::ilft, Large::ihfbit, &ioff);
    zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 0; j < idim_max_buf; j++)*/
  return;
}/*std::complex<double> X_child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
  When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::C::Half::general_int_MPIsingle(
  long int i_int,//!<[in] Interaction ID
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1,//!<[in] v0 = H v1
  long int i_max,
  long int* list_1, 
  long int* list_2_1, 
  long int* list_2_2
){

  mltply::Spin::C::Half::X_general_int_MPIsingle(
    (int)Def::InterAll_OffDiagonal[i_int][0], (int)Def::InterAll_OffDiagonal[i_int][1], 
    (int)Def::InterAll_OffDiagonal[i_int][3], (int)Def::InterAll_OffDiagonal[i_int][4],
    (int)Def::InterAll_OffDiagonal[i_int][5], (int)Def::InterAll_OffDiagonal[i_int][7],
    Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1, i_max, list_1, list_2_1, list_2_2);
}/*void child_general_int_spin_MPIsingle*/
/*
@brief General interaction term of canonical spin system.
site 3 is in the inter process region
*/
void mltply::Spin::C::Half::X_general_int_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  std::complex<double> tmp_J,//!<[in] Copupling constatnt
  int nstate,
  std::complex<double> **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  std::complex<double> **tmp_v1,//!<[in] Vector to be producted
  long int i_max,
  long int* list_1, 
  long int* list_2_1, 
  long int* list_2_2
) {
  int mask2, state2, origin;
  long int mask1, idim_max_buf, j, ioff, state1, jreal, state1check;
  std::complex<double> Jint;
  int one = 1;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)Def::Tpow[org_isite3];
  origin = MP::myrank ^ mask2;
  state2 = (origin & mask2) / mask2;

  if (state2 == org_ispin4) {
    state1check = (long int) org_ispin2;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (long int) org_ispin1;
    Jint = conj(tmp_J);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return;

  idim_max_buf = wrapperMPI::SendRecv_i(origin, i_max);
  wrapperMPI::SendRecv_iv(origin, i_max, idim_max_buf, list_1, List::buf);
  wrapperMPI::SendRecv_cv(origin, i_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);
  /*
  Index in the intra PE
  */
  mask1 = Def::Tpow[org_isite1];
  //printf("debug1 %ld\n", idim_max_buf);

#pragma omp parallel for default(none) private(j, ioff, jreal, state1) \
shared(idim_max_buf, Jint, mask1, state1check, org_isite1, \
Large::irght, Large::ilft, Large::ihfbit, \
list_2_1, list_2_2, List::buf, Wave::v1buf, tmp_v1, tmp_v0,nstate,one)
  for (j = 0; j < idim_max_buf; j++) {
    jreal = List::buf[j];

    state1 = (jreal & mask1) / mask1;
    if (state1 == state1check) {
      GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
        Large::irght, Large::ilft, Large::ihfbit, &ioff);
      zaxpy_(&nstate, &Jint, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
    }
  }
}/*std::complex<double> X_child_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::GC::Half::general_int_MPIdouble(
  long int i_int,//!<[in] Interaction ID
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  if (Def::InterAll_OffDiagonal[i_int][1] == Def::InterAll_OffDiagonal[i_int][3] &&
      Def::InterAll_OffDiagonal[i_int][5] != Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::Half::CisAisCjuAjv_MPIdouble(i_int, nstate, tmp_v0, tmp_v1);
  }
  else if (Def::InterAll_OffDiagonal[i_int][1] != Def::InterAll_OffDiagonal[i_int][3] &&
           Def::InterAll_OffDiagonal[i_int][5] == Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::Half::CisAitCjuAju_MPIdouble(i_int, nstate, tmp_v0, tmp_v1);
  }
  else {
    mltply::Spin::GC::Half::CisAitCiuAiv_MPIdouble(i_int, nstate, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::GC::Half::general_int_MPIsingle(
  long int i_int,//!<[in] Interaction ID
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  if (Def::InterAll_OffDiagonal[i_int][1] == Def::InterAll_OffDiagonal[i_int][3] &&
      Def::InterAll_OffDiagonal[i_int][5] != Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::Half::X_CisAisCjuAjv_MPIsingle(i_int, nstate, tmp_v0, tmp_v1);
  }
  else if (Def::InterAll_OffDiagonal[i_int][1] != Def::InterAll_OffDiagonal[i_int][3] &&
           Def::InterAll_OffDiagonal[i_int][5] == Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::Half::X_CisAitCjuAju_MPIsingle(i_int, nstate, tmp_v0, tmp_v1);
  }
  else {
    mltply::Spin::GC::Half::X_CisAitCiuAiv_MPIsingle(i_int, nstate, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::GC::General::general_int_MPIdouble(
  long int i_int,//!<[in] Interaction ID
  int nstate, std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  if (Def::InterAll_OffDiagonal[i_int][1] == Def::InterAll_OffDiagonal[i_int][3] &&
      Def::InterAll_OffDiagonal[i_int][5] != Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::General::X_CisAisCjuAjv_MPIdouble(
      Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
      Def::InterAll_OffDiagonal[i_int][4], Def::InterAll_OffDiagonal[i_int][5],
      Def::InterAll_OffDiagonal[i_int][7], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
  }
  else if (Def::InterAll_OffDiagonal[i_int][1] != Def::InterAll_OffDiagonal[i_int][3] &&
           Def::InterAll_OffDiagonal[i_int][5] == Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::General::X_CisAitCjuAju_MPIdouble(
      Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
      Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
      Def::InterAll_OffDiagonal[i_int][5], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
  }
  else {
    mltply::Spin::GC::General::X_CisAitCjuAjv_MPIdouble(
      Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1], 
      Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
      Def::InterAll_OffDiagonal[i_int][5], Def::InterAll_OffDiagonal[i_int][7],
      Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::GC::General::general_int_MPIsingle(
  long int i_int,//!<[in] Interaction ID
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){

  if (Def::InterAll_OffDiagonal[i_int][1] == Def::InterAll_OffDiagonal[i_int][3] &&
      Def::InterAll_OffDiagonal[i_int][5] != Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::General::X_CisAisCjuAjv_MPIsingle(
      Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
      Def::InterAll_OffDiagonal[i_int][4], Def::InterAll_OffDiagonal[i_int][5],
      Def::InterAll_OffDiagonal[i_int][7], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
  }
  else if (Def::InterAll_OffDiagonal[i_int][1] != Def::InterAll_OffDiagonal[i_int][3] &&
           Def::InterAll_OffDiagonal[i_int][5] == Def::InterAll_OffDiagonal[i_int][7]) {
    mltply::Spin::GC::General::X_CisAitCjuAju_MPIsingle(
      Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1], 
      Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
      Def::InterAll_OffDiagonal[i_int][5], Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
    }
  else {
    mltply::Spin::GC::General::X_CisAitCjuAjv_MPIsingle(
      Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1], 
      Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
      Def::InterAll_OffDiagonal[i_int][5], Def::InterAll_OffDiagonal[i_int][7],
      Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1);
  }

}/*void GC_child_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::C::General::general_int_MPIdouble(
  long int i_int,//!<[in] Interaction ID
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1,//!<[in] v0 = H v1
  long int i_max,
  long int* list_1, 
  long int* list_2_1,
  long int* list_2_2
){
    mltply::Spin::C::General::X_CisAitCjuAjv_MPIdouble(
      Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1], 
      Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
      Def::InterAll_OffDiagonal[i_int][5], Def::InterAll_OffDiagonal[i_int][7],
      Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1,
      i_max, list_1, list_2_1, list_2_2);

}/*void GC_child_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Spin::C::General::general_int_MPIsingle(
  long int i_int,//!<[in] Interaction ID
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1,//!<[in] v0 = H v1
  long int i_max, 
  long int* list_1, 
  long int* list_2_1, 
  long int* list_2_2
){

  mltply::Spin::C::General::X_CisAitCjuAjv_MPIsingle(
    Def::InterAll_OffDiagonal[i_int][0], Def::InterAll_OffDiagonal[i_int][1],
    Def::InterAll_OffDiagonal[i_int][3], Def::InterAll_OffDiagonal[i_int][4],
    Def::InterAll_OffDiagonal[i_int][5], Def::InterAll_OffDiagonal[i_int][7],
    Def::ParaInterAll_OffDiagonal[i_int], nstate, tmp_v0, tmp_v1, 
    i_max, list_1, list_2_1, list_2_2);

}/*void GC_child_general_int_spin_MPIsingle*/
