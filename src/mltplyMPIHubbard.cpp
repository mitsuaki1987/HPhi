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
@brief Functions for Hubbard Hamiltonian + MPI
*/
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "mltplyCommon.hpp"
#include "mltplyMPIHubbard.hpp"
#include "global.hpp"
/**
@brief Hopping term in Hubbard + GC
When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Hubbard::GC::general_hopp_MPIdouble
(
 long int itrans,//!<[in] Transfer ID
 int nstate, 
 std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
 std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  mltply::Hubbard::GC::X_general_hopp_MPIdouble(
    Def::EDGeneralTransfer[itrans][0], Def::EDGeneralTransfer[itrans][1],
    Def::EDGeneralTransfer[itrans][2], Def::EDGeneralTransfer[itrans][3],
    Def::EDParaGeneralTransfer[itrans], nstate, tmp_v0, tmp_v1);
}/*void mltply::Hubbard::GC::general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard + GC
When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
@return fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
*/
void mltply::Hubbard::GC::X_general_hopp_MPIdouble(
  int org_isite1,//!<[in] @f$i_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin1,//!<[in] @f$\sigma_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_isite2,//!<[in] @f$i_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin2,//!<[in] @f$\sigma_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  std::complex<double> tmp_trans,//!<[in] Transfer @f$t@f$
  int nstate, 
  std::complex<double> **tmp_v0,//!< [out] Result v0 = H v1
  std::complex<double> **tmp_v1 //!< [in] v0 = H v1
) {
  int mask1, mask2, state1, state2, origin, bitdiff, Fsgn;
  std::complex<double> trans;

  mask1 = (int)Def::Tpow[2 * org_isite1 + org_ispin1];
  mask2 = (int)Def::Tpow[2 * org_isite2 + org_ispin2];
  if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
  else bitdiff = mask1 - mask2 * 2;
  origin = MP::myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

  SgnBit((long int) (origin & bitdiff), &Fsgn); // Fermion sign

  if (state1 == 0 && state2 == mask2) {
    trans = -(double)Fsgn * tmp_trans;
  }/*if (state1 == 0 && state2 == mask2)*/
  else if (state1 == mask1 && state2 == 0) {
    trans = -(double)Fsgn * conj(tmp_trans);
    if (Large::mode == M_CORR || Large::mode == M_CALCSPEC) trans = 0.0;
  }/*if (state1 == mask1 && state2 == 0)*/
  else return;

  wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);

  zaxpy_long(Check::idim_max*nstate, trans, &Wave::v1buf[0][0], &tmp_v0[0][0]);
}/*void mltply::Hubbard::GC::general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard + MPI
When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
@return fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
*/
void mltply::Hubbard::C::X_CisAjt_MPIdouble(
  int org_isite1,//!<[in] @f$i_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin1,//!<[in] @f$\sigma_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_isite2,//!<[in] @f$i_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin2,//!<[in] @f$\sigma_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  std::complex<double> tmp_trans,//!<[in] Transfer @f$t@f$
  int nstate,
  std::complex<double> **tmp_v0,//!< [out] Result v0 = H v1
  std::complex<double> **tmp_v1//!< [in] v0 = H v1
) {
  int mask1, mask2, state1, state2, origin, bitdiff, Fsgn;
  long int idim_max_buf, j, ioff;
  std::complex<double> trans;
  int one = 1;

  mask1 = (int) Def::Tpow[2 * org_isite1 + org_ispin1];
  mask2 = (int) Def::Tpow[2 * org_isite2 + org_ispin2];
  if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
  else bitdiff = mask1 - mask2 * 2;
  origin = MP::myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

  SgnBit((long int) (origin & bitdiff), &Fsgn); // Fermion sign

  if (state1 == 0 && state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
  }/*if (state1 == 0 && state2 == mask2)*/
  else if (state1 == mask1 && state2 == 0) {
    trans = -(double) Fsgn * conj(tmp_trans);
    if (Large::mode == M_CORR|| Large::mode == M_CALCSPEC) {
      trans = 0;
    }
  }/*if (state1 == mask1 && state2 == 0)*/
  else return;

  idim_max_buf = wrapperMPI::SendRecv_i(origin, Check::idim_maxOrg);
  wrapperMPI::SendRecv_iv(origin, Check::idim_maxOrg, idim_max_buf, List::c1_org, List::c1buf_org);
  wrapperMPI::SendRecv_cv(origin, Check::idim_maxOrg*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);
  
#pragma omp parallel for default(none) private(j, ioff) \
shared(idim_max_buf, trans,List::c2_1,List::c2_2,List::c1buf_org,Wave::v1buf, tmp_v0, \
nstate,one,Large::irght, Large::ilft, Large::ihfbit)
  for (j = 0; j < idim_max_buf; j++) {
    GetOffComp(List::c2_1, List::c2_2, List::c1buf_org[j],
               Large::irght, Large::ilft, Large::ihfbit, &ioff);
    zaxpy_(&nstate, &trans, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 0; j < idim_max_buf; j++)*/
}/*void child_CisAjt_MPIdouble*/
/**
@brief Hopping term in Hubbard + GC
 When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Hubbard::GC::general_hopp_MPIsingle(
  long int itrans,//!<[in] Transfer ID
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  mltply::Hubbard::GC::X_general_hopp_MPIsingle(
    Def::EDGeneralTransfer[itrans][0], Def::EDGeneralTransfer[itrans][1],
    Def::EDGeneralTransfer[itrans][2], Def::EDGeneralTransfer[itrans][3],
    Def::EDParaGeneralTransfer[itrans], nstate, tmp_v0, tmp_v1       );
}/*void mltply::Hubbard::GC::general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard + GC
 When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltply::Hubbard::GC::X_general_hopp_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Hopping integral
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
) {
  int mask2, state1, state2, origin, bit2diff, Fsgn;
  long int j, mask1, state1check, bit1diff, ioff;
  std::complex<double> trans, dmv;
  int one = 1;
  /*
    Prepare index in the inter PE
  */
  mask2 = (int) Def::Tpow[2 * org_isite2 + org_ispin2];
  bit2diff = mask2 - 1;
  origin = MP::myrank ^ mask2;
  state2 = origin & mask2;

  SgnBit((long int) (origin & bit2diff), &Fsgn); // Fermion sign

  wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, Check::idim_max*nstate,
    &tmp_v1[0][0], &Wave::v1buf[0][0]);

  /*
    Index in the intra PE
  */
  mask1 = Def::Tpow[2 * org_isite1 + org_ispin1];

  if (state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
    state1check = 0;
  }/*if (state2 == mask2)*/
  else if (state2 == 0) {
    state1check = mask1;
    trans = -(double) Fsgn * conj(tmp_trans);
    if (Large::mode == M_CORR|| Large::mode == M_CALCSPEC) trans = 0;
  }/*if (state2 != mask2)*/
  else return;

  bit1diff = Def::Tpow[2 * Def::Nsite - 1] * 2 - mask1 * 2;

#pragma omp parallel default(none) private(j,dmv,state1,Fsgn,ioff) \
shared(Check::idim_max,trans,mask1,state1check,bit1diff,Wave::v1buf,tmp_v1,tmp_v0,nstate,one)
  {
#pragma omp for
    for (j = 0; j < Check::idim_max; j++) {

      state1 = j & mask1;

      if (state1 == state1check) {

        SgnBit(j & bit1diff, &Fsgn);
        ioff = j ^ mask1;

        dmv = (double)Fsgn * trans;
        zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
      }/*if (state1 == state1check)*/
    }/*for (j = 0; j < idim_max_buf; j++)*/

  }/*End of parallel region*/
}/*void mltply::Hubbard::GC::general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Hubbard::C::general_hopp_MPIdouble(
  long int itrans,//!<[in] Transfer ID
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  mltply::Hubbard::C::X_general_hopp_MPIdouble( 
    Def::EDGeneralTransfer[itrans][0], Def::EDGeneralTransfer[itrans][1],
    Def::EDGeneralTransfer[itrans][2], Def::EDGeneralTransfer[itrans][3],
    Def::EDParaGeneralTransfer[itrans], nstate, tmp_v0, tmp_v1);
}/*void mltply::Hubbard::C::general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Hubbard::C::X_general_hopp_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Hopping integral
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
) {
  int mask1, mask2, state1, state2, origin, bitdiff, Fsgn;
  long int idim_max_buf, j, ioff;
  std::complex<double> trans;
  int one = 1;

  mask1 = (int) Def::Tpow[2 * org_isite1 + org_ispin1];
  mask2 = (int) Def::Tpow[2 * org_isite2 + org_ispin2];

  if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
  else bitdiff = mask1 - mask2 * 2;
  origin = MP::myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

  SgnBit((long int) (origin & bitdiff), &Fsgn); // Fermion sign

  if (state1 == 0 && state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
  }
  else if (state1 == mask1 && state2 == 0) {
    trans = -(double) Fsgn * conj(tmp_trans);
    if (Large::mode == M_CORR|| Large::mode == M_CALCSPEC) trans = 0;
  }
  else return;

  idim_max_buf = wrapperMPI::SendRecv_i(origin, Check::idim_max);
  wrapperMPI::SendRecv_iv(origin, Check::idim_max, idim_max_buf, List::c1, List::c1buf);
  wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate,
    &tmp_v1[0][0], &Wave::v1buf[0][0]);

#pragma omp parallel default(none) private(j,ioff) \
shared(idim_max_buf,trans,Large::irght, Large::ilft, Large::ihfbit, \
List::c2_1,List::c2_2,List::c1buf,Wave::v1buf,tmp_v1,tmp_v0,nstate,one)
  {
#pragma omp for
    for (j = 0; j < idim_max_buf; j++) {
      GetOffComp(List::c2_1, List::c2_2, List::c1buf[j],
                 Large::irght, Large::ilft, Large::ihfbit, &ioff);
      zaxpy_(&nstate, &trans, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
    }/*for (j = 0; j < idim_max_buf; j++)*/
  }/*End of parallel region*/
}/*void mltply::Hubbard::C::general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Hubbard::C::general_hopp_MPIsingle(
  long int itrans,//!<[in] Transfer ID
  int nstate,
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  mltply::Hubbard::C::X_general_hopp_MPIsingle(
    Def::EDGeneralTransfer[itrans][0], Def::EDGeneralTransfer[itrans][1],
    Def::EDGeneralTransfer[itrans][2], Def::EDGeneralTransfer[itrans][3],
    Def::EDParaGeneralTransfer[itrans], nstate, tmp_v0, tmp_v1);
}/*void mltply::Hubbard::C::general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Hubbard::C::X_general_hopp_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Hopping integral
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
) {
  int mask2, state2, origin, bit2diff, Fsgn;
  long int mask1, state1, idim_max_buf, j, state1check, bit1diff, ioff, jreal;
  std::complex<double> trans, dmv;
  int one = 1;
  /*
    Prepare index in the inter PE
  */
  mask2 = (int)Def::Tpow[2 * org_isite2+org_ispin2];
  bit2diff = mask2 - 1;
  origin = MP::myrank ^ mask2;

  state2 = origin & mask2;

  SgnBit((long int) (origin & bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = wrapperMPI::SendRecv_i(origin, Check::idim_max);
  wrapperMPI::SendRecv_iv(origin, Check::idim_max, idim_max_buf, List::c1, List::c1buf);
  wrapperMPI::SendRecv_cv(origin, Check::idim_max*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);
  /*
    Index in the intra PE
  */
  mask1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  if (state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
    state1check = 0;
  }
  else if (state2 == 0) {
    state1check = mask1;
    trans = -(double) Fsgn * conj(tmp_trans);
    if (Large::mode == M_CORR|| Large::mode == M_CALCSPEC) {
      trans = 0;
    }
  }
  else return;

  bit1diff = Def::Tpow[2 * Def::Nsite - 1] * 2 - mask1 * 2;

#pragma omp parallel default(none) private(j,dmv,Fsgn,ioff,jreal,state1) \
shared(idim_max_buf,trans,mask1,state1check,bit1diff,MP::myrank,Large::irght, Large::ilft, Large::ihfbit, \
List::c1,List::c2_1,List::c2_2,List::c1buf,Wave::v1buf,tmp_v1,tmp_v0,nstate,one)
  {
#pragma omp for
    for (j = 0; j < idim_max_buf; j++) {

      jreal = List::c1buf[j];
      state1 = jreal & mask1;

      if (state1 == state1check) {
        SgnBit(jreal & bit1diff,&Fsgn);
        GetOffComp(List::c2_1, List::c2_2, jreal ^ mask1,
            Large::irght, Large::ilft, Large::ihfbit, &ioff);

        dmv = (double)Fsgn * trans;
        zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
      }/*if (state1 == state1check)*/
    }/*for (j = 0; j < idim_max_buf; j++)*/
  }/*End of parallel region*/
}/*std::complex<double> mltply::Hubbard::C::general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
  When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void mltply::Hubbard::C::X_CisAjt_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  std::complex<double> tmp_trans,//!<[in] Hopping integral
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1//!<[in] v0 = H v1
){
  int mask2, state2, origin, bit2diff, Fsgn;
  long int mask1, state1, idim_max_buf, j, state1check, bit1diff, ioff, jreal;
  std::complex<double> trans, dmv;
  int one = 1;
  /*
    Prepare index in the inter PE
  */
  mask2 = (int)Def::Tpow[2 * org_isite2+org_ispin2];
  bit2diff = mask2 - 1;
  origin = MP::myrank ^ mask2;

  state2 = origin & mask2;

  SgnBit((long int) (origin & bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = wrapperMPI::SendRecv_i(origin, Check::idim_maxOrg);
  wrapperMPI::SendRecv_iv(origin, Check::idim_maxOrg, idim_max_buf, List::c1_org, List::c1buf_org);
  wrapperMPI::SendRecv_cv(origin, Check::idim_maxOrg*nstate, idim_max_buf*nstate, &tmp_v1[0][0], &Wave::v1buf[0][0]);
  /*
    Index in the intra PE
  */
  mask1 = Def::Tpow[2 * org_isite1 + org_ispin1];
  if (state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
    state1check = 0;
  }
  else if (state2 == 0) {
    state1check = mask1;
    trans = -(double) Fsgn * conj(tmp_trans);
  }
  else return;

  bit1diff = Def::Tpow[2 * Def::Nsite - 1] * 2 - mask1 * 2;

#pragma omp parallel for default(none) private(j,dmv,Fsgn,ioff,jreal,state1) \
shared(idim_max_buf,trans,mask1,state1check,bit1diff,List::c2_1,List::c2_2,List::c1buf_org, \
List::c1,Wave::v1buf, tmp_v0,nstate,one,Large::irght, Large::ilft, Large::ihfbit)
  for (j = 0; j < idim_max_buf; j++) {
    jreal = List::c1buf_org[j];
    state1 = jreal & mask1;
    if (state1 == state1check) {
      SgnBit(jreal & bit1diff, &Fsgn);
      GetOffComp(List::c2_1, List::c2_2, jreal ^ mask1,
        Large::irght, Large::ilft, Large::ihfbit, &ioff);
      if (ioff != 0) {
        dmv = (double)Fsgn * trans;
        zaxpy_(&nstate, &dmv, &Wave::v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
      }/*if(ioff !=0)*/
    }/*if (state1 == state1check)*/
  }/*for (j = 0; j < idim_max_buf; j++)*/
}/*std::complex<double> mltply::Hubbard::C::general_hopp_MPIsingle*/
