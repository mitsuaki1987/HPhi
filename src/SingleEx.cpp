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
@brief Functions to compute singly excited state
*/
#include "bitcalc.hpp"
#include "SingleEx.hpp"
#include "SingleExHubbard.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
/**
@brief Calculation of single excited state
Target System: Hubbard, Kondo
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetExcitedState::Single::main(
  int nstate, 
  std::complex<double> **tmp_v0,//!<[out] Result v0 = H v1
  std::complex<double> **tmp_v1,//!<[in] v0 = H v1
  int iEx
) {
  int iret = 0;
  //tmp_v0
  if (Def::NSingleExcitationOperator == 0) return TRUE;

  switch (Def::iCalcModel) {
  case DC::HubbardGC:
    iret = GetExcitedState::Single::Hubbard::GC(nstate, tmp_v0, tmp_v1, iEx);
    break;

  case DC::KondoGC:
  case DC::Hubbard:
  case DC::Kondo:
    iret = GetExcitedState::Single::Hubbard::C(nstate, tmp_v0, tmp_v1, iEx);
    break;

  case DC::Spin:
  case DC::SpinGC:
    iret = FALSE;
    break;

  default:
    iret = FALSE;
    break;
  }/*switch (Def::iCalcModel)*/
  return iret;
}/*int GetSingleExcitedState*/
