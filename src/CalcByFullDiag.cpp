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
#include "CalcByFullDiag.hpp"
#include "input.hpp"
#include "wrapperMPI.hpp"
#include "CalcTime.hpp"
#include "mltplyCommon.hpp"
#include "lapack_diag.hpp"
#include "phys.hpp"
#include "output.hpp"
#include "mltply.hpp"
#include "global.hpp"

/**
@brief Parent function for FullDiag mode
@retval TRUE(=1) normally finished.
@retval FALSE(=0) abnormally finished.
*/
int CalcByFullDiag::main()
{
  int iret = 0;
  long int idim;

  fprintf(MP::STDOUT, "%s", "######  Start: Setting Hamiltonian.  ######\n\n");
  StartTimer(5100);
  if (Def::iInputHam == FALSE) {
    zclear(Check::idim_max * Check::idim_max, &Wave::v0[0][0]);
    zclear(Check::idim_max * Check::idim_max, &Wave::v1[0][0]);
    for (idim = 0; idim < Check::idim_max; idim++) Wave::v1[idim][idim] = 1.0;
    mltply::main(Check::idim_max, Wave::v0, Wave::v1, Check::idim_max, List::a1, List::a2_1, List::a2_2, List::Diagonal);
  }
  else if (Def::iInputHam == TRUE) {
    fprintf(MP::STDOUT, "%s", "######  Start: Input Hamiltonian.  ######\n\n");
    inputHam();
    fprintf(MP::STDOUT, "%s", "######  End  : Input Hamiltonian.  ######\n\n");
  }
  StopTimer(5100);
  fprintf(MP::STDOUT, "%s", "######  End  : Setting Hamiltonian.  ######\n\n");
  if (iret != 0) return FALSE;


  if (Def::iOutputHam == TRUE) {
    fprintf(MP::STDOUT, "%s", "######  Start: Output Hamiltonian.  ######\n\n");
    StartTimer(5500);
    iret = outputHam();
    StopTimer(5500);
    fprintf(MP::STDOUT, "%s", "######  End  : Output Hamiltonian.  ######\n\n");
    if (iret != 0) return FALSE;
    return TRUE;
  }

  fprintf(MP::STDOUT, "%s", "######  Start: Diagonalization.  ######\n\n");
  StartTimer(5200);
  iret = lapack_diag(Phys::energy);
  StopTimer(5200);
  fprintf(MP::STDOUT, "%s", "######  End  : Diagonalization.  ######\n\n");
  if (iret != 0) return FALSE;

  Def::St = 0;
  fprintf(MP::STDOUT, "%s", "######  Start: Calc Expected value.  ######\n\n");
  StartTimer(5300);
  phys(Check::idim_max);
  StopTimer(5300);
  fprintf(MP::STDOUT, "%s", "######  End  : Calc Expected value.  ######\n\n");

  StartTimer(5400);
  iret = output();
  StopTimer(5400);
  fprintf(MP::STDOUT, "%s", "######  Finish Calculation.  ######\n");
  if (iret != 0) return FALSE;

  return TRUE;
}
