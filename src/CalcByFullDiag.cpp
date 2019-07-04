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
int CalcByFullDiag::CalcByFullDiag(
)
{
  int iret=0;
  long int idim;

  fprintf(stdoutMPI, "%s", "######  Start: Setting Hamiltonian.  ######\n\n");
  StartTimer(5100);
  if(Def::iInputHam==FALSE){
    zclear((Check::idim_max + 1)*Check::idim_max, &v0[0][0]);
    zclear((Check::idim_max + 1)*Check::idim_max, &v1[0][0]);
    for (idim = 1; idim <= Check::idim_max; idim++) v1[idim][idim-1] = 1.0;
    mltply(Check::idim_max, v0, v1);
  }
  else if(Def::iInputHam==TRUE){
    fprintf(stdoutMPI, "%s", "######  Start: Input Hamiltonian.  ######\n\n");
    inputHam();
    fprintf(stdoutMPI, "%s", "######  End  : Input Hamiltonian.  ######\n\n");
  }
  StopTimer(5100);
  fprintf(stdoutMPI, "%s", "######  End  : Setting Hamiltonian.  ######\n\n");
  if(iret != 0) return FALSE;


  if(Def::iOutputHam == TRUE){
    fprintf(stdoutMPI, "%s", "######  Start: Output Hamiltonian.  ######\n\n");
    StartTimer(5500);
    iret=outputHam();
    StopTimer(5500);
    fprintf(stdoutMPI, "%s", "######  End  : Output Hamiltonian.  ######\n\n");
    if(iret != 0) return FALSE;
    return TRUE;
  }

  fprintf(stdoutMPI, "%s", "######  Start: Diagonalization.  ######\n\n");
  StartTimer(5200);
  iret=lapack_diag();
  StopTimer(5200);
  fprintf(stdoutMPI, "%s", "######  End  : Diagonalization.  ######\n\n");
  if(iret != 0) return FALSE;

  Def::St=0;
  fprintf(stdoutMPI, "%s", "######  Start: Calc Expected value.  ######\n\n");
  StartTimer(5300);
  phys(Check::idim_max);
  StopTimer(5300);
  fprintf(stdoutMPI, "%s", "######  End  : Calc Expected value.  ######\n\n");

  StartTimer(5400);
  iret=output();
  StopTimer(5400);  
  fprintf(stdoutMPI, "%s", "######  Finish Calculation.  ######\n");
  if(iret != 0) return FALSE;

  return TRUE;
}
