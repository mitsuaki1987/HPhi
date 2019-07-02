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
#include "output_list.hpp"
#include "FileIO.hpp"
#include "wrapperMPI.hpp"
#include "global.hpp"
#include "DefCommon.hpp"

/**
 * @file
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @brief  Output list_1 for canonical ensembles obtained in sz.cpp.
 *
 */


/** 
 * @brief Output list_1 for canonical ensembles.
 * 
 * @param X [in] Struct to get the dimension of the target Hilbert space and the information for making the output file name.
 * 
 * @retval -1 fail to output the list.
 * @retval 0  success to output the list.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int output_list(){
  
  FILE *fp;
  char sdt[D_FileNameMax];
  int i,i_max;
  
  fprintf(stdoutMPI, "%s", "  Start: output list. \n");
  i_max=Check::idim_max;
  switch(Def::iCalcModel){
  case HubbardGC:
  case Hubbard:
  case Spin:
  case SpinGC:
    sprintf(sdt, "ListForModel_Ns%d_Nup%dNdown%d.dat", Def::Nsite,Def::Nup,Def::Ndown);
  break;
  case Kondo:
  case KondoGC:
    sprintf(sdt, "ListForKondo_Ns%d_Ncond%d", Def::Nsite,Def::Ne);
    break;
  default:
    return -1;
  }
  if(childfopenMPI(sdt,"w",&fp)!=0){
    return -1;
  }
  for(i=1;i<=i_max;i++){
    fprintf(fp," %lu \n",list_1[i]);
  }
  fclose(fp);
  
  fprintf(stdoutMPI, "%s", "  End  : output list. \n");
  return 0;
}
