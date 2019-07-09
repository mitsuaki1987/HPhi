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
#include "output.hpp"
#include "FileIO.hpp"
#include "wrapperMPI.hpp"
#include "DefCommon.hpp"
#include "global.hpp"
#include <cstring>

/**
\brief output function for FullDiag mode
\retval 0 normally finished.
\retval -1 abnormally finished.
*/
int output() {

  FILE *fp;
  char sdt[D_FileNameMax];
  long int i, i_max;
  i_max = Check::idim_max;

  if (Def::iCalcType == DC::FullDiag) {
    switch (Def::iCalcModel) {
      case DC::Spin:
      case DC::Hubbard:
      case DC::Kondo:
        sprintf(sdt, "%s_phys_Nup%d_Ndown%d.dat", Def::CDataFileHead, Def::Nup, Def::Ndown);
        break;
      case DC::SpinGC:
      case DC::HubbardGC:
      case DC::KondoGC:
        sprintf(sdt, "%s_phys.dat", Def::CDataFileHead);
        break;
      default:
        break;
    }
    if (childfopenMPI(sdt, "w", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "  <H>         <N>        <Sz>       <S2>       <D> \n");
    for (i = 0; i < i_max; i++) {
      fprintf(fp, " %10lf %10lf %10lf %10lf %10lf\n", Phys::energy[i], Phys::num_up[i]+Phys::num_down[i], Phys::Sz[i],
              Phys::s2[i], Phys::doublon[i]);
    }
    fclose(fp);
  }
  else{
    fprintf(MP::STDOUT, "Error: output function is used only for FullDiag mode.");
    return -1;
  }

  return 0;
}
/**
\brief output Hamiltonian only used for FullDiag mode
\note global: [in] Ham
\retval 0 normally finished.
\retval -1 abnormally finished.
*/
int outputHam(){
  long int i=0;
  long int j=0;
  long int imax = Check::idim_max;
  long int ihermite=0;
  char cHeader[256];
  FILE *fp;
  char sdt[D_FileNameMax];

#pragma omp parallel for default(none) reduction(+:ihermite) private(i, j) \
shared(Wave::v0,imax)
  for (i = 1; i <= imax; i++) {
    for (j = 1; j <= i; j++) {
      if (abs(Wave::v0[i][j]) > 1.0e-13) {
        ihermite += 1;
      }
    }
  }

  strcpy(cHeader, "%%%%MatrixMarket matrix coordinate complex hermitian\n");
  sprintf(sdt,"%s_Ham.dat", Def::CDataFileHead);
  if(childfopenMPI(sdt,"w",&fp)!=0){
    return -1;
  }
  fprintf(fp, "%s", cHeader);
  fprintf(fp, "%ld %ld %ld \n", imax, imax, ihermite);
  for (i=1; i<=imax; i++){
    for (j=1; j<=i; j++){
      if(abs(Wave::v0[i][j])>1.0e-13){
        fprintf(fp, "%ld %ld %lf %lf\n",i,j,real(Wave::v0[i][j]),imag(Wave::v0[i][j]));
      }
    }
  }
  fclose(fp);
  return 0;
}
