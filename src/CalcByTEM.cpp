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
#include "readdef.hpp"
#include "FirstMultiply.hpp"
#include "Multiply.hpp"
#include "diagonalcalc.hpp"
#include "expec_energy_flct.hpp"
#include "expec_cisajs.hpp"
#include "expec_cisajscktaltdc.hpp"
#include "CalcByTEM.hpp"
#include "FileIO.hpp"
#include "wrapperMPI.hpp"
#include "HPhiTrans.hpp"
#include "common/setmemory.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include "log.hpp"
#include <iostream>
/**@file
 * @brief  File to define functions to calculate expected values by Time evolution method.
*/
 /**\brief 
 Set transfer integrals at timeidx-th time
 */
void MakeTEDTransfer( 
  const int timeidx/**<index of time*/
) {
  int i, j;
  //Clear values
  for (i = 0; i < Def::NTETransferMax; i++) {
    for (j = 0; j < 4; j++) {
      Def::EDGeneralTransfer[i + Def::EDNTransfer][j] = 0;
    }
    Def::EDParaGeneralTransfer[i + Def::EDNTransfer] = 0.0;
  }

  //Input values
  for (i = 0; i < Def::NTETransfer[timeidx]; i++) {
    for (j = 0; j < 4; j++) {
      Def::EDGeneralTransfer[i + Def::EDNTransfer][j] = Def::TETransfer[timeidx][i][j];
    }
    Def::EDParaGeneralTransfer[i + Def::EDNTransfer] = Def::ParaTETransfer[timeidx][i];
  }
  Def::EDNTransfer += Def::NTETransfer[timeidx];
}
/**\brief 
Set interall interactions at timeidx-th time
*/
void MakeTEDInterAll( 
  const int timeidx/**<index of time*/
) {
  int i, j;
  //Clear values
  for (i = 0; i < Def::NTEInterAllMax; i++) {
    for (j = 0; j < 8; j++) {
      Def::InterAll_OffDiagonal[i + Def::NInterAll_OffDiagonal][j] = 0;
    }
    Def::ParaInterAll_OffDiagonal[i + Def::NInterAll_OffDiagonal] = 0.0;
  }
  //Input values
  for (i = 0; i < Def::NTEInterAllOffDiagonal[timeidx]; i++) {
    for (j = 0; j < 8; j++) {
      Def::InterAll_OffDiagonal[i + Def::NInterAll_OffDiagonal][j] = Def::TEInterAllOffDiagonal[timeidx][i][j];
    }
    Def::ParaInterAll_OffDiagonal[i + Def::NInterAll_OffDiagonal] = Def::ParaTEInterAllOffDiagonal[timeidx][i];
  }
  Def::NInterAll_OffDiagonal += Def::NTEInterAllOffDiagonal[timeidx];
}
/** 
 * @brief main function of time evolution calculation
 * 
 * @return 0 normally finished
 * @return -1 unnormally finished
 *
 * @author Kota Ido (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CalcByTEM(
  const int ExpecInterval/**<[in] interval to output expected values*/
) {
  size_t byte_size;
  char *defname;
  char sdt[D_FileNameMax];
  char sdt_phys[D_FileNameMax];
  char sdt_norm[D_FileNameMax];
  char sdt_flct[D_FileNameMax];
  int rand_i = 0;
  int step_initial = 0;
  long int i_max = 0;
  FILE *fp;
  double Time = Param::Tinit;
  double dt = ((Def::NLaser == 0) ? 0.0 : Param::TimeSlice);
  std::complex<double> **v2;  /**< Ttemporary vector for time evolution calculation, @f$ v2 = H*v1 = H^coef |psi(t)>@f$.*/

  global_norm = d_1d_allocate(1);

  if (Def::NTETimeSteps < Def::Lanczos_max) {
    fprintf(stdoutMPI, "Error: NTETimeSteps must be larger than Lanczos_max.\n");
    return -1;
  }
  step_spin = ExpecInterval;
  Def::St = 0;
  fprintf(stdoutMPI, "%s", "######  Start: TimeEvolution by Taylor Expansion.  ######\n\n");
  if (Def::iInputEigenVec == FALSE) {
    fprintf(stderr, "Error: A file of Inputvector is not inputted.\n");
    return -1;
  }
  else {
    //input v1
    fprintf(stdoutMPI, "%s", "An Initial Vector is inputted.\n");
    TimeKeeper("%s_TimeKeeper.dat", "Reading an input Eigenvector starts: %s", "a");
    GetFileNameByKW(KWSpectrumVec, &defname);
    strcat(defname, "_rank_%d.dat");
    sprintf(sdt, defname, myrank);
    childfopenALL(sdt, "rb", &fp);
    if (fp == NULL) {
      fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
      fclose(fp);
      exitMPI(-1);
    }
    byte_size = fread(&step_initial, sizeof(int), 1, fp);
    byte_size = fread(&i_max, sizeof(long int), 1, fp);
    if (i_max != Check::idim_max) {
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      fclose(fp);
      printf("byte_size : %d\n", (int)byte_size);
      exitMPI(-1);
    }
    byte_size = fread(&v1[0][0], sizeof(std::complex<double>), Check::idim_max + 1, fp);
    fclose(fp);
    if (Def::iReStart == RESTART_NOT || Def::iReStart == RESTART_OUT) {
      step_initial = 0;
    }
  }

  sprintf(sdt_phys, "%s", "SS.dat");
  if (childfopenMPI(sdt_phys, "w", &fp) != 0) {
    return -1;
  }
  fprintf(fp, "%s", " # time, energy, phys_var, phys_doublon, phys_num, step_i\n");
  fclose(fp);

  sprintf(sdt_norm, "%s", "Norm.dat");
  if (childfopenMPI(sdt_norm, "w", &fp) != 0) {
    return -1;
  }
  fprintf(fp, "%s", " # time, norm, step_i \n");
  fclose(fp);

  sprintf(sdt_flct, "%s", "Flct.dat");
  if (childfopenMPI(sdt_flct, "w", &fp) != 0) {
    return -1;
  }
  fprintf(fp, "%s", " # time, N, N^2, D, D^2, Sz, Sz^2, step_i \n");
  fclose(fp);


  int iInterAllOffDiagonal_org = Def::NInterAll_OffDiagonal;
  int iTransfer_org = Def::EDNTransfer;
  v2 = cd_2d_allocate(Check::idim_max + 1, 1);
  for (step_i = step_initial; step_i < Def::Lanczos_max; step_i++) {
    Def::istep = step_i;

    //Reset total number of interactions (changed in MakeTED***function.)
    Def::EDNTransfer = iTransfer_org;
    Def::NInterAll_OffDiagonal = iInterAllOffDiagonal_org;

    if (step_i % (Def::Lanczos_max / 10) == 0) {
      fprintf(stdoutMPI, "    step_i/total_step = %d/%d \n", step_i, Def::Lanczos_max);
    }

    if (Def::NLaser != 0) {
      TransferWithPeierls(Time);
    }
    else {
      // common procedure
      Time = Def::TETime[step_i];
      if (step_i == 0) dt = 0.0;
      else {
        dt = Def::TETime[step_i] - Def::TETime[step_i - 1];
      }
      Param::TimeSlice = dt;

      // Set interactions
      if (Def::NTETransferMax != 0 && Def::NTEInterAllMax != 0) {
        fprintf(stdoutMPI,
          "Error: Time Evolution mode does not support TEOneBody and TETwoBody interactions at the same time. \n");
        return -1;
      }
      else if (Def::NTETransferMax > 0) { //One-Body type
        MakeTEDTransfer(step_i);
      }
      else if (Def::NTEInterAllMax > 0) { //Two-Body type
        MakeTEDInterAll(step_i);
      }
      //[e] Yoshimi
    }

    if (step_i == step_initial) {
      TimeKeeperWithStep("%s_Time_TE_Step.dat", "step %d:TE begins: %s", "w", step_i);
    }
    else {
      TimeKeeperWithStep("%s_Time_TE_Step.dat", "step %d:TE begins: %s", "a", step_i);
    }
    MultiplyForTEM(v2);
    //Add Diagonal Parts
    //Multiply Diagonal
    expec_energy_flct(1, v0, v1);

    if (Def::NLaser > 0) Time += dt;
    if (childfopenMPI(sdt_phys, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n",
            Time, Phys::energy[0], Phys::var[0],
            Phys::doublon[0], Phys::num[0], step_i);
    fclose(fp);

    if (childfopenMPI(sdt_norm, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf %.16lf %d\n", Time, global_norm[0], step_i);
    fclose(fp);

    if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", 
      Time, Phys::num[0], Phys::num2[0], Phys::doublon[0],
            Phys::doublon2[0], Phys::Sz[0], Phys::Sz2[0], step_i);
    fclose(fp);

    if (step_i % step_spin == 0) {
      expec_cisajs(1, v2, v1);
      expec_cisajscktaltdc(1, v2, v1);
    }
    if (Def::iOutputEigenVec == TRUE) {
      if (step_i % Param::OutputInterval == 0) {
        sprintf(sdt, "%s_eigenvec_%d_rank_%d.dat", Def::CDataFileHead, step_i, myrank);
        if (childfopenALL(sdt, "wb", &fp) != 0) {
          fclose(fp);
          exitMPI(-1);
        }
        fwrite(&step_i, sizeof(step_i), 1, fp);
        fwrite(&Check::idim_max, sizeof(long int), 1, fp);
        fwrite(&v1[0][0], sizeof(std::complex<double>), Check::idim_max + 1, fp);
        fclose(fp);
      }
    }
  }/*for (step_i = step_initial; step_i < Def::Lanczos_max; step_i++)*/
  free_cd_2d_allocate(v2);

  if (Def::iOutputEigenVec == TRUE) {
    sprintf(sdt, "%s_eigenvec_%d_rank_%d.dat", Def::CDataFileHead, rand_i, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      fclose(fp);
      exitMPI(-1);
    }
    fwrite(&step_i, sizeof(step_i), 1, fp);
    fwrite(&Check::idim_max, sizeof(long int), 1, fp);
    fwrite(&v1[0][0], sizeof(std::complex<double>), Check::idim_max + 1, fp);
    fclose(fp);
  }

  fprintf(stdoutMPI, "%s", "######  End  : TimeEvolution by Taylor Expansion.  ######\n\n");
  return 0;
}
