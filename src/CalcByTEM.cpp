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
#include "Common.hpp"
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
#include <iostream>
/**
 * @file   CalcByTEM.c
 *
 * @brief  File to define functions to calculate expected values by Time evolution method.
 */
 /// \brief Set transfer integrals at timeidx-th time
 /// \param X struct for getting information of transfer integrals
 /// \param timeidx index of time
void MakeTEDTransfer(struct BindStruct *X, const int timeidx) {
  int i, j;
  //Clear values
  for (i = 0; i < X->Def.NTETransferMax; i++) {
    for (j = 0; j < 4; j++) {
      X->Def.EDGeneralTransfer[i + X->Def.EDNTransfer][j] = 0;
    }
    X->Def.EDParaGeneralTransfer[i + X->Def.EDNTransfer] = 0.0;
  }

  //Input values
  for (i = 0; i < X->Def.NTETransfer[timeidx]; i++) {
    for (j = 0; j < 4; j++) {
      X->Def.EDGeneralTransfer[i + X->Def.EDNTransfer][j] = X->Def.TETransfer[timeidx][i][j];
    }
    X->Def.EDParaGeneralTransfer[i + X->Def.EDNTransfer] = X->Def.ParaTETransfer[timeidx][i];
  }
  X->Def.EDNTransfer += X->Def.NTETransfer[timeidx];
}
/// \brief Set interall interactions at timeidx-th time
/// \param X struct for getting information of interall interactions
/// \param timeidx index of time
void MakeTEDInterAll(struct BindStruct *X, const int timeidx) {
  int i, j;
  //Clear values
  for (i = 0; i < X->Def.NTEInterAllMax; i++) {
    for (j = 0; j < 8; j++) {
      X->Def.InterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal][j] = 0;
    }
    X->Def.ParaInterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal] = 0.0;
  }
  //Input values
  for (i = 0; i < X->Def.NTEInterAllOffDiagonal[timeidx]; i++) {
    for (j = 0; j < 8; j++) {
      X->Def.InterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal][j] = X->Def.TEInterAllOffDiagonal[timeidx][i][j];
    }
    X->Def.ParaInterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal] = X->Def.ParaTEInterAllOffDiagonal[timeidx][i];
  }
  X->Def.NInterAll_OffDiagonal += X->Def.NTEInterAllOffDiagonal[timeidx];
}
/** 
 * @brief main function of time evolution calculation
 * 
 * @param ExpecInterval interval to output expected values
 * @param X struct to get information of calculations.
 * @return 0 normally finished
 * @return -1 unnormally finished
 *
 * @author Kota Ido (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CalcByTEM(
  const int ExpecInterval,
  struct EDMainCalStruct *X
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
  double Time = X->Bind.Def.Param.Tinit;
  double dt = ((X->Bind.Def.NLaser == 0) ? 0.0 : X->Bind.Def.Param.TimeSlice);
  std::complex<double> **v2;  /**< Ttemporary vector for time evolution calculation, @f$ v2 = H*v1 = H^coef |psi(t)>@f$.*/

  global_norm = d_1d_allocate(1);

  if (X->Bind.Def.NTETimeSteps < X->Bind.Def.Lanczos_max) {
    fprintf(stdoutMPI, "Error: NTETimeSteps must be larger than Lanczos_max.\n");
    return -1;
  }
  step_spin = ExpecInterval;
  X->Bind.Def.St = 0;
  fprintf(stdoutMPI, "%s", "######  Start: TimeEvolution by Taylor Expansion.  ######\n\n");
  if (X->Bind.Def.iInputEigenVec == FALSE) {
    fprintf(stderr, "Error: A file of Inputvector is not inputted.\n");
    return -1;
  }
  else {
    //input v1
    fprintf(stdoutMPI, "%s", "An Initial Vector is inputted.\n");
    TimeKeeper(&(X->Bind), "%s_TimeKeeper.dat", "Reading an input Eigenvector starts: %s", "a");
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
    if (i_max != X->Bind.Check.idim_max) {
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      fclose(fp);
      printf("byte_size : %d\n", (int)byte_size);
      exitMPI(-1);
    }
    fread(&v1[0][0], sizeof(std::complex<double>), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);
    if (X->Bind.Def.iReStart == RESTART_NOT || X->Bind.Def.iReStart == RESTART_OUT) {
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


  int iInterAllOffDiagonal_org = X->Bind.Def.NInterAll_OffDiagonal;
  int iTransfer_org = X->Bind.Def.EDNTransfer;
  v2 = cd_2d_allocate(X->Bind.Check.idim_max + 1, 1);
  for (step_i = step_initial; step_i < X->Bind.Def.Lanczos_max; step_i++) {
    X->Bind.Def.istep = step_i;

    //Reset total number of interactions (changed in MakeTED***function.)
    X->Bind.Def.EDNTransfer = iTransfer_org;
    X->Bind.Def.NInterAll_OffDiagonal = iInterAllOffDiagonal_org;

    if (step_i % (X->Bind.Def.Lanczos_max / 10) == 0) {
      fprintf(stdoutMPI, "    step_i/total_step = %d/%d \n", step_i, X->Bind.Def.Lanczos_max);
    }

    if (X->Bind.Def.NLaser != 0) {
      TransferWithPeierls(&(X->Bind), Time);
    }
    else {
      // common procedure
      Time = X->Bind.Def.TETime[step_i];
      if (step_i == 0) dt = 0.0;
      else {
        dt = X->Bind.Def.TETime[step_i] - X->Bind.Def.TETime[step_i - 1];
      }
      X->Bind.Def.Param.TimeSlice = dt;

      // Set interactions
      if (X->Bind.Def.NTETransferMax != 0 && X->Bind.Def.NTEInterAllMax != 0) {
        fprintf(stdoutMPI,
          "Error: Time Evoluation mode does not support TEOneBody and TETwoBody interactions at the same time. \n");
        return -1;
      }
      else if (X->Bind.Def.NTETransferMax > 0) { //One-Body type
        MakeTEDTransfer(&(X->Bind), step_i);
      }
      else if (X->Bind.Def.NTEInterAllMax > 0) { //Two-Body type
        MakeTEDInterAll(&(X->Bind), step_i);
      }
      //[e] Yoshimi
    }

    if (step_i == step_initial) {
      TimeKeeperWithStep(&(X->Bind), "%s_Time_TE_Step.dat", "step %d:TE begins: %s", "w", step_i);
    }
    else {
      TimeKeeperWithStep(&(X->Bind), "%s_Time_TE_Step.dat", "step %d:TE begins: %s", "a", step_i);
    }
    MultiplyForTEM(&(X->Bind), v2);
    //Add Diagonal Parts
    //Multiply Diagonal
    expec_energy_flct(&(X->Bind), 1, v0, v1);

    if (X->Bind.Def.NLaser > 0) Time += dt;
    if (childfopenMPI(sdt_phys, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n",
            Time, X->Bind.Phys.energy[0], X->Bind.Phys.var[0],
            X->Bind.Phys.doublon[0], X->Bind.Phys.num[0], step_i);
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
      Time, X->Bind.Phys.num[0], X->Bind.Phys.num2[0], X->Bind.Phys.doublon[0],
            X->Bind.Phys.doublon2[0], X->Bind.Phys.Sz[0], X->Bind.Phys.Sz2[0], step_i);
    fclose(fp);

    if (step_i % step_spin == 0) {
      expec_cisajs(&(X->Bind), 1, v2, v1);
      expec_cisajscktaltdc(&(X->Bind), 1, v2, v1);
    }
    if (X->Bind.Def.iOutputEigenVec == TRUE) {
      if (step_i % X->Bind.Def.Param.OutputInterval == 0) {
        sprintf(sdt, "%s_eigenvec_%d_rank_%d.dat", X->Bind.Def.CDataFileHead, step_i, myrank);
        if (childfopenALL(sdt, "wb", &fp) != 0) {
          fclose(fp);
          exitMPI(-1);
        }
        fwrite(&step_i, sizeof(step_i), 1, fp);
        fwrite(&X->Bind.Check.idim_max, sizeof(long int), 1, fp);
        fwrite(&v1[0][0], sizeof(std::complex<double>), X->Bind.Check.idim_max + 1, fp);
        fclose(fp);
      }
    }
  }/*for (step_i = step_initial; step_i < X->Bind.Def.Lanczos_max; step_i++)*/
  free_cd_2d_allocate(v2);

  if (X->Bind.Def.iOutputEigenVec == TRUE) {
    sprintf(sdt, "%s_eigenvec_%d_rank_%d.dat", X->Bind.Def.CDataFileHead, rand_i, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      fclose(fp);
      exitMPI(-1);
    }
    fwrite(&step_i, sizeof(step_i), 1, fp);
    fwrite(&X->Bind.Check.idim_max, sizeof(long int), 1, fp);
    fwrite(&v1[0][0], sizeof(std::complex<double>), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);
  }

  fprintf(stdoutMPI, "%s", "######  End  : TimeEvolution by Taylor Expansion.  ######\n\n");
  return 0;
}
