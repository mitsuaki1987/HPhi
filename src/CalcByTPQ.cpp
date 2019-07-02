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
#include "FirstMultiply.hpp"
#include "Multiply.hpp"
#include "expec_energy_flct.hpp"
#include "expec_cisajs.hpp"
#include "expec_cisajscktaltdc.hpp"
#include "CalcByTPQ.hpp"
#include "FileIO.hpp"
#include "wrapperMPI.hpp"
#include "CalcTime.hpp"
#include "common/setmemory.hpp"
#include "mltplyCommon.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include "log.hpp"
/**
 * @file
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @brief  File for givinvg functions of TPQ method
 */
/** 
 * @brief A main function to calculate physical quqntities by TPQ method
 *
 * @param [in] NumAve  Number of samples
 * @param [in] ExpecInterval interval steps between the steps to calculate physical quantities
 * @param [in,out] X CalcStruct list for getting and giving calculation information
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 */
int CalcByTPQ(
  const int NumAve,
  const int ExpecInterval
) {
  char sdt[D_FileNameMax];
  char **sdt_phys, **sdt_norm, **sdt_flct;
  int rand_i, iret;
  long int i_max;
  int step_iO = 0;
  FILE *fp;
  double *inv_temp, Ns;
  size_t byte_size;

  inv_temp = d_1d_allocate(NumAve);

  step_spin = ExpecInterval;
  Def::St = 0;
  fprintf(stdoutMPI, "%s", "######  Start: TPQCalculation.  ######\n\n");
  global_norm = d_1d_allocate(NumAve);
  global_1st_norm = d_1d_allocate(NumAve);

  //for rand_i =0, rand_i<NumAve
  sdt_phys = (char**)malloc(sizeof(char*)*NumAve);
  sdt_norm = (char**)malloc(sizeof(char*)*NumAve);
  sdt_flct = (char**)malloc(sizeof(char*)*NumAve);
  for (rand_i = 0; rand_i < NumAve; rand_i++) {
    sdt_phys[rand_i] = (char*)malloc(sizeof(char)*D_FileNameMax);
    sdt_norm[rand_i] = (char*)malloc(sizeof(char)*D_FileNameMax);
    sdt_flct[rand_i] = (char*)malloc(sizeof(char)*D_FileNameMax);
    sprintf(sdt_phys[rand_i], "SS_rand%d.dat", rand_i);
    sprintf(sdt_norm[rand_i], "Norm_rand%d.dat", rand_i);
    sprintf(sdt_flct[rand_i], "Flct_rand%d.dat", rand_i);
  }
  Ns = 1.0 * Def::NsiteMPI;
  fprintf(stdoutMPI, "  rand_i / rand_max  =  %d / %d\n", 1, NumAve);
  iret = 0;

  //Make or Read initial vector
  if (Def::iReStart == RESTART_INOUT || Def::iReStart == RESTART_IN) {
    StartTimer(3600);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector starts: %s\n", "a", 0, step_i);
    fprintf(stdoutMPI, "%s", "  Start:  Input vector.\n");
    sprintf(sdt, "tmpvec_set%d_rank_%d.dat", rand_i, myrank);
    childfopenALL(sdt, "rb", &fp);
    if (fp == NULL) {
      fprintf(stdout, "A file of Inputvector does not exist.\n");
      fprintf(stdout, "Start to calculate in normal procedure.\n");
      iret = 1;
    }
    byte_size = fread(&step_i, sizeof(step_i), 1, fp);
    byte_size = fread(&i_max, sizeof(long int), 1, fp);
    if (i_max != Check::idim_max) {
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      exitMPI(-1);
    }
    byte_size = fread(v0, sizeof(std::complex<double>), (Check::idim_max + 1)*NumAve, fp);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector finishes: %s\n", "a", 0, step_i);
    fprintf(stdoutMPI, "%s", "  End  :  Input vector.\n");
    fclose(fp);
    StopTimer(3600);
    Def::istep = step_i;
    StartTimer(3200);
    iret = expec_energy_flct(NumAve, v0, v1);
    StopTimer(3200);
    if (iret != 0) return -1;

    step_iO = step_i - 1;
    if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
  }/*if (Def::iReStart == RESTART_INOUT || Def::iReStart == RESTART_IN)*/

  if (Def::iReStart == RESTART_NOT || Def::iReStart == RESTART_OUT || iret == 1) {
    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      if (childfopenMPI(sdt_phys[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", " # inv_tmp, energy, phys_var, phys_doublon, phys_num, step_i\n");
        fclose(fp);
      }
      else return -1;
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", " # inv_temp, global_norm, global_1st_norm, step_i \n");
        fclose(fp); 
      }
      else return -1;
      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", " # inv_temp, N, N^2, D, D^2, Sz, Sz^2, step_i \n");
        fclose(fp);
      }
      else return -1;
    }
    StopTimer(3600);

    step_i = 0;

    StartTimer(3100);
    if (rand_i == 0) {
      TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "set %d step %d:TPQ begins: %s", "w", 0, step_i);
    }
    else {
      TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "set %d step %d:TPQ begins: %s", "a", 0, step_i);
    }
    /**@brief
    Initialize v1 and compute v0 = H*v1
    */
    FirstMultiply();
    StopTimer(3100);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = 0.0;
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::energy[rand_i], Phys::var[rand_i],
          Phys::doublon[rand_i], Phys::num[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], global_1st_norm[rand_i], global_1st_norm[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
    }
    /**@brief
    Compute v1=0, and compute v0 = H*v1
    */
    StartTimer(3200);
    iret = expec_energy_flct(NumAve, v0, v1); //v0 = H*v1
    StopTimer(3200);
    if (iret != 0) return -1;
    step_i += 1;
    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = (2.0 / Ns) / (LargeValue - Phys::energy[rand_i] / Ns);
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::energy[rand_i], Phys::var[rand_i],
          Phys::doublon[rand_i], Phys::num[rand_i], step_i);
        fclose(fp);
      }
      else return -1;     
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], global_norm[rand_i], global_1st_norm[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::num[rand_i], Phys::num2[rand_i], 
          Phys::doublon[rand_i], Phys::doublon2[rand_i], 
          Phys::Sz[rand_i], Phys::Sz2[rand_i], step_i);
        fclose(fp);
      }
      else return -1;      
    }/*for (rand_i = 0; rand_i < NumAve; rand_i++)*/
    StopTimer(3600);
    step_i += 1;
    Def::istep = step_i;
    step_iO = 0;
  }/*if (Def::iReStart == RESTART_NOT || Def::iReStart == RESTART_OUT || iret == 1)*/

  for (step_i = Def::istep; step_i < Def::Lanczos_max; step_i++) {
    Def::istep = step_i;
    if (step_i % ((Def::Lanczos_max - step_iO) / 10) == 0) {
      fprintf(stdoutMPI, "    step_i/total_step = %d/%d \n", step_i, Def::Lanczos_max);
    }
    Def::istep = step_i;
    StartTimer(3600);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "set %d step %d:TPQ begins: %s", "a", 0, step_i);
    StopTimer(3600);
    StartTimer(3500);
    Multiply();
    StopTimer(3500);

    if (step_i%step_spin == 0) {
      StartTimer(3300);
      iret = expec_cisajs(NumAve, v1, v0);
      StopTimer(3300);
      if (iret != 0) return -1;

      StartTimer(3400);
      iret = expec_cisajscktaltdc(NumAve, v1, v0);
      StopTimer(3400);
      if (iret != 0) return -1;
    }

    StartTimer(3200);
    iret = expec_energy_flct(NumAve, v0, v1);
    StopTimer(3200);
    if (iret != 0) return -1;

    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = (2.0*step_i / Ns) / (LargeValue - Phys::energy[rand_i] / Ns);
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::energy[rand_i], Phys::var[rand_i],
          Phys::doublon[rand_i], Phys::num[rand_i], step_i);
        // for
        fclose(fp);
      }
      else return FALSE;

      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], global_norm[rand_i], global_1st_norm[rand_i], step_i);
        fclose(fp);
      }
      else return FALSE;

      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::num[rand_i], Phys::num2[rand_i],
          Phys::doublon[rand_i], Phys::doublon2[rand_i],
          Phys::Sz[rand_i], Phys::Sz2[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
    }/*for (rand_i = 0; rand_i < NumAve; rand_i++)*/
    StopTimer(3600);
  }/*for (step_i = Def::istep; step_i < Def::Lanczos_max; step_i++)*/

  if (Def::iReStart == RESTART_OUT || Def::iReStart == RESTART_INOUT) {
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector starts: %s\n", "a", 0, step_i);
    fprintf(stdoutMPI, "%s", "  Start:  Output vector.\n");
    sprintf(sdt, "tmpvec_set%d_rank_%d.dat", 0, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      exitMPI(-1);
    }
    fwrite(&step_i, sizeof(step_i), 1, fp);
    fwrite(&Check::idim_max, sizeof(Check::idim_max), 1, fp);
    fwrite(v1, sizeof(std::complex<double>), (Check::idim_max + 1)*NumAve, fp);
    fclose(fp);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector finishes: %s\n", "a", 0, step_i);
    fprintf(stdoutMPI, "%s", "  End  :  Output vector.\n");
  }

  fprintf(stdoutMPI, "%s", "######  End  : TPQCalculation.  ######\n\n");


  free_d_1d_allocate(inv_temp);

  for (rand_i = 0; rand_i < NumAve; rand_i++) {
    free(sdt_phys[rand_i]);
    free(sdt_norm[rand_i]);
    free(sdt_flct[rand_i]);
  }
  free(sdt_phys);
  free(sdt_norm);
  free(sdt_flct);

  return TRUE;
}
