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
#include "CalcSpectrum.hpp"
#include "CalcSpectrumByBiCG.hpp"
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
  const int ExpecInterval,
  int NdcSpectrum
) {
  char sdt[D_FileNameMax];
  char **sdt_phys, **sdt_norm, **sdt_flct;
  int rand_i, iret, iomega;
  long int i_max;
  int step_iO = 0;
  FILE *fp;
  double *inv_temp, Ns;
  size_t byte_size;
  std::complex<double>*** dcSpectrum;
  std::complex<double>** dcomega;

  inv_temp = d_1d_allocate(NumAve);
  if (Def::iFlgCalcSpec != DC::CALCSPEC_NOT) {
    dcomega = cd_2d_allocate(Def::k_exct, Def::iNOmega);
    dcSpectrum = cd_3d_allocate(Def::k_exct, Def::iNOmega, NdcSpectrum);
  }/*if (Def::iFlgCalcSpec != DC::CALCSPEC_NOT)*/

  Step::step_spin = ExpecInterval;
  Def::St = 0;
  fprintf(MP::STDOUT, "%s", "######  Start: TPQCalculation.  ######\n\n");
  Step::global_norm = d_1d_allocate(NumAve);
  Step::global_1st_norm = d_1d_allocate(NumAve);

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
  fprintf(MP::STDOUT, "  rand_i / rand_max  =  %d / %d\n", 1, NumAve);
  iret = 0;

  //Make or Read initial vector
  if (Def::iReStart == DC::RESTART_INOUT || Def::iReStart == DC::RESTART_IN) {
    StartTimer(3600);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector starts: %s\n", "a", 0, Step::step_i);
    fprintf(MP::STDOUT, "%s", "  Start:  Input vector.\n");
    sprintf(sdt, "tmpvec_set%d_rank_%d.dat", rand_i, MP::myrank);
    childfopenALL(sdt, "rb", &fp);
    if (fp == NULL) {
      fprintf(stdout, "A file of Inputvector does not exist.\n");
      fprintf(stdout, "Start to calculate in normal procedure.\n");
      iret = 1;
    }
    byte_size = fread(&Step::step_i, sizeof(Step::step_i), 1, fp);
    byte_size = fread(&i_max, sizeof(long int), 1, fp);
    if (i_max != Check::idim_max) {
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      wrapperMPI::Exit(-1);
    }
    byte_size = fread(Wave::v0, sizeof(std::complex<double>), Check::idim_max*NumAve, fp);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector finishes: %s\n", "a", 0, Step::step_i);
    fprintf(MP::STDOUT, "%s", "  End  :  Input vector.\n");
    fclose(fp);
    StopTimer(3600);
    Def::istep = Step::step_i;
    StartTimer(3200);
    iret = expec::energy_flct::main(NumAve, Wave::v0, Wave::v1);
    StopTimer(3200);
    if (iret != 0) return -1;

    step_iO = Step::step_i - 1;
    if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
  }/*if (Def::iReStart == DC::RESTART_INOUT || Def::iReStart == DC::RESTART_IN)*/

  if (Def::iReStart == DC::RESTART_NOT || Def::iReStart == DC::RESTART_OUT || iret == 1) {
    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      if (childfopenMPI(sdt_phys[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", " # inv_tmp, energy, phys_var, phys_doublon, phys_num, Step::step_i\n");
        fclose(fp);
      }
      else return -1;
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", " # inv_temp, Step::global_norm, Step::global_1st_norm, Step::step_i \n");
        fclose(fp); 
      }
      else return -1;
      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", " # inv_temp, N, N^2, D, D^2, Sz, Sz^2, Step::step_i \n");
        fclose(fp);
      }
      else return -1;
    }
    StopTimer(3600);

    Step::step_i = 0;

    StartTimer(3100);
    if (rand_i == 0) {
      TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "set %d step %d:TPQ begins: %s", "w", 0, Step::step_i);
    }
    else {
      TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "set %d step %d:TPQ begins: %s", "a", 0, Step::step_i);
    }
    /**@brief
    Initialize Wave::v1 and compute v0 = H*Wave::v1
    */
    FirstMultiply();
    StopTimer(3100);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = 0.0;
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::energy[rand_i], Phys::var[rand_i],
          Phys::doublon[rand_i], Phys::num[rand_i], Step::step_i);
        fclose(fp);
      }
      else return -1;
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Step::global_1st_norm[rand_i], Step::global_1st_norm[rand_i], Step::step_i);
        fclose(fp);
      }
      else return -1;
    }
    /**@brief
    Compute Wave::v1=0, and compute v0 = H*Wave::v1
    */
    StartTimer(3200);
    iret = expec::energy_flct::main(NumAve, Wave::v0, Wave::v1); //v0 = H*Wave::v1
    StopTimer(3200);
    if (iret != 0) return -1;
    Step::step_i += 1;
    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = (2.0 / Ns) / (Step::LargeValue - Phys::energy[rand_i] / Ns);
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::energy[rand_i], Phys::var[rand_i],
          Phys::doublon[rand_i], Phys::num[rand_i], Step::step_i);
        fclose(fp);
      }
      else return -1;     
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Step::global_norm[rand_i], Step::global_1st_norm[rand_i], Step::step_i);
        fclose(fp);
      }
      else return -1;
      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::num[rand_i], Phys::num2[rand_i], 
          Phys::doublon[rand_i], Phys::doublon2[rand_i], 
          Phys::Sz[rand_i], Phys::Sz2[rand_i], Step::step_i);
        fclose(fp);
      }
      else return -1;      
    }/*for (rand_i = 0; rand_i < NumAve; rand_i++)*/
    StopTimer(3600);
    Step::step_i += 1;
    Def::istep = Step::step_i;
    step_iO = 0;
  }/*if (Def::iReStart == DC::RESTART_NOT || Def::iReStart == DC::RESTART_OUT || iret == 1)*/

  for (Step::step_i = Def::istep; Step::step_i < Def::Lanczos_max; Step::step_i++) {
    Def::istep = Step::step_i;
    if (Step::step_i % ((Def::Lanczos_max - step_iO) / 10) == 0) {
      fprintf(MP::STDOUT, "    Step::step_i/total_step = %d/%d \n", Step::step_i, Def::Lanczos_max);
    }
    Def::istep = Step::step_i;
    StartTimer(3600);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "set %d step %d:TPQ begins: %s", "a", 0, Step::step_i);
    StopTimer(3600);
    StartTimer(3500);
    Multiply();
    StopTimer(3500);

    if (Step::step_i%Step::step_spin == 0) {
      StartTimer(3300);
      iret = expec::cisajs::main(NumAve, Wave::v1, Wave::v0);
      StopTimer(3300);
      if (iret != 0) return -1;

      StartTimer(3400);
      iret = expec::cisajscktalt::main(NumAve, Wave::v1, Wave::v0);
      StopTimer(3400);
      if (iret != 0) return -1;
    }/*if (Step::step_i%Step::step_spin == 0)*/

    StartTimer(3200);
    iret = expec::energy_flct::main(NumAve, Wave::v0, Wave::v1);
    StopTimer(3200);
    if (iret != 0) return -1;
    /*
     For spectrum calculation with mTPQ+sBiCG
    */
    if (Def::iFlgCalcSpec != DC::CALCSPEC_NOT && Step::step_i % Step::step_spin == 0) {
      for (rand_i = 0; rand_i < NumAve; rand_i++) {
        for (iomega = 0; iomega < Def::iNOmega; iomega++) {
          dcomega[rand_i][iomega] = Phys::energy[rand_i] + Def::dcOmegaMin
            + (Def::dcOmegaMax - Def::dcOmegaMin) / (std::complex<double>)Def::iNOmega * (std::complex<double>)iomega;
        }
      }/*for (rand_i = 0; rand_i < NumAve; rand_i++)*/
      Spectrum::BiCG(NumAve, Def::iNOmega, NdcSpectrum, dcSpectrum, dcomega, Wave::v1);
      Spectrum::OutputSpectrum(NumAve, Def::iNOmega, NdcSpectrum, dcSpectrum, dcomega);
    }/*if (Def::iFlgCalcSpec != DC::CALCSPEC_NOT)*/
    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = (2.0*Step::step_i / Ns) / (Step::LargeValue - Phys::energy[rand_i] / Ns);
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::energy[rand_i], Phys::var[rand_i],
          Phys::doublon[rand_i], Phys::num[rand_i], Step::step_i);
        // for
        fclose(fp);
      }
      else return FALSE;

      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Step::global_norm[rand_i], Step::global_1st_norm[rand_i], Step::step_i);
        fclose(fp);
      }
      else return FALSE;

      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], Phys::num[rand_i], Phys::num2[rand_i],
          Phys::doublon[rand_i], Phys::doublon2[rand_i],
          Phys::Sz[rand_i], Phys::Sz2[rand_i], Step::step_i);
        fclose(fp);
      }
      else return -1;
    }/*for (rand_i = 0; rand_i < NumAve; rand_i++)*/
    StopTimer(3600);
  }/*for (Step::step_i = Def::istep; Step::step_i < Def::Lanczos_max; Step::step_i++)*/

  if (Def::iReStart == DC::RESTART_OUT || Def::iReStart == DC::RESTART_INOUT) {
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector starts: %s\n", "a", 0, Step::step_i);
    fprintf(MP::STDOUT, "%s", "  Start:  Output vector.\n");
    sprintf(sdt, "tmpvec_set%d_rank_%d.dat", 0, MP::myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      wrapperMPI::Exit(-1);
    }
    fwrite(&Step::step_i, sizeof(Step::step_i), 1, fp);
    fwrite(&Check::idim_max, sizeof(Check::idim_max), 1, fp);
    fwrite(Wave::v1, sizeof(std::complex<double>), Check::idim_max*NumAve, fp);
    fclose(fp);
    TimeKeeperWithRandAndStep("%s_Time_TPQ_Step.dat", "  set %d step %d:output vector finishes: %s\n", "a", 0, Step::step_i);
    fprintf(MP::STDOUT, "%s", "  End  :  Output vector.\n");
  }

  fprintf(MP::STDOUT, "%s", "######  End  : TPQCalculation.  ######\n\n");


  free_d_1d_allocate(inv_temp);

  for (rand_i = 0; rand_i < NumAve; rand_i++) {
    free(sdt_phys[rand_i]);
    free(sdt_norm[rand_i]);
    free(sdt_flct[rand_i]);
  }
  free(sdt_phys);
  free(sdt_norm);
  free(sdt_flct);

  if (Def::iFlgCalcSpec != DC::CALCSPEC_NOT) {
    free_cd_3d_allocate(dcSpectrum);
    free_cd_2d_allocate(dcomega);
  }/*if (Def::iFlgCalcSpec != DC::CALCSPEC_NOT)*/

  return TRUE;
}
