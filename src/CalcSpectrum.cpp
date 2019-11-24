/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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

/**
 * @file
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @brief  File for givinvg functions of calculating spectrum
 *
 *
 */

#include "mltply.hpp"
#include "CalcSpectrum.hpp"
#include "CalcSpectrumByBiCG.hpp"
#include "CalcSpectrumByFullDiag.hpp"
#include "CalcTime.hpp"
#include "SingleEx.hpp"
#include "PairEx.hpp"
#include "wrapperMPI.hpp"
#include "FileIO.hpp"
#include "common/setmemory.hpp"
#include "readdef.hpp"
#include "xsetmem.hpp"
#include "sz.hpp"
#include "check.hpp"
#include "diagonalcalc.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include "log.hpp"
#include <iostream>

/**
@brief Set target frequencies
Output: dcOmegaMax, dcOmegaMin
@retval FALSE Fail to set frequencies.
retval TRUE Success to set frequencies.
*/
int SetOmega()
{
  FILE *fp;
  char sdt[D_FileNameMax], ctmp[256];
  int istp = 4;
  double E1, E2, E3, E4, Emax;
  long int iline_countMax = 2;
  long int iline_count = 2;


  if (Def::iFlgSpecOmegaMax == TRUE && Def::iFlgSpecOmegaMin == TRUE) {
    return TRUE;
  }
  else {
    if (Def::iCalcType == DC::Lanczos || Def::iCalcType == DC::FullDiag) {
      sprintf(sdt, "%s_Lanczos_Step.dat", Def::CDataFileHead);
      childfopenMPI(sdt, "r", &fp);
      if (fp == NULL) {
        fprintf(MP::STDOUT, "Error: xx_Lanczos_Step.dat does not exist.\n");
        return FALSE;
      }
      wrapperMPI::Fgets(ctmp, 256, fp); //1st line is skipped
      wrapperMPI::Fgets(ctmp, 256, fp); //2nd line is skipped
      while (wrapperMPI::Fgets(ctmp, 256, fp) != NULL) {
        iline_count++;
      }
      iline_countMax = iline_count;
      iline_count = 2;
      rewind(fp);
      wrapperMPI::Fgets(ctmp, 256, fp); //1st line is skipped
      wrapperMPI::Fgets(ctmp, 256, fp); //2nd line is skipped

      while (wrapperMPI::Fgets(ctmp, 256, fp) != NULL) {
        sscanf(ctmp, "stp=%d %lf %lf %lf %lf %lf\n",
          &istp,
          &E1,
          &E2,
          &E3,
          &E4,
          &Emax);
        iline_count++;
        if (iline_count == iline_countMax) break;
      }
      fclose(fp);
      if (istp < 4) {
        fprintf(MP::STDOUT, "Error: Lanczos step must be greater than 4 for using spectrum calculation.\n");
        return FALSE;
      }
    }/*if (Def::iCalcType == DC::Lanczos || Def::iCalcType == DC::FullDiag)*/
    else
    {
      sprintf(sdt, "%s_energy.dat", Def::CDataFileHead);
      childfopenMPI(sdt, "r", &fp);
      if (fp == NULL) {
        fprintf(MP::STDOUT, "Error: xx_energy.dat does not exist.\n");
        return FALSE;
      }/*if (fp == NULL)*/
      wrapperMPI::Fgets(ctmp, 256, fp); //1st line is skipped
      wrapperMPI::Fgets(ctmp, 256, fp); //1st line is skipped
      sscanf(ctmp, "  Energy  %lf \n", &E1);
      Emax = Step::LargeValue;
    }/**/
    //Read Lanczos_Step
    if (Def::iFlgSpecOmegaMax == FALSE) {
      Def::dcOmegaMax = Emax * (double)Def::Nsite;
    }
    if (Def::iFlgSpecOmegaMin == FALSE) {
      Def::dcOmegaMin = E1;
    }
  }/*Omegamax and omegamin is not specified in modpara*/

  return TRUE;
}
/**
@brief Make the lists for the excited state; List::a1, List::a2_1 and List::a2_2 (for canonical ensemble).
The original lists before the excitation are given by List::cxxx_org
Output: iCalcModel (From HubbardNConserved to Hubbard), {Ne, Nup, Ndown, Nsite, Total2Sz} (update for MPI)
@param iFlgListModifed [out] If the list is modified due to the excitation, the value becomes TRUE(1), otherwise FALSE(0).
@retval -1 fail to make lists.
@retval 0  sucsess to make lists.
*/
int MakeExcitedList(
  int *iFlgListModifed
) {
  int Ne, Nup, Ndown, Total2Sz;
  long int j;
  *iFlgListModifed = FALSE;
  //To Get Original space

  if (Def::NNSingleExcitationOperator > 0) {
    switch (Def::iCalcModel) {
    case DC::HubbardGC:
      break;
    case DC::HubbardNConserved:
    case DC::KondoGC:
    case DC::Hubbard:
    case DC::Kondo:
      *iFlgListModifed = TRUE;
      break;
    case DC::Spin:
    case DC::SpinGC:
      return FALSE;
    }
  }
  else if (Def::NNPairExcitationOperator > 0) {
    switch (Def::iCalcModel) {
    case DC::HubbardGC:
    case DC::SpinGC:
    case DC::HubbardNConserved:
      break;
    case DC::KondoGC:
    case DC::Hubbard:
    case DC::Kondo:
    case DC::Spin:
      if (Def::PairExcitationOperator[0][0][1] != Def::PairExcitationOperator[0][0][3]) {
        *iFlgListModifed = TRUE;
      }
      break;
    }
  }
  else {
    return FALSE;
  }

  if (*iFlgListModifed == TRUE) {

    Ne = Def::NeMPI;
    Nup = Def::NupMPI;
    Ndown = Def::NdownMPI;
    Total2Sz = Def::Total2SzMPI;
    if (Def::NNSingleExcitationOperator > 0) {
      switch (Def::iCalcModel) {
      case DC::HubbardGC:
        break;
      case DC::HubbardNConserved:
        if (Def::SingleExcitationOperator[0][0][2] == 1) { //cis
          Ne = Def::NeMPI + 1;
        }
        else {
          Ne = Def::NeMPI - 1;
        }
        break;
      case DC::KondoGC:
      case DC::Hubbard:
      case DC::Kondo:
        if (Def::SingleExcitationOperator[0][0][2] == 1) { //cis
          Ne = Def::NeMPI + 1;
          if (Def::SingleExcitationOperator[0][0][1] == 0) {//up
            Nup = Def::NupMPI + 1;
            Ndown = Def::NdownMPI;
          }
          else {//down
            Nup = Def::NupMPI;
            Ndown = Def::NdownMPI + 1;
          }
        }
        else {//ajt
          Ne = Def::NeMPI - 1;
          if (Def::SingleExcitationOperator[0][0][1] == 0) {//up
            Nup = Def::NupMPI - 1;
            Ndown = Def::NdownMPI;

          }
          else {//down
            Nup = Def::NupMPI;
            Ndown = Def::NdownMPI - 1;
          }
        }
        break;
      case DC::Spin:
      case DC::SpinGC:
        return FALSE;
      }
    }
    else if (Def::NNPairExcitationOperator > 0) {
      Ne = Def::NeMPI;
      switch (Def::iCalcModel) {
      case DC::HubbardGC:
      case DC::SpinGC:
      case DC::HubbardNConserved:
        break;
      case DC::KondoGC:
      case DC::Hubbard:
      case DC::Kondo:
        if (Def::PairExcitationOperator[0][0][1] != Def::PairExcitationOperator[0][0][3]) {
          if (Def::PairExcitationOperator[0][0][1] == 0) {//up
            Nup = Def::NupMPI + 1;
            Ndown = Def::NdownMPI - 1;
          }
          else {//down
            Nup = Def::NupMPI - 1;
            Ndown = Def::NdownMPI + 1;
          }
        }
        break;
      case DC::Spin:
        if (Def::PairExcitationOperator[0][0][1] != Def::PairExcitationOperator[0][0][3]) {
          if (Def::iFlgGeneralSpin == FALSE) {
            if (Def::PairExcitationOperator[0][0][1] == 0) {//down
              Nup = Def::NupMPI - 1;
              Ndown = Def::NdownMPI + 1;
            }
            else {//up
              Nup = Def::NupMPI + 1;
              Ndown = Def::NdownMPI - 1;
            }
          }
          else {//for general spin
            Total2Sz = Def::Total2SzMPI + 2 * (Def::PairExcitationOperator[0][0][1] - Def::PairExcitationOperator[0][0][3]);
          }
        }
        break;
      }
    }
    else {
      return FALSE;
    }
    //Update Infomation
    Def::Nsite = Def::NsiteMPI;

    if (check(&Ne, &Nup, &Ndown, &Total2Sz, &Check::idim_maxs) == MPIFALSE) {
      wrapperMPI::Finalize();
      return FALSE;
    }
    if (GetlistSize() == TRUE) {
      List::b1 = li_1d_allocate(Check::idim_maxs);
      List::b2_1 = li_1d_allocate(Large::SizeOflist_2_1);
      List::b2_2 = li_1d_allocate(Large::SizeOflist_2_2);
      if (List::b1 == NULL
        || List::b2_1 == NULL
        || List::b2_2 == NULL
        )
      {
        return -1;
      }
      for (j = 0; j < Large::SizeOflist_2_1; j++) {
        List::b2_1[j] = 0;
      }
      for (j = 0; j < Large::SizeOflist_2_2; j++) {
        List::b2_2[j] = 0;
      }
    }/*if (GetlistSize() == TRUE)*/
    List::Diagonals = d_1d_allocate(Check::idim_maxs);
    if (sz(List::b1, List::b2_1, List::b2_2, Ne, Nup, Ndown, Total2Sz, Check::idim_maxs) != 0) {
      return FALSE;
    }
    /*
     MPI buffer is common for list_a and list_b
    */
#ifdef __MPI
    long int MAXidim_max, MAXidim_maxs;
    MAXidim_max = wrapperMPI::Max_li(Check::idim_max);
    MAXidim_maxs = wrapperMPI::Max_li(Check::idim_maxs);
    if (MAXidim_max < MAXidim_maxs) {
      free_cd_2d_allocate(Wave::v1buf);
      free_li_1d_allocate(List::buf);
      Wave::v1buf = cd_2d_allocate(MAXidim_maxs, Def::k_exct);
      List::buf = li_1d_allocate(MAXidim_maxs);
    }
#endif // MPI
  }
  else {
    Check::idim_maxs = Check::idim_max;
    List::b1 = List::a1;
    List::b2_1 = List::a2_1;
    List::b2_2 = List::a2_2;
    List::Diagonals = List::Diagonal;
  }
  Wave::v0 = cd_2d_allocate(Check::idim_maxs, Def::k_exct);
  Wave::v1 = cd_2d_allocate(Check::idim_maxs, Def::k_exct);

  if (Def::iCalcModel == DC::HubbardNConserved) {
    Def::iCalcModel = DC::Hubbard;
  }

  return TRUE;
}
/**
\brief Output spectrum.
\param Nomega [in] A total number of discrete frequencies.
\param dcSpectrum [in] Array of spectrum.
\param dcomega [in] Array of discrete frequencies.
\retval FALSE Fail to open the output file.
\retval TRUE Success to output the spectrum.
*/
void OutputSpectrum(
  int nstate,
  int Nomega,
  int NdcSpectrum,
  std::complex<double> ***dcSpectrum,
  std::complex<double> **dcomega)
{
  FILE *fp;
  char sdt[D_FileNameMax];
  int iomega, idcSpectrum, istate;

  for (istate = 0; istate < nstate; istate++) {
    sprintf(sdt, "%s_DynamicalGreen.dat", Def::CDataFileHead);
    if (childfopenMPI(sdt, "w", &fp) != 0) wrapperMPI::Exit(-1);

    for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
      for (iomega = 0; iomega < Nomega; iomega++) {
        fprintf(fp, "%.10lf %.10lf %.10lf %.10lf \n",
          real(dcomega[istate][iomega] - Def::dcOmegaOrg - Phys::energy[istate]),
          imag(dcomega[istate][iomega] - Def::dcOmegaOrg - Phys::energy[istate]),
          real(dcSpectrum[istate][iomega][idcSpectrum]),
          imag(dcSpectrum[istate][iomega][idcSpectrum]));
      }/*for (i = 0; i < Nomega; i++)*/
      fprintf(fp, "\n");
    }
    fclose(fp);
  }/*for (istate = 0; istate < nstate; istate++)*/
}/*int OutputSpectrum*/
/**
\brief Parent function to calculate the excited state.
\param tmp_v0 [out] Result @f$ v_0 = H_{ex} v_1 @f$.
\param tmp_v1 [in] The original state before excitation  @f$ v_1 @f$.
\retval FALSE Fail to calculate the excited state.
\retval TRUE Success to calculate the excited state.
*/
int GetExcitedState
(
  
  int nstate,
  std::complex<double> **tmp_v0,
  std::complex<double> **tmp_v1,
  int iEx
)
{
  if (Def::NNSingleExcitationOperator > 0) {
    if (GetSingleExcitedState(nstate, tmp_v0, tmp_v1, iEx) != TRUE) {
      return FALSE;
    }
  }
  else if (Def::NNPairExcitationOperator > 0) {
    if (GetPairExcitedState(nstate, tmp_v0, tmp_v1, iEx) != TRUE) {
      return FALSE;
    }
  }
  return TRUE;
}
/**
 * @brief A main function to calculate spectrum.
 *
 * input: iFlgSpecOmegaOrg, dcOmegaMax, dcOmegaMin, iNOmega etc.\n
 * output: dcOmegaOrg, iFlagListModified.
 *
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 *
 */
void CalcSpectrum()
{
  char sdt[D_FileNameMax];
  char *defname;
  long int i;
  long int i_max = 0;
  int i_stp, NdcSpectrum, istate;
  int iFlagListModified = FALSE;
  FILE *fp;
  std::complex<double> **v1Org; /**< Input vector to calculate spectrum function.*/

  //ToDo: Nomega should be given as a parameter
  int Nomega;
  std::complex<double> OmegaMax, OmegaMin;
  std::complex<double> ***dcSpectrum;
  std::complex<double> **dcomega;
  size_t byte_size;

  v1Org = cd_2d_allocate(Check::idim_max, Def::k_exct);
  for (i = 0; i < Check::idim_max; i++) 
    for (istate = 0; istate < Def::k_exct; istate++)
      v1Org[i][istate] = Wave::v1[i][istate];
  free_cd_2d_allocate(Wave::v0);
  free_cd_2d_allocate(Wave::v1);

  //set omega
  if (SetOmega() != TRUE) {
    fprintf(stderr, "Error: Fail to set Omega.\n");
    wrapperMPI::Exit(-1);
  }
  else {
    if (Def::iFlgSpecOmegaOrg == FALSE) {
      Def::dcOmegaOrg = std::complex<double>(0.0,1.0) * (Def::dcOmegaMax - Def::dcOmegaMin) / (double)Def::iNOmega;
    }
  }
  /*
   Set & malloc omega grid
  */
  Nomega = Def::iNOmega;
  dcomega = cd_2d_allocate(Def::k_exct, Nomega);
  OmegaMax = Def::dcOmegaMax + Def::dcOmegaOrg;
  OmegaMin = Def::dcOmegaMin + Def::dcOmegaOrg;
  for (istate = 0; istate < Def::k_exct; istate++) {
    for (i = 0; i < Nomega; i++) {
      dcomega[istate][i] = Phys::energy[istate] + OmegaMin
        + (OmegaMax - OmegaMin) / (std::complex<double>)Nomega * (std::complex<double>)i;
    }
  }

  fprintf(MP::STDOUT, "\nFrequency range:\n");
  fprintf(MP::STDOUT, "  Omega Max. : %15.5e %15.5e\n", real(OmegaMax), imag(OmegaMax));
  fprintf(MP::STDOUT, "  Omega Min. : %15.5e %15.5e\n", real(OmegaMin), imag(OmegaMin));
  fprintf(MP::STDOUT, "  Num. of Omega : %d\n", Nomega);

  if (Def::NNSingleExcitationOperator == 0) {
    if (Def::NNPairExcitationOperator == 0) {
      fprintf(stderr, "Error: Any excitation operators are not defined.\n");
      wrapperMPI::Exit(-1);
    }
    else {
      NdcSpectrum = Def::NNPairExcitationOperator - 1;
    }
  }
  else if (Def::NNPairExcitationOperator == 0) {
    NdcSpectrum = Def::NNSingleExcitationOperator - 1;
  }
  else {
    fprintf(stderr, "Error: Both single and pair excitation operators exist.\n");
    wrapperMPI::Exit(-1);
  }
  dcSpectrum = cd_3d_allocate(Def::k_exct, Nomega, NdcSpectrum);

  //Make New Lists
  if (MakeExcitedList(&iFlagListModified) == FALSE) {
    wrapperMPI::Exit(-1);
  }
  Def::iFlagListModified = iFlagListModified;
  Wave::v0 = cd_2d_allocate(Check::idim_maxs, Def::k_exct);
  Wave::v1 = cd_2d_allocate(Check::idim_maxs, Def::k_exct);

  //Make excited state
  StartTimer(6100);
  if (Def::iFlgCalcSpec == DC::RECALC_NOT ||
    Def::iFlgCalcSpec == DC::RECALC_OUTPUT_TMComponents_VEC ||
    (Def::iFlgCalcSpec == DC::RECALC_INOUT_TMComponents_VEC && Def::iCalcType == DC::CG)) {
    v1Org = cd_2d_allocate(Check::idim_max, 1);
    //input eigen vector
    StartTimer(6101);
    fprintf(MP::STDOUT, "  Start: An Eigenvector is inputted in CalcSpectrum.\n");
    TimeKeeper("%s_TimeKeeper.dat", "Reading an input Eigenvector starts: %s", "a");
    GetFileNameByKW(KWSpectrumVec, &defname);
    strcat(defname, "_rank_%d.dat");
    //    sprintf(sdt, "%s_eigenvec_%d_rank_%d.dat", Def::CDataFileHead, Def::k_exct - 1, MP::myrank);
    sprintf(sdt, defname, MP::myrank);
    childfopenALL(sdt, "rb", &fp);

    if (fp == NULL) {
      fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
      wrapperMPI::Exit(-1);
    }

    byte_size = fread(&i_stp, sizeof(i_stp), 1, fp);
    Large::itr = i_stp; //For TPQ
    byte_size = fread(&i_max, sizeof(i_max), 1, fp);
    if (i_max != Check::idim_max) {
      fprintf(stderr, "Error: MP::myrank=%d, i_max=%ld\n", MP::myrank, i_max);
      fprintf(stderr, "Error: A file of Input vector is incorrect.\n");
      wrapperMPI::Exit(-1);
    }
    byte_size = fread(&v1Org[0][0], sizeof(std::complex<double>), i_max + 1, fp);
    fclose(fp);
    StopTimer(6101);
    if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
  }
  StopTimer(6100);

  diagonalcalc(Check::idim_maxs, List::Diagonals, List::b1);

  fprintf(MP::STDOUT, "  Start: Calculating a spectrum.\n\n");
  TimeKeeper("%s_TimeKeeper.dat", "Calculating a spectrum starts: %s", "a");
  StartTimer(6200);
  switch (Def::iCalcType) {
  case DC::CG:
    CalcSpectrumByBiCG(Def::k_exct, Wave::v0, Wave::v1, Nomega, NdcSpectrum, dcSpectrum, dcomega, v1Org);
    OutputSpectrum(Def::k_exct, Nomega, NdcSpectrum, dcSpectrum, dcomega);
    break;
  case DC::FullDiag:
    CalcSpectrumByFullDiag(Nomega, NdcSpectrum, dcSpectrum, dcomega, v1Org);
    OutputSpectrum(Check::idim_max, Nomega, NdcSpectrum, dcSpectrum, dcomega);
    break;
  default:
    break;
  }
  StopTimer(6200);

  fprintf(MP::STDOUT, "  End:  Calculating a spectrum.\n\n");
  TimeKeeper("%s_TimeKeeper.dat", "Calculating a spectrum finishes: %s", "a");
  free_cd_3d_allocate(dcSpectrum);
  free_cd_2d_allocate(dcomega);
}/*int CalcSpectrum*/
