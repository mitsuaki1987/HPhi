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
/*-------------------------------------------------------------
 *[ver.2008.11.4]
 *  Read Definition files
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 *-------------------------------------------------------------*/
/**
 * @file
 * 
 * @brief  File to define functions of reading input files
 */
#include "readdef.hpp"
#include <cctype>
#include "wrapperMPI.hpp"
#include "common/setmemory.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include <iostream>
/**
 @brief Keyword List in NameListFile.
 */
static char cKWListOfFileNameList[][D_CharTmpReadDef]={
  "CalcMod",
  "ModPara",
  "LocSpin",
  "Trans",
  "CoulombIntra",
  "CoulombInter",
  "Hund",
  "PairHop",
  "Exchange",
  "InterAll",
  "OneBodyG",
  "TwoBodyG",
  "PairLift",
  "Ising",
  "Boost",
  "SingleExcitation",
  "PairExcitation",
  "SpectrumVec",
  "Laser",
  "TEOneBody",
  "TETwoBody"
};

int D_iKWNumDef = sizeof(cKWListOfFileNameList)/sizeof(cKWListOfFileNameList[0]);
double eps_CheckImag0;/**> epsilon for checking values of one-body and two-body interactions.*/

/**
 * File Name List in NameListFile.
 **/
static char **cFileNameListFile;

int CheckInterAllCondition(
        int iCalcModel,
        int Nsite,
        int iFlgGeneralSpin,
        int *iLocSpin,
        int isite1, int isigma1,
        int isite2, int isigma2,
        int isite3, int isigma3,
        int isite4, int isigma4
);

int InputInterAllInfo(
        int *icnt_interall,
        int **iInterAllInfo,
        std::complex<double> *cInterAllValue,
        int isite1, int isigma1,
        int isite2, int isigma2,
        int isite3, int isigma3,
        int isite4, int isigma4,
        double re_value, double im_value
);


/**
 * @brief Error Function of reading def files.
 * @param[in] defname name of def file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileError(
                     const char *defname
                     ){
  fprintf(MP::STDOUT, "Error: %s (Broken file or Not exist)\n", defname);
  return (-1);
}

/**
 * @brief Function of Validating value.
 * @param[in] icheckValue value to validate.
 * @param[in] ilowestValue lowest value which icheckValue can be set.
 * @param[in] iHighestValue heighest value which icheckValue can be set.
 * @retval 0 value is correct.
 * @retval -1 value is incorrect.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int ValidateValue(
                  const int icheckValue,
                  const int ilowestValue, 
                  const int iHighestValue
                  ){

  if(icheckValue < ilowestValue || icheckValue > iHighestValue){
    return(-1);
  }
  return 0;
}

/**
 * @brief Function of Checking keyword in NameList file.
 * @param[in] cKW keyword candidate
 * @param[in] cKWList Reffercnce of keyword List
 * @param[in] iSizeOfKWidx number of keyword
 * @param[out] iKWidx index of keyword
 * @retval 0 keyword is correct.
 * @retval -1 keyword is incorrect.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckKW(
            const char* cKW,
            char  cKWList[][D_CharTmpReadDef],
            int iSizeOfKWidx,
            int* iKWidx
            ){
  *iKWidx=-1;
  int itmpKWidx;
  int iret=-1;
  for(itmpKWidx=0; itmpKWidx<iSizeOfKWidx; itmpKWidx++){
    if(strcmp(cKW,"")==0){
      break;
    }
    else if(CheckWords(cKW, cKWList[itmpKWidx])==0){
      iret=0;
      *iKWidx=itmpKWidx;
    }
  }
  return iret;
}

/**
 * @brief Function of Getting keyword and it's variable from characters.
 * @param[in] ctmpLine characters including keyword and it's variable 
 * @param[out] ctmp keyword
 * @param[out] itmp variable for a keyword
 * @retval 0 keyword and it's variable are obtained.
 * @retval 1 ctmpLine is a comment line.
 * @retval -1 format of ctmpLine is incorrect.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int GetKWWithIdx(
                 char *ctmpLine,
                 char *ctmp,
                 int *itmp
                 )
{
  char *ctmpRead;
  char *cerror;
  char csplit[] = " ,.\t\n";
  if(*ctmpLine=='\n') return 1;
  ctmpRead = strtok(ctmpLine, csplit);
  if(strncmp(ctmpRead, "=", 1)==0 || strncmp(ctmpRead, "#", 1)==0 || ctmpRead==NULL){
    return 1;
  }
  strcpy(ctmp, ctmpRead);
    
  ctmpRead = strtok( NULL, csplit );
  *itmp = strtol(ctmpRead, &cerror, 0);
  //if ctmpRead is not integer type
  if(*cerror != '\0'){
    fprintf(MP::STDOUT, "Error: incorrect format= %s. \n", cerror);
    return(-1);
  }

  ctmpRead = strtok( NULL, csplit );
  if(ctmpRead != NULL){
    fprintf(MP::STDOUT, "Error: incorrect format= %s. \n", ctmpRead);
    return(-1);
  }
    
  return 0;
}

/**
 * @brief Function of Reading calcmod file.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int ReadcalcmodFile(
  const char* defname/**<[in] file name to read.*/
)
{
  FILE *fp;
  int itmp, iret;
  char ctmpLine[D_CharTmpReadDef+D_CharKWDMAX];
  char ctmp[D_CharKWDMAX];
  Def::iCalcType=0;
  Def::iFlgFiniteTemperature=0;
  Def::iCalcModel=0;
  Def::iOutputMode=0;
  Def::iCalcEigenVec=0;
  Def::iInitialVecType=0;
  Def::iOutputEigenVec=0;
  Def::iInputEigenVec=0;
  Def::iOutputHam=0;
  Def::iInputHam=0;
  Def::iOutputExVec = 0;
  Def::iFlgCalcSpec=0;
  Def::iReStart=0;
  Def::iFlgMPI=0;
  Def::iFlgScaLAPACK=0;
#ifdef _MAGMA
  Def::iNGPU=2;
#else
  Def::iNGPU=0;
#endif
  /*=======================================================================*/
  fp = wrapperMPI::Fopen(defname, "r");
  if(fp==NULL) return ReadDefFileError(defname);
  /* read Parameters from calcmod.def*/
  while (wrapperMPI::Fgets(ctmpLine, D_CharTmpReadDef + D_CharKWDMAX, fp) != NULL) {
    if( (iret=GetKWWithIdx(ctmpLine, ctmp, &itmp)) !=0){
      if(iret==1) continue;
      return(-1);
    }   
    if(CheckWords(ctmp, "CalcType")==0){
      Def::iCalcType=itmp;
      if (Def::iCalcType == DC::Lanczos) {
        fprintf(MP::STDOUT, "  LOBPCG is used alternative to Lanczos.\n");
        Def::iCalcType = DC::CG;
      }
    }
    else if(CheckWords(ctmp, "FlgFiniteTemperature")==0){
      Def::iFlgFiniteTemperature = itmp;
    }
    else if(CheckWords(ctmp, "CalcModel")==0){
      Def::iCalcModel=itmp;
    }
    else if(CheckWords(ctmp, "OutputMode")==0){
      Def::iOutputMode=itmp;
    }
    else if(CheckWords(ctmp, "CalcEigenVec")==0){
      Def::iCalcEigenVec=itmp;
    }
    else if(CheckWords(ctmp, "InitialVecType")==0){
      Def::iInitialVecType=itmp;
    }
    else if(CheckWords(ctmp, "OutputEigenVec")==0 || CheckWords(ctmp, "OEV")==0){
      Def::iOutputEigenVec=itmp;
    }
    else if(CheckWords(ctmp, "InputEigenVec")==0 || CheckWords(ctmp, "IEV")==0){
      Def::iInputEigenVec=itmp;
    }
    else if(CheckWords(ctmp, "OutputHam")==0){
      Def::iOutputHam=itmp;
    }
    else if(CheckWords(ctmp, "InputHam")==0){
      Def::iInputHam=itmp;
    }
    else if(CheckWords(ctmp, "OutputExcitedVec")==0|| CheckWords(ctmp, "OutputExVec")==0){
      Def::iOutputExVec=itmp;
    }
    else if(CheckWords(ctmp, "CalcSpec")==0 || CheckWords(ctmp, "CalcSpectrum")==0){
      Def::iFlgCalcSpec=itmp;
    }
    else if(CheckWords(ctmp, "ReStart")==0){
      Def::iReStart=itmp;
    }
    else if(CheckWords(ctmp, "NGPU")==0){
        Def::iNGPU=itmp;
    }
    else if(CheckWords(ctmp, "ScaLAPACK")==0){
#ifdef _SCALAPACK
      Def::iFlgScaLAPACK=itmp;
#endif
    }
    else{
      fprintf(MP::STDOUT, "Error: In %s, wrong parameter name:%s \n", defname, ctmp);
      return(-1);
    }
  }
  fclose(fp);
    
  /* Check values*/
  if(ValidateValue(Def::iCalcModel, 0, DC::NUM_CALCMODEL-1)){
    fprintf(MP::STDOUT, "Error in %s\n CalcType: 0: Lanczos Method, 1: Thermal Pure Quantum State Method, 2: Full Diagonalization Method, 3: Calculation Spectrum mode.\n", defname);
    return (-1);
  }
  if(ValidateValue(Def::iCalcType, 0, DC::NUM_CALCTYPE-1)){
    fprintf(MP::STDOUT, "Error in %s\n CalcType: 0: Lanczos Method, 1: Thermal Pure Quantum State Method, 2: Full Diagonalization Method, 3: Calculation Spectrum mode.\n", defname);
    return (-1);
  }
  if(ValidateValue(Def::iOutputMode, 0, DC::NUM_OUTPUTMODE-1)){
    fprintf(MP::STDOUT, "Error in %s\n OutputMode: \n 0: calc one body green function and two body green functions,\n 1: calc one body green function and two body green functions and correlatinos for charge and spin.\n", defname);
    return (-1);
  }
  
  if(ValidateValue(Def::iCalcEigenVec, -1, DC::NUM_CALCEIGENVEC-1)){
    fprintf(MP::STDOUT, "Error in %s\n CalcEigenVec: \n 0: Lanczos+CG method,\n 1: Lanczos method.\n", defname);
    return (-1);
  }
  
  if(ValidateValue(Def::iInitialVecType, 0, DC::NUM_SETINITAILVEC-1)){
    fprintf(MP::STDOUT, "Error in %s\n InitialVecType: \n 0: complex type,\n 1: real type.\n", defname);
    return (-1);
  }

  if(ValidateValue(Def::iOutputHam, 0, DC::NUM_OUTPUTHAM-1)){
    fprintf(MP::STDOUT, "Error in %s\n OutputHam: \n 0: not output Hamiltonian,\n 1: output Hamiltonian.\n", defname);
    return (-1);
  }
  if(ValidateValue(Def::iInputHam, 0, DC::NUM_INPUTHAM-1)){
    fprintf(MP::STDOUT, "Error in %s\n InputHam: 0: not input Hamiltonian,\n 1: input Hamiltonian.\n", defname);
    return (-1);
  }
  if(Def::iInputHam == 1 && Def::iOutputHam==1){
    fprintf(MP::STDOUT,
            "Error in %s\n OutputHam=1 and InputHam=1.\n", defname);
    return (-1);
  }
  if(ValidateValue(Def::iReStart, 0, DC::NUM_RESTART-1)){
    fprintf(MP::STDOUT, "Error in %s Restart: \n 0: not restart (default).\n 1: output a restart vector.\n 2: input a restart vector and output a new restart vector.\n 3: input a restart vector.\n", defname);
    return (-1);
  }
  if(Def::iNGPU < 0){
    fprintf(MP::STDOUT, "Error in %s\n NGPU: NGPU must be greater than 0.\n", defname);
    return (-1);
  }
  if(ValidateValue(Def::iFlgScaLAPACK, 0, 1)){
    fprintf(MP::STDOUT, "Error in %s\n NGPU: NGPU must be greater than 0.\n", defname);
    return (-1);
  }

  /* In the case of Full Diagonalization method(iCalcType=2)*/
  if(Def::iCalcType==2 && ValidateValue(Def::iFlgFiniteTemperature, 0, 1)){
    fprintf(MP::STDOUT,
            "Error in %s\n FlgFiniteTemperature: Finite Temperature, 1: Zero Temperature.\n",
            defname);
    return (-1);
  }

  if(Def::iCalcType !=2 && Def::iOutputHam ==TRUE) {
    fprintf(MP::STDOUT,
            "Error in %s\n OutputHam is only defined for FullDiag mode, CalcType=2.\n",
            defname);
    return (-1);
  }

  return 0;
}

/**
 * @brief Function of Fitting FileName
 * @param[in]  cFileListNameFile file for getting names of input files.
 * @param[out] cFileNameList arrays for getting names of input files.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int GetFileName(
                const char* cFileListNameFile,
                char **cFileNameList
                )
{
  FILE *fplist;
  int itmpKWidx=-1;
  char ctmpFileName[D_FileNameMaxReadDef];
  char ctmpKW[D_CharTmpReadDef], ctmp2[256];
  int i;
  for(i=0; i< D_iKWNumDef; i++){
    strcpy(cFileNameList[i],"");
  }

  fplist = wrapperMPI::Fopen(cFileListNameFile, "r");
  if(fplist==NULL) return ReadDefFileError(cFileListNameFile);

  while(wrapperMPI::Fgets(ctmp2, 256, fplist) != NULL){
    memset(ctmpKW, '\0', strlen(ctmpKW));
    memset(ctmpFileName, '\0', strlen(ctmpFileName));
    sscanf(ctmp2,"%s %s\n", ctmpKW, ctmpFileName);

    if(strncmp(ctmpKW, "#", 1)==0 || *ctmp2=='\n' || (strcmp(ctmpKW, "")&&strcmp(ctmpFileName,""))==0){
      continue;
    }
    else if(strcmp(ctmpKW, "")*strcmp(ctmpFileName, "")==0){
      fprintf(MP::STDOUT,
              "Error: keyword and filename must be set as a pair in %s.\n",
              cFileListNameFile);
      fclose(fplist);
      return(-1);
    }
    /*!< Check KW */
    if( CheckKW(ctmpKW, cKWListOfFileNameList, D_iKWNumDef, &itmpKWidx)!=0 ){
      fprintf(MP::STDOUT, "Error: Wrong keywords %s in %s.\n", ctmpKW, cFileListNameFile);
      fprintf(MP::STDOUT, "%s", "Choose Keywords as follows: \n");
      for(i=0; i<D_iKWNumDef;i++){
        fprintf(MP::STDOUT, "%s \n", cKWListOfFileNameList[i]);
      }
      fclose(fplist);
      return(-1);
    }
    /*!< Check cFileNameList to prevent from double registering the file name */
    if(strcmp(cFileNameList[itmpKWidx], "") !=0){
      fprintf(MP::STDOUT, "Error: Same keywords exist in %s.\n", cFileListNameFile);
      fclose(fplist);
      return(-1);
    }

    /*!< Copy FileName */
    strcpy(cFileNameList[itmpKWidx], ctmpFileName);
  }
  fclose(fplist);  
  return 0;
}

/** 
 * @brief  Function of reading information about "ModPara" file and total number of parameters from other def files.
 *
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileNInt(
  char *xNameListFile/**<[in]  List of Input File names*/
)
{
  FILE *fp;
  char defname[D_FileNameMaxReadDef];
  char ctmp[D_CharTmpReadDef], ctmp2[256];
  int i, itmp;
  int iline = 0;
  Def::nvec = 0;
  Def::iFlgSpecOmegaMax = FALSE;
  Def::iFlgSpecOmegaMin = FALSE;
  Def::iFlgSpecOmegaOrg = FALSE;
  Def::iNOmega = 1000;
  Def::NCond = 0;
  Def::iFlgSzConserved = FALSE;
  Def::dcOmegaOrg = 0;
  int iReadNCond = FALSE;
  Boost::flgBoost = FALSE;
  InitializeInteractionNum();
  Step::NumAve = 1;
  Param::ExpecInterval = 1;
  cFileNameListFile = (char**)malloc(sizeof(char*)*D_iKWNumDef);
  for (i = 0; i < D_iKWNumDef; i++)
    cFileNameListFile[i] = (char*)malloc(sizeof(char)*D_CharTmpReadDef);

  fprintf(MP::STDOUT, "  Read File %s.\n", xNameListFile);
  if (GetFileName(xNameListFile, cFileNameListFile) != 0) {
    return(-1);
  }

  /*=======================================================================*/
  int iKWidx = 0;
  //Check the existence of Essensial Files.
  Def::READ = 0;
  Def::WRITE = 0;

  for (iKWidx = 0; iKWidx < D_iKWNumDef; iKWidx++) {
    strcpy(defname, cFileNameListFile[iKWidx]);
    if (strcmp(defname, "") == 0) {
      switch (iKWidx) {
      case KWCalcMod:
      case KWModPara:
      case KWLocSpin:
        fprintf(MP::STDOUT, "Error: Need to make a def file for %s.\n", cKWListOfFileNameList[iKWidx]);
        return(-1);
      default:
        break;
      }
    }
  }

  for (iKWidx = 0; iKWidx < D_iKWNumDef; iKWidx++) {
    strcpy(defname, cFileNameListFile[iKWidx]);

    if (strcmp(defname, "") == 0) continue;
    if (iKWidx == KWSpectrumVec) {
      continue;
    }
    fprintf(MP::STDOUT, "  Read File %s for %s.\n", defname, cKWListOfFileNameList[iKWidx]);
    fp = wrapperMPI::Fopen(defname, "r");
    if (fp == NULL) return ReadDefFileError(defname);
    switch (iKWidx) {
    case KWCalcMod:
      /* Read calcmod.def---------------------------------------*/
      if (ReadcalcmodFile(defname) != 0) {
        fclose(fp);
        return ReadDefFileError(defname);
      }
      break;

    case KWModPara:
      /* Read modpara.def---------------------------------------*/
      //TODO: add error procedure here when parameters are not enough.
      //! Read Header (5 lines).
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &itmp); //2
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //3
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //4
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //5
  //! Read header name for files about data
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %s\n", ctmp, Def::CDataFileHead); //6
  //! Read header name for files about parameters
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %s\n", ctmp, Def::CParaFileHead); //7
  //! Read header (1 line).
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);   //8
      double dtmp, dtmp2;
      Def::read_hacker = 1;
      //! Read lines.
      while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
        if (*ctmp2 == '\n') continue;
        sscanf(ctmp2, "%s %lf %lf\n", ctmp, &dtmp, &dtmp2);
        if (CheckWords(ctmp, "Nsite") == 0) {
          Def::Nsite = (int)dtmp;
        }
        else if (CheckWords(ctmp, "Nup") == 0) {
          Def::Nup = (int)dtmp;
        }
        else if (CheckWords(ctmp, "Ndown") == 0) {
          Def::Ndown = (int)dtmp;
          Def::Total2Sz = Def::Nup - Def::Ndown;
        }
        else if (CheckWords(ctmp, "2Sz") == 0) {
          Def::Total2Sz = (int)dtmp;
          Def::iFlgSzConserved = TRUE;
        }
        else if (CheckWords(ctmp, "Ncond") == 0) {
          if ((int)dtmp < 0) {
            fprintf(MP::STDOUT, "Error in %s\n Ncond must be greater than 0.\n", defname);
            return (-1);
          }
          Def::NCond = (int)dtmp;
          iReadNCond = TRUE;
        }
        else if (CheckWords(ctmp, "Lanczos_max") == 0) {
          Def::Lanczos_max = (int)dtmp;
        }
        else if (CheckWords(ctmp, "initial_iv") == 0) {
          Def::initial_iv = (int)dtmp;
        }
        else if (CheckWords(ctmp, "nvec") == 0) {
          Def::nvec = (int)dtmp;
        }
        else if (CheckWords(ctmp, "exct") == 0) {
          Def::k_exct = (int)dtmp;
        }
        else if (CheckWords(ctmp, "LanczosEps") == 0) {
          Def::LanczosEps = (int)dtmp;
        }
        else if (CheckWords(ctmp, "LanczosTarget") == 0) {
          Def::LanczosTarget = (int)dtmp;
        }
        else if (CheckWords(ctmp, "LargeValue") == 0) {
          Step::LargeValue = dtmp;
        }
        else if (CheckWords(ctmp, "NumAve") == 0) {
          Step::NumAve = (int)dtmp;
        }
        else if (strcmp(ctmp, "TimeSlice") == 0) {
          Param::TimeSlice = dtmp;
        }
        else if (strcmp(ctmp, "ExpandCoef") == 0) {
          Param::ExpandCoef = (int)dtmp;
        }
        else if (strcmp(ctmp, "OutputInterval") == 0) {
          Param::OutputInterval = (int)dtmp;
        }
        else if (CheckWords(ctmp, "ExpecInterval") == 0) {
          Param::ExpecInterval = (int)dtmp;
        }
        else if (strcmp(ctmp, "Tinit") == 0) {
          Param::Tinit = dtmp;
        }
        else if (CheckWords(ctmp, "CalcHS") == 0) {
          Def::read_hacker = (int)dtmp;
        }
        else if (CheckWords(ctmp, "OmegaMax") == 0) {
          Def::dcOmegaMax = std::complex<double>(dtmp, dtmp2);
          Def::iFlgSpecOmegaMax = TRUE;
        }
        else if (CheckWords(ctmp, "OmegaMin") == 0) {
          Def::dcOmegaMin = std::complex<double>(dtmp, dtmp2);
          Def::iFlgSpecOmegaMin = TRUE;
        }
        else if (CheckWords(ctmp, "OmegaIm") == 0) {
          Def::dcOmegaOrg += std::complex<double>(0.0, dtmp);
          Def::iFlgSpecOmegaOrg = TRUE;
        }
        else if (CheckWords(ctmp, "OmegaOrg") == 0) {
          Def::dcOmegaOrg += std::complex<double>(dtmp, dtmp2);
          Def::iFlgSpecOmegaOrg = TRUE;
        }
        else if (CheckWords(ctmp, "NOmega") == 0) {
          Def::iNOmega = (int)dtmp;
        }
        else if (CheckWords(ctmp, "TargetTPQRand") == 0) {
          Def::irand = (int)dtmp;
        }
        else {
          return (-1);
        }
      }
      break;

    case KWLocSpin:
      // Read locspn.def
      Def::iFlgGeneralSpin = FALSE;
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NLocSpn));
      break;
    case KWTrans:
      // Read transfer.def
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NTransfer));
      break;
    case KWCoulombIntra:
      /* Read coulombintra.def----------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NCoulombIntra));
      break;
    case KWCoulombInter:
      /* Read coulombinter.def----------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NCoulombInter));
      break;
    case KWHund:
      /* Read hund.def------------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NHundCoupling));
      break;
    case KWPairHop:
      /* Read pairhop.def---------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NPairHopping));
      Def::NPairHopping *= 2;
      break;
    case KWExchange:
      /* Read exchange.def--------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NExchangeCoupling));
      break;
    case KWIsing:
      /* Read ising.def--------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NIsingCoupling));
      break;
    case KWPairLift:
      /* Read exchange.def--------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NPairLiftCoupling));
      break;
    case KWInterAll:
      /* Read InterAll.def--------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NInterAll));
      break;
    case KWOneBodyG:
      /* Read cisajs.def----------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NCisAjt));
      break;
    case KWTwoBodyG:
      /* Read cisajscktaltdc.def--------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NCisAjtCkuAlvDC));
      break;
    case KWLaser:
      /* Read laser.def--------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NLaser));
      break;

    case KWTEOneBody: {
      if (Def::iCalcType != DC::TimeEvolution) break;
      /* Read TEOnebody.def--------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NTETimeSteps));
      wrapperMPI::Fgets(ctmp2, 256, fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      int iTETransMax = 0;
      if (Def::NTETimeSteps > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%lf %d \n", &dtmp, &itmp);
          for (i = 0; i < itmp; ++i) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
          }
          if (iTETransMax < itmp) iTETransMax = itmp;
        }
      }
      Def::NTETransferMax = iTETransMax;
      break;
    }
    case KWTETwoBody: {
      if (Def::iCalcType != DC::TimeEvolution) break;
      /* Read TETwobody.def--------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NTETimeSteps));
      wrapperMPI::Fgets(ctmp2, 256, fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      int iTEInterAllMax = 0;
      if (Def::NTETimeSteps > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%lf %d \n", &dtmp, &itmp);
          for (i = 0; i < itmp; ++i) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
          }
          if (iTEInterAllMax < itmp) iTEInterAllMax = itmp;
        }
      }
      Def::NTEInterAllMax = iTEInterAllMax;
      break;
    }

    case KWBoost:
      /* Read boost.def--------------------------------*/
      Boost::NumarrayJ = 0;
      Boost::W0 = 0;
      Boost::R0 = 0;
      Boost::num_pivot = 0;
      Boost::ishift_nspin = 0;
      Boost::flgBoost = TRUE;
      //first line is skipped
      wrapperMPI::Fgets(ctmp2, 256, fp);
      //read numarrayJ
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%d\n", &(Boost::NumarrayJ));
      //skipp arrayJ
      for (iline = 0; iline < Boost::NumarrayJ * 3; iline++) {
        wrapperMPI::Fgets(ctmp2, 256, fp);
      }
      //read W0 R0 num_pivot ishift_nspin
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%ld %ld %ld %ld\n", &(Boost::W0), &(Boost::R0), &(Boost::num_pivot),
        &(Boost::ishift_nspin));

      break;

    case KWSingleExcitation:
      /* Read singleexcitation.def----------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NNSingleExcitationOperator));
      break;

    case KWPairExcitation:
      /* Read pairexcitation.def----------------------------------------*/
      wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%s %d\n", ctmp, &(Def::NNPairExcitationOperator));
      break;

    default:
      fprintf(MP::STDOUT, "%s", "Error: incorrect file.\n");
      fclose(fp);
      return (-1);
    }
    /*=======================================================================*/
    fclose(fp);
  }

  //Sz, Ncond
  switch (Def::iCalcModel) {
  case DC::Spin:
  case DC::Hubbard:
  case DC::Kondo:
  case DC::SpinlessFermion:

    if (iReadNCond == TRUE) {
      if (Def::iCalcModel == DC::Spin) {
        fprintf(MP::STDOUT, "For Spin, Ncond should not be defined.\n");
        return(-1);
      }
      else {
        if (Def::iFlgSzConserved == TRUE) {
          if (Def::iCalcModel == DC::SpinlessFermion) {
            fprintf(MP::STDOUT, "  Warning: For Spinless fermion, 2Sz should not be defined.\n");
            Def::Ne = Def::NCond;
            Def::Nup = Def::NCond;
            Def::Ndown = 0;
            break;
          }
          Def::Nup = Def::NLocSpn + Def::NCond + Def::Total2Sz;
          Def::Ndown = Def::NLocSpn + Def::NCond - Def::Total2Sz;
          Def::Nup /= 2;
          Def::Ndown /= 2;
          Def::Ne = Def::Nup + Def::Ndown;
        }
        else {
          if (Def::iCalcModel == DC::Hubbard) {
            Def::Ne = Def::NCond;
            if (Def::Ne < 1) {
              fprintf(MP::STDOUT, "Ncond is incorrect.\n");
              return(-1);
            }
            Def::iCalcModel = DC::HubbardNConserved;
          }
          else if (Def::iCalcModel == DC::SpinlessFermion) {
            Def::Ne = Def::NCond;
            Def::Nup = Def::NCond;
            Def::Ndown = 0;
          }
          else {
            fprintf(MP::STDOUT, " 2Sz is not defined.\n");
            return(-1);
          }
        }
      }
    }
    else if (iReadNCond == FALSE && Def::iFlgSzConserved == TRUE) {
      if (Def::iCalcModel != DC::Spin) {
        fprintf(MP::STDOUT, " NCond is not defined.\n");
        return(-1);
      }
      Def::Nup = Def::NLocSpn + Def::Total2Sz;
      Def::Ndown = Def::NLocSpn - Def::Total2Sz;
      Def::Nup /= 2;
      Def::Ndown /= 2;
    }
    else {
      if (Def::Nup == 0 && Def::Ndown == 0) {
        if (Def::iCalcModel == DC::Spin) {
          fprintf(MP::STDOUT, " 2Sz is not defined.\n");
          return(-1);
        }
        else {
          fprintf(MP::STDOUT, " NCond is not defined.\n");
          return(-1);
        }
      }
    }

    if (Def::iCalcModel == DC::Spin) {
      Def::Ne = Def::Nup;
    }
    else {
      if (Def::Ne == 0) {
        Def::Ne = Def::Nup + Def::Ndown;
      }
      if (Def::NLocSpn > Def::Ne) {
        fprintf(MP::STDOUT, "%s", "Error: Ne=Nup+Ndown must be (Ne >= NLocalSpin).\n");
        fprintf(MP::STDOUT, "NLocalSpin=%d, Ne=%d\n", Def::NLocSpn, Def::Ne);
        return(-1);
      }
    }
    break;
  case DC::SpinGC:
  case DC::KondoGC:
  case DC::HubbardGC:
  case DC::SpinlessFermionGC:
    if (iReadNCond == TRUE || Def::iFlgSzConserved == TRUE) {
      fprintf(MP::STDOUT, "\n  Warning: For GC, both Ncond and 2Sz should not be defined.\n");
      //return(-1);
    }
    break;
  default:
    break;
  }

  /* Check values (Positive)*/
  if (Def::Nsite <= 0) {// Nsite must be positve
    fprintf(MP::STDOUT, "Error in %s\n Nsite must be positive value.\n", defname);
    return (-1);
  }
  if (Def::Lanczos_max <= 0) {// Lanczos_max must be positive
    fprintf(MP::STDOUT, "Error in %s\n Lanczos_max must be positive value.\n", defname);
    return (-1);
  }
  if (Def::LanczosEps <= 0) {// Lanczos_eps must be positive
    fprintf(MP::STDOUT, "Error in %s\n Lanczos_eps must be positive value.\n", defname);
    return (-1);
  }
  if (Step::NumAve <= 0) { // Average number must be positive
    fprintf(MP::STDOUT, "Error in %s\n Step::NumAve must be positive value.\n", defname);
    return (-1);
  }
  if (Param::ExpecInterval <= 0) {// Interval to calculate expected values must be positive
    fprintf(MP::STDOUT, "Error in %s\n ExpecInterval must be positive value.\n", defname);
    return (-1);
  }
  if (Def::nvec == 0) {
    Def::nvec = Def::Lanczos_max;
  }

  if (Def::nvec < Def::k_exct) {
    Def::nvec = Def::k_exct;
  }
  if (Def::LanczosTarget < Def::k_exct) Def::LanczosTarget = Def::k_exct;

  if (ValidateValue(Def::k_exct, 1, Def::nvec)) {
    fprintf(MP::STDOUT, "Error in %s \n exct=%d must be greater than 1 and smaller than nvec=%d.\n", 
      defname, Def::k_exct, Def::nvec);
    return (-1);
  }

  if (Def::k_exct > Def::LanczosTarget) {
    fprintf(MP::STDOUT, "Error in %s\n LanczosTarget=%d must be greater than exct=%d.\n", defname, Def::LanczosTarget, Def::k_exct);
    return (-1);
  }

  Def::fidx = 0;
  Def::NeMPI = Def::Ne;
  Def::NupMPI = Def::Nup;
  Def::NdownMPI = Def::Ndown;
  Def::NupOrg = Def::Nup;
  Def::NdownOrg = Def::Ndown;
  return 0;
}
///
/// \brief function of checking hermite conditions about interall interactions
/// \param InterAll arrays of information of interall interactions
/// \param ParaInterAll arrays of values of interall interactions
/// \param InterAllOffDiagonal arrays of information of off-diagonal part of interall interactions
/// \param ParaInterAllOffDiagonal arrays of values of off-diagonal part of interall interactions
/// \param NInterAllOffDiagonal total number of off-diagonal part of interall interactions
/// \param iCalcModel Target Model defined in CalcMod file (ex. Spin, SpinGC etc.)
/// \retval 0 Hermite condition is satisfied
/// \retval -1 Hermite condition is not satisfied
/// @version 0.2
/// @details rearray a InterAll_OffDiagonal array to satisfy a condition of hermite conjugation between 2*i and 2*i+1 components.
///
/// @version 0.1
/// @author Takahiro Misawa (The University of Tokyo)
/// @author Kazuyoshi Yoshimi (The University of Tokyo)
static int CheckInterAllHermite
(
  int **InterAll,
  std::complex<double>* ParaInterAll,
  int **InterAllOffDiagonal,
  std::complex<double>*ParaInterAllOffDiagonal,
  const int NInterAllOffDiagonal,
  const int iCalcModel
) {
  int i, j, icntincorrect, itmpret;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  int itmpsite1, itmpsite2, itmpsite3, itmpsite4;
  int itmpsigma1, itmpsigma2, itmpsigma3, itmpsigma4;
  int itmpIdx, icntHermite;
  int icheckHermiteCount = FALSE;
  std::complex<double> ddiff_intall;
  icntincorrect = 0;
  icntHermite = 0;
  for (i = 0; i < NInterAllOffDiagonal; i++) {
    itmpret = 0;
    isite1 = InterAllOffDiagonal[i][0];
    isigma1 = InterAllOffDiagonal[i][1];
    isite2 = InterAllOffDiagonal[i][2];
    isigma2 = InterAllOffDiagonal[i][3];
    isite3 = InterAllOffDiagonal[i][4];
    isigma3 = InterAllOffDiagonal[i][5];
    isite4 = InterAllOffDiagonal[i][6];
    isigma4 = InterAllOffDiagonal[i][7];
    icheckHermiteCount = FALSE;

    for (j = 0; j < NInterAllOffDiagonal; j++) {
      itmpsite1 = InterAllOffDiagonal[j][0];
      itmpsigma1 = InterAllOffDiagonal[j][1];
      itmpsite2 = InterAllOffDiagonal[j][2];
      itmpsigma2 = InterAllOffDiagonal[j][3];
      itmpsite3 = InterAllOffDiagonal[j][4];
      itmpsigma3 = InterAllOffDiagonal[j][5];
      itmpsite4 = InterAllOffDiagonal[j][6];
      itmpsigma4 = InterAllOffDiagonal[j][7];

      if (isite1 == itmpsite4 && isite2 == itmpsite3 && isite3 == itmpsite2 && isite4 == itmpsite1) {
        if (isigma1 == itmpsigma4 && isigma2 == itmpsigma3 && isigma3 == itmpsigma2 && isigma4 == itmpsigma1) {
          ddiff_intall = abs(ParaInterAllOffDiagonal[i] - conj(ParaInterAllOffDiagonal[j]));

          if (abs(ddiff_intall) < eps_CheckImag0) {
            itmpret = 1;
            if (icheckHermiteCount == FALSE) {
              icheckHermiteCount = TRUE; //for not double counting
              if (i <= j) {
                if (2 * icntHermite >= NInterAllOffDiagonal) {
                  fprintf(MP::STDOUT, "Elements of InterAll are incorrect.\n");
                  return (-1);
                }

                for (itmpIdx = 0; itmpIdx < 8; itmpIdx++) {
                  InterAll[2 * icntHermite][itmpIdx] = InterAllOffDiagonal[i][itmpIdx];
                  InterAll[2 * icntHermite + 1][itmpIdx] = InterAllOffDiagonal[j][itmpIdx];
                }

                ParaInterAll[2 * icntHermite] = ParaInterAllOffDiagonal[i];
                ParaInterAll[2 * icntHermite + 1] = ParaInterAllOffDiagonal[j];
                icntHermite++;
              }
              break;
            }
          }
        }
      }
      else if (isite1 == itmpsite2 && isite2 == itmpsite1 && isite3 == itmpsite4 &&
        isite4 == itmpsite3) {      //for spin and Kondo
        if (iCalcModel == DC::Kondo || iCalcModel == DC::KondoGC || iCalcModel == DC::Spin || iCalcModel == DC::SpinGC) {
          if (isigma1 == itmpsigma2 && isigma2 == itmpsigma1 && isigma3 == itmpsigma4 && isigma4 == itmpsigma3) {
            ddiff_intall = ParaInterAllOffDiagonal[i] - conj(ParaInterAllOffDiagonal[j]);
            if (abs(ddiff_intall) < eps_CheckImag0) {
              itmpret = 1;
              if (icheckHermiteCount == FALSE) {
                icheckHermiteCount = TRUE; // for not double-counting
                if (i <= j) {
                  if (2 * icntHermite >= NInterAllOffDiagonal) {
                    fprintf(MP::STDOUT, "Elements of InterAll are incorrect.\n");
                    return (-1);
                  }
                  for (itmpIdx = 0; itmpIdx < 8; itmpIdx++) {
                    InterAll[2 * icntHermite][itmpIdx] = InterAllOffDiagonal[i][itmpIdx];
                  }
                  for (itmpIdx = 0; itmpIdx < 4; itmpIdx++) {
                    InterAll[2 * icntHermite + 1][2 * itmpIdx] = InterAllOffDiagonal[i][6 -
                      2 *
                      itmpIdx];
                    InterAll[2 * icntHermite + 1][2 * itmpIdx + 1] = InterAllOffDiagonal[i][7 - 2 *
                      itmpIdx];

                  }
                  ParaInterAll[2 * icntHermite] = ParaInterAllOffDiagonal[i];
                  ParaInterAll[2 * icntHermite + 1] = ParaInterAllOffDiagonal[j];
                  icntHermite++;
                }
                break;
              }
            }
          }
        }
      }
    }
    //if counterpart for satisfying hermite conjugate does not exist.
    if (itmpret != 1) {
      fprintf(MP::STDOUT, "Error: NonHermite (i, spni, j, spnj, k, spnk, l, spnl) = (%d, %d, %d, %d, %d, %d, %d, %d), InterAll_re= %lf, InterAll_im= %lf . \n", isite1, isigma1, isite2, isigma2, isite3, isigma3, isite4, isigma4,
        real(ParaInterAllOffDiagonal[i]), imag(ParaInterAllOffDiagonal[i]));
      icntincorrect++;
    }
  }

  if (icntincorrect != 0) {
    return (-1);
  }

  for (i = 0; i < NInterAllOffDiagonal; i++) {
    for (itmpIdx = 0; itmpIdx < 8; itmpIdx++) {
      InterAllOffDiagonal[i][itmpIdx] = InterAll[i][itmpIdx];
    }
    ParaInterAllOffDiagonal[i] = ParaInterAll[i];
  }

  return 0;
}/*CheckInterAllHermite*/
/// \brief function of getting diagonal components
/// \param InterAll  arrays of information of interall interactions
/// \param ParaInterAll arrays of values of interall interactions
/// \param NInterAll total number of interall interactions
/// \param InterAllDiagonal arrays of information of diagonal part of interall interactions
/// \param ParaInterAllDiagonal arrays of values of diagonal part of interall interactions
/// \param InterAllOffDiagonal arrays of information of off-diagonal part of interall interactions
/// \param ParaInterAllOffDiagonal arrays of values of off-diagonal part of interall interactions
/// \param Chemi arrays of the site of chemical potential
/// \param SpinChemi arrays of the spin of chemical potential
/// \param ParaChemi arrays of the value of chemical potential
/// \param NChemi total number of chemical potential
/// \param iCalcModel Target Model defined in CalcMod file (ex. Spin, SpinGC etc.)
/// \retval 0 succeed to get diagonal interactions.
/// \retval -1 format of interall interactions is incorrect.
/// \version 2.1
/// \author Kazuyoshi Yoshimi (The University of Tokyo)
int GetDiagonalInterAll
(
  int **InterAll,
  std::complex<double> *ParaInterAll,
  const int NInterAll,
  int **InterAllDiagonal,
  double *ParaInterAllDiagonal,
  int **InterAllOffDiagonal,
  std::complex<double> *ParaInterAllOffDiagonal,
  int *Chemi,
  int *SpinChemi,
  double *ParaChemi,
  int *NChemi,
  const int iCalcModel
)
{
  int i, icnt_diagonal, icnt_offdiagonal, tmp_i;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  int iret = 0;
  icnt_diagonal = 0;
  icnt_offdiagonal = 0;

  for (i = 0; i < NInterAll; i++) {
    isite1 = InterAll[i][0];
    isigma1 = InterAll[i][1];
    isite2 = InterAll[i][2];
    isigma2 = InterAll[i][3];
    isite3 = InterAll[i][4];
    isigma3 = InterAll[i][5];
    isite4 = InterAll[i][6];
    isigma4 = InterAll[i][7];

    //Get Diagonal term
    if (isite1 == isite2 && isite3 == isite4 &&
      isigma1 == isigma2 && isigma3 == isigma4)
    {
      InterAllDiagonal[icnt_diagonal][0] = isite1;
      InterAllDiagonal[icnt_diagonal][1] = isigma1;
      InterAllDiagonal[icnt_diagonal][2] = isite3;
      InterAllDiagonal[icnt_diagonal][3] = isigma3;
      ParaInterAllDiagonal[icnt_diagonal] = real(ParaInterAll[i]);
      icnt_diagonal++;
      continue;
    }
    else if (isite1 == isite4 && isite2 == isite3 &&
      isigma1 == isigma4 && isigma2 == isigma3)
    {
      InterAllDiagonal[icnt_diagonal][0] = isite1;
      InterAllDiagonal[icnt_diagonal][1] = isigma1;
      InterAllDiagonal[icnt_diagonal][2] = isite2;
      InterAllDiagonal[icnt_diagonal][3] = isigma2;
      ParaInterAllDiagonal[icnt_diagonal] = -real(ParaInterAll[i]);
      Chemi[*NChemi] = isite1;
      SpinChemi[*NChemi] = isigma1;
      //transfer integral has minus sign for default setting
      ParaChemi[*NChemi] = -real(ParaInterAll[i]);
      icnt_diagonal++;
      *NChemi += 1;
      continue;
    }
    else {
      //Get Off-Diagonal term
      switch (iCalcModel) {
      case DC::Hubbard:
      case DC::HubbardNConserved:
      case DC::Kondo:
      case DC::KondoGC:
      case DC::HubbardGC:
        if (isigma1 == isigma2 && isigma3 == isigma4) {
          for (tmp_i = 0; tmp_i < 8; tmp_i++) {
            InterAllOffDiagonal[icnt_offdiagonal][tmp_i] = InterAll[i][tmp_i];
          }
          ParaInterAllOffDiagonal[icnt_offdiagonal] = ParaInterAll[i];
        }
        else if (isigma1 == isigma4 && isigma2 == isigma3) {
          InterAllOffDiagonal[icnt_offdiagonal][0] = isite1;
          InterAllOffDiagonal[icnt_offdiagonal][1] = isigma1;
          InterAllOffDiagonal[icnt_offdiagonal][2] = isite4;
          InterAllOffDiagonal[icnt_offdiagonal][3] = isigma1;
          InterAllOffDiagonal[icnt_offdiagonal][4] = isite3;
          InterAllOffDiagonal[icnt_offdiagonal][5] = isigma2;
          InterAllOffDiagonal[icnt_offdiagonal][6] = isite2;
          InterAllOffDiagonal[icnt_offdiagonal][7] = isigma2;
          ParaInterAllOffDiagonal[icnt_offdiagonal] = -ParaInterAll[i];
        }
        else {
          // Sz symmetry is assumed
          if (iCalcModel == DC::Hubbard || iCalcModel == DC::Kondo) {
            fprintf(MP::STDOUT, 
              "Error: This operator breaks Sz Symmetry (i, spni, j, spnj, k, spnk, l, spnl) = (%d, %d, %d, %d, %d, %d, %d, %d), InterAll_re= %lf, InterAll_im= %lf . \n",
              isite1,
              isigma1,
              isite2,
              isigma2,
              isite3,
              isigma3,
              isite4,
              isigma4,
              real(ParaInterAll[i]),
              imag(ParaInterAll[i])
            );
            iret = -1;
          }
          else {
            for (tmp_i = 0; tmp_i < 8; tmp_i++) {
              InterAllOffDiagonal[icnt_offdiagonal][tmp_i] = InterAll[i][tmp_i];
            }
            ParaInterAllOffDiagonal[icnt_offdiagonal] = ParaInterAll[i];
          }
        }
        break;
      case DC::Spin:
      case DC::SpinGC:
        if (isite1 == isite2 && isite3 == isite4) {
          for (tmp_i = 0; tmp_i < 8; tmp_i++) {
            InterAllOffDiagonal[icnt_offdiagonal][tmp_i] = InterAll[i][tmp_i];
          }
          ParaInterAllOffDiagonal[icnt_offdiagonal] = ParaInterAll[i];
        }
        break;
      default:
        return(-1);
      }
      if (iret != -1) {
        icnt_offdiagonal++;
      }
    }

    if (iret != 0) {
      return(-1);
    }
  }

  return 0;
}/*GetDiagonalInterAll*/
/**
 * @brief Check Hermite for TETransfer integrals.
 * @param[in] X Define List for getting transfer integrals.
 * @param[in] NTETransfer total number of transfer integrals
 * @param[in] idx index for time step.
 * @retval 0 Hermite.
 * @retval -1 NonHermite.
 **/
static int CheckTETransferHermite
(
  
  const int NTETransfer,
  const int idx
)
{
  int i, j;
  int isite1, isite2;
  int isigma1, isigma2;
  int itmpsite1, itmpsite2;
  int itmpsigma1, itmpsigma2;
  int itmperrsite1, itmperrsite2;
  int itmperrsigma1, itmperrsigma2;
  std::complex<double> dcerrTrans;
  int icheckHermiteCount;
  int iCount = 0;

  std::complex<double> ddiff_trans;
  int itmpIdx, icntHermite, icntchemi;
  icntHermite = 0;
  icntchemi = 0;

  int** tmp_TETransfer = i_2d_allocate(NTETransfer, 4);
  std::complex<double>*tmp_paraTETransfer = (std::complex<double>*)malloc((NTETransfer) * sizeof(std::complex<double>));

  //copy
  for (i = 0; i < NTETransfer; i++) {
    for (j = 0; j < 4; j++) {
      tmp_TETransfer[i][j] = Def::TETransfer[idx][i][j];
      Def::TETransfer[idx][i][j] = 0;
    }
    tmp_paraTETransfer[i] = Def::ParaTETransfer[idx][i];
    Def::ParaTETransfer[idx][i] = 0.0;
  }

  for (i = 0; i < NTETransfer; i++) {
    isite1 = tmp_TETransfer[i][0];
    isigma1 = tmp_TETransfer[i][1];
    isite2 = tmp_TETransfer[i][2];
    isigma2 = tmp_TETransfer[i][3];
    icheckHermiteCount = FALSE;
    for (j = 0; j < NTETransfer; j++) {
      itmpsite1 = tmp_TETransfer[j][0];
      itmpsigma1 = tmp_TETransfer[j][1];
      itmpsite2 = tmp_TETransfer[j][2];
      itmpsigma2 = tmp_TETransfer[j][3];
      if (isite1 == itmpsite2 && isite2 == itmpsite1) {
        if (isigma1 == itmpsigma2 && isigma2 == itmpsigma1) {

          ddiff_trans = tmp_paraTETransfer[i] - conj(tmp_paraTETransfer[j]);
          if (abs(ddiff_trans) > eps_CheckImag0) {
            itmperrsite1 = itmpsite1;
            itmperrsigma1 = itmpsigma1;
            itmperrsite2 = itmpsite2;
            itmperrsigma2 = itmpsigma2;
            dcerrTrans = tmp_paraTETransfer[j];
            fprintf(MP::STDOUT, 
              "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n", 
              isite1, isigma1, isite2, isigma2, real(tmp_paraTETransfer[i]), imag(tmp_paraTETransfer[i]));
            fprintf(MP::STDOUT,
              "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n",
              itmperrsite1, itmperrsigma1, itmperrsite2, itmperrsigma2, real(dcerrTrans), imag(dcerrTrans));
            iCount++;
          }
          else {
            if (icheckHermiteCount == FALSE) {
              if (i <= j) {
                if (2 * icntHermite >= NTETransfer) {
                  fprintf(stderr, "Elements of Transfers are incorrect.\n");
                  return(-1);
                }
                if (isite1 != isite2 || isigma1 != isigma2) {
                  for (itmpIdx = 0; itmpIdx < 4; itmpIdx++) {
                    Def::TETransfer[idx][2 * icntHermite][itmpIdx] = tmp_TETransfer[i][itmpIdx];
                    Def::TETransfer[idx][2 * icntHermite + 1][itmpIdx] = tmp_TETransfer[j][itmpIdx];
                  }
                  Def::ParaTETransfer[idx][2 * icntHermite] = tmp_paraTETransfer[i];
                  Def::ParaTETransfer[idx][2 * icntHermite + 1] = tmp_paraTETransfer[j];
                  icntHermite++;
                }
                else {
                  Def::TETransferDiagonal[idx][icntchemi][0] = tmp_TETransfer[i][0];
                  Def::TETransferDiagonal[idx][icntchemi][1] = tmp_TETransfer[i][1];
                  Def::ParaTETransferDiagonal[idx][icntchemi] = real(tmp_paraTETransfer[i]);
                  icntchemi += 1;
                }
              }
              icheckHermiteCount = TRUE;
            }
          }
        }
      }
    }

    //if counterpart for satisfying hermite conjugate does not exist.
    if (icheckHermiteCount == FALSE) {
      fprintf(MP::STDOUT,
        "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n",
        isite1, isigma1, isite2, isigma2, 
        real(tmp_paraTETransfer[i]), imag(tmp_paraTETransfer[i]));
      iCount++;
      //fprintf(MP::STDOUT, "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n", itmperrsite1, itmperrsigma1, itmperrsite2, itmperrsigma2, real(dcerrTrans), imag(dcerrTrans));
      //return(-1);
    }
  }

  if (iCount != 0) {
    return -1;
  }

  Def::NTETransfer[idx] = 2 * icntHermite;
  Def::NTETransferDiagonal[idx] = icntchemi;


  free_i_2d_allocate(tmp_TETransfer);
  free(tmp_paraTETransfer);
  return 0;
}/*CheckTETransferHermite*/
/**
 * @brief function of reading def files to get keyword index
 * 
 * @param X define list to get and put informations for calcuation
 * @param xBoost list to get and put informations for Boost (CMA) calcuation
 *
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileIdxPara()
{
  FILE *fp;
  char defname[D_FileNameMaxReadDef];
  char ctmp[D_CharTmpReadDef], ctmp2[256];

  int i, idx, itype;
  int xitmp[8];
  int iKWidx = 0;
  int iboolLoc = 0;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  double dvalue_re, dvalue_im;
  double dArrayValue_re[3];
  int icnt_diagonal = 0;
  int ieps_CheckImag0 = -12;
  eps_CheckImag0 = pow(10.0, ieps_CheckImag0);
  int iline = 0;
  int ilineIn = 0;
  int ilineIn2 = 0;
  int itmp = 0;
  int icnt_trans = 0;
  int iflg_trans = 0;
  int icnt_interall = 0;

  int iloop = 0;

  for (iKWidx = KWLocSpin; iKWidx < D_iKWNumDef; iKWidx++) {
    strcpy(defname, cFileNameListFile[iKWidx]);
    if (strcmp(defname, "") == 0 || iKWidx == KWSpectrumVec) continue;
    fprintf(MP::STDOUT, "  Read File %s.\n", defname);
    fp = wrapperMPI::Fopen(defname, "r");
    if (fp == NULL) return ReadDefFileError(defname);
    if (iKWidx != KWBoost) {
      for (i = 0; i < IgnoreLinesInDef; i++) wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
    }

    idx = 0;
    /*=======================================================================*/
    switch (iKWidx) {
    case KWLocSpin:
      /* Read locspn.def----------------------------------------*/
      while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
        if (idx == Def::Nsite) {
          fclose(fp);
          return ReadDefFileError(defname);
        }

        sscanf(ctmp2, "%d %d\n", &(xitmp[0]), &(xitmp[1]));
        Def::LocSpn[xitmp[0]] = xitmp[1];
        Def::SiteToBit[xitmp[0]] = (Def::LocSpn[xitmp[0]] + 1);//2S+1
        if (CheckSite(xitmp[0], Def::Nsite) != 0) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
        idx++;
      }
      if (CheckLocSpin() == FALSE) {
        fclose(fp);
        return ReadDefFileError(defname);
      }

      break;

    case KWTrans:
      /* transfer.def--------------------------------------*/
      if (Def::NTransfer > 0) {
        icnt_trans = 0;
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL)
        {
          if (idx == Def::NTransfer) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %d %d %lf %lf\n",
            &isite1,
            &isigma1,
            &isite2,
            &isigma2,
            &dvalue_re,
            &dvalue_im
          );

          if (CheckPairSite(isite1, isite2, Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          if (isite1 == isite2 && isigma1 == isigma2) {
            if (fabs(dvalue_im) > eps_CheckImag0) {
              //NonHermite
              fprintf(MP::STDOUT,
                "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n",
                isite1, isigma1, isite2, isigma2, dvalue_re, dvalue_im);
              fclose(fp);
              return ReadDefFileError(defname);
            }
          }

          if (Def::iCalcModel == DC::Spin) {
            if (isite1 != isite2) {
              iboolLoc = 1;
              fprintf(MP::STDOUT, "Warning: Site component of (i, j) =(%d, %d) is ignored.\n",
                isite1, isite2);
            }
          }
          else if (Def::iCalcModel == DC::Kondo) {
            if (Def::LocSpn[isite1] != ITINERANT || Def::LocSpn[isite2] != ITINERANT) {
              if (isite1 != isite2) {
                iboolLoc = 1;
                fprintf(MP::STDOUT, "Error: Site component of (i, j) =(%d, %d) is incorrect.\n", isite1, isite2);
              }
            }
          }
          else if (Def::iCalcModel == DC::SpinlessFermion || Def::iCalcModel == DC::SpinlessFermionGC) {
            if (isigma1 != 0 || isigma2 != 0) {
              //Not allowed
              fprintf(stderr,
                "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n",
                isite1, isigma1, isite2, isigma2, dvalue_re, dvalue_im);
              fclose(fp);
              return ReadDefFileError(defname);
            }
          }

          iflg_trans = 0;
          for (i = 0; i < icnt_trans; i++) {
            if (isite1 == Def::GeneralTransfer[i][0] && isite2 == Def::GeneralTransfer[i][2]
              && isigma1 == Def::GeneralTransfer[i][1] && isigma2 == Def::GeneralTransfer[i][3])
            {
              Def::ParaGeneralTransfer[i] += std::complex<double>(dvalue_re, dvalue_im);
              iflg_trans = 1;
              continue;
            }
          }

          if (iflg_trans == 0) {
            Def::GeneralTransfer[icnt_trans][0] = isite1;
            Def::GeneralTransfer[icnt_trans][1] = isigma1;
            Def::GeneralTransfer[icnt_trans][2] = isite2;
            Def::GeneralTransfer[icnt_trans][3] = isigma2;
            Def::ParaGeneralTransfer[icnt_trans] = std::complex<double>(dvalue_re, dvalue_im);
            icnt_trans++;
          }
          idx++;
        }

        if (iboolLoc == 1) {
          fclose(fp);
          return(-1);
        }
      }

      Def::NTransfer = icnt_trans;

      if (CheckSpinIndexForTrans() == FALSE) {
        fclose(fp);
        return(-1);
      }

      if (CheckTransferHermite() != 0) {
        fprintf(MP::STDOUT, "%s", "Error: NonHermite Pair exists in transfer integral. \n");
        fclose(fp);
        return(-1);
      }
      break;

    case KWCoulombIntra:
      /*coulombintra.def----------------------------------*/
      if (Def::NCoulombIntra > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NCoulombIntra) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %lf\n",
            &(Def::CoulombIntra[idx][0]),
            &(Def::ParaCoulombIntra[idx])
          );

          if (CheckSite(Def::CoulombIntra[idx][0], Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          idx++;
        }
      }
      break;

    case KWCoulombInter:
      /*coulombinter.def----------------------------------*/
      if (Def::NCoulombInter > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NCoulombInter) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
            &(Def::CoulombInter[idx][0]),
            &(Def::CoulombInter[idx][1]),
            &(Def::ParaCoulombInter[idx])
          );

          if (CheckPairSite(Def::CoulombInter[idx][0], Def::CoulombInter[idx][1], Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;

    case KWHund:
      /*hund.def------------------------------------------*/
      if (Def::NHundCoupling > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL)
        {
          if (idx == Def::NHundCoupling) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
            &(Def::HundCoupling[idx][0]),
            &(Def::HundCoupling[idx][1]),
            &(Def::ParaHundCoupling[idx])
          );

          if (CheckPairSite(Def::HundCoupling[idx][0], Def::HundCoupling[idx][1], Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;
    case KWPairHop:
      /*pairhop.def---------------------------------------*/
      if (Def::iCalcModel == DC::Spin || Def::iCalcModel == DC::SpinGC) {
        fprintf(MP::STDOUT, "PairHop is not active in Spin and SpinGC.\n");
        return(-1);
      }

      if (Def::NPairHopping > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NPairHopping / 2) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %d %lf\n",
            &(Def::PairHopping[2 * idx][0]),
            &(Def::PairHopping[2 * idx][1]),
            &(Def::ParaPairHopping[2 * idx])
          );

          if (CheckPairSite(Def::PairHopping[2 * idx][0], Def::PairHopping[2 * idx][1], Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          Def::PairHopping[2 * idx + 1][0] = Def::PairHopping[2 * idx][1];
          Def::PairHopping[2 * idx + 1][1] = Def::PairHopping[2 * idx][0];
          Def::ParaPairHopping[2 * idx + 1] = Def::ParaPairHopping[2 * idx];
          idx++;
        }
      }
      break;

    case KWExchange:
      /*exchange.def--------------------------------------*/
      if (Def::NExchangeCoupling > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NExchangeCoupling) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
            &(Def::ExchangeCoupling[idx][0]),
            &(Def::ExchangeCoupling[idx][1]),
            &(Def::ParaExchangeCoupling[idx])
          );

          if (CheckPairSite(Def::ExchangeCoupling[idx][0], Def::ExchangeCoupling[idx][1], Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;

    case KWIsing:
      /*ising.def--------------------------------------*/
      if (Def::NIsingCoupling > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NIsingCoupling) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
            &isite1,
            &isite2,
            &dvalue_re
          );

          if (CheckPairSite(isite1, isite2, Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          //input into exchange couplings
          Def::HundCoupling[Def::NHundCoupling + idx][0] = isite1;
          Def::HundCoupling[Def::NHundCoupling + idx][1] = isite2;
          Def::ParaHundCoupling[Def::NHundCoupling + idx] = -dvalue_re / 2.0;
          //input into inter Coulomb
          Def::CoulombInter[Def::NCoulombInter + idx][0] = isite1;
          Def::CoulombInter[Def::NCoulombInter + idx][1] = isite2;
          Def::ParaCoulombInter[Def::NCoulombInter + idx] = -dvalue_re / 4.0;
          idx++;
        }
      }
      break;

    case KWPairLift:
      /*pairlift.def--------------------------------------*/
      if (Def::NPairLiftCoupling > 0) {
        if (Def::iCalcModel != DC::SpinGC) {
          fprintf(MP::STDOUT, "PairLift is active only in SpinGC.\n");
          return(-1);
        }
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL)
        {
          if (idx == Def::NPairLiftCoupling) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
            &(Def::PairLiftCoupling[idx][0]),
            &(Def::PairLiftCoupling[idx][1]),
            &(Def::ParaPairLiftCoupling[idx])
          );

          if (CheckPairSite(Def::PairLiftCoupling[idx][0], Def::PairLiftCoupling[idx][1], Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;

    case KWInterAll:
      /*interall.def---------------------------------------*/
      Def::NInterAll_Diagonal = 0;
      Def::NInterAll_OffDiagonal = 0;
      if (Def::NInterAll > 0) {
        icnt_interall = 0;
        icnt_diagonal = 0;
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NInterAll) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %d %d %d %d %d %d %d %lf %lf\n",
            &isite1,
            &isigma1,
            &isite2,
            &isigma2,
            &isite3,
            &isigma3,
            &isite4,
            &isigma4,
            &dvalue_re,
            &dvalue_im
          );

          if (CheckInterAllCondition(Def::iCalcModel, Def::Nsite, Def::iFlgGeneralSpin, Def::LocSpn,
            isite1, isigma1, isite2, isigma2,
            isite3, isigma3, isite4, isigma4) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          if (InputInterAllInfo(&icnt_interall,
            Def::InterAll,
            Def::ParaInterAll,
            isite1, isigma1,
            isite2, isigma2,
            isite3, isigma3,
            isite4, isigma4,
            dvalue_re, dvalue_im
          ) != 0) {
            icnt_diagonal += 1;
          }
          idx++;
        }
      }

      Def::NInterAll = icnt_interall;
      Def::NInterAll_Diagonal = icnt_diagonal;
      Def::NInterAll_OffDiagonal = Def::NInterAll - Def::NInterAll_Diagonal;

      if (GetDiagonalInterAll(
        Def::InterAll, Def::ParaInterAll, Def::NInterAll,
        Def::InterAll_Diagonal, Def::ParaInterAll_Diagonal,
        Def::InterAll_OffDiagonal, Def::ParaInterAll_OffDiagonal,
        Def::EDChemi, Def::EDSpinChemi, Def::EDParaChemi, &Def::EDNChemi,
        Def::iCalcModel
      ) != 0) {
        fclose(fp);
        return(-1);
      }

      if (CheckInterAllHermite(
        Def::InterAll, Def::ParaInterAll,
        Def::InterAll_OffDiagonal, Def::ParaInterAll_OffDiagonal,
        Def::NInterAll_OffDiagonal, Def::iCalcModel
      ) != 0) {
        fprintf(MP::STDOUT, "%s", "Error: NonHermite Pair exists in InterAll. \n");
        fclose(fp);
        return (-1);
      }
      break;

    case KWOneBodyG:
      /*cisajs.def----------------------------------------*/
      if (Def::NCisAjt > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NCisAjt) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %d %d %d\n",
            &isite1,
            &isigma1,
            &isite2,
            &isigma2);

          if (Def::iCalcModel == DC::Spin) {
            if (isite1 != isite2) {
              fprintf(MP::STDOUT, "Warning: Site component of (i, j) =(%d, %d) is ignored.\n",
                isite1, isite2);
              Def::NCisAjt--;
              continue;
            }
          }

          Def::CisAjt[idx][0] = isite1;
          Def::CisAjt[idx][1] = isigma1;
          Def::CisAjt[idx][2] = isite2;
          Def::CisAjt[idx][3] = isigma2;

          if (CheckPairSite(isite1, isite2, Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;

    case KWTwoBodyG:
      /*cisajscktaltdc.def--------------------------------*/
      if (Def::NCisAjtCkuAlvDC > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          if (idx == Def::NCisAjtCkuAlvDC) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %d %d %d %d %d %d\n",
            &isite1,
            &isigma1,
            &isite2,
            &isigma2,
            &isite3,
            &isigma3,
            &isite4,
            &isigma4
          );

          if (Def::iCalcModel == DC::Spin || Def::iCalcModel == DC::SpinGC) {
            if (CheckFormatForSpinInt(isite1, isite2, isite3, isite4) != 0) {
              wrapperMPI::Exit(-1);
              //Def::NCisAjtCkuAlvDC--;
              //continue;
            }
          }


          Def::CisAjtCkuAlvDC[idx][0] = isite1;
          Def::CisAjtCkuAlvDC[idx][1] = isigma1;
          Def::CisAjtCkuAlvDC[idx][2] = isite2;
          Def::CisAjtCkuAlvDC[idx][3] = isigma2;
          Def::CisAjtCkuAlvDC[idx][4] = isite3;
          Def::CisAjtCkuAlvDC[idx][5] = isigma3;
          Def::CisAjtCkuAlvDC[idx][6] = isite4;
          Def::CisAjtCkuAlvDC[idx][7] = isigma4;

          if (CheckQuadSite(isite1, isite2, isite3, isite4, Def::Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          idx++;
        }
      }
      break;

    case KWLaser:
      //printf("KWLaser\n");
      /*laser.def----------------------------------*/
      if (Def::NLaser > 0) {
        //printf("Read Start\n");
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%s %lf\n", ctmp, &(Def::ParaLaser[idx]));
          //printf("[%d]:%f\n",idx,Def::ParaLaser[idx]);
          idx++;
        }
        if (idx != Def::NLaser) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
      }
      break;

    case KWTEOneBody:
      if (Def::NTETimeSteps > 0) {
        idx = 0;
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%lf %d\n", &(Def::TETime[idx]), &(Def::NTETransfer[idx]));
          for (i = 0; i < Def::NTETransfer[idx]; ++i) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
            sscanf(ctmp2, "%d %d %d %d %lf %lf\n",
              &isite1,
              &isigma1,
              &isite2,
              &isigma2,
              &dvalue_re,
              &dvalue_im
            );
            Def::TETransfer[idx][i][0] = isite1;
            Def::TETransfer[idx][i][1] = isigma1;
            Def::TETransfer[idx][i][2] = isite2;
            Def::TETransfer[idx][i][3] = isigma2;
            Def::ParaTETransfer[idx][i] = std::complex<double>(dvalue_re, dvalue_im);
          }
          //check Transfer Hermite
          if (CheckTETransferHermite(Def::NTETransfer[idx], idx) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          idx++;
        }
        if (idx != Def::NTETimeSteps) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
      }
      break;

    case KWTETwoBody:
      if (Def::NTETimeSteps > 0) {
        idx = 0;
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%lf %d\n", &(Def::TETime[idx]), &(Def::NTEInterAll[idx]));
          icnt_interall = 0;
          icnt_diagonal = 0;
          for (i = 0; i < Def::NTEInterAll[idx]; ++i) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
            sscanf(ctmp2, "%d %d %d %d %d %d %d %d %lf %lf\n",
              &isite1,
              &isigma1,
              &isite2,
              &isigma2,
              &isite3,
              &isigma3,
              &isite4,
              &isigma4,
              &dvalue_re,
              &dvalue_im
            );
            if (CheckInterAllCondition(Def::iCalcModel, Def::Nsite, Def::iFlgGeneralSpin, Def::LocSpn,
              isite1, isigma1, isite2, isigma2,
              isite3, isigma3, isite4, isigma4) != 0) {
              fclose(fp);
              return ReadDefFileError(defname);
            }
            if (InputInterAllInfo(&icnt_interall,
              Def::TEInterAll[idx],
              Def::ParaTEInterAll[idx],
              isite1, isigma1,
              isite2, isigma2,
              isite3, isigma3,
              isite4, isigma4,
              dvalue_re, dvalue_im
            ) != 0) {
              icnt_diagonal += 1;
            }
          }

          Def::NTEInterAll[idx] = icnt_interall;
          Def::NTEInterAllDiagonal[idx] = icnt_diagonal;
          Def::NTEInterAllOffDiagonal[idx] = icnt_interall - icnt_diagonal;
          //Diagonal -> OffDiagonal -> search pair -> hermite
          if (GetDiagonalInterAll(Def::TEInterAll[idx], Def::ParaTEInterAll[idx], 
            Def::NTEInterAll[idx], Def::TEInterAllDiagonal[idx], Def::ParaTEInterAllDiagonal[idx],
            Def::TEInterAllOffDiagonal[idx], Def::ParaTEInterAllOffDiagonal[idx], Def::TEChemi[idx], 
            Def::SpinTEChemi[idx], Def::ParaTEChemi[idx], &Def::NTEChemi[idx], Def::iCalcModel) != 0)
          {
            fclose(fp);
            return (-1);
          }

          if (CheckInterAllHermite(
            Def::TEInterAll[idx], Def::ParaTEInterAll[idx],
            Def::TEInterAllOffDiagonal[idx], Def::ParaTEInterAllOffDiagonal[idx],
            Def::NTEInterAllOffDiagonal[idx], Def::iCalcModel
          ) != 0) {
            fprintf(MP::STDOUT, "%s", "Error: NonHermite Pair exists in InterAll. \n");
            fclose(fp);
            return (-1);
          }
          idx++;
        }

        if (idx != Def::NTETimeSteps) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
      }
      break;

    case KWBoost:
      /* boost.def--------------------------------*/
      //input magnetic field
      wrapperMPI::Fgets(ctmp2, 256, fp);
      sscanf(ctmp2, "%lf %lf %lf\n",
        &dArrayValue_re[0],
        &dArrayValue_re[1],
        &dArrayValue_re[2]);
      for (iline = 0; iline < 3; iline++) {
        Boost::vecB[iline] = dArrayValue_re[iline];
      }

      //this line is skipped;
      wrapperMPI::Fgets(ctmp2, 256, fp);

      //input arrayJ
      if (Boost::NumarrayJ > 0) {
        for (iline = 0; iline < Boost::NumarrayJ; iline++) {
          for (ilineIn = 0; ilineIn < 3; ilineIn++) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
            sscanf(ctmp2, "%lf %lf %lf\n",
              &dArrayValue_re[0],
              &dArrayValue_re[1],
              &dArrayValue_re[2]);
            for (ilineIn2 = 0; ilineIn2 < 3; ilineIn2++) {
              Boost::arrayJ[iline][ilineIn][ilineIn2] = dArrayValue_re[ilineIn2];
            }
          }
        }
      }

      //this line is skipped;
      wrapperMPI::Fgets(ctmp2, 256, fp);

      //read list_6spin_star
      if (Boost::num_pivot > 0) {
        for (iline = 0; iline < Boost::num_pivot; iline++) {
          //input
          wrapperMPI::Fgets(ctmp2, 256, fp);
          sscanf(ctmp2, "%d %d %d %d %d %d %d\n",
            &Boost::list_6spin_star[iline][0],
            &Boost::list_6spin_star[iline][1],
            &Boost::list_6spin_star[iline][2],
            &Boost::list_6spin_star[iline][3],
            &Boost::list_6spin_star[iline][4],
            &Boost::list_6spin_star[iline][5],
            &Boost::list_6spin_star[iline][6]
          );
          //copy
          for (iloop = 0; iloop < Boost::R0; iloop++) {
            for (itmp = 0; itmp < 7; itmp++) {
              Boost::list_6spin_star[iloop* Boost::num_pivot + iline][itmp]
                = Boost::list_6spin_star[iline][itmp];
            }
          }
        }
      }

      //read list_6spin_pair
      if (Boost::num_pivot > 0) {
        for (iline = 0; iline < Boost::num_pivot; iline++) {
          //input
          for (ilineIn2 = 0; ilineIn2 < Boost::list_6spin_star[iline][0]; ilineIn2++) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
            sscanf(ctmp2, "%d %d %d %d %d %d %d\n",
              &Boost::list_6spin_pair[iline][0][ilineIn2],
              &Boost::list_6spin_pair[iline][1][ilineIn2],
              &Boost::list_6spin_pair[iline][2][ilineIn2],
              &Boost::list_6spin_pair[iline][3][ilineIn2],
              &Boost::list_6spin_pair[iline][4][ilineIn2],
              &Boost::list_6spin_pair[iline][5][ilineIn2],
              &Boost::list_6spin_pair[iline][6][ilineIn2]
            );

            //copy
            for (iloop = 0; iloop < Boost::R0; iloop++) {
              for (itmp = 0; itmp < 7; itmp++) {
                Boost::list_6spin_pair[iloop* Boost::num_pivot + iline][itmp][ilineIn2]
                  = Boost::list_6spin_pair[iline][itmp][ilineIn2];
              }
            }
          }
        }

      }

      break;

    case KWSingleExcitation:
      /*singleexcitation.def----------------------------------------*/
      if (Def::NNSingleExcitationOperator > 0) {
        if (Def::iCalcModel == DC::Spin || Def::iCalcModel == DC::SpinGC) {
          fprintf(stderr, "SingleExcitation is not allowed for spin system.\n");
          fclose(fp);
          return ReadDefFileError(defname);
        }
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%d\n", &Def::NSingleExcitationOperator[idx]);
          Def::SingleExcitationOperator[idx] = (int**)malloc(sizeof(int*)*Def::NSingleExcitationOperator[idx]);
          Def::ParaSingleExcitationOperator[idx] = (std::complex<double>*)malloc(
            sizeof(std::complex<double>)*Def::NSingleExcitationOperator[idx]);
          for (i = 0; i < Def::NSingleExcitationOperator[idx]; ++i) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
            sscanf(ctmp2, "%d %d %d %lf %lf\n",
              &isite1,
              &isigma1,
              &itype,
              &dvalue_re,
              &dvalue_im
            );

            if (CheckSite(isite1, Def::Nsite) != 0) {
              fclose(fp);
              return ReadDefFileError(defname);
            }

            Def::SingleExcitationOperator[idx][i] = (int*)malloc(sizeof(int) * 3);
            Def::SingleExcitationOperator[idx][i][0] = isite1;
            Def::SingleExcitationOperator[idx][i][1] = isigma1;
            Def::SingleExcitationOperator[idx][i][2] = itype;
            Def::ParaSingleExcitationOperator[idx][i] 
              = std::complex<double>(dvalue_re,dvalue_im);
          }/*for (i = 0; i < Def::NSingleExcitationOperator[idx]; ++i)*/
          idx++;
        }
        if (idx != Def::NNSingleExcitationOperator) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
      }
      break;

    case KWPairExcitation:
      /*pairexcitation.def----------------------------------------*/
      if (Def::NNPairExcitationOperator > 0) {
        while (wrapperMPI::Fgets(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%d\n", &Def::NPairExcitationOperator[idx]);
          Def::PairExcitationOperator[idx] = (int**)malloc(sizeof(int*)*Def::NPairExcitationOperator[idx]);
          Def::ParaPairExcitationOperator[idx] = (std::complex<double>*)malloc(
            sizeof(std::complex<double>)*Def::NPairExcitationOperator[idx]);
          for (i = 0; i < Def::NPairExcitationOperator[idx]; ++i) {
            wrapperMPI::Fgets(ctmp2, 256, fp);
            sscanf(ctmp2, "%d %d %d %d %d %lf %lf\n",
              &isite1,
              &isigma1,
              &isite2,
              &isigma2,
              &itype,
              &dvalue_re,
              &dvalue_im
            );
            if (CheckPairSite(isite1, isite2, Def::Nsite) != 0) {
              fclose(fp);
              return ReadDefFileError(defname);
            }

            Def::PairExcitationOperator[idx][i] = (int*)malloc(sizeof(int) * 5);
            if (itype == 1) {
              Def::PairExcitationOperator[idx][i][0] = isite1;
              Def::PairExcitationOperator[idx][i][1] = isigma1;
              Def::PairExcitationOperator[idx][i][2] = isite2;
              Def::PairExcitationOperator[idx][i][3] = isigma2;
              Def::PairExcitationOperator[idx][i][4] = itype;
              Def::ParaPairExcitationOperator[idx][i] = 
                std::complex<double>(dvalue_re, dvalue_im);
            }
            else {
              Def::PairExcitationOperator[idx][i][0] = isite2;
              Def::PairExcitationOperator[idx][i][1] = isigma2;
              Def::PairExcitationOperator[idx][i][2] = isite1;
              Def::PairExcitationOperator[idx][i][3] = isigma1;
              Def::PairExcitationOperator[idx][i][4] = itype;
              Def::ParaPairExcitationOperator[idx][i] = 
                std::complex<double>(-dvalue_re, -dvalue_im);
            }
          }/*for (i = 0; i < Def::NPairExcitationOperator[idx]; ++i)*/
          idx++;
        }
        if (idx != Def::NNPairExcitationOperator) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
      }
      break;

    default:
      break;
    }
    fclose(fp);

    switch (iKWidx) {
    case KWCoulombIntra:
    case KWCoulombInter:
    case KWHund:
    case KWPairHop:
    case KWExchange:
    case KWIsing:
    case KWPairLift:
      if (Def::iFlgGeneralSpin == TRUE) {
        fprintf(MP::STDOUT, "%s", 
          "Error: Use only InterAll for setteing interactions for general spin.\n");
        return(-1);
      }
      break;
    default:
      break;
    }
  }

  ResetInteractionNum();
  /*=======================================================================*/
  return 0;
}
/**
 * @brief Check Site Number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckSite(
  const int iSite/**<[in] site number*/,
  const int iMaxNum/**<[in] Max site number*/
)
{
  if (iSite >= iMaxNum) return(-1);
  return 0;
}
/**
 * @brief Check Site Number for a pair -> (siteA, siteB).
 * @param[in] iSite1 a site number on a site A.
 * @param[in] iSite2 a site number on a site B.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckPairSite(
                  const int iSite1,
                  const int iSite2,
                  const int iMaxNum
                  )
{
  if(CheckSite(iSite1, iMaxNum)!=0){
    return(-1);
  }
  if(CheckSite(iSite2, iMaxNum)!=0){
    return(-1);
  }
  return 0;
}

/**
 * @brief Check Site Number for a quad -> (siteA, siteB, siteC, siteD).
 * @param[in] iSite1 a site number on site A.
 * @param[in] iSite2 a site number on site B.
 * @param[in] iSite3 a site number on site C.
 * @param[in] iSite4 a site number on site D.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckQuadSite(
                  const int iSite1,
                  const int iSite2,
                  const int iSite3,
                  const int iSite4,
                  const int iMaxNum
                  )
{
  if(CheckPairSite(iSite1, iSite2, iMaxNum)!=0){
    return(-1);
  }
  if(CheckPairSite(iSite3, iSite4, iMaxNum)!=0){
    return(-1);
  }
  return 0;
}

/**
 * @brief Check Hermite for Transfer integrals.
 * @param[in] X Define List for getting transfer integrals.
 * @retval 0 Hermite.
 * @retval -1 NonHermite.
 * @version 0.2
 * @details rearray a GeneralTransfer array to satisfy a condition of hermite conjugation between 2*i and 2*i+1 components.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckTransferHermite
(
  
)
{
  int i, j;
  int isite1, isite2;
  int isigma1, isigma2;
  int itmpsite1, itmpsite2;
  int itmpsigma1, itmpsigma2;
  int itmperrsite1, itmperrsite2;
  int itmperrsigma1, itmperrsigma2;
  std::complex<double> dcerrTrans;
  int icheckHermiteCount = FALSE;
  int iCount = 0;

  std::complex<double> ddiff_trans;
  int itmpIdx, icntHermite, icntchemi;
  icntHermite = 0;
  icntchemi = 0;

  for (i = 0; i < Def::NTransfer; i++) {
    isite1 = Def::GeneralTransfer[i][0];
    isigma1 = Def::GeneralTransfer[i][1];
    isite2 = Def::GeneralTransfer[i][2];
    isigma2 = Def::GeneralTransfer[i][3];
    icheckHermiteCount = FALSE;
    for (j = 0; j < Def::NTransfer; j++) {
      itmpsite1 = Def::GeneralTransfer[j][0];
      itmpsigma1 = Def::GeneralTransfer[j][1];
      itmpsite2 = Def::GeneralTransfer[j][2];
      itmpsigma2 = Def::GeneralTransfer[j][3];
      if (isite1 == itmpsite2 && isite2 == itmpsite1) {
        if (isigma1 == itmpsigma2 && isigma2 == itmpsigma1) {

          ddiff_trans = Def::ParaGeneralTransfer[i] - conj(Def::ParaGeneralTransfer[j]);
          if (abs(ddiff_trans) > eps_CheckImag0) {
            itmperrsite1 = itmpsite1;
            itmperrsigma1 = itmpsigma1;
            itmperrsite2 = itmpsite2;
            itmperrsigma2 = itmpsigma2;
            dcerrTrans = Def::ParaGeneralTransfer[j];
            fprintf(MP::STDOUT, "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n",
              isite1, isigma1, isite2, isigma2, real(Def::ParaGeneralTransfer[i]), imag(Def::ParaGeneralTransfer[i]));
            fprintf(MP::STDOUT, "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n", 
              itmperrsite1, itmperrsigma1, itmperrsite2, itmperrsigma2, real(dcerrTrans), imag(dcerrTrans));
            iCount++;
          }
          else {
            if (icheckHermiteCount == FALSE) {
              if (i <= j) {
                if (2 * icntHermite >= Def::NTransfer) {
                  fprintf(stderr, "Elements of Transfers are incorrect.\n");
                  return(-1);
                }
                if (isite1 != isite2 || isigma1 != isigma2) {
                  for (itmpIdx = 0; itmpIdx < 4; itmpIdx++) {
                    Def::EDGeneralTransfer[2 * icntHermite][itmpIdx] = Def::GeneralTransfer[i][itmpIdx];
                    Def::EDGeneralTransfer[2 * icntHermite + 1][itmpIdx] = Def::GeneralTransfer[j][itmpIdx];
                  }
                  Def::EDParaGeneralTransfer[2 * icntHermite] = Def::ParaGeneralTransfer[i];
                  Def::EDParaGeneralTransfer[2 * icntHermite + 1] = Def::ParaGeneralTransfer[j];
                  icntHermite++;
                }
                else {
                  Def::EDChemi[icntchemi] = Def::GeneralTransfer[i][0];
                  Def::EDSpinChemi[icntchemi] = Def::GeneralTransfer[i][1];
                  Def::EDParaChemi[icntchemi] = real(Def::ParaGeneralTransfer[i]);
                  icntchemi += 1;
                }
              }
              icheckHermiteCount = TRUE;
            }
          }
        }

      }
    }

    //if counterpart for satisfying hermite conjugate does not exist.
    if (icheckHermiteCount == FALSE) {
      fprintf(MP::STDOUT, "Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n", isite1, isigma1, isite2, isigma2, real(Def::ParaGeneralTransfer[i]), imag(Def::ParaGeneralTransfer[i]));
      iCount++;
    }
  }

  if (iCount != 0) {
    return -1;
  }
  Def::EDNTransfer = 2 * icntHermite;
  Def::EDNChemi = icntchemi;

  //To realize ido-san's result
  for (i = 0; i < Def::EDNTransfer; i++) {
    for (itmpIdx = 0; itmpIdx < 4; itmpIdx++) {
      Def::GeneralTransfer[i][itmpIdx] = Def::EDGeneralTransfer[i][itmpIdx];
    }
    Def::ParaGeneralTransfer[i] = Def::EDParaGeneralTransfer[i];
  }
  return 0;
}
/** 
 * @brief function of judging a type of define files.
 * 
 * @param[in] argc argument count
 * @param[in] argv argument vector 
 * @param[out] mode a number to show a type of a define file
 * 
 * @retval 0 format is correct
 * @retval -1 format is incorrect
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void JudgeDefType
(
 const int argc,
 char *argv[],
 int *mode
 )
{
   int ver_maj =
#include "version_major.hpp"
;
   int ver_min =
#include "version_miner.hpp"
;
   int ver_pat =
#include "version_patch.hpp"
;

  if(argc == 3 && 
     (CheckWords(argv[1], "-e") == 0 ||
      CheckWords(argv[1], "--Expert") == 0)){
    *mode=EXPERT_MODE;
  }
  else if (argc == 3 && 
           (CheckWords(argv[1], "-s") ==0 ||
            CheckWords(argv[1], "--Standard") == 0 )){
    *mode=STANDARD_MODE;
  }
  else if (argc == 3 && 
           (CheckWords(argv[1], "-sdry") == 0 ||
            CheckWords(argv[1], "-s-dry") == 0)
           ){
    *mode = STANDARD_DRY_MODE;
  }
  else if (argc >= 2 &&
           (CheckWords(argv[1], "-v") == 0
            || CheckWords(argv[1], "--version") == 0)
           ) {
    fprintf(MP::STDOUT, "\nHPhi version %d.%d.%d \n\n", ver_maj, ver_min, ver_pat);
    exit(-1);
  }
  else{
    fprintf(MP::STDOUT, "\n[Usage] \n");
    fprintf(MP::STDOUT, "* Expert mode \n");
    fprintf(MP::STDOUT, "   $ HPhi -e {namelist_file} \n");
    fprintf(MP::STDOUT, "* Standard mode \n");
    fprintf(MP::STDOUT, "   $ HPhi -s {input_file} \n");
    fprintf(MP::STDOUT, "* Standard DRY mode \n");
    fprintf(MP::STDOUT, "   $ HPhi -sdry {input_file} \n");
    fprintf(MP::STDOUT, "   In this mode, Hphi stops after it generats expert input files. \n");
    fprintf(MP::STDOUT, "* Print the version \n");
    fprintf(MP::STDOUT, "   $ HPhi -v \n\n");
    exit(-1);
  }
}

/** 
 * @brief function of checking format of spin interactions
 * 
 * @param[in] site1 a site number on site1.
 * @param[in] site2 a site number on site2.
 * @param[in] site3 a site number on site3.
 * @param[in] site4 a site number on site4.
 * 
 * @retval 0 format is correct
 * @retval -1 format is incorrect
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckFormatForSpinInt
(
 const int site1,
 const int site2,
 const int site3,
 const int site4
 ){
  if(site1==site2 && site3==site4){
    return 0;
  }

  fprintf(MP::STDOUT,
          "Warning: Site component of (i, j, k, l) =(%d, %d, %d, %d) is not correct; i=j and k=l must be satisfied. \n",
          site1, site2, site3, site4);
  return(-1);

}


/// \brief function of checking format of Kondo interactions
/// \param isite1 a site number on site1
/// \param isite2 a site number on site2
/// \param isite3 a site number on site3
/// \param isite4 a site number on site4
/// \param iLocInfo An array with the value of S at each site.
/// \retval  0 format is correct
/// \retval  -1 format is incorrect
/// \version 0.1
/// \author Takahiro Misawa (The University of Tokyo)
/// \author Kazuyoshi Yoshimi (The University of Tokyo)
int CheckFormatForKondoInt
        (
                const int isite1, const int isite2,
                const int isite3, const int isite4,
                int* iLocInfo
        )
{
  if (iLocInfo[isite1] != ITINERANT || iLocInfo[isite2] != ITINERANT) {
    if (isite1 != isite2) {
      fprintf(MP::STDOUT, "Error: Site component of (i, j, k, l) =(%d, %d, %d, %d) is incorrect.\n", isite1, isite2, isite3, isite4);
      return -1;
    }
  }
  if (iLocInfo[isite3] != ITINERANT || iLocInfo[isite4] != ITINERANT) {
    if (isite3 != isite4) {
      fprintf(MP::STDOUT, "Error: Site component of (i, j, k, l) =(%d, %d, %d, %d) is incorrect.\n", isite1, isite2, isite3, isite4);
      return -1;
    }
  }
  return 0;
}

/** 
 * @brief function to set convergence factors
 * 
 * @param[in] X Define list to get Lanczos eps.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void SetConvergenceFactor()
{
  //In future, convergence facator can be set by a def file.
  Def::eps_Lanczos     = pow(10,-Def::LanczosEps);
}
/** 
 * @brief function of checking indexies of localized spin
 * 
 * @param [inout] X Define list to get and put information of localized spin
 * 
 * @return TURE Indecies of localizes spin is correct
 * @return FALSE Indecies of localizes spin is incorrect
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
int CheckLocSpin()
{

  int i=0;
  switch(Def::iCalcModel){
  case DC::Hubbard:
  case DC::HubbardNConserved:
  case DC::HubbardGC:
  case DC::SpinlessFermion:
  case DC::SpinlessFermionGC:
    for(i=0; i<Def::Nsite; i++){
      if(Def::LocSpn[i]!=ITINERANT){
        return FALSE;
      }
    }
    break;

  case DC::Kondo:
  case DC::KondoGC:
    for(i=0; i<Def::Nsite; i++){
      if(Def::LocSpn[i]>LOCSPIN){
        Def::iFlgGeneralSpin=TRUE;
      }
      else if(Def::LocSpn[i]<ITINERANT){
        return FALSE;
      }
    }
    break;

  case DC::Spin:
  case DC::SpinGC:
    for(i=0; i<Def::Nsite; i++){
      if(Def::LocSpn[i]>LOCSPIN){
        Def::iFlgGeneralSpin=TRUE;
      }
      else if(Def::LocSpn[i]<LOCSPIN){
        return FALSE;
      }
    }
    break;
  
  default:
    return FALSE;
    //break;
  }

  if(CheckTotal2Sz() != TRUE){
    return FALSE;
  }
  return TRUE;
}  
/** 
 * 
 * @brief function of resetting number of interactions
 * 
 * @param[out] X Define list to add number of ising coulomnb interactions
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
void ResetInteractionNum()
{
  Def::NHundCoupling += Def::NIsingCoupling;
  Def::NCoulombInter += Def::NIsingCoupling;
}

/** 
 * @brief function of initializing interactions
 * 
 * @param[out] X Define list to initialize number of interactions
 * @version 0.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
void InitializeInteractionNum()
{
  Def::NTransfer=0;
  Def::NCoulombIntra=0;
  Def::NCoulombInter=0;
  Def::NHundCoupling = 0;
  Def::NExchangeCoupling = 0;
  Def::NPairHopping = 0;
  Def::NIsingCoupling=0;
  Def::NPairLiftCoupling=0;
  Def::NInterAll=0;
  Def::NInterAll_Diagonal = 0;
  Def::NInterAll_OffDiagonal = 0;
  Def::NCisAjt=0;
  Def::NCisAjtCkuAlvDC=0;
  Def::NNSingleExcitationOperator=0;
  Def::NNPairExcitationOperator=0;
  //[s] Time Evolution
  Def::NTETimeSteps=0;
  Def::NLaser=0;
  Def::NTEInterAll=0;
  Def::NTETransfer=0;
  //[e] Time Evolution

}


///
/// \brief function of checking spin index for all interactions
/// \param isite1 a site number on site1
/// \param isigma1 a spin index on site1
/// \param isite2 a site number on site2
/// \param isigma2 a spin index on site2
/// \param isite3 a site number on site3
/// \param isigma3 a spin index on site3
/// \param isite4 a site number on site4
/// \param isigma4 a spin index on site4
/// \param iLocInfo An array with the value of S at each site.
/// \retval  TRUE spin index is correct
/// \retval  FALSE spin index is incorrect
/// \version 0.2
/// \author Kazuyoshi Yoshimi (The University of Tokyo)
/// \author Takahiro Misawa (The University of Tokyo)
int CheckGeneralSpinIndexForInterAll
(
        const int isite1, const int isigma1,
        const int isite2, const int isigma2,
        const int isite3, const int isigma3,
        const int isite4, const int isigma4,
        int* iLocInfo
 )
{
   if( isigma1 > iLocInfo[isite1] || isigma2 >iLocInfo[isite2]
         ||isigma3 > iLocInfo[isite3] || isigma4 >iLocInfo[isite4]){
        fprintf(MP::STDOUT, "%s", "Error: Spin index is incorrect for interactions defined in InterAll file.\n");
        return FALSE;
    }
  return TRUE;
}

/** 
 * @brief function of checking spin index for transfers
 * 
 * @param[in] X Define list to get informations of transfers
 * @retval TRUE spin index is correct
 * @retval FALSE spin index is incorrect
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
int CheckSpinIndexForTrans()
{
  int i = 0;
  int isite1, isite2;
  int isigma1, isigma2;
  if (Def::iFlgGeneralSpin == TRUE) {
    for (i = 0; i < Def::NTransfer; i++) {
      isite1 = Def::GeneralTransfer[i][0];
      isigma1 = Def::GeneralTransfer[i][1];
      isite2 = Def::GeneralTransfer[i][2];
      isigma2 = Def::GeneralTransfer[i][3];
      if (isigma1 > Def::LocSpn[isite1] || isigma2 > Def::LocSpn[isite2]) {
        fprintf(MP::STDOUT, "%s", "Error: Spin index is incorrect for transfers defined in Trans file.\n");
        return FALSE;
      }
    }
  }
  return TRUE;
}
/** 
 * @brief function of checking an input data of total2Sz
 * 
 * @param[in] X Define list to get informations of transfers
 * @retval TRUE spin index is correct
 * @retval FALSE spin index is incorrect
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
int CheckTotal2Sz()
{
  if(Def::iFlgSzConserved==TRUE && Def::iFlgGeneralSpin==FALSE){
    int tmp_Nup=Def::NLocSpn+Def::NCond+Def::Total2Sz;
    int tmp_Ndown=Def::NLocSpn+Def::NCond-Def::Total2Sz;
    if(tmp_Nup%2 != 0 && tmp_Ndown%2 !=0){
      printf("Nup=%d, Ndown=%d\n",Def::Nup,Def::Ndown);
      fprintf(MP::STDOUT, "2Sz is incorrect.\n");
      return FALSE;
    }
  }
  return TRUE;
}

/**
 *
 * @brief function of checking whether ctmp is same as cKeyWord or not
 *
 * @param[in] ctmp A word to be checked whether it matches the registerd keyword or not.
 * @param[in] cKeyWord Registered keyword name
 * @return 0 ctmp is same as cKeyWord
 *
 * @version 1.1.0
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckWords(
  const char* ctmp,
  const char* cKeyWord
)
{
  int i = 0;
  char ctmp_small[256] = { 0 };
  char cKW_small[256] = { 0 };
  int NcKeyWord, Nctmp;
  NcKeyWord = strlen(cKeyWord);
  strncpy(cKW_small, cKeyWord, NcKeyWord);

  for (i = 0; i < NcKeyWord; i++) {
    cKW_small[i] = tolower(cKW_small[i]);
  }
  Nctmp = strlen(ctmp);
  strncpy(ctmp_small, ctmp, Nctmp);
  for (i = 0; i < Nctmp; i++) {
    ctmp_small[i] = tolower(ctmp_small[i]);
  }
  if (Nctmp < NcKeyWord) Nctmp = NcKeyWord;
  return(strncmp(ctmp_small, cKW_small, Nctmp));
}
///
/// \brief function of getting file name labeled by the keyword
/// \param iKWidx index of keyword
/// \param FileName filename
/// \retval 0 normally finished getting file name.
/// \retval -1 unnormally finished getting file name.
int GetFileNameByKW(
        int iKWidx,
        char **FileName
){
  if(cFileNameListFile == NULL){
    return -1;
  }
  *FileName=cFileNameListFile[iKWidx];
  return 0;
}


/**
 * @brief Check InterAll condition.
 * @param[in] iCalcModel Target Model defined in CalcMod file (ex. Spin, SpinGC etc.).
 * @param[in] Nsite  A total number of site.
 * @param[in] iFlgGeneralSpin  Flag for general spin (TRUE: General Spin, FALSE: Spin-1/2).
 * @param[in] iLocInfo An array with the value of S at each site
 * @param[in] isite1 a site number on the site A.
 * @param[in] isigma1 a spin index on the site A.
 * @param[in] isite2 a site number on the site B.
 * @param[in] isigma2 a spin index on the site B.
 * @param[in] isite3 a site number on the site C.
 * @param[in] isigma3 a spin index on the site C.
 * @param[in] isite4 a site number on the site D.
 * @param[in] isigma4 a spin index on the site D.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 2.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckInterAllCondition(
        int iCalcModel,
        int Nsite,
        int iFlgGeneralSpin,
        int *iLocInfo,
        int isite1, int isigma1,
        int isite2, int isigma2,
        int isite3, int isigma3,
        int isite4, int isigma4
){
  if(CheckQuadSite(isite1, isite2, isite3, isite4, Nsite) !=0){
    fprintf(stderr, "%s", "Error: Site index of InterAll is incorrect.\n");
    return(-1);
  }

  if(iCalcModel == DC::Spin || iCalcModel == DC::SpinGC){
    if(CheckFormatForSpinInt(isite1, isite2, isite3, isite4)!=0){
      fprintf(stderr, "%s", "Error: Spin index of InterAll is incorrect.\n");
      return(-1);
    }
  }
  else if(iCalcModel == DC::SpinlessFermion || iCalcModel== DC::SpinlessFermionGC){
    if(isigma1 !=0 || isigma2 != 0 || isigma3 != 0 || isigma4 !=0){
      fprintf(stderr, "%s", "Error: Spin index of InterAll is incorrect.\n");
      return -1;
    }
  }
  else if(iCalcModel == DC::Kondo){
    if(CheckFormatForKondoInt(isite1, isite2, isite3, isite4, iLocInfo)!=0){
      return -1;
    }
  }

  if(iFlgGeneralSpin ==TRUE) {
    if(CheckGeneralSpinIndexForInterAll(isite1, isigma1, isite2, isigma2, isite3, isigma3, isite4, isigma4, iLocInfo)!=TRUE){
      return -1;
    }
  }
  return 0;
}
/**
\brief Input InterAll Interactions (Operators of the same kinds are grouped together).
\return 0 Interaction is off-diagonal
\return 1 Interaction is diagonal
*/
int InputInterAllInfo(
  int* icnt_interall,/**<total number of interall interactions*/
  int** iInterAllInfo,/**<arrays of information of interall interactions*/
  std::complex<double>* cInterAllValue,/**<arrays of values of interall interactions*/
  int isite1,/**<[in] site number on the site 1*/
  int isigma1,/**<[in] spin index on the site 1*/
  int isite2,/**<[in] site number on the site 2*/
  int isigma2,/**<[in] spin index on the site 2*/
  int isite3,/**<[in] site number on the site 3*/
  int isigma3,/**<[in] spin index on the site 3*/
  int isite4, /**<[in] site number on the site 4*/
  int isigma4,/**<[in] spin index on the site 4*/
  double dvalue_re, /**<[in]*/
  double dvalue_im/**<[in]*/
) {
  int i = 0;
  int iflg_interall = 0;
  //Collect and sum same components of InterAll interactions
  for (i = 0; i < *icnt_interall; i++) {
    if (isite1 == iInterAllInfo[i][0] && isite2 == iInterAllInfo[i][2] &&
        isite3 == iInterAllInfo[i][4] && isite4 == iInterAllInfo[i][6] &&
        isigma1 == iInterAllInfo[i][1] && isigma2 == iInterAllInfo[i][3] &&
        isigma3 == iInterAllInfo[i][5] && isigma4 == iInterAllInfo[i][7]) {
      cInterAllValue[i] += std::complex<double>(dvalue_re, dvalue_im);
      iflg_interall = 1;
      return 0;
    }
  }

  //Input all InterAll interactions
  if (iflg_interall == 0) {
    iInterAllInfo[*icnt_interall][0] = isite1;
    iInterAllInfo[*icnt_interall][1] = isigma1;
    iInterAllInfo[*icnt_interall][2] = isite2;
    iInterAllInfo[*icnt_interall][3] = isigma2;
    iInterAllInfo[*icnt_interall][4] = isite3;
    iInterAllInfo[*icnt_interall][5] = isigma3;
    iInterAllInfo[*icnt_interall][6] = isite4;
    iInterAllInfo[*icnt_interall][7] = isigma4;
    cInterAllValue[*icnt_interall] = std::complex<double>(dvalue_re, dvalue_im);
    *icnt_interall+=1;
    //Check Diagonal part or not
    if (isite1 == isite2 && isite3 == isite4 &&
        isigma1 == isigma2 && isigma3 == isigma4) { //normal diagonal part
      return 1;
    } else if (isite1 == isite4 && isite2 == isite3 &&
               isigma1 == isigma4 && isigma2 == isigma3) { //hund term
      return 1;
    }
  }
  return 0;
}
/**
@page page_addexpert Add new input-file for Expert mode
When you add a new input file to _namelist file_,
the following procedures must be needed.
In the following, we add the keyword "Test" as an example.

1. Add a new keyword to the end of @c cKWListOfFileNameList in @c readdef.cpp.
```
static char cKWListOfFileNameList[][D_CharTmpReadDef]
={
  "CalcMod",
  "ModPara",
  "LocSpin",
  "Trans",
  "CoulombIntra",
  "CoulombInter",
  "Hund",
  "PairHop",
  "Exchange",
  "InterAll",
  "OneBodyG",
  "TwoBodyG",
  "PairLift",
  "Ising",
  "Boost",
  "SingleExcitation",
  "PairExcitation",
  "SpectrumVec",
  "Laser",
  "TEOneBody",
  "TETwoBody"
  "Test"
}
```
Here, `` D_CharTmpReadDef `` is set as `` 200 `` in readdef.h.
If the the character number of added keyword exceeds `` 200 ``, please change the value.

2. Define the index of keyword (such as KWCalcMod, KWModPara...) in @c readdef.h.
```
#define KWCalcMod 0
#define KWModPara 1
#define KWLocSpin 2
#define KWTrans 3
#define KWCoulombIntra 4
#define KWCoulombInter 5
#define KWHund 6
#define KWPairHop 7
#define KWExchange 8
#define KWInterAll 9
#define KWOneBodyG 10
#define KWTwoBodyG 11
#define KWPairLift 12
#define KWIsing 13
#define KWBoost 14
#define KWSingleExcitation 15
#define KWPairExcitation 16
#define KWSpectrumVec 17
#define KWLaser 18
#define KWTEOneBody 19
#define KWTETwoBody 20
#define KWTest 21
```
The defined value must be same as the index of cKWListOfFileNameList
to get the name of keyword, i.e. cKWListOfFileNameList[KWTest] = "Test".

3. Add procedure of reading the file in @c ReadDefFileNInt function in @c readdef.cpp.
 ```
 for(iKWidx=0; iKWidx< D_iKWNumDef; iKWidx++) {
    strcpy(defname, cFileNameListFile[iKWidx]);

    if (strcmp(defname, "") == 0) continue;
    if(iKWidx==KWSpectrumVec){
      continue;
    }
    fprintf(MP::STDOUT, "  Read File %s for %s.\n", defname, cKWListOfFileNameList[iKWidx]);
    fp = wrapperMPI::Fopen(defname, "r");
    if (fp == NULL) return ReadDefFileError(defname);
    switch (iKWidx) {
    case KWTest:
        wrapperMPI::Fgets(...); //Add the procedure to read-line here.
    }
 ```
@sa ReadDefFileNInt

4. Use @c InitializeInteractionNum function to initialize variables.
 ```
 void InitializeInteractionNum
(
 
 )
{
  Def::NTransfer=0;
  Def::NCoulombIntra=0;
  Def::NCoulombInter=0;
  Def::NIsingCoupling=0;
  Def::NPairLiftCoupling=0;
  Def::NInterAll=0;
  Def::NCisAjt=0;
  Def::NCisAjtCkuAlvDC=0;
  Def::NSingleExcitationOperator=0;
  Def::NPairExcitationOperator=0;
  //[s] Time Evolution
  Def::NTETimeSteps=0;
  Def::NLaser=0;
  Def::NTEInterAll=0;
  Def::NTETransfer=0;
  //[e] Time Evolution
  Def::NTest = 0;
}
 ```
    @sa  InitializeInteractionNum
5. The memories of arrays are stored by xsetmem::def function in @c xsetmem.cpp.
   @sa  xsetmem::def
 **/

/**
@page page_addmodpara Add new parameter into modpara

You can set a value of parameters with a new keyword in ``modpara`` file by following way.

- Define a new variable corresponding to the above parameter in @c global.hpp file.

- The value with the keyword are read by `` ReadDefFileNInt `` function in @c readdef.cpp.

  In the following, we describe the detail of the flow of reading the parameter.

  To read the parameter, the switch statement where ``iKWidx`` matches ``KWModPara`` is used.
  The detail of the reading flow in this function are described as follows.

1. The first eight lines are header (not touch!).
  ```
        //! Read Header (5 lines).
       wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //1
       wrapperMPI::Fgets(ctmp2, 256, fp);
       sscanf(ctmp2, "%s %d\n", ctmp, &itmp); //2
       wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //3
       wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //4
       wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //5
       //! Read header name for files about data
       wrapperMPI::Fgets(ctmp2, 256, fp);
       sscanf(ctmp2, "%s %s\n", ctmp, Def::CDataFileHead); //6
        //! Read header name for files about parameters
       wrapperMPI::Fgets(ctmp2, 256, fp);
       sscanf(ctmp2, "%s %s\n", ctmp, Def::CParaFileHead); //7
       //! Read header (1 line).
       wrapperMPI::Fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);   //8
  ```

2. Each line is read by ``` wrapperMPI::Fgets(ctmp2, 256, fp) ``` function.

3. The line is divided into keyword and number by using ``CheckWords`` function.

   For example, when you add new key word "NTest", you can get the value as follows:
   ```
        if (CheckWords(ctmp, "NTest") == 0) {
                Def::NTest = (int) dtmp;
              }
   ```

@sa ReadDefFileNInt, CheckWords
*/

/**
@page page_addcalcmod Add new calculation mode into calcmod

You can set a new calculation mode with a new keyword in ``calcmod`` file by following way.

- Define a new variable corresponding to the new calculation mode in @c global.hpp file.

- The value with the keyword are read by `` ReadcalcmodFile `` function in @c readdef.cpp.

In the following, we describe the detail of the flow of setting the calculation mode.

1. Set initial value at the beginning of ReadcalcmodFile function.
  ```
  Def::iCalcType=0;
  Def::iFlgFiniteTemperature=0;
  Def::iCalcModel=0;
  Def::iOutputMode=0;
  Def::iCalcEigenVec=0;
  Def::iInitialVecType=0;
  Def::iOutputEigenVec=0;
  Def::iInputEigenVec=0;
  Def::iOutputHam=0;
  Def::iInputHam=0;
  Def::iFlgCalcSpec=0;
  Def::iReStart=0;
  Def::iFlgMPI=0;
  ```

2. Each line is read by ``` wrapperMPI::Fgets ``` function. 
   ``GetKWWithIdx`` function reads ctmp = keyword, itmp=index.
  ```
   while( wrapperMPI::Fgets(ctmpLine, D_CharTmpReadDef+D_CharKWDMAfp)!=NULL ){
    if( (iret=GetKWWithIdx(ctmpLine, ctmp, &itmp)) !=0){
      if(iret==1) continue;
      return(-1);
    }   
    if(CheckWords(ctmp, "CalcType")==0){
      Def::iCalcType=itmp;
    }
   ...
   }
  ```


3. The line is divided into keyword and number by using ``CheckWords`` function.

   For example, when you add new key word "NTest", you can get the value as follows:
   ```
        if (CheckWords(ctmp, "NTest") == 0) {
                Def::NTest = (int) dtmp;
              }
   ```

@sa ReadcalcmodFile, CheckWords
*/
