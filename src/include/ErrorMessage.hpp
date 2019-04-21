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

#ifndef HPHI_ERRORMESSAGE_H
#define HPHI_ERRORMESSAGE_H

//! Error Message in HPhiMain.c
extern int iErrCodeMem;

extern const char *cErrOutput;
extern const char *cErrDefFile;
extern const char *cErrnvec;

//! Error Message in HPhiTrans.c
extern const char *cErrLargeMem;


//! Error Message in readdef.c
extern const char *cErrReadDefFile;
extern const char *cErrDefFileFormat;
extern const char *cErrNLoc;
extern const char *cErrDefFileParam;
extern const char *cErrCalcType;
extern const char *cErrOutputMode;
extern const char *cErrCalcModel;
extern const char *cErrCalcEigenVec;
extern const char *cErrSetIniVec;
extern const char *cErrOutputHam;
extern const char *cErrInputHam;
extern const char *cErrRestart;
extern const char *cErrCUDA;
extern const char *cErrScaLAPACK;

extern const char *cErrKW;
extern const char *cErrKW_ShowList;
extern const char *cErrKW_Same;
extern const char *cErrKW_InCorPair;

extern const char *cErrNsite;
extern const char *cErrNcond;
extern const char *cErrNumAve;
extern const char *cErrExpecInterval;
extern const char *cErrLanczos_max;
extern const char *cErrLanczos_eps;
extern const char *cErrLanczosTarget;
extern const char *cErrLanczosExct;

extern const char *cErrMakeDef;
extern const char *cErrIncorrectDef;
extern const char *cErrNonHermiteTrans;
extern const char *cErrNonHermiteTransForAll;
extern const char *cErrNonHermiteInterAll;
extern const char *cErrNonConservedInterAll;
extern const char *cErrNonHermiteInterAllForAll;
extern const char *cErrIncorrectFormatForKondoInt;
extern const char *cErrIncorrectFormatForKondoTrans;
extern const char *cErrIncorrectFormatInter;
extern const char *cErrIncorrectSpinIndexForInter;
extern const char *cErrIncorrectSpinIndexForTrans;

extern const char *cErrIncorrectFormatForSpinTrans;
extern const char *cWarningIncorrectFormatForSpin2;


//! Error Message in CheckMPI.c
extern const char *cErrNProcNumberHubbard;
extern const char *cErrNProcNumberSpin;
extern const char *cErrNProcNumberGneralSpin;
extern const char *cErrNProcNumber;
extern const char *cErrNProcNumberSet;

//! Error Message in diagonal calc.c
extern const char *cErrNoModel;
extern const char *cErrNoHilbertSpace;

//! Error Message in bitcalc.c
extern const char *cErrSiteNumber;

//! Error Message in mltiply.c


//! Error Message in FileIO.c
extern const char *cErrFIOpen;

//! Error Message in sz.c
extern const char* cErrSz;
extern const char* cErrSz_NoFile;
extern const char* cErrSz_NoFile_Show;

#endif /* HPHI_ERRORMESSAGE_H */
