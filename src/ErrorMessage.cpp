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

#include "ErrorMessage.hpp"

//! Error Message in HPhiMain.c
int iErrCodeMem=1;

const char *cErrOutput="Error: Fail to make output folder in current directory. \n";
const char *cErrDefFile="Error: Definition files(*.def) are incomplete.\n";
const char *cErrnvec="Error: nvec should be smaller than exct are incorrect.\n";

//! Error Message in HPhiTrans.c
const char *cErrLargeMem="Error: Fail for memory allocation. see error code %d in user manual. \n";


//! Error Message in readdef.c
const char *cErrReadDefFile="Error: %s (Broken file or Not exist)\n";
const char *cErrDefFileFormat="Error: incorrect format= %s. \n";
const char *cErrNLoc ="Error: Ne=Nup+Ndown must be (Ne >= NLocalSpin).\n";
const char *cErrDefFileParam="Error: In %s, wrong parameter name:%s \n";
const char *cErrCalcType="Error in %s\n CalcType: 0: Lanczos Method, 1: Thermal Pure Quantum State Method, 2: Full Diagonalization Method, 3: Calculation Spectrum mode.\n";
const char *cErrOutputMode="Error in %s\n OutputMode: \n 0: calc one body green function and two body green functions,\n 1: calc one body green function and two body green functions and correlatinos for charge and spin.\n";
const char *cErrCalcEigenVec="Error in %s\n CalcEigenVec: \n 0: Lanczos+CG method,\n 1: Lanczos method.\n";
const char *cErrOutputHam="Error in %s\n OutputHam: \n 0: not output Hamiltonian,\n 1: output Hamiltonian.\n";
const char *cErrInputHam="Error in %s\n InputHam: \n 0: not input Hamiltonian,\n 1: input Hamiltonian.\n";
const char *cErrCalcModel="Error in %s\n CalcModel: \n 0: Hubbard, 1: Spin, 2: Kondo, 3: HubbardGC, 4: SpinGC, 5:KondoGC.\n";
const char *cErrSetIniVec="Error in %s\n InitialVecType: \n 0: complex type,\n 1: real type.\n";
const char *cErrRestart="Error in %s\n Restart: \n 0: not restart (default).\n 1: output a restart vector.\n 2: input a restart vector and output a new restart vector.\n 3: input a restart vector.\n";
const char *cErrCUDA="Error in %s\n NGPU: NGPU must be greater than 0.\n";
const char *cErrScaLAPACK="Error in %s\n ScaLAPACK: \n 0: Use LAPACK for FullDiag mode,\n 1: Use ScaLAPACK for FullDiag mode.\n";

const char *cErrNcond= "Error in %s\n Ncond must be greater than 0.\n ";
const char *cErrNsite= "Error in %s\n Nsite must be positive value.\n ";
const char *cErrNumAve= "Error in %s\n NumAve must be positive value.\n ";
const char *cErrExpecInterval= "Error in %s\n ExpecInterval must be positive value.\n ";
const char *cErrLanczos_max="Error in %s\n Lanczos_max must be positive value.\n ";
const char *cErrLanczos_eps="Error in %s\n Lanczos_eps must be positive value.\n ";
const char *cErrLanczosExct="Error in %s\n exct=%d must be greater than 1 and smaller than nvec=%d.\n ";
const char *cErrLanczosTarget="Error in %s\n LanczosTarget=%d must be greater than exct=%d.\n ";


const char *cErrKW="Error: Wrong keywords '%s' in %s.\n";
const char *cErrKW_ShowList="Choose Keywords as follows: \n";
const char *cErrKW_Same="Error: Same keywords exist in %s.\n";
const char *cErrKW_InCorPair="Error: keyword and filename must be set as a pair in %s.\n";

const char *cErrMakeDef="Error: Need to make a def file for %s.\n";
const char *cErrIncorrectDef="Error: incorrect file.\n";
const char *cErrNonHermiteTrans="Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n";
const char *cErrNonHermiteTransForAll= "Error: NonHermite Pair exists in transfer integral. \n";
const char *cErrNonHermiteInterAll="Error: NonHermite (i, spni, j, spnj, k, spnk, l, spnl) = (%d, %d, %d, %d, %d, %d, %d, %d), InterAll_re= %lf, InterAll_im= %lf . \n";
const char *cErrNonConservedInterAll="Error: This operator breaks Sz Symmetry (i, spni, j, spnj, k, spnk, l, spnl) = (%d, %d, %d, %d, %d, %d, %d, %d), InterAll_re= %lf, InterAll_im= %lf . \n";
const char *cErrNonHermiteInterAllForAll= "Error: NonHermite Pair exists in InterAll. \n";
const char *cErrIncorrectFormatForKondoInt= "Error: Site component of (i, j, k, l) =(%d, %d, %d, %d) is incorrect.\n";
const char *cErrIncorrectFormatForKondoTrans= "Error: Site component of (i, j) =(%d, %d) is incorrect.\n";


const char *cErrIncorrectFormatForSpinTrans= "Error: Site component of (i, j) =(%d, %d) is incorrect.\n";
const char *cWarningIncorrectFormatForSpin2= "Warning: Site component of (i, j) =(%d, %d) is ignored.\n";
const char *cErrIncorrectFormatInter= "Error: Use only InterAll for setteing interactions for general spin.\n";
const char *cErrIncorrectSpinIndexForInter="Error: Spin index is incorrect for interactions defined in InterAll file.\n";
const char *cErrIncorrectSpinIndexForTrans="Error: Spin index is incorrect for transfers defined in Trans file.\n";

//! Error Message in CheckMPI.c
const char *cErrNProcNumberHubbard = "Error ! The number of PROCESS should be 4-exponent !\n";
const char *cErrNProcNumberSpin = "Error ! The number of PROCESS should be 2-exponent !\n";
const char *cErrNProcNumberGneralSpin = "Error ! The number of PROCESS is wrong !\n";
const char *cErrNProcNumber = "        The number of PROCESS : %d\n";
const char *cErrNProcNumberSet = "        Set the number of PROCESS as %d or %d.\n";


//! Error Message in diagonal calc.c
const char *cErrNoModel ="Error: CalcModel %d is incorrect.\n";
const char *cErrNoHilbertSpace = "Error: There is no target of Hilbert space.\n";

//! Error Message in bitcalc.c
const char *cErrSiteNumber = "Error: Total Site Number is incorrect.\n";


//! Error Message in mltiply.c


//! Error Message in FileIO.c
const char *cErrFIOpen ="FileOpenError: %s.\n";

//! Error Message in sz.c
const char* cErrSz="Error: in sz. \n";
const char* cErrSz_NoFile="No file. Please set READ=0.\n";
const char* cErrSz_NoFile_Show=" %s does not exist. \n";
