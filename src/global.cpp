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

#include "global.hpp"

std::complex<double> I(0.0, 1.0);
std::complex<double> **v0;  /**< A vector after multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
std::complex<double> **v1;  /**< A vector before multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
std::complex<double> **v1buf; /**< A temporary vector for MPI. */

double *list_Diagonal; /**< list for diagonal components.*/
long unsigned int *list_1; /**< list of getting real-space configuration for canonical state*/
long unsigned int *list_1buf;/**< list of getting real-space configuration for canonical state across processes*/
long unsigned int *list_2_1;/**< list to get index of list_1*/
long unsigned int *list_2_2;/**< list to get index of list_1*/

/*[s] For Spectrum */
long unsigned int *list_1_org; /**< list of getting real-space configuration for canonical state before excitation*/
long unsigned int *list_1buf_org;/**< list of getting real-space configuration for canonical state before excitation across processes*/
long unsigned int *list_2_1_org;/**< list to get index of list_1_org*/
long unsigned int *list_2_2_org;/**< list to get index of list_1_org*/
/*[e] For Spectrum */

/*[s] For Lanczos */
int     initial_mode;/**< mode to get initial state (0: use same random generator for MPI, 1: use each random generator for MPI)*/
/*[e] For Lanczos */

/*[s] For TPQ*/
double LargeValue;/**< constant value l for TPQ calculation.*/
int    NumAve;/**< Average number for TPQ calculation*/
int step_i;/**< step for TPQ calculation*/
double *global_norm;/**< norm before normalization for TPQ calculation*/
double *global_1st_norm;/**< 1-st norm for TPQ calculation*/
int step_spin;/**< output step for TE calculation.*/
/*[e] For TPQ*/

/*[s] For All Diagonalization*/
#ifdef _SCALAPACK
std::complex<double> *Z_vec; /**> distributed matrix of eigen vector*/
int descZ_vec[9]; /*descriptor for Z_vec*/
#endif
/*[e] For All Diagonalization*/

/*OutputPath Name*/
const char* cParentOutputFolder = "./output/";

//For TimeKeep
const char* cFileNameTimeKeep="%s_TimeKeeper.dat";
const char* cFileNameSzTimeKeep="%s_sz_TimeKeeper.dat";

//For Check
const char* cFileNameCheckCoulombIntra="CHECK_CoulombIntra.dat";
const char* cFileNameCheckChemi="CHECK_Chemi.dat";
const char* cFileNameCheckInterU="CHECK_INTER_U.dat";
const char* cFileNameCheckHund="CHECK_Hund.dat";
const char* cFileNameCheckInterAll="CHECK_InterAll.dat";
const char* cFileNameCheckMemory="CHECK_Memory.dat";
const char* cFileNameCheckSdim="CHECK_Sdim.dat";

//For EDTrans
const char* cFileNameWarningOnTransfer="WarningOnTransfer.dat";

//For Lanczos
const char* cFileNameLanczosStep="%s_Lanczos_Step.dat";
const char* cFileNameEnergy_Lanczos= "%s_energy.dat";
const char* cFileNameEigenvalue_Lanczos= "Eigenvalue.dat";
const char* cFileNameEnergy_CG="%s_energy.dat";
const char* cFileName1BGreen_Lanczos="%s_cisajs.dat";
const char* cFileName1BGreen_CG="%s_cisajs.dat";
const char* cFileName2BGreen_Lanczos="%s_cisajscktalt.dat";
const char* cFileName2BGreen_CG="%s_cisajscktalt.dat";
const char* cFileNameTimeEV_CG="Time_EigenVector.dat";
const char* cFileNameListModel="ListForModel_Ns%d_Nup%dNdown%d.dat";
const char* cFileNameOutputEigen="%s_eigenvec_%d_rank_%d.dat";
const char* cFileNameInputEigen="%s_eigenvec_%d_rank_%d.dat";
const char* cFileNameCalcDynamicalGreen="%s_DynamicalGreen.dat";
const char* cFileNameTridiagonalMatrixComponents="%s_TMComponents.dat";


//For TPQ
const char* cFileNameSSRand="SS_rand%d.dat";
const char* cFileNameTPQStep="%s_Time_TPQ_Step.dat";
const char* cFileNameNormRand="Norm_rand%d.dat";
const char* cFileNameFlctRand="Flct_rand%d.dat";
const char* cFileName1BGreen_TPQ="%s_cisajs_set%dstep%d.dat";
const char* cFileName2BGreen_TPQ="%s_cisajscktalt_set%dstep%d.dat";
const char* cFileName1BGreen_TE="%s_cisajs_step%d.dat";
const char* cFileName2BGreen_TE="%s_cisajscktalt_step%d.dat";
const char* cFileNameOutputVector="tmpvec_set%d_rank_%d.dat";
const char* cFileNameInputVector="tmpvec_set%d_rank_%d.dat";

//Fot Time evolution
const char* cFileNameTEStep="%s_Time_TE_Step.dat";
const char* cFileNameSS="SS.dat";
const char* cFileNameNorm="Norm.dat";
const char* cFileNameFlct="Flct.dat";

//For FullDiag
const char* cFileNamePhys_FullDiag="%s_phys_Nup%d_Ndown%d.dat";
const char* cFileNamePhys_FullDiag_GC="%s_phys.dat";
const char* cFileName1BGreen_FullDiag="%s_cisajs_eigen%d.dat";
const char* cFileName2BGreen_FullDiag="%s_cisajscktalt_eigen%d.dat";
const char* cFileNamePhys_FullDiag_Ham="%s_Ham.dat";

//For Spectrum
const char* cFileNameOutputRestartVec="%s_recalcvec_rank_%d.dat";
const char* cFileNameOutputExcitedVec="%s_excitedvec_rank_%d.dat";
//For Error
const char* cFileNameErrorSz="Err_sz.dat";

//For Timer
double *Timer; /**> The procedure execution time.*/
double *TimerStart;/**> Timer when the procedure starts.*/

/********************************************************************/
/********************************************************************/
double eps; /**> epsilon used in getting spectrum by Lanczos method and Lanczos eigenvector by CG method.*/
double eps_CG;/**> epsilon used in getting Lanczos eigenvector by CG method.*/
double eps_Lanczos;/**> epsilon used in LOBPCG, BiCG and Lanczos eigen value.*/
double eps_Energy;/**> epsilon for energy*/
double eps_CheckImag0;/**> epsilon for checking values of one-body and two-body interactions.*/

/*
 Variables for the MPI parallelism
*/
int nproc;//!< Number of processors, defined in InitializeMPI()
int myrank;//!< Process ID, defined in InitializeMPI()
int nthreads;//!< Number of Threads, defined in InitializeMPI()
FILE *stdoutMPI;/**<@brief File pointer to the standard output
                defined in InitializeMPI()*/
