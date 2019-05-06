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
long int *list_1; /**< list of getting real-space configuration for canonical state*/
long int *list_1buf;/**< list of getting real-space configuration for canonical state across processes*/
long int *list_2_1;/**< list to get index of list_1*/
long int *list_2_2;/**< list to get index of list_1*/

/*[s] For Spectrum */
long int *list_1_org; /**< list of getting real-space configuration for canonical state before excitation*/
long int *list_1buf_org;/**< list of getting real-space configuration for canonical state before excitation across processes*/
long int *list_2_1_org;/**< list to get index of list_1_org*/
long int *list_2_2_org;/**< list to get index of list_1_org*/
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


/**@page page_variable Global variables and Data structure
  In HPhi, global variables are used. List of them can be found in global.h
  Sometimes, we pass variables to the function as
  @code
    func(&(X.Bind.Def))
  @endcode
 This C-structure is defined in struct.hpp.
*/

