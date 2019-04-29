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

#ifndef HPHI_GLOBAL_H
#define HPHI_GLOBAL_H

#include <complex>
#include <cmath>
#include <cstdio>
#define D_FileNameMax 256
#define MPIFALSE -1
#define FALSE 0
#define TRUE 1

/**
 * Kind of electron. For ver.0.1, ITINERANT=1. From ver.0.2, ITINERANT=0.
 **/
#define ITINERANT 0
#define LOCSPIN 1

extern std::complex<double> I;
extern std::complex<double> **v0;  /**< A vector after multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
extern std::complex<double> **v1;  /**< A vector before multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
extern std::complex<double> **v1buf; /**< A temporary vector for MPI. */

extern double *list_Diagonal; /**< list for diagonal components.*/
extern long int *list_1; /**< list of getting real-space configuration for canonical state*/
extern long int *list_1buf;/**< list of getting real-space configuration for canonical state across processes*/
extern long int *list_2_1;/**< list to get index of list_1*/
extern long int *list_2_2;/**< list to get index of list_1*/

/*[s] For Spectrum */
extern long int *list_1_org; /**< list of getting real-space configuration for canonical state before excitation*/
extern long int *list_1buf_org;/**< list of getting real-space configuration for canonical state before excitation across processes*/
extern long int *list_2_1_org;/**< list to get index of list_1_org*/
extern long int *list_2_2_org;/**< list to get index of list_1_org*/
/*[e] For Spectrum */

/*[s] For Lanczos */
extern int     initial_mode;/**< mode to get initial state (0: use same random generator for MPI, 1: use each random generator for MPI)*/
/*[e] For Lanczos */

/*[s] For TPQ*/
extern double LargeValue;/**< constant value l for TPQ calculation.*/
extern int    NumAve;/**< Average number for TPQ calculation*/
extern int step_i;/**< step for TPQ calculation*/
extern double *global_norm;/**< norm before normalization for TPQ calculation*/
extern double *global_1st_norm;/**< 1-st norm for TPQ calculation*/
extern int step_spin;/**< output step for TE calculation.*/
/*[e] For TPQ*/

/*[s] For All Diagonalization*/
#ifdef _SCALAPACK
extern std::complex<double> *Z_vec; /**> distributed matrix of eigen vector*/
extern int descZ_vec[9]; /*descriptor for Z_vec*/
#endif
/*[e] For All Diagonalization*/

//For Timer
extern double *Timer; /**> The procedure execution time.*/
extern double *TimerStart;/**> Timer when the procedure starts.*/

/********************************************************************/
/********************************************************************/
extern double eps; /**> epsilon used in getting spectrum by Lanczos method and Lanczos eigenvector by CG method.*/
extern double eps_CG;/**> epsilon used in getting Lanczos eigenvector by CG method.*/
extern double eps_Lanczos;/**> epsilon used in LOBPCG, BiCG and Lanczos eigen value.*/
extern double eps_Energy;/**> epsilon for energy*/
extern double eps_CheckImag0;/**> epsilon for checking values of one-body and two-body interactions.*/

/*
 Variables for the MPI parallelism
*/
extern int nproc;//!< Number of processors, defined in InitializeMPI()
extern int myrank;//!< Process ID, defined in InitializeMPI()
extern int nthreads;//!< Number of Threads, defined in InitializeMPI()
extern FILE *stdoutMPI;/**<@brief File pointer to the standard output
                defined in InitializeMPI()*/

#endif /* HPHI_GLOBAL_H */

/**
@page page_variable Global variables and Data structure

In HPhi, global variables are used. List of them can be found in global.h

Sometimes, we pass variables to the function as
@code
func(&(X.Bind.Def))
@endcode
This C-structure is defined in struct.h.
*/

