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

namespace Param {
  extern double Tinit;
  extern double TimeSlice;
  extern int OutputInterval;
  extern int ExpandCoef;
  extern int ExpecInterval;
}/*namespace Param*/
/**
@brief Definision of system (Hamiltonian) etc.
*/
namespace Def {
  extern char* CDataFileHead;
  extern char* CParaFileHead;
  extern int nvec;/**<@brief Read from Calcmod in readdef.h*/
  extern int k_exct;/**<@brief Read from Calcmod in readdef.h*/
  extern int LanczosEps;/**<@brief log(10 base) of the convergence threshold.
                 Read from Calcmod in readdef.h*/
  extern int LanczosTarget;/**<@brief Which eigenstate is used to check convergence.
                    Read from Calcmod in readdef.h.*/
  extern int read_hacker;/**<@brief Whether use an efficient method (=1) in sz.c or not (=0)*/
  extern int READ;/**<@brief It is ALWAYS 0 ???*/
  extern int WRITE;/**<@brief It is ALWAYS 0 ???*/

  extern int Nsite;/**<@brief Number of sites in the INTRA process region*/
  extern int NsiteMPI;/**<@brief Total number of sites, differ from DefineList::Nsite*/
  extern int Nup;/**<@brief Number of spin-up electrons in this process. */
  extern int Ndown;/**<@brief Number of spin-down electrons in this process. */
  extern int NupMPI;/**<@brief Total number of spin-up electrons across processes.
                      Deffer from DefineList::Nup. Read from modpara in readdef.h*/
  extern int NdownMPI;/**<@brief Total number of spin-down electrons across processes.
                      Deffer from DefineList::Ndown. Read from modpara in readdef.h*/
  extern int NupOrg;/**<@brief Number of spin-up electrons before exitation. Used only in
                      the spectrum calculation. Read from modpara in readdef.h*/
  extern int NdownOrg;/**<@brief Number of spin-down electrons before exitation. Used only in
                      the spectrum calculation. Read from modpara in readdef.h*/

  extern int Total2Sz;/**<@brief Total @f$2S_z@f$ in this process.*/
  extern int Total2SzMPI;/**<@brief Total @f$2S_z@f$ across processes.*/
  extern int Ne;/**<@brief Number of electrons in this process.*/
  extern int NeMPI;/**<@brief Total number of electrons across process.
                     Differ from DefineList::Ne .*/
  extern int Lanczos_max;/**<@brief Maximum number of iterations.*/
  extern int Lanczos_restart;/**<@brief Number of iterations performed in the restart computation.*/
  extern long int initial_iv;/**<@brief Seed of random number for initial guesss of wavefunctions.*/

  extern int istep;/**<@brief Index of TPQ step ???*/
  extern int irand;/**<@brief Input keyword TargetTPQRand ???*/
  extern int St;/**<@brief 0 or 1, but it affects nothing.*/

  extern int* LocSpn;/**<@brief [DefineList::NLocSpn] Flag (and size) of the local spin.
              malloc in setmem_def().*/
  extern int NLocSpn;/**<@brief Number of local spins*/
  extern int NCond;/**<@brief Number of itinerant electrons*/
  extern int iFlgGeneralSpin;/**<@brief Flag for the general (Sz/=1/2) spin*/
  extern int iFlgSzConserved;/**<@brief Flag whether Sz is conserved.*/

  extern int fidx;/**<@brief Always 0, it is not used ???*/
  extern long int* Tpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$
                          malloc in setmem_def().*/
  extern long int* OrgTpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$
                             malloc in setmem_def().*/
  extern long int* SiteToBit;/**<@brief [DefineList::NsiteMPI] Similar to DefineList::Tpow.
                      For general spin.*/

  extern int EDNChemi;/**<@brief Number of on-site term.*/
  extern int* EDChemi;/**<@brief [DefineList::Nsite] Chemical potential. malloc in setmem_def().*/
  extern int* EDSpinChemi;/**<@brief [DefineList::Nsite]*/
  extern double* EDParaChemi;/**<@brief [DefineList::Nsite] On-site potential parameter.
                      malloc in setmem_def().*/

                      //[s] Transfer
  extern int NTransfer;/**<@brief Number of transfer integrals obtained by a def file.*/
  extern int EDNTransfer;/**<@brief Number of transfer integrals for calculation. */
  extern int** GeneralTransfer;/**<@brief Index of transfer integrals obtained by a def file.
                        malloc in setmem_def().\n
                        Data Format [DefineList::NTransfer][4]:
                        0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  extern int** EDGeneralTransfer;/**<@brief Index of transfer integrals for calculation.
                          malloc in setmem_def().\n
                          Data Format [DefineList::NTransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  extern std::complex<double>* ParaGeneralTransfer;/**<@brief Value of general transfer integrals by a def file.
                                      malloc in setmem_def().\n
                                      Data Format [DefineList::NTransfer].*/
  extern std::complex<double>* EDParaGeneralTransfer;/**<@brief Value of general transfer integrals  by a def file.
                                        malloc in setmem_def().\n
                                        Data Format [DefineList::NTransfer].*/
                                        //[e] Transfer

  extern int NCoulombIntra;/**< Number of on-site Coulomb interaction*/
  extern int** CoulombIntra;/**< [DefineList::NCoulombIntra][1] Index of on-site coulomb interaction.
                     malloc in setmem_def().*/
  extern double* ParaCoulombIntra;/**< [DefineList::NCoulombIntra] Coupling constant of on-site
                           Coulomb interaction. malloc in setmem_def().*/

  extern int NCoulombInter;/**<@brief Number of off-site Coulomb interaction*/
  extern int** CoulombInter;/**< [DefineList::NCoulombInter][2] Index of off-site coulomb interaction.
                     malloc in setmem_def().*/
  extern double* ParaCoulombInter;/**<@brief [DefineList::NCoulombInter]Coupling constant of off-site
                           Coulomb interaction. malloc in setmem_def().*/

  extern int NHundCoupling;/**<@brief Number of Hund coupling*/
  extern int** HundCoupling;/**<@brief [DefineList::NHundCoupling][2] Index of Hund coupling.
                     malloc in setmem_def().*/
  extern double* ParaHundCoupling;/**<@brief [DefineList::NHundCoupling] Hund coupling constant.
                           malloc in setmem_def().*/

  extern int NPairHopping;/**<@brief Number of pair-hopping term*/
  extern int** PairHopping;/**<@brief [DefineList::NPairHopping][2] Index of pair-hopping.
                    malloc in setmem_def().*/
  extern double* ParaPairHopping;/**<@brief [DefineList::NPairHopping] Coupling constant of
                          pair-hopping term. malloc in setmem_def().*/

  extern int NExchangeCoupling;/**<@brief Number of exchange term*/
  extern int** ExchangeCoupling;/**<@brief [DefineList::NExchangeCoupling][2] Index of exchange term.
                         malloc in setmem_def().*/
  extern double* ParaExchangeCoupling;/**<@brief [DefineList::NExchangeCoupling] Coupling constant of
                               exchange term. malloc in setmem_def().*/

  extern int NIsingCoupling;/**<@brief Number of Ising term.*/

  extern int NPairLiftCoupling;/**<@brief Number of pair-lift term*/
  extern int** PairLiftCoupling;/**<@brief [DefineList::NPairHopping][2] Index of pair-lift term.
                         malloc in setmem_def().*/
  extern double* ParaPairLiftCoupling;/**<@brief [DefineList::NPairHopping] Coupling constant of
                               pair-lift term. malloc in setmem_def().*/

                               //[s] For InterAll
  extern int** InterAll;/**<@brief [DefineList::NinterAll][8] Interacted quartet*/
  extern int** InterAll_OffDiagonal;/**<@brief [DefineList::NinterAll_OffDiagonal][8] Interacted quartet*/
  extern int** InterAll_Diagonal;/**<@brief [DefineList::NinterAll_Diagonal][4] Interacted quartet*/
  extern int NInterAll;/**<@brief Total Number of Interacted quartet*/
  extern int NInterAll_Diagonal;/**<@brief Number of interall term (diagonal)*/
  extern int NInterAll_OffDiagonal;/**<@brief Number of interall term (off-diagonal)*/
  extern std::complex<double>* ParaInterAll;/**<@brief [DefineList::NInterAll] Coupling constant of
                               inter-all term. malloc in setmem_def().*/
  extern double* ParaInterAll_Diagonal;/**<@brief [DefineList::NInterAll_Diagonal] Coupling constant of
                               diagonal inter-all term. malloc in setmem_def().*/
  extern std::complex<double>* ParaInterAll_OffDiagonal;/**<@brief [DefineList::NInterAll_OffDiagonal] Coupling constant of
                               off-diagonal inter-all term. malloc in setmem_def().*/
                               //[e] For InterAll

  extern int** CisAjt;/**<@brief [DefineList::NCisAjt][4] Indices of one-body correlation function. malloc in setmem_def().*/
  extern int NCisAjt;/**<@brief Number of indices of two-body correlation function.*/

  extern int** CisAjtCkuAlvDC;/**<@brief [DefineList::NCisAjtCkuAlvDC][4] Indices of two-body correlation function. malloc in setmem_def().*/
  extern int NCisAjtCkuAlvDC;/**<@brief Number of indices of two-body correlation function.*/

  extern int*** SingleExcitationOperator;/**<@brief [DefineList::NSingleExcitationOperator][3]
                                 Indices of single excitaion operator for spectrum. malloc in setmem_def().*/
  extern int NNSingleExcitationOperator;/**<@brief Number of single excitaion operator for spectrum.*/
  extern int* NSingleExcitationOperator;/**<@brief Number of single excitaion operator for spectrum.*/
  extern std::complex<double>** ParaSingleExcitationOperator;/**<@brief [DefineList::NSingleExcitationOperator]
              Coefficient of single excitaion operator for spectrum. malloc in setmem_def().*/

  extern int*** PairExcitationOperator;/**<@brief [DefineList::NPairExcitationOperator][5]
                               Indices of pair excitaion operator for spectrum. malloc in setmem_def().*/
  extern int NNPairExcitationOperator;/**<@brief Number of pair excitaion operator for spectrum.*/
  extern int* NPairExcitationOperator;/**<@brief Number of pair excitaion operator for spectrum.*/
  extern std::complex<double>** ParaPairExcitationOperator;/**<@brief [DefineList::NPairExcitationOperator]
                           Coefficient of pair excitaion operator for spectrum. malloc in setmem_def().*/

  extern int iCalcType;/**<@brief Switch for calculation type. 0:Lanczos, 1:TPQCalc, 2:FullDiag.*/
  extern int iCalcEigenVec;/**<@brief Switch for method to calculate eigenvectors.
                    0:Lanczos+CG, 1: Lanczos. default value is set as 0 in readdef.c*/
  extern int iInitialVecType;/**<@brief Switch for type of inital vectors.
                      0:complex type, 1: real type. default value is set as 0 in readdef.c*/
  extern int iFlgFiniteTemperature;/**<@brief ???*/
  extern int iCalcModel;/**<@brief Switch for model. 0:Hubbard, 1:Spin, 2:Kondo,
                 3:HubbardGC, 4:SpinGC, 5:KondoGC, 6:HubbardNConserved*/
  extern int iOutputMode;/**<@brief Switch for output mode. 0: OneBodyG and TwoBodyG.
                  1: OneBodyG and TwoBodyG and correlations for charge and spin.*/
  extern int iOutputEigenVec;/**<@brief ASwitch for outputting an eigenvector. 0: no output, 1:output.*/
  extern int iInputEigenVec;/**<@brief Switch for reading an eigenvector. 0: no input, 1:input*/
  extern int iOutputHam;/**<brief Switch for outputting a Hamiltonian. 0: no output, 1:output*/
  extern int iInputHam;/**<brief Switch for reading a Hamiltonian. 0: no input, 1:input*/
  extern int iOutputExVec; /**<brief Switch for outputting an excited vector. 0: no output, 1:output*/

    //[s] For Spectrum
  extern std::complex<double> dcOmegaMax;/**<@brief Upper limit of the frequency for the spectrum.*/
  extern std::complex<double> dcOmegaMin;/**<@brief Lower limit of the frequency for the spectrum.*/
  extern std::complex<double> dcOmegaOrg;/**<@brief Origin limit of the frequency for the spectrum.*/
  extern int iNOmega;/**<@brief Number of frequencies for spectrum.*/
  extern int iFlgSpecOmegaMax;/**<@brief Whether DefineList::dcOmegaMax is input or not.*/
  extern int iFlgSpecOmegaMin;/**<@brief Whether DefineList::dcOmegaMin is input or not.*/
  extern int iFlgSpecOmegaOrg;/**<@brief Whether DefineList::dcOmegaOrg is input or not.*/
  extern int iFlgCalcSpec;/**<@brief Input parameter CalcSpec in teh CalcMod file.*/
  extern int iFlagListModified;/**<@brief When the Hilbert space of excited state differs from the original one.*/

    //[e] For Spectrum

  extern int iReStart;/**< An integer for restarting output a Hamiltonian.
     - 0: not restart
     - 1:restart (output restart vector),
     - 2: restart (input and output restart vector) */
  extern int iFlgMPI;/**<@brief MPI mode
    - 0: butterfly
    - 1: Parallel Interaction [to be supported]
    */

  extern int iNGPU;/**<@brief GPU mode ( only for FullDiag )
  - 0: Use lapack
  - >0: Use GPU
  */

  extern int iFlgScaLAPACK;/**<@brief ScaLAPACK mode ( only for FullDiag )
  - 0: Use lapack
  - 1: Use ScaLAPACK
  */

  //[s] For Time Evolution

  //Information of Time
  extern int NTETimeSteps;
  extern double* TETime;

  //[s]For Ido-san version
  extern int NLaser;
  extern double* ParaLaser;
  //[e]For Ido-san version

  //Information of Transfer integrals
  extern int NTETransferMax;
  extern int* NTETransfer;        /**< Number of time-dependent transfer integrals for Time Evolution.\n
             Data Format [NTE]*/
  extern int* NTETransferDiagonal;        /**< Number of time-dependent transfer integrals for Time Evolution.\n
             Data Format [NTE]*/
  extern int*** TETransfer;      /**< Index of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  extern int*** TETransferDiagonal;      /**< Index of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer][2]: 0->site number i, 1-> spin index on i. */
  extern std::complex<double>** ParaTETransfer;  /**< Value of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */
  extern double** ParaTETransferDiagonal;  /**< Value of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */

             //Two-body part
  extern int NTEInterAllMax;
  extern int* NTEInterAll;        /**< Number of time-dependent InterAll for Time Evolution.\n
             Data Format [NTE]*/
  extern int* NTEInterAllOffDiagonal;        /**< Number of off-diagonal part of time-dependent InterAll for Time Evolution.\n
             Data Format [NTE]*/

  extern int* NTEInterAllDiagonal;        /**< Number of diagonal part of time-dependent InterAll for Time Evolution.\n
             Data Format [NTE]*/
  extern int*** TEInterAll;      /**< Index of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][NTEInterAll][8]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j.
             4->site number k, 5-> spin index on k, 6-> site number l, 7-> spin index on l.*/
  extern int*** TEInterAllOffDiagonal;      /**< Index of off-diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][NTEInterAll][8]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j.
             4->site number k, 5-> spin index on k, 6-> site number l, 7-> spin index on l.*/
  extern int*** TEInterAllDiagonal;      /**< Index of diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][NTEInterAll][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  extern std::complex<double>** ParaTEInterAll;  /**< Value of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */
  extern std::complex<double>** ParaTEInterAllOffDiagonal;  /**< Value of off-diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */

  extern double** ParaTEInterAllDiagonal;  /**< Value of diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */
  extern int** TEChemi;    /**< [NTE][Nsite] */
  extern int* NTEChemi;   /**< [NTE] */
  extern int** SpinTEChemi;  /**< [NTE][Nsite] */
  extern double** ParaTEChemi;  /**< [NTE][Nsite] */
  //[e] For Time Evolution
};/*namespace Def */
/**
@brief Size of the Hilbert space
*/
namespace Check {
  extern long int idim_max;/**<@brief The dimension of the Hilbert space of this process.*/
  extern long int idim_maxMPI;/**<@brief The total dimension across process.*/
  extern long int idim_maxOrg;/**<@brief The local Hilbert-space dimention of original state for the spectrum.*/
  extern long int idim_maxMPIOrg;/**<@brief The global Hilbert-space dimention of original state for the spectrum.*/
  extern long int sdim;/**<@brief Dimension for Ogata-Lin ???*/
  extern double max_mem;/**<@brief Estimated memory size.*/
};/*namespace Check*/
/**
@brief For Matrix-Vector product
*/
namespace Large {
  extern int itr;/**<@brief Iteration number.*/
  extern long int iv;/**<@brief Used for initializing vector.*/
  extern long int i_max;/**<@brief Length of eigenvector*/
  extern long int SizeOflist_2_1;/**<@brief Size of ::list_2_1*/
  extern long int SizeOflist_2_2;/**<@brief Size of ::list_2_2*/
  extern long int SizeOflistjb;/**<@brief Used for computing Sz.*/

  extern std::complex<double> tmp_trans;/**<@brief Hopping parameter.*/
  extern std::complex<double> tmp_J;/**<@brief Coupling constant*/

  extern long int is1_up;/**<@brief Mask used in the bit oeration.*/
  extern long int is1_down;/**<@brief Mask used in the bit oeration.*/
  extern long int is2_up;/**<@brief Mask used in the bit oeration.*/
  extern long int is2_down;/**<@brief Mask used in the bit oeration.*/

  extern int mode;/**<@brief multiply or expectation value.*/
  extern double sgn;/**<@brief Not used ???*/
  extern long int is1_spin;/**<@brief Mask used in the bit oeration.*/
  extern long int is2_spin;/**<@brief Mask used in the bit oeration.*/
  extern long int is3_spin;/**<@brief Mask used in the bit oeration.*/
  extern long int is4_spin;/**<@brief Mask used in the bit oeration.*/
  extern int isite1;/**<@brief Is it realy used ???*/
  extern int isite2;/**<@brief Is it realy used ???*/
  extern int isite3;/**<@brief Is it realy used ???*/
  extern int isite4;/**<@brief Is it realy used ???*/

  extern long int A_spin;/**<@brief Mask used in the bit oeration.*/
  extern long int B_spin;/**<@brief Mask used in the bit oeration.*/
  extern long int irght;/**<@brief Used for Ogata-Lin ???*/
  extern long int ilft;/**<@brief Used for Ogata-Lin ???*/
  extern long int ihfbit;/**<@brief Used for Ogata-Lin ???*/
  extern long int isA_spin;/**<@brief Mask used in the bit oeration.*/
  extern long int isB_spin;/**<@brief Mask used in the bit oeration.*/
  extern std::complex<double> tmp_V;/**<@brief Coupling constant*/
};/*namespace Large*/
/**
@brief Physical quantities (Expectation value)
*/
namespace Phys {
  //double energy,doublon;
  extern double* energy;/**<@brief Expectation value of the total energy.*/
  extern double* doublon;/**<@brief Expectation value of the Doublon*/
  extern double* doublon2;/**<@brief Expectation value of the Square of doublon*/
  extern double* num;/**<@brief Expectation value of the Number of electrons*/
  extern double* num2;/**<@brief Expectation value of the quare of the number of electrons*/
  extern double* Sz;/**<@brief Expectation value of the Total Sz*/
  extern double* Sz2;/**<@brief Expectation value of the Square of total Sz*/
  extern double* num_up;/**<@brief Expectation value of the number of up-spin electtrons*/
  extern double* num_down;/**<@brief Expectation value of the number of down-spin electtrons*/
  extern double* s2;/**<@brief Expectation value of the square of the total S.*/
    /*[s] For TPQ*/
  extern double* var;/**<@brief Expectation value of the Energy variance.*/
    /*[e] For TPQ*/

  extern double* spin_real_cor;/**<@brief Malloc, but Not used ???*/
  extern double* charge_real_cor;/**<@brief Malloc, but Not used ???*/
  extern double* loc_spin_z;/**<@brief Malloc, but Not used ???*/
  extern double Target_energy;/**<@brief Is it really used ???*/
  extern double Target_CG_energy;/**<@brief Taget energy of CG-inversed iteration (NOT LOBCG) method.*/
};/*namespace Phys*/
/**
@brief For Boost
*/
namespace Boost {
  extern int flgBoost;/**<@brief Flag whether use CMA algorithm.*/
  extern long int R0;
  extern long int W0;
  extern long int num_pivot;
  extern long int ishift_nspin;
  extern int NumarrayJ;/**<@brief */
  extern std::complex<double>*** arrayJ;/**<@brief */
  extern std::complex<double> vecB[3];/**<@brief */
  extern int** list_6spin_star;/**<@brief */
  extern int*** list_6spin_pair;/**<@brief */
};/*namespace BoostList*/


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
