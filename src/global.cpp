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


namespace Param {
  double Tinit;
  double TimeSlice;
  int OutputInterval = 0;
  int ExpandCoef = 0;
  int ExpecInterval = 0;
}/*namespace Param*/
/**
@brief Definision of system (Hamiltonian) etc.
*/
namespace Def {
  char* CDataFileHead;/**<@brief Read from Calcmod in readdef.h.
                      Header of output file such as Green's function*/
  char* CParaFileHead;/**<@brief Read from Calcmod in readdef.h.
                      It is not used. Just for the compatibility to mVMC*/
  int nvec = 0;/**<@brief Read from Calcmod in readdef.h*/
  int k_exct = 0;/**<@brief Read from Calcmod in readdef.h*/
  int LanczosEps = 0;/**<@brief log(10 base) of the convergence threshold.
                 Read from Calcmod in readdef.h*/
  int LanczosTarget = 0;/**<@brief Which eigenstate is used to check convergence.
                    Read from Calcmod in readdef.h.*/
  int read_hacker = 0;/**<@brief Whether use an efficient method (=1) in sz.c or not (=0)*/
  int READ = 0;/**<@brief It is ALWAYS 0 ???*/
  int WRITE = 0;/**<@brief It is ALWAYS 0 ???*/

  int Nsite = 0;/**<@brief Number of sites in the INTRA process region*/
  int NsiteMPI = 0;/**<@brief Total number of sites, differ from DefineList::Nsite*/
  int Nup = 0;/**<@brief Number of spin-up electrons in this process. */
  int Ndown = 0;/**<@brief Number of spin-down electrons in this process. */
  int NupMPI = 0;/**<@brief Total number of spin-up electrons across processes.
                      Deffer from DefineList::Nup. Read from modpara in readdef.h*/
  int NdownMPI = 0;/**<@brief Total number of spin-down electrons across processes.
                      Deffer from DefineList::Ndown. Read from modpara in readdef.h*/
  int NupOrg = 0;/**<@brief Number of spin-up electrons before exitation. Used only in
                      the spectrum calculation. Read from modpara in readdef.h*/
  int NdownOrg = 0;/**<@brief Number of spin-down electrons before exitation. Used only in
                      the spectrum calculation. Read from modpara in readdef.h*/

  int Total2Sz = 0;/**<@brief Total @f$2S_z@f$ in this process.*/
  int Total2SzMPI = 0;/**<@brief Total @f$2S_z@f$ across processes.*/
  int Ne = 0;/**<@brief Number of electrons in this process.*/
  int NeMPI = 0;/**<@brief Total number of electrons across process.
                     Differ from DefineList::Ne .*/
  int Lanczos_max = 0;/**<@brief Maximum number of iterations.*/
  int Lanczos_restart = 0;/**<@brief Number of iterations performed in the restart computation.*/
  long int initial_iv = 0;/**<@brief Seed of random number for initial guesss of wavefunctions.*/

  int istep = 0;/**<@brief Index of TPQ step ???*/
  int irand = 0;/**<@brief Input keyword TargetTPQRand ???*/
  int St = 0;/**<@brief 0 or 1, but it affects nothing.*/

  int* LocSpn;/**<@brief [DefineList::NLocSpn] Flag (and size) of the local spin.
              malloc in setmem_def().*/
  int NLocSpn = 0;/**<@brief Number of local spins*/
  int NCond = 0;/**<@brief Number of itinerant electrons*/
  int iFlgGeneralSpin = 0;/**<@brief Flag for the general (Sz/=1/2) spin*/
  int iFlgSzConserved = 0;/**<@brief Flag whether Sz is conserved.*/

  int fidx = 0;/**<@brief Always 0, it is not used ???*/
  long int* Tpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$
                          malloc in setmem_def().*/
  long int* OrgTpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$
                             malloc in setmem_def().*/
  long int* SiteToBit;/**<@brief [DefineList::NsiteMPI] Similar to DefineList::Tpow.
                      For general spin.*/

  int EDNChemi = 0;/**<@brief Number of on-site term.*/
  int* EDChemi;/**<@brief [DefineList::Nsite] Chemical potential. malloc in setmem_def().*/
  int* EDSpinChemi;/**<@brief [DefineList::Nsite]*/
  double* EDParaChemi;/**<@brief [DefineList::Nsite] On-site potential parameter.
                      malloc in setmem_def().*/

                      //[s] Transfer
  int NTransfer = 0;/**<@brief Number of transfer integrals obtained by a def file.*/
  int EDNTransfer = 0;/**<@brief Number of transfer integrals for calculation. */
  int** GeneralTransfer;/**<@brief Index of transfer integrals obtained by a def file.
                        malloc in setmem_def().\n
                        Data Format [DefineList::NTransfer][4]:
                        0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  int** EDGeneralTransfer;/**<@brief Index of transfer integrals for calculation.
                          malloc in setmem_def().\n
                          Data Format [DefineList::NTransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  std::complex<double>* ParaGeneralTransfer;/**<@brief Value of general transfer integrals by a def file.
                                      malloc in setmem_def().\n
                                      Data Format [DefineList::NTransfer].*/
  std::complex<double>* EDParaGeneralTransfer;/**<@brief Value of general transfer integrals  by a def file.
                                        malloc in setmem_def().\n
                                        Data Format [DefineList::NTransfer].*/
                                        //[e] Transfer

  int NCoulombIntra = 0;/**< Number of on-site Coulomb interaction*/
  int** CoulombIntra;/**< [DefineList::NCoulombIntra][1] Index of on-site coulomb interaction.
                     malloc in setmem_def().*/
  double* ParaCoulombIntra;/**< [DefineList::NCoulombIntra] Coupling constant of on-site
                           Coulomb interaction. malloc in setmem_def().*/

  int NCoulombInter = 0;/**<@brief Number of off-site Coulomb interaction*/
  int** CoulombInter;/**< [DefineList::NCoulombInter][2] Index of off-site coulomb interaction.
                     malloc in setmem_def().*/
  double* ParaCoulombInter;/**<@brief [DefineList::NCoulombInter]Coupling constant of off-site
                           Coulomb interaction. malloc in setmem_def().*/

  int NHundCoupling = 0;/**<@brief Number of Hund coupling*/
  int** HundCoupling;/**<@brief [DefineList::NHundCoupling][2] Index of Hund coupling.
                     malloc in setmem_def().*/
  double* ParaHundCoupling;/**<@brief [DefineList::NHundCoupling] Hund coupling constant.
                           malloc in setmem_def().*/

  int NPairHopping = 0;/**<@brief Number of pair-hopping term*/
  int** PairHopping;/**<@brief [DefineList::NPairHopping][2] Index of pair-hopping.
                    malloc in setmem_def().*/
  double* ParaPairHopping;/**<@brief [DefineList::NPairHopping] Coupling constant of
                          pair-hopping term. malloc in setmem_def().*/

  int NExchangeCoupling = 0;/**<@brief Number of exchange term*/
  int** ExchangeCoupling;/**<@brief [DefineList::NExchangeCoupling][2] Index of exchange term.
                         malloc in setmem_def().*/
  double* ParaExchangeCoupling;/**<@brief [DefineList::NExchangeCoupling] Coupling constant of
                               exchange term. malloc in setmem_def().*/

  int NIsingCoupling = 0;/**<@brief Number of Ising term.*/

  int NPairLiftCoupling = 0;/**<@brief Number of pair-lift term*/
  int** PairLiftCoupling;/**<@brief [DefineList::NPairHopping][2] Index of pair-lift term.
                         malloc in setmem_def().*/
  double* ParaPairLiftCoupling;/**<@brief [DefineList::NPairHopping] Coupling constant of
                               pair-lift term. malloc in setmem_def().*/

                               //[s] For InterAll
  int** InterAll;/**<@brief [DefineList::NinterAll][8] Interacted quartet*/
  int** InterAll_OffDiagonal;/**<@brief [DefineList::NinterAll_OffDiagonal][8] Interacted quartet*/
  int** InterAll_Diagonal;/**<@brief [DefineList::NinterAll_Diagonal][4] Interacted quartet*/
  int NInterAll = 0;/**<@brief Total Number of Interacted quartet*/
  int NInterAll_Diagonal = 0;/**<@brief Number of interall term (diagonal)*/
  int NInterAll_OffDiagonal = 0;/**<@brief Number of interall term (off-diagonal)*/
  std::complex<double>* ParaInterAll;/**<@brief [DefineList::NInterAll] Coupling constant of
                               inter-all term. malloc in setmem_def().*/
  double* ParaInterAll_Diagonal;/**<@brief [DefineList::NInterAll_Diagonal] Coupling constant of
                               diagonal inter-all term. malloc in setmem_def().*/
  std::complex<double>* ParaInterAll_OffDiagonal;/**<@brief [DefineList::NInterAll_OffDiagonal] Coupling constant of
                               off-diagonal inter-all term. malloc in setmem_def().*/
                               //[e] For InterAll

  int** CisAjt;/**<@brief [DefineList::NCisAjt][4] Indices of one-body correlation function. malloc in setmem_def().*/
  int NCisAjt = 0;/**<@brief Number of indices of two-body correlation function.*/

  int** CisAjtCkuAlvDC;/**<@brief [DefineList::NCisAjtCkuAlvDC][4] Indices of two-body correlation function. malloc in setmem_def().*/
  int NCisAjtCkuAlvDC = 0;/**<@brief Number of indices of two-body correlation function.*/

  int*** SingleExcitationOperator;/**<@brief [DefineList::NSingleExcitationOperator][3]
                                 Indices of single excitaion operator for spectrum. malloc in setmem_def().*/
  int NNSingleExcitationOperator = 0;/**<@brief Number of single excitaion operator for spectrum.*/
  int* NSingleExcitationOperator;/**<@brief Number of single excitaion operator for spectrum.*/
  std::complex<double>** ParaSingleExcitationOperator;/**<@brief [DefineList::NSingleExcitationOperator]
              Coefficient of single excitaion operator for spectrum. malloc in setmem_def().*/

  int*** PairExcitationOperator;/**<@brief [DefineList::NPairExcitationOperator][5]
                               Indices of pair excitaion operator for spectrum. malloc in setmem_def().*/
  int NNPairExcitationOperator = 0;/**<@brief Number of pair excitaion operator for spectrum.*/
  int* NPairExcitationOperator;/**<@brief Number of pair excitaion operator for spectrum.*/
  std::complex<double>** ParaPairExcitationOperator;/**<@brief [DefineList::NPairExcitationOperator]
                           Coefficient of pair excitaion operator for spectrum. malloc in setmem_def().*/

  int iCalcType = 0;/**<@brief Switch for calculation type. 0:Lanczos, 1:TPQCalc, 2:FullDiag.*/
  int iCalcEigenVec = 0;/**<@brief Switch for method to calculate eigenvectors.
                    0:Lanczos+CG, 1: Lanczos. default value is set as 0 in readdef.c*/
  int iInitialVecType = 0;/**<@brief Switch for type of inital vectors.
                      0:complex type, 1: real type. default value is set as 0 in readdef.c*/
  int iFlgFiniteTemperature = 0;/**<@brief ???*/
  int iCalcModel = 0;/**<@brief Switch for model. 0:Hubbard, 1:Spin, 2:Kondo,
                 3:HubbardGC, 4:SpinGC, 5:KondoGC, 6:HubbardNConserved*/
  int iOutputMode = 0;/**<@brief Switch for output mode. 0: OneBodyG and TwoBodyG.
                  1: OneBodyG and TwoBodyG and correlations for charge and spin.*/
  int iOutputEigenVec = 0;/**<@brief ASwitch for outputting an eigenvector. 0: no output, 1:output.*/
  int iInputEigenVec = 0;/**<@brief Switch for reading an eigenvector. 0: no input, 1:input*/
  int iOutputHam = 0;/**<brief Switch for outputting a Hamiltonian. 0: no output, 1:output*/
  int iInputHam = 0;/**<brief Switch for reading a Hamiltonian. 0: no input, 1:input*/
  int iOutputExVec = 0; /**<brief Switch for outputting an excited vector. 0: no output, 1:output*/

    //[s] For Spectrum
  std::complex<double> dcOmegaMax;/**<@brief Upper limit of the frequency for the spectrum.*/
  std::complex<double> dcOmegaMin;/**<@brief Lower limit of the frequency for the spectrum.*/
  std::complex<double> dcOmegaOrg;/**<@brief Origin limit of the frequency for the spectrum.*/
  int iNOmega = 0;/**<@brief Number of frequencies for spectrum.*/
  int iFlgSpecOmegaMax = 0;/**<@brief Whether DefineList::dcOmegaMax is input or not.*/
  int iFlgSpecOmegaMin = 0;/**<@brief Whether DefineList::dcOmegaMin is input or not.*/
  int iFlgSpecOmegaOrg = 0;/**<@brief Whether DefineList::dcOmegaOrg is input or not.*/
  int iFlgCalcSpec = 0;/**<@brief Input parameter CalcSpec in teh CalcMod file.*/
  int iFlagListModified = 0;/**<@brief When the Hilbert space of excited state differs from the original one.*/

    //[e] For Spectrum

  int iReStart = 0;/**< An integer for restarting output a Hamiltonian.
     - 0: not restart
     - 1:restart (output restart vector),
     - 2: restart (input and output restart vector) */
  int iFlgMPI = 0;/**<@brief MPI mode
    - 0: butterfly
    - 1: Parallel Interaction [to be supported]
    */

  int iNGPU = 0;/**<@brief GPU mode ( only for FullDiag )
  - 0: Use lapack
  - >0: Use GPU
  */

  int iFlgScaLAPACK = 0;/**<@brief ScaLAPACK mode ( only for FullDiag )
  - 0: Use lapack
  - 1: Use ScaLAPACK
  */

  //[s] For Time Evolution

  //Information of Time
  int NTETimeSteps = 0;
  double* TETime;

  //[s]For Ido-san version
  int NLaser = 0;
  double* ParaLaser;
  //[e]For Ido-san version

  //Information of Transfer integrals
  int NTETransferMax = 0;
  int* NTETransfer;        /**< Number of time-dependent transfer integrals for Time Evolution.\n
             Data Format [NTE]*/
  int* NTETransferDiagonal;        /**< Number of time-dependent transfer integrals for Time Evolution.\n
             Data Format [NTE]*/
  int*** TETransfer;      /**< Index of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  int*** TETransferDiagonal;      /**< Index of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer][2]: 0->site number i, 1-> spin index on i. */
  std::complex<double>** ParaTETransfer;  /**< Value of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */
  double** ParaTETransferDiagonal;  /**< Value of time-dependent transfer integrals for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */

             //Two-body part
  int NTEInterAllMax = 0;
  int* NTEInterAll;        /**< Number of time-dependent InterAll for Time Evolution.\n
             Data Format [NTE]*/
  int* NTEInterAllOffDiagonal;        /**< Number of off-diagonal part of time-dependent InterAll for Time Evolution.\n
             Data Format [NTE]*/

  int* NTEInterAllDiagonal;        /**< Number of diagonal part of time-dependent InterAll for Time Evolution.\n
             Data Format [NTE]*/
  int*** TEInterAll;      /**< Index of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][NTEInterAll][8]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j.
             4->site number k, 5-> spin index on k, 6-> site number l, 7-> spin index on l.*/
  int*** TEInterAllOffDiagonal;      /**< Index of off-diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][NTEInterAll][8]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j.
             4->site number k, 5-> spin index on k, 6-> site number l, 7-> spin index on l.*/
  int*** TEInterAllDiagonal;      /**< Index of diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][NTEInterAll][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  std::complex<double>** ParaTEInterAll;  /**< Value of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */
  std::complex<double>** ParaTEInterAllOffDiagonal;  /**< Value of off-diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */

  double** ParaTEInterAllDiagonal;  /**< Value of diagonal part of time-dependent InterAll for Time Evolution. \n
             Data Format [NTE][Ntransfer]. */
  int** TEChemi;    /**< [NTE][Nsite] */
  int* NTEChemi;   /**< [NTE] */
  int** SpinTEChemi;  /**< [NTE][Nsite] */
  double** ParaTEChemi;  /**< [NTE][Nsite] */
  //[e] For Time Evolution
};/*namespace Def */
/**
@brief Size of the Hilbert space
*/
namespace Check {
  long int idim_max = 0;/**<@brief The dimension of the Hilbert space of this process.*/
  long int idim_maxMPI = 0;/**<@brief The total dimension across process.*/
  long int idim_maxOrg = 0;/**<@brief The local Hilbert-space dimention of original state for the spectrum.*/
  long int idim_maxMPIOrg = 0;/**<@brief The global Hilbert-space dimention of original state for the spectrum.*/
  long int sdim = 0;/**<@brief Dimension for Ogata-Lin ???*/
  double max_mem;/**<@brief Estimated memory size.*/
};/*namespace Check*/
/**
@brief For Matrix-Vector product
*/
namespace Large {
  int itr = 0;/**<@brief Iteration number.*/
  long int iv = 0;/**<@brief Used for initializing vector.*/
  long int i_max = 0;/**<@brief Length of eigenvector*/
  long int SizeOflist_2_1 = 0;/**<@brief Size of ::list_2_1*/
  long int SizeOflist_2_2 = 0;/**<@brief Size of ::list_2_2*/
  long int SizeOflistjb = 0;/**<@brief Used for computing Sz.*/

  std::complex<double> tmp_trans;/**<@brief Hopping parameter.*/
  std::complex<double> tmp_J;/**<@brief Coupling constant*/

  long int is1_up = 0;/**<@brief Mask used in the bit oeration.*/
  long int is1_down = 0;/**<@brief Mask used in the bit oeration.*/
  long int is2_up = 0;/**<@brief Mask used in the bit oeration.*/
  long int is2_down = 0;/**<@brief Mask used in the bit oeration.*/

  int mode = 0;/**<@brief multiply or expectation value.*/
  double sgn;/**<@brief Not used ???*/
  long int is1_spin = 0;/**<@brief Mask used in the bit oeration.*/
  long int is2_spin = 0;/**<@brief Mask used in the bit oeration.*/
  long int is3_spin = 0;/**<@brief Mask used in the bit oeration.*/
  long int is4_spin = 0;/**<@brief Mask used in the bit oeration.*/
  int isite1 = 0;/**<@brief Is it realy used ???*/
  int isite2 = 0;/**<@brief Is it realy used ???*/
  int isite3 = 0;/**<@brief Is it realy used ???*/
  int isite4 = 0;/**<@brief Is it realy used ???*/

  long int A_spin = 0;/**<@brief Mask used in the bit oeration.*/
  long int B_spin = 0;/**<@brief Mask used in the bit oeration.*/
  long int irght = 0;/**<@brief Used for Ogata-Lin ???*/
  long int ilft = 0;/**<@brief Used for Ogata-Lin ???*/
  long int ihfbit = 0;/**<@brief Used for Ogata-Lin ???*/
  long int isA_spin = 0;/**<@brief Mask used in the bit oeration.*/
  long int isB_spin = 0;/**<@brief Mask used in the bit oeration.*/
  std::complex<double> tmp_V;/**<@brief Coupling constant*/
};/*namespace Large*/
/**
@brief Physical quantities (Expectation value)
*/
namespace Phys {
  //double energy,doublon;
  double* energy;/**<@brief Expectation value of the total energy.*/
  double* doublon;/**<@brief Expectation value of the Doublon*/
  double* doublon2;/**<@brief Expectation value of the Square of doublon*/
  double* num;/**<@brief Expectation value of the Number of electrons*/
  double* num2;/**<@brief Expectation value of the quare of the number of electrons*/
  double* Sz;/**<@brief Expectation value of the Total Sz*/
  double* Sz2;/**<@brief Expectation value of the Square of total Sz*/
  double* num_up;/**<@brief Expectation value of the number of up-spin electtrons*/
  double* num_down;/**<@brief Expectation value of the number of down-spin electtrons*/
  double* s2;/**<@brief Expectation value of the square of the total S.*/
    /*[s] For TPQ*/
  double* var;/**<@brief Expectation value of the Energy variance.*/
    /*[e] For TPQ*/

  double* spin_real_cor;/**<@brief Malloc, but Not used ???*/
  double* charge_real_cor;/**<@brief Malloc, but Not used ???*/
  double* loc_spin_z;/**<@brief Malloc, but Not used ???*/
  double Target_energy;/**<@brief Is it really used ???*/
  double Target_CG_energy;/**<@brief Taget energy of CG-inversed iteration (NOT LOBCG) method.*/
};/*namespace Phys*/
/**
@brief For Boost
*/
namespace Boost {
  int flgBoost = 0;/**<@brief Flag whether use CMA algorithm.*/
  long int R0 = 0;
  long int W0 = 0;
  long int num_pivot = 0;
  long int ishift_nspin = 0;
  int NumarrayJ = 0;/**<@brief */
  std::complex<double>*** arrayJ;/**<@brief */
  std::complex<double> vecB[3];/**<@brief */
  int** list_6spin_star;/**<@brief */
  int*** list_6spin_pair;/**<@brief */
};/*namespace BoostList*/

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
*/

