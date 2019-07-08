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
  extern int nvec;
  extern int k_exct;/**<@brief Read from Calcmod in readdef.h*/
  extern int LanczosEps;
  extern int LanczosTarget;
  extern int read_hacker;
  extern int READ;
  extern int WRITE;
  extern int Nsite;/**<@brief Number of sites in the INTRA process region*/
  extern int NsiteMPI;/**<@brief Total number of sites, differ from DefineList::Nsite*/
  extern int Nup;/**<@brief Number of spin-up electrons in this process. */
  extern int Ndown;/**<@brief Number of spin-down electrons in this process. */
  extern int NupMPI;
  extern int NdownMPI;
  extern int NupOrg;
  extern int NdownOrg;
  extern int Total2Sz;/**<@brief Total @f$2S_z@f$ in this process.*/
  extern int Total2SzMPI;/**<@brief Total @f$2S_z@f$ across processes.*/
  extern int Ne;/**<@brief Number of electrons in this process.*/
  extern int NeMPI;
  extern int Lanczos_max;
  extern double eps_Lanczos;
  extern int Lanczos_restart;
  extern long int initial_iv;
  extern int istep;/**<@brief Index of TPQ step ???*/
  extern int irand;/**<@brief Input keyword TargetTPQRand ???*/
  extern int St;
  extern int* LocSpn;
  extern int NLocSpn;/**<@brief Number of local spins*/
  extern int NCond;/**<@brief Number of itinerant electrons*/
  extern int iFlgGeneralSpin;/**<@brief Flag for the general (Sz/=1/2) spin*/
  extern int iFlgSzConserved;
  extern int fidx;/**<@brief Always 0, it is not used ???*/
  extern long int* Tpow;
  extern long int* OrgTpow;
  extern long int* SiteToBit;
  extern int EDNChemi;/**<@brief Number of on-site term.*/
  extern int* EDChemi;
  extern int* EDSpinChemi;/**<@brief [DefineList::Nsite]*/
  extern double* EDParaChemi; 
  extern int NTransfer;
  extern int EDNTransfer;
  extern int** GeneralTransfer;  
  extern int** EDGeneralTransfer; 
  extern std::complex<double>* ParaGeneralTransfer; 
  extern std::complex<double>* EDParaGeneralTransfer;  
  extern int NCoulombIntra;
  extern int** CoulombIntra;
  extern double* ParaCoulombIntra;
  extern int NCoulombInter;
  extern int** CoulombInter;
  extern double* ParaCoulombInter;
  extern int NHundCoupling;
  extern int** HundCoupling;
  extern double* ParaHundCoupling;
  extern int NPairHopping;
  extern int** PairHopping;
  extern double* ParaPairHopping;
  extern int NExchangeCoupling;
  extern int** ExchangeCoupling;
  extern double* ParaExchangeCoupling;
  extern int NIsingCoupling;
  extern int NPairLiftCoupling;
  extern int** PairLiftCoupling;
  extern double* ParaPairLiftCoupling;
  extern int** InterAll;
  extern int** InterAll_OffDiagonal;
  extern int** InterAll_Diagonal;
  extern int NInterAll;
  extern int NInterAll_Diagonal;
  extern int NInterAll_OffDiagonal;
  extern std::complex<double>* ParaInterAll;
  extern double* ParaInterAll_Diagonal;
  extern std::complex<double>* ParaInterAll_OffDiagonal;
  extern int** CisAjt;
  extern int NCisAjt;
  extern int** CisAjtCkuAlvDC;
  extern int NCisAjtCkuAlvDC;
  extern int*** SingleExcitationOperator;
  extern int NNSingleExcitationOperator;
  extern int* NSingleExcitationOperator;
  extern std::complex<double>** ParaSingleExcitationOperator;
  extern int*** PairExcitationOperator;
  extern int NNPairExcitationOperator;
  extern int* NPairExcitationOperator;
  extern std::complex<double>** ParaPairExcitationOperator;
  extern int iCalcType;
  extern int iCalcEigenVec;
  extern int iInitialVecType;
  extern int iFlgFiniteTemperature;
  extern int iCalcModel;
  extern int iOutputMode;
  extern int iOutputEigenVec;
  extern int iInputEigenVec;
  extern int iOutputHam;
  extern int iInputHam;
  extern int iOutputExVec; 
  extern std::complex<double> dcOmegaMax;
  extern std::complex<double> dcOmegaMin;
  extern std::complex<double> dcOmegaOrg;
  extern int iNOmega;
  extern int iFlgSpecOmegaMax;
  extern int iFlgSpecOmegaMin;
  extern int iFlgSpecOmegaOrg;
  extern int iFlgCalcSpec;
  extern int iFlagListModified;
  extern int iReStart;
  extern int iFlgMPI;
  extern int iNGPU;
  extern int iFlgScaLAPACK;
  extern int NTETimeSteps;
  extern double* TETime;
  extern int NLaser;
  extern double* ParaLaser;
  extern int NTETransferMax;
  extern int* NTETransfer;
  extern int* NTETransferDiagonal;
  extern int*** TETransfer;
  extern int*** TETransferDiagonal;
  extern std::complex<double>** ParaTETransfer;
  extern double** ParaTETransferDiagonal; 
  extern int NTEInterAllMax;
  extern int* NTEInterAll; 
  extern int* NTEInterAllOffDiagonal; 
  extern int* NTEInterAllDiagonal;
  extern int*** TEInterAll;  
  extern int*** TEInterAllOffDiagonal;
  extern int*** TEInterAllDiagonal; 
  extern std::complex<double>** ParaTEInterAll;
  extern std::complex<double>** ParaTEInterAllOffDiagonal; 
  extern double** ParaTEInterAllDiagonal; 
  extern int** TEChemi;
  extern int* NTEChemi; 
  extern int** SpinTEChemi; 
  extern double** ParaTEChemi; 
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

namespace Wave {
  extern std::complex<double>** v0;  /**< A vector after multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
  extern std::complex<double>** v1;  /**< A vector before multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
  extern std::complex<double>** v1buf; /**< A temporary vector for MPI. */
}

namespace List {
  extern double* Diagonal; /**< list for diagonal components.*/
  extern long int* c1; /**< list of getting real-space configuration for canonical state*/
  extern long int* c1buf;/**< list of getting real-space configuration for canonical state across processes*/
  extern long int* c2_1;/**< list to get index of list_1*/
  extern long int* c2_2;/**< list to get index of list_1*/

  /*[s] For Spectrum */
  extern long int* c1_org; /**< list of getting real-space configuration for canonical state before excitation*/
  extern long int* c1buf_org;/**< list of getting real-space configuration for canonical state before excitation across processes*/
  extern long int* c2_1_org;/**< list to get index of list_1_org*/
  extern long int* c2_2_org;/**< list to get index of list_1_org*/
  /*[e] For Spectrum */
}

namespace Step {
  extern double LargeValue;
  extern int    NumAve;
  extern int step_i;
  extern double* global_norm;
  extern double* global_1st_norm;
  extern int step_spin;/**< output step for TE calculation.*/
}

/*[s] For All Diagonalization*/
#ifdef _SCALAPACK
extern std::complex<double> *Z_vec; /**> distributed matrix of eigen vector*/
extern int descZ_vec[9]; /*descriptor for Z_vec*/
#endif
/*
 Variables for the MPI parallelism
*/
namespace MP {
  extern int nproc;//!< Number of processors, defined in InitializeMPI()
  extern int myrank;//!< Process ID, defined in InitializeMPI()
  extern int nthreads;//!< Number of Threads, defined in InitializeMPI()
  extern FILE* STDOUT;/**<@brief File pointer to the standard output
                  defined in InitializeMPI()*/
}

#endif /* HPHI_GLOBAL_H */
