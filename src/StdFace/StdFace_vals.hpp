/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**@file
@brief Variables used in the Standard mode.
These variables are passed as a pointer of the structure(StdIntList).
*/
#include <complex>
#include <vector>

namespace StdI {
  /*
   Initial (undefined)
  */
  extern const int NaN_i;
  extern const double pi;
  /*
  Parameters for LATTICE
  */
  extern char lattice[256];
  extern double a;
  extern double length[3];
  extern int W;
  extern int L;
  extern int Height;
  extern double direct[3][3];
  extern int box[3][3];
  extern int rbox[3][3];
  extern int NCell;
  extern int **Cell;
  extern int NsiteUC;
  extern double **tau;
  /*
  Parameters for MODEL
  */
  extern char model[256];
  extern double mu;
  extern std::complex<double> t;
  extern std::complex<double> tp;
  extern std::complex<double> t0;
  extern std::complex<double> t0p;
  extern std::complex<double> t0pp;
  extern std::complex<double> t1;
  extern std::complex<double> t1p;
  extern std::complex<double> t1pp;/**<@brief Anisotropic hopping (3rd), input parameter*/
  extern std::complex<double> t2;/**<@brief Anisotropic hopping (1st), input parameter*/
  extern std::complex<double> t2p;/**<@brief Anisotropic hopping (2nd), input parameter*/
  extern std::complex<double> t2pp;/**<@brief Anisotropic hopping (3rd), input parameter*/
  extern std::complex<double> tpp;/**<@brief 3rd-nearest hopping, input parameter*/
  extern double U;/**<@brief On-site Coulomb potential, input parameter*/
  extern double V;/**<@brief Off-site Coulomb potential (1st), input parameter*/
  extern double Vp;/**<@brief Off-site Coulomb potential (2nd), input parameter*/
  extern double V0;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  extern double V0p;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  extern double V0pp;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  extern double V1;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  extern double V1p;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  extern double V1pp;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  extern double V2;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  extern double V2p;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  extern double V2pp;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  extern double Vpp;/**<@brief Off-site Coulomb potential (3rd), input parameter*/
  /**/
  extern double JAll;
  extern double JpAll;
  extern double J0All;
  extern double J0pAll;
  extern double J0ppAll;
  extern double J1All;
  extern double J1pAll;
  extern double J1ppAll;
  extern double J2All;
  extern double J2pAll;
  extern double J2ppAll;
  extern double JppAll;
  extern double J[3][3];
  extern double Jp[3][3];
  extern double J0[3][3];
  extern double J0p[3][3];
  extern double J0pp[3][3];
  extern double J1[3][3];
  extern double J1p[3][3];
  extern double J1pp[3][3];
  extern double J2[3][3];
  extern double J2p[3][3];
  extern double J2pp[3][3];
  extern double Jpp[3][3];
  extern double D[3][3];
  extern double h;/**<@brief Longitudinal magnetic field, input parameter.*/
  extern double Gamma;/**<@brief Transvars magnetic field, input parameter.*/
  extern double K;
  /*
   Phase for the boundary
  */
  extern const double pi180;/**<@brief @f$\pi/180@f$, set in StdFace_ResetVals().*/
  extern double phase[3];/**<@brief Boundary phase, input parameter phase0, etc.*/
  extern std::complex<double> ExpPhase[3];/**<@brief @f$\exp(i \pi {\rm phase}/180)@f$.*/
  extern int AntiPeriod[3];
  /*
   Transfer, Interaction, Locspin
  */
  extern int nsite;
  extern int *locspinflag;
  extern std::vector<std::vector<int> > transindx;
  extern std::vector<std::complex<double> > trans;
  extern int Lintr;
  extern std::vector<std::vector<int> > intrindx;
  extern std::vector<std::complex<double> > intr;
  extern int LCintra;
  extern std::vector<int> CintraIndx;
  extern std::vector<double> Cintra;
  extern int LCinter;
  extern std::vector<std::vector<int> > CinterIndx;
  extern std::vector<double> Cinter;
  extern int LHund;
  extern std::vector<std::vector<int> > HundIndx;
  extern std::vector<double>Hund;
  extern int LEx;
  extern std::vector<std::vector<int> > ExIndx;
  extern std::vector<double>Ex;
  extern int LPairLift;
  extern std::vector<std::vector<int> > PLIndx;
  extern std::vector<double>PairLift;
  extern int LPairHopp;
  extern std::vector<std::vector<int> >PHIndx;
  extern std::vector<double>PairHopp;
  extern int lBoost;
  /*
   Calculation conditions
  */
  extern int lGC;
  extern int nelec;
  extern int S2;
  extern char outputmode[256];
  extern char CDataFileHead[256];
  extern int Sz2;
  extern int ioutputmode;
  /*
   Wannier90 mode
  */
  extern double cutoff_t;
  extern double cutoff_u;
  extern double cutoff_j;
  extern double cutoff_length_t;
  extern double cutoff_length_U;
  extern double cutoff_length_J;
  extern int cutoff_tR[3];
  extern int cutoff_UR[3];
  extern int cutoff_JR[3];
  extern int double_counting;
#if defined(_HPhi)
  /*
  HPhi modpara
  */
  extern char method[256];/**<@brief The name of method, input from file.*/
  extern char Restart[256];/**<@brief The name of restart mode, input from file.*/
  extern char InitialVecType[256];/**<@brief The name of initialguess-type, input from file.*/
  extern char EigenVecIO[256];/**<@brief The name of I/O mode for eigenvector, input from file*/
  extern int FlgTemp;/**<@brief */
  extern int Lanczos_max;/**<@brief The maxixmum number of iterations, input from file*/
  extern int initial_iv; /**<@brief the number for generating random number, input from file.*/
  extern int nvec;/**<@brief */
  extern int exct;/**<@brief The number of eigenvectors to be computed. input from file*/
  extern int LanczosEps;/**<@brief Convergence threshold for the Lanczos method.*/
  extern int LanczosTarget;/**<@brief Which eigenvector is used for the convergence check.*/
  extern int NumAve;/**<@brief Number of trials for TPQ calculation.*/
  extern int ExpecInterval;
  extern double LargeValue;/**<@brief The shift parameter for the TPQ calculation.*/
  /*
  Boost
  */
  extern int ***list_6spin_pair;/**<@brief */
  extern int **list_6spin_star;/**<@brief */
  extern int num_pivot;/**<@brief */
  extern int ishift_nspin;/**<@brief */
  /*
  Spectrum
  */
  extern char CalcSpec[256];/**<@brief The name of mode for spectrum, input from file.*/
  extern char SpectrumType[256];/**<@brief The type of mode for spectrum, input from file.*/
  extern int Nomega;/**<@brief Number of frequencies, input from file.*/
  extern double OmegaMax;/**<@brief Maximum of frequency for spectrum, input from file.*/
  extern double OmegaMin;/**<@brief Minimum of frequency for spectrum, input from file.*/
  extern double OmegaIm;/**<@brief Imaginary part of frequency.*/
  extern double SpectrumQ[3];/**<@brief wavenumver (q-vector) in fractional coordinate*/
  extern int SpectrumBody;
  /*
  Time evolution
  */
  extern double dt;/**<@brief Time step*/
  extern double tshift;/**<@brief Shift of time-step of laser*/
  extern double tdump;/**<@brief Time scale of dumping*/
  extern double freq;/**<@brief Frequency of laser*/
  extern double Uquench;/**<@brief Quenched on-site potential*/
  extern double VecPot[3];/**<@brief Vector potential*/
  extern char PumpType[256];/**<@brief The type of pump*/
  extern int PumpBody;
  extern int *npump;/**<@brief [StdI::nt] Number of transfer, counted in each lattice file.*/
  extern std::vector<std::vector<std::vector<int> > > pumpindx;
  extern std::vector<std::vector<std::complex<double> > > pump;
  extern double **At;/**<@brief [StdI::nt][3] Vector potential.*/
  extern int ExpandCoef;/**<@brief The number of Hamiltonian-vector operation for the time-evolution*/
#elif defined(_mVMC)
  /*mVMC modpara*/
  extern char CParaFileHead[256];
  extern int NVMCCalMode;
  extern int NLanczosMode;/**<@brief Power Lanczos(=1), input from file*/
  extern int NDataIdxStart;/**<@brief Start index of trials, input from file.*/
  extern int NDataQtySmp;/**<@brief Number of trials, input from file.*/
  extern int NSPGaussLeg;
  extern int NMPTrans;/**<@brief Number of translation symmetry*/
  extern int NSROptItrStep;/**<@brief Number of iterations for stocastic reconfiguration*/
  extern int NSROptItrSmp;/**<@brief Number of steps for sampling*/
  extern int NSROptFixSmp;/**<@brief */
  extern double DSROptRedCut;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  extern double DSROptStaDel;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  extern double DSROptStepDt;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  extern int NVMCWarmUp;/**<@brief */
  extern int NVMCInterval;/**<@brief */
  extern int NVMCSample;/**<@brief */
  extern int NExUpdatePath;/**<@brief */
  extern int RndSeed;/**<@brief */
  extern int NSplitSize;/**<@brief */
  extern int NSPStot;/**<@brief */
  extern int NStore;/**<@brief */
  extern int NSRCG;/**<@brief */
  extern int ComplexType;/**<@brief */
  /*
   Sub-lattice
  */
  extern int Lsub;/**<@brief Sublattice*/
  extern int Wsub;/**<@brief Sublattice*/
  extern int Hsub;/**<@brief Sublattice*/
  extern int NCellsub;/**<@brief Number of cells in a sublattice*/
  extern int boxsub[3][3];/**<@brief Sublattice*/
  extern int rboxsub[3][3];/**<@brief Sublattice*/
  /*
   2-body part of the trial wavefunction
  */
  extern int **Orb;/**<@brief [StdI::nsite][StdI::nsite] Orbital index*/
  extern int **AntiOrb;/**<@brief [StdI::nsite][StdI::nsite] Anti-periodic switch*/
  extern int NOrb;/**<@brief Number of independent orbital index*/
  extern int NSym;
#endif
};
