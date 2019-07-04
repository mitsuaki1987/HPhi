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

namespace StdI {
  /*
   Initial (undefined)
  */
  extern int NaN_i;/**<@brief It is used for initializing input parameter. 
            This means that a parameter wich is not specified in input file.
            Set in StdFace_ResetVals().*/
  extern double pi;/**<@brief @f$\pi=3.14...@f$*/
  /*
  Parameters for LATTICE
  */
  extern char lattice[256];/**<@brief Name of lattice. Input parameter.*/
  extern double a; /**<@brief The lattice constant. Input parameter.*/
  extern double length[3];/**<@brief Anisotropic lattice constant, 
                   input parameter wlength, llength, hlength.*/
  extern int W;/**<@brief Number of sites along the 1st axis, input parameter.*/
  extern int L;/**<@brief Number of sites along the 2nd axis, input parameter.*/
  extern int Height;/**<@brief Number of sites along the 3rd axis, input parameter.*/
  extern double direct[3][3];/**<@brief The unit direct lattice vector. 
                      Set in StdFace_InitSite().*/
  extern int box[3][3];/**<@brief The shape of the super-cell. Input parameter
                a0W, a0L, a0H, etc. or defined from StdI::W, etc. in 
                StdFace_InitSite().*/
  extern int rbox[3][3];/**<@brief The inversion of StdI::box.
                 Set in StdFace_InitSite().*/
  extern int NCell;/**<@brief The number of the unit cell in the super-cell
            (determinant of StdI::box). Set in StdFace_InitSite().*/
  extern int **Cell;/**<@brief [StdIntList][3] The cell position in the fractional 
             coordinate. Malloc and Set in StdFace_InitSite().*/
  extern int NsiteUC;/**<@brief Number of sites in the unit cell. Defined in the
              beginning of each lattice function*/
  extern double **tau;/**<@brief Cell-internal site position in the fractional 
               coordinate. Defined in the beginning of each lattice function*/
  /*
  Parameters for MODEL
  */
  extern char model[256];/**<@brief Name of model, input parameter*/
  extern double mu;/**<@brief Chemical potential, input parameter*/
  extern std::complex<double> t;/**<@brief Nearest-neighbor hopping, input parameter*/
  extern std::complex<double> tp;/**<@brief 2nd-nearest hopping, input parameter*/
  extern std::complex<double> t0;/**<@brief Anisotropic hopping (1st), input parameter*/
  extern std::complex<double> t0p;/**<@brief Anisotropic hopping (2nd), input parameter*/
  extern std::complex<double> t0pp;/**<@brief Anisotropic hopping (3rd), input parameter*/
  extern std::complex<double> t1;/**<@brief Anisotropic hopping (1st), input parameter*/
  extern std::complex<double> t1p;/**<@brief Anisotropic hopping (2nd), input parameter*/
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
  extern double JAll;/**<@brief Isotropic, diagonal spin coupling (1st Near.), 
              input parameter J.*/
  extern double JpAll;/**<@brief Isotropic, diagonal spin coupling (2nd Near), 
               input parameter Jp.*/
  extern double J0All;/**<@brief Anisotropic, diagonal spin coupling (1st Near), 
               input parameter J0.*/
  extern double J0pAll;/**<@brief Anisotropic, diagonal spin coupling (2nd Near), 
               input parameter J0'.*/
  extern double J0ppAll;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J0''.*/
  extern double J1All;/**<@brief Anisotropic, diagonal spin coupling (1st Near),
               input parameter J1.*/
  extern double J1pAll;/**<@brief Anisotropic, diagonal spin coupling (2nd Near), 
               input parameter J1'.*/
  extern double J1ppAll;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J1''.*/
  extern double J2All;/**<@brief Anisotropic, diagonal spin coupling (1st Near),
               input parameter J2.*/
  extern double J2pAll;/**<@brief Anisotropic, diagonal spin coupling (2nd Near), 
               input parameter J2'.*/
  extern double J2ppAll;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J2''.*/
  extern double JppAll;/**<@brief Isotropic, diagonal spin coupling (3rd Near),
               input parameter J''.*/
  extern double J[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                 (1st Near.), input parameter Jx, Jy, Jz, Jxy, etc.*/
  extern double Jp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J'x, J'y, J'z, J'xy, etc.*/
  extern double J0[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                 (1st Near.), input parameter J0x, J0y, J0z, J0xy, etc. 
                 or set in StdFace_InputSpinNN().*/
  extern double J0p[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J0'x, J0'y, J0'z, J0'xy, etc. 
                   or set in StdFace_InputSpin().*/
  extern double J0pp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J0''x, J0''y, J0''z, J0''xy, etc.
                   or set in StdFace_InputSpin().*/
  extern double J1[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                  (1st Near.), input parameter J1x, J1y, J1z, J1xy, etc. 
                  or set in StdFace_InputSpinNN().*/
  extern double J1p[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J1'x, J1'y, J1'z, J1'xy, etc. 
                   or set in StdFace_InputSpin().*/
  extern double J1pp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J1''x, J1''y, J1''z, J1''xy, etc.
                   or set in StdFace_InputSpin().*/
  extern double J2[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                  (1st Near.), input parameter J2x, J2y, J2z, J2xy, etc. 
                  or set in StdFace_InputSpinNN().*/
  extern double J2p[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J2'x, J2'y, J2'z, J2'xy, etc. 
                   or set in StdFace_InputSpin().*/
  extern double J2pp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J2''x, J2''y, J2''z, J2''xy, etc.
                   or set in StdFace_InputSpin().*/
  extern double Jpp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J''x, J''y, J''z, J''xy, etc.*/
  extern double D[3][3];/**<@brief Coefficient for @f${\hat S}_{i z} {\hat S}_{i z}@f$
                 input parameter D. Only D[2][2] is used.*/
  extern double h;/**<@brief Longitudinal magnetic field, input parameter.*/
  extern double Gamma;/**<@brief Transvars magnetic field, input parameter.*/
  extern double K;/**<@brief 4-spin term. Not used.*/
  /*
   Phase for the boundary
  */
  extern double pi180;/**<@brief @f$\pi/180@f$, set in StdFace_ResetVals().*/
  extern double phase[3];/**<@brief Boundary phase, input parameter phase0, etc.*/
  extern std::complex<double> ExpPhase[3];/**<@brief @f$\exp(i \pi {\rm phase}/180)@f$.*/
  extern int AntiPeriod[3];/**<@brief If corresponding StdI::phase = 180,
                    it becomes 1.*/
  /*
   Transfer, Interaction, Locspin
  */
  extern int nsite;/**<@brief Number of sites, set in the each lattice file.*/
  extern int *locspinflag;/**<@brief [StdI::nsite] LocSpin in Expert mode, 
                   malloc and set in each lattice file.*/
  extern int ntrans;/**<@brief Number of transfer, counted in each lattice file.*/
  extern int **transindx;/**<@brief [StdI::ntrans][4] Site/spin indices of 
                  one-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_trans().*/
  extern std::complex<double> *trans;/**<@brief [StdI::ntrans] Coefficient of 
                  one-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_trans().*/
  extern int nintr;/**<@brief Number of InterAll, counted in each lattice file.*/
  extern int Lintr;/**<@brief Print interall.def or not, set in PrintInteractions().*/
  extern int **intrindx;/**<@brief [StdI::nintr][8] Site/spin indices of 
                  two-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern std::complex<double> *intr;/**<@brief [StdI::nintr] Coefficient of general
                  two-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern int NCintra;/**<@brief Number of intra-site Coulomb interaction,
              counted in each lattice file.*/
  extern int LCintra;/**<@brief Print coulombintra.def or not, set in PrintInteractions().*/
  extern int **CintraIndx;/**<@brief [StdI::NCintra][1] Site indices of 
                  intra-site Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern double *Cintra;/**<@brief [StdI::NCintra] Coefficient of intra-site
                  Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern int NCinter;/**<@brief Number of inter-site Coulomb interaction,
              counted in each lattice file.*/
  extern int LCinter;/**<@brief Print coulombinter.def or not, set in PrintInteractions().*/
  extern int **CinterIndx;/**<@brief [StdI::NCinter][2] Site indices of 
                  inter-site Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern double *Cinter;/**<@brief [StdI::NCinter] Coefficient of inter-site
                  Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern int NHund;/**<@brief Number of Hund term, counted in each lattice file.*/
  extern int LHund;/**<@brief Print hund.def or not, set in PrintInteractions().*/
  extern int **HundIndx;/**<@brief [StdI::NHund][2] Site indices of 
                  Hund term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern double *Hund;/**<@brief [StdI::NHund] Coefficient of Hund term, 
               malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  extern int NEx;/**<@brief Number of exchange term, counted in each lattice file.*/
  extern int LEx;/**<@brief Print exchange.def or not, set in PrintInteractions().*/
  extern int **ExIndx;/**<@brief [StdI::NEx][2] Site indices of 
                  exchange term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern double *Ex;/**<@brief [StdI::NEx] Coefficient of exchange term, 
               malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  extern int NPairLift;/**<@brief Number of pair-lift term, counted in each lattice file.*/
  extern int LPairLift;/**<@brief Print pairlift.def or not, set in PrintInteractions().*/
  extern int **PLIndx;/**<@brief [StdI::NPairLift][2] Site indices of 
                  pair-lift term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  extern double *PairLift;/**<@brief [StdI::NPairLift] Coefficient of 
                   pair-lift term, malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  extern int NPairHopp;/**<@brief Number of pair-hopping term, counted in each lattice file.*/
  extern int LPairHopp;/**<@brief Print pairhopp.def or not, set in PrintInteractions().*/
  extern int **PHIndx;/**<@brief [StdI::NPairLift][2] Site indices of
               pair-hopping term, malloc in StdFace_MallocInteractions()
               and set in StdFace_intr().*/
  extern double *PairHopp;/**<@brief [StdI::NPairLift] Coefficient of
                   pair-hopping term, malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  extern int lBoost;
  /*
   Calculation conditions
  */
  extern int lGC;/**<@brief Switch for computing Grandcanonical ensemble(== 1).
          Setted in StdFace_main() after all keywords are read.*/
  extern int nelec;/**<@brief Number of electrons, input from file.*/
  extern int S2;/**<@brief Total spin |S| of a local spin, input from file.*/
  extern char outputmode[256];/**<@brief Select amount of correlation function,
                       input from file.*/
  extern char CDataFileHead[256];/**<@brief Header of the output files.
                          Input from file*/
  extern int Sz2;/**<@brief Total Sz, input from file.*/
  extern int ioutputmode;/**<@brief Switch associated to StdI::outputmode*/
  /*
   Wannier90 mode
  */
  extern double cutoff_t;/**<@brief Cutoof for the hopping in wannier90, input from file*/
  extern double cutoff_u;/**<@brief Cutoof for the Coulomb in wannier90, input from file*/
  extern double cutoff_j;/**<@brief Cutoof for the Hund in wannier90, input from file*/
  extern double cutoff_length_t; /**<@brief Cutoof for R in wannier90, input from file.*/
  extern double cutoff_length_U; /**<@brief Cutoof for R in wannier90, input from file.*/
  extern double cutoff_length_J; /**<@brief Cutoof for R in wannier90, input from file.*/
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
  extern int ExpecInterval;/**<@brief Interval for the iteration when the expectation 
                    value is computed.*/
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
  extern int SpectrumBody;/**<@brief one- or two-body excitation, defined from
                   StdI::SpectrumType*/
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
  extern int PumpBody;/**<@brief one- or two-body pumping, defined from
                   StdI::PumpType*/
  extern int *npump;/**<@brief [StdI::nt] Number of transfer, counted in each lattice file.*/
  extern int ***pumpindx;/**<@brief [StdI::nt][StdI::npump][4] Site/spin indices of
                  one-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_trans().*/
  extern std::complex<double> **pump;/**<@brief [StdI::nt][StdI::npump] Coefficient of
                        one-body term, malloc in StdFace_MallocInteractions()
                        and set in StdFace_trans().*/
  extern double **At;/**<@brief [StdI::nt][3] Vector potential.*/
  extern int ExpandCoef;/**<@brief The number of Hamiltonian-vector operation for the time-evolution*/
#elif defined(_mVMC)
  /*mVMC modpara*/
  extern char CParaFileHead[256];/**<@brief Header of the optimized wavefunction,
                          input from file*/
  extern int NVMCCalMode;/**<@brief Optimization(=0) or compute correlation
                  function(=1), input from file.*/
  extern int NLanczosMode;/**<@brief Power Lanczos(=1), input from file*/
  extern int NDataIdxStart;/**<@brief Start index of trials, input from file.*/
  extern int NDataQtySmp;/**<@brief Number of trials, input from file.*/
  extern int NSPGaussLeg;/**<@brief Number of Gauss-Legendre points for spin projection,
                  input from file.*/
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
  extern int NSym;/**<@brief Number of translation symmetries, 
           Defined from the number of cells in the sub-lattice.*/
#endif
};
