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
@brief Read Input file and write files for Expert mode.
       Initialize variables.
       Check parameters.

The following lattices are supported:
- 1D Chain : StdFace::Chain()
- 1D Ladder : StdFace::Ladder()
- 2D Tetragonal : StdFace::Tetragonal()
- 2D Triangular : StdFace::Triangular()
- 2D Honeycomb : StdFace::Honeycomb()
- 2D Kagome : StdFace::Kagome()
- 3D Simple Orthorhombic : StdFace::Orthorhombic()
- 3D Face Centered Orthorhombic : StdFace::FCOrtho()
- 3D Pyrochlore : StdFace::Pyrochlore()

*/
#include "StdFace_main.hpp"
#include "StdFace_vals.hpp"
#include "StdFace_ModelUtil.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <iostream>

namespace StdI {
  /*
   Initial (undefined)
  */
  const int NaN_i = 2147483647;/**<@brief It is used for initializing input parameter.
            This means that a parameter wich is not specified in input file.
            Set in StdFace::ResetVals().*/
  const double pi = acos(-1.0);/**<@brief @f$\pi=3.14...@f$*/
  const double NaN_d = 0.0 / 0.0;
  /*
  Parameters for LATTICE
  */
  char lattice[256] = "****\0";/**<@brief Name of lattice. Input parameter.*/
  double a = NaN_d; /**<@brief The lattice constant. Input parameter.*/
  double length[3] = { NaN_d ,NaN_d ,NaN_d };
  /**<@brief Anisotropic lattice constant,
                   input parameter wlength, llength, hlength.*/
  int W = NaN_i;/**<@brief Number of sites along the 1st axis, input parameter.*/
  int L = NaN_i;/**<@brief Number of sites along the 2nd axis, input parameter.*/
  int Height = NaN_i;/**<@brief Number of sites along the 3rd axis, input parameter.*/
  double direct[3][3] = {{ NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d }};
  /**<@brief The unit direct lattice vector.
                      Set in StdFace::InitSite().*/
  int box[3][3] = { { NaN_i ,NaN_i ,NaN_i },
                    { NaN_i, NaN_i, NaN_i },
                    { NaN_i, NaN_i, NaN_i } };
  /**<@brief The shape of the super-cell. Input parameter
                a0W, a0L, a0H, etc. or defined from StdI::W, etc. in
                StdFace::InitSite().*/
  int rbox[3][3];/**<@brief The inversion of StdI::box.
                 Set in StdFace::InitSite().*/
  int NCell;/**<@brief The number of the unit cell in the super-cell
            (determinant of StdI::box). Set in StdFace::InitSite().*/
  int** Cell;/**<@brief [StdIntList][3] The cell position in the fractional
             coordinate. Malloc and Set in StdFace::InitSite().*/
  int NsiteUC;/**<@brief Number of sites in the unit cell. Defined in the
              beginning of each lattice function*/
  double** tau;/**<@brief Cell-internal site position in the fractional
               coordinate. Defined in the beginning of each lattice function*/
               /*
               Parameters for MODEL
               */
  char model[256] = "****\0";/**<@brief Name of model, input parameter*/
  double mu = NaN_d;/**<@brief Chemical potential, input parameter*/
  std::complex<double> t = NaN_d;/**<@brief Nearest-neighbor hopping, input parameter*/
  std::complex<double> tp = NaN_d;/**<@brief 2nd-nearest hopping, input parameter*/
  std::complex<double> t0 = NaN_d;/**<@brief Anisotropic hopping (1st), input parameter*/
  std::complex<double> t0p = NaN_d;/**<@brief Anisotropic hopping (2nd), input parameter*/
  std::complex<double> t0pp = NaN_d;/**<@brief Anisotropic hopping (3rd), input parameter*/
  std::complex<double> t1 = NaN_d;/**<@brief Anisotropic hopping (1st), input parameter*/
  std::complex<double> t1p = NaN_d;/**<@brief Anisotropic hopping (2nd), input parameter*/
  std::complex<double> t1pp = NaN_d;/**<@brief Anisotropic hopping (3rd), input parameter*/
  std::complex<double> t2 = NaN_d;/**<@brief Anisotropic hopping (1st), input parameter*/
  std::complex<double> t2p = NaN_d;/**<@brief Anisotropic hopping (2nd), input parameter*/
  std::complex<double> t2pp = NaN_d;/**<@brief Anisotropic hopping (3rd), input parameter*/
  std::complex<double> tpp = NaN_d;/**<@brief 3rd-nearest hopping, input parameter*/
  double U = NaN_d;/**<@brief On-site Coulomb potential, input parameter*/
  double V = NaN_d;/**<@brief Off-site Coulomb potential (1st), input parameter*/
  double Vp = NaN_d;/**<@brief Off-site Coulomb potential (2nd), input parameter*/
  double V0 = NaN_d;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  double V0p = NaN_d;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  double V0pp = NaN_d;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  double V1 = NaN_d;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  double V1p = NaN_d;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  double V1pp = NaN_d;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  double V2 = NaN_d;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  double V2p = NaN_d;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  double V2pp = NaN_d;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  double Vpp = NaN_d;/**<@brief Off-site Coulomb potential (3rd), input parameter*/
  /**/
  double JAll = NaN_d;/**<@brief Isotropic, diagonal spin coupling (1st Near.),
              input parameter J.*/
  double JpAll = NaN_d;/**<@brief Isotropic, diagonal spin coupling (2nd Near),
               input parameter Jp.*/
  double J0All = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (1st Near),
               input parameter J0.*/
  double J0pAll = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (2nd Near),
               input parameter J0'.*/
  double J0ppAll = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J0''.*/
  double J1All = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (1st Near),
               input parameter J1.*/
  double J1pAll = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (2nd Near),
               input parameter J1'.*/
  double J1ppAll = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J1''.*/
  double J2All = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (1st Near),
               input parameter J2.*/
  double J2pAll = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (2nd Near),
               input parameter J2'.*/
  double J2ppAll = NaN_d;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J2''.*/
  double JppAll = NaN_d;/**<@brief Isotropic, diagonal spin coupling (3rd Near),
               input parameter J''.*/
  double J[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                     { NaN_d, NaN_d, NaN_d },
                     { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                 (1st Near.), input parameter Jx, Jy, Jz, Jxy, etc.*/
  double Jp[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J'x, J'y, J'z, J'xy, etc.*/
  double J0[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                 (1st Near.), input parameter J0x, J0y, J0z, J0xy, etc.
                 or set in StdFace::InputSpinNN().*/
  double J0p[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J0'x, J0'y, J0'z, J0'xy, etc.
                   or set in StdFace::InputSpin().*/
  double J0pp[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J0''x, J0''y, J0''z, J0''xy, etc.
                   or set in StdFace::InputSpin().*/
  double J1[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                      { NaN_d, NaN_d, NaN_d },
                      { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                  (1st Near.), input parameter J1x, J1y, J1z, J1xy, etc.
                  or set in StdFace::InputSpinNN().*/
  double J1p[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J1'x, J1'y, J1'z, J1'xy, etc.
                   or set in StdFace::InputSpin().*/
  double J1pp[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J1''x, J1''y, J1''z, J1''xy, etc.
                   or set in StdFace::InputSpin().*/
  double J2[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                  (1st Near.), input parameter J2x, J2y, J2z, J2xy, etc.
                  or set in StdFace::InputSpinNN().*/
  double J2p[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J2'x, J2'y, J2'z, J2'xy, etc.
                   or set in StdFace::InputSpin().*/
  double J2pp[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J2''x, J2''y, J2''z, J2''xy, etc.
                   or set in StdFace::InputSpin().*/
  double Jpp[3][3] = { { NaN_d ,NaN_d ,NaN_d },
                         { NaN_d, NaN_d, NaN_d },
                         { NaN_d, NaN_d, NaN_d } };
  /**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J''x, J''y, J''z, J''xy, etc.*/
  double D[3][3] = { {  0.0 , 0.0 , 0.0 },
                     {  0.0,  0.0, 0.0 },
                     {  0.0,  0.0, NaN_d } };
  /**<@brief Coefficient for @f${\hat S}_{i z} {\hat S}_{i z}@f$
                 input parameter D. Only D[2][2] is used.*/
  double h = NaN_d;/**<@brief Longitudinal magnetic field, input parameter.*/
  double Gamma = NaN_d;/**<@brief Transvars magnetic field, input parameter.*/
  double K = NaN_d;/**<@brief 4-spin term. Not used.*/
  /*
   Phase for the boundary
  */
  const double pi180 = StdI::pi / 180.0;/**<@brief @f$\pi/180@f$, set in StdFace::ResetVals().*/
  double phase[3] = { NaN_d, NaN_d, NaN_d };/**<@brief Boundary phase, input parameter phase0, etc.*/
  std::complex<double> ExpPhase[3];/**<@brief @f$\exp(i \pi {\rm phase}/180)@f$.*/
  int AntiPeriod[3];/**<@brief If corresponding StdI::phase = 180,
                    it becomes 1.*/
                    /*
                     Transfer, Interaction, Locspin
                    */
  int nsite;/**<@brief Number of sites, set in the each lattice file.*/
  int* locspinflag;/**<@brief [StdI::nsite] LocSpin in Expert mode,
                   malloc and set in each lattice file.*/
  std::vector<std::vector<int> > transindx;/**<@brief [StdI::ntrans][4] Site/spin indices of
                  one-body term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::trans().*/
  std::vector<std::complex<double> > trans;/**<@brief [StdI::ntrans] Coefficient of
                  one-body term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::trans().*/
  int Lintr;/**<@brief Print interall.def or not, set in PrintInteractions().*/
  std::vector<std::vector<int> > intrindx;/**<@brief [StdI::nintr][8] Site/spin indices of
                  two-body term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  std::vector<std::complex<double> > intr;/**<@brief [StdI::nintr] Coefficient of general
                  two-body term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  int LCintra;/**<@brief Print coulombintra.def or not, set in PrintInteractions().*/
  std::vector<int> CintraIndx;/**<@brief [StdI::NCintra][1] Site indices of
                  intra-site Coulomb term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  std::vector<double> Cintra;/**<@brief [StdI::NCintra] Coefficient of intra-site
                  Coulomb term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  int LCinter;/**<@brief Print coulombinter.def or not, set in PrintInteractions().*/
  std::vector<std::vector<int> > CinterIndx;/**<@brief [StdI::NCinter][2] Site indices of
                  inter-site Coulomb term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  std::vector<double> Cinter;/**<@brief [StdI::NCinter] Coefficient of inter-site
                  Coulomb term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  int LHund;/**<@brief Print hund.def or not, set in PrintInteractions().*/
  std::vector<std::vector<int> > HundIndx;/**<@brief [StdI::NHund][2] Site indices of
                  Hund term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  std::vector<double> Hund;/**<@brief [StdI::NHund] Coefficient of Hund term,
               malloc in StdFace::MallocInteractions()
                   and set in StdFace::intr().*/
  int LEx;/**<@brief Print exchange.def or not, set in PrintInteractions().*/
  std::vector<std::vector<int> > ExIndx;/**<@brief [StdI::NEx][2] Site indices of
                  exchange term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  std::vector<double> Ex;/**<@brief [StdI::NEx] Coefficient of exchange term,
               malloc in StdFace::MallocInteractions()
                   and set in StdFace::intr().*/
  int LPairLift;/**<@brief Print pairlift.def or not, set in PrintInteractions().*/
  std::vector<std::vector<int> > PLIndx;/**<@brief [StdI::NPairLift][2] Site indices of
                  pair-lift term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::intr().*/
  std::vector<double> PairLift;/**<@brief [StdI::NPairLift] Coefficient of
                   pair-lift term, malloc in StdFace::MallocInteractions()
                   and set in StdFace::intr().*/
  int LPairHopp;/**<@brief Print pairhopp.def or not, set in PrintInteractions().*/
  std::vector<std::vector<int> > PHIndx;/**<@brief [StdI::NPairLift][2] Site indices of
               pair-hopping term, malloc in StdFace::MallocInteractions()
               and set in StdFace::intr().*/
  std::vector<double> PairHopp;/**<@brief [StdI::NPairLift] Coefficient of
                   pair-hopping term, malloc in StdFace::MallocInteractions()
                   and set in StdFace::intr().*/
  int lBoost;
  /*
   Calculation conditions
  */
  int lGC;/**<@brief Switch for computing Grandcanonical ensemble(== 1).
          Setted in StdFace::main() after all keywords are read.*/
  int nelec = NaN_i;/**<@brief Number of electrons, input from file.*/
  int S2 = NaN_i;/**<@brief Total spin |S| of a local spin, input from file.*/
  char outputmode[256] = "****\0";/**<@brief Select amount of correlation function,
                       input from file.*/
  char CDataFileHead[256] = "****\0";/**<@brief Header of the output files.
                          Input from file*/
  int Sz2 = NaN_i;/**<@brief Total Sz, input from file.*/
  int ioutputmode;/**<@brief Switch associated to StdI::outputmode*/
  /*
   Wannier90 mode
  */
  double cutoff_t = NaN_d;/**<@brief Cutoof for the hopping in wannier90, input from file*/
  double cutoff_u = NaN_d;/**<@brief Cutoof for the Coulomb in wannier90, input from file*/
  double cutoff_j = NaN_d;/**<@brief Cutoof for the Hund in wannier90, input from file*/
  double cutoff_length_t = NaN_d; /**<@brief Cutoof for R in wannier90, input from file.*/
  double cutoff_length_U = NaN_d; /**<@brief Cutoof for R in wannier90, input from file.*/
  double cutoff_length_J = NaN_d; /**<@brief Cutoof for R in wannier90, input from file.*/
  int cutoff_tR[3] = { NaN_i ,NaN_i ,NaN_i };
  int cutoff_UR[3] = { NaN_i ,NaN_i ,NaN_i };
  int cutoff_JR[3] = { NaN_i ,NaN_i ,NaN_i };
  int double_counting = NaN_i;
#if defined(_HPhi)
  /*
  HPhi modpara
  */
  char method[256] = "****\0";/**<@brief The name of method, input from file.*/
  char Restart[256] = "****\0";/**<@brief The name of restart mode, input from file.*/
  char InitialVecType[256] = "****\0";/**<@brief The name of initialguess-type, input from file.*/
  char EigenVecIO[256] = "****\0";/**<@brief The name of I/O mode for eigenvector, input from file*/
  int FlgTemp = 1;/**<@brief */
  int Lanczos_max = NaN_i;/**<@brief The maxixmum number of iterations, input from file*/
  int initial_iv = NaN_i; /**<@brief the number for generating random number, input from file.*/
  int nvec = NaN_i;/**<@brief */
  int exct = NaN_i;/**<@brief The number of eigenvectors to be computed. input from file*/
  int LanczosEps = NaN_i;/**<@brief Convergence threshold for the Lanczos method.*/
  int LanczosTarget = NaN_i;/**<@brief Which eigenvector is used for the convergence check.*/
  int NumAve = NaN_i;/**<@brief Number of trials for TPQ calculation.*/
  int ExpecInterval = NaN_i;/**<@brief Interval for the iteration when the expectation
                    value is computed.*/
  double LargeValue = NaN_d;/**<@brief The shift parameter for the TPQ calculation.*/
  /*
  Boost
  */
  int*** list_6spin_pair;/**<@brief */
  int** list_6spin_star;/**<@brief */
  int num_pivot;/**<@brief */
  int ishift_nspin;/**<@brief */
  /*
  Spectrum
  */
  char CalcSpec[256] = "****\0";/**<@brief The name of mode for spectrum, input from file.*/
  char SpectrumType[256] = "****\0";/**<@brief The type of mode for spectrum, input from file.*/
  int Nomega= NaN_i;/**<@brief Number of frequencies, input from file.*/
  double OmegaMax = NaN_d;/**<@brief Maximum of frequency for spectrum, input from file.*/
  double OmegaMin = NaN_d;/**<@brief Minimum of frequency for spectrum, input from file.*/
  double OmegaIm = NaN_d;/**<@brief Imaginary part of frequency.*/
  double SpectrumQ[3] = { NaN_d ,NaN_d ,NaN_d };
  /**<@brief wavenumver (q-vector) in fractional coordinate*/
  int SpectrumBody;/**<@brief one- or two-body excitation, defined from
                   StdI::SpectrumType*/
                   /*
                   Time evolution
                   */
  double dt = NaN_d;/**<@brief Time step*/
  double tshift = NaN_d;/**<@brief Shift of time-step of laser*/
  double tdump = NaN_d;/**<@brief Time scale of dumping*/
  double freq = NaN_d;/**<@brief Frequency of laser*/
  double Uquench = NaN_d;/**<@brief Quenched on-site potential*/
  double VecPot[3] = { NaN_d ,NaN_d ,NaN_d };/**<@brief Vector potential*/
  char PumpType[256] = "****\0";/**<@brief The type of pump*/
  int PumpBody;/**<@brief one- or two-body pumping, defined from
                   StdI::PumpType*/
  int* npump;/**<@brief [StdI::nt] Number of transfer, counted in each lattice file.*/
  std::vector<std::vector<std::vector<int> > > pumpindx;/**<@brief [StdI::nt][StdI::npump][4] Site/spin indices of
                  one-body term, malloc in StdFace::MallocInteractions()
                  and set in StdFace::trans().*/
  std::vector<std::vector<std::complex<double> > > pump;/**<@brief [StdI::nt][StdI::npump] Coefficient of
                        one-body term, malloc in StdFace::MallocInteractions()
                        and set in StdFace::trans().*/
  double** At;/**<@brief [StdI::nt][3] Vector potential.*/
  int ExpandCoef = NaN_i;/**<@brief The number of Hamiltonian-vector operation for the time-evolution*/
#elif defined(_mVMC)
  /*mVMC modpara*/
  char CParaFileHead[256] = "****\0";/**<@brief Header of the optimized wavefunction,
                          input from file*/
  int NVMCCalMode = NaN_i;/**<@brief Optimization(=0) or compute correlation
                  function(=1), input from file.*/
  int NLanczosMode = NaN_i;/**<@brief Power Lanczos(=1), input from file*/
  int NDataIdxStart = NaN_i;/**<@brief Start index of trials, input from file.*/
  int NDataQtySmp = NaN_i;/**<@brief Number of trials, input from file.*/
  int NSPGaussLeg = NaN_i;/**<@brief Number of Gauss-Legendre points for spin projection,
                  input from file.*/
  int NMPTrans = NaN_i;/**<@brief Number of translation symmetry*/
  int NSROptItrStep = NaN_i;/**<@brief Number of iterations for stocastic reconfiguration*/
  int NSROptItrSmp = NaN_i;/**<@brief Number of steps for sampling*/
  int NSROptFixSmp = NaN_i;/**<@brief */
  double DSROptRedCut = NaN_d;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  double DSROptStaDel = NaN_d;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  double DSROptStepDt = NaN_d;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  int NVMCWarmUp = NaN_i;/**<@brief */
  int NVMCInterval = NaN_i;/**<@brief */
  int NVMCSample = NaN_i;/**<@brief */
  int NExUpdatePath = NaN_i;/**<@brief */
  int RndSeed = NaN_i;/**<@brief */
  int NSplitSize = NaN_i;/**<@brief */
  int NSPStot = NaN_i;/**<@brief */
  int NStore = NaN_i;/**<@brief */
  int NSRCG = NaN_i;/**<@brief */
  int ComplexType = NaN_i;/**<@brief */
  /*
   Sub-lattice
  */
  int Lsub = NaN_i;/**<@brief Sublattice*/
  int Wsub = NaN_i;/**<@brief Sublattice*/
  int Hsub = NaN_i;/**<@brief Sublattice*/
  int NCellsub = NaN_i;/**<@brief Number of cells in a sublattice*/
  int boxsub[3][3] = { { NaN_i ,NaN_i ,NaN_i },
                       { NaN_i ,NaN_i ,NaN_i },
                       { NaN_i ,NaN_i ,NaN_i } };/**<@brief Sublattice*/
  int rboxsub[3][3];/**<@brief Sublattice*/
  /*
   2-body part of the trial wavefunction
  */
  int** Orb;/**<@brief [StdI::nsite][StdI::nsite] Orbital index*/
  int** AntiOrb;/**<@brief [StdI::nsite][StdI::nsite] Anti-periodic switch*/
  int NOrb;/**<@brief Number of independent orbital index*/
  int NSym;/**<@brief Number of translation symmetries,
           Defined from the number of cells in the sub-lattice.*/
#endif
};

#if defined(_HPhi)
/**
@brief Set Largevalue (StdI::LargeValue) for TPQ.
       Sum absolute-value of all one- and two- body terms.
*/
void StdFace::LargeValue() {
  unsigned int ktrans, kintr;
  double LargeValue0;

  LargeValue0 = 0.0;
  for (ktrans = 0; ktrans < StdI::trans.size(); ktrans++) {
    LargeValue0 += std::abs(StdI::trans.at(ktrans));
  }
  for (kintr = 0; kintr < StdI::intr.size(); kintr++) {
    LargeValue0 += std::abs(StdI::intr.at(kintr));
  }
  for (kintr = 0; kintr < StdI::Cintra.size(); kintr++) {
    LargeValue0 += std::abs(StdI::Cintra.at(kintr));
  }
  for (kintr = 0; kintr < StdI::Cinter.size(); kintr++) {
    LargeValue0 += std::abs(StdI::Cinter.at(kintr));
  }
  for (kintr = 0; kintr < StdI::Ex.size(); kintr++) {
    LargeValue0 += 2.0 * std::abs(StdI::Ex.at(kintr));
  }
  for (kintr = 0; kintr < StdI::PairLift.size(); kintr++) {
    LargeValue0 += 2.0 * std::abs(StdI::PairLift.at(kintr));
  }
  for (kintr = 0; kintr < StdI::Hund.size(); kintr++) {
    LargeValue0 += 2.0 * std::abs(StdI::Hund.at(kintr));
  }
  LargeValue0 /= (double)StdI::nsite;
  StdFace::PrintVal_d("LargeValue", &StdI::LargeValue, LargeValue0);
}/*static void StdFace::LargeValue*/
/**
@brief Print calcmod.def
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintCalcMod()
{
  FILE *fp;
  int iCalcType, iCalcModel, iRestart, iCalcSpec, 
    iCalcEigenvec, iInitialVecTpye, InputEigenVec, OutputEigenVec;
  /*
  First, check all parameters and exit if invalid parameters
  */
  fprintf(stdout, "\n  @ CalcMod\n\n");
  /*
   Method
  */
  iCalcEigenvec = 0;
  if (strcmp(StdI::method, "****") == 0){
    fprintf(stdout, "ERROR ! Method is NOT specified !\n");
    StdFace::exit(-1);
  }
  else if (strcmp(StdI::method, "lanczos") == 0) iCalcType = 0;
  else if (strcmp(StdI::method, "lanczosenergy") == 0) { 
    iCalcType = 0; 
    iCalcEigenvec = 1;
  }
  else if (strcmp(StdI::method, "tpq") == 0) iCalcType = 1;
  else if (strcmp(StdI::method, "fulldiag") == 0 ) iCalcType = 2;
  else if (strcmp(StdI::method, "cg") == 0) iCalcType = 3;
  else if (strcmp(StdI::method, "timeevolution") == 0) iCalcType = 4;
  else{
    fprintf(stdout, "\n ERROR ! Unsupported Solver : %s\n", StdI::method);
    StdFace::exit(-1);
  }/*if (strcmp(StdI::method, METHODS) != 0*/
  if (iCalcType != 4) StdI::PumpBody = 0;
  /*
   Model
  */
  if (strcmp(StdI::model, "hubbard") == 0) {
    if (StdI::lGC == 0)iCalcModel = 0;
    else iCalcModel = 3;
  }/*if (strcmp(StdI::model, "hubbard") == 0)*/
  else if (strcmp(StdI::model, "spin") == 0) {
    if (StdI::lGC == 0)iCalcModel = 1;
    else iCalcModel = 4;
  }/*if (strcmp(StdI::model, "spin") == 0)*/
  else if (strcmp(StdI::model, "kondo") == 0) {
    if (StdI::lGC == 0)iCalcModel = 2;
    else iCalcModel = 5;
  }/*if (strcmp(StdI::model, "kondo") == 0)*/
  /*
  Restart
  */
  if (strcmp(StdI::Restart, "****") == 0) {
    strcpy(StdI::Restart, "none\0");
    fprintf(stdout, "          Restart = none        ######  DEFAULT VALUE IS USED  ######\n");
    iRestart = 0;
  }/*if (strcmp(StdI::Restart, "****") == 0)*/
  else {
    fprintf(stdout, "          Restart = %s\n", StdI::Restart);
    if (strcmp(StdI::Restart, "none") == 0) iRestart = 0;
    else if (strcmp(StdI::Restart, "restart_out") == 0 ||
             strcmp(StdI::Restart, "save") == 0) iRestart = 1;
    else if (strcmp(StdI::Restart, "restartsave") == 0 ||
             strcmp(StdI::Restart, "restart")     == 0) iRestart = 2;
    else if (strcmp(StdI::Restart, "restart_in") == 0) iRestart = 3;
    else {
      fprintf(stdout, "\n ERROR ! Restart Mode : %s\n", StdI::Restart);
      StdFace::exit(-1);
    }
  }/*if (strcmp(StdI::Restart, "****") != 0)*/
  /*
  InitialVecType
  */
  if (strcmp(StdI::InitialVecType, "****") == 0) {
    strcpy(StdI::InitialVecType, "c\0");
    fprintf(stdout, "   InitialVecType = c           ######  DEFAULT VALUE IS USED  ######\n");
    iInitialVecTpye = 0;
  }/*if (strcmp(StdI::InitialVecType, "****") == 0)*/
  else {
    fprintf(stdout, "   InitialVecType = %s\n", StdI::InitialVecType);
    if (strcmp(StdI::InitialVecType, "c") == 0) iInitialVecTpye = 0;
    else if (strcmp(StdI::InitialVecType, "r") == 0) iInitialVecTpye = 1;
    else {
      fprintf(stdout, "\n ERROR ! Restart Mode : %s\n", StdI::Restart);
      StdFace::exit(-1);
    }
  }/*if (strcmp(StdI::InitialVecType, "****") != 0)*/
  /*
  EigenVecIO
  */
  InputEigenVec = 0;
  OutputEigenVec = 0;
  if (strcmp(StdI::EigenVecIO, "****") == 0) {
    strcpy(StdI::EigenVecIO, "none\0");
    fprintf(stdout, "       EigenVecIO = none        ######  DEFAULT VALUE IS USED  ######\n");
  }/*if (strcmp(StdI::EigenVecIO, "****") == 0)*/
  else {
    fprintf(stdout, "       EigenVecIO = %s\n", StdI::EigenVecIO);
    if (strcmp(StdI::EigenVecIO, "none") == 0) InputEigenVec = 0;
    else if (strcmp(StdI::EigenVecIO, "in") == 0) InputEigenVec = 1;
    else if (strcmp(StdI::EigenVecIO, "out") == 0) OutputEigenVec = 1;
    else if (strcmp(StdI::EigenVecIO, "inout") == 0) {
      InputEigenVec = 1;
      OutputEigenVec = 1;
    }/*if (strcmp(StdI::EigenVecIO, "inout") == 0)*/
    else {
      fprintf(stdout, "\n ERROR ! EigenVecIO Mode : %s\n", StdI::Restart);
      StdFace::exit(-1);
    }
  }/*if (strcmp(StdI::EigenVecIO, "****") != 0)*/
  if (strcmp(StdI::method, "timeevolution") == 0) InputEigenVec = 1;
  /*
  CalcSpec
  */
  if (strcmp(StdI::CalcSpec, "****") == 0) {
    strcpy(StdI::CalcSpec, "none\0");
    fprintf(stdout, "         CalcSpec = none        ######  DEFAULT VALUE IS USED  ######\n");
    iCalcSpec = 0;
  }/*if (strcmp(StdI::CalcSpec, "****") == 0)*/
  else {
    fprintf(stdout, "         CalcSpec = %s\n", StdI::CalcSpec);
    if (strcmp(StdI::CalcSpec, "none") == 0) iCalcSpec = 0;
    else if (strcmp(StdI::CalcSpec, "normal") == 0) iCalcSpec = 1;
    else if (strcmp(StdI::CalcSpec, "noiteration") == 0) iCalcSpec = 2;
    else if (strcmp(StdI::CalcSpec, "restart_out") == 0) iCalcSpec = 3;
    else if (strcmp(StdI::CalcSpec, "restart_in") == 0) iCalcSpec = 4;
    else if (strcmp(StdI::CalcSpec, "restartsave") == 0 ||
             strcmp(StdI::CalcSpec, "restart")     == 0) iCalcSpec = 5;
    else if (strcmp(StdI::CalcSpec, "scratch") == 0) iCalcSpec = 6;
    else {
      fprintf(stdout, "\n ERROR ! CalcSpec : %s\n", StdI::CalcSpec);
      StdFace::exit(-1);
    }
  }/*if (strcmp(StdI::CalcSpec, "****") != 0)*/

  fp = fopen("calcmod.def", "w");
  fprintf(fp, "#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution\n");
  fprintf(fp, "#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC\n");
  fprintf(fp, "#Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart\n");
  fprintf(fp, "#CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save, 6:Scratch\n");
  fprintf(fp, "CalcType %3d\n", iCalcType);
  fprintf(fp, "CalcModel %3d\n", iCalcModel);
  fprintf(fp, "ReStart %3d\n", iRestart);
  fprintf(fp, "CalcSpec %3d\n", iCalcSpec);
  fprintf(fp, "CalcEigenVec %3d\n", iCalcEigenvec);
  fprintf(fp, "InitialVecType %3d\n", iInitialVecTpye);
  fprintf(fp, "InputEigenVec %3d\n", InputEigenVec);
  fprintf(fp, "OutputEigenVec %3d\n", OutputEigenVec);
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "     calcmod.def is written.\n\n");
}/*static void PrintCalcMod*/
/**
@brief Print single.def or pair.def
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintExcitation() {
  FILE *fp;
  int NumOp, **spin, isite, ispin, icell, itau, iEx, lR;
  double *coef, Cphase, S, Sz;
  double *fourier_r, *fourier_i;

  if (strcmp(StdI::model, "spin") == 0 && StdI::S2 > 1) {
    coef = (double *)malloc(sizeof(double) * (StdI::S2 + 1));
    spin = (int **)malloc(sizeof(int*) * (StdI::S2 + 1));
    for (ispin = 0; ispin < StdI::S2 + 1; ispin++) spin[ispin] = (int *)malloc(sizeof(int) * 2);
  }
  else {
    coef = (double *)malloc(sizeof(double) * 2);
    spin = (int **)malloc(sizeof(int*) * 2);
    for (ispin = 0; ispin < 2; ispin++) spin[ispin] = (int *)malloc(sizeof(int) * 2);
  }

  fourier_r = (double *)malloc(sizeof(double) * StdI::nsite);
  fourier_i = (double *)malloc(sizeof(double) * StdI::nsite);
  
  fprintf(stdout, "\n  @ Spectrum\n\n");

  StdFace::PrintVal_d("SpectrumQW", &StdI::SpectrumQ[0], 0.0);
  StdFace::PrintVal_d("SpectrumQL", &StdI::SpectrumQ[1], 0.0);
  StdFace::PrintVal_d("SpectrumQH", &StdI::SpectrumQ[2], 0.0);

  if (strcmp(StdI::SpectrumType, "****") == 0) {
    strcpy(StdI::SpectrumType, "szsz\0");
    fprintf(stdout, "     SpectrumType = szsz        ######  DEFAULT VALUE IS USED  ######\n");
    if (strcmp(StdI::model, "spin") == 0) {
      NumOp = StdI::S2 + 1;
      for (ispin = 0; ispin <= StdI::S2; ispin++) {
        Sz = (double)ispin - (double)StdI::S2 * 0.5;
        coef[ispin] = Sz;
        spin[ispin][0] = ispin;
        spin[ispin][1] = ispin;
      }
    }
    else {
      NumOp = 2;
      coef[0] = 0.5;
      coef[1] = -0.5;
      spin[0][0] = 0;
      spin[0][1] = 0;
      spin[1][0] = 1;
      spin[1][1] = 1;
    }
    StdI::SpectrumBody = 2;
    lR = 0;
  }
  else {
    fprintf(stdout, "     SpectrumType = %s\n", StdI::SpectrumType);
    if (strcmp(StdI::SpectrumType, "szsz") == 0 ||
        strcmp(StdI::SpectrumType, "szsz_r") == 0) {
      if (strcmp(StdI::model, "spin") == 0) {
        NumOp = StdI::S2 + 1;
        for (ispin = 0; ispin <= StdI::S2; ispin++) {
          Sz = (double)ispin - (double)StdI::S2 * 0.5;
          coef[ispin] = Sz;
          spin[ispin][0] = ispin;
          spin[ispin][1] = ispin;
        }
      }
      else {
        NumOp = 2;
        coef[0] = 0.5;
        coef[1] = -0.5;
        spin[0][0] = 0;
        spin[0][1] = 0;
        spin[1][0] = 1;
        spin[1][1] = 1;
      }
      if (strcmp(StdI::SpectrumType, "szsz") == 0) lR = 0;
      else lR = 1;
      StdI::SpectrumBody = 2;
    }
    else if (strcmp(StdI::SpectrumType, "s+s-") == 0 ||
             strcmp(StdI::SpectrumType, "s+s-_r") == 0) {
      if (strcmp(StdI::model, "spin") == 0 && StdI::S2 > 1) {
        NumOp = StdI::S2;
        S = (double)StdI::S2 * 0.5;
        for (ispin = 0; ispin < StdI::S2; ispin++) {
          Sz = (double)ispin - (double)StdI::S2 * 0.5;
          coef[ispin] = sqrt(S*(S + 1.0) - Sz*(Sz + 1.0));
          spin[ispin][0] = ispin;
          spin[ispin][1] = ispin + 1;
        }
      }
      else {
        NumOp = 1;
        coef[0] = 1.0;
        spin[0][0] = 0;
        spin[0][1] = 1;
      }
      if (strcmp(StdI::SpectrumType, "s+s-") == 0) lR = 0;
      else lR = 1;
      StdI::SpectrumBody = 2;
    }
    else if (strcmp(StdI::SpectrumType, "density") == 0 ||
             strcmp(StdI::SpectrumType, "density_r") == 0) {
      NumOp = 2;
      coef[0] = 1.0;
      coef[1] = 1.0;
      spin[0][0] = 0;
      spin[0][1] = 0;
      spin[1][0] = 1;
      spin[1][1] = 1;
      if (strcmp(StdI::SpectrumType, "density") == 0) lR = 0;
      else lR = 1;
      StdI::SpectrumBody = 2;
    }
    else if (strcmp(StdI::SpectrumType, "up") == 0 ||
             strcmp(StdI::SpectrumType, "up_r") == 0) {
      NumOp = 1;
      coef[0] = 1.0;
      spin[0][0] = 0;
      if (strcmp(StdI::SpectrumType, "up") == 0) lR = 0;
      else lR = 1;
      StdI::SpectrumBody = 1;
    }
    else if (strcmp(StdI::SpectrumType, "down") == 0 ||
             strcmp(StdI::SpectrumType, "down_r") == 0) {
      NumOp = 1;
      coef[0] = 1.0;
      spin[0][0] = 1;
      if (strcmp(StdI::SpectrumType, "down") == 0) lR = 0;
      else lR = 1;
      StdI::SpectrumBody = 1;
    }
    else {
      fprintf(stdout, "\n ERROR ! SpectrumType : %s\n", StdI::SpectrumType);
      StdFace::exit(-1);
    }
  }

  isite = 0;
  for (icell = 0; icell < StdI::NCell; icell++) {
    for (itau = 0; itau < StdI::NsiteUC; itau++) {
      Cphase = (StdI::Cell[icell][0] + StdI::tau[itau][0])*StdI::SpectrumQ[0]
             + (StdI::Cell[icell][1] + StdI::tau[itau][1])*StdI::SpectrumQ[1]
             + (StdI::Cell[icell][2] + StdI::tau[itau][2])*StdI::SpectrumQ[2];
      fourier_r[isite] = cos(2.0*StdI::pi*Cphase);
      fourier_i[isite] = sin(2.0*StdI::pi*Cphase);
      isite += 1;
    }
  }
  if (strcmp(StdI::model, "kondo") == 0) {
    for (isite = 0; isite < StdI::nsite / 2; isite++) {
      fourier_r[isite + StdI::nsite / 2] = fourier_r[isite];
      fourier_i[isite + StdI::nsite / 2] = fourier_i[isite];
    }/*for (isite = 0; isite < StdI::nsite; isite++)*/
  }/*if (strcmp(StdI::model, "kondo") == 0)*/

  if (StdI::SpectrumBody == 1) {
    fp = fopen("single.def", "w");
    fprintf(fp, "=============================================\n");
    if (lR == 0) fprintf(fp, "NSingle %d\n", 2);
    else fprintf(fp, "NSingle %d\n", 1+ StdI::nsite);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "============== Single Excitation ============\n");
    fprintf(fp, "=============================================\n");
    if (lR == 0) {
      if (strcmp(StdI::model, "kondo") == 0) {
        for (iEx = 0; iEx < 2; iEx++) {
          fprintf(fp, "%d\n", StdI::nsite / 2 * NumOp);
          for (isite = StdI::nsite / 2; isite < StdI::nsite; isite++) {
            fprintf(fp, "%d %d 0 %25.15f %25.15f\n", isite, spin[0][0],
              fourier_r[isite] * coef[0], fourier_i[isite] * coef[0]);
          }/*for (isite = 0; isite < StdI::nsite; isite++)*/
        }/*for (iEx = 0; iEx < 2; iEx++)*/
      }/*if (strcmp(StdI::model, "kondo") == 0)*/
      else {
        for (iEx = 0; iEx < 2; iEx++) {
          fprintf(fp, "%d\n", StdI::nsite * NumOp);
          for (isite = 0; isite < StdI::nsite; isite++) {
            fprintf(fp, "%d %d 0 %25.15f %25.15f\n", isite, spin[0][0],
              fourier_r[isite] * coef[0], fourier_i[isite] * coef[0]);
          }/*for (isite = 0; isite < StdI::nsite; isite++)*/
        }/*for (iEx = 0; iEx < 2; iEx++)*/
      }
    }/*if (lR == 0)*/
    else {
      if (strcmp(StdI::model, "kondo") == 0) {
        fprintf(fp, "%d\n", NumOp);
        fprintf(fp, "%d %d 0 %25.15f 0.0\n", StdI::nsite / 2, spin[0][0], coef[0]);
        for (isite = StdI::nsite / 2; isite < StdI::nsite; isite++) {
          fprintf(fp, "%d\n", NumOp);
          fprintf(fp, "%d %d 0 %25.15f 0.0\n", isite, spin[0][0], coef[0]);
        }
      }
      else {
        fprintf(fp, "%d\n", NumOp);
        fprintf(fp, "%d %d 0 %25.15f 0.0\n", 0, spin[0][0], coef[0]);
        for (isite = 0; isite < StdI::nsite; isite++) {
          fprintf(fp, "%d\n", NumOp);
          fprintf(fp, "%d %d 0 %25.15f 0.0\n", isite, spin[0][0], coef[0]);
        }
      }
    }/*if (lR != 0)*/
    fprintf(stdout, "      single.def is written.\n\n");
  }/*if (StdI::SpectrumBody == 1)*/
  else {
    fp = fopen("pair.def", "w");
    fprintf(fp, "=============================================\n");
    if (lR == 0) fprintf(fp, "NPair %d\n", 2);
    else fprintf(fp, "NSingle %d\n", 1 + StdI::nsite);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "=============== Pair Excitation =============\n");
    fprintf(fp, "=============================================\n");
    if (lR == 0) {
      for (iEx = 0; iEx < 2; iEx++) {
        fprintf(fp, "%d\n", StdI::nsite * NumOp);
        for (isite = 0; isite < StdI::nsite; isite++) {
          for (ispin = 0; ispin < NumOp; ispin++) {
            fprintf(fp, "%d %d %d %d 1 %25.15f %25.15f\n",
              isite, spin[ispin][0], isite, spin[ispin][1],
              fourier_r[isite] * coef[ispin], fourier_i[isite] * coef[ispin]);
          }
        }
      }/*for (iEx = 0; iEx < 2; iEx++)*/
    }/*if (lR == 0)*/
    else {
      fprintf(fp, "%d\n", NumOp);
      for (ispin = 0; ispin < NumOp; ispin++) {
        fprintf(fp, "%d %d %d %d 1 %25.15f 0.0\n",
          0, spin[ispin][0], 0, spin[ispin][1], coef[ispin]);
      }
      for (isite = 0; isite < StdI::nsite; isite++) {
        fprintf(fp, "%d\n", NumOp);
        for (ispin = 0; ispin < NumOp; ispin++) {
          fprintf(fp, "%d %d %d %d 1 %25.15f 0.0\n",
            isite, spin[ispin][0], isite, spin[ispin][1], coef[ispin]);
        }
      }
    }/*if (lR != 0)*/
    fprintf(stdout, "        pair.def is written.\n\n");
  }/*if (StdI::SpectrumBody == 2)*/
  fflush(fp);
  fclose(fp);

  free(fourier_r);
  free(fourier_i);
  if (strcmp(StdI::model, "spin") == 0) 
    for (ispin = 0; ispin < StdI::S2 + 1; ispin++) free(spin[ispin]);
  else 
    for (ispin = 0; ispin < 2; ispin++) free(spin[ispin]);
  free(spin);
  free(coef);

}/*static void PrintExcitation()*/
/*
@brief Compute vectorpotential
*/
static void VectorPotential() {
  FILE *fp;
  int it, ii;
  double time;
  double **Et;

  fprintf(stdout, "\n  @ Time-evolution\n\n");

  StdFace::PrintVal_d("VecPotW", &StdI::VecPot[0], 0.0);
  StdFace::PrintVal_d("VecPotL", &StdI::VecPot[1], 0.0);
  StdFace::PrintVal_d("VecPotH", &StdI::VecPot[2], 0.0);
  StdFace::PrintVal_i("Lanczos_max", &StdI::Lanczos_max, 1000);
  StdFace::PrintVal_d("dt", &StdI::dt, 0.1);
  StdFace::PrintVal_d("freq", &StdI::freq, 0.1);
  StdFace::PrintVal_d("tshift", &StdI::tshift, 0.0);
  StdFace::PrintVal_d("tdump", &StdI::tdump, 0.1);
  StdFace::PrintVal_d("Uquench", &StdI::Uquench, 0.0);
  StdFace::PrintVal_i("ExpandCoef", &StdI::ExpandCoef, 10);
  StdI::At = (double **)malloc(sizeof(double*) * StdI::Lanczos_max);
  Et = (double **)malloc(sizeof(double*) * StdI::Lanczos_max);
  for (it = 0; it < StdI::Lanczos_max; it++) {
    StdI::At[it] = (double *)malloc(sizeof(double) * 3);
    Et[it] = (double *)malloc(sizeof(double) * 3);
  }

  if (strcmp(StdI::PumpType, "****") == 0) {
    strcpy(StdI::PumpType, "quench\0");
    fprintf(stdout, "     PumpType = quench        ######  DEFAULT VALUE IS USED  ######\n");
    StdI::PumpBody = 2;
  }/*if (strcmp(StdI::PumpType, "****")*/
  else {
    fprintf(stdout, "     PumpType = %s\n", StdI::PumpType);
    if (strcmp(StdI::PumpType, "quench") == 0) {
      StdI::PumpBody = 2;
    }/*if (strcmp(StdI::PumpType, "quench")*/
    else if (strcmp(StdI::PumpType, "pulselaser") == 0) {
      for (it = 0; it < StdI::Lanczos_max; it++) {
        time = StdI::dt*(double)it;
        for (ii = 0; ii < 3; ii++) {
          StdI::At[it][ii] = StdI::VecPot[ii] * cos(StdI::freq*(time - StdI::tshift))
            * exp(-0.5* (time - StdI::tshift)*(time - StdI::tshift) / (StdI::tdump*StdI::tdump));
          Et[it][ii] = -StdI::VecPot[ii]
            * (
            (StdI::tshift - time) / (StdI::tdump*StdI::tdump) * cos(StdI::freq*(time - StdI::tshift))
              - StdI::freq* sin(StdI::freq*(time - StdI::tshift))
              )
            * exp(-0.5* (time - StdI::tshift)*(time - StdI::tshift) / (StdI::tdump*StdI::tdump));
        }
      }/*for (it = 0; it < StdI::Lanczos_max; it++)*/
      StdI::PumpBody = 1;
    }/*if (strcmp(StdI::PumpType, "pulselaser") == 0)*/
    else if (strcmp(StdI::PumpType, "aclaser") == 0) {
      for (it = 0; it < StdI::Lanczos_max; it++) {
        time = StdI::dt*(double)it;
        for (ii = 0; ii < 3; ii++) {
          StdI::At[it][ii] = StdI::VecPot[ii] * sin(StdI::freq*(time - StdI::tshift));
          Et[it][ii] = StdI::VecPot[ii] * cos(StdI::freq*(time - StdI::tshift)) * StdI::freq;
        }
      }/*for (it = 0; it < StdI::Lanczos_max; it++)*/
      StdI::PumpBody = 1;
    }/*if (strcmp(StdI::PumpType, "aclaser") == 0)*/
    else if (strcmp(StdI::PumpType, "dclaser") == 0) {
      for (it = 0; it < StdI::Lanczos_max; it++) {
        time = StdI::dt*(double)it;
        for (ii = 0; ii < 3; ii++) {
          StdI::At[it][ii] = StdI::VecPot[ii] * time;
          Et[it][ii] = -StdI::VecPot[ii];
        }
      }/*for (it = 0; it < StdI::Lanczos_max; it++)*/
      StdI::PumpBody = 1;
    }/* if (strcmp(StdI::PumpType, "dclaser") == 0)*/
    else {
      fprintf(stdout, "\n ERROR ! PumpType : %s\n", StdI::PumpType);
      StdFace::exit(-1);
    }
  }/*if (! strcmp(StdI::PumpType, "****"))*/

  if (StdI::PumpBody == 1) {
    fp = fopen("potential.dat", "w");
    fprintf(fp, "# Time A_W A_L A_H E_W E_L E_H\n");
    for (it = 0; it < StdI::Lanczos_max; it++) {
      time = StdI::dt*(double)it;
      fprintf(fp, "%f %f %f %f %f %f %f\n",
        time, StdI::At[it][0], StdI::At[it][1], StdI::At[it][2], Et[it][0], Et[it][1], Et[it][2]);
    }
    fflush(fp);
    fclose(fp);

    StdI::pumpindx.resize(StdI::Lanczos_max);
    StdI::pump.resize(StdI::Lanczos_max);
  }/*if (StdI::PumpBody == 1)*/

  for (it = 0; it < StdI::Lanczos_max; it++) free(Et[it]);
  free(Et);
}/*static void VectorPotential()*/
/**
@brief Print single.def or pair.def
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintPump() {
  FILE *fp;
  unsigned int it, isite, ipump, jpump, npump0;

  if (StdI::PumpBody == 1) {

    fp = fopen("teone.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "AllTimeStep %d\n", StdI::Lanczos_max);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "=========  OneBody Time Evolution  ==========\n");
    fprintf(fp, "=============================================\n");
    for (it = 0; it < StdI::Lanczos_max; it++) {
      /*
      Sum equivalent pumping
      */
      for (ipump = 0; ipump < StdI::pump.at(it).size(); ipump++) {
        for (jpump = ipump + 1; jpump < StdI::pump.at(it).size(); jpump++) {
          if (StdI::pumpindx.at(it).at(ipump).at(0) == StdI::pumpindx.at(it).at(jpump).at(0)
            && StdI::pumpindx.at(it).at(ipump).at(1) == StdI::pumpindx.at(it).at(jpump).at(1)
            && StdI::pumpindx.at(it).at(ipump).at(2) == StdI::pumpindx.at(it).at(jpump).at(2)
            && StdI::pumpindx.at(it).at(ipump).at(3) == StdI::pumpindx.at(it).at(jpump).at(3)) {
            StdI::pump.at(it).at(ipump) = StdI::pump.at(it).at(ipump) + StdI::pump.at(it).at(jpump);
            StdI::pump.at(it).at(jpump) = 0.0;
          }
        }/*for (ktrans = jtrans + 1; ktrans < StdI::ntrans; ktrans++)*/
      }/*for (jtrans = 0; jtrans < StdI::ntrans; jtrans++)*/
      /*
      Count the number of finite pumping
      */
      npump0 = 0;
      for (ipump = 0; ipump < StdI::pump.at(it).size(); ipump++)
        if (std::abs(StdI::pump.at(it).at(ipump)) > 0.000001) npump0 += 1;

      fprintf(fp, "%f  %d\n", StdI::dt*(double)it, npump0);
      for (ipump = 0; ipump < StdI::pump.at(it).size(); ipump++) {

        if (std::abs(StdI::pump.at(it).at(ipump)) <= 0.000001) continue;

        fprintf(fp, "%5d %5d %5d %5d %25.15f %25.15f\n",
          StdI::pumpindx.at(it).at(ipump).at(0), StdI::pumpindx.at(it).at(ipump).at(1),
          StdI::pumpindx.at(it).at(ipump).at(2), StdI::pumpindx.at(it).at(ipump).at(3),
          real(StdI::pump.at(it).at(ipump)), imag(StdI::pump.at(it).at(ipump)));
      }/*for (itrans = 0; itrans < StdI::ntrans; itrans++)*/
    }/*for (it = 0; it < StdI::Lanczos_max; it++)*/
    fprintf(stdout, "      teone.def is written.\n\n");
  }
  else {
    fp = fopen("tetwo.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "AllTimeStep %d\n", StdI::Lanczos_max);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "========== TwoBody Time Evolution ===========\n");
    fprintf(fp, "=============================================\n");
    for (it = 0; it < StdI::Lanczos_max; it++) {
      fprintf(fp, "%f  %d\n", StdI::dt*(double)it, StdI::nsite);
      for (isite = 0; isite < StdI::nsite; isite++) {
        fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d %25.15f  %25.15f\n",
          isite, 0, isite, 0, isite, 1, isite, 1, StdI::Uquench, 0.0);
      }/*for (isite = 0; isite < StdI::nsite; isite++)*/
    }/*for (it = 0; it < StdI::Lanczos_max; it++)*/
    fprintf(stdout, "        tetwo.def is written.\n\n");
  }
  fflush(fp);
  fclose(fp);
}/*tatic void PrintPump*/
#elif defined(_mVMC)
/**
@brief Output Anti-parallel orbital index
Free StdI::Orb
*/
static void PrintOrb() {
  FILE *fp;
  int isite, jsite, iOrb;

  fp = fopen("orbitalidx.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NOrbitalIdx %10d\n", StdI::NOrb);
  fprintf(fp, "ComplexType %10d\n", StdI::ComplexType);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "=============================================\n");

  for (isite = 0; isite < StdI::nsite; isite++) {
    for (jsite = 0; jsite < StdI::nsite; jsite++) {
      if (StdI::AntiPeriod[0] == 1 || StdI::AntiPeriod[1] == 1 || StdI::AntiPeriod[2] == 1) {
        fprintf(fp, "%5d  %5d  %5d  %5d\n", isite, jsite, StdI::Orb[isite][jsite], StdI::AntiOrb[isite][jsite]);
      }
      else {
        fprintf(fp, "%5d  %5d  %5d\n", isite, jsite, StdI::Orb[isite][jsite]);
      }
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI::nsite; isite++)*/

  for (iOrb = 0; iOrb < StdI::NOrb; iOrb++)
    fprintf(fp, "%5d  %5d\n", iOrb, 1);

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    orbitalidx.def is written.\n");

  for (isite = 0; isite < StdI::nsite; isite++) free(StdI::Orb[isite]);
  free(StdI::Orb);
}/*void PrintOrb*/
/**
@brief Output parallel orbitalIdx
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintOrbPara() {
  FILE *fp;
  int isite, jsite, NOrbGC, iOrbGC, isite1, jsite1, iorb;
  int **OrbGC, **AntiOrbGC;
  /**@brief
  (1) Copy from anti-parallel orbital index
  */
  OrbGC = (int **)malloc(sizeof(int*) * StdI::nsite);
  AntiOrbGC = (int **)malloc(sizeof(int*) * StdI::nsite);
  for (isite = 0; isite < StdI::nsite; isite++) {
    OrbGC[isite] = (int *)malloc(sizeof(int) * StdI::nsite);
    AntiOrbGC[isite] = (int *)malloc(sizeof(int) * StdI::nsite);
    for (jsite = 0; jsite < StdI::nsite; jsite++) {
      OrbGC[isite][jsite] = StdI::Orb[isite][jsite];
      AntiOrbGC[isite][jsite] = StdI::AntiOrb[isite][jsite];
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI::nsite; isite++)*/
  /**@brief
  (2) Symmetrize
  */
  for (iorb = 0; iorb < StdI::NOrb; iorb++) {
    for (isite = 0; isite < StdI::nsite; isite++) {
      for (jsite = 0; jsite < StdI::nsite; jsite++) {
        if (OrbGC[isite][jsite] == iorb) {
          OrbGC[jsite][isite] = OrbGC[isite][jsite];
        }
      }/*for (jsite = 0; jsite < isite; jsite++)*/
    }/*for (isite = 0; isite < StdI::nsite; isite++)*/
  }/*for (iorb = 0; iorb < StdI::NOrb; iorb++)*/
   /**/
  NOrbGC = 0;
  for (isite = 0; isite < StdI::nsite; isite++) {
    for (jsite = 0; jsite < isite; jsite++) {
      if (OrbGC[isite][jsite] >= 0) {
        iOrbGC = OrbGC[isite][jsite];
        NOrbGC -= 1;
        for (isite1 = 0; isite1 < StdI::nsite; isite1++) {
          for (jsite1 = 0; jsite1 < StdI::nsite; jsite1++) {
            if (OrbGC[isite1][jsite1] == iOrbGC)
              OrbGC[isite1][jsite1] = NOrbGC;
          }/*for (jsite1 = 0; jsite1 < StdI::nsite; jsite1++)*/
        }/*for (isite1 = 0; isite1 < StdI::nsite; isite1++)*/
      }/*if (OrbGC[isite][jsite] >= 0)*/
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI::nsite; isite++)*/
   /**/
  NOrbGC = -NOrbGC;
  for (isite = 0; isite < StdI::nsite; isite++) {
    for (jsite = 0; jsite < StdI::nsite; jsite++) {
      OrbGC[isite][jsite] = -1 - OrbGC[isite][jsite];
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI::nsite; isite++)*/

  fp = fopen("orbitalidxpara.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NOrbitalIdx %10d\n", NOrbGC);
  fprintf(fp, "ComplexType %10d\n", StdI::ComplexType);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "=============================================\n");

  for (isite = 0; isite < StdI::nsite; isite++) {
    for (jsite = 0; jsite < StdI::nsite; jsite++) {
      if (isite >= jsite) continue;
      if (StdI::AntiPeriod[0] == 1 || StdI::AntiPeriod[1] == 1 || StdI::AntiPeriod[2] == 1)
        fprintf(fp, "%5d  %5d  %5d  %5d\n", isite, jsite, OrbGC[isite][jsite], AntiOrbGC[isite][jsite]);
      else
        fprintf(fp, "%5d  %5d  %5d\n", isite, jsite, OrbGC[isite][jsite]);
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI::nsite; isite++)*/

  for (iOrbGC = 0; iOrbGC < NOrbGC; iOrbGC++)
    fprintf(fp, "%5d  %5d\n", iOrbGC, 1);

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    orbitalidxpara.def is written.\n");

  for (isite = 0; isite < StdI::nsite; isite++) {
    free(OrbGC[isite]);
    free(AntiOrbGC[isite]);
  }
  free(OrbGC);
  free(AntiOrbGC);
}/*static void PrintOrbPara*/
/**
@brief Output .def file for Gutzwiller
*/
static void PrintGutzwiller()
{
  FILE *fp;
  int iCell, isite, jsite, NGutzwiller, iGutz;
  int *Gutz;

  Gutz = (int *)malloc(sizeof(int) * StdI::nsite);

  if (std::abs(StdI::NMPTrans) == 1 || StdI::NMPTrans == StdI::NaN_i) {
    if (strcmp(StdI::model, "hubbard") == 0) NGutzwiller = 0;
    else NGutzwiller = -1;

    for (isite = 0; isite < StdI::nsite; isite++) Gutz[isite] = StdI::Orb[isite][isite];

    for (isite = 0; isite < StdI::nsite; isite++) {
      /*
      For Local spin
      */
      if (StdI::locspinflag[isite] != 0) {
        Gutz[isite] = -1;
        continue;
      }
      /**/
      if (Gutz[isite] >= 0) {
        iGutz = Gutz[isite];
        NGutzwiller -= 1;
        for (jsite = 0; jsite < StdI::nsite; jsite++) {
          if (Gutz[jsite] == iGutz)
            Gutz[jsite] = NGutzwiller;
        }/*for (jsite = 0; jsite < StdI::nsite; jsite++)*/
      }/*if (Gutz[isite] >= 0)*/
    }/*for (isite = 0; isite < StdI::nsite; isite++)*/
     /**/
    NGutzwiller = -NGutzwiller;
    for (isite = 0; isite < StdI::nsite; isite++) {
      Gutz[isite] = -1 - Gutz[isite];
    }/*for (isite = 0; isite < StdI::nsite; isite++)*/
  }/*if (std::abs(StdI::NMPTrans) == 1)*/
  else {
    if (strcmp(StdI::model, "hubbard") == 0) NGutzwiller = StdI::NsiteUC;
    else if (strcmp(StdI::model, "spin") == 0) NGutzwiller = 1;
    else NGutzwiller = StdI::NsiteUC + 1;

    for (iCell = 0; iCell < StdI::NCell; iCell++) {
      for (isite = 0; isite < StdI::NsiteUC; isite++) {
        if (strcmp(StdI::model, "hubbard") == 0)
          Gutz[isite + StdI::NsiteUC*iCell] = isite;
        else if (strcmp(StdI::model, "spin") == 0)
          Gutz[isite + StdI::NsiteUC*iCell] = 0;
        else {
          Gutz[isite + StdI::NsiteUC*iCell] = 0;
          Gutz[isite + StdI::NsiteUC*(iCell + StdI::NCell)] = isite + 1;
        }
      }/*for (isite = 0; isite < StdI::NsiteUC; isite++)*/
    }/*for (iCell = 0; iCell < StdI::NCell; iCell++)*/
  }/*if (std::abs(StdI::NMPTrans) != 1)*/

  fp = fopen("gutzwilleridx.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NGutzwillerIdx %10d\n", NGutzwiller);
  fprintf(fp, "ComplexType %10d\n", 0);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "=============================================\n");

  for (isite = 0; isite < StdI::nsite; isite++)
    fprintf(fp, "%5d  %5d\n", isite, Gutz[isite]);

  for (iGutz = 0; iGutz < NGutzwiller; iGutz++) {
    if (strcmp(StdI::model, "hubbard") == 0 || iGutz > 0)
      fprintf(fp, "%5d  %5d\n", iGutz, 1);
    else
      fprintf(fp, "%5d  %5d\n", iGutz, 0);
  }/*for (iGutz = 0; iGutz < NGutzwiller; iGutz++)*/
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    gutzwilleridx.def is written.\n");

  free(Gutz);
}/*static void PrintGutzwiller*/
#endif
/*
@brief Make all characters lower
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void Text2Lower(char *value //!<[inout] @brief Keyword or value
){
  char value2;
  int valuelen, ii;

  valuelen = strlen(value);
  for (ii = 0; ii < valuelen; ii++) {
    value2 = tolower(value[ii]);
    value[ii] = value2;
  }
}/*static void Text2Lower*/
/**
@brief Remove : space etc. from keyword and value in an iput file
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void TrimSpaceQuote(char *value //!<[inout] @brief Keyword or value
){
  char value2[256];
  int valuelen, valuelen2, ii;

  valuelen = strlen(value);
  valuelen2 = 0;
  for (ii = 0; ii < valuelen; ii++){
    if (value[ii] != ' ' &&
      value[ii] != ':' &&
      value[ii] != ';' &&
      value[ii] != '\"' &&
      value[ii] != '\b' &&
      value[ii] != '\\' &&
      value[ii] != '\v' &&
      value[ii] != '\n' &&
      value[ii] != '\0'){
      value2[valuelen2] = value[ii];
      valuelen2++;
    }
  }

  strncpy(value, value2, valuelen2);
  value[valuelen2] = '\0';

}/*static void TrimSpaceQuote*/
/**
@brief Store an input value into the valiable (string)
 If duplicated, HPhi will stop.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StoreWithCheckDup_s(
  char *keyword,//!<[in] keyword read from the input file
  char *valuestring,//!<[in] value read from the input file
  char *value//!<[out]
)
{
  if (strcmp(value, "****") != 0){
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace::exit(-1);
  }
  else{
    strcpy(value, valuestring);
  }
}/*static void StoreWithCheckDup_s*/
/**
@brief Store an input value into the valiable (string) 
Force string lower. If duplicated, HPhi will stop.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StoreWithCheckDup_sl(
  char *keyword,//!<[in] keyword read from the input file
  char *valuestring,//!<[in] value read from the input file
  char *value//!<[out]
)
{
  if (strcmp(value, "****") != 0) {
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace::exit(-1);
  }
  else {
    strcpy(value, valuestring);
    Text2Lower(value);
  }
}/*static void StoreWithCheckDup_sl*/
/**
@brief Store an input value into the valiable (integer)
If duplicated, HPhi will stop.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StoreWithCheckDup_i(
  char *keyword,//!<[in] keyword read from the input file
  char *valuestring,//!<[in] value read from the input file
  int *value//!<[out]
)
{
  int NaN_i = 2147483647;

  if (*value != NaN_i){
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace::exit(-1);
  }
  else{
    sscanf(valuestring, "%d", value);
  }
}/*static void StoreWithCheckDup_i*/
/**
@brief Store an input value into the valiable (double)
If duplicated, HPhi will stop.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StoreWithCheckDup_d(
  char *keyword,//!<[in] keyword read from the input file
  char *valuestring,//!<[in] value read from the input file
  double *value//!<[out]
)
{
  if (std::isnan(*value) == 0){
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace::exit(-1);
  }
  else{
    sscanf(valuestring, "%lf", value);
  }
}/*static void StoreWithCheckDup_d*/
/**
@brief Store an input value into the valiable (Double complex)
      If duplicated, HPhi will stop.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StoreWithCheckDup_c(
  char *keyword,//!<[in] keyword read from the input file
  char *valuestring,//!<[in] value read from the input file
  std::complex<double> *value//!<[out]
)
{
  int num;
  char *valuestring_r, *valuestring_i;
  double value_r, value_i;

  if (std::isnan(real(*value)) == 0) {
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace::exit(-1);
  }
  else {

    if (valuestring[0] == ',') {
      valuestring_r = NULL;
      valuestring_i = strtok(valuestring, ",");
    }
    else {
      valuestring_r = strtok(valuestring, ",");
      valuestring_i = strtok(NULL, ",");
    }
    
    if (valuestring_r == NULL) {
      *value = 0.0;
    }
    else {
      num = sscanf(valuestring_r, "%lf", &value_r);
      if (num == 1) *value = value_r;
      else *value = 0.0;
    }

    if (valuestring_i == NULL) {
      *value += std::complex<double>(0.0, 0.0);
    }
    else {
        num = sscanf(valuestring_i, "%lf", &value_i);
      if (num == 1) *value += std::complex<double>(0.0, value_i);
      else *value += std::complex<double>(0.0, 0.0);
    }
  }
}/*static void StoreWithCheckDup_c*/
/**
@brief Print the locspin file
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintLocSpin() {
  FILE *fp;
  int isite, nlocspin;

  nlocspin = 0;
  for (isite = 0; isite < StdI::nsite; isite++)
    if (StdI::locspinflag[isite] != 0) nlocspin = nlocspin + 1;

  fp = fopen("locspn.def", "w");
  fprintf(fp, "================================ \n");
  fprintf(fp, "NlocalSpin %5d  \n", nlocspin);
  fprintf(fp, "================================ \n");
  fprintf(fp, "========i_0LocSpn_1IteElc ====== \n");
  fprintf(fp, "================================ \n");

  for (isite = 0; isite < StdI::nsite; isite++)
    fprintf(fp, "%5d  %5d\n", isite, StdI::locspinflag[isite]);

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    locspn.def is written.\n");
}/*static void PrintLocSpin*/
/**
@brief Print the transfer file
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintTrans(){
  FILE *fp;
  unsigned int jtrans, ktrans, ntrans0;

  for (jtrans = 0; jtrans < StdI::trans.size(); jtrans++){
    for (ktrans = jtrans + 1; ktrans < StdI::trans.size(); ktrans++){
      if (StdI::transindx.at(jtrans).at(0) == StdI::transindx.at(ktrans).at(0)
        && StdI::transindx.at(jtrans).at(1) == StdI::transindx.at(ktrans).at(1)
        && StdI::transindx.at(jtrans).at(2) == StdI::transindx.at(ktrans).at(2)
        && StdI::transindx.at(jtrans).at(3) == StdI::transindx.at(ktrans).at(3)){
        StdI::trans.at(jtrans) = StdI::trans.at(jtrans) + StdI::trans.at(ktrans);
        StdI::trans.at(ktrans) = 0.0;
      }
    }/*for (ktrans = jtrans + 1; ktrans < StdI::ntrans; ktrans++)*/
  }/*for (jtrans = 0; jtrans < StdI::ntrans; jtrans++)*/

  ntrans0 = 0;
  for (ktrans = 0; ktrans < StdI::trans.size(); ktrans++){
    if (std::abs(StdI::trans.at(ktrans)) > 0.000001) ntrans0 += 1;
  }

  fp = fopen("trans.def", "w");
  fprintf(fp, "======================== \n");
  fprintf(fp, "NTransfer %7d  \n", ntrans0);
  fprintf(fp, "======================== \n");
  fprintf(fp, "========i_j_s_tijs====== \n");
  fprintf(fp, "======================== \n");

  ntrans0 = 0;
  for (ktrans = 0; ktrans < StdI::trans.size(); ktrans++) {
    if (std::abs(StdI::trans.at(ktrans)) > 0.000001)
      fprintf(fp, "%5d %5d %5d %5d %25.15f %25.15f\n",
        StdI::transindx.at(ktrans).at(0), StdI::transindx.at(ktrans).at(1),
        StdI::transindx.at(ktrans).at(2), StdI::transindx.at(ktrans).at(3),
        real(StdI::trans.at(ktrans)), imag(StdI::trans.at(ktrans)));
  }

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "      trans.def is written.\n");
}/*static void PrintTrans*/
/**
@brief Print namelist.def  
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintNamelist(){
  FILE *fp;

  fp = fopen("namelist.def", "w");
  fprintf(                         fp, "         ModPara  modpara.def\n");
  fprintf(                         fp, "         LocSpin  locspn.def\n");
  fprintf(                         fp, "           Trans  trans.def\n");
  if (StdI::LCintra == 1) fprintf( fp, "    CoulombIntra  coulombintra.def\n");
  if (StdI::LCinter == 1) fprintf( fp, "    CoulombInter  coulombinter.def\n");
  if (StdI::LHund == 1)fprintf(    fp, "            Hund  hund.def\n");
  if (StdI::LEx == 1)fprintf(      fp, "        Exchange  exchange.def\n");
  if (StdI::LPairLift == 1)fprintf(fp, "        PairLift  pairlift.def\n");
  if (StdI::LPairHopp == 1)fprintf(fp, "         PairHop  pairhopp.def\n");
  if (StdI::Lintr == 1)fprintf(    fp, "        InterAll  interall.def\n");
  if (StdI::ioutputmode != 0) {
    fprintf(                       fp, "        OneBodyG  greenone.def\n");
    fprintf(                       fp, "        TwoBodyG  greentwo.def\n");
  }
#if defined(_HPhi)
  fprintf(                         fp, "         CalcMod  calcmod.def\n");
  if(StdI::SpectrumBody == 1) 
    fprintf(                       fp, "SingleExcitation  single.def\n");
  else fprintf(                    fp, "  PairExcitation  pair.def\n");
  if (strcmp(StdI::method, "timeevolution") == 0) {
    if (StdI::PumpBody == 1)
      fprintf(fp, "       TEOneBody  teone.def\n");
    else if (StdI::PumpBody == 2)
      fprintf(fp, "       TETwoBody  tetwo.def\n");
  }/*if (strcmp(StdI::method, "timeevolution") == 0)*/
  fprintf(                         fp, "     SpectrumVec  %s_eigenvec_0\n",
                                   StdI::CDataFileHead);
  if (StdI::lBoost == 1) fprintf(  fp, "           Boost  boost.def\n");
#elif defined(_mVMC)
  fprintf(                         fp, "      Gutzwiller  gutzwilleridx.def\n");
  fprintf(                         fp, "         Jastrow  jastrowidx.def\n");
  fprintf(                         fp, "         Orbital  orbitalidx.def\n");
  if (StdI::lGC == 1 || (StdI::Sz2 != 0 && StdI::Sz2 != StdI::NaN_i))
    fprintf(fp, " OrbitalParallel  orbitalidxpara.def\n");
  fprintf(                         fp, "        TransSym  qptransidx.def\n");
  if(strcmp(StdI::lattice, "wannier90") == 0)
    fprintf(fp, "        Initial  initial.def\n");
#endif
  
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    namelist.def is written.\n");
}/*static void PrintNamelist*/
/**
@brief Print modpara.def
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintModPara()
{
  FILE *fp;

  fp = fopen("modpara.def", "w");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Model_Parameters   0\n");
  fprintf(fp, "--------------------\n");
#if defined(_HPhi)
  fprintf(fp, "HPhi_Cal_Parameters\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "CDataFileHead  %s\n", StdI::CDataFileHead);
  fprintf(fp, "CParaFileHead  zqp\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Nsite          %-5d\n", StdI::nsite);
  if (StdI::Sz2 != StdI::NaN_i) fprintf(fp, "2Sz            %-5d\n", StdI::Sz2);
  if (StdI::nelec != StdI::NaN_i) fprintf(fp, "Ncond          %-5d\n", StdI::nelec);
  fprintf(fp, "Lanczos_max    %-5d\n", StdI::Lanczos_max);
  fprintf(fp, "initial_iv     %-5d\n", StdI::initial_iv);
  if(StdI::nvec != StdI::NaN_i) fprintf(fp, "nvec           %-5d\n", StdI::nvec);
  fprintf(fp, "exct           %-5d\n", StdI::exct);
  fprintf(fp, "LanczosEps     %-5d\n", StdI::LanczosEps);
  fprintf(fp, "LanczosTarget  %-5d\n", StdI::LanczosTarget);
  fprintf(fp, "LargeValue     %-25.15e\n", StdI::LargeValue);
  fprintf(fp, "NumAve         %-5d\n", StdI::NumAve);
  fprintf(fp, "ExpecInterval  %-5d\n", StdI::ExpecInterval);
  fprintf(fp, "NOmega         %-5d\n", StdI::Nomega);
  fprintf(fp, "OmegaMax       %-25.15e %-25.15e\n", StdI::OmegaMax, StdI::OmegaIm);
  fprintf(fp, "OmegaMin       %-25.15e %-25.15e\n", StdI::OmegaMin, StdI::OmegaIm);
  fprintf(fp, "OmegaOrg       0.0 0.0\n");
  if (strcmp(StdI::method, "timeevolution") == 0)
    fprintf(fp, "ExpandCoef     %-5d\n", StdI::ExpandCoef);
#elif defined(_mVMC)
  fprintf(fp, "VMC_Cal_Parameters\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "CDataFileHead  %s\n", StdI::CDataFileHead);
  fprintf(fp, "CParaFileHead  %s\n", StdI::CParaFileHead);
  fprintf(fp, "--------------------\n");
  fprintf(fp, "NVMCCalMode    %d\n", StdI::NVMCCalMode);
  /*fprintf(fp, "NLanczosMode   %d\n", StdI::NLanczosMode);*/
  fprintf(fp, "--------------------\n");
  fprintf(fp, "NDataIdxStart  %d\n", StdI::NDataIdxStart);
  fprintf(fp, "NDataQtySmp    %d\n", StdI::NDataQtySmp);
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Nsite          %d\n", StdI::nsite);
  fprintf(fp, "Ncond          %-5d\n", StdI::nelec);
  if (StdI::Sz2 != StdI::NaN_i)
    fprintf(fp, "2Sz            %d\n", StdI::Sz2);
  if (StdI::NSPGaussLeg != StdI::NaN_i)
    fprintf(fp, "NSPGaussLeg    %d\n", StdI::NSPGaussLeg);
  if (StdI::NSPStot != StdI::NaN_i)
    fprintf(fp, "NSPStot        %d\n", StdI::NSPStot);
  fprintf(fp, "NMPTrans       %d\n", StdI::NMPTrans);
  fprintf(fp, "NSROptItrStep  %d\n", StdI::NSROptItrStep);
  fprintf(fp, "NSROptItrSmp   %d\n", StdI::NSROptItrSmp);
  fprintf(fp, "DSROptRedCut   %.10f\n", StdI::DSROptRedCut);
  fprintf(fp, "DSROptStaDel   %.10f\n", StdI::DSROptStaDel);
  fprintf(fp, "DSROptStepDt   %.10f\n", StdI::DSROptStepDt);
  fprintf(fp, "NVMCWarmUp     %d\n", StdI::NVMCWarmUp);
  fprintf(fp, "NVMCInterval   %d\n", StdI::NVMCInterval);
  fprintf(fp, "NVMCSample     %d\n", StdI::NVMCSample);
  fprintf(fp, "NExUpdatePath  %d\n", StdI::NExUpdatePath);
  fprintf(fp, "RndSeed        %d\n", StdI::RndSeed);
  fprintf(fp, "NSplitSize     %d\n", StdI::NSplitSize);
  fprintf(fp, "NStore         %d\n", StdI::NStore);
  fprintf(fp, "NSRCG          %d\n", StdI::NSRCG);
#endif

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "     modpara.def is written.\n");
}/*static void PrintModPara*/
/**
@brief Print greenone.def
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void Print1Green()
{
  FILE *fp;
  int ngreen, igreen, store, xkondo;
  int isite, jsite, ispin, jspin, SiMax, SjMax, isite_k;
  int **greenindx;
  /*
   Set Indices of correlation functions
  */
  ngreen = 0;
  if (StdI::ioutputmode != 0) {
    for (store = 0; store < 2; store++) {

      if (store == 1) {
        greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
        for (igreen = 0; igreen < ngreen; igreen++) {
          greenindx[igreen] = (int *)malloc(sizeof(int) * 4);
        }
        ngreen = 0;
      }/*if (store == 1)*/

      if (strcmp(StdI::model, "kondo") == 0) xkondo = 2;
      else xkondo = 1;

      if (StdI::ioutputmode == 1) {
        for (isite = 0; isite < StdI::NsiteUC*xkondo; isite++) {

          if (isite >= StdI::NsiteUC) isite_k = isite + StdI::nsite / 2;
          else isite_k = isite;

          if (StdI::locspinflag[isite_k] == 0) SiMax = 1;
          else SiMax = StdI::locspinflag[isite_k];

          for (ispin = 0; ispin <= SiMax; ispin++) {
            for (jsite = 0; jsite < StdI::nsite; jsite++) {

              if (StdI::locspinflag[jsite] == 0) SjMax = 1;
              else SjMax = StdI::locspinflag[jsite];

              for (jspin = 0; jspin <= SjMax; jspin++) {

                if (isite_k != jsite &&
                  (StdI::locspinflag[isite_k] != 0 && StdI::locspinflag[jsite] != 0)) continue;

                if (ispin == jspin){
                  if (store == 1) {
                    greenindx[ngreen][0] = isite_k;
                    greenindx[ngreen][1] = ispin;
                    greenindx[ngreen][2] = jsite;
                    greenindx[ngreen][3] = jspin;
                  }
                  ngreen++;
                }

              }/*for (jspin = 0; jspin <= SjMax; jspin++)*/
            }/*for (jsite = 0; jsite < StdI::nsite; jsite++)*/
          }/*for (ispin = 0; ispin <= SiMax; ispin++)*/
        }/*for (isite = 0; isite < StdI::nsite; isite++)*/
      }/*if (StdI::ioutputmode == 1)*/
      else {
        for (isite = 0; isite < StdI::nsite; isite++) {

          if (StdI::locspinflag[isite] == 0) SiMax = 1;
          else SiMax = StdI::locspinflag[isite];

          for (ispin = 0; ispin <= SiMax; ispin++) {
            for (jsite = 0; jsite < StdI::nsite; jsite++) {

              if (StdI::locspinflag[jsite] == 0) SjMax = 1;
              else SjMax = StdI::locspinflag[jsite];

              for (jspin = 0; jspin <= SjMax; jspin++) {

                if (isite != jsite &&
                  (StdI::locspinflag[isite] != 0 && StdI::locspinflag[jsite] != 0)) continue;

                if (store == 1) {
                  greenindx[ngreen][0] = isite;
                  greenindx[ngreen][1] = ispin;
                  greenindx[ngreen][2] = jsite;
                  greenindx[ngreen][3] = jspin;
                }
                ngreen++;

              }/*for (jspin = 0; jspin <= SjMax; jspin++)*/
            }/*for (jsite = 0; jsite < StdI::nsite; jsite++)*/
          }/*for (ispin = 0; ispin <= SiMax; ispin++)*/
        }/*for (isite = 0; isite < StdI::nsite; isite++)*/
      }/*if (StdI::ioutputmode == 2)*/
    }/*if (StdI::ioutputmode != 0)*/

    fp = fopen("greenone.def", "w");
    fprintf(fp, "===============================\n");
    fprintf(fp, "NCisAjs %10d\n", ngreen);
    fprintf(fp, "===============================\n");
    fprintf(fp, "======== Green functions ======\n");
    fprintf(fp, "===============================\n");
    for (igreen = 0; igreen < ngreen; igreen++) {
    fprintf(fp,   "%5d %5d %5d %5d\n",
      greenindx[igreen][0], greenindx[igreen][1], greenindx[igreen][2], greenindx[igreen][3]);
    }
    fflush(fp);
    fclose(fp);

    fprintf(stdout, "    greenone.def is written.\n");

    for (igreen = 0; igreen < ngreen; igreen++) {
      free(greenindx[igreen]);
    }
    free(greenindx);

  }/*if (StdI::ioutputmode != 0) */
}/*static void Print1Green*/
/**
@brief Print greentwo.def
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void Print2Green() {
  FILE *fp;
  int ngreen, store, igreen, xkondo;
  int site1, site2, site3, site4, site1k;
  int spin1, spin2, spin3, spin4;
  int S1Max, S2Max, S3Max, S4Max;
  int **greenindx;
  /*
   Set Indices of correlation functions
  */
  ngreen = 0;
  if (StdI::ioutputmode == 1) {
    for (store = 0; store < 2; store++) {

      if (store == 1) {
        greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
        for (igreen = 0; igreen < ngreen; igreen++)
          greenindx[igreen] = (int *)malloc(sizeof(int) * 8);
        ngreen = 0;
      }/*if (store == 1)*/

      if (strcmp(StdI::model, "kondo") == 0) xkondo = 2;
      else xkondo = 1;

      for (site1 = 0; site1 < StdI::NsiteUC*xkondo; site1++) {

        if (site1 >= StdI::NsiteUC) site1k = site1 + StdI::nsite / 2;
        else site1k = site1;

        if (StdI::locspinflag[site1k] == 0) S1Max = 1;
        else S1Max = StdI::locspinflag[site1k];
        for (spin1 = 0; spin1 <= S1Max; spin1++) {
          for (spin2 = 0; spin2 <= S1Max; spin2++) {

            for (site3 = 0; site3 < StdI::nsite; site3++) {

              if (StdI::locspinflag[site3] == 0) S3Max = 1;
              else S3Max = StdI::locspinflag[site3];
              for (spin3 = 0; spin3 <= S3Max; spin3++) {
                for (spin4 = 0; spin4 <= S3Max; spin4++) {

                  if (spin1 - spin2 + spin3 - spin4 == 0) {
                    if (store == 1) {
#if defined(_mVMC)
                      if (spin1 != spin2 || spin3 != spin4)
                      {
                        greenindx[ngreen][0] = site1k;
                        greenindx[ngreen][1] = spin1;
                        greenindx[ngreen][2] = site3;
                        greenindx[ngreen][3] = spin4;
                        greenindx[ngreen][4] = site3;
                        greenindx[ngreen][5] = spin3;
                        greenindx[ngreen][6] = site1k;
                        greenindx[ngreen][7] = spin2;
                      }
                      else
#endif
                      {
                        greenindx[ngreen][0] = site1k;
                        greenindx[ngreen][1] = spin1;
                        greenindx[ngreen][2] = site1k;
                        greenindx[ngreen][3] = spin2;
                        greenindx[ngreen][4] = site3;
                        greenindx[ngreen][5] = spin3;
                        greenindx[ngreen][6] = site3;
                        greenindx[ngreen][7] = spin4;
                      }
                    }/*if (store == 1)*/
                    ngreen++;
                  }/*if (spin1 - spin2 + spin3 - spin4 == 0)*/

                }/*for (spin4 = 0; spin4 <= S3Max; spin4++)*/
              }/*for (spin3 = 0; spin3 <= S3Max; spin3++*/
            }/*for (site3 = 0; site3 < StdI::nsite; site3++)*/
          }/*for (spin2 = 0; spin2 <= S1Max; spin2++)*/
        }/*for (spin1 = 0; spin1 <= S1Max; spin1++)*/
      }/*for (site1 = 0; site1 < StdI::nsite; site1++)*/

    }/*for (store = 0; store < 2; store++)*/
  }/*if (StdI::ioutputmode == 1)*/
  else if (StdI::ioutputmode == 2) {
    for (store = 0; store < 2; store++) {

      if (store == 1) {
        greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
        for (igreen = 0; igreen < ngreen; igreen++)
          greenindx[igreen] = (int *)malloc(sizeof(int) * 8);
        ngreen = 0;
      }/*if (store == 1)*/

      for (site1 = 0; site1 < StdI::nsite; site1++) {

        if (StdI::locspinflag[site1] == 0) S1Max = 1;
        else S1Max = StdI::locspinflag[site1];
        for (spin1 = 0; spin1 <= S1Max; spin1++) {

          for (site2 = 0; site2 < StdI::nsite; site2++) {

            if (StdI::locspinflag[site1] != 0 && StdI::locspinflag[site2] != 0
              && site1 != site2) continue;

            if (StdI::locspinflag[site2] == 0) S2Max = 1;
            else S2Max = StdI::locspinflag[site2];
            for (spin2 = 0; spin2 <= S2Max; spin2++) {

              for (site3 = 0; site3 < StdI::nsite; site3++) {

                if (StdI::locspinflag[site3] == 0) S3Max = 1;
                else S3Max = StdI::locspinflag[site3];
                for (spin3 = 0; spin3 <= S3Max; spin3++) {

                  for (site4 = 0; site4 < StdI::nsite; site4++) {

                    if (StdI::locspinflag[site3] != 0 && StdI::locspinflag[site4] != 0
                      && site3 != site4) continue;

                    if (StdI::locspinflag[site4] == 0) S4Max = 1;
                    else S4Max = StdI::locspinflag[site4];
                    for (spin4 = 0; spin4 <= S4Max; spin4++) {

                      if (store == 1) {
                        greenindx[ngreen][0] = site1;
                        greenindx[ngreen][1] = spin1;
                        greenindx[ngreen][2] = site2;
                        greenindx[ngreen][3] = spin2;
                        greenindx[ngreen][4] = site3;
                        greenindx[ngreen][5] = spin3;
                        greenindx[ngreen][6] = site4;
                        greenindx[ngreen][7] = spin4;
                      }/*if (store == 1)*/
                      ngreen++;

                    }/*for (spin4 = 0; spin4 <= S4Max; spin4++)*/
                  }/*for (site4 = 0; site4 < StdI::nsite; site4++)*/
                }/*for (spin3 = 0; spin3 <= S3Max; spin3++*/
              }/*for (site3 = 0; site3 < StdI::nsite; site3++)*/
            }/*for (spin2 = 0; spin2 <= S2Max; spin2++)*/
          }/*for (site2 = 0; site2 < StdI::nsite; site2++)*/
        }/*for (spin1 = 0; spin1 <= S1Max; spin1++)*/
      }/*for (site1 = 0; site1 < StdI::nsite; site1++)*/

    }/*for (store = 0; store < 2; store++)*/
  }/*if (StdI::ioutputmode == 2)*/
  if (StdI::ioutputmode != 0) {
    fp = fopen("greentwo.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NCisAjsCktAltDC %10d\n", ngreen);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "======== Green functions for Sq AND Nq ======\n");
    fprintf(fp, "=============================================\n");
    for (igreen = 0; igreen < ngreen; igreen++) {
      fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d\n",
        greenindx[igreen][0], greenindx[igreen][1], greenindx[igreen][2], greenindx[igreen][3],
        greenindx[igreen][4], greenindx[igreen][5], greenindx[igreen][6], greenindx[igreen][7]);
    }
    fflush(fp);
    fclose(fp);

    fprintf(stdout, "    greentwo.def is written.\n");

    for (igreen = 0; igreen < ngreen; igreen++) {
      free(greenindx[igreen]);
    }
    free(greenindx);
  }/*if (StdI::ioutputmode != 0)*/
}/*static void Print2Green()*/
/**
@brief Stop HPhi if unsupported model is read 
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void UnsupportedSystem(
  char *model,//!<[in]
  char *lattice//!<[in]
)
{
  fprintf(stdout, "\nSorry, specified combination, \n");
  fprintf(stdout, "    MODEL : %s  \n", model);
  fprintf(stdout, "  LATTICE : %s, \n", lattice);
  fprintf(stdout, "is unsupported in the STANDARD MODE...\n");
  fprintf(stdout, "Please use the EXPART MODE, or write a NEW FUNCTION and post us.\n");
  StdFace::exit(-1);
}/*static void UnsupportedSystem*/
/**
@brief Verify outputmode
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void CheckOutputMode()
{
  /*
  Form for Correlation function
  */
  if (strcmp(StdI::outputmode, "non") == 0
    || strcmp(StdI::outputmode, "none") == 0
    || strcmp(StdI::outputmode, "off") == 0) {
    StdI::ioutputmode = 0;
    fprintf(stdout, "      ioutputmode = %-10d\n", StdI::ioutputmode);
  }
  else if (strcmp(StdI::outputmode, "cor") == 0
    || strcmp(StdI::outputmode, "corr") == 0
    || strcmp(StdI::outputmode, "correlation") == 0) {
    StdI::ioutputmode = 1;
    fprintf(stdout, "      ioutputmode = %-10d\n", StdI::ioutputmode);
  }
  else if (strcmp(StdI::outputmode, "****") == 0) {
    StdI::ioutputmode = 1;
    fprintf(stdout, "      ioutputmode = %-10d  ######  DEFAULT VALUE IS USED  ######\n", StdI::ioutputmode);
  }
  else if (strcmp(StdI::outputmode, "raw") == 0
    || strcmp(StdI::outputmode, "all") == 0
    || strcmp(StdI::outputmode, "full") == 0) {
    StdI::ioutputmode = 2;
    fprintf(stdout, "      ioutputmode = %-10d\n", StdI::ioutputmode);
  }
  else{
    fprintf(stdout, "\n ERROR ! Unsupported OutPutMode : %s\n", StdI::outputmode);
    StdFace::exit(-1);
  }
}/*static void CheckOutputMode*/
/**
@brief Summary numerical parameter check the combination of
 the number of sites, total spin, the number of electrons
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void CheckModPara()
{

  /**/
#if defined(_HPhi)
  StdFace::PrintVal_i("Lanczos_max", &StdI::Lanczos_max, 2000);
  StdFace::PrintVal_i("initial_iv", &StdI::initial_iv, -1);
  /*StdFace::PrintVal_i("nvec", &StdI::nvec, 1);*/
  StdFace::PrintVal_i("exct", &StdI::exct, 1);
  StdFace::PrintVal_i("LanczosEps", &StdI::LanczosEps, 14);
  StdFace::PrintVal_i("LanczosTarget", &StdI::LanczosTarget, 2);
  if(StdI::LanczosTarget < StdI::exct) StdI::LanczosTarget = StdI::exct;
  StdFace::PrintVal_i("NumAve", &StdI::NumAve, 5);
  StdFace::PrintVal_i("ExpecInterval", &StdI::ExpecInterval, 20);
  StdFace::PrintVal_i("NOmega", &StdI::Nomega, 200);
  StdFace::PrintVal_d("OmegaMax", &StdI::OmegaMax, StdI::LargeValue*StdI::nsite);
  StdFace::PrintVal_d("OmegaMin", &StdI::OmegaMin, -StdI::LargeValue*StdI::nsite);
  StdFace::PrintVal_d("OmegaIm", &StdI::OmegaIm, 0.01* (int)StdI::LargeValue);
#elif defined(_mVMC)
  if (strcmp(StdI::CParaFileHead, "****") == 0) {
    strcpy(StdI::CParaFileHead, "zqp\0");
    fprintf(stdout, "    CParaFileHead = %-12s######  DEFAULT VALUE IS USED  ######\n", StdI::CParaFileHead);
  }
  else fprintf(stdout, "    CParaFileHead = %-s\n", StdI::CParaFileHead);
  
  StdFace::PrintVal_i("NVMCCalMode", &StdI::NVMCCalMode, 0);
  StdFace::PrintVal_i("NLanczosMode", &StdI::NLanczosMode, 0);
  StdFace::PrintVal_i("NDataIdxStart", &StdI::NDataIdxStart, 1);

  if (StdI::NVMCCalMode == 0) StdFace::NotUsed_i("NDataQtySmp", StdI::NDataQtySmp);
  /*else*/StdFace::PrintVal_i("NDataQtySmp", &StdI::NDataQtySmp, 1);

  if (StdI::lGC == 0 && (StdI::Sz2 == 0 || StdI::Sz2 == StdI::NaN_i)) {
    StdFace::PrintVal_i("NSPGaussLeg", &StdI::NSPGaussLeg, 8);
    StdFace::PrintVal_i("NSPStot", &StdI::NSPStot, 0);
  }
  else {
    StdFace::NotUsed_i("NSPGaussLeg", StdI::NSPGaussLeg);
    StdFace::NotUsed_i("NSPStot", StdI::NSPStot);
  }
 
  if (StdI::AntiPeriod[0] == 1 || StdI::AntiPeriod[1] == 1 || StdI::AntiPeriod[2] == 2)
    StdFace::PrintVal_i("NMPTrans", &StdI::NMPTrans, -1);
  else StdFace::PrintVal_i("NMPTrans", &StdI::NMPTrans, 1);

  StdFace::PrintVal_i("NSROptItrStep", &StdI::NSROptItrStep, 1000);
  
  if (StdI::NVMCCalMode == 1) StdFace::NotUsed_i("NSROptItrSmp", StdI::NSROptItrSmp);
  /*else*/ StdFace::PrintVal_i("NSROptItrSmp", &StdI::NSROptItrSmp, StdI::NSROptItrStep/10);

  StdFace::PrintVal_i("NVMCWarmUp", &StdI::NVMCWarmUp, 10);
  StdFace::PrintVal_i("NVMCInterval", &StdI::NVMCInterval, 1);
  StdFace::PrintVal_i("NVMCSample", &StdI::NVMCSample, 1000);

  if (strcmp(StdI::model, "hubbard") == 0) StdI::NExUpdatePath = 0;
  else if (strcmp(StdI::model, "spin") == 0) StdI::NExUpdatePath = 2;
  else if (strcmp(StdI::model, "kondo") == 0) { 
    if(StdI::lGC==0) StdI::NExUpdatePath = 1; 
    else StdI::NExUpdatePath = 3;
  }
  fprintf(stdout, "  %15s = %-10d\n", "NExUpdatePath", StdI::NExUpdatePath);

  StdFace::PrintVal_i("RndSeed", &StdI::RndSeed, 123456789);
  StdFace::PrintVal_i("NSplitSize", &StdI::NSplitSize, 1);
  StdFace::PrintVal_i("NStore", &StdI::NStore, 1);
  StdFace::PrintVal_i("NSRCG", &StdI::NSRCG, 0);

  StdFace::PrintVal_d("DSROptRedCut", &StdI::DSROptRedCut, 0.001);
  StdFace::PrintVal_d("DSROptStaDel", &StdI::DSROptStaDel, 0.02);
  StdFace::PrintVal_d("DSROptStepDt", &StdI::DSROptStepDt, 0.02);
#endif
  /*
   (Un)Conserved variables (Number of electrons, total Sz)
  */
  if (strcmp(StdI::model, "hubbard") == 0){
#if defined(_HPhi)
    if (StdI::lGC == 0) StdFace::RequiredVal_i("nelec", StdI::nelec);
    else {
      StdFace::NotUsed_i("nelec", StdI::nelec);
      StdFace::NotUsed_i("2Sz", StdI::Sz2);
    }
#else
    StdFace::RequiredVal_i("nelec", StdI::nelec);
    if (StdI::lGC == 0) StdFace::PrintVal_i("2Sz", &StdI::Sz2, 0);
    else StdFace::NotUsed_i("2Sz", StdI::Sz2);
#endif
  }
  else if (strcmp(StdI::model, "spin") == 0) {
    StdFace::NotUsed_i("nelec", StdI::nelec);
#if defined(_mVMC)
    StdI::nelec = 0;
#endif
    if (StdI::lGC == 0) StdFace::RequiredVal_i("2Sz", StdI::Sz2);
    else StdFace::NotUsed_i("2Sz", StdI::Sz2);
  }/*else if (strcmp(StdI::model, "spin") == 0)*/
  else if (strcmp(StdI::model, "kondo") == 0) {
#if defined(_HPhi)
    if (StdI::lGC == 0) StdFace::RequiredVal_i("nelec", StdI::nelec);
    else {
      StdFace::NotUsed_i("nelec", StdI::nelec);
      StdFace::NotUsed_i("2Sz", StdI::Sz2);
    }
#else
    StdFace::RequiredVal_i("nelec", StdI::nelec);
    if (StdI::lGC == 0) StdFace::PrintVal_i("2Sz", &StdI::Sz2, 0);
    else StdFace::NotUsed_i("2Sz", StdI::Sz2);
#endif
  }/*else if (strcmp(StdI::model, "kondo") == 0)*/
}/*static void CheckModPara*/
/**
@brief Output .def file for Specific interaction
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintInteractions()
{
  FILE *fp;
  unsigned int nintr0, kintr, jintr;
  /*
   Coulomb INTRA
  */
  for (kintr = 0; kintr < StdI::Cintra.size(); kintr++) {
    for (jintr = kintr + 1; jintr < StdI::Cintra.size(); jintr++)
      if(StdI::CintraIndx.at(jintr) == StdI::CintraIndx.at(kintr))
      {
        StdI::Cintra.at(kintr) += StdI::Cintra.at(jintr);
        StdI::Cintra.at(jintr) = 0.0;
      }
  }
  nintr0 = 0;
  for (kintr = 0; kintr < StdI::Cintra.size(); kintr++) {
    if (fabs(StdI::Cintra.at(kintr)) > 0.000001) nintr0 += 1;
  }
  if (nintr0 == 0 || StdI::lBoost == 1) StdI::LCintra = 0;
  else StdI::LCintra = 1;

  if (StdI::LCintra == 1) {
    fp = fopen("coulombintra.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NCoulombIntra %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "================== CoulombIntra ================\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI::Cintra.size(); kintr++) {
      if (fabs(StdI::Cintra.at(kintr)) > 0.000001)
        fprintf(fp, "%5d %25.15f\n",
          StdI::CintraIndx.at(kintr), StdI::Cintra.at(kintr));
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    coulombintra.def is written.\n");
  }/*if (StdI::LCintra == 1)*/
  /*
  Coulomb INTER
  */
  for (kintr = 0; kintr < StdI::Cinter.size(); kintr++) {
    for (jintr = kintr + 1; jintr < StdI::Cinter.size(); jintr++)
      if (
        (    StdI::CinterIndx.at(jintr).at(0) == StdI::CinterIndx.at(kintr).at(0)
          && StdI::CinterIndx.at(jintr).at(1) == StdI::CinterIndx.at(kintr).at(1))
        ||
        (    StdI::CinterIndx.at(jintr).at(0) == StdI::CinterIndx.at(kintr).at(1)
          && StdI::CinterIndx.at(jintr).at(1) == StdI::CinterIndx.at(kintr).at(0))
        )
      {
        StdI::Cinter.at(kintr) += StdI::Cinter.at(jintr);
        StdI::Cinter.at(jintr) = 0.0;
      }
  }/*for (kintr = 0; kintr < StdI::NCinter; kintr++)*/
  nintr0 = 0;
  for (kintr = 0; kintr < StdI::Cinter.size(); kintr++) {
    if (fabs(StdI::Cinter.at(kintr)) > 0.000001) nintr0 += 1;
  }
  if (nintr0 == 0 || StdI::lBoost == 1) StdI::LCinter = 0;
  else StdI::LCinter = 1;

  if (StdI::LCinter == 1) {
    fp = fopen("coulombinter.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NCoulombInter %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "================== CoulombInter ================\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI::Cinter.size(); kintr++) {
      if (fabs(StdI::Cinter.at(kintr)) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI::CinterIndx.at(kintr).at(0), StdI::CinterIndx.at(kintr).at(1), StdI::Cinter.at(kintr));
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    coulombinter.def is written.\n");
  }/*if (StdI::LCinter == 1)*/
  /*
  Hund
  */
  for (kintr = 0; kintr < StdI::Hund.size(); kintr++) {
    for (jintr = kintr + 1; jintr < StdI::Hund.size(); jintr++)
      if (
        (StdI::HundIndx.at(jintr).at(0) == StdI::HundIndx.at(kintr).at(0)
          && StdI::HundIndx.at(jintr).at(1) == StdI::HundIndx.at(kintr).at(1))
        ||
        (StdI::HundIndx.at(jintr).at(0) == StdI::HundIndx.at(kintr).at(1)
          && StdI::HundIndx.at(jintr).at(1) == StdI::HundIndx.at(kintr).at(0))
        )
      {
        StdI::Hund.at(kintr) += StdI::Hund.at(jintr);
        StdI::Hund.at(jintr) = 0.0;
      }
  }/*for (kintr = 0; kintr < StdI::NHund; kintr++)*/
  nintr0 = 0;
  for (kintr = 0; kintr < StdI::Hund.size(); kintr++) {
    if (fabs(StdI::Hund.at(kintr)) > 0.000001) nintr0 += 1;
  }
  if (nintr0 == 0 || StdI::lBoost == 1) StdI::LHund = 0;
  else StdI::LHund = 1;

  if (StdI::LHund == 1) {
    fp = fopen("hund.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NHund %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "=============== Hund coupling ===============\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI::Hund.size(); kintr++) {
      if (fabs(StdI::Hund.at(kintr)) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI::HundIndx.at(kintr).at(0), StdI::HundIndx.at(kintr).at(1), StdI::Hund.at(kintr));
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    hund.def is written.\n");
  }/*if (StdI::LHund == 1)*/
  /*
  Exchange
  */
  for (kintr = 0; kintr < StdI::Ex.size(); kintr++) {
    for (jintr = kintr + 1; jintr < StdI::Ex.size(); jintr++)
      if (
        (StdI::ExIndx.at(jintr).at(0) == StdI::ExIndx.at(kintr).at(0)
          && StdI::ExIndx.at(jintr).at(1) == StdI::ExIndx.at(kintr).at(1))
        ||
        (StdI::ExIndx.at(jintr).at(0) == StdI::ExIndx.at(kintr).at(1)
          && StdI::ExIndx.at(jintr).at(1) == StdI::ExIndx.at(kintr).at(0))
        )
      {
        StdI::Ex.at(kintr) += StdI::Ex.at(jintr);
        StdI::Ex.at(jintr) = 0.0;
      }
  }/*for (kintr = 0; kintr < StdI::NEx; kintr++)*/
  nintr0 = 0;
  for (kintr = 0; kintr < StdI::Ex.size(); kintr++) {
    if (fabs(StdI::Ex.at(kintr)) > 0.000001) nintr0 += 1;
  }
  if (nintr0 == 0 || StdI::lBoost == 1) StdI::LEx = 0;
  else StdI::LEx = 1;

  if (StdI::LEx == 1) {
    fp = fopen("exchange.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NExchange %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "====== ExchangeCoupling coupling ============\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI::Ex.size(); kintr++) {
      if (fabs(StdI::Ex.at(kintr)) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI::ExIndx.at(kintr).at(0), StdI::ExIndx.at(kintr).at(1), StdI::Ex.at(kintr));
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    exchange.def is written.\n");
  }
  /*
    PairLift
  */
  for (kintr = 0; kintr < StdI::PairLift.size(); kintr++) {
    for (jintr = kintr + 1; jintr < StdI::PairLift.size(); jintr++)
      if (
        (StdI::PLIndx.at(jintr).at(0) == StdI::PLIndx.at(kintr).at(0)
          && StdI::PLIndx.at(jintr).at(1) == StdI::PLIndx.at(kintr).at(1))
        ||
        (StdI::PLIndx.at(jintr).at(0) == StdI::PLIndx.at(kintr).at(1)
          && StdI::PLIndx.at(jintr).at(1) == StdI::PLIndx.at(kintr).at(0))
        )
      {
        StdI::PairLift.at(kintr) += StdI::PairLift.at(jintr);
        StdI::PairLift.at(jintr) = 0.0;
      }
  }/*for (kintr = 0; kintr < StdI::NPairLift; kintr++)*/
  nintr0 = 0;
  for (kintr = 0; kintr < StdI::PairLift.size(); kintr++) {
    if (fabs(StdI::PairLift.at(kintr)) > 0.000001) nintr0 += 1;
  }
  if (nintr0 == 0 || StdI::lBoost == 1) StdI::LPairLift = 0;
  else StdI::LPairLift = 1;

  if (StdI::LPairLift == 1) {
    fp = fopen("pairlift.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NPairLift %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "====== Pair-Lift term ============\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI::PairLift.size(); kintr++) {
      if (fabs(StdI::PairLift.at(kintr)) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI::PLIndx.at(kintr).at(0), StdI::PLIndx.at(kintr).at(1), StdI::PairLift.at(kintr));
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    pairlift.def is written.\n");
  }
  /*
  PairHopp
  */
  for (kintr = 0; kintr < StdI::PairHopp.size(); kintr++) {
    for (jintr = kintr + 1; jintr < StdI::PairHopp.size(); jintr++)
      if (
        (StdI::PHIndx.at(jintr).at(0) == StdI::PHIndx.at(kintr).at(0)
          && StdI::PHIndx.at(jintr).at(1) == StdI::PHIndx.at(kintr).at(1))
        ||
        (StdI::PHIndx.at(jintr).at(0) == StdI::PHIndx.at(kintr).at(1)
          && StdI::PHIndx.at(jintr).at(1) == StdI::PHIndx.at(kintr).at(0))
        )
      {
        StdI::PairHopp.at(kintr) += StdI::PairHopp.at(jintr);
        StdI::PairHopp.at(jintr) = 0.0;
      }
  }/*for (kintr = 0; kintr < StdI::NPairHopp; kintr++)*/
  nintr0 = 0;
  for (kintr = 0; kintr < StdI::PairHopp.size(); kintr++) {
    if (fabs(StdI::PairHopp.at(kintr)) > 0.000001) nintr0 += 1;
  }
  if (nintr0 == 0 || StdI::lBoost == 1) StdI::LPairHopp = 0;
  else StdI::LPairHopp = 1;

  if (StdI::LPairHopp == 1) {
    fp = fopen("pairhopp.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NPairHopp %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "====== Pair-Hopping term ============\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI::PairHopp.size(); kintr++) {
      if (fabs(StdI::PairHopp.at(kintr)) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI::PHIndx.at(kintr).at(0), StdI::PHIndx.at(kintr).at(1), StdI::PairHopp.at(kintr));
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    pairhopp.def is written.\n");
  }
  /*
   InterAll
  */
  //
  // Merge equivalent terms
  //
  for (jintr = 0; jintr < StdI::intr.size(); jintr++) {
    for (kintr = jintr + 1; kintr < StdI::intr.size(); kintr++) {
      if (
        (    StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(kintr).at(0)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(kintr).at(1)
          && StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(kintr).at(2)
          && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(kintr).at(3)
          && StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(kintr).at(4)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(kintr).at(5)
          && StdI::intrindx.at(jintr).at(6) == StdI::intrindx.at(kintr).at(6)
          && StdI::intrindx.at(jintr).at(7) == StdI::intrindx.at(kintr).at(7))
        ||
        (    StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(kintr).at(4)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(kintr).at(5)
          && StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(kintr).at(6)
          && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(kintr).at(7)
          && StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(kintr).at(0)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(kintr).at(1)
          && StdI::intrindx.at(jintr).at(6) == StdI::intrindx.at(kintr).at(2)
          && StdI::intrindx.at(jintr).at(7) == StdI::intrindx.at(kintr).at(3)
          && !(StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(jintr).at(6)
            && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(jintr).at(7))
          && !(StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(jintr).at(4)
            && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(jintr).at(5)))
        ) {
        StdI::intr.at(jintr) = StdI::intr.at(jintr) + StdI::intr.at(kintr);
        StdI::intr.at(kintr) = 0.0;
      }
      else if (
        (    StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(kintr).at(4)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(kintr).at(5)
          && StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(kintr).at(2)
          && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(kintr).at(3)
          && StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(kintr).at(0)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(kintr).at(1)
          && StdI::intrindx.at(jintr).at(6) == StdI::intrindx.at(kintr).at(6)
          && StdI::intrindx.at(jintr).at(7) == StdI::intrindx.at(kintr).at(7)
          && !(StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(jintr).at(0)
            && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(jintr).at(1))
          && !(StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(jintr).at(4)
            && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(jintr).at(5))
          )
        ||
        (    StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(kintr).at(0)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(kintr).at(1)
          && StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(kintr).at(6)
          && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(kintr).at(7)
          && StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(kintr).at(4)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(kintr).at(5)
          && StdI::intrindx.at(jintr).at(6) == StdI::intrindx.at(kintr).at(2)
          && StdI::intrindx.at(jintr).at(7) == StdI::intrindx.at(kintr).at(3)
          && !(StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(jintr).at(2)
            && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(jintr).at(3))
          && !(StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(jintr).at(6)
            && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(jintr).at(7))
          )
        ) {
        StdI::intr.at(jintr) = StdI::intr.at(jintr) - StdI::intr.at(kintr);
        StdI::intr.at(kintr) = 0.0;
      }
    }/*for (kintr = jintr + 1; kintr < StdI::nintr; kintr++)*/
  }/*for (jintr = 0; jintr < StdI::nintr; jintr++)*/
  //
  // Force Hermite term as
  // (c1+ c2 c3+ c4)+ = c4+ c3 c2+ c1
  //
  for (jintr = 0; jintr < StdI::intr.size(); jintr++) {
    for (kintr = jintr + 1; kintr < StdI::intr.size(); kintr++) {
      if ( StdI::intrindx.at(jintr).at(6) == StdI::intrindx.at(kintr).at(4)
        && StdI::intrindx.at(jintr).at(7) == StdI::intrindx.at(kintr).at(5)
        && StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(kintr).at(6)
        && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(kintr).at(7)
        && StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(kintr).at(0)
        && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(kintr).at(1)
        && StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(kintr).at(2)
        && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(kintr).at(3)
        && !(StdI::intrindx.at(kintr).at(0) == StdI::intrindx.at(kintr).at(6)
          && StdI::intrindx.at(kintr).at(1) == StdI::intrindx.at(kintr).at(7))
        && !(StdI::intrindx.at(kintr).at(2) == StdI::intrindx.at(kintr).at(4)
          && StdI::intrindx.at(kintr).at(3) == StdI::intrindx.at(kintr).at(5))
        ) {
        StdI::intrindx.at(kintr).at(0) = StdI::intrindx.at(jintr).at(6);
        StdI::intrindx.at(kintr).at(1) = StdI::intrindx.at(jintr).at(7);
        StdI::intrindx.at(kintr).at(2) = StdI::intrindx.at(jintr).at(4);
        StdI::intrindx.at(kintr).at(3) = StdI::intrindx.at(jintr).at(5);
        StdI::intrindx.at(kintr).at(4) = StdI::intrindx.at(jintr).at(2);
        StdI::intrindx.at(kintr).at(5) = StdI::intrindx.at(jintr).at(3);
        StdI::intrindx.at(kintr).at(6) = StdI::intrindx.at(jintr).at(0);
        StdI::intrindx.at(kintr).at(7) = StdI::intrindx.at(jintr).at(1);
      }
      else if (
        (    StdI::intrindx.at(jintr).at(6) == StdI::intrindx.at(kintr).at(4)
          && StdI::intrindx.at(jintr).at(7) == StdI::intrindx.at(kintr).at(5)
          && StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(kintr).at(2)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(kintr).at(3)
          && StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(kintr).at(0)
          && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(kintr).at(1)
          && StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(kintr).at(6)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(kintr).at(7)
          && !(StdI::intrindx.at(kintr).at(2) == StdI::intrindx.at(kintr).at(0)
            && StdI::intrindx.at(kintr).at(3) == StdI::intrindx.at(kintr).at(1))
          && !(StdI::intrindx.at(kintr).at(2) == StdI::intrindx.at(kintr).at(4)
            && StdI::intrindx.at(kintr).at(3) == StdI::intrindx.at(kintr).at(5)))
        ||
        (    StdI::intrindx.at(jintr).at(6) == StdI::intrindx.at(kintr).at(0)
          && StdI::intrindx.at(jintr).at(7) == StdI::intrindx.at(kintr).at(1)
          && StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(kintr).at(6)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(kintr).at(7)
          && StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(kintr).at(4)
          && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(kintr).at(5)
          && StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(kintr).at(2)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(kintr).at(3)
          && !(StdI::intrindx.at(kintr).at(4) == StdI::intrindx.at(kintr).at(2)
            && StdI::intrindx.at(kintr).at(5) == StdI::intrindx.at(kintr).at(3))
          && !(StdI::intrindx.at(kintr).at(4) == StdI::intrindx.at(kintr).at(6)
            && StdI::intrindx.at(kintr).at(5) == StdI::intrindx.at(kintr).at(7)))
        ) {
        StdI::intrindx.at(kintr).at(0) = StdI::intrindx.at(jintr).at(6);
        StdI::intrindx.at(kintr).at(1) = StdI::intrindx.at(jintr).at(7);
        StdI::intrindx.at(kintr).at(2) = StdI::intrindx.at(jintr).at(4);
        StdI::intrindx.at(kintr).at(3) = StdI::intrindx.at(jintr).at(5);
        StdI::intrindx.at(kintr).at(4) = StdI::intrindx.at(jintr).at(2);
        StdI::intrindx.at(kintr).at(5) = StdI::intrindx.at(jintr).at(3);
        StdI::intrindx.at(kintr).at(6) = StdI::intrindx.at(jintr).at(0);
        StdI::intrindx.at(kintr).at(7) = StdI::intrindx.at(jintr).at(1);

        StdI::intr.at(kintr) = -StdI::intr.at(kintr);
      }
    }/*for (kintr = jintr + 1; kintr < StdI::nintr; kintr++)*/
  }/*for (jintr = 0; jintr < StdI::nintr; jintr++)*/

  for (jintr = 0; jintr < StdI::intr.size(); jintr++) {

    if (
      (StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(jintr).at(4)
        && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(jintr).at(5)) ||
        (StdI::intrindx.at(jintr).at(2) == StdI::intrindx.at(jintr).at(6)
          && StdI::intrindx.at(jintr).at(3) == StdI::intrindx.at(jintr).at(7))
      ) {

      if (!(
        (StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(jintr).at(2)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(jintr).at(3))
        ||
        (StdI::intrindx.at(jintr).at(0) == StdI::intrindx.at(jintr).at(6)
          && StdI::intrindx.at(jintr).at(1) == StdI::intrindx.at(jintr).at(7))
        ||
        (StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(jintr).at(2)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(jintr).at(3))
        ||
        (StdI::intrindx.at(jintr).at(4) == StdI::intrindx.at(jintr).at(6)
          && StdI::intrindx.at(jintr).at(5) == StdI::intrindx.at(jintr).at(7))
        ))
        StdI::intr.at(jintr) = 0.0;
    }
  }/*for (jintr = 0; jintr < StdI::nintr; jintr++)*/

  nintr0 = 0;
  for (kintr = 0; kintr < StdI::intr.size(); kintr++) {
    if (std::abs(StdI::intr.at(kintr)) > 0.000001) nintr0 += 1;
  }
  if (nintr0 == 0 || StdI::lBoost == 1) StdI::Lintr = 0;
  else StdI::Lintr = 1;

  if (StdI::Lintr == 1) {
    fp = fopen("interall.def", "w");
    fprintf(fp, "====================== \n");
    fprintf(fp, "NInterAll %7d  \n", nintr0);
    fprintf(fp, "====================== \n");
    fprintf(fp, "========zInterAll===== \n");
    fprintf(fp, "====================== \n");

    if (StdI::lBoost == 0) {
      nintr0 = 0;
      for (kintr = 0; kintr < StdI::intr.size(); kintr++) {
        if (std::abs(StdI::intr.at(kintr)) > 0.000001)
          fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d %25.15f  %25.15f\n",
            StdI::intrindx.at(kintr).at(0), StdI::intrindx.at(kintr).at(1),
            StdI::intrindx.at(kintr).at(2), StdI::intrindx.at(kintr).at(3),
            StdI::intrindx.at(kintr).at(4), StdI::intrindx.at(kintr).at(5),
            StdI::intrindx.at(kintr).at(6), StdI::intrindx.at(kintr).at(7),
            real(StdI::intr.at(kintr)), imag(StdI::intr.at(kintr)));
      }/*for (kintr = 0; kintr < StdI::nintr; kintr++)*/
    }/* if (StdI::lBoost == 0)*/

    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    interall.def is written.\n");
  }
}/*static void PrintInteractions*/
/**
@brief Main routine for the standard mode
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace::main(
  char *fname//!<[in] Input file name for the standard mode
)
{
  FILE *fp;
  char ctmpline[256];
  char *keyword, *value;

  fprintf(stdout, "\n######  Input Parameter of Standard Intarface  ######\n");
  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stdout, "\n  ERROR !  Cannot open input file %s !\n\n", fname);
    StdFace::exit(-1);
  }
  else {
    fprintf(stdout, "\n  Open Standard-Mode Inputfile %s \n\n", fname);
  }

  while (fgets(ctmpline, 256, fp) != NULL) {

    TrimSpaceQuote(ctmpline);
    if (strncmp(ctmpline, "//", 2) == 0) {
      fprintf(stdout, "  Skipping a line.\n");
      continue;
    }
    else if (ctmpline[0] == '\0') {
      fprintf(stdout, "  Skipping a line.\n");
      continue;
    }
    keyword = strtok(ctmpline, "=");
    value = strtok(NULL, "=");
    if (value == NULL) {
      fprintf(stdout, "\n  ERROR !  \"=\" is NOT found !\n\n");
      StdFace::exit(-1);
    }
    Text2Lower(keyword);
    fprintf(stdout, "  KEYWORD : %-20s | VALUE : %s \n", keyword, value);

    if (strcmp(keyword, "a") == 0) StoreWithCheckDup_d(keyword, value, &StdI::a);
    else if (strcmp(keyword, "a0h") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[0][2]);
    else if (strcmp(keyword, "a0l") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[0][1]);
    else if (strcmp(keyword, "a0w") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[0][0]);
    else if (strcmp(keyword, "a1h") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[1][2]);
    else if (strcmp(keyword, "a1l") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[1][1]);
    else if (strcmp(keyword, "a1w") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[1][0]);
    else if (strcmp(keyword, "a2h") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[2][2]);
    else if (strcmp(keyword, "a2l") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[2][1]);
    else if (strcmp(keyword, "a2w") == 0) StoreWithCheckDup_i(keyword, value, &StdI::box[2][0]);
    else if (strcmp(keyword, "cutoff_j") == 0) StoreWithCheckDup_d(keyword, value, &StdI::cutoff_j);
    else if (strcmp(keyword, "cutoff_jh") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_JR[2]);
    else if (strcmp(keyword, "cutoff_jl") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_JR[1]);
    else if (strcmp(keyword, "cutoff_jw") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_JR[0]);
    else if (strcmp(keyword, "cutoff_length_j") == 0) StoreWithCheckDup_d(keyword, value, &StdI::cutoff_length_J);
    else if (strcmp(keyword, "cutoff_length_u") == 0) StoreWithCheckDup_d(keyword, value, &StdI::cutoff_length_U);
    else if (strcmp(keyword, "cutoff_length_t") == 0) StoreWithCheckDup_d(keyword, value, &StdI::cutoff_length_t);
    else if (strcmp(keyword, "cutoff_t") == 0) StoreWithCheckDup_d(keyword, value, &StdI::cutoff_t);
    else if (strcmp(keyword, "cutoff_th") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_tR[2]);
    else if (strcmp(keyword, "cutoff_tl") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_tR[1]);
    else if (strcmp(keyword, "cutoff_tw") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_tR[0]);
    else if (strcmp(keyword, "cutoff_u") == 0) StoreWithCheckDup_d(keyword, value, &StdI::cutoff_u);
    else if (strcmp(keyword, "cutoff_uh") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_UR[2]);
    else if (strcmp(keyword, "cutoff_ul") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_UR[1]);
    else if (strcmp(keyword, "cutoff_uw") == 0) StoreWithCheckDup_i(keyword, value, &StdI::cutoff_UR[0]);
    else if (strcmp(keyword, "d") == 0) StoreWithCheckDup_d(keyword, value, &StdI::D[2][2]);
    else if (strcmp(keyword, "doublecounting") == 0) StoreWithCheckDup_i(keyword, value, &StdI::double_counting);
    else if (strcmp(keyword, "gamma") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Gamma);
    else if (strcmp(keyword, "h") == 0) StoreWithCheckDup_d(keyword, value, &StdI::h);
    else if (strcmp(keyword, "height") == 0) StoreWithCheckDup_i(keyword, value, &StdI::Height);
    else if (strcmp(keyword, "hlength") == 0) StoreWithCheckDup_d(keyword, value, &StdI::length[2]);
    else if (strcmp(keyword, "hx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[2][0]);
    else if (strcmp(keyword, "hy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[2][1]);
    else if (strcmp(keyword, "hz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[2][2]);
    else if (strcmp(keyword, "j") == 0) StoreWithCheckDup_d(keyword, value, &StdI::JAll);
    else if (strcmp(keyword, "jx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[0][0]);
    else if (strcmp(keyword, "jxy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[0][1]);
    else if (strcmp(keyword, "jxz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[0][2]);
    else if (strcmp(keyword, "jy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[1][1]);
    else if (strcmp(keyword, "jyx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[1][0]);
    else if (strcmp(keyword, "jyz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[1][2]);
    else if (strcmp(keyword, "jz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[2][2]);
    else if (strcmp(keyword, "jzx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[2][0]);
    else if (strcmp(keyword, "jzy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J[2][1]);
    else if (strcmp(keyword, "j0") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0All);
    else if (strcmp(keyword, "j0x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[0][0]);
    else if (strcmp(keyword, "j0xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[0][1]);
    else if (strcmp(keyword, "j0xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[0][2]);
    else if (strcmp(keyword, "j0y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[1][1]);
    else if (strcmp(keyword, "j0yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[1][0]);
    else if (strcmp(keyword, "j0yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[1][2]);
    else if (strcmp(keyword, "j0z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[2][2]);
    else if (strcmp(keyword, "j0zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[2][0]);
    else if (strcmp(keyword, "j0zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0[2][1]);
    else if (strcmp(keyword, "j0'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pAll);
    else if (strcmp(keyword, "j0'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[0][0]);
    else if (strcmp(keyword, "j0'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[0][1]);
    else if (strcmp(keyword, "j0'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[0][2]);
    else if (strcmp(keyword, "j0'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[1][1]);
    else if (strcmp(keyword, "j0'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[1][0]);
    else if (strcmp(keyword, "j0'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[1][2]);
    else if (strcmp(keyword, "j0'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[2][2]);
    else if (strcmp(keyword, "j0'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[2][0]);
    else if (strcmp(keyword, "j0'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0p[2][1]);
    else if (strcmp(keyword, "j0''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0ppAll);
    else if (strcmp(keyword, "j0''x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[0][0]);
    else if (strcmp(keyword, "j0''xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[0][1]);
    else if (strcmp(keyword, "j0''xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[0][2]);
    else if (strcmp(keyword, "j0''y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[1][1]);
    else if (strcmp(keyword, "j0''yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[1][0]);
    else if (strcmp(keyword, "j0''yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[1][2]);
    else if (strcmp(keyword, "j0''z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[2][2]);
    else if (strcmp(keyword, "j0''zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[2][0]);
    else if (strcmp(keyword, "j0''zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J0pp[2][1]);
    else if (strcmp(keyword, "j1") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1All);
    else if (strcmp(keyword, "j1x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[0][0]);
    else if (strcmp(keyword, "j1xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[0][1]);
    else if (strcmp(keyword, "j1xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[0][2]);
    else if (strcmp(keyword, "j1y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[1][1]);
    else if (strcmp(keyword, "j1yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[1][0]);
    else if (strcmp(keyword, "j1yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[1][2]);
    else if (strcmp(keyword, "j1z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[2][2]);
    else if (strcmp(keyword, "j1zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[2][0]);
    else if (strcmp(keyword, "j1zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1[2][1]);
    else if (strcmp(keyword, "j1'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pAll);
    else if (strcmp(keyword, "j1'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[0][0]);
    else if (strcmp(keyword, "j1'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[0][1]);
    else if (strcmp(keyword, "j1'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[0][2]);
    else if (strcmp(keyword, "j1'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[1][1]);
    else if (strcmp(keyword, "j1'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[1][0]);
    else if (strcmp(keyword, "j1'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[1][2]);
    else if (strcmp(keyword, "j1'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[2][2]);
    else if (strcmp(keyword, "j1'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[2][0]);
    else if (strcmp(keyword, "j1'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1p[2][1]);
    else if (strcmp(keyword, "j1''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1ppAll);
    else if (strcmp(keyword, "j1''x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[0][0]);
    else if (strcmp(keyword, "j1''xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[0][1]);
    else if (strcmp(keyword, "j1''xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[0][2]);
    else if (strcmp(keyword, "j1''y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[1][1]);
    else if (strcmp(keyword, "j1''yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[1][0]);
    else if (strcmp(keyword, "j1''yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[1][2]);
    else if (strcmp(keyword, "j1''z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[2][2]);
    else if (strcmp(keyword, "j1''zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[2][0]);
    else if (strcmp(keyword, "j1''zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J1pp[2][1]);
    else if (strcmp(keyword, "j2") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2All);
    else if (strcmp(keyword, "j2x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[0][0]);
    else if (strcmp(keyword, "j2xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[0][1]);
    else if (strcmp(keyword, "j2xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[0][2]);
    else if (strcmp(keyword, "j2y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[1][1]);
    else if (strcmp(keyword, "j2yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[1][0]);
    else if (strcmp(keyword, "j2yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[1][2]);
    else if (strcmp(keyword, "j2z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[2][2]);
    else if (strcmp(keyword, "j2zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[2][0]);
    else if (strcmp(keyword, "j2zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2[2][1]);
    else if (strcmp(keyword, "j2'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pAll);
    else if (strcmp(keyword, "j2'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[0][0]);
    else if (strcmp(keyword, "j2'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[0][1]);
    else if (strcmp(keyword, "j2'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[0][2]);
    else if (strcmp(keyword, "j2'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[1][1]);
    else if (strcmp(keyword, "j2'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[1][0]);
    else if (strcmp(keyword, "j2'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[1][2]);
    else if (strcmp(keyword, "j2'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[2][2]);
    else if (strcmp(keyword, "j2'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[2][0]);
    else if (strcmp(keyword, "j2'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2p[2][1]);
    else if (strcmp(keyword, "j2''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2ppAll);
    else if (strcmp(keyword, "j2''x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[0][0]);
    else if (strcmp(keyword, "j2''xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[0][1]);
    else if (strcmp(keyword, "j2''xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[0][2]);
    else if (strcmp(keyword, "j2''y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[1][1]);
    else if (strcmp(keyword, "j2''yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[1][0]);
    else if (strcmp(keyword, "j2''yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[1][2]);
    else if (strcmp(keyword, "j2''z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[2][2]);
    else if (strcmp(keyword, "j2''zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[2][0]);
    else if (strcmp(keyword, "j2''zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::J2pp[2][1]);
    else if (strcmp(keyword, "j'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::JpAll);
    else if (strcmp(keyword, "j'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[0][0]);
    else if (strcmp(keyword, "j'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[0][1]);
    else if (strcmp(keyword, "j'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[0][2]);
    else if (strcmp(keyword, "j'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[1][1]);
    else if (strcmp(keyword, "j'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[1][0]);
    else if (strcmp(keyword, "j'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[1][2]);
    else if (strcmp(keyword, "j'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[2][2]);
    else if (strcmp(keyword, "j'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[2][0]);
    else if (strcmp(keyword, "j'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jp[2][1]);
    else if (strcmp(keyword, "j''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::JppAll);
    else if (strcmp(keyword, "j''x") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[0][0]);
    else if (strcmp(keyword, "j''xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[0][1]);
    else if (strcmp(keyword, "j''xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[0][2]);
    else if (strcmp(keyword, "j''y") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[1][1]);
    else if (strcmp(keyword, "j''yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[1][0]);
    else if (strcmp(keyword, "j''yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[1][2]);
    else if (strcmp(keyword, "j''z") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[2][2]);
    else if (strcmp(keyword, "j''zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[2][0]);
    else if (strcmp(keyword, "j''zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Jpp[2][1]);
    else if (strcmp(keyword, "k") == 0) StoreWithCheckDup_d(keyword, value, &StdI::K);
    else if (strcmp(keyword, "l") == 0) StoreWithCheckDup_i(keyword, value, &StdI::L);
    else if (strcmp(keyword, "lattice") == 0) StoreWithCheckDup_sl(keyword, value, StdI::lattice);
    else if (strcmp(keyword, "llength") == 0) StoreWithCheckDup_d(keyword, value, &StdI::length[1]);
    else if (strcmp(keyword, "lx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[1][0]);
    else if (strcmp(keyword, "ly") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[1][1]);
    else if (strcmp(keyword, "lz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[1][2]);
    else if (strcmp(keyword, "model") == 0) StoreWithCheckDup_sl(keyword, value, StdI::model);
    else if (strcmp(keyword, "mu") == 0) StoreWithCheckDup_d(keyword, value, &StdI::mu);
    else if (strcmp(keyword, "nelec") == 0) StoreWithCheckDup_i(keyword, value, &StdI::nelec);
    else if (strcmp(keyword, "outputmode") == 0) StoreWithCheckDup_sl(keyword, value, StdI::outputmode);
    else if (strcmp(keyword, "phase0") == 0) StoreWithCheckDup_d(keyword, value, &StdI::phase[0]);
    else if (strcmp(keyword, "phase1") == 0) StoreWithCheckDup_d(keyword, value, &StdI::phase[1]);
    else if (strcmp(keyword, "phase2") == 0) StoreWithCheckDup_d(keyword, value, &StdI::phase[2]);
    else if (strcmp(keyword, "t") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t);
    else if (strcmp(keyword, "t0") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t0);
    else if (strcmp(keyword, "t0'") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t0p);
    else if (strcmp(keyword, "t0''") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t0pp);
    else if (strcmp(keyword, "t1") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t1);
    else if (strcmp(keyword, "t1'") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t1p);
    else if (strcmp(keyword, "t1''") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t1pp);
    else if (strcmp(keyword, "t2") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t2);
    else if (strcmp(keyword, "t2'") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t2p);
    else if (strcmp(keyword, "t2''") == 0) StoreWithCheckDup_c(keyword, value, &StdI::t2pp);
    else if (strcmp(keyword, "t'") == 0) StoreWithCheckDup_c(keyword, value, &StdI::tp);
    else if (strcmp(keyword, "t''") == 0) StoreWithCheckDup_c(keyword, value, &StdI::tpp);
    else if (strcmp(keyword, "u") == 0) StoreWithCheckDup_d(keyword, value, &StdI::U);
    else if (strcmp(keyword, "v") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V);
    else if (strcmp(keyword, "v0") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V0);
    else if (strcmp(keyword, "v0'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V0p);
    else if (strcmp(keyword, "v0''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V0pp);
    else if (strcmp(keyword, "v1") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V1);
    else if (strcmp(keyword, "v1'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V1p);
    else if (strcmp(keyword, "v1''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V1pp);
    else if (strcmp(keyword, "v2") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V2);
    else if (strcmp(keyword, "v2'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V2p);
    else if (strcmp(keyword, "v2''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::V2pp);
    else if (strcmp(keyword, "v'") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Vp);
    else if (strcmp(keyword, "v''") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Vpp);
    else if (strcmp(keyword, "w") == 0) StoreWithCheckDup_i(keyword, value, &StdI::W);
    else if (strcmp(keyword, "wlength") == 0) StoreWithCheckDup_d(keyword, value, &StdI::length[0]);
    else if (strcmp(keyword, "wx") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[0][0]);
    else if (strcmp(keyword, "wy") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[0][1]);
    else if (strcmp(keyword, "wz") == 0) StoreWithCheckDup_d(keyword, value, &StdI::direct[0][2]);
    else if (strcmp(keyword, "2sz") == 0) StoreWithCheckDup_i(keyword, value, &StdI::Sz2);

#if defined(_HPhi)
    else if (strcmp(keyword, "calcspec") == 0) StoreWithCheckDup_sl(keyword, value, StdI::CalcSpec);
    else if (strcmp(keyword, "exct") == 0) StoreWithCheckDup_i(keyword, value, &StdI::exct);
    else if (strcmp(keyword, "eigenvecio") == 0) StoreWithCheckDup_sl(keyword, value, StdI::EigenVecIO);
    else if (strcmp(keyword, "expandcoef") == 0) StoreWithCheckDup_i(keyword, value, &StdI::ExpandCoef);
    else if (strcmp(keyword, "expecinterval") == 0) StoreWithCheckDup_i(keyword, value, &StdI::ExpecInterval);
    else if (strcmp(keyword, "cdatafilehead") == 0) StoreWithCheckDup_s(keyword, value, StdI::CDataFileHead);
    else if (strcmp(keyword, "dt") == 0) StoreWithCheckDup_d(keyword, value, &StdI::dt);
    else if (strcmp(keyword, "flgtemp") == 0) StoreWithCheckDup_i(keyword, value, &StdI::FlgTemp);
    else if (strcmp(keyword, "freq") == 0) StoreWithCheckDup_d(keyword, value, &StdI::freq);
    else if (strcmp(keyword, "initialvectype") == 0) StoreWithCheckDup_sl(keyword, value, StdI::InitialVecType);
    else if (strcmp(keyword, "initial_iv") == 0) StoreWithCheckDup_i(keyword, value, &StdI::initial_iv);
    else if (strcmp(keyword, "lanczoseps") == 0) StoreWithCheckDup_i(keyword, value, &StdI::LanczosEps);
    else if (strcmp(keyword, "lanczostarget") == 0) StoreWithCheckDup_i(keyword, value, &StdI::LanczosTarget);
    else if (strcmp(keyword, "lanczos_max") == 0) StoreWithCheckDup_i(keyword, value, &StdI::Lanczos_max);
    else if (strcmp(keyword, "largevalue") == 0) StoreWithCheckDup_d(keyword, value, &StdI::LargeValue);
    else if (strcmp(keyword, "method") == 0) StoreWithCheckDup_sl(keyword, value, StdI::method);
    else if (strcmp(keyword, "nomega") == 0) StoreWithCheckDup_i(keyword, value, &StdI::Nomega);
    else if (strcmp(keyword, "numave") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NumAve);
    else if (strcmp(keyword, "nvec") == 0) StoreWithCheckDup_i(keyword, value, &StdI::nvec);
    else if (strcmp(keyword, "omegamax") == 0) StoreWithCheckDup_d(keyword, value, &StdI::OmegaMax);
    else if (strcmp(keyword, "omegamin") == 0) StoreWithCheckDup_d(keyword, value, &StdI::OmegaMin);
    else if (strcmp(keyword, "omegaim") == 0) StoreWithCheckDup_d(keyword, value, &StdI::OmegaIm);
    else if (strcmp(keyword, "pumptype") == 0) StoreWithCheckDup_sl(keyword, value, StdI::PumpType);
    else if (strcmp(keyword, "restart") == 0) StoreWithCheckDup_sl(keyword, value, StdI::Restart);
    else if (strcmp(keyword, "spectrumqh") == 0) StoreWithCheckDup_d(keyword, value, &StdI::SpectrumQ[2]);
    else if (strcmp(keyword, "spectrumql") == 0) StoreWithCheckDup_d(keyword, value, &StdI::SpectrumQ[1]);
    else if (strcmp(keyword, "spectrumqw") == 0) StoreWithCheckDup_d(keyword, value, &StdI::SpectrumQ[0]);
    else if (strcmp(keyword, "spectrumtype") == 0) StoreWithCheckDup_sl(keyword, value, StdI::SpectrumType);
    else if (strcmp(keyword, "tdump") == 0) StoreWithCheckDup_d(keyword, value, &StdI::tdump);
    else if (strcmp(keyword, "tshift") == 0) StoreWithCheckDup_d(keyword, value, &StdI::tshift);
    else if (strcmp(keyword, "uquench") == 0) StoreWithCheckDup_d(keyword, value, &StdI::Uquench);
    else if (strcmp(keyword, "vecpoth") == 0) StoreWithCheckDup_d(keyword, value, &StdI::VecPot[2]);
    else if (strcmp(keyword, "vecpotl") == 0) StoreWithCheckDup_d(keyword, value, &StdI::VecPot[1]);
    else if (strcmp(keyword, "vecpotw") == 0) StoreWithCheckDup_d(keyword, value, &StdI::VecPot[0]);
    else if (strcmp(keyword, "2s") == 0) StoreWithCheckDup_i(keyword, value, &StdI::S2);
#elif defined(_mVMC)
    else if (strcmp(keyword, "a0hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[0][2]);
    else if (strcmp(keyword, "a0lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[0][1]);
    else if (strcmp(keyword, "a0wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[0][0]);
    else if (strcmp(keyword, "a1hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[1][2]);
    else if (strcmp(keyword, "a1lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[1][1]);
    else if (strcmp(keyword, "a1wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[1][0]);
    else if (strcmp(keyword, "a2hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[2][2]);
    else if (strcmp(keyword, "a2lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[2][1]);
    else if (strcmp(keyword, "a2wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::boxsub[2][0]);
    else if (strcmp(keyword, "complextype") == 0) StoreWithCheckDup_i(keyword, value, &StdI::ComplexType);
    else if (strcmp(keyword, "cparafilehead") == 0) StoreWithCheckDup_s(keyword, value, StdI::CParaFileHead);
    else if (strcmp(keyword, "dsroptredcut") == 0) StoreWithCheckDup_d(keyword, value, &StdI::DSROptRedCut);
    else if (strcmp(keyword, "dsroptstadel") == 0) StoreWithCheckDup_d(keyword, value, &StdI::DSROptStaDel);
    else if (strcmp(keyword, "dsroptstepdt") == 0) StoreWithCheckDup_d(keyword, value, &StdI::DSROptStepDt);
    else if (strcmp(keyword, "hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::Hsub);
    else if (strcmp(keyword, "lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::Lsub);
    else if (strcmp(keyword, "nvmccalmode") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NVMCCalMode);
    else if (strcmp(keyword, "ndataidxstart") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NDataIdxStart);
    else if (strcmp(keyword, "ndataqtysmp") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NDataQtySmp);
    else if (strcmp(keyword, "nlanczosmode") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NLanczosMode);
    else if (strcmp(keyword, "nmptrans") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NMPTrans);
    else if (strcmp(keyword, "nspgaussleg") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NSPGaussLeg);
    else if (strcmp(keyword, "nsplitsize") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NSplitSize);
    else if (strcmp(keyword, "nspstot") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NSPStot);
    else if (strcmp(keyword, "nsroptitrsmp") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NSROptItrSmp);
    else if (strcmp(keyword, "nsroptitrstep") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NSROptItrStep);
    else if (strcmp(keyword, "nstore") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NStore);
    else if (strcmp(keyword, "nsrcg") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NSRCG);
    else if (strcmp(keyword, "nvmcinterval") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NVMCInterval);
    else if (strcmp(keyword, "nvmcsample") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NVMCSample);
    else if (strcmp(keyword, "nvmcwarmup") == 0) StoreWithCheckDup_i(keyword, value, &StdI::NVMCWarmUp);
    else if (strcmp(keyword, "rndseed") == 0) StoreWithCheckDup_i(keyword, value, &StdI::RndSeed);
    else if (strcmp(keyword, "wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI::Wsub);
#endif
    else {
      fprintf(stdout, "ERROR ! Unsupported Keyword in Standard mode!\n");
      StdFace::exit(-1);
    }
  }
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Construct Model  #######\n");
  fprintf(stdout, "\n");
  /*
  Check the model
  */
  if (strcmp(StdI::CDataFileHead, "****") == 0) {
    strcpy(StdI::CDataFileHead, "zvo\0");
    fprintf(stdout, "    CDataFileHead = %-12s######  DEFAULT VALUE IS USED  ######\n", StdI::CDataFileHead);
  }
  else fprintf(stdout, "    CDataFileHead = %-s\n", StdI::CDataFileHead);
  /**/
  StdI::lGC = 0;
  StdI::lBoost = 0;
  if (strcmp(StdI::model, "fermionhubbard") == 0
    || strcmp(StdI::model, "hubbard") == 0)
    strcpy(StdI::model, "hubbard\0");
  else if(strcmp(StdI::model, "fermionhubbardgc") == 0
    || strcmp(StdI::model, "hubbardgc") == 0) {
    strcpy(StdI::model, "hubbard\0");
    StdI::lGC = 1;
  }
  else if (strcmp(StdI::model, "spin") == 0)
    strcpy(StdI::model, "spin\0");
  else if (strcmp(StdI::model, "spingc") == 0) {
    strcpy(StdI::model, "spin\0");
    StdI::lGC = 1;
  }
#if defined(_HPhi)
  else if(strcmp(StdI::model, "spingcboost") == 0 ||
    strcmp(StdI::model, "spingccma") == 0) {
    strcpy(StdI::model, "spin\0");
    StdI::lGC = 1;
    StdI::lBoost = 1;
  }
#endif
  else if (strcmp(StdI::model, "kondolattice") == 0
    || strcmp(StdI::model, "kondo") == 0) {
    strcpy(StdI::model, "kondo\0");
  }
  else if(strcmp(StdI::model, "kondolatticegc") == 0
    || strcmp(StdI::model, "kondogc") == 0) {
    strcpy(StdI::model, "kondo\0");
    StdI::lGC = 1;
  }
  else UnsupportedSystem(StdI::model, StdI::lattice);
#if defined(_HPhi)
  /*
  Check the method
  */
  if (strcmp(StdI::method, "direct") == 0
    || strcmp(StdI::method, "alldiag") == 0)
    strcpy(StdI::method, "fulldiag\0");
  else if (strcmp(StdI::method, "te") == 0
    || strcmp(StdI::method, "time-evolution") == 0) {
    strcpy(StdI::method, "timeevolution\0");
  }
  /*
  Compute vector potential and electrical field
  */
  if (strcmp(StdI::method, "timeevolution") == 0) VectorPotential();
#endif
  /*>>
  Generate Hamiltonian definition files
  */
  if (strcmp(StdI::lattice, "chain") == 0
    || strcmp(StdI::lattice, "chainlattice") == 0) StdFace::Chain();
  else if (strcmp(StdI::lattice, "face-centeredorthorhombic") == 0
    || strcmp(StdI::lattice, "fcorthorhombic") == 0
    || strcmp(StdI::lattice, "fco") == 0) StdFace::FCOrtho();
  else if (strcmp(StdI::lattice, "honeycomb") == 0
    || strcmp(StdI::lattice, "honeycomblattice") == 0) StdFace::Honeycomb();
  else if (strcmp(StdI::lattice, "kagome") == 0
    || strcmp(StdI::lattice, "kagomelattice") == 0) StdFace::Kagome();
  else if (strcmp(StdI::lattice, "ladder") == 0
    || strcmp(StdI::lattice, "ladderlattice") == 0) StdFace::Ladder();
  else if (strcmp(StdI::lattice, "orthorhombic") == 0
    || strcmp(StdI::lattice, "simpleorthorhombic") == 0) StdFace::Orthorhombic();
  else if (strcmp(StdI::lattice, "pyrochlore") == 0) StdFace::Pyrochlore();
  else if (strcmp(StdI::lattice, "tetragonal") == 0
    || strcmp(StdI::lattice, "tetragonallattice") == 0
    || strcmp(StdI::lattice, "square") == 0
    || strcmp(StdI::lattice, "squarelattice") == 0) StdFace::Tetragonal();
  else if (strcmp(StdI::lattice, "triangular") == 0
    || strcmp(StdI::lattice, "triangularlattice") == 0) StdFace::Triangular();
  else if (strcmp(StdI::lattice, "wannier90") == 0) StdFace::Wannier90();
  else UnsupportedSystem(StdI::model, StdI::lattice);//<<
  /**/
#if defined(_HPhi)
  StdFace::LargeValue();
  /*
  Generate Hamiltonian for Boost
  */
  if (StdI::lBoost == 1) {
    if (strcmp(StdI::lattice, "chain") == 0
      || strcmp(StdI::lattice, "chainlattice") == 0) StdFace::Chain_Boost();
    else if (strcmp(StdI::lattice, "honeycomb") == 0
      || strcmp(StdI::lattice, "honeycomblattice") == 0) StdFace::Honeycomb_Boost();
    else if (strcmp(StdI::lattice, "kagome") == 0
      || strcmp(StdI::lattice, "kagomelattice") == 0) StdFace::Kagome_Boost();
    else if (strcmp(StdI::lattice, "ladder") == 0
      || strcmp(StdI::lattice, "ladderlattice") == 0) StdFace::Ladder_Boost();
    else UnsupportedSystem(StdI::model, StdI::lattice);
  }
#endif
  /**/
  fprintf(stdout, "\n");
  fprintf(stdout, "######  Print Expert input files  ######\n");
  fprintf(stdout, "\n");
  PrintLocSpin();
  PrintTrans();
  PrintInteractions();
  CheckModPara();
  PrintModPara(); 
#if defined(_HPhi)
  PrintExcitation();
  if (strcmp(StdI::method, "timeevolution") == 0) PrintPump();
  PrintCalcMod();
#elif defined(_mVMC)

  if(StdI::lGC == 0 && (StdI::Sz2 == 0 || StdI::Sz2 == StdI::NaN_i)) 
    StdFace::PrintVal_i("ComplexType", &StdI::ComplexType, 0);
  else StdFace::PrintVal_i("ComplexType", &StdI::ComplexType, 1);

  StdFace::generate_orb();
  StdFace::Proj();
  PrintJastrow();
  if(StdI::lGC == 1 || (StdI::Sz2 != 0 && StdI::Sz2 != StdI::NaN_i) )
    PrintOrbPara();
  PrintGutzwiller();
  PrintOrb();
#endif
  CheckOutputMode();
  Print1Green();
  Print2Green();
  PrintNamelist();

  fprintf(stdout, "\n######  Input files are generated.  ######\n\n");
}/*void StdFace::main*/
/**
@page page_addstandard Add new lattice model into Standard mode

@section sec_stan_proc Overall procedure

If you want to create a new lattice file, the following procedures are needed.

-# Copy one of lattice files such as Kagome.cpp 
   (Probably the most similar one) and rename it.
-# @ref sec_lattice
-# Add the function in the header file, StdFace_ModelUtil.hpp.
-# Add entry at
   @dontinclude StdFace_main.cpp
   @skip StdFace\_main
   @until StdIntList
   :
   @skip >>
   @until <<
.
<HR> 
@section sec_lattice Modify lattice model file

To create a new lattice file, please modify the following part
(Kagome.cpp as an example):

@dontinclude Kagome.cpp
Define function as
@skip StdFace\_Kagome(
@until {
Lattice parameters are used only in geometry.dat and lattice.gp
@skip StdFace\_PrintVal\_d
@until Ly
These are unit lattice vectors.\n
Just call this function to initialize all lattice related parameters
@skipline StdFace\_InitSite
where "2" indicates 2D.
@skip tau
@until tau\[2\]\[0\]
These are the fractional coordinates of internal sites.
Then set parameters of Hamiltonian
@skip StdFace\_NotUsed\_J
@until @@
to determine the default values of them and unused parameters.
For more details, please see the description of each function.
Then compute the upper limit of the number of Transfer & Interaction and malloc them.
@skip >>
@until <<
Please estimate the number of bonds per site.
@skipline kCell
In this loop, the parameters of Hamiltonian are computed & stored.
The local term is computed as follows:
@skip >>
@until <<
Probably, it is not necessary to modify this part.
The non-local term is as follows:
@skip >>
@until <<
For more details, please see each function.

@page page_addstandardval Add new input variable into Standard mode

We add new input variable in Standard mode through the following procedure:

@section sec_parse_standard Parse the input file

The input file for Standared mode is read in StdFace::main().
In that function, the keyword value pair is found as follows:

@dontinclude StdFace_main.cpp
@skip (fgets(ctmpline
@until fclose

We have to add new variable (new_val in this case) as
@code{C}
else if (strcmp(keyword, "new_val") == 0) StoreWithCheckDup_i(keyword, value, &StdI::new_val);
@endcode
where StoreWithCheckDup_i() is for the integer variable;
for other type, please refer the above link.

@section sec_share_standard If it should be shared

If the inputted variable should be shared among routines in Standard mode,
we have to add it to the list in StdFace_vals.hpp.

Also, the variable should be intialized before it is read.
This initiallization is performed in the function StdFace::ResetVals().
We have to initialize new variable in this function as:
@code{C}
StdI::new_val = NaN_d;
\endcode
for the float,
@code{C}
StdI::new_val = NaN_i;
\endcode
for the integer, and
@code{C}
strcpy(StdI::new_val, "****\0");
\endcode
for the string.
*/
