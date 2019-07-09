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

namespace DC {
  /*!< CalcType */

  enum {
    Lanczos = 0, /*!< CalcType is Exact Diagonalization method.*/
    TPQCalc, /*!< CalcType is TPQ calculation.*/
    FullDiag, /*!< CalcType is Full Diagonalization method.*/
    CG, /*!< CalcType is CG method */
    TimeEvolution, /*!< CalcType is Time Evolution method*/
    NUM_CALCTYPE
  };

/*!< CalcModel */
  enum {
    Hubbard = 0, /*!< CalcModel is Hubbard model.*/
    Spin, /*!< CalcModel is Spin system.*/
    Kondo, /*!< CalcModel is Kondo model.*/
    HubbardGC, /*!< CalcModel is GrandCanonical Hubbard model.*/
    SpinGC, /*!< CalcModel is GrandCanonical Spin system.*/
    KondoGC, /*!< CalcModel is GrandCanonical Kondo model.*/
    HubbardNConserved, /*!< CalcModel is Hubbard model under particle number conserved.
                       This symmetry is automatically introduced by not defining
                       2Sz in a modpara file.*/
    SpinlessFermion, /*!< CalcModel is GrandCanonical Spinless fermion model.*/
    SpinlessFermionGC, /*!< CalcModel is GrandCanonical Spinless fermionGC model.*/
    NUM_CALCMODEL/*!< Number of model types defined by CalcModel in calcmodfile. 
      Note: HubbardNConserved is not explicitly defined in calcmod file and thus not counted. SpinlessFermion and SpinlessFermionGC are not yet supported*/
  };

/*!< OutputMode */
  enum {
    RAWMODE = 0, /*!< calc one body green function and two body green functions.*/
    CORRMODE,/*!< calc one body green function and two body green functions
              and correlatinos for charge and spin.*/
    NUM_OUTPUTHAM, /*!< Number of output Hamiltonian mode */
    NUM_OUTPUTMODE /*!< Number of output mode.*/
  };

/*!< InputMode */
  enum {
    NUM_INPUTHAM = 2 /*!< Number of input Hamiltonian mode */
  };

/*!< CalcEigenVector */
  enum {
    CALCVEC_NOT = -1, /*!< eigenvector is not calculated*/
    CALCVEC_LANCZOSCG, /*!< Lanczos + CG method*/
    CALCVEC_LANCZOS, /*!< Lanczos method*/
    NUM_CALCEIGENVEC /*!< Number of calculating eigenvector mode.*/
  };
  enum {
    NUM_SETINITAILVEC = 2 /*!< Number of setting type of initial vectors.*/
  };
/*!< CalcSpectrum */
  enum {
    CALCSPEC_NOT,
    RECALC_NOT,
    RECALC_FROM_TMComponents,
    RECALC_OUTPUT_TMComponents_VEC,
    RECALC_FROM_TMComponents_VEC,
    RECALC_INOUT_TMComponents_VEC,
    CALCSPEC_SCRATCH
  };

/*!< ReStartVector */
  enum {
    RESTART_NOT = 0,
    RESTART_OUT,
    RESTART_INOUT,
    RESTART_IN,
    NUM_RESTART 
  };
}
