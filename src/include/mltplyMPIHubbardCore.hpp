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
#include <complex>

namespace mltply {
  namespace Hubbard {
    void X_CisAisCjtAjt_MPI(int org_isite1, int org_ispin1,
      int org_isite3, int org_ispin3, std::complex<double> tmp_V, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAjtCkuAlv_MPI(int isite1, int isigma1,
      int isite2, int isigma2, int isite3, int isigma3, int isite4, int isigma4,
      std::complex<double> tmp_V, int nstate, 
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAjtCkuAku_MPI(int isite1, int isigma1,
      int isite2, int isigma2, int isite3, int isigma3,
      std::complex<double> tmp_V, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAisCjtAku_MPI(int isite1, int isigma1, int isite3, int isigma3,
      int isite4, int isigma4, std::complex<double> tmp_V, int nstate, 
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAis_MPI(int org_isite1, int org_ispin1,
      std::complex<double> tmp_V, int nstate, 
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_Cis_MPI(int org_isite, int org_ispin,
      std::complex<double> tmp_trans, int nstate, 
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
      long int idim_max, long int* Tpow,
      long int _irght, long int _ilft, long int _ihfbit);
    void X_Ajt_MPI(int org_isite, int org_ispin,
      std::complex<double> tmp_trans, int nstate, 
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
      long int idim_max, long int* Tpow,
      long int _irght, long int _ilft, long int _ihfbit);
  }
  namespace HubbardGC {
    void X_CisAisCjtAjt_MPI(int org_isite1, int org_ispin1,
      int org_isite3, int org_ispin3,
      std::complex<double> tmp_V, int nstate, 
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAjtCkuAlv_MPI(int isite1, int isigma1, int isite2, int isigma2, 
      int isite3, int isigma3, int isite4, int isigma4,
      std::complex<double> tmp_V, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAjtCkuAku_MPI(int isite1, int isigma1,
      int isite2, int isigma2, int isite3, int isigma3,
      std::complex<double> tmp_V, int nstate, 
      std::complex<double>** tmp_v0,
      std::complex<double>** tmp_v1
    );
    void X_CisAisCjtAku_MPI(int isite1, int isigma1,
      int isite3, int isigma3, int isite4, int isigma4,
      std::complex<double> tmp_V, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAis_MPI(int org_isite1, int org_ispin1,
      std::complex<double> tmp_V, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_CisAjt_MPI(int org_isite1, int org_ispin1,
      int org_isite2, int org_ispin2,
      std::complex<double> tmp_V, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1);
    void X_Cis_MPI(int org_isite, int org_ispin,
      std::complex<double> tmp_trans, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
      long int idim_max, long int* Tpow);
    void X_Ajt_MPI(int org_isite, int org_ispin,
      std::complex<double> tmp_trans, int nstate,
      std::complex<double>** tmp_v0, std::complex<double>** tmp_v1,
      long int idim_max, long int* Tpow);
  }
}
int CheckPE(int isite);
int CheckBit_Cis(long int is1_spin, long int orgbit, long int* offbit);
int CheckBit_Ajt(long int is1_spin, long int orgbit, long int* offbit);
int CheckBit_InterAllPE(int isite1, int isigma1, int isite2, int isigma2,
  int isite3, int isigma3, int isite4, int isigma4, long int orgbit, long int* offbit);

int CheckBit_PairPE(int isite1, int isigma1, int isite3, int isigma3, long int orgbit);
int GetSgnInterAll(long int isite1, long int isite2, long int isite3, long int isite4,
  int* Fsgn, long int orgbit, long int* offbit);
