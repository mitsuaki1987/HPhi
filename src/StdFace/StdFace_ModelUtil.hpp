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
#include <complex>
#include <cstdio>

namespace StdFace {
  void exit(int errorcode);

  void intr(std::complex<double> intr0,
    int site1, int spin1, int site2, int spin2,
    int site3, int spin3, int site4, int spin4);

  void Hopping(std::complex<double> trans0, int isite, int jsite, double* dR);
  void trans(std::complex<double> trans0, int isite, int ispin, int jsite, int jspin);
  void HubbardLocal(double mu0, double h0,
    double Gamma0, double U0, int isite);
  void MagField(int S2, double h, double Gamma, int isite);

  void Coulomb(double V, int isite, int jsite);
  void GeneralJ(double J[3][3],
    int Si2, int Sj2, int isite, int jsite);

  void PrintVal_d(const char* valname, double* val, double val0);
  void PrintVal_dd(char* valname, double* val, double val0, double val1);
  void PrintVal_c(const char* valname, std::complex<double>* val, std::complex<double> val0);
  void PrintVal_i(const char* valname, int* val, int val0);

  void NotUsed_d(const char* valname, double val);
  void NotUsed_i(const char* valname, int val);
  void NotUsed_c(const char* valname, std::complex<double> val);
  void NotUsed_J(const char* valname, double JAll, double J[3][3]);

  void RequiredVal_i(const char* valname, int val);
  void InputSpinNN(double J[3][3], double JAll, double J0[3][3], double J0All, const char* J0name);
  void InputSpin(double Jp[3][3], double JpAll, const char* Jpname);
  void InputCoulombV(double V, double* V0, const char* V0name);
  void InputHopp(std::complex<double> t, std::complex<double>* t0, const char* t0name);

  void InitSite(FILE* fp, int dim);
  void SetLabel(FILE* fp,
    int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC,
    int* isite, int* jsite, int connect, std::complex<double>* Cphase, double* dR);
  void PrintGeometry();
  void FindSite(int iW, int iL, int iH, int diW, int diL, int diH,
    int isiteUC, int jsiteUC,
    int* isite, int* jsite, std::complex<double>* Cphase, double* dR);
  void PrintXSF();
  void FoldSite(int iCellV[3], int nBox[3], int iCellV_fold[3]);

  void Tetragonal();
  void Chain();
  void Ladder();
  void Triangular();
  void Honeycomb();
  void Kagome();
  void Orthorhombic();
  void FCOrtho();
  void Pyrochlore();
  void Wannier90();

#if defined(_HPhi)
  void Chain_Boost();
  void Ladder_Boost();
  void Honeycomb_Boost();
  void Kagome_Boost();
#elif defined(_mVMC)
  void generate_orb();
  void Proj();
  void PrintJastrow();
#endif
}
