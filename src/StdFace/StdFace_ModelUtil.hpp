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

void StdFace_exit(int errorcode);

void StdFace_intr(std::complex<double> intr0,
  int site1, int spin1, int site2, int spin2,
  int site3, int spin3, int site4, int spin4);

void StdFace_Hopping(std::complex<double> trans0, int isite, int jsite, double *dR);
void StdFace_trans(std::complex<double> trans0,int isite,int ispin,int jsite,int jspin);
void StdFace_HubbardLocal(double mu0, double h0,
  double Gamma0, double U0, int isite);
void StdFace_MagField(int S2, double h, double Gamma, int isite);

void StdFace_Coulomb(double V, int isite, int jsite);
void StdFace_GeneralJ(double J[3][3],
  int Si2, int Sj2, int isite, int jsite);

void StdFace_PrintVal_d(const char* valname, double *val, double val0);
void StdFace_PrintVal_dd(char* valname, double *val, double val0, double val1);
void StdFace_PrintVal_c(const char* valname, std::complex<double> *val, std::complex<double> val0);
void StdFace_PrintVal_i(const char* valname, int *val, int val0);

void StdFace_NotUsed_d(const char* valname, double val);
void StdFace_NotUsed_i(const char* valname, int val);
void StdFace_NotUsed_c(const char* valname, std::complex<double> val);
void StdFace_NotUsed_J(const char* valname, double JAll, double J[3][3]);

void StdFace_RequiredVal_i(const char* valname, int val);
void StdFace_InputSpinNN(double J[3][3], double JAll, double J0[3][3], double J0All, const char *J0name);
void StdFace_InputSpin(double Jp[3][3], double JpAll, const char *Jpname);
void StdFace_InputCoulombV(double V, double *V0, const char *V0name);
void StdFace_InputHopp(std::complex<double> t, std::complex<double> *t0, const char *t0name);

void StdFace_InitSite(FILE *fp, int dim);
void StdFace_SetLabel(FILE *fp,
  int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC,
  int *isite, int *jsite, int connect, std::complex<double> *Cphase, double *dR);
void StdFace_PrintGeometry();
void StdFace_FindSite(int iW, int iL, int iH, int diW, int diL, int diH,
  int isiteUC, int jsiteUC,
  int *isite, int *jsite, std::complex<double> *Cphase, double *dR);
void StdFace_PrintXSF();

void StdFace_Tetragonal();
void StdFace_Chain();
void StdFace_Ladder();
void StdFace_Triangular();
void StdFace_Honeycomb();
void StdFace_Kagome();
void StdFace_Orthorhombic();
void StdFace_FCOrtho();
void StdFace_Pyrochlore();
void StdFace_Wannier90();

#if defined(_HPhi)
void StdFace_Chain_Boost();
void StdFace_Ladder_Boost();
void StdFace_Honeycomb_Boost();
void StdFace_Kagome_Boost();
#elif defined(_mVMC)
void StdFace_generate_orb();
void StdFace_Proj();
void PrintJastrow();
#endif
