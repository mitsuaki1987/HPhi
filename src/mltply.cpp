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
/**
 * @file
 *
 * @brief  Multiplying the wavefunction by the Hamiltonian. @f$ H v_1@f$.
 *
 * @version 0.2
 * @details add function to treat the case of generalspin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
 // Define Mode for mltply
 // complex version
#include "bitcalc.hpp"
#include "mltply.hpp"
#include "mltplySpin.hpp"
#include "mltplyHubbard.hpp"
#include "wrapperMPI.hpp"
#include "CalcTime.hpp"
#include "mltplyCommon.hpp"
#include "diagonalcalc.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
/**
 * @brief Parent function of multiplying the wavefunction by the Hamiltonian. @f$ H v_1@f$.\n
 * First, the calculation of diagonal term is done by using the list @f$ \verb|list_diaognal| @f$. \n
 * Next, the calculation of off-diagonal term is done.\n
 * @note If @f$ \verb|mode| @f$ in BindStruct X is @f$ \verb|M_CORR| @f$, the wave function is not updated. The expected values are only calculated.\n
 * Otherwise, the wavefunction @f$ v_0 @f$ is updated as @f$ v_0 += H v_1@f$.
 *
 * @param X [in] Struct for getting the information of the operators.
 * @param tmp_v0 [in, out]
 * @param tmp_v1 [in]
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int mltply( int nstate, std::complex<double> **tmp_v0,std::complex<double> **tmp_v1) {
  int one = 1;
  long int j=0;
  long int irght=0;
  long int ilft=0;
  long int ihfbit=0;
  std::complex<double> dmv;


  long int i_max;

  StartTimer(1);
  i_max = Check::idim_max;

  if (Def::iFlgGeneralSpin == FALSE) {
    if (GetSplitBitByModel(Def::Nsite, Def::iCalcModel, &irght, &ilft, &ihfbit) != 0) {
      return -1;
    }
  }
  else {
    if (Def::iCalcModel == Spin) {
      if (GetSplitBitForGeneralSpin(Def::Nsite, &ihfbit, Def::SiteToBit) != 0) {
        return -1;
      }
    }
  }  
 
  Large::i_max = i_max;
  Large::irght = irght;
  Large::ilft = ilft;
  Large::ihfbit = ihfbit;
  Large::mode = M_MLTPLY;

  StartTimer(100);
#pragma omp parallel for default(none) private(dmv) \
  firstprivate(i_max) shared(tmp_v0, tmp_v1, list_Diagonal,one,nstate)
  for (j = 1; j <= i_max; j++) {
    dmv = list_Diagonal[j];
    zaxpy_(&nstate, &dmv, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }
  StopTimer(100);
  if (Def::iCalcType == TimeEvolution) diagonalcalcForTE(step_i, &tmp_v0[0][0], &tmp_v1[0][0]);
  
  switch (Def::iCalcModel) {
  case HubbardGC:
    mltplyHubbardGC(nstate, tmp_v0, tmp_v1);
    break;
      
  case KondoGC:
  case Hubbard:
  case Kondo:
    mltplyHubbard(nstate, tmp_v0, tmp_v1);
    break;
      
  case Spin:
    mltplySpin(nstate, tmp_v0, tmp_v1);
    break;
      
  case SpinGC:
    mltplySpinGC(nstate, tmp_v0, tmp_v1);
    break;
      
  default:
    return -1;
  }
  
  StopTimer(1);
  return 0;
}
/**
@brief Wrapper of zaxpy.
*/
void zaxpy_long(
  long int n, 
  std::complex<double> a, 
  std::complex<double> *x, 
  std::complex<double> *y
) {
  long int i;

#pragma omp parallel for default(none) private(i) shared(n, a, x, y)
  for (i = 0; i < n; i++) 
    y[i] += a * x[i];
}
/**
@brief clear std::complex<double> array.
*/
void zclear(
  long int n,
  std::complex<double> *x
) {
  long int i;

#pragma omp parallel for default(none) private(i) shared(n, x)
  for (i = 0; i < n; i++) 
    x[i] = 0.0;
}
