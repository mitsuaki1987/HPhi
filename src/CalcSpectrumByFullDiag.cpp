/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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
/**@file
@brief Functions to perform spectrum calculations with the
full-diagonalization method.
*/
#include <complex>
#include "global.hpp"
#include <ctime>
#include "struct.hpp"
#include "lapack_diag.hpp"
#include "mltply.hpp"
#include "mltplyCommon.hpp"
#include "CalcTime.hpp"
#include "common/setmemory.hpp"
#include "CalcSpectrum.hpp"
/**
@brief Compute the Green function with the Lehmann representation and FD
@f[
G(z) = \sum_n \frac{|\langle n|c|0\rangle|^2}{z - E_n}
@f]
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int CalcSpectrumByFullDiag(
  struct EDMainCalStruct *X,//!<[inout]
  int Nomega,//!<[in] Number of frequencies
  int NdcSpectrum,
  std::complex<double> **dcSpectrum,//!<[out] [Nomega] Spectrum
  std::complex<double> *dcomega,//!<[in] [Nomega] Frequency
  std::complex<double> **v1Org
)
{
  int idim, jdim, iomega, idim_max_int, idcSpectrum;
  std::complex<double> **vR, **vL, vRv, vLv, *vLvvRv;
  /**
  <ul>
  <li>Generate fully stored Hamiltonian. Because ::v0 & ::v1 are overwritten,
  copy ::v0 into ::vg.</li>
  */
  idim_max_int = (int)X->Bind.Check.idim_max;
  vR = cd_2d_allocate(idim_max_int, 1);
  vL = cd_2d_allocate(idim_max_int, 1);
  vLvvRv = cd_1d_allocate(idim_max_int);

  StartTimer(6301);
  zclear((X->Bind.Check.idim_max + 1)*(X->Bind.Check.idim_max + 1), &v0[0][0]);
  zclear((X->Bind.Check.idim_max + 1)*(X->Bind.Check.idim_max + 1), &v1[0][0]);
  for (idim = 1; idim <= X->Bind.Check.idim_max; idim++) v1[idim][idim] = 1.0;
  mltply(&(X->Bind), X->Bind.Check.idim_max, v0, v1);
  StopTimer(6301);
  /**
  <li>::v0 becomes eigenvalues in lapack_diag(), and
   ::v1 becomes eigenvectors</li>
  */
  StartTimer(6302);
  lapack_diag(&(X->Bind));
  StopTimer(6302);
  /**
  <li>Compute @f$|\langle n|c|0\rangle|^2@f$ for all @f$n@f$ and store them into ::v1,
  where @f$c|0\rangle@f$ is ::vg.</li>
  */
  zclear(X->Bind.Check.idim_max, &vR[1][0]);
  GetExcitedState(&(X->Bind), 1, vR, v1Org, 0);
  for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
    StartTimer(6303);
    zclear(X->Bind.Check.idim_max, &vL[1][0]);
    GetExcitedState(&(X->Bind), 1, vL, v1Org, idcSpectrum + 1);
    for (idim = 0; idim < idim_max_int; idim++) {
      vRv = 0.0;
      vLv = 0.0;
      for (jdim = 0; jdim < idim_max_int; jdim++) {
        vRv += conj(v1[jdim][idim]) * vR[jdim][1];
        vLv += conj(v1[jdim][idim]) * vL[jdim][1];
      }
      vLvvRv[idim] = conj(vLv) * vRv;
    }/*for (idim = 0; idim < idim_max_int; idim++)*/
    StopTimer(6303);
    /**
    <li>Compute spectrum
    @f[
    \sum_n \frac{|\langle n|c|0\rangle|^2}{z - E_n}
    @f]
    </li>
    </ul>
    */
    StartTimer(6304);
    for (iomega = 0; iomega < Nomega; iomega++) {
      dcSpectrum[iomega][idcSpectrum] = 0.0;
      for (idim = 0; idim < idim_max_int; idim++) {
        dcSpectrum[iomega][idcSpectrum] += vLvvRv[idim] / (dcomega[iomega] - X->Bind.Phys.energy[idim]);
      }/*for (idim = 0; idim < idim_max_int; idim++)*/
    }/*for (iomega = 0; iomega < Nomega; iomega++)*/
    StopTimer(6304);
  }/*for (idcSpectrum = 1; idcSpectrum < NdcSpectrum; idcSpectrum++)*/
  free_cd_2d_allocate(vL);
  free_cd_2d_allocate(vR);
  free_cd_1d_allocate(vLvvRv);
  return TRUE;
}/*CalcSpectrumByFullDiag*/

