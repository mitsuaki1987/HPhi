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
#include "global.hpp"
#include "lapack_diag.hpp"
#include "mltply.hpp"
#include "mltplyCommon.hpp"
#include "CalcTime.hpp"
#include "common/setmemory.hpp"
#include "CalcSpectrum.hpp"
#include <ctime>
#include <complex>
/**
@brief Compute the Green function with the Lehmann representation and FD
@f[
G(z) = \sum_n \frac{|\langle n|c|0\rangle|^2}{z - E_n}
@f]
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void CalcSpectrumByFullDiag(
  int Nomega,//!<[in] Number of frequencies
  int NdcSpectrum,
  std::complex<double> ***dcSpectrum,//!<[out] [Nomega] Spectrum
  std::complex<double> **dcomega,//!<[in] [Nomega] Frequency
  std::complex<double> **v1Org
)
{
  int kdim, idims, jdims, iomega, idcSpectrum;
  double* energy;
  std::complex<double> **vR, **vL, vRv, vLv, **vLvvRv, **HamS, **VecS;
  /**
  <ul>
  <li>Generate fully stored Hamiltonian. Because HamS & VecS are overwritten,
  copy ::HamS into ::vg.</li>
  */
  energy = d_1d_allocate(Check::idim_maxs);
  vR = cd_2d_allocate(Check::idim_maxs, Check::idim_max);
  vL = cd_2d_allocate(Check::idim_maxs, Check::idim_max);
  vLvvRv = cd_2d_allocate(Check::idim_max, Check::idim_maxs);
  HamS = cd_2d_allocate(Check::idim_maxs, Check::idim_maxs);
  VecS = cd_2d_allocate(Check::idim_maxs, Check::idim_maxs);

  StartTimer(6301);
  zclear(Check::idim_maxs*Check::idim_maxs, &HamS[0][0]);
  zclear(Check::idim_maxs*Check::idim_maxs, &VecS[0][0]);
  for (idims = 0; idims < Check::idim_maxs; idims++) VecS[idims][idims] = 1.0;
  mltply::main(Check::idim_maxs, HamS, VecS, Check::idim_maxs,
               List::b1, List::b2_1, List::b2_2, List::Diagonals);
  StopTimer(6301);
  /**
  <li>::HamS becomes eigenvalues in lapack_diag(), and
   ::VecS becomes eigenvectors</li>
  */
  StartTimer(6302);
  lapack_diag(HamS, energy, VecS);
  StopTimer(6302);
  /**
  <li>Compute @f$|\langle n|c|0\rangle|^2@f$ for all @f$n@f$ and store them into ::VecS,
  where @f$c|0\rangle@f$ is ::vg.</li>
  */
  zclear(Check::idim_maxs * Check::idim_max, &vR[0][0]);
  GetExcitedState::main(Check::idim_max, vR, v1Org, 0);
  for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
    StartTimer(6303);
    zclear(Check::idim_maxs * Check::idim_max, &vL[0][0]);
    GetExcitedState::main(Check::idim_max, vL, v1Org, idcSpectrum + 1);
    for (kdim = 0; kdim < Check::idim_max; kdim++) {
      for (idims = 0; idims < Check::idim_maxs; idims++) {
        vRv = 0.0;
        vLv = 0.0;
        for (jdims = 0; jdims < Check::idim_maxs; jdims++) {
          vRv += conj(VecS[jdims][idims]) * vR[jdims][kdim];
          vLv += conj(VecS[jdims][idims]) * vL[jdims][kdim];
        }
        vLvvRv[kdim][idims] = conj(vLv) * vRv;
      }/*for (idim = 0; idim < idim_max_int; idim++)*/
    }/*for (kdim = 0; kdim < Check::idim_max; kdim++)*/
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
      for (kdim = 0; kdim < Check::idim_max; kdim++) {
        dcSpectrum[kdim][iomega][idcSpectrum] = 0.0;
        for (idims = 0; idims < Check::idim_maxs; idims++) {
          dcSpectrum[kdim][iomega][idcSpectrum] += 
            vLvvRv[kdim][idims] / (dcomega[kdim][iomega] - energy[idims]);
        }/*for (idim = 0; idim < idim_max_int; idim++)*/
      }
    }/*for (iomega = 0; iomega < Nomega; iomega++)*/
    StopTimer(6304);
  }/*for (idcSpectrum = 1; idcSpectrum < NdcSpectrum; idcSpectrum++)*/
  free_cd_2d_allocate(vL);
  free_cd_2d_allocate(vR);
  free_cd_2d_allocate(vLvvRv);
  free_cd_2d_allocate(HamS);
  free_cd_2d_allocate(VecS);
  free_d_1d_allocate(energy);
}/*CalcSpectrumByFullDiag*/

