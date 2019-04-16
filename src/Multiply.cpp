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

#include "Common.hpp"
#include "diagonalcalc.hpp"
#include "Multiply.hpp"
#include "wrapperMPI.hpp"
#include "mltply.hpp"

/**
 * @file   Multiply.c
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Kota Ido (The University of Tokyo)
 *
 * @brief  File for giving multiplying functions to the wave vectors for TPQ and TE method
 *
 */

/**
 * @brief Function of calculating the i-th step norm as @f[ N^{(i)} = |L v_1^{(i)}-v_0^{(i)}/N_s | @f]
 * and update the i+1-th wave vector as @f[ v_0^{(i+1)} = \frac{L v_1^{(i)}-v_0^{(i)}/N_s}{N^{(i)}} @f] for TPQ calculation
 * @param X [in] data list for calculation (idim_max and NsiteMPI)
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int Multiply
(
  struct BindStruct *X
)
{
  long int i, i_max;
  double Ns;
  int rand_i;

  i_max = X->Check.idim_max;
  Ns = 1.0*X->Def.NsiteMPI;
  // mltply is in expec_energy.c v0=H*v1
#pragma omp parallel for default(none) private(i,rand_i)  \
  shared(v0, v1,NumAve) firstprivate(i_max, Ns, LargeValue)
  for (i = 1; i <= i_max; i++) {
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      v0[i][rand_i] = LargeValue * v1[i][rand_i] - v0[i][rand_i] / Ns;  //v0=(l-H/Ns)*v1
    }
  }
  NormMPI_dv(i_max, NumAve, v0, global_norm);
#pragma omp parallel for default(none) private(i,rand_i) \
shared(v0,NumAve,global_norm) firstprivate(i_max)
  for (i = 1; i <= i_max; i++) 
    for (rand_i = 0; rand_i < NumAve; rand_i++)
      v0[i][rand_i] = v0[i][rand_i] / global_norm[rand_i];
  return 0;
}
/**
 * @brief  Function of multiplying Hamiltonian for Time Evolution.
 *
 * Make @f$ |v_0 \rangle = |\psi(t+dt) \rangle @f$ from @f$ |v_1 \rangle = | \psi(t) \rangle  @f$ and @f$ |v_0 \rangle = H |\psi(t) \rangle @f$.
 * @param X [in] data list for calculation (idim_max and TimeSlice)
 *
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 */
int MultiplyForTEM
(
  struct BindStruct *X,
  std::complex<double> **v2
)
{
  long int i, i_max;
  int coef;
  double dnorm = 0.0;
  std::complex<double> tmp1 = 1.0;
  std::complex<double> tmp2 = 0.0;
  double dt = X->Def.Param.TimeSlice;

  //Make |v0> = |psi(t+dt)> from |v1> = |psi(t)> and |v0> = H |psi(t)>
  i_max = X->Check.idim_max;
  // mltply is in expec_energy.c v0=H*v1
  if (dt < pow(10.0, -14)) {
#pragma omp parallel for default(none) private(i) \
shared(I,v0, v1, v2) firstprivate(i_max, dt, tmp2)
    for (i = 1; i <= i_max; i++) {
      tmp2 = v0[i][0];
      v0[i][0] = v1[i][0];  //v0=(1-i*dt*H)*v1
      v1[i][0] = tmp2;
      v2[i][0] = 0.0 + I * 0.0;
    }
    mltply(X, 1, v2, v1);
  }
  else {
    tmp1 *= -I * dt;
#pragma omp parallel for default(none) private(i) \
shared(v0, v1, v2,I) firstprivate(i_max, dt, tmp1, tmp2)
    for (i = 1; i <= i_max; i++) {
      tmp2 = v0[i][0];
      v0[i][0] = v1[i][0] + tmp1 * tmp2;  //v0=(1-i*dt*H)*v1
      v1[i][0] = tmp2;
      v2[i][0] = 0.0 + I * 0.0;
    }
    for (coef = 2; coef <= X->Def.Param.ExpandCoef; coef++) {
      tmp1 *= -I * dt / (std::complex<double>)coef;
      //v2 = H*v1 = H^coef |psi(t)>
      mltply(X, 1, v2, v1);

#pragma omp parallel for default(none) private(i) shared(I, v0, v1, v2) \
firstprivate(i_max, tmp1, myrank)
      for (i = 1; i <= i_max; i++) {
        v0[i][0] += tmp1 * v2[i][0];
        v1[i][0] = v2[i][0];
        v2[i][0] = 0.0 + I * 0.0;
      }
    }
  }
  dnorm = 0.0;
#pragma omp parallel for default(none) reduction(+: dnorm) private(i) shared(v0) \
firstprivate(i_max, dt)
  for (i = 1; i <= i_max; i++) {
    dnorm += real(conj(v0[i][0])*v0[i][0]);
  }
  dnorm = SumMPI_d(dnorm);
  dnorm = sqrt(dnorm);
  global_norm[0] = dnorm;
#pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max, dnorm)
  for (i = 1; i <= i_max; i++) {
    v0[i][0] = v0[i][0] / dnorm;
  }
  return 0;
}
