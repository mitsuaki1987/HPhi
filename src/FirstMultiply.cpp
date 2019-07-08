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
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @version 0.1
 * @brief  Multiplication @f$ v_0 = H v_1 @f$ at the first step for TPQ mode (@f$ v_1 @f$ is the random or inputted vector).
 *
 */

#include "FirstMultiply.hpp"
#include "expec_energy_flct.hpp"
#include "common/setmemory.hpp"
#include "wrapperMPI.hpp"
#include "CalcTime.hpp"
#include "mltplyCommon.hpp"
#include "expec_cisajs.hpp"
#include "expec_cisajscktaltdc.hpp"
#include "dSFMT.hpp"
#include "global.hpp"
#include "log.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
///
/// \brief Multiplication @f$ v_0 = H v_1 @f$ at the first step for TPQ mode (@f$ v_1 @f$ is the random or inputted vector).
/// \param rand_i [in] A rundom number seed for giving the initial vector @f$ v_1 @f$.
/// \param X [in] Struct to get information of the vector @f$ v_1 @f$ for the first step calculation.
/// \retval -1 fail the multiplication @f$ v_0 = H v_1 @f$.
/// \retval 0 succeed the multiplication @f$ v_0 = H v_1 @f$.
/// \version 0.1
/// \author Takahiro Misawa (The University of Tokyo)
/// \author Kazuyoshi Yoshimi (The University of Tokyo)
int FirstMultiply() {

  long int i, i_max;
  double dnorm;
  double Ns;
  long int u_long_i;
  dsfmt_t dsfmt;
  int mythread, rand_i, iret;

  Ns = 1.0*Def::NsiteMPI;
  i_max = Check::idim_max;

  for (rand_i = 0; rand_i < Step::NumAve; rand_i++) {
#pragma omp parallel default(none) private(i, mythread, u_long_i, dsfmt) \
shared(Wave::v0,Wave::v1,MP::nthreads,MP::myrank,rand_i,MP::STDOUT,i_max,Def::initial_iv,Def::iInitialVecType)
  {
#pragma omp for
    for (i = 1; i <= i_max; i++) {
      Wave::v0[i][rand_i] = 0.0;
      Wave::v1[i][rand_i] = 0.0;
    }
    /*
    Initialise MT
    */
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
    u_long_i = 123432 + (rand_i + 1)*labs(Def::initial_iv) + mythread + MP::nthreads * MP::myrank;
    dsfmt_init_gen_rand(&dsfmt, u_long_i);

    if (Def::iInitialVecType == 0) {

      StartTimer(3101);
#pragma omp for
      for (i = 1; i <= i_max; i++)
        Wave::v1[i][rand_i] = std::complex<double>
        (2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5),
          2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5));
      }/*if (Def::iInitialVecType == 0)*/
      else {
#pragma omp for
        for (i = 1; i <= i_max; i++)
          Wave::v1[i][rand_i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
      }
      StopTimer(3101);
    }/*#pragma omp parallel*/
    /*
    Normalize v
    */
  dnorm = 0.0;
#pragma omp parallel for default(none) private(i) \
shared(Wave::v1, i_max, rand_i) reduction(+:dnorm)
      for (i = 1; i <= i_max; i++)
        dnorm += real(conj(Wave::v1[i][rand_i])*Wave::v1[i][rand_i]);
    dnorm = SumMPI_d(dnorm);
    dnorm = sqrt(dnorm);
    Step::global_1st_norm[rand_i] = dnorm;
#pragma omp parallel for default(none) private(i) shared(Wave::v1,rand_i,i_max, dnorm)
    for (i = 1; i <= i_max; i++) Wave::v1[i][rand_i] = Wave::v1[i][rand_i] / dnorm;
  }/*for (rand_i = 0; rand_i < Step::NumAve; rand_i++)*/

  TimeKeeperWithRandAndStep("%s_TimeKeeper.dat", 
    "set %d step %d:TPQ begins: %s", "a", rand_i, Step::step_i);
  /**@brief
Compute expectation value at infinite temperature
*/
  Def::istep = 0;
  StartTimer(3300);
  iret = expec_cisajs(Step::NumAve, Wave::v0, Wave::v1);
  StopTimer(3300);
  if (iret != 0) return -1;

  StartTimer(3400);
  iret = expec_cisajscktaltdc(Step::NumAve, Wave::v0, Wave::v1);
  StopTimer(3400);
  if (iret != 0) return -1;

#pragma omp parallel for default(none) private(i,rand_i) shared(Wave::v0,Wave::v1,i_max,Step::NumAve)
  for (i = 1; i <= i_max; i++) 
    for (rand_i = 0; rand_i < Step::NumAve; rand_i++) Wave::v0[i][rand_i] = Wave::v1[i][rand_i];
  StartTimer(3102);
  if(expec_energy_flct(Step::NumAve, Wave::v0, Wave::v1) !=0){
    StopTimer(3102);
    return -1;
  }
  StopTimer(3102);

  for (rand_i = 0; rand_i < Step::NumAve; rand_i++) {
#pragma omp parallel for default(none) private(i) \
shared(Wave::v0, Wave::v1,rand_i, i_max, Ns, Step::LargeValue, MP::myrank)
    for (i = 1; i <= i_max; i++) {
      Wave::v0[i][rand_i] = Step::LargeValue * Wave::v1[i][rand_i] - Wave::v0[i][rand_i] / Ns;
    }

    dnorm = 0.0;
#pragma omp parallel for default(none) private(i) shared(Wave::v0,rand_i,i_max) reduction(+:dnorm)
    for (i = 1; i <= i_max; i++)
      dnorm += real(conj(Wave::v0[i][rand_i])*Wave::v0[i][rand_i]);
    dnorm = SumMPI_d(dnorm);
    dnorm = sqrt(dnorm);
    Step::global_norm[rand_i] = dnorm;
#pragma omp parallel for default(none) private(i) shared(Wave::v0,rand_i,i_max, dnorm)
    for (i = 1; i <= i_max; i++) Wave::v0[i][rand_i] = Wave::v0[i][rand_i] / dnorm;
  }/*for (rand_i = 0; rand_i < Step::NumAve; rand_i++)*/
  TimeKeeperWithRandAndStep("%s_TimeKeeper.dat", "set %d step %d:TPQ finishes: %s", "a", rand_i, Step::step_i);
  return 0;
}
