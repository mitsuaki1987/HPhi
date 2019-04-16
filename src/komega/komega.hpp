/*
ISSP Math Library - A library for solving linear systems in materials science
Copyright (C) 2016 Mitsuaki Kawamura

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

For more details, See ‘COPYING.LESSER’ in the root directory of this library.
*/

#include <complex>
#pragma once

extern "C" {
  extern void komega_bicg_init(int *ndim, int *nl, int *nz, std::complex<double> *x,
    std::complex<double> *z, int *itermax, double *threshold, int *comm);
  extern void komega_bicg_restart(int *ndim, int *nl, int *nz, std::complex<double> *x,
    std::complex<double> *z, int *itermax, double *threshold,
    int *status, int *iter_old, std::complex<double> *v2,
    std::complex<double> *v12, std::complex<double> *v4, std::complex<double> *v14,
    std::complex<double> *alpha_save, std::complex<double> *beta_save,
    std::complex<double> *z_seed, std::complex<double> *r_l_save, int *comm);
  extern void komega_bicg_update(std::complex<double> *v12, std::complex<double> *v2,
    std::complex<double> *v14, std::complex<double> *v4,
    std::complex<double> *x, std::complex<double> *r_l, int *status);
  extern void komega_bicg_getcoef(std::complex<double> *alpha_save, std::complex<double> *beta_save,
    std::complex<double> *z_seed, std::complex<double> *r_l_save);
  extern void komega_bicg_getvec(std::complex<double> *r_old, std::complex<double> *r_tilde_old);
  extern void komega_bicg_getresidual(double *res);
  extern void komega_bicg_finalize();
}
