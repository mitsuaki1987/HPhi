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
#pragma once
#include "Common.hpp"

long int Binomial(
         int n,
         int k,
         long int **comb,
         int Nsite
         );

int sz(
       struct BindStruct *X,
       long int *list_1_,
       long int *list_2_1_,
       long int *list_2_2_
       );

int Read_sz
(
 struct BindStruct *X,
 const long int irght,
 const long int ilft,
 const long int ihfbit,
 long int *i_max
 );
