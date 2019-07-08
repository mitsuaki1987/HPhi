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
/*-------------------------------------------------------------*/
/**@file
@brief Print logo mark and version number
*/
#include <cstdio>
#include "global.hpp"
/**
@brief Print logo mark and version number
*/
void splash(){
   int ver_maj =
#include "version_major.hpp"
;
   int ver_min =
#include "version_miner.hpp"
;
   int ver_pat =
#include "version_patch.hpp"
;

  fprintf(MP::STDOUT, "                                                                \n");
  fprintf(MP::STDOUT, "      ,ammmmmmmmmmmmmmb,,        Welcome to the                 \n");
  fprintf(MP::STDOUT, "    ,@@` dm          mb  ===m                                    \n");
  fprintf(MP::STDOUT, "  ,@@` d@@@@@@@@@@@@@@@@b Pm,    @@          @@        @@       \n");
  fprintf(MP::STDOUT, " d@  d@@@ @@@ @@@@@@ @@@@b ~@a   @@          @@     @@@@@@@@    \n");
  fprintf(MP::STDOUT, "d@   @@@@ ^^^ @@@@ m m @@@   @,  @@          @@   @@@  @@  @@@  \n");
  fprintf(MP::STDOUT, "@    @@@@_@@@_@@@@mm mm@@@   @|  @@mmmmmmmmmm@@  @@    @@    @@ \n");
  fprintf(MP::STDOUT, "P@    9@@@@@@@@@@@@@@@@@P    @~  @@@@@@@@@@@@@@  @@    @@    @@ \n");
  fprintf(MP::STDOUT, " @@      ~~9@@@@@@PPP~      @P   @@          @@   @@@  @@  @@@  \n");
  fprintf(MP::STDOUT, "  ~@@b      @@@@@@@      ,@@~    @@          @@     @@@@@@@@    \n");
  fprintf(MP::STDOUT, "    ~@@@m,,@@@@@@@@@  ,m@~`      @@          @@        @@       \n");
  fprintf(MP::STDOUT, "        ~~9@@@@@@@@@  ~                                         \n");
  fprintf(MP::STDOUT, "           9@P~~~9@P             Version %d.%d.%d    \n", ver_maj, ver_min, ver_pat);
  fprintf(MP::STDOUT, "                                                                \n");

}/*void splash()*/
