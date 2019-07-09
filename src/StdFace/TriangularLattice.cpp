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
/**@file
@brief Standard mode for the triangular lattice
*/
#include "StdFace_vals.hpp"
#include "StdFace_ModelUtil.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <cstring>

/**
@brief Setup a Hamiltonian for the Triangular lattice
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace::Triangular()
{
  int isite, jsite, kCell;
  int iL, iW;
  FILE *fp;
  std::complex<double> Cphase;
  double dR[3];

  /**@brief
  (1) Compute the shape of the super-cell and sites in the super-cell
  */
  fp = fopen("lattice.gp", "w");
  /**/
  StdI::NsiteUC = 1;
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace::PrintVal_d("a", &StdI::a, 1.0);
  StdFace::PrintVal_d("Wlength", &StdI::length[0], StdI::a);
  StdFace::PrintVal_d("Llength", &StdI::length[1], StdI::a);
  StdFace::PrintVal_d("Wx", &StdI::direct[0][0], StdI::length[0]);
  StdFace::PrintVal_d("Wy", &StdI::direct[0][1], 0.0);
  StdFace::PrintVal_d("Lx", &StdI::direct[1][0], StdI::length[1] * 0.5);
  StdFace::PrintVal_d("Ly", &StdI::direct[1][1], StdI::length[1] * 0.5 * sqrt(3.0));
  
  StdFace::PrintVal_d("phase0", &StdI::phase[0], 0.0);
  StdFace::PrintVal_d("phase1", &StdI::phase[1], 0.0);
  /**/
  StdFace::InitSite(fp, 2);
  StdI::tau[0][0] = 0.0; StdI::tau[0][1] = 0.0; StdI::tau[0][2] = 0.0;
  /**@brief
  (2) check & store parameters of Hamiltonian
  */
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  StdFace::NotUsed_d("K", StdI::K);
  StdFace::PrintVal_d("h", &StdI::h, 0.0);
  StdFace::PrintVal_d("Gamma", &StdI::Gamma, 0.0);
  /**/
  if (strcmp(StdI::model, "spin") == 0 ) {
    StdFace::PrintVal_i("2S", &StdI::S2, 1);
    StdFace::PrintVal_d("D", &StdI::D[2][2], 0.0);
    StdFace::InputSpinNN(StdI::J, StdI::JAll, StdI::J0, StdI::J0All, "J0");
    StdFace::InputSpinNN(StdI::J, StdI::JAll, StdI::J1, StdI::J1All, "J1");
    StdFace::InputSpinNN(StdI::J, StdI::JAll, StdI::J2, StdI::J2All, "J2");
    StdFace::InputSpinNN(StdI::Jp, StdI::JpAll, StdI::J0p, StdI::J0pAll, "J0'");
    StdFace::InputSpinNN(StdI::Jp, StdI::JpAll, StdI::J1p, StdI::J1pAll, "J1'");
    StdFace::InputSpinNN(StdI::Jp, StdI::JpAll, StdI::J2p, StdI::J2pAll, "J2'");
    StdFace::InputSpinNN(StdI::Jpp, StdI::JppAll, StdI::J0pp, StdI::J0ppAll, "J0'");
    StdFace::InputSpinNN(StdI::Jpp, StdI::JppAll, StdI::J1pp, StdI::J1ppAll, "J1'");
    StdFace::InputSpinNN(StdI::Jpp, StdI::JppAll, StdI::J2pp, StdI::J2ppAll, "J2'");
    /**/
    StdFace::NotUsed_d("mu", StdI::mu);
    StdFace::NotUsed_d("U", StdI::U);
    StdFace::NotUsed_c("t", StdI::t);
    StdFace::NotUsed_c("t0", StdI::t0);
    StdFace::NotUsed_c("t1", StdI::t1);
    StdFace::NotUsed_c("t2", StdI::t2);
    StdFace::NotUsed_c("t'", StdI::tp);
    StdFace::NotUsed_c("t0'", StdI::t0p);
    StdFace::NotUsed_c("t1'", StdI::t1p);
    StdFace::NotUsed_c("t2'", StdI::t2p);
    StdFace::NotUsed_c("t''", StdI::tp);
    StdFace::NotUsed_c("t0''", StdI::t0pp);
    StdFace::NotUsed_c("t1''", StdI::t1pp);
    StdFace::NotUsed_c("t2''", StdI::t2pp);
    StdFace::NotUsed_d("V", StdI::V);
    StdFace::NotUsed_d("V0", StdI::V0);
    StdFace::NotUsed_d("V1", StdI::V1);
    StdFace::NotUsed_d("V2", StdI::V2);
    StdFace::NotUsed_d("V'", StdI::Vp);
    StdFace::NotUsed_d("V0'", StdI::V0p);
    StdFace::NotUsed_d("V1'", StdI::V1p);
    StdFace::NotUsed_d("V2'", StdI::V2p);
    StdFace::NotUsed_d("V''", StdI::Vpp);
    StdFace::NotUsed_d("V0''", StdI::V0pp);
    StdFace::NotUsed_d("V1''", StdI::V1pp);
    StdFace::NotUsed_d("V2''", StdI::V2pp);
  }/*if (strcmp(StdI::model, "spin") == 0 )*/
  else {
    StdFace::PrintVal_d("mu", &StdI::mu, 0.0);
    StdFace::PrintVal_d("U", &StdI::U, 0.0);
    StdFace::InputHopp(StdI::t, &StdI::t0, "t0");
    StdFace::InputHopp(StdI::t, &StdI::t1, "t1");
    StdFace::InputHopp(StdI::t, &StdI::t2, "t2");
    StdFace::InputHopp(StdI::tp, &StdI::t0p, "t0'");
    StdFace::InputHopp(StdI::tp, &StdI::t1p, "t1'");
    StdFace::InputHopp(StdI::tp, &StdI::t2p, "t2'");
    StdFace::InputHopp(StdI::tpp, &StdI::t0pp, "t0''");
    StdFace::InputHopp(StdI::tpp, &StdI::t1pp, "t1''");
    StdFace::InputHopp(StdI::tpp, &StdI::t2pp, "t2''");
    StdFace::InputCoulombV(StdI::V, &StdI::V0, "V0");
    StdFace::InputCoulombV(StdI::V, &StdI::V1, "V1");
    StdFace::InputCoulombV(StdI::V, &StdI::V2, "V2");
    StdFace::InputCoulombV(StdI::Vp, &StdI::V0p, "V0'");
    StdFace::InputCoulombV(StdI::Vp, &StdI::V1p, "V1'");
    StdFace::InputCoulombV(StdI::Vp, &StdI::V2p, "V2'");
    StdFace::InputCoulombV(StdI::Vpp, &StdI::V0pp, "V0''");
    StdFace::InputCoulombV(StdI::Vpp, &StdI::V1pp, "V1''");
    StdFace::InputCoulombV(StdI::Vpp, &StdI::V2pp, "V2''");
    /**/
    StdFace::NotUsed_J("J0", StdI::J0All, StdI::J0);
    StdFace::NotUsed_J("J1", StdI::J1All, StdI::J1);
    StdFace::NotUsed_J("J2", StdI::J2All, StdI::J2);
    StdFace::NotUsed_J("J0'", StdI::J0pAll, StdI::J0p);
    StdFace::NotUsed_J("J1'", StdI::J1pAll, StdI::J1p);
    StdFace::NotUsed_J("J2'", StdI::J2pAll, StdI::J2p);
    StdFace::NotUsed_J("J0''", StdI::J0ppAll, StdI::J0pp);
    StdFace::NotUsed_J("J1''", StdI::J1ppAll, StdI::J1pp);
    StdFace::NotUsed_J("J2''", StdI::J2ppAll, StdI::J2pp);
    StdFace::NotUsed_d("D", StdI::D[2][2]);

    if (strcmp(StdI::model, "hubbard") == 0 ) {
      StdFace::NotUsed_i("2S", StdI::S2);
      StdFace::NotUsed_J("J", StdI::JAll, StdI::J);
    }/*if (strcmp(StdI::model, "hubbard") == 0 )*/
    else {
      StdFace::PrintVal_i("2S", &StdI::S2, 1);
      StdFace::InputSpin(StdI::J, StdI::JAll, "J");
    }/*if (model != "hubbard")*/

  }/*if (model != "spin")*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");
  /**@brief
  (3) Set local spin flag (StdI::locspinflag) and
  the number of sites (StdI::nsite)
  */
  StdI::nsite = StdI::NsiteUC * StdI::NCell;
  if (strcmp(StdI::model, "kondo") == 0 ) StdI::nsite *= 2;
  StdI::locspinflag = (int *)malloc(sizeof(int) * StdI::nsite);
  /**/
  if (strcmp(StdI::model, "spin") == 0 )
    for (isite = 0; isite < StdI::nsite; isite++) StdI::locspinflag[isite] = StdI::S2;
  else if (strcmp(StdI::model, "hubbard") == 0 )
    for (isite = 0; isite < StdI::nsite; isite++) StdI::locspinflag[isite] = 0;
  else
    for (iL = 0; iL < StdI::nsite / 2; iL++) {
      StdI::locspinflag[iL] = StdI::S2;
      StdI::locspinflag[iL + StdI::nsite / 2] = 0;
    }
  /**@brief
  (5) Set Transfer & Interaction
  */
  for (kCell = 0; kCell < StdI::NCell; kCell++) {
    /**/
    iW = StdI::Cell[kCell][0];
    iL = StdI::Cell[kCell][1];
    /*
    Local term
    */
    isite = kCell;
    if (strcmp(StdI::model, "kondo") == 0 ) isite += StdI::NCell;
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::MagField(StdI::S2, -StdI::h, -StdI::Gamma, isite);
      StdFace::GeneralJ(StdI::D, StdI::S2, StdI::S2, isite, isite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::HubbardLocal(StdI::mu, -StdI::h, -StdI::Gamma, StdI::U, isite);
      if (strcmp(StdI::model, "kondo") == 0 ) {
        jsite = kCell;
        StdFace::GeneralJ(StdI::J, 1, StdI::S2, isite, jsite);
        StdFace::MagField(StdI::S2, -StdI::h, -StdI::Gamma, jsite);
      }/*if (strcmp(StdI::model, "kondo") == 0 )*/
    }
    /*
    Nearest neighbor along W
    */
    StdFace::SetLabel(fp, iW, iL, 1, 0, 0, 0, &isite, &jsite, 1, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::GeneralJ(StdI::J0, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t0, isite, jsite, dR);
      StdFace::Coulomb(StdI::V0, isite, jsite);
    }
    /*
    Nearest neighbor along L
    */
    StdFace::SetLabel(fp, iW, iL, 0, 1, 0, 0, &isite, &jsite, 1, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::GeneralJ(StdI::J1, StdI::S2, StdI::S2, isite, jsite);
    }
    else {
      StdFace::Hopping(Cphase * StdI::t1, isite, jsite, dR);
      StdFace::Coulomb(StdI::V1, isite, jsite);
    }
    /*
    Nearest neighbor along W - L
    */
    StdFace::SetLabel(fp, iW, iL, 1, - 1, 0, 0, &isite, &jsite, 1, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::GeneralJ(StdI::J2, StdI::S2, StdI::S2, isite, jsite);
    }
    else {
      StdFace::Hopping(Cphase * StdI::t2, isite, jsite, dR);
      StdFace::Coulomb(StdI::V2, isite, jsite);
    }
    /*
    Second nearest neighbor 2W - L
    */
    StdFace::SetLabel(fp, iW, iL, 2, - 1, 0, 0, &isite, &jsite, 2, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::GeneralJ(StdI::J1p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t1p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V1p, isite, jsite);
    }
    /*
    Second nearest neighbor W+L
    */
    StdFace::SetLabel(fp, iW, iL, 1, 1, 0, 0, &isite, &jsite, 2, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::GeneralJ(StdI::J2p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t2p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V2p, isite, jsite);
    }/*if (model != "spin")*/
     /*
     Second nearest neighbor -W+2L
     */
    StdFace::SetLabel(fp, iW, iL, - 1, 2, 0, 0, &isite, &jsite, 2, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::GeneralJ(StdI::J0p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t0p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V0p, isite, jsite);
    }/*if (model != "spin")*/
    /*
    Third neighbor along 2W
    */
    StdFace::SetLabel(fp, iW, iL, 2, 0, 0, 0, &isite, &jsite, 3, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J0pp, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t0pp, isite, jsite, dR);
      StdFace::Coulomb(StdI::V0pp, isite, jsite);
    }
    /*
    Third neighbor along 2L
    */
    StdFace::SetLabel(fp, iW, iL, 0, 2, 0, 0, &isite, &jsite, 3, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J1pp, StdI::S2, StdI::S2, isite, jsite);
    }
    else {
      StdFace::Hopping(Cphase * StdI::t1pp, isite, jsite, dR);
      StdFace::Coulomb(StdI::V1pp, isite, jsite);
    }
    /*
    Nearest neighbor along 2W - 2L
    */
    StdFace::SetLabel(fp, iW, iL, 2, -2, 0, 0, &isite, &jsite, 3, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J2pp, StdI::S2, StdI::S2, isite, jsite);
    }
    else {
      StdFace::Hopping(Cphase * StdI::t2pp, isite, jsite, dR);
      StdFace::Coulomb(StdI::V2pp, isite, jsite);
    }
  }/*for (kCell = 0; kCell < StdI::NCell; kCell++)*/

  fprintf(fp, "plot \'-\' w d lc 7\n0.0 0.0\nend\npause -1\n");
  fclose(fp);
  StdFace::PrintGeometry();
}

