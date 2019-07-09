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
@brief Standard mode for the pyrochlore lattice
*/
#include "StdFace_vals.hpp"
#include "StdFace_ModelUtil.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <cstring>

/**
@brief Setup a Hamiltonian for the Pyrochlore structure
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace::Pyrochlore()
{
  int isite, jsite, isiteUC;
  int iL, iW, iH, kCell;
  FILE *fp;
  std::complex<double> Cphase;
  double dR[3];

  /**@brief
  (1) Compute the shape of the super-cell and sites in the super-cell
  */
  fp = fopen("lattice.xsf", "w");
  /**/
  StdI::NsiteUC = 4;
  /**/
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace::PrintVal_d("a", &StdI::a, 1.0);
  StdFace::PrintVal_d("Wlength", &StdI::length[0], StdI::a);
  StdFace::PrintVal_d("Llength", &StdI::length[1], StdI::a);
  StdFace::PrintVal_d("Hlength", &StdI::length[2], StdI::a);
  StdFace::PrintVal_d("Wx", &StdI::direct[0][0], 0.0);
  StdFace::PrintVal_d("Wy", &StdI::direct[0][1], 0.5 * StdI::length[0]);
  StdFace::PrintVal_d("Wz", &StdI::direct[0][2], 0.5 * StdI::length[0]);
  StdFace::PrintVal_d("Lx", &StdI::direct[1][0], 0.5 * StdI::length[1]);
  StdFace::PrintVal_d("Ly", &StdI::direct[1][1], 0.0);
  StdFace::PrintVal_d("Lz", &StdI::direct[1][2], 0.5 * StdI::length[1]);
  StdFace::PrintVal_d("Hx", &StdI::direct[2][0], 0.5 * StdI::length[2]);
  StdFace::PrintVal_d("Hy", &StdI::direct[2][1], 0.5 * StdI::length[2]);
  StdFace::PrintVal_d("Hz", &StdI::direct[2][2], 0.0);

  StdFace::PrintVal_d("phase0", &StdI::phase[0], 0.0);
  StdFace::PrintVal_d("phase1", &StdI::phase[1], 0.0);
  StdFace::PrintVal_d("phase2", &StdI::phase[2], 0.0);
  /**/
  StdFace::InitSite(fp, 3);
  StdI::tau[0][0] = 0.0; StdI::tau[0][1] = 0.0; ; StdI::tau[0][2] = 0.0;
  StdI::tau[1][0] = 0.5; StdI::tau[1][1] = 0.0; ; StdI::tau[1][2] = 0.0;
  StdI::tau[2][0] = 0.0; StdI::tau[2][1] = 0.5; ; StdI::tau[2][2] = 0.0;
  StdI::tau[3][0] = 0.0; StdI::tau[3][1] = 0.0; ; StdI::tau[3][2] = 0.5;
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
    StdFace::InputSpinNN(StdI::J, StdI::JAll, StdI::J0p, StdI::J0pAll, "J0'");
    StdFace::InputSpinNN(StdI::J, StdI::JAll, StdI::J1p, StdI::J1pAll, "J1'");
    StdFace::InputSpinNN(StdI::J, StdI::JAll, StdI::J2p, StdI::J2pAll, "J2'");
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
    StdFace::NotUsed_c("t''", StdI::tpp);
    StdFace::NotUsed_d("V", StdI::V);
    StdFace::NotUsed_d("V0", StdI::V0);
    StdFace::NotUsed_d("V1", StdI::V1);
    StdFace::NotUsed_d("V'", StdI::Vp);
  }/*if (strcmp(StdI::model, "spin") == 0 )*/
  else {
    StdFace::PrintVal_d("mu", &StdI::mu, 0.0);
    StdFace::PrintVal_d("U", &StdI::U, 0.0);
    StdFace::InputHopp(StdI::t, &StdI::t0, "t0");
    StdFace::InputHopp(StdI::t, &StdI::t1, "t1");
    StdFace::InputHopp(StdI::t, &StdI::t2, "t2");
    StdFace::InputHopp(StdI::t, &StdI::t0p, "t0'");
    StdFace::InputHopp(StdI::t, &StdI::t1p, "t1'");
    StdFace::InputHopp(StdI::t, &StdI::t2p, "t2'");
    StdFace::InputCoulombV(StdI::V, &StdI::V0, "V0");
    StdFace::InputCoulombV(StdI::V, &StdI::V1, "V1");
    StdFace::InputCoulombV(StdI::V, &StdI::V2, "V2");
    StdFace::InputCoulombV(StdI::V, &StdI::V0p, "V0'");
    StdFace::InputCoulombV(StdI::V, &StdI::V1p, "V1'");
    StdFace::InputCoulombV(StdI::V, &StdI::V2p, "V2'");
    /**/
    StdFace::NotUsed_J("J0", StdI::J0All, StdI::J0);
    StdFace::NotUsed_J("J1", StdI::J1All, StdI::J1);
    StdFace::NotUsed_J("J2", StdI::J2All, StdI::J2);
    StdFace::NotUsed_J("J0'", StdI::J0pAll, StdI::J0p);
    StdFace::NotUsed_J("J1'", StdI::J1pAll, StdI::J1p);
    StdFace::NotUsed_J("J2'", StdI::J2pAll, StdI::J2p);
    StdFace::NotUsed_J("J''", StdI::JppAll, StdI::Jpp);
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
  if(strcmp(StdI::model, "spin") == 0 )
    for (isite = 0; isite < StdI::nsite; isite++) StdI::locspinflag[isite] = StdI::S2;
  else if(strcmp(StdI::model, "hubbard") == 0 )
    for (isite = 0; isite < StdI::nsite; isite++) StdI::locspinflag[isite] = 0;
  else 
    for (iL = 0; iL < StdI::nsite / 2; iL++) {
      StdI::locspinflag[iL] = StdI::S2;
      StdI::locspinflag[iL + StdI::nsite / 2] = 0;
    }
  /**@brief
  (5) Set Transfer & Interaction
  */
  for (kCell = 0; kCell < StdI::NCell; kCell++){
    /**/
    iW = StdI::Cell[kCell][0];
    iL = StdI::Cell[kCell][1];
    iH = StdI::Cell[kCell][2];
    /*
     (1) Local term
    */
    isite = StdI::NsiteUC * kCell;
    if (strcmp(StdI::model, "kondo") == 0) isite += StdI::nsite / 2;
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      for (isiteUC = 0; isiteUC < StdI::NsiteUC; isiteUC++) {
        StdFace::MagField(StdI::S2, -StdI::h, -StdI::Gamma, isite + isiteUC);
        StdFace::GeneralJ(StdI::D, StdI::S2, StdI::S2, isite + isiteUC, isite + isiteUC);
      }/*for (isiteUC = 0; isiteUC < StdI::NsiteUC; isiteUC++)*/
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      for (isiteUC = 0; isiteUC < StdI::NsiteUC; isiteUC++) {
        StdFace::HubbardLocal(StdI::mu, -StdI::h, -StdI::Gamma, StdI::U, isite + isiteUC);
      }/*for (isiteUC = 0; isiteUC < StdI::NsiteUC; isiteUC++)*/
      /**/
      if (strcmp(StdI::model, "kondo") == 0) {
        jsite = StdI::NsiteUC * kCell;
        for (isiteUC = 0; isiteUC < StdI::NsiteUC; isiteUC++) {
          StdFace::GeneralJ(StdI::J, 1, StdI::S2, isite + 3, jsite + isiteUC);
          StdFace::MagField(StdI::S2, -StdI::h, -StdI::Gamma, jsite + isiteUC);
        }/*for (isiteUC = 0; isiteUC < StdI::NsiteUC; isiteUC++)*/
      }/*if (strcmp(StdI::model, "kondo") == 0 )*/
    }/*if (strcmp(StdI::model, "spin") != 0 )*/
    /*
     (2) Intra-Cell along W
    */
    StdFace::FindSite(iW, iL, iH, 0, 0, 0, 0, 1, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0 ) {
      StdFace::GeneralJ(StdI::J0, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t0, isite, jsite, dR);
      StdFace::Coulomb(StdI::V0, isite, jsite);
    }
    /*
     (3) Intra-Cell along L
    */
    StdFace::FindSite(iW, iL, iH, 0, 0, 0, 0, 2, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J1, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t1, isite, jsite, dR);
      StdFace::Coulomb(StdI::V1, isite, jsite);
    }
    /*
     (4) Intra-Cell along H
    */
    StdFace::FindSite(iW, iL, iH, 0, 0, 0, 0, 3, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J2, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t2, isite, jsite, dR);
      StdFace::Coulomb(StdI::V2, isite, jsite);
    }
    /*
     (5) Intra-Cell along L-H
    */
    StdFace::FindSite(iW, iL, iH, 0, 0, 0, 2, 3, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J0p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t0p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V0p, isite, jsite);
    }
    /*
     (6) Intra-Cell along H-W
    */
    StdFace::FindSite(iW, iL, iH, 0, 0, 0, 3, 1, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J1p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t1p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V1p, isite, jsite);
    }
    /*
     (7) Intra-Cell along W-L
    */
    StdFace::FindSite(iW, iL, iH, 0, 0, 0, 1, 2, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J2p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t2p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V2p, isite, jsite);
    }
    /*
     (8) Inter-Cell along W
    */
    StdFace::FindSite(iW, iL, iH, 1, 0, 0, 1, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J0, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t0, isite, jsite, dR);
      StdFace::Coulomb(StdI::V0, isite, jsite);
    }
    /*
     (9) Inter-Cell along L
    */
    StdFace::FindSite(iW, iL, iH, 0, 1, 0, 2, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J1, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t1, isite, jsite, dR);
      StdFace::Coulomb(StdI::V1, isite, jsite);
    }
    /*
     (10) Inter-Cell along H
    */
    StdFace::FindSite(iW, iL, iH, 0, 0, 1, 3, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J2, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t2, isite, jsite, dR);
      StdFace::Coulomb(StdI::V2, isite, jsite);
    }
    /*
    (11) Inter-Cell along L-H
    */
    StdFace::FindSite(iW, iL, iH, 0, -1, 1, 3, 2, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J0p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t0p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V0p, isite, jsite);
    }
    /*
    (12) Intra-Cell along H-W
    */
    StdFace::FindSite(iW, iL, iH, 1, 0, -1, 1, 3, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J1p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t1p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V1p, isite, jsite);
    }
    /*
    (13) Intra-Cell along W-L
    */
    StdFace::FindSite(iW, iL, iH, -1, 1, 0, 2, 1, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI::model, "spin") == 0) {
      StdFace::GeneralJ(StdI::J2p, StdI::S2, StdI::S2, isite, jsite);
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      StdFace::Hopping(Cphase * StdI::t2p, isite, jsite, dR);
      StdFace::Coulomb(StdI::V2p, isite, jsite);
    }
  }/*for (kCell = 0; kCell < StdI::NCell; kCell++)*/

  fclose(fp);
  StdFace::PrintXSF();
  StdFace::PrintGeometry();
}/*void StdFace::Pyrochlore*/
