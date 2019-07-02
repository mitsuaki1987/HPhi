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
/**@file
@brief Compute total number of electrons, spins
*/
#include "Common.hpp"
#include "wrapperMPI.hpp"
/**
@brief Define the number of sites in each PE (DefineList.Nsite).
 Reduce the number of electrons (DefineList.Ne), 
 total Sz (DefineList.Total2Sz) by them in the inter process region 
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int CheckMPI()
{
  int isite;
  int NDimInterPE, SmallDim, SpinNum, ishift;
  int ipivot, isiteMax, isiteMax0;

  /**@brief
  Branch for each model
  <ul>
  */
  Def::NsiteMPI = Def::Nsite;
  Def::Total2SzMPI = Def::Total2Sz;
  switch (Def::iCalcModel) {
  case HubbardGC: /****************************************************/
  case Hubbard:
  case HubbardNConserved:
  case Kondo:
  case KondoGC:

    /**@brief
     <li> For Hubbard & Kondo
     Define local dimension DefineList::Nsite</li>
     <ul>
    */
    NDimInterPE = 1;
    for (isite = Def::NsiteMPI; isite > 0; isite--) {
      if (NDimInterPE == nproc) {
        Def::Nsite = isite;
        break;
      } /*if (NDimInterPE == nproc)*/
      NDimInterPE *= 4;
    } /*for (isite = NsiteMPI; isite > 0; isite--)*/

    if (isite == 0) {
      fprintf(stdoutMPI, "%s", "Error ! The number of PROCESS should be 4-exponent !\n");
      fprintf(stdoutMPI, "        The number of PROCESS : %d\n", nproc);
      NDimInterPE = 1;
      int ismallNproc = 1;
      int ilargeNproc = 1;
      for (isite = Def::NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE > nproc) {
          ilargeNproc = NDimInterPE;
          if (isite > 1)
            ismallNproc = NDimInterPE / 4;
          break;
        }/*if (NDimInterPE > nproc)*/
        NDimInterPE *= 4;
      }/*for (isite = Def::NsiteMPI; isite > 0; isite--)*/
      fprintf(stdoutMPI, "        Set the number of PROCESS as %d or %d.\n", ismallNproc, ilargeNproc);
      return FALSE;
      //return FALSE;
    } /*if (isite == 0)*/

    switch (Def::iCalcModel) /*2 (inner)*/ {

    case Hubbard:
      /**@brief
      <li>For canonical Hubbard
      DefineList::Nup, DefineList::Ndown, and DefineList::Ne should be
      differerent in each PE.</li>
      */
      SmallDim = myrank;
      for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
        SpinNum = SmallDim % 4;
        SmallDim /= 4;
        if (SpinNum == 1 /*01*/) {
          Def::Nup -= 1;
          Def::Ne -= 1;
        }
        else if (SpinNum == 2 /*10*/) {
          Def::Ndown -= 1;
          Def::Ne -= 1;
        }
        else if (SpinNum == 3 /*11*/) {
          Def::Nup -= 1;
          Def::Ndown -= 1;
          Def::Ne -= 2;
        }
      } /*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/

      break;/*case Hubbard:*/

    case HubbardNConserved:
      /**@brief
      <li>For N-conserved canonical Hubbard
      DefineList::Ne should be differerent in each PE.</li>
      */
      SmallDim = myrank;
      for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
        SpinNum = SmallDim % 4;
        SmallDim /= 4;
        if (SpinNum == 1 /*01*/ || SpinNum == 2 /*10*/) Def::Ne -= 1;
        else if (SpinNum == 3 /*11*/) Def::Ne -= 2;
      } /*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/

      break; /*case HubbardNConserved:*/

    case KondoGC:
    case Kondo:
      /**@brief
      <li>For canonical Kondo system
      DefineList::Nup, DefineList::Ndown, and DefineList::Ne should be
      differerent in each PE.</li>
      */
      for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)
        if (Def::LocSpn[isite] != ITINERANT) Def::NLocSpn -= 1;

      if (Def::iCalcModel == Kondo) {
        SmallDim = myrank;
        for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
          SpinNum = SmallDim % 4;
          SmallDim /= 4;
          if (Def::LocSpn[isite] == ITINERANT) {
            if (SpinNum == 1 /*01*/) {
              Def::Nup -= 1;
              Def::Ne -= 1;
            }
            else if (SpinNum == 2 /*10*/) {
              Def::Ndown -= 1;
              Def::Ne -= 1;
            }
            else if (SpinNum == 3 /*11*/) {
              Def::Nup -= 1;
              Def::Ndown -= 1;
              Def::Ne -= 2;
            }
          }
          else {
            fprintf(stdoutMPI, "\n Stop because local spin in the inter process region\n");
            return FALSE;
          }
        }/*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/
      } /*if (Def::iCalcModel == Kondo)*/
      else {
        Def::Nup = 0;
        Def::Ndown = 0;
        Def::Ne = 0;
      }

      break; /*case KondoGC, Kondo*/

    case HubbardGC:
      Def::Nup = 0;
      Def::Ndown = 0;
      Def::Ne = 0;
      Def::Total2Sz = 0;
      break;
    } /*switch (Def::iCalcModel) 2(inner)*/

    break; /*case HubbardGC, Hubbard, HubbardNConserved, Kondo, KondoGC:*/
    /**@brief</ul>*/
  case SpinGC:/********************************************************/
  case Spin:

    if (Def::iFlgGeneralSpin == FALSE) {
      /**@brief
      <li> For 1/2 Spin system,
      define local dimension DefineList::Nsite</li>
      */
      NDimInterPE = 1;
      for (isite = Def::NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          Def::Nsite = isite;
          break;
        }/*if (NDimInterPE == nproc)*/
        NDimInterPE *= 2;
      }/*for (isite = Def::NsiteMPI; isite > 0; isite--)*/

      if (isite == 0) {
        fprintf(stdoutMPI, "%s", "Error ! The number of PROCESS should be 2-exponent !\n");
        fprintf(stdoutMPI, "        The number of PROCESS : %d\n", nproc);
        NDimInterPE = 1;
        int ismallNproc = 1;
        int ilargeNproc = 1;
        for (isite = Def::NsiteMPI; isite > 0; isite--) {
          if (NDimInterPE > nproc) {
            ilargeNproc = NDimInterPE;
            if (isite > 1)
              ismallNproc = NDimInterPE / 2;
            break;
          }/*if (NDimInterPE > nproc)*/
          NDimInterPE *= 2;
        }/*for (isite = Def::NsiteMPI; isite > 0; isite--)*/
        fprintf(stdoutMPI, "        Set the number of PROCESS as %d or %d.\n", ismallNproc, ilargeNproc);
        return FALSE;
      }/*if (isite == 0)*/

      if (Def::iCalcModel == Spin) {
        /*Def::NeMPI = Def::Ne;*/

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
          SpinNum = SmallDim % 2;
          SmallDim /= 2;
          if (SpinNum == 0) {
            Def::Ndown -= 1;
          }
          else {
            Def::Ne -= 1;
            Def::Nup -= 1;
          }
        }/*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/
      }/*if (Def::iCalcModel == Spin)*/

    } /*if (Def::iFlgGeneralSpin == FALSE)*/
    else {/* General Spin */
      /**@brief
      <li> For general Spin system,
      define local dimension DefineList::Nsite</li>
      */
      NDimInterPE = 1;
      for (isite = Def::NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          Def::Nsite = isite;
          break;
        }/*if (NDimInterPE == nproc)*/
        NDimInterPE *= Def::SiteToBit[isite - 1];
      }/*for (isite = Def::NsiteMPI; isite > 0; isite--)*/

      if (isite == 0) {
        fprintf(stdoutMPI, "%s", "Error ! The number of PROCESS is wrong !\n");
        fprintf(stdoutMPI, "        The number of PROCESS : %d\n", nproc);
        NDimInterPE = 1;
        int ismallNproc = 1;
        int ilargeNproc = 1;
        for (isite = Def::NsiteMPI; isite > 0; isite--) {
          if (NDimInterPE > nproc) {
            ilargeNproc = NDimInterPE;
            if (isite > 1)
              ismallNproc = NDimInterPE / Def::SiteToBit[isite - 2];
            break;
          }/*if (NDimInterPE > nproc)*/
          NDimInterPE *= Def::SiteToBit[isite - 1];
        }/*for (isite = Def::NsiteMPI; isite > 0; isite--)*/
        fprintf(stdoutMPI, "        Set the number of PROCESS as %d or %d.\n", ismallNproc, ilargeNproc);
        return FALSE;
      }/*if (isite == 0)*/

      if (Def::iCalcModel == Spin) {
        Def::Total2SzMPI = Def::Total2Sz;

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
          SpinNum = SmallDim % Def::SiteToBit[isite];
          SmallDim /= Def::SiteToBit[isite];

          Def::Total2Sz += Def::SiteToBit[isite] - 1 - 2 * SpinNum;
        }/*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/
      }/*if (Def::iCalcModel == Spin)*/
    }/*if (Def::iFlgGeneralSpin == TRUE)*/

     /**@brief</ul>*/
    break; /*case SpinGC, Spin*/

  default:
    fprintf(stdoutMPI, "Error ! Wrong model !\n");
    return FALSE;
  }/*switch (Def::iCalcModel)*/

  /**@brief
   Check the number of processes for Boost
  */
  if (Boost::flgBoost == 1) {
    isiteMax = Boost::W0;
    ishift = 0;
    for (ipivot = 0; ipivot < Boost::num_pivot; ipivot++) {
      isiteMax0 = Boost::list_6spin_star[ipivot][1]
        + Boost::list_6spin_star[ipivot][2]
        + Boost::list_6spin_star[ipivot][3]
        + Boost::list_6spin_star[ipivot][4]
        + Boost::list_6spin_star[ipivot][5];
      if (ishift > 1) isiteMax0 = Def::NsiteMPI - isiteMax0 - 1 - ishift;
      else isiteMax0 = Def::NsiteMPI - isiteMax0 - 2;
      if (isiteMax0 < isiteMax) isiteMax = isiteMax0;
      if (Boost::list_6spin_star[ipivot][6] == 1) ishift += Boost::ishift_nspin;
    }/*for (ipivot = 0; ipivot < Boost::num_pivot; ipivot++)*/

    NDimInterPE = 1;
    for (isite = 0; isite < isiteMax; isite++) NDimInterPE *= 2;

    if (NDimInterPE < nproc) {
      fprintf(stderr, "\n Error ! in ReadDefFileIdxPara.\n");
      fprintf(stderr, "Too many MPI processes ! It should be <= %d. \n\n", NDimInterPE);
      exitMPI(-1);
    }/*if (NDimInterPE < nproc)*/
  }/*if (Boost::flgBoost == 1)*/

  return TRUE;
}/*void CheckMPI*/
/**
@brief Print infomation of MPI parallelization
Modify Definelist::Tpow in the inter process region
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void CheckMPI_Summary()
{
  int iproc, SmallDim, SpinNum, Nelec;
  int isite;
  long int idimMPI;

  if (Def::iFlgScaLAPACK == 0) {
    fprintf(stdoutMPI, "\n\n######  MPI site separation summary  ######\n\n");
    fprintf(stdoutMPI, "  INTRA process site\n");
    fprintf(stdoutMPI, "    Site    Bit\n");
    for (isite = 0; isite < Def::Nsite; isite++) {
      switch (Def::iCalcModel) {
      case HubbardGC:
      case Hubbard:
      case HubbardNConserved:
      case Kondo:
      case KondoGC:

        fprintf(stdoutMPI, "    %4d    %4d\n", isite, 4);
        break;

      case Spin:
      case SpinGC:

        if (Def::iFlgGeneralSpin == FALSE) {
          fprintf(stdoutMPI, "    %4d    %4d\n", isite, 2);
        }/*if (Def::iFlgGeneralSpin == FALSE)*/
        else {
          fprintf(stdoutMPI, "    %4d    %4ld\n", isite, Def::SiteToBit[isite]);
        }/*if (Def::iFlgGeneralSpin == TRUE)*/

        break;

      } /*switch (Def::iCalcModel)*/
    } /*for (isite = 0; isite < Def::Nsite; isite++)*/

    fprintf(stdoutMPI, "\n  INTER process site\n");
    fprintf(stdoutMPI, "    Site    Bit\n");
    for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
      switch (Def::iCalcModel) {
      case HubbardGC:
      case Hubbard:
      case HubbardNConserved:
      case Kondo:
      case KondoGC:

        fprintf(stdoutMPI, "    %4d    %4d\n", isite, 4);
        break;

      case Spin:
      case SpinGC:

        if (Def::iFlgGeneralSpin == FALSE) {
          fprintf(stdoutMPI, "    %4d    %4d\n", isite, 2);
        }/*if (Def::iFlgGeneralSpin == FALSE) */
        else {
          fprintf(stdoutMPI, "    %4d    %4ld\n", isite, Def::SiteToBit[isite]);
        }/*if (Def::iFlgGeneralSpin == TRUE) */

        break;

      }/*switch (Def::iCalcModel)*/
    }/*for (isite = Def::Nsite; isite < NsiteMPI; isite++)*/

    fprintf(stdoutMPI, "\n  Process element info\n");
    fprintf(stdoutMPI, "    Process       Dimension   Nup  Ndown  Nelec  Total2Sz   State\n");

    for (iproc = 0; iproc < nproc; iproc++) {

      fprintf(stdoutMPI, "    %7d", iproc);

      if (myrank == iproc) idimMPI = Check::idim_max;
      else idimMPI = 0;
      fprintf(stdoutMPI, " %15ld", SumMPI_li(idimMPI));

      if (myrank == iproc) Nelec = Def::Nup;
      else Nelec = 0;
      fprintf(stdoutMPI, "  %4d", SumMPI_i(Nelec));

      if (myrank == iproc) Nelec = Def::Ndown;
      else Nelec = 0;
      fprintf(stdoutMPI, "  %5d", SumMPI_i(Nelec));

      if (myrank == iproc) {
        Nelec = Def::Ne; //Def::Nup
        if (Def::iCalcModel == Spin || Def::iCalcModel == SpinGC) Nelec += Def::Ndown;
      }
      else Nelec = 0;

      fprintf(stdoutMPI, "  %5d", SumMPI_i(Nelec));

      if (myrank == iproc) Nelec = Def::Total2Sz;
      else Nelec = 0;
      fprintf(stdoutMPI, "  %8d   ", SumMPI_i(Nelec));
      /**@brief
       Print the configuration in the inter process region of each PE
       as a binary (excepting general spin) format.
      */
      switch (Def::iCalcModel) {
      case HubbardGC: /****************************************************/
      case Hubbard:
      case HubbardNConserved:
      case Kondo:
      case KondoGC:

        SmallDim = iproc;
        for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
          SpinNum = SmallDim % 4;
          SmallDim /= 4;
          if (SpinNum == 0) fprintf(stdoutMPI, "00");
          else if (SpinNum == 1) fprintf(stdoutMPI, "01");
          else if (SpinNum == 2) fprintf(stdoutMPI, "10");
          else if (SpinNum == 3) fprintf(stdoutMPI, "11");
        } /*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/

        break;

      case Spin:
      case SpinGC:

        SmallDim = iproc;
        if (Def::iFlgGeneralSpin == FALSE) {
          for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
            SpinNum = SmallDim % 2;
            SmallDim /= 2;
            fprintf(stdoutMPI, "%1d", SpinNum);
          }/*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/
        }/*if (Def::iFlgGeneralSpin == FALSE)*/
        else {
          SmallDim = iproc;
          for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++) {
            SpinNum = SmallDim % (int)Def::SiteToBit[isite];
            SmallDim /= Def::SiteToBit[isite];
            fprintf(stdoutMPI, "%1d", SpinNum);
          }/*for (isite = Def::Nsite; isite < Def::NsiteMPI; isite++)*/
        }/*if (Def::iFlgGeneralSpin == TRUE)*/

        break;

      }/*switch (Def::iCalcModel)*/
      fprintf(stdoutMPI, "\n");
    }/*for (iproc = 0; iproc < nproc; iproc++)*/

    Check::idim_maxMPI = SumMPI_li(Check::idim_max);
    fprintf(stdoutMPI, "\n   Total dimension : %ld\n\n", Check::idim_maxMPI);
    if (Check::idim_maxMPI < 1) {
      fprintf(stdoutMPI, "ERROR! Total dimension < 1\n");
      exitMPI(-1);
    }
  }
  else {
    fprintf(stdoutMPI, "\n   Total dimension : %ld\n\n", Check::idim_max);
  }

  /**@brief
    Reset DefineList::Tpow[DefNsite], DefineList::Tpow[DefNsite + 1] ...
    as inter process space
    For Hubbard & Kondo system, define DefineList::OrgTpow which is not
    affected by the number of processes.
  */
  switch (Def::iCalcModel) {
  case HubbardGC: /****************************************************/
  case Hubbard:
  case HubbardNConserved:
  case Kondo:
  case KondoGC:

    Def::Tpow[2 * Def::Nsite] = 1;
    for (isite = 2 * Def::Nsite + 1; isite < 2 * Def::NsiteMPI; isite++)
      Def::Tpow[isite] = Def::Tpow[isite - 1] * 2;

    Def::OrgTpow[0] = 1;
    for (isite = 1; isite < 2 * Def::NsiteMPI; isite++)
      Def::OrgTpow[isite] = Def::OrgTpow[isite - 1] * 2;

    break;

  case SpinGC:/********************************************************/
  case Spin:

    if (Def::iFlgGeneralSpin == FALSE) {

      Def::Tpow[Def::Nsite] = 1;
      for (isite = Def::Nsite + 1; isite < Def::NsiteMPI; isite++)
        Def::Tpow[isite] = Def::Tpow[isite - 1] * 2;

    }/*if (Def::iFlgGeneralSpin == FALSE)*/
    else {

      Def::Tpow[Def::Nsite] = 1;
      for (isite = Def::Nsite + 1; isite < Def::NsiteMPI; isite++)
        Def::Tpow[isite] = Def::Tpow[isite - 1] * Def::SiteToBit[isite - 1];

    }/*if (Def::iFlgGeneralSpin == TRUE)*/
    break;
  } /*switch (Def::iCalcModel)*/
}/*void CheckMPI_Summary*/
