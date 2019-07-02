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

#include "common/setmemory.hpp"
#include "xsetmem.hpp"
#include "wrapperMPI.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include <iostream>
/**
 * @file
 *
 * @brief  Set size of memories to be needed for calculation.
 * @version 2.0
 * @version 1.2
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */

///
/// \brief Set size of memories headers of output files.
/// Output: CDataFileHead, CParaFileHead
/// \version 0.1
void setmem_HEAD()
{
  Def::CDataFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
  Def::CParaFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
}

///
/// \brief Set size of memories for Def and Phys in BindStruct.
/// \param X [in,out] BindStruct to get information of Def and Phys structs.
/// \param xBoost [in,out] Struct for Boost mode.
/// \version 0.1
void setmem_def()
{
  Def::Tpow = li_1d_allocate(2 * Def::Nsite + 2);
  Def::OrgTpow = li_1d_allocate(2 * Def::Nsite + 2);
  Def::SiteToBit = li_1d_allocate(Def::Nsite + 1);
  Def::LocSpn = i_1d_allocate(Def::Nsite);
  Phys::spin_real_cor = d_1d_allocate(Def::Nsite * Def::Nsite);
  Phys::charge_real_cor = d_1d_allocate(Def::Nsite * Def::Nsite);
  Phys::loc_spin_z = d_1d_allocate(Def::Nsite * Def::Nsite);
  Def::EDChemi = i_1d_allocate(Def::NInterAll + Def::NTransfer);
  Def::EDSpinChemi = i_1d_allocate(Def::NInterAll + Def::NTransfer);
  Def::EDParaChemi = d_1d_allocate(Def::NInterAll + Def::NTransfer);
  Def::GeneralTransfer = i_2d_allocate(Def::NTransfer, 4);
  Def::ParaGeneralTransfer = cd_1d_allocate(Def::NTransfer);

  if (Def::iCalcType == TimeEvolution) {
    Def::EDGeneralTransfer = i_2d_allocate(Def::NTransfer + Def::NTETransferMax, 4);
    Def::EDParaGeneralTransfer = cd_1d_allocate(Def::NTransfer + Def::NTETransferMax);
  } else {
    Def::EDGeneralTransfer = i_2d_allocate(Def::NTransfer, 4);
    Def::EDParaGeneralTransfer = cd_1d_allocate(Def::NTransfer);
  }

  Def::CoulombIntra = i_2d_allocate(Def::NCoulombIntra, 1);
  Def::ParaCoulombIntra = d_1d_allocate(Def::NCoulombIntra);
  Def::CoulombInter = i_2d_allocate(Def::NCoulombInter + Def::NIsingCoupling, 2);
  Def::ParaCoulombInter = d_1d_allocate(Def::NCoulombInter + Def::NIsingCoupling);
  Def::HundCoupling = i_2d_allocate(Def::NHundCoupling + Def::NIsingCoupling, 2);
  Def::ParaHundCoupling = d_1d_allocate(Def::NHundCoupling + Def::NIsingCoupling);
  Def::PairHopping = i_2d_allocate(Def::NPairHopping, 2);
  Def::ParaPairHopping = d_1d_allocate(Def::NPairHopping);
  Def::ExchangeCoupling = i_2d_allocate(Def::NExchangeCoupling, 2);
  Def::ParaExchangeCoupling = d_1d_allocate(Def::NExchangeCoupling);
  Def::PairLiftCoupling = i_2d_allocate(Def::NPairLiftCoupling, 2);
  Def::ParaPairLiftCoupling = d_1d_allocate(Def::NPairLiftCoupling);

  Def::InterAll = i_2d_allocate(Def::NInterAll, 8);
  Def::ParaInterAll = cd_1d_allocate(Def::NInterAll);

  Def::CisAjt = i_2d_allocate(Def::NCisAjt, 4);
  Def::CisAjtCkuAlvDC = i_2d_allocate(Def::NCisAjtCkuAlvDC, 8);

  Def::NSingleExcitationOperator = i_1d_allocate(Def::NNSingleExcitationOperator);
  Def::SingleExcitationOperator = (int***)malloc(sizeof(int**)*Def::NNSingleExcitationOperator);
  Def::ParaSingleExcitationOperator = (std::complex<double>**)malloc(
    sizeof(std::complex<double>*)*Def::NNSingleExcitationOperator);
  Def::NPairExcitationOperator = i_1d_allocate(Def::NNPairExcitationOperator);
  Def::PairExcitationOperator = (int***)malloc(sizeof(int**)*Def::NNPairExcitationOperator);
  Def::ParaPairExcitationOperator = (std::complex<double>**)malloc(
    sizeof(std::complex<double>*)*Def::NNPairExcitationOperator);

  Def::ParaLaser = d_1d_allocate(Def::NLaser);

  Boost::list_6spin_star = i_2d_allocate(Boost::R0 * Boost::num_pivot, 7);
  Boost::list_6spin_pair = i_3d_allocate(Boost::R0 * Boost::num_pivot, 7, 15);
  Boost::arrayJ = cd_3d_allocate(Boost::NumarrayJ, 3, 3);

  int NInterAllSet;
  NInterAllSet = (Def::iCalcType == TimeEvolution) ? Def::NInterAll + Def::NTEInterAllMax : Def::NInterAll;
  Def::InterAll_OffDiagonal = i_2d_allocate(NInterAllSet, 8);
  Def::ParaInterAll_OffDiagonal = cd_1d_allocate(NInterAllSet);
  Def::InterAll_Diagonal = i_2d_allocate(NInterAllSet, 4);
  Def::ParaInterAll_Diagonal = d_1d_allocate(NInterAllSet);

  if (Def::iCalcType == TimeEvolution) {
    Def::TETime = d_1d_allocate(Def::NTETimeSteps);
    //Time-dependent Transfer
    Def::NTETransfer = i_1d_allocate(Def::NTETimeSteps);
    Def::NTETransferDiagonal = i_1d_allocate(Def::NTETimeSteps);
    Def::TETransfer = i_3d_allocate(Def::NTETimeSteps, Def::NTETransferMax, 4);
    Def::TETransferDiagonal = i_3d_allocate(Def::NTETimeSteps, Def::NTETransferMax, 2);
    Def::ParaTETransfer = cd_2d_allocate(Def::NTETimeSteps, Def::NTETransferMax);
    Def::ParaTETransferDiagonal = d_2d_allocate(Def::NTETimeSteps, Def::NTETransferMax);
    //Time-dependent InterAll
    Def::NTEInterAll = i_1d_allocate(Def::NTETimeSteps);
    Def::NTEInterAllDiagonal = i_1d_allocate(Def::NTETimeSteps);
    Def::TEInterAll = i_3d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax, 8);
    Def::TEInterAllDiagonal = i_3d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax, 4);
    Def::ParaTEInterAll = cd_2d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax);
    Def::ParaTEInterAllDiagonal = d_2d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax);
    Def::NTEInterAllOffDiagonal = i_1d_allocate(Def::NTETimeSteps);
    Def::TEInterAllOffDiagonal = i_3d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax, 8);
    Def::ParaTEInterAllOffDiagonal = cd_2d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax);
    //Time-dependent Chemi generated by InterAll diagonal components
    Def::NTEChemi = i_1d_allocate(Def::NTETimeSteps);
    Def::TEChemi = i_2d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax);
    Def::SpinTEChemi = i_2d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax);
    Def::ParaTEChemi = d_2d_allocate(Def::NTETimeSteps, Def::NTEInterAllMax);
  }
}
///
/// \brief Set size of memories for vectors(vg, v0, v1, v2, vec, alpha, beta), lists (list_1, list_2_1, list_2_2, list_Diagonal) and Phys(BindStruct.PhysList) struct in the case of Full Diag mode.
/// \retval -1 Fail to set memories.
/// \retval 0 Normal to set memories.
/// \version 0.1
int setmem_large()
{
  int nstate;

  if (GetlistSize() == TRUE) {
    list_1 = li_1d_allocate(Check::idim_max + 1);
    list_2_1 = li_1d_allocate(Large::SizeOflist_2_1);
    list_2_2 = li_1d_allocate(Large::SizeOflist_2_2);
    if (list_1 == NULL
      || list_2_1 == NULL
      || list_2_2 == NULL
      ) {
      return -1;
    }
  }

  list_Diagonal = d_1d_allocate(Check::idim_max + 1);

  if (Def::iCalcType == FullDiag) {
    nstate = Check::idim_max;
  }
  else if (Def::iCalcType == CG) {
    nstate = Def::k_exct;
  }
  else if (Def::iCalcType == TPQCalc) {
    nstate = NumAve;
  }
  else {
    nstate = 1;
  }
  v0 = cd_2d_allocate(Check::idim_max + 1, nstate);
  v1 = cd_2d_allocate(Check::idim_max + 1, nstate);
#ifdef __MPI
  long int MAXidim_max;
  MAXidim_max = MaxMPI_li(Check::idim_max);
  if (GetlistSize() == TRUE) list_1buf = li_1d_allocate(MAXidim_max + 1);
  v1buf = cd_2d_allocate(MAXidim_max + 1, nstate);
#else
  if (Def::iCalcType == CG) v1buf = cd_2d_allocate(Check::idim_max + 1, nstate);
#endif // MPI

  Phys::num_down = d_1d_allocate(nstate);
  Phys::num_up = d_1d_allocate(nstate);
  Phys::num = d_1d_allocate(nstate);
  Phys::num2 = d_1d_allocate(nstate);
  Phys::energy = d_1d_allocate(nstate);
  Phys::var = d_1d_allocate(nstate);
  Phys::doublon = d_1d_allocate(nstate);
  Phys::doublon2 = d_1d_allocate(nstate);
  Phys::Sz = d_1d_allocate(nstate);
  Phys::Sz2 = d_1d_allocate(nstate);
  Phys::s2 = d_1d_allocate(nstate);

  fprintf(stdoutMPI, "%s", "\n######  LARGE ALLOCATE FINISH !  ######\n\n");
  return 0;
}
///
/// \brief Set size of lists for the canonical ensemble.
/// Input: DefineList.iFlgGeneralSpin, DefineList.iCalcModel, DefineList.Nsite, CheckList.sdim, DefineList.Tpow, DefineList.SiteToBit\n
/// Output: LargeList.SizeOflist_2_1, LargeList.SizeOflist_2_2, LargeList.SizeOflistjb
/// \retval TRUE: Normally finished
/// \retval FALSE: Unnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetlistSize()
{
  switch (Def::iCalcModel) {
  case Spin:
  case Hubbard:
  case HubbardNConserved:
  case Kondo:
  case KondoGC:
    if (Def::iFlgGeneralSpin == FALSE) {
      if (Def::iCalcModel == Spin && Def::Nsite % 2 == 1) {
        Large::SizeOflist_2_1 = Check::sdim * 2 + 2;
      }
      else {
        Large::SizeOflist_2_1 = Check::sdim + 2;
      }
      Large::SizeOflist_2_2 = Check::sdim + 2;
      Large::SizeOflistjb = Check::sdim + 2;
    }
    else {//for spin-canonical general spin
      Large::SizeOflist_2_1 = Check::sdim + 2;
      Large::SizeOflist_2_2 =
        Def::Tpow[Def::Nsite - 1] * Def::SiteToBit[Def::Nsite - 1] / Check::sdim + 2;
      Large::SizeOflistjb =
        Def::Tpow[Def::Nsite - 1] * Def::SiteToBit[Def::Nsite - 1] / Check::sdim + 2;
    }
    break;
  default:
    Large::SizeOflistjb = 1;
    return FALSE;
  }
  return TRUE;
}
/**
@page page_setmem Malloc vectors

 To set memory, we use ```?malloc%``` function defined in @c mfmemmory.h,
 where ```?``` indicates the type of the array and ```%``` means the dimension.

 For example, ```char_malloc2(N1, N2)``` function sets the size of memories N1@f$ \times @f$ N2 characters to two dimensional array X.

 To set memories to global arrays, we prepare two functions, setmem_def() and setmem_large() functions.

 - setmem_def()

    In this function, the memories of the arrays which do not have large memory are stored.

    Arrays for defining interactions and correlation functions are mainly defined.

 - setmem_large()

    In this function, the memories of the arrays which have large memory are stored.

    Arrays for defining Hamiltonian and vectors are mainly defined.

@sa setmem_def(), setmem_large()
*/
