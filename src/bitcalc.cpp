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
#include "bitcalc.hpp"
#include "wrapperMPI.hpp"
#include "DefCommon.hpp"
#include "global.hpp"

/**
@file
@version 0.1, 0.2
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
@brief  File for giving functions of treating bits on the target of Hilbert space.
*/

/** 
@brief function of getting right, left and half bits corresponding to a original Hilbert space.
@version 0.1
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
@return 
*/
int GetSplitBit(
  const int Nsite, //!<[in] total number of sites
  long int* irght, //!<[out] bit to split original Hilbert space into right space @f$2^{(Ns+2)/2}-1@f$
  long int* ilft, //!<[out] bit to split original Hilbert space into left space 
  long int* ihfbit//!<[out] half bit to split original Hilbert space @f$2^{(Ns+2)/2}@f$
) {
  if (Nsite < 1) {
    fprintf(stderr, "%s", "Error: Total Site Number is incorrect.\n");
    return -1;
  }
  *ihfbit = 1;
  *ihfbit = (*ihfbit << (long int)((Nsite + 1) / 2));
  *irght = *ihfbit - 1;
  *ilft = 1;
  *ilft = (*ilft << (long int)Nsite) - 1;
  *ilft = *ilft ^ *irght;
  return 0;
}
/** 
@brief function of splitting original bit into right and left  spaces.
@version 0.1
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
@return 
*/
int GetSplitBitByModel(
  const int Nsite,  //!<[in] a total number of sites 
  const int iCalcModel,  //!<[in] Calc model defined in CalcMode file
  long int* irght,  //!<[out] bit to split original  space into right space
  long int* ilft,  //!<[out] bit to split original  space into left space
  long int* ihfbit  //!<[out] half bit to split original  space
)
{
  int tmpNsite = Nsite;
  switch (iCalcModel) {
  case DC::HubbardGC:
  case DC::KondoGC:
  case DC::HubbardNConserved:
  case DC::Hubbard:
  case DC::Kondo:
    tmpNsite *= 2;
    break;
  case DC::Spin:
  case DC::SpinGC:
    break;
  default:
    fprintf(stderr, "Error: CalcModel %d is incorrect.\n", iCalcModel);
    return -1;
  }

  if (GetSplitBit(tmpNsite, irght, ilft, ihfbit) != 0) {
    return -1;
  }
  return 0;
}
/** 
@brief function of getting right, left and half bits corresponding to a original  space.
@retval 0 normally finished
@retval -1 unnormally finished
@version 0.2
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
int GetSplitBitForGeneralSpin(
  const int Nsite,  //!<[in] total number of sites
  long int* ihfbit, //!<[out] bit to split original  space
  const long int* SiteToBit  //!<[in]
) {
  int isite = 0;
  long int isqrtMaxDim = 1;
  long int tmpbit = 1;

  if (Nsite < 1) {
    fprintf(stderr, "%s", "Error: Total Site Number is incorrect.\n");
    return -1;
  }

  for (isite = 1; isite <= Nsite; isite++) {
    isqrtMaxDim *= SiteToBit[isite - 1];
  }
  isqrtMaxDim = (long int)sqrt(isqrtMaxDim);

  for (isite = 1; isite <= Nsite; isite++) {
    tmpbit *= SiteToBit[isite - 1];
    if (tmpbit >= isqrtMaxDim) break;
  }
  *ihfbit = tmpbit;
  return 0;
}
/** 
@brief function of splitting a original bit to right and left spaces
@version 0.1 
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
void SplitBit(
  const long int ibit,  //!<[in] original bit
  const long int irght,  //!<[in] bit to split original  space into right space
  const long int ilft,  //!<[in] bit to split original  space into left space
  const long int ihfbit,  //!<[in] half bit to split original  space
  long int* isplited_Bit_right,  //!<[out] splitted bit reflected on right space
  long int* isplited_Bit_left  //!<[out] splitted bit reflected on left space
)
{
  *isplited_Bit_right = ibit & irght;
  *isplited_Bit_left = ibit & ilft;
  *isplited_Bit_left = *isplited_Bit_left / ihfbit;
}
/** 
@brief function of getting off-diagonal component
@version 0.1 
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
int GetOffComp(
  long int* _list_2_1,  //!<[in] list to right space
  long int* _list_2_2, //!<[in] list to left space
  long int _ibit,  //!<[in] original bit 
  const long int _irght, //!<[in] bit to split original  space into right space
  const long int _ilft,  //!<[in] bit to split original  space into left space
  const long int _ihfbit,  //!<[in] half bit to split original  space
  long int* _ioffComp  //!<[out] off diagonal component
)
{
  long int ia, ib;
  SplitBit(_ibit, _irght, _ilft, _ihfbit, &ia, &ib);

  if (_list_2_1[ia] * _list_2_2[ib] == 0) {
    *_ioffComp = 0;
    return FALSE;
  }
  *_ioffComp = _list_2_1[ia] - 1;
  *_ioffComp += _list_2_2[ib] - 1;

  return TRUE;
}
/** 
@brief function of getting off-diagonal component for general spin
@retval FALSE off-diagonal component does not exist
@retval TRUE off-diagonal component exists
@version 0.2
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
int GetOffCompGeneralSpin(
  const long int org_ibit,  //!<[in] original bit
  const int org_isite,  //!<[in] target site 
  const int org_ispin, //!<[in] target spin to delete.
  const int off_ispin, //!<[in] target spin to create.
  long  int* _ioffComp,  //!<[out] generated bit 
  const long int* SiteToBit,  //!<[in]t List for getting bit at a site
  const long int* Tpow  //!<[in] List for getting total bit at a site before
)
{
  if (off_ispin > SiteToBit[org_isite] - 1 ||
    off_ispin<0 ||
    org_ispin>SiteToBit[org_isite] - 1 ||
    org_ispin < 0) {
    *_ioffComp = 0;
    return FALSE;
  }
  if (BitCheckGeneral(org_ibit, org_isite, org_ispin, SiteToBit, Tpow) == FALSE) {
    *_ioffComp = 0;
    return FALSE;
  }

  //delete org_ispin and create off_ispin
  long int tmp_off = 0;
  tmp_off = (long int)(off_ispin - org_ispin);
  tmp_off *= Tpow[org_isite];
  tmp_off += org_ibit;
  *_ioffComp = tmp_off;
  return TRUE;
}
/** 
@brief function of converting component to list_1
@version 0.2
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
int ConvertToList1GeneralSpin(
  const long int org_ibit,  //!<[in] original bit
  const long int ihlfbit,    //!<[in] split bit for general spin
  long int* _ilist1Comp,     //!<[out] component converted to list_1
  long int *list_2_1,
  long int *list_2_2
)
{
  long int ia, ib;
  ia = org_ibit % ihlfbit;
  ib = org_ibit / ihlfbit;
  if (list_2_1[ia] * list_2_2[ib] == 0) {
    *_ilist1Comp = 0;
    return FALSE;
  }
  *_ilist1Comp = list_2_1[ia] + list_2_2[ib] - 2;
  return TRUE;
}
/**
@brief function of getting fermion signs (for 32bit)
@version 0.1
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
void SgnBit_old(
  const long int org_bit,  //!<[in] an original bit
  int* sgn  //!<[out] fermion sign 
)
{
  long int bit;

  bit = org_bit ^ (org_bit >> 1);
  bit = (bit ^ (bit >> 2)) & 0x11111111;
  bit = bit * 0x11111111;
  *sgn = 1 - 2 * ((bit >> 28) & 1); // sgn = pm 1
}
// for 64 bit
/** 
@brief function of getting fermion sign (64 bit)
@version 0.1
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
void SgnBit(
  const long int org_bit,  //!<[in] an original bit
  int* sgn  //!<[out] sgn fermion sign 
)
{
  long int bit;

  bit = org_bit ^ (org_bit >> 1);
  bit = bit ^ (bit >> 2);
  bit = bit ^ (bit >> 4);
  bit = bit ^ (bit >> 8);
  bit = bit ^ (bit >> 16);
  bit = bit ^ (bit >> 32);
  *sgn = 1 - 2 * (bit & 1); // sgn = pm 1
}
/** 
@brief bit check function
@retval 1
@retval 0 
@version 0.1
@author Takahiro Misawa (The University of Tokyo) 
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
int BitCheck(
  const long int org_bit,  //!<[in] original bit to check
  const long int target_bit  //!<[in] target bit to check
)
{
  return  (org_bit >> target_bit) & 1;
  // (org_bit & (2^target_bit))/2^target_bit
}
/** 
 * 
 * @brief bit check function for general spin
 * @retval 0 bit does not exists
 * @retval 1 bit exists
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int BitCheckGeneral(
  const long int org_bit,  //!<[in] original bit to check
  const int org_isite,  //!<[in] site index (org_isite >= 1)
  const int target_ispin,  //!<[in] target spin to check 
  const long int* SiteToBit,  //!<[in] List for getting bit at a site
  const long int* Tpow  //!<[in] List for getting total bit at a site before
)
{
  if (GetBitGeneral(org_isite, org_bit, SiteToBit, Tpow) != target_ispin) {
    return FALSE;
  }
  return TRUE;
}
/**
@brief get bit at a site for general spin
@return bit at a site
@version 0.2
@author Kazuyoshi Yoshimi (The University of Tokyo) 
*/
int GetBitGeneral(
  const int isite,  //!<[in] site index (isite >= 1)
  const long int org_bit,  //!<[in] original bit to check 
  const long int* SiteToBit,  //!<[in] List for getting bit at a site
  const long int* Tpow  //!<[in] List for getting total bit at a site before
)
{
  long int tmp_bit = (org_bit / Tpow[isite]) % SiteToBit[isite];
  return (tmp_bit);
}
/**
 @brief get 2sz at a site for general spin
 @return 2sz at isite
 @version 0.2
 @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetLocal2Sz
(
  const int isite,  //!<[in] site index (isite >= 1)
  const long int org_bit,  //!<[in] original bit to check 
  const long int* SiteToBit,  //!<[in] List for getting bit at a site
  const long int* Tpow  //!<[in] List for getting total bit at a site before
)
{
  int TwiceSz = 0;
  int bitAtSite = 0;
  //get bit
  bitAtSite = GetBitGeneral(isite, org_bit, SiteToBit, Tpow);
  TwiceSz = -(SiteToBit[isite] - 1) + 2 * bitAtSite; //-2S^{total}_i+2Sz_i
  return TwiceSz;
}
/** 
@brief "finding the next higher number after a given number that has the same number of 1-bits" 
This method is introduced in S.H. Warren, Hacker$B!G(Bs Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
@param x 
@version 2.0
@author Takahiro Misawa (The University of Tokyo) 
*/
long int snoob(long int x) {
  long int smallest, ripple, ones;
  smallest = x & (-x);
  ripple = x + smallest;
  ones = x ^ ripple;
  ones = (ones >> 2) / smallest;
  return   ripple | ones;
}
/** 
@brief calculating number of 1-bits in x (32 bit)
This method is introduced in S.H. Warren, Hacker$B!G(Bs Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
@version 2.0
@author Takahiro Misawa (The University of Tokyo) 
*/
int pop(int x) {
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x + (x >> 4)) & 0x0F0F0F0F;
  x = x + (x >> 8);
  x = x + (x >> 16);
  return  x & 0x0000003F;
}
