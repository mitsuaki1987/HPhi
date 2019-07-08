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

int GetSplitBit(const int Nsite, long int* irght, long int* ilft, long int* ihfbit);
int GetSplitBitByModel(const int Nsite, const int iCalcModel, long int* irght,
  long int* ilft, long int* ihfbit);
int GetSplitBitForGeneralSpin(const int Nsite, long int* ihfbit, const long int* SiteToBit);
void SplitBit(const long int ibit, const long int irght, const long int ilft,
  const long int ihfbit, long int* isplited_Bit_right, long int* isplited_Bit_left);
int GetOffComp(long int* _list_2_1, long int* _list_2_2, long int _ibit,
  const long int _irght, const long int _ilft, const long int _ihfbit, long int* _ioffComp);
void SgnBit_old(const long int bit, int* sgn);
void SgnBit(const long int bit, int* sgn);
int BitCheck(const long int org_bit, const long int target_bit);
int BitCheckGeneral(const long int org_bit, const int org_isite, const int target_ispin,
  const long int* SiteToBit, const long int* TPow);
int GetBitGeneral(const int isite, const long int org_bit,
  const long int* SiteToBit, const long int* TPow);
int GetOffCompGeneralSpin(const long int org_ibit, const int org_isite, const int org_ispin,
  const int off_ispin, long  int* _ioffComp, const long int* SiteToBit, const long int* TPow);
int GetLocal2Sz(const int isite, const long int org_bit, 
  const long int* SiteToBit, const long int* Tpow);
int ConvertToList1GeneralSpin(const long int org_ibit,
  const long int ihlfbit, long int* _ilist1Comp);
long int snoob(long int x);
int pop(int x);
