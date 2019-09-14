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
#include "common/setmemory.hpp"
#include "FileIO.hpp"
#include "sz.hpp"
#include "wrapperMPI.hpp"
#include "xsetmem.hpp"
#include "log.hpp"
#include "global.hpp"
#include "DefCommon.hpp"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @file
 * 
 * @brief  Generating Hilbert spaces
 * 
 * @version 0.2
 * @details 
 *
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */

 /**
  * @brief calculating restricted Hilbert space for Kondo-GC systems
  *
  * @param[in] ib   upper half bit of i
  * @param[in] ihfbit 2^(Ns/2)
  * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia)
  * @param[out] list_2_1_  list_2_1_[ib] = jb
  * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
  * @param[in] list_jb_   list_jb_[ib]  = jb
  *
  * @return
  * @author Takahiro Misawa (The University of Tokyo)
  */
int child_omp_sz_Kondo_hacker(
  long int ib,
  long int ihfbit,
  long int *list_1_,
  long int *list_2_1_,
  long int *list_2_2_,
  long int *list_jb_
)
{
  long int j;
  unsigned long int i;
  unsigned long int ia, ja, jb;
  unsigned long int div_down, div_up;
  unsigned long int num_up, num_down;
  unsigned long int tmp_num_up, tmp_num_down;
  unsigned int icheck_loc;

  jb = list_jb_[ib];
  i = ib * ihfbit;

  num_up = 0;
  num_down = 0;
  icheck_loc = 1;
  for (j = Def::Nsite / 2; j < Def::Nsite; j++) {
    div_up = i & Def::Tpow[2 * j];
    div_up = div_up / Def::Tpow[2 * j];
    div_down = i & Def::Tpow[2 * j + 1];
    div_down = div_down / Def::Tpow[2 * j + 1];

    if (Def::LocSpn[j] == ITINERANT) {
      num_up += div_up;
      num_down += div_down;
    }
    else {
      num_up += div_up;
      num_down += div_down;
      if (Def::Nsite % 2 == 1 && j == (Def::Nsite / 2)) {
        icheck_loc = icheck_loc;
      }
      else {
        icheck_loc = icheck_loc * (div_up^div_down);// exclude doubly occupied site
      }
    }
  }
  //[s] get ja  
  ja = 0;
  tmp_num_up = num_up;
  tmp_num_down = num_down;
  if (icheck_loc == 1) {
    //for(ia=0;ia<Check::sdim;ia++){
    ia = Def::Tpow[Def::Nup + Def::Ndown - tmp_num_up - tmp_num_down] - 1;
    //ia = 1;
    //if(ia < Check::sdim && ia!=0){
    //ia = snoob(ia);
    while (ia < (unsigned long)Check::sdim && ia != 0) {
      // for(ia=0;ia<Check::sdim;ia++){
          //[s] proceed ja
      i = ia;
      num_up = tmp_num_up;
      num_down = tmp_num_down;
      icheck_loc = 1;
      for (j = 0; j < (Def::Nsite + 1) / 2; j++) {
        div_up = i & Def::Tpow[2 * j];
        div_up = div_up / Def::Tpow[2 * j];
        div_down = i & Def::Tpow[2 * j + 1];
        div_down = div_down / Def::Tpow[2 * j + 1];

        if (Def::LocSpn[j] == ITINERANT) {
          num_up += div_up;
          num_down += div_down;
        }
        else {
          num_up += div_up;
          num_down += div_down;
          if (Def::Nsite % 2 == 1 && j == (Def::Nsite / 2)) {
            icheck_loc = icheck_loc;
          }
          else {
            icheck_loc = icheck_loc * (div_up^div_down);// exclude doubllly ocupited site
          }
        }
      }

      if (icheck_loc == 1 && Def::LocSpn[Def::Nsite / 2] != ITINERANT && Def::Nsite % 2 == 1) {
        div_up = ia & Def::Tpow[Def::Nsite - 1];
        div_up = div_up / Def::Tpow[Def::Nsite - 1];
        div_down = (ib*ihfbit) & Def::Tpow[Def::Nsite];
        div_down = div_down / Def::Tpow[Def::Nsite];
        icheck_loc = icheck_loc * (div_up^div_down);
      }

      if (num_up == (unsigned long)Def::Nup 
        && num_down == (unsigned long)Def::Ndown && icheck_loc == 1) {
        //printf("ia=%ud ja=%ud \n",ia,ja);
        list_1_[ja + jb] = ia + ib * ihfbit;
        list_2_1_[ia] = ja + 1;
        list_2_2_[ib] = jb + 1;
        ja += 1;
      }
      ia = snoob(ia);
      //[e] proceed ja
      //ia+=1;
    //}
    }
  }
  //[e] get ja
  return ja;
}
/** 
 * 
 * @file
 * @brief  calculating binomial coefficients
 * 
 * @param[in] n      n for @f$_nC_k = \frac{n!}{(n-k)!k!}@f$
 * @param[in] k      k for @f$_nC_k = \frac{n!}{(n-k)!k!}@f$
 * @param[out] comb   binomial coefficients @f$_nC_k@f$
 * @param[in] Nsite  # of sites
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
long int Binomial(int n,int k,long int **comb,int Nsite){
  // nCk, Nsite=max(n)
  int tmp_i,tmp_j;

  if(n==0 && k==0){
    return 1;
  } 
  else if(n<0 || k<0 || n<k){
    return 0;
  }
  
  for(tmp_i=0;tmp_i<=Nsite;tmp_i++){
    for(tmp_j=0;tmp_j<=Nsite;tmp_j++){
      comb[tmp_i][tmp_j] = 0;
    }
  }

  comb[0][0] = 1;
  comb[1][0] = 1;
  comb[1][1] = 1;
  for(tmp_i=2;tmp_i<=n;tmp_i++){
    for(tmp_j=0;tmp_j<=tmp_i;tmp_j++){
      if(tmp_j==0){
        comb[tmp_i][tmp_j] = 1;
      }else if(tmp_j==tmp_i){
        comb[tmp_i][tmp_j] = 1;
      }else{
        comb[tmp_i][tmp_j] = comb[tmp_i-1][tmp_j-1]+comb[tmp_i-1][tmp_j];
      }
    }
  }
  return comb[n][k];
}

/** 
 * @brief calculating restricted Hilbert space for Hubbard systems
 * 
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_omp_sz(
  long int ib,    //!<[in]
  long int ihfbit, //!<[in]
  long int* list_1_, //!<[out]
  long int* list_2_1_,//!<[out]
  long int* list_2_2_,//!<[out]
  long int* list_jb_ //!<[in]
)
{
  long int i,j; 
  long int ia,ja,jb;
  long int div_down, div_up;
  long int num_up,num_down;
  long int tmp_num_up,tmp_num_down;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  for(j=0;j< Def::Nsite ;j++){
    div_up    = i & Def::Tpow[2*j];
    div_up    = div_up/Def::Tpow[2*j];
    div_down  = i & Def::Tpow[2*j+1];
    div_down  = div_down/Def::Tpow[2*j+1];
    num_up += div_up;
    num_down += div_down;
  }
  
  ja=0;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;

  if(Def::iCalcModel== DC::Hubbard){
    for(ia=0;ia<Check::sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      
      for(j=0;j<Def::Nsite;j++){
        div_up    = i & Def::Tpow[2*j];
        div_up    = div_up/Def::Tpow[2*j];
        div_down  = i & Def::Tpow[2*j+1];
        div_down  = div_down/Def::Tpow[2*j+1];
        num_up += div_up;
        num_down += div_down;
      }
      if(num_up == Def::Nup && num_down == Def::Ndown){
        list_1_[ja+jb]=ia+ib*ihfbit;
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        ja+=1;
      } 
    }
  }
  else if(Def::iCalcModel== DC::HubbardNConserved){
    for(ia=0;ia<Check::sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      
      for(j=0;j<Def::Nsite;j++){
        div_up    = i & Def::Tpow[2*j];
        div_up    = div_up/Def::Tpow[2*j];
        div_down  = i & Def::Tpow[2*j+1];
        div_down  = div_down/Def::Tpow[2*j+1];
        num_up += div_up;
        num_down += div_down;
      }
      if( (num_up+num_down) == Def::Ne){
        list_1_[ja+jb]=ia+ib*ihfbit;
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        ja+=1;
      } 
    }  
  }
  return ja; 
}
/** 
 * @brief efficient version of calculating restricted Hilbert space for Hubbard systems  using snoob
 * details of snoob is found in S.H. Warren, Hacker’s Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int child_omp_sz_hacker(long int ib,
  long int ihfbit,
  long int *list_1_,
  long int *list_2_1_,
  long int *list_2_2_,
  long int *list_jb_
)
{
  long int j;
  unsigned long int i;
  unsigned long int ia, ja, jb;
  unsigned long int div_down, div_up;
  unsigned long int num_up, num_down;
  unsigned long int tmp_num_up, tmp_num_down;

  jb = list_jb_[ib];
  i = ib * ihfbit;

  num_up = 0;
  num_down = 0;
  for (j = 0; j < Def::Nsite; j++) {
    div_up = i & Def::Tpow[2 * j];
    div_up = div_up / Def::Tpow[2 * j];
    div_down = i & Def::Tpow[2 * j + 1];
    div_down = div_down / Def::Tpow[2 * j + 1];
    num_up += div_up;
    num_down += div_down;
  }

  ja = 0;
  tmp_num_up = num_up;
  tmp_num_down = num_down;

  if (Def::iCalcModel == DC::Hubbard) {
    if (tmp_num_up <= (unsigned long)Def::Nup
      && tmp_num_down <= (unsigned long)Def::Ndown) { //do not exceed Nup and Ndown
      ia = Def::Tpow[Def::Nup + Def::Ndown - tmp_num_up - tmp_num_down] - 1;
      if (ia < (unsigned long)Check::sdim) {
        num_up = tmp_num_up;
        num_down = tmp_num_down;
        for (j = 0; j < Def::Nsite; j++) {
          div_up = ia & Def::Tpow[2 * j];
          div_up = div_up / Def::Tpow[2 * j];
          div_down = ia & Def::Tpow[2 * j + 1];
          div_down = div_down / Def::Tpow[2 * j + 1];
          num_up += div_up;
          num_down += div_down;
        }
        if (num_up == (unsigned long)Def::Nup && num_down == (unsigned long)Def::Ndown) {
          list_1_[ja + jb] = ia + ib * ihfbit;
          list_2_1_[ia] = ja + 1;
          list_2_2_[ib] = jb + 1;
          ja += 1;
        }
        if (ia != 0) {
          ia = snoob(ia);
          while (ia < (unsigned long)Check::sdim) {
            num_up = tmp_num_up;
            num_down = tmp_num_down;
            for (j = 0; j < Def::Nsite; j++) {
              div_up = ia & Def::Tpow[2 * j];
              div_up = div_up / Def::Tpow[2 * j];
              div_down = ia & Def::Tpow[2 * j + 1];
              div_down = div_down / Def::Tpow[2 * j + 1];
              num_up += div_up;
              num_down += div_down;
            }
            if (num_up == (unsigned long)Def::Nup && num_down == (unsigned long)Def::Ndown) {
              list_1_[ja + jb] = ia + ib * ihfbit;
              list_2_1_[ia] = ja + 1;
              list_2_2_[ib] = jb + 1;
              ja += 1;
            }
            ia = snoob(ia);
          }
        }
      }
    }
  }
  else if (Def::iCalcModel == DC::HubbardNConserved) {
    if (tmp_num_up + tmp_num_down <= (unsigned long)Def::Ne) { //do not exceed Ne
      ia = Def::Tpow[Def::Ne - tmp_num_up - tmp_num_down] - 1;
      if (ia < (unsigned long)Check::sdim) {
        list_1_[ja + jb] = ia + ib * ihfbit;
        list_2_1_[ia] = ja + 1;
        list_2_2_[ib] = jb + 1;
        ja += 1;
        if (ia != 0) {
          ia = snoob(ia);
          while (ia < (unsigned long)Check::sdim) {
            list_1_[ja + jb] = ia + ib * ihfbit;
            list_2_1_[ia] = ja + 1;
            list_2_2_[ib] = jb + 1;
            ja += 1;
            ia = snoob(ia);
          }
        }
      }
    }
  }
  return ja;
}
/** 
 * @brief calculating restricted Hilbert space for Kondo systems
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int child_omp_sz_Kondo(
  long int ib,        //[in]
  long int ihfbit,    //[in]
  long int* list_1_,  //[out]
  long int* list_2_1_,//[out]
  long int* list_2_2_,//[out]
  long int* list_jb_  //[in]
)
{
  long int i,j; 
  long int ia,ja,jb;
  long int div_down, div_up;
  long int num_up,num_down;
  long int tmp_num_up,tmp_num_down;
  int icheck_loc;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  icheck_loc=1;
  for(j=Def::Nsite/2; j< Def::Nsite ;j++){
    div_up    = i & Def::Tpow[2*j];
    div_up    = div_up/Def::Tpow[2*j];
    div_down  = i & Def::Tpow[2*j+1];
    div_down  = div_down/Def::Tpow[2*j+1];

    if(Def::LocSpn[j] == ITINERANT){
      num_up   += div_up;        
      num_down += div_down;  
    }else{    
      num_up   += div_up;        
      num_down += div_down;
      if(Def::Nsite%2==1 && j==(Def::Nsite/2)){
        icheck_loc= icheck_loc;
      }
      else{
        icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
      }
    }
  }
  
  ja = 0;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;
  if(icheck_loc ==1){
    for(ia=0;ia<Check::sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      icheck_loc=1;
      for(j=0;j<(Def::Nsite+1)/2;j++){
        div_up    = i & Def::Tpow[2*j];
        div_up    = div_up/Def::Tpow[2*j];
        div_down  = i & Def::Tpow[2*j+1];
        div_down  = div_down/Def::Tpow[2*j+1];

        if(Def::LocSpn[j] ==  ITINERANT){
          num_up   += div_up;        
          num_down += div_down;  
        }else{    
          num_up   += div_up;        
          num_down += div_down;  
          if(Def::Nsite%2==1 && j==(Def::Nsite/2)){
            icheck_loc= icheck_loc;
          }
          else{
            icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
          }
        }
      }
      
      if(icheck_loc == 1 && Def::LocSpn[Def::Nsite/2] != ITINERANT && Def::Nsite%2==1){
        div_up    = ia & Def::Tpow[Def::Nsite-1];
        div_up    = div_up/Def::Tpow[Def::Nsite-1];
        div_down  = (ib*ihfbit) & Def::Tpow[Def::Nsite];
        div_down  = div_down/Def::Tpow[Def::Nsite];
        icheck_loc= icheck_loc*(div_up^div_down);
      }
      
      if(num_up == Def::Nup && num_down == Def::Ndown && icheck_loc==1){
        list_1_[ja+jb]=ia+ib*ihfbit;
        /*
        list_2_1_[ia]=ja;
        list_2_2_[ib]=jb;
         */
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        //printf("DEBUG: rank=%d, list_1[%d]=%d, list_2_1_[%d]=%d, list_2_2_[%d]=%d\n", MP::myrank, ja+jb, list_1_[ja+jb], ia, list_2_1[ia], ib, list_2_2[ib]);
        ja+=1;
      }
    }
  }
  return ja; 
}
/** 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_omp_sz_KondoGC(
  long int ib,  //!<[in]
  long int ihfbit,//!<[in]
  long int* list_1_, //!<[out]
  long int* list_2_1_,//!<[out]
  long int* list_2_2_,//!<[out]
  long int* list_jb_//!<[in]
)
{
  long int i,j; 
  long int ia,ja,jb;
  long int div_down, div_up;
  int icheck_loc;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
  icheck_loc=1;
  for(j=Def::Nsite/2; j< Def::Nsite ;j++){
    div_up    = i & Def::Tpow[2*j];
    div_up    = div_up/Def::Tpow[2*j];
    div_down  = i & Def::Tpow[2*j+1];
    div_down  = div_down/Def::Tpow[2*j+1];
    if(Def::LocSpn[j] !=  ITINERANT){
      if(Def::Nsite%2==1 && j==(Def::Nsite/2)){
        icheck_loc= icheck_loc;
      }
      else{
        icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
      }
    }
  }

  ja=0;
  if(icheck_loc ==1){
    for(ia=0;ia<Check::sdim;ia++){
      i=ia;
      icheck_loc =1;
      for(j=0;j<(Def::Nsite+1)/2;j++){
        div_up    = i & Def::Tpow[2*j];
        div_up    = div_up/Def::Tpow[2*j];
        div_down  = i & Def::Tpow[2*j+1];
        div_down  = div_down/Def::Tpow[2*j+1];
        if(Def::LocSpn[j] !=  ITINERANT){
          if(Def::Nsite%2==1 && j==(Def::Nsite/2)){
            icheck_loc= icheck_loc;
          }
          else{
            icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
          }
        }
      }

      if(icheck_loc == 1 && Def::LocSpn[Def::Nsite/2] != ITINERANT && Def::Nsite%2==1){
        div_up    = ia & Def::Tpow[Def::Nsite-1];
        div_up    = div_up/Def::Tpow[Def::Nsite-1];
        div_down  = (ib*ihfbit) & Def::Tpow[Def::Nsite];
        div_down  = div_down/Def::Tpow[Def::Nsite];
        icheck_loc= icheck_loc*(div_up^div_down);
      }
      
      if(icheck_loc==1){
        list_1_[ja+jb]=ia+ib*ihfbit;
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        ja+=1;
      }
    }
  }
  return ja; 
}

/** 
 * @brief calculating restricted Hilbert space for spin-1/2 systems
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] N ???
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia)
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int child_omp_sz_spin(
  long int ib,
  long int ihfbit,
  int N,
  long int *list_1_,
  long int *list_2_1_,
  long int *list_2_2_,
  long int *list_jb_
)
{
  long int i, j, div;
  long int ia, ja, jb;
  long int num_up;
  int tmp_num_up;

  jb = list_jb_[ib];
  i = ib * ihfbit;
  num_up = 0;
  for (j = 0; j < N; j++) {
    div = i & Def::Tpow[j];
    div = div / Def::Tpow[j];
    num_up += div;
  }
  ja = 0;
  tmp_num_up = num_up;

  for (ia = 0; ia < ihfbit; ia++) {
    i = ia;
    num_up = tmp_num_up;
    for (j = 0; j < N; j++) {
      div = i & Def::Tpow[j];
      div = div / Def::Tpow[j];
      num_up += div;
    }

    if (num_up == Def::Ne) {
      list_1_[ja + jb] = ia + ib * ihfbit;
      list_2_1_[ia] = ja + 1;
      list_2_2_[ib] = jb + 1;
      ja += 1;
    }
  }
  return ja;
}
/** 
 * @brief efficient version of calculating restricted Hilbert space for spin-1/2 systems 
 * details of snoob is found in S.H. Warren, Hacker’s Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] N ???
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia)
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int child_omp_sz_spin_hacker(
  long int ib,
  long int ihfbit,
  int N,
  long int *list_1_,
  long int *list_2_1_,
  long int *list_2_2_,
  long int *list_jb_
)
{
  long int j;
  unsigned long int i, div;
  unsigned long int ia, ja, jb;
  unsigned long int num_up;
  unsigned int tmp_num_up;

  jb = list_jb_[ib];
  i = ib * ihfbit;
  num_up = 0;
  for (j = 0; j < N; j++) {
    div = i & Def::Tpow[j];
    div = div / Def::Tpow[j];
    num_up += div;
  }
  ja = 0;
  tmp_num_up = num_up;

  // using hacker's delight
  if (tmp_num_up <= (unsigned long)Def::Ne
    && ((unsigned long)Def::Ne - tmp_num_up) <= (unsigned long)Def::Nsite - 1) { // do not exceed Ne
    ia = Def::Tpow[Def::Ne - tmp_num_up] - 1;
    if (ia < (unsigned long)ihfbit) {          // do not exceed Ne
      list_1_[ja + jb] = ia + ib * ihfbit;
      list_2_1_[ia] = ja + 1;
      list_2_2_[ib] = jb + 1;
      ja += 1;

      if (ia != 0) {
        ia = snoob(ia);
        while (ia < (unsigned long)ihfbit) {
          //fprintf(MP::STDOUT, " X: ia= %ld ia=%ld \n", ia,ia);
          list_1_[ja + jb] = ia + ib * ihfbit;
          list_2_1_[ia] = ja + 1;
          list_2_2_[ib] = jb + 1;
          ja += 1;
          ia = snoob(ia);
        }
      }
    }
  }
  return ja;
}
/** 
 * @brief calculating restricted Hilbert space for general spin systems (S>1/2)
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[out] list_2_1_Sz_  
 * @param[out] list_2_2_Sz_  
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int child_omp_sz_GeneralSpin(
  long int ib,
  long int ihfbit,
  long int* list_1_,
  long int* list_2_1_,
  long int* list_2_2_,
  long int* list_2_1_Sz_,
  long int* list_2_2_Sz_,
  long int* list_jb_
)
{
  long int ia,ja,jb;  
  int list_2_2_Sz_ib=0;
  int tmp_2Sz=0;
  jb = list_jb_[ib];
  list_2_2_Sz_ib =list_2_2_Sz_[ib];
  ja = 0;
  for(ia=0;ia<ihfbit;ia++){
    tmp_2Sz=list_2_1_Sz_[ia]+list_2_2_Sz_ib;
    if(tmp_2Sz == Def::Total2Sz){
      list_1_[ja+jb]=ia+ib*ihfbit;
      list_2_1_[ia]=ja+1;
      list_2_2_[ib]=jb+1;
      ja+=1;
    } 
  }
  return ja; 
}

/** 
 * @brief reading the list of the restricted Hilbert space
 * 
 * @param[in] irght 
 * @param[in] ilft 
 * @param[in] ihfbit 
 * @param[in] i_max 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int Read_sz
(
 const long int irght,
 const long int ilft,
 const long int ihfbit,
 long int *i_max
 )
{
  FILE *fp,*fp_err;
  char sdt[D_FileNameMax];
  char buf[D_FileNameMax];
    
  long int icnt=0; 
  long int ia,ib;
  long int ja=0;
  long int jb=0;
  long int ibpatn=0;
  long int dam; 

  TimeKeeper("%s_sz_TimeKeeper.dat","READ = 1: read starts: %s", "a");
  TimeKeeper("%s_TimeKeeper.dat","READ = 1: read starts: %s", "a");

  switch(Def::iCalcModel){
  case DC::Hubbard:
  case DC::HubbardGC:
  case DC::Spin:
  case DC::SpinGC:
    sprintf(sdt,"ListForModel_Ns%d_Nup%dNdown%d.dat", Def::Nsite, Def::Nup, Def::Ndown);
    break;
  case DC::Kondo:
    sprintf(sdt,"ListForKondo_Ns%d_Ncond%d.dat",Def::Nsite,Def::Ne);
    break;
  }
  if(childfopenMPI(sdt,"r", &fp)!=0){
    wrapperMPI::Exit(-1);
  }  

  if(fp == NULL){
    if(childfopenMPI("Err_sz.dat","a",&fp_err)!=0){
      wrapperMPI::Exit(-1);
    }
    fprintf(fp_err, "%s", "No file. Please set READ=0.\n");
    fprintf(stderr, "%s", "No file. Please set READ=0.\n");
    fprintf(fp_err, " %s does not exist. \n",sdt);
    fprintf(stderr, " %s does not exist. \n", sdt);
    fclose(fp_err);
  }else{
    while(NULL != wrapperMPI::Fgets(buf,sizeof(buf),fp)){  
      dam=atol(buf);  
      List::c1[icnt]=dam;
            
      ia= dam & irght;
      ib= dam & ilft;
      ib=ib/ihfbit; 
            
      if(ib==ibpatn){
        ja=ja+1;
      }else{
        ibpatn=ib;
        ja=1;
        jb=icnt-1;
      }
            
      List::c2_1[ia]=ja;
      List::c2_2[ib]=jb;
      icnt+=1;
                
    }
    fclose(fp);
    *i_max=icnt-1;
  }

  TimeKeeper("%s_sz_TimeKeeper.dat", "READ = 1: read finishes: %s", "a");
  TimeKeeper("%s_TimeKeeper.dat", "READ = 1: read finishes: %s", "a");

  return 0;
}
/**
 *
 * @brief generating Hilbert space
 *
 * @param[out] list_1_   list_1[icnt] = i (index of full Hilbert space) : icnt = index in the restricted Hilbert space
 * @param[out] list_2_1_ icnt=list_2_1[]+list_2_2[]
 * @param[out] list_2_2_
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int sz
(
  long int *list_1_,
  long int *list_2_1_,
  long int *list_2_2_
)
{
  FILE *fp, *fp_err;
  char sdt[D_FileNameMax], sdt_err[D_FileNameMax];
  long int *HilbertNumToSz;
  long int i, icnt;
  long int ib, jb;

  long int j;
  long int div;
  long int num_up, num_down;
  long int irght, ilft, ihfbit;

  //*[s] for omp parall
  int  all_up, all_down, tmp_res, num_threads;
  long int tmp_1, tmp_2, tmp_3;
  long int **comb;
  //*[e] for omp parall

  // [s] for Kondo
  int N_all_up, N_all_down;
  int all_loc;
  long int num_loc, div_down;
  int num_loc_up;
  int icheck_loc;
  int ihfSpinDown = 0;
  // [e] for Kondo

  long int i_max = 0;
  double idim = 0.0;
  long int div_up;

  // [s] for general spin
  long int *list_2_1_Sz;
  long int *list_2_2_Sz;
  if (Def::iFlgGeneralSpin == TRUE) {
    list_2_1_Sz = li_1d_allocate(Check::sdim + 2);
    list_2_2_Sz = li_1d_allocate((Def::Tpow[Def::Nsite - 1] * Def::SiteToBit[Def::Nsite - 1] / Check::sdim) + 2);
    for (j = 0; j < Check::sdim + 2; j++) {
      list_2_1_Sz[j] = 0;
    }
    for (j = 0; j < (Def::Tpow[Def::Nsite - 1] * Def::SiteToBit[Def::Nsite - 1] / Check::sdim) + 2; j++) {
      list_2_2_Sz[j] = 0;
    }
  }
  // [e] for general spin

  long int *list_jb;
  list_jb = li_1d_allocate(Large::SizeOflistjb);
  for (i = 0; i < Large::SizeOflistjb; i++) {
    list_jb[i] = 0;
  }

  //hacker
  int hacker;
  long int tmp_i, tmp_j, tmp_pow, max_tmp_i;
  long int ia, ja;
  long int ibpatn = 0;
  //hacker

  int iSpnup, iMinup, iAllup;
  int N2 = 0;
  int N = 0;
  fprintf(MP::STDOUT, "%s", "  Start: Calculate HilbertNum for fixed Sz. \n");
  TimeKeeper("%s_sz_TimeKeeper.dat", "initial sz : %s", "w");
  TimeKeeper("%s_TimeKeeper.dat", "initial sz : %s", "a");

  if (Check::idim_max != 0) {
    switch (Def::iCalcModel) {
    case DC::HubbardGC:
    case DC::HubbardNConserved:
    case DC::Hubbard:
      N2 = 2 * Def::Nsite;
      idim = pow(2.0, N2);
      break;
    case DC::KondoGC:
    case DC::Kondo:
      N2 = 2 * Def::Nsite;
      N = Def::Nsite;
      idim = pow(2.0, N2);
      for (j = 0; j < N; j++) {
        fprintf(MP::STDOUT, "  j  =  %ld loc %d \n", j, Def::LocSpn[j]);
      }
      break;
    case DC::SpinGC:
    case DC::Spin:
      N = Def::Nsite;
      if (Def::iFlgGeneralSpin == FALSE) {
        idim = pow(2.0, N);
      }
      else {
        idim = 1;
        for (j = 0; j < N; j++) {
          idim *= Def::SiteToBit[j];
        }
      }
      break;
    }
    comb = li_2d_allocate(Def::Nsite + 1, Def::Nsite + 1);
    i_max = Check::idim_max;

    switch (Def::iCalcModel) {
    case DC::HubbardNConserved:
    case DC::Hubbard:
    case DC::KondoGC:
    case DC::Kondo:
    case DC::Spin:
      if (Def::iFlgGeneralSpin == FALSE) {
        if (GetSplitBitByModel(Def::Nsite, Def::iCalcModel, &irght, &ilft, &ihfbit) != 0) {
          wrapperMPI::Exit(-1);
        }
        Large::irght = irght;
        Large::ilft = ilft;
        Large::ihfbit = ihfbit;
        //fprintf(MP::STDOUT, "idim=%lf irght=%ld ilft=%ld ihfbit=%ld \n",idim,irght,ilft,ihfbit);
      }
      else {
        ihfbit = Check::sdim;
        //fprintf(MP::STDOUT, "idim=%lf ihfbit=%ld \n",idim, ihfbit);
      }
      break;
    default:
      break;
    }

    icnt = 1;
    jb = 0;

    if (Def::READ == 1) {
      if (Read_sz(irght, ilft, ihfbit, &i_max) != 0) {
        wrapperMPI::Exit(-1);
      }
    }
    else {
      sprintf(sdt, "%s_sz_TimeKeeper.dat", Def::CDataFileHead);
#ifdef _OPENMP
      num_threads = omp_get_max_threads();
#else
      num_threads = 1;
#endif
      childfopenMPI(sdt, "a", &fp);
      fprintf(fp, "num_threads==%d\n", num_threads);
      fclose(fp);

      //*[s] omp parallel

      TimeKeeper("%s_sz_TimeKeeper.dat", "omp parallel sz starts: %s", "a");
      TimeKeeper("%s_TimeKeeper.dat", "omp parallel sz starts: %s", "a");
      switch (Def::iCalcModel) {
      case DC::HubbardGC:
        icnt = Def::Tpow[2 * Def::Nsite - 1] * 2 + 0;/*Tpow[2*Def::Nsit]=1*/
        break;

      case DC::SpinGC:
        if (Def::iFlgGeneralSpin == FALSE) {
          icnt = Def::Tpow[Def::Nsite - 1] * 2 + 0;/*Tpow[Def::Nsit]=1*/
        }
        else {
          icnt = Def::Tpow[Def::Nsite - 1] * Def::SiteToBit[Def::Nsite - 1];
        }
        break;

      case DC::KondoGC:
        // this part can not be parallelized
        jb = 0;
        num_loc = 0;
        for (j = Def::Nsite / 2; j < Def::Nsite; j++) { // counting # of localized spins
          if (Def::LocSpn[j] != ITINERANT) { // //ITINERANT ==0 -> itinerant
            num_loc += 1;
          }
        }

        for (ib = 0; ib < Check::sdim; ib++) {
          list_jb[ib] = jb;
          i = ib * ihfbit;
          icheck_loc = 1;
          for (j = (Def::Nsite + 1) / 2; j < Def::Nsite; j++) {
            div_up = i & Def::Tpow[2 * j];
            div_up = div_up / Def::Tpow[2 * j];
            div_down = i & Def::Tpow[2 * j + 1];
            div_down = div_down / Def::Tpow[2 * j + 1];
            if (Def::LocSpn[j] != ITINERANT) {
              if (Def::Nsite % 2 == 1 && j == (Def::Nsite / 2)) {
                icheck_loc = icheck_loc;
              }
              else {
                icheck_loc = icheck_loc * (div_up^div_down);// exclude doubllly ocupited site
              }
            }
          }
          if (icheck_loc == 1) {
            if (Def::Nsite % 2 == 1 && Def::LocSpn[Def::Nsite / 2] != ITINERANT) {
              jb += Def::Tpow[Def::Nsite - 1 - (Def::NLocSpn - num_loc)];
            }
            else {
              jb += Def::Tpow[Def::Nsite - (Def::NLocSpn - num_loc)];
            }
          }
        }

        icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_jb, ihfbit, N2, Check::sdim)
        for (ib = 0; ib < Check::sdim; ib++) {
          icnt += child_omp_sz_KondoGC(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_jb);
        }
        break;

      case DC::Hubbard:
        hacker = Def::read_hacker;
        if (hacker == 0) {
          // this part can not be parallelized
          jb = 0;
          for (ib = 0; ib < Check::sdim; ib++) { // sdim = 2^(N/2)
            list_jb[ib] = jb;
            i = ib * ihfbit;
            //[s] counting # of up and down electrons
            num_up = 0;
            for (j = 0; j <= N2 - 2; j += 2) { // even -> up spin
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_up += div;
            }
            num_down = 0;
            for (j = 1; j <= N2 - 1; j += 2) { // odd -> down spin
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_down += div;
            }
            //[e] counting # of up and down electrons
            tmp_res = Def::Nsite % 2; // even Ns-> 0, odd Ns -> 1
            all_up = (Def::Nsite + tmp_res) / 2;
            all_down = (Def::Nsite - tmp_res) / 2;

            tmp_1 = Binomial(all_up, Def::Nup - num_up, comb, all_up);
            tmp_2 = Binomial(all_down, Def::Ndown - num_down, comb, all_down);
            jb += tmp_1 * tmp_2;
          }

          //#pragma omp barrier
          TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
          TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

          icnt = 0;
          for (ib = 0; ib < Check::sdim; ib++) {
            icnt += child_omp_sz(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_jb);
          }
          break;
        }
        else if (hacker == 1) {
          // this part can not be parallelized
          jb = 0;

          for (ib = 0; ib < Check::sdim; ib++) {
            list_jb[ib] = jb;

            i = ib * ihfbit;
            num_up = 0;
            for (j = 0; j <= N2 - 2; j += 2) {
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_up += div;
            }
            num_down = 0;
            for (j = 1; j <= N2 - 1; j += 2) {
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_down += div;
            }

            tmp_res = Def::Nsite % 2; // even Ns-> 0, odd Ns -> 1
            all_up = (Def::Nsite + tmp_res) / 2;
            all_down = (Def::Nsite - tmp_res) / 2;

            tmp_1 = Binomial(all_up, Def::Nup - num_up, comb, all_up);
            tmp_2 = Binomial(all_down, Def::Ndown - num_down, comb, all_down);
            jb += tmp_1 * tmp_2;
          }

          //#pragma omp barrier
          TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
          TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_jb, ihfbit, Check::sdim)
          for (ib = 0; ib < Check::sdim; ib++) {
            icnt += child_omp_sz_hacker(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_jb);
          }
          break;
        }
        else {
          fprintf(stderr, "Error: CalcHS in ModPara file must be 0 or 1 for Hubbard model.");
          return -1;
        }

      case DC::HubbardNConserved:
        hacker = Def::read_hacker;
        if (hacker == 0) {
          // this part can not be parallelized
          jb = 0;
          iSpnup = 0;
          iMinup = 0;
          iAllup = Def::Ne;
          if (Def::Ne > Def::Nsite) {
            iMinup = Def::Ne - Def::Nsite;
            iAllup = Def::Nsite;
          }
          for (ib = 0; ib < Check::sdim; ib++) {
            list_jb[ib] = jb;
            i = ib * ihfbit;
            num_up = 0;
            for (j = 0; j <= N2 - 2; j += 2) {
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_up += div;
            }
            num_down = 0;
            for (j = 1; j <= N2 - 1; j += 2) {
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_down += div;
            }
            tmp_res = Def::Nsite % 2; // even Ns-> 0, odd Ns -> 1
            all_up = (Def::Nsite + tmp_res) / 2;
            all_down = (Def::Nsite - tmp_res) / 2;

            for (iSpnup = iMinup; iSpnup <= iAllup; iSpnup++) {
              tmp_1 = Binomial(all_up, iSpnup - num_up, comb, all_up);
              tmp_2 = Binomial(all_down, Def::Ne - iSpnup - num_down, comb, all_down);
              jb += tmp_1 * tmp_2;
            }
          }
          //#pragma omp barrier
          TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
          TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_jb, ihfbit, N2, Check::sdim) 
          for (ib = 0; ib < Check::sdim; ib++) {
            icnt += child_omp_sz(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_jb);
          }
          break;
        }
        else if (hacker == 1) {
          // this part can not be parallelized
          jb = 0;
          iSpnup = 0;
          iMinup = 0;
          iAllup = Def::Ne;
          if (Def::Ne > Def::Nsite) {
            iMinup = Def::Ne - Def::Nsite;
            iAllup = Def::Nsite;
          }
          for (ib = 0; ib < Check::sdim; ib++) {
            list_jb[ib] = jb;
            i = ib * ihfbit;
            num_up = 0;
            for (j = 0; j <= N2 - 2; j += 2) {
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_up += div;
            }
            num_down = 0;
            for (j = 1; j <= N2 - 1; j += 2) {
              div = i & Def::Tpow[j];
              div = div / Def::Tpow[j];
              num_down += div;
            }
            tmp_res = Def::Nsite % 2; // even Ns-> 0, odd Ns -> 1
            all_up = (Def::Nsite + tmp_res) / 2;
            all_down = (Def::Nsite - tmp_res) / 2;

            for (iSpnup = iMinup; iSpnup <= iAllup; iSpnup++) {
              tmp_1 = Binomial(all_up, iSpnup - num_up, comb, all_up);
              tmp_2 = Binomial(all_down, Def::Ne - iSpnup - num_down, comb, all_down);
              jb += tmp_1 * tmp_2;
            }
          }
          //#pragma omp barrier
          TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
          TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_jb, ihfbit, N2, Check::sdim) 
          for (ib = 0; ib < Check::sdim; ib++) {
            icnt += child_omp_sz_hacker(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_jb);
          }

          break;
        }
        else {
          fprintf(stderr, "Error: CalcHS in ModPara file must be 0 or 1 for Hubbard model.");
          return -1;
        }

      case DC::Kondo:
        // this part can not be parallelized
        N_all_up = Def::Nup;
        N_all_down = Def::Ndown;
        fprintf(MP::STDOUT, "  N_all_up = %d N_all_down = %d \n", N_all_up, N_all_down);

        jb = 0;
        num_loc = 0;
        for (j = Def::Nsite / 2; j < Def::Nsite; j++) {// counting localized # of spins
          if (Def::LocSpn[j] != ITINERANT) {
            num_loc += 1;
          }
        }

        for (ib = 0; ib < Check::sdim; ib++) { //sdim = 2^(N/2)
          list_jb[ib] = jb;
          i = ib * ihfbit; // ihfbit=pow(2,((Nsite+1)/2))
          num_up = 0;
          num_down = 0;
          icheck_loc = 1;

          for (j = Def::Nsite / 2; j < Def::Nsite; j++) {
            div_up = i & Def::Tpow[2 * j];
            div_up = div_up / Def::Tpow[2 * j];
            div_down = i & Def::Tpow[2 * j + 1];
            div_down = div_down / Def::Tpow[2 * j + 1];
            if (Def::LocSpn[j] == ITINERANT) {
              num_up += div_up;
              num_down += div_down;
            }
            else {
              num_up += div_up;
              num_down += div_down;
              if (Def::Nsite % 2 == 1 && j == (Def::Nsite / 2)) { // odd site
                icheck_loc = icheck_loc;
                ihfSpinDown = div_down;
                if (div_down == 0) {
                  num_up += 1;
                }
              }
              else {
                icheck_loc = icheck_loc * (div_up^div_down);// exclude empty or doubly occupied site
              }
            }
          }

          if (icheck_loc == 1) { // itinerant of local spins without holon or doublon
            tmp_res = Def::Nsite % 2; // even Ns-> 0, odd Ns -> 1
            all_loc = Def::NLocSpn - num_loc; // # of local spins
            all_up = (Def::Nsite + tmp_res) / 2 - all_loc;
            all_down = (Def::Nsite - tmp_res) / 2 - all_loc;
            if (Def::Nsite % 2 == 1 && Def::LocSpn[Def::Nsite / 2] != ITINERANT) {
              all_up = (Def::Nsite) / 2 - all_loc;
              all_down = (Def::Nsite) / 2 - all_loc;
            }

            for (num_loc_up = 0; num_loc_up <= all_loc; num_loc_up++) {
              tmp_1 = Binomial(all_loc, num_loc_up, comb, all_loc);
              if (Def::Nsite % 2 == 1 && Def::LocSpn[Def::Nsite / 2] != ITINERANT) {
                if (ihfSpinDown != 0) {
                  tmp_2 = Binomial(all_up, Def::Nup - num_up - num_loc_up, comb, all_up);
                  tmp_3 = Binomial(all_down, Def::Ndown - num_down - (all_loc - num_loc_up), comb, all_down);
                }
                else {
                  tmp_2 = Binomial(all_up, Def::Nup - num_up - num_loc_up, comb, all_up);
                  tmp_3 = Binomial(all_down, Def::Ndown - num_down - (all_loc - num_loc_up), comb, all_down);
                }
              }
              else {
                tmp_2 = Binomial(all_up, Def::Nup - num_up - num_loc_up, comb, all_up);
                tmp_3 = Binomial(all_down, Def::Ndown - num_down - (all_loc - num_loc_up), comb, all_down);
              }
              jb += tmp_1 * tmp_2*tmp_3;
            }
          }

        }
        //#pragma omp barrier
        TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
        TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

        hacker = Def::read_hacker;
        if (hacker == 0) {
          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_jb, ihfbit, N2, Check::sdim)
          for (ib = 0; ib < Check::sdim; ib++) {
            icnt += child_omp_sz_Kondo(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_jb);
          }
        }
        else if (hacker == 1) {
          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_jb, ihfbit, N2, Check::sdim)
          for (ib = 0; ib < Check::sdim; ib++) {
            icnt += child_omp_sz_Kondo_hacker(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_jb);
          }
        }
        break;

      case DC::Spin:
        // this part can not be parallelized
        if (Def::iFlgGeneralSpin == FALSE) {
          hacker = Def::read_hacker;
          //printf(" rank=%d:Ne=%ld ihfbit=%ld sdim=%ld\n", MP::myrank,Def::Ne,ihfbit,Check::sdim);
          // using hacker's delight only + no open mp 
          if (hacker == -1) {
            icnt = 1;
            tmp_pow = 1;
            tmp_i = 0;
            jb = 0;
            ja = 0;
            while (tmp_pow < Def::Tpow[Def::Ne]) {
              tmp_i += tmp_pow;
              tmp_pow = tmp_pow * 2;
            }
            //printf("DEBUG: %ld %ld %ld %ld\n",tmp_i,Check::sdim,Def::Tpow[Def::Ne],Def::Nsite);
            if (Def::Nsite % 2 == 0) {
              max_tmp_i = Check::sdim*Check::sdim;
            }
            else {
              max_tmp_i = Check::sdim*Check::sdim * 2 - 1;
            }
            while (tmp_i < max_tmp_i) {
              list_1_[icnt] = tmp_i;

              ia = tmp_i & irght;
              ib = tmp_i & ilft;
              ib = ib / ihfbit;
              if (ib == ibpatn) {
                ja = ja + 1;
              }
              else {
                ibpatn = ib;
                ja = 1;
                jb = icnt - 1;
              }

              list_2_1_[ia] = ja + 1;
              list_2_2_[ib] = jb + 1;
              tmp_j = snoob(tmp_i);
              tmp_i = tmp_j;
              icnt += 1;
            }
            icnt = icnt - 1;
            // old version + hacker's delight
          }
          else if (hacker == 1) {
            jb = 0;
            for (ib = 0; ib < Check::sdim; ib++) {
              list_jb[ib] = jb;
              i = ib * ihfbit;
              num_up = 0;
              for (j = 0; j < N; j++) {
                div_up = i & Def::Tpow[j];
                div_up = div_up / Def::Tpow[j];
                num_up += div_up;
              }
              all_up = (Def::Nsite + 1) / 2;
              tmp_1 = Binomial(all_up, Def::Ne - num_up, comb, all_up);
              jb += tmp_1;
            }
            //#pragma omp barrier
            TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
            TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

            icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(ihfbit, N, list_1_, list_2_1_, list_2_2_, list_jb, Check::sdim)
            for (ib = 0; ib < Check::sdim; ib++) {
              icnt += child_omp_sz_spin_hacker(ib, ihfbit, N, list_1_, list_2_1_, list_2_2_, list_jb);
            }
            //printf(" rank=%d ib=%ld:Ne=%d icnt=%ld :idim_max=%ld N=%d\n", MP::myrank,ib,Def::Ne,icnt,Check::idim_max,N);
            // old version
          }
          else if (hacker == 0) {
            jb = 0;
            for (ib = 0; ib < Check::sdim; ib++) {
              list_jb[ib] = jb;
              i = ib * ihfbit;
              num_up = 0;
              for (j = 0; j < N; j++) {
                div_up = i & Def::Tpow[j];
                div_up = div_up / Def::Tpow[j];
                num_up += div_up;
              }
              all_up = (Def::Nsite + 1) / 2;
              tmp_1 = Binomial(all_up, Def::Ne - num_up, comb, all_up);
              jb += tmp_1;
            }
            //#pragma omp barrier
            TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
            TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

            icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_jb, ihfbit, N, Check::sdim)
            for (ib = 0; ib < Check::sdim; ib++) {
              icnt += child_omp_sz_spin(ib, ihfbit, N, list_1_, list_2_1_, list_2_2_, list_jb);
            }
          }
          else {
            fprintf(stderr, "Error: CalcHS in ModPara file must be -1 or 0 or 1 for Spin model.");
            return -1;
          }
        }
        else {
          int Max2Sz = 0;
          int irghtsite = 1;
          long int itmpSize = 1;
          int i2Sz = 0;
          for (j = 0; j < Def::Nsite; j++) {
            itmpSize *= Def::SiteToBit[j];
            if (itmpSize == ihfbit) {
              break;
            }
            irghtsite++;
          }
          for (j = 0; j < Def::Nsite; j++) {
            Max2Sz += Def::LocSpn[j];
          }

          HilbertNumToSz = li_1d_allocate(2 * Max2Sz + 1);
          for (ib = 0; ib < 2 * Max2Sz + 1; ib++) {
            HilbertNumToSz[ib] = 0;
          }

          for (ib = 0; ib < ihfbit; ib++) {
            i2Sz = 0;
            for (j = 0; j < irghtsite; j++) {
              i2Sz += GetLocal2Sz(j, ib, Def::SiteToBit, Def::Tpow);
            }
            list_2_1_Sz[ib] = i2Sz;
            HilbertNumToSz[i2Sz + Max2Sz]++;
          }
          jb = 0;
          long int ilftdim = (Def::Tpow[Def::Nsite - 1] * Def::SiteToBit[Def::Nsite - 1]) / ihfbit;
          for (ib = 0; ib < ilftdim; ib++) {
            list_jb[ib] = jb;
            i2Sz = 0;
            for (j = irghtsite; j < N; j++) {
              i2Sz += GetLocal2Sz(j, ib*ihfbit, Def::SiteToBit, Def::Tpow);
            }
            list_2_2_Sz[ib] = i2Sz;
            if ((Def::Total2Sz - i2Sz + (int)Max2Sz) >= 0 && (Def::Total2Sz - i2Sz) <= (int)Max2Sz) {
              jb += HilbertNumToSz[Def::Total2Sz - i2Sz + Max2Sz];
            }
          }

          TimeKeeper("%s_sz_TimeKeeper.dat", "mid omp parallel sz : %s", "a");
          TimeKeeper("%s_TimeKeeper.dat", "mid omp parallel sz : %s", "a");

          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) \
shared(list_1_, list_2_1_, list_2_2_, list_2_1_Sz, list_2_2_Sz,list_jb, ilftdim, ihfbit)
          for (ib = 0; ib < ilftdim; ib++) {
            icnt += child_omp_sz_GeneralSpin(ib, ihfbit, list_1_, list_2_1_, list_2_2_, list_2_1_Sz, list_2_2_Sz, list_jb);
          }
          free_li_1d_allocate(HilbertNumToSz);
        }

        break;
      default:
        wrapperMPI::Exit(-1);
      }
      i_max = icnt;
      //fprintf(MP::STDOUT, "Debug: Xicnt=%ld \n",icnt);
      TimeKeeper("%s_sz_TimeKeeper.dat", "omp parallel sz finishes: %s", "a");
      TimeKeeper("%s_TimeKeeper.dat", "omp parallel sz finishes: %s", "a");
    }

    if (Def::iFlgCalcSpec == DC::CALCSPEC_NOT) {
      if (Def::iCalcModel == DC::HubbardNConserved) {
        Def::iCalcModel = DC::Hubbard;
      }
    }

    //Error message
    //i_max=i_max+1;
    if (i_max != Check::idim_max) {
      fprintf(stderr, "%s", "Error: in sz. \n");
      fprintf(stderr, "imax = %ld, Check.idim_max=%ld \n", i_max, Check::idim_max);
      strcpy(sdt_err, "Err_sz.dat");
      if (childfopenMPI(sdt_err, "a", &fp_err) != 0) {
        wrapperMPI::Exit(-1);
      }
      fprintf(fp_err, "%s", "Caution!!  Error in sz !!!! idim_max is not correct \n");
      fclose(fp_err);
      wrapperMPI::Exit(-1);
    }

    free_li_2d_allocate(comb);
  }
  fprintf(MP::STDOUT, "%s", "  End  : Calculate HilbertNum for fixed Sz. \n\n");

  free(list_jb);
  if (Def::iFlgGeneralSpin == TRUE) {
    free(list_2_1_Sz);
    free(list_2_2_Sz);
  }
  return 0;
}
