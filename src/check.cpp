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
#include "sz.hpp"
#include "FileIO.hpp"
#include "common/setmemory.hpp"
#include "check.hpp"
#include "wrapperMPI.hpp"
#include "CheckMPI.hpp"
#include "global.hpp"
#include "DefCommon.hpp"

/**
 * @file
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for giving a function of calculating size of Hilbert space.
 * 
 */


/** 
 * @brief A program to check size of dimension for Hilbert-space.
 * 
 * @param[in,out] X  Common data set used in HPhi.
 * 
 * @retval TRUE normal termination
 * @retval FALSE abnormal termination
 * @retval MPIFALSE CheckMPI abnormal termination
 * @version 0.2
 * @details add function of calculating Hilbert space for canonical ensemble.
 *  
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int check(){
    
  FILE *fp;
  long int i,tmp_sdim;
  int NLocSpn, NCond;
  int Nup, Ndown;
  long int u_tmp;
  long int tmp;
  long int Ns,comb_1,comb_2,comb_3,comb_sum, comb_up, comb_down;
  int u_loc;
  long int **comb;    
  long int idimmax=0;
  long int idim=0;
  long int isite=0;
  int tmp_sz=0;
  int iMinup=0;
  if(Def::iCalcModel ==Spin ||Def::iCalcModel ==SpinGC )
  {
    Def::Ne=Def::Nup;
  }

  int iAllup=Def::Ne;

  if(Def::iFlgScaLAPACK == 0) {
    /*
      Set Site number per MPI process
    */
    if (CheckMPI() != TRUE) {
      return MPIFALSE;
    }
  }
  else{
    Def::NsiteMPI = Def::Nsite;
    Def::Total2SzMPI = Def::Total2Sz;
  }

  Ns = Def::Nsite;

  comb = li_2d_allocate(Ns+1,Ns+1);

  //idim_max
  switch(Def::iCalcModel){
  case HubbardGC:
    //comb_sum = 2^(2*Ns)=4^Ns
    comb_sum = 1;
    for(i=0;i<2*Def::Nsite;i++){
      comb_sum= 2*comb_sum;     
    }
    break;
  case SpinGC:
    //comb_sum = 2^(Ns)
    comb_sum = 1;
    if(Def::iFlgGeneralSpin ==FALSE){
      for(i=0;i<Def::Nsite;i++){
        comb_sum= 2*comb_sum;     
      }
    }
    else{
      for(i=0; i<Def::Nsite;i++){
        comb_sum=comb_sum*Def::SiteToBit[i];
      }
    }
    break;

  case Hubbard:
    comb_up= Binomial(Ns, Def::Nup, comb, Ns);
    comb_down= Binomial(Ns, Def::Ndown, comb, Ns);
    comb_sum=comb_up*comb_down;
    break;

  case HubbardNConserved:
    comb_sum=0;
    if(Def::Ne > Def::Nsite){
      iMinup = Def::Ne-Def::Nsite;
      iAllup = Def::Nsite;
    }

    for(i=iMinup; i<= iAllup; i++){
      comb_up= Binomial(Ns, i, comb, Ns);
      comb_down= Binomial(Ns, Def::Ne-i, comb, Ns);
      comb_sum +=comb_up*comb_down;
    }
    break;
    
  case Kondo:
    Nup     = Def::Nup;
    Ndown   = Def::Ndown;
    NCond   = Def::Nsite-Def::NLocSpn;
    NLocSpn = Def::NLocSpn;
    comb_sum = 0;
    for(u_loc=0;u_loc<=Def::Nup;u_loc++){
      comb_1     = Binomial(NLocSpn,u_loc,comb,Ns);
      comb_2     = Binomial(NCond,Nup-u_loc,comb,Ns);
      comb_3     = Binomial(NCond,Ndown+u_loc-NLocSpn,comb,Ns);
      comb_sum  += comb_1*comb_2*comb_3;
    }
    break;
  case KondoGC:
    comb_sum = 1;
    NCond   = Def::Nsite-Def::NLocSpn;
    NLocSpn = Def::NLocSpn;
    //4^Nc*2^Ns
    for(u_loc=0;u_loc <(2*NCond+NLocSpn); u_loc++){
      comb_sum= 2*comb_sum;     
    }
    break;
  case Spin:

    if(Def::iFlgGeneralSpin ==FALSE){
      if(Def::Nup+Def::Ndown != Def::Nsite){
        fprintf(stderr, " 2Sz is incorrect.\n");
        return FALSE;
      }
      //comb_sum= Binomial(Ns, Def::Ne, comb, Ns);
      comb_sum= Binomial(Ns, Def::Nup, comb, Ns);
    }
    else{
      idimmax = 1;
      Def::Tpow[0]=idimmax;
      for(isite=0; isite<Def::Nsite;isite++){
        idimmax=idimmax*Def::SiteToBit[isite];
        Def::Tpow[isite+1]=idimmax;
      }
      comb_sum=0;
#pragma omp parallel for default(none) reduction(+:comb_sum) private(tmp_sz, isite) \
shared(idimmax, Def::Nsite, Def::Tpow, Def::SiteToBit, Def::Total2Sz) 
      for (idim = 0; idim < idimmax; idim++) {
        tmp_sz = 0;
        for (isite = 0; isite < Def::Nsite; isite++) {
          tmp_sz += GetLocal2Sz(isite + 1, idim, Def::SiteToBit, Def::Tpow);
        }
        if(tmp_sz == Def::Total2Sz){
          comb_sum +=1;
        }
      }
    }
    
    break;
  default:
    fprintf(stderr, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
    free_li_2d_allocate(comb);
    return FALSE;
  }  

  //fprintf(stdoutMPI, "Debug: comb_sum= %ld \n",comb_sum);

  Check::idim_max = comb_sum;
  switch(Def::iCalcType) {
    case CG:
      switch (Def::iCalcModel) {
        case Hubbard:
        case HubbardNConserved:
        case Kondo:
        case KondoGC:
        case Spin:
          Check::max_mem = (7 * Def::k_exct + 1.5) * Check::idim_max * 16.0 / (pow(10, 9));
          break;
        case HubbardGC:
        case SpinGC:
          Check::max_mem = (7 * Def::k_exct + 1.0) * Check::idim_max * 16.0 / (pow(10, 9));
          break;
      }
      break;
    case TPQCalc:
      switch (Def::iCalcModel) {
        case Hubbard:
        case HubbardNConserved:
        case Kondo:
        case KondoGC:
        case Spin:
          if (Def::iFlgCalcSpec != CALCSPEC_NOT) {
            Check::max_mem = NumAve * 3 * Check::idim_max * 16.0 / (pow(10, 9));
          } else {
            Check::max_mem = 4.5 * Check::idim_max * 16.0 / (pow(10, 9));
          }
          break;
        case HubbardGC:
        case SpinGC:
          if (Def::iFlgCalcSpec != CALCSPEC_NOT) {
            Check::max_mem = NumAve * 3 * Check::idim_max * 16.0 / (pow(10, 9));
          } else {
            Check::max_mem = 3.5 * Check::idim_max * 16.0 / (pow(10, 9));
          }
          break;
      }
      break;
    case FullDiag:
      Check::max_mem = Check::idim_max * 8.0 * Check::idim_max * 8.0 / (pow(10, 9));
      break;
    case TimeEvolution:
      Check::max_mem = (4 + 2 + 1) * Check::idim_max * 16.0 / (pow(10, 9));
      break;
    default:
      return FALSE;
      //break;
  }

  //fprintf(stdoutMPI, "  MAX DIMENSION idim_max=%ld \n",Check::idim_max);
  //fprintf(stdoutMPI, "  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",Check::max_mem);
  long int li_dim_max=MaxMPI_li(Check::idim_max);
  fprintf(stdoutMPI, "  MAX DIMENSION idim_max=%ld \n",li_dim_max);
  double dmax_mem=MaxMPI_d(Check::max_mem);
  fprintf(stdoutMPI, "  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",dmax_mem);
  if(childfopenMPI("CHECK_Memory.dat","w", &fp)!=0){
    free_li_2d_allocate(comb);
    return FALSE;
  }
  fprintf(fp,"  MAX DIMENSION idim_max=%ld \n", li_dim_max);
  fprintf(fp,"  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n", dmax_mem);

  
  /*
  fprintf(fp,"  MAX DIMENSION idim_max=%ld \n",Check::idim_max);
  fprintf(fp,"  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",Check::max_mem);
  */
  fclose(fp);

  //sdim 
  tmp=1;
  tmp_sdim=1;

  switch(Def::iCalcModel){
  case HubbardGC:
  case KondoGC:
  case HubbardNConserved:
  case Hubbard:
  case Kondo:
    while(tmp <= Def::Nsite){
      tmp_sdim=tmp_sdim*2;
      tmp+=1;
    }
    break;
  case Spin:
  case SpinGC:
    if(Def::iFlgGeneralSpin==FALSE){ 
      while(tmp <= Def::Nsite/2){
        tmp_sdim=tmp_sdim*2;
        tmp+=1;
      }
    }
    else{
      GetSplitBitForGeneralSpin(Def::Nsite, &tmp_sdim, Def::SiteToBit);
    }
    break;
  default:
    fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
    free_li_2d_allocate(comb);
    return FALSE;
  }  
  Check::sdim=tmp_sdim;
  
  if(childfopenMPI("CHECK_Sdim.dat","w", &fp)!=0){
    free_li_2d_allocate(comb);
    return FALSE;
  }

  switch(Def::iCalcModel){
  case HubbardGC:
  case KondoGC:
  case HubbardNConserved:
  case Hubbard:
  case Kondo:
    //fprintf(stdoutMPI, "sdim=%ld =2^%d\n",Check::sdim,Def::Nsite);
    fprintf(fp,"sdim=%ld =2^%d\n",Check::sdim,Def::Nsite);
    break;
  case Spin:
  case SpinGC:
    if(Def::iFlgGeneralSpin==FALSE){
      //fprintf(stdoutMPI, "sdim=%ld =2^%d\n",Check::sdim,Def::Nsite/2);
      fprintf(fp,"sdim=%ld =2^%d\n",Check::sdim,Def::Nsite/2);
    }
    break;
  default:
    break;
  }

  free_li_2d_allocate(comb);

  u_tmp=1;
  Def::Tpow[0]=u_tmp;
  switch(Def::iCalcModel){
  case HubbardGC:
  case KondoGC:
    for(i=1;i<=2*Def::Nsite;i++){
      u_tmp=u_tmp*2;
      Def::Tpow[i]=u_tmp;
      fprintf(fp,"%ld %ld \n",i,u_tmp);
    }
    break;
  case HubbardNConserved:
  case Hubbard:
  case Kondo:
    for(i=1;i<=2*Def::Nsite-1;i++){
      u_tmp=u_tmp*2;
      Def::Tpow[i]=u_tmp;
      fprintf(fp,"%ld %ld \n",i,u_tmp);
    }
    break;
 case SpinGC:
   if(Def::iFlgGeneralSpin==FALSE){
     for(i=1;i<=Def::Nsite;i++){
       u_tmp=u_tmp*2;
       Def::Tpow[i]=u_tmp;
       fprintf(fp,"%ld %ld \n",i,u_tmp);
     }
   }
   else{
     Def::Tpow[0]=u_tmp;
     fprintf(fp,"%d %ld \n", 0, u_tmp);
     for(i=1;i<Def::Nsite;i++){
       u_tmp=u_tmp*Def::SiteToBit[i-1];
       Def::Tpow[i]=u_tmp;
       fprintf(fp,"%ld %ld \n",i,u_tmp);
     }
   }
   break;
 case Spin:
   if(Def::iFlgGeneralSpin==FALSE){
     for(i=1;i<=Def::Nsite-1;i++){
       u_tmp=u_tmp*2;
       Def::Tpow[i]=u_tmp;
       fprintf(fp,"%ld %ld \n",i,u_tmp);
     }
   }
   else{
     for(i=0;i<Def::Nsite;i++){
       fprintf(fp,"%ld %ld \n",i,Def::Tpow[i]);
     }
   }     
    break;
  default:
    fprintf(stdoutMPI, "Error: CalcModel %d is incorrect.\n", Def::iCalcModel);
    free_li_2d_allocate(comb);
    return FALSE;
  }  
  fclose(fp);
  /*
    Print MPI-site information and Modify Tpow 
    in the inter process region.
  */
  CheckMPI_Summary();
  
  return TRUE;
}    
