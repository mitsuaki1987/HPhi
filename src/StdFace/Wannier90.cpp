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
@brief Standard mode for wannier90
*/
#include "StdFace_vals.hpp"
#include "StdFace_ModelUtil.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <cstring>
/**
@brief Read Geometry file for wannier90
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void geometry_W90()
{
  int isite, ii, ierr;
  char filename[265];
  FILE *fp;

  sprintf(filename, "%s_geom.dat", StdI::CDataFileHead);
  fprintf(stdout, "    Wannier90 Geometry file = %s\n", filename);

  fp = fopen(filename, "r");
  /**@brief
   Direct lattice vector StdI::direct
  */
  for (ii = 0; ii < 3; ii++)
    ierr = fscanf(fp, "%lf%lf%lf", &StdI::direct[ii][0], &StdI::direct[ii][1], &StdI::direct[ii][2]);
  if(ierr == EOF) printf("%d\n", ierr);
  /**@brief
   Intrinsic site position StdI::tau and its number StdI::NsiteUC
  */
  for (isite = 0; isite < StdI::NsiteUC; isite++) free(StdI::tau[isite]);
  free(StdI::tau);
  ierr = fscanf(fp, "%d", &StdI::NsiteUC);
  fprintf(stdout, "    Number of Correlated Sites = %d\n", StdI::NsiteUC);

  StdI::tau = (double **)malloc(sizeof(double*) * StdI::NsiteUC);
  for (ii = 0; ii < StdI::NsiteUC; ii++) StdI::tau[ii] = (double *)malloc(sizeof(double) * 3);

  for (isite = 0; isite < StdI::NsiteUC; isite++)
    ierr = fscanf(fp, "%lf%lf%lf", &StdI::tau[isite][0], &StdI::tau[isite][1], &StdI::tau[isite][2]);
  fclose(fp);

  printf("    Direct lattice vectors:\n");
  for (ii = 0; ii < 3; ii++) printf("      %10.5f %10.5f %10.5f\n",
    StdI::direct[ii][0], StdI::direct[ii][1], StdI::direct[ii][2]);
  printf("    Wannier centres:\n");
  for (isite = 0; isite < StdI::NsiteUC; isite++) printf("      %10.5f %10.5f %10.5f\n",
    StdI::tau[isite][0], StdI::tau[isite][1], StdI::tau[isite][2]);

}/*static void geometry_W90() */
/**
 @brief Read Wannier90 hamiltonian file (*_hr)
 @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void read_W90(
  char *filename,//!<[in] Input file name
  double cutoff,//!<[in] Threshold for the Hamiltonian 
  int *cutoff_R,
  double cutoff_length,
  int itUJ,
  int *NtUJ,
  int ***tUJindx,//!<[out] R, band index of matrix element
  std::complex<double> **tUJ//!<[out] Matrix element
)
{
  FILE *fp;
  int ierr, nWan, nWSC, iWSC, jWSC, iWan, jWan, iWan0, jWan0, ii, jj;
  double dtmp[2], dR[3], length;
  char ctmp[256], *ctmp2;
  std::complex<double> ***Mat_tot;
  int **indx_tot;
  /*
  Header part
  */
  fp = fopen(filename, "r");
  ctmp2 = fgets(ctmp, 256, fp);
  ierr = fscanf(fp, "%d", &nWan);
  if (ierr == EOF) printf("%d %s\n", ierr, ctmp2);
  ierr = fscanf(fp, "%d", &nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    ierr = fscanf(fp, "%d", &ii);
  }
  /*
  Malloc Matrix elements and their indices
  */
  Mat_tot = (std::complex<double> ***)malloc(sizeof(std::complex<double> **) * nWSC);
  indx_tot = (int **)malloc(sizeof(int*) * nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    Mat_tot[iWSC] = (std::complex<double> **)malloc(sizeof(std::complex<double> *) * nWan);
    indx_tot[iWSC] = (int *)malloc(sizeof(int) * 3);
    for (iWan = 0; iWan < nWan; iWan++) {
      Mat_tot[iWSC][iWan] = (std::complex<double> *)malloc(sizeof(std::complex<double>) * nWan);
    }
  }
  /*
  Read body
  */
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      for (jWan = 0; jWan < nWan; jWan++) {
        ierr = fscanf(fp, "%d%d%d%d%d%lf%lf",
          &indx_tot[iWSC][0], &indx_tot[iWSC][1], &indx_tot[iWSC][2],
          &iWan0, &jWan0,
          &dtmp[0], &dtmp[1]);
        /*
         Compute Euclid length
        */
        for (ii = 0; ii < 3; ii++) {
          dR[ii] = 0.0;
          for (jj = 0; jj < 3; jj++)
            dR[ii] += StdI::direct[jj][ii] * (StdI::tau[jWan][jj] - StdI::tau[iWan][jj] + indx_tot[iWSC][jj]);
        }
        length = sqrt(dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2]);
        if (length > cutoff_length && cutoff_length > 0.0) {
          dtmp[0] = 0.0;
          dtmp[1] = 0.0;
        }
        if (abs(indx_tot[iWSC][0]) > cutoff_R[0] ||
          abs(indx_tot[iWSC][1]) > cutoff_R[1] ||
          abs(indx_tot[iWSC][2]) > cutoff_R[2]) {
          dtmp[0] = 0.0;
          dtmp[1] = 0.0;
        }
        if (iWan0 <= StdI::NsiteUC && jWan0 <= StdI::NsiteUC)
          Mat_tot[iWSC][iWan0 - 1][jWan0 - 1] = std::complex<double>(dtmp[0], dtmp[1]);
      }
    }
    /**@brief
    (1) Apply inversion symmetry and delete duplication
    */
    for (jWSC = 0; jWSC < iWSC; jWSC++) {
      if (
        indx_tot[iWSC][0] == -indx_tot[jWSC][0] &&
        indx_tot[iWSC][1] == -indx_tot[jWSC][1] &&
        indx_tot[iWSC][2] == -indx_tot[jWSC][2]
        )
        for (iWan = 0; iWan < StdI::NsiteUC; iWan++) {
          for (jWan = 0; jWan < StdI::NsiteUC; jWan++) {
            Mat_tot[iWSC][iWan][jWan] = 0.0;
          }
        }
    }/*for (jWSC = 0; jWSC < iWSC; jWSC++)*/
    if (indx_tot[iWSC][0] == 0 &&
      indx_tot[iWSC][1] == 0 &&
      indx_tot[iWSC][2] == 0)
      for (iWan = 0; iWan < StdI::NsiteUC; iWan++) {
        for (jWan = 0; jWan < iWan; jWan++) {
          Mat_tot[iWSC][iWan][jWan] = 0.0;
        }
      }
  }/*for (iWSC = 0; iWSC < nWSC; iWSC++)*/
  fclose(fp);
  /**@brief
  (3-1)  Compute the number of terms lerger than cut-off.
  */
  fprintf(stdout, "\n      EFFECTIVE terms:\n");
  fprintf(stdout, "           R0   R1   R2 band_i band_f Hamiltonian\n");
  NtUJ[itUJ] = 0;
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < StdI::NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI::NsiteUC; jWan++) {
        if (cutoff < abs(Mat_tot[iWSC][iWan][jWan])) {
          fprintf(stdout, "        %5d%5d%5d%5d%5d%12.6f%12.6f\n",
            indx_tot[iWSC][0], indx_tot[iWSC][1], indx_tot[iWSC][2], iWan, jWan,
            real(Mat_tot[iWSC][iWan][jWan]), imag(Mat_tot[iWSC][iWan][jWan]));
          NtUJ[itUJ] += 1;
        }
      }
    }
  }
  fprintf(stdout, "      Total number of EFFECTIVE term = %d\n", NtUJ[itUJ]);
  tUJ[itUJ] = (std::complex<double> *)malloc(sizeof(std::complex<double>) * NtUJ[itUJ]);
  tUJindx[itUJ] = (int **)malloc(sizeof(int*) * NtUJ[itUJ]);
  for (ii = 0; ii < NtUJ[itUJ]; ii++) tUJindx[itUJ][ii] = (int *)malloc(sizeof(int) * 5);
  /**@brief
  Then Store to the hopping Integeral and
  its site index.
  */
  NtUJ[itUJ] = 0;
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < StdI::NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI::NsiteUC; jWan++) {
        if (cutoff < abs(Mat_tot[iWSC][iWan][jWan])) {
          for (ii = 0; ii < 3; ii++) tUJindx[itUJ][NtUJ[itUJ]][ii] = indx_tot[iWSC][ii];
          tUJindx[itUJ][NtUJ[itUJ]][3] = iWan;
          tUJindx[itUJ][NtUJ[itUJ]][4] = jWan;
          tUJ[itUJ][NtUJ[itUJ]] = Mat_tot[iWSC][iWan][jWan];
          NtUJ[itUJ] += 1;
        }
      }/*for (jWan = 0; jWan < StdI::NsiteUC; jWan++)*/
    }/*for (iWan = 0; iWan < StdI::NsiteUC; iWan++)*/
  }/*for (iWSC = 0; iWSC < nWSC; iWSC++)*/

  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      free(Mat_tot[iWSC][iWan]);
    }
    free(Mat_tot[iWSC]);
    free(indx_tot[iWSC]);
  }
  free(Mat_tot);
  free(indx_tot);
}/*static int read_W90(char *model)*/
/**
 @brief Read RESPACK Density-matrix file (*_dr.dat)
 @author Mitsuaki Kawamura (The University of Tokyo)
*/
static std::complex<double>***** read_density_matrix(
  char *filename//!<[in] Input file name
)
{
  FILE *fp;
  int ierr, nWan, nWSC, iWSC, iWan, jWan, iWan0, jWan0, ii, Rmax[3], Rmin[3], NR[3], i0, i1, i2;
  double dtmp[2];
  char ctmp[256], *ctmp2;
  std::complex<double> ***Mat_tot, *****DenMat;
  int **indx_tot;
  /*
  Header part
  */
  fp = fopen(filename, "r");
  ctmp2 = fgets(ctmp, 256, fp);
  ierr = fscanf(fp, "%d", &nWan);
  if (ierr == EOF) printf("%d %s\n", ierr, ctmp2);
  ierr = fscanf(fp, "%d", &nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    ierr = fscanf(fp, "%d", &ii);
  }
  /*
  Malloc Matrix elements and their indices
  */
  Mat_tot = (std::complex<double> ***)malloc(sizeof(std::complex<double> **) * nWSC);
  indx_tot = (int **)malloc(sizeof(int*) * nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    Mat_tot[iWSC] = (std::complex<double> **)malloc(sizeof(std::complex<double> *) * nWan);
    indx_tot[iWSC] = (int *)malloc(sizeof(int) * 3);
    for (iWan = 0; iWan < nWan; iWan++) {
      Mat_tot[iWSC][iWan] = (std::complex<double> *)malloc(sizeof(std::complex<double>) * nWan);
    }
  }
  /*
  Read body
  */
  for (ii = 0; ii < 3; ii++) {
    Rmin[ii] = 0;
    Rmax[ii] = 0;
  }
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      for (jWan = 0; jWan < nWan; jWan++) {
        ierr = fscanf(fp, "%d%d%d%d%d%lf%lf",
          &indx_tot[iWSC][0], &indx_tot[iWSC][1], &indx_tot[iWSC][2],
          &iWan0, &jWan0,
          &dtmp[0], &dtmp[1]);
        if (iWan0 <= StdI::NsiteUC && jWan0 <= StdI::NsiteUC)
          Mat_tot[iWSC][iWan0 - 1][jWan0 - 1] = std::complex<double>(dtmp[0], dtmp[1]);
        for (ii = 0; ii < 3; ii++) {
          if (indx_tot[iWSC][ii] < Rmin[ii]) Rmin[ii] = indx_tot[iWSC][ii];
          if (indx_tot[iWSC][ii] > Rmax[ii]) Rmax[ii] = indx_tot[iWSC][ii];
        }
      }
    }
  }/*for (iWSC = 0; iWSC < nWSC; iWSC++)*/
  fclose(fp);
  for (ii = 0; ii < 3; ii++) NR[ii] = Rmax[ii] - Rmin[ii] + 1;
  fprintf(stdout, "      Minimum R : %d %d %d\n", Rmin[0], Rmin[1], Rmin[2]);
  fprintf(stdout, "      Maximum R : %d %d %d\n", Rmax[0], Rmax[1], Rmax[2]);
  fprintf(stdout, "      Numver of R : %d %d %d\n", NR[0], NR[1], NR[2]);

  DenMat = (std::complex<double> *****)malloc(sizeof(std::complex<double> ****)*NR[0]);
  DenMat = DenMat - Rmin[0]; // base shift
  for (i0 = Rmin[0]; i0 <= Rmax[0]; i0++) {
    DenMat[i0] = (std::complex<double> ****)malloc(sizeof(std::complex<double> ***)*NR[1]);
    DenMat[i0] = DenMat[i0] - Rmin[1]; // base shift
    for (i1 = Rmin[1]; i1 <= Rmax[1]; i1++) {
      DenMat[i0][i1] = (std::complex<double> ***)malloc(sizeof(std::complex<double> **)*NR[2]);
      DenMat[i0][i1] = DenMat[i0][i1] - Rmin[2]; // base shift
      for (i2 = Rmin[2]; i2 <= Rmax[2]; i2++) {
        DenMat[i0][i1][i2] = (std::complex<double> **)malloc(sizeof(std::complex<double> *) * StdI::NsiteUC);
        for (iWan = 0; iWan < StdI::NsiteUC; iWan++) {
          DenMat[i0][i1][i2][iWan] = (std::complex<double> *)malloc(sizeof(std::complex<double>) * StdI::NsiteUC);
          for (jWan = 0; jWan < StdI::NsiteUC; jWan++)DenMat[i0][i1][i2][iWan][jWan] = 0.0;
        }
      }
    }
  }
  for (iWSC = 0; iWSC < nWSC; iWSC++) 
    for (iWan = 0; iWan < nWan; iWan++)
      for (jWan = 0; jWan < nWan; jWan++)
        DenMat[indx_tot[iWSC][0]][indx_tot[iWSC][1]][indx_tot[iWSC][2]][iWan][jWan]
        = Mat_tot[iWSC][iWan][jWan];

  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      free(Mat_tot[iWSC][iWan]);
    }
    free(Mat_tot[iWSC]);
    free(indx_tot[iWSC]);
  }
  free(Mat_tot);
  free(indx_tot);

  return DenMat;
}/*static int read_W90(char *model)*/
/**
@brief Print the initial guess of UHF
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintUHFinitial(
  int *NtUJ,
  std::complex<double> *****DenMat,
  int ***tUJindx
) {
  FILE *fp;
  int iW, iL, iH, kCell, it, isite, jsite ,NIniGuess, ispin;
  double dR[3];
  std::complex<double> Cphase, **IniGuess;

  IniGuess = (std::complex<double> **)malloc(sizeof(std::complex<double>*) * StdI::nsite);
  for (isite = 0; isite < StdI::nsite; isite++) 
    IniGuess[isite] = (std::complex<double> *)malloc(sizeof(std::complex<double>) * StdI::nsite);

  for (kCell = 0; kCell < StdI::NCell; kCell++) {
    /**/
    iW = StdI::Cell[kCell][0];
    iL = StdI::Cell[kCell][1];
    iH = StdI::Cell[kCell][2];
    /*
     Diagonal term
    */
    for (isite = 0; isite < StdI::NsiteUC; isite++) {
      jsite = isite + StdI::NsiteUC*kCell;
      IniGuess[jsite][jsite] = DenMat[0][0][0][isite][isite];
    }
    /*
    Coulomb integral (U)
    */
    for (it = 0; it < NtUJ[1]; it++) {
      StdFace::FindSite(iW, iL, iH,
        tUJindx[1][it][0], tUJindx[1][it][1], tUJindx[1][it][2],
        tUJindx[1][it][3], tUJindx[1][it][4], &isite, &jsite, &Cphase, dR);
      IniGuess[isite][jsite] = DenMat[tUJindx[1][it][0]][tUJindx[1][it][1]][tUJindx[1][it][2]]
                                     [tUJindx[1][it][3]][tUJindx[1][it][4]];
      IniGuess[jsite][isite] = conj(DenMat[tUJindx[1][it][0]][tUJindx[1][it][1]][tUJindx[1][it][2]]
                                          [tUJindx[1][it][3]][tUJindx[1][it][4]]);
    }/*for (it = 0; it < NtUJ[1]; it++)*/
    /*
    Wxchange integral (J)
    */
    for (it = 0; it < NtUJ[2]; it++) {
      StdFace::FindSite(iW, iL, iH,
        tUJindx[2][it][0], tUJindx[2][it][1], tUJindx[2][it][2],
        tUJindx[2][it][3], tUJindx[2][it][4], &isite, &jsite, &Cphase, dR);
      IniGuess[isite][jsite] = DenMat[tUJindx[2][it][0]][tUJindx[2][it][1]][tUJindx[2][it][2]]
                                     [tUJindx[2][it][3]][tUJindx[2][it][4]];
      IniGuess[jsite][isite] = conj(DenMat[tUJindx[2][it][0]][tUJindx[2][it][1]][tUJindx[2][it][2]]
                                          [tUJindx[2][it][3]][tUJindx[2][it][4]]);
    }/*for (it = 0; it < NtUJ[1]; it++)*/
  }/*for (kCell = 0; kCell < StdI::NCell; kCell++)*/

  NIniGuess = 0;
  for (isite = 0; isite < StdI::nsite; isite++) {
    for (jsite = 0; jsite < StdI::nsite; jsite++) {
      if (abs(IniGuess[isite][jsite]) > 1.0e-6)NIniGuess += 1;
    }/*for (jsite = 0; jsite < StdI::nsite; jsite++)*/
  }/*for (isite = 0; isite < StdI::nsite; isite++)*/

  fp = fopen("initial.def", "w");
  fprintf(fp, "======================== \n");
  fprintf(fp, "NInitialGuess %7d  \n", NIniGuess*2);
  fprintf(fp, "======================== \n");
  fprintf(fp, "========i_j_s_tijs====== \n");
  fprintf(fp, "======================== \n");

  for (isite = 0; isite < StdI::nsite; isite++) {
    for (jsite = 0; jsite < StdI::nsite; jsite++) {
      if (abs(IniGuess[isite][jsite]) > 1.0e-6) 
        for (ispin = 0; ispin < 2; ispin++) {
          fprintf(fp, "%5d %5d %5d %5d %25.15f %25.15f\n",
            jsite, ispin, isite, ispin,
            0.5*real(IniGuess[isite][jsite]),
            0.5*imag(IniGuess[isite][jsite]));
        }
    }/*for (jsite = 0; jsite < StdI::nsite; jsite++)*/
  }/*for (isite = 0; isite < StdI::nsite; isite++)*/

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "      initial.def is written.\n");

  for (isite = 0; isite < StdI::nsite; isite++)
    free(IniGuess[isite]);
  free(IniGuess);
}/*static void PrintTrans*/
/**
@brief Setup a Hamiltonian for the Wannier90 *_hr.dat
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace::Wannier90()
{
  int isite, jsite, ispin;
  int iL, iW, iH, kCell, it, ii;
  double Jtmp[3][3] = { {0.0} };
  FILE *fp;
  std::complex<double> Cphase, DenMat0;
  double dR[3], *Uspin;
  int NtUJ[3];
  std::complex<double> **tUJ, *****DenMat;
  int ***tUJindx;
  char filename[263];
  /**@brief
  (1) Compute the shape of the super-cell and sites in the super-cell
  */
  fp = fopen("lattice.xsf", "w");
  /**/
  StdFace::PrintVal_d("phase0", &StdI::phase[0], 0.0);
  StdFace::PrintVal_d("phase1", &StdI::phase[1], 0.0);
  StdFace::PrintVal_d("phase2", &StdI::phase[2], 0.0);
  StdI::NsiteUC = 1;
  StdFace::InitSite(fp, 3);
  fprintf(stdout, "\n  @ Wannier90 Geometry \n\n");
  geometry_W90();
  StdFace::PrintVal_i("DoubleCounting", &StdI::double_counting, 1);

  tUJ = (std::complex<double> **)malloc(sizeof(std::complex<double>*) * 3);
  tUJindx = (int ***)malloc(sizeof(int**) * 3);

  /*
  Read Hopping
  */
  fprintf(stdout, "\n  @ Wannier90 hopping \n\n");
  StdFace::PrintVal_d("cutoff_t", &StdI::cutoff_t, 1.0e-8);
  StdFace::PrintVal_d("cutoff_length_t", &StdI::cutoff_length_t, -1.0);
  sprintf(filename, "%s_hr.dat", StdI::CDataFileHead);
  read_W90(filename,
    StdI::cutoff_t, StdI::cutoff_tR, StdI::cutoff_length_t,
     0, NtUJ, tUJindx, tUJ);
  /*
  Read Coulomb
  */
  fprintf(stdout, "\n  @ Wannier90 Coulomb \n\n");
  StdFace::PrintVal_d("cutoff_u", &StdI::cutoff_u, 1.0e-8);
  StdFace::PrintVal_d("cutoff_length_U", &StdI::cutoff_length_U, -1.0);
  sprintf(filename, "%s_ur.dat", StdI::CDataFileHead);
  read_W90(filename,
    StdI::cutoff_u, StdI::cutoff_UR, StdI::cutoff_length_U, 
    1, NtUJ, tUJindx, tUJ);
  /*
  Read Hund
  */
  fprintf(stdout, "\n  @ Wannier90 Hund \n\n");
  StdFace::PrintVal_d("cutoff_j", &StdI::cutoff_j, 1.0e-8);
  StdFace::PrintVal_d("cutoff_length_J", &StdI::cutoff_length_J, -1.0);
  sprintf(filename, "%s_jr.dat", StdI::CDataFileHead);
  read_W90(filename,
    StdI::cutoff_j, StdI::cutoff_JR, StdI::cutoff_length_J, 
    2, NtUJ, tUJindx, tUJ);
  /*
  Read Density matrix
  */
  fprintf(stdout, "\n  @ Wannier90 Density-matrix \n\n");
  sprintf(filename, "%s_dr.dat", StdI::CDataFileHead);
  DenMat = read_density_matrix(filename);
  /**@brief
  (2) check & store parameters of Hamiltonian
  */
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  StdFace::NotUsed_d("K", StdI::K);
  StdFace::PrintVal_d("h", &StdI::h, 0.0);
  StdFace::PrintVal_d("Gamma", &StdI::Gamma, 0.0);
  StdFace::NotUsed_d("U", StdI::U);
  /**/
  if (strcmp(StdI::model, "spin") == 0 ) {
    StdFace::PrintVal_i("2S", &StdI::S2, 1);
  }/*if (strcmp(StdI::model, "spin") == 0 )*/
  else if (strcmp(StdI::model, "hubbard") == 0) {
    StdFace::PrintVal_d("mu", &StdI::mu, 0.0);
  }
  else{
    printf("wannier + Kondo is not available !\n");
    StdFace::exit(-1);
  }/*if (model != "spin")*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");
  /**@brief
  (3) Set local spin flag (StdI::locspinflag) and 
  the number of sites (StdI::nsite)
  */
  StdI::nsite = StdI::NsiteUC * StdI::NCell;
  StdI::locspinflag = (int *)malloc(sizeof(int) * StdI::nsite);
  /**/
  if(strcmp(StdI::model, "spin") == 0 )
    for (isite = 0; isite < StdI::nsite; isite++) StdI::locspinflag[isite] = StdI::S2;
  else if(strcmp(StdI::model, "hubbard") == 0 )
    for (isite = 0; isite < StdI::nsite; isite++) StdI::locspinflag[isite] = 0;
  /**@brief
  (4.5) For spin system, compute super exchange interaction.
  */
  if (strcmp(StdI::model, "spin") == 0) {
    Uspin = (double *)malloc(sizeof(double) * StdI::NsiteUC);
    for (it = 0; it < NtUJ[1]; it++)
      if (tUJindx[1][it][0] == 0 && tUJindx[1][it][1] == 0 && tUJindx[1][it][2] == 0
        && tUJindx[1][it][3] == tUJindx[1][it][4])     
        Uspin[tUJindx[1][it][3]] = real(tUJ[1][it]);
  }/*if (strcmp(StdI::model, "spin") == 0)*/
  /**@brief
  (5) Set Transfer & Interaction
  */
  for (kCell = 0; kCell < StdI::NCell; kCell++){
    /**/
    iW = StdI::Cell[kCell][0];
    iL = StdI::Cell[kCell][1];
    iH = StdI::Cell[kCell][2];
    /*
     Local term 1
    */
    if (strcmp(StdI::model, "spin") == 0) {
      for (isite = StdI::NsiteUC*kCell; isite < StdI::NsiteUC*(kCell + 1); isite++) {
        StdFace::MagField(StdI::S2, -StdI::h, -StdI::Gamma, isite);
      }
    }/*if (strcmp(StdI::model, "spin") == 0 )*/
    else {
      for (isite = StdI::NsiteUC*kCell; isite < StdI::NsiteUC*(kCell + 1); isite++) {
        StdFace::HubbardLocal(StdI::mu, -StdI::h, -StdI::Gamma, 0.0, isite);
      }
    }/*if (strcmp(StdI::model, "spin") != 0 )*/
    /*
     Hopping
    */
    for (it = 0; it < NtUJ[0]; it++) {
      /*
       Local term
      */
      if (tUJindx[0][it][0] == 0 && tUJindx[0][it][1] == 0 && tUJindx[0][it][2] == 0
        && tUJindx[0][it][3] == tUJindx[0][it][4])
      {
        if (strcmp(StdI::model, "hubbard") == 0) {
          isite = StdI::NsiteUC*kCell + tUJindx[0][it][3];
          for (ispin = 0; ispin < 2; ispin++) {
            StdI::trans.push_back(-tUJ[0][it]);
            StdI::transindx.push_back(std::vector<int>{isite, ispin, isite, ispin});
          }/*for (ispin = 0; ispin < 2; ispin++)*/
        }/*if (strcmp(StdI::model, "hubbrad") == 0 )*/
      }/*Local term*/
      else {
        /*
         Non-local term
        */
        StdFace::FindSite(iW, iL, iH,
          tUJindx[0][it][0], tUJindx[0][it][1], tUJindx[0][it][2],
          tUJindx[0][it][3], tUJindx[0][it][4], &isite, &jsite, &Cphase, dR);
        if (strcmp(StdI::model, "spin") == 0) {
          for (ii = 0; ii < 3; ii++) 
            Jtmp[ii][ii] = 2.0 * real(tUJ[0][it] * conj(tUJ[0][it]))
            * (1.0 / Uspin[tUJindx[0][it][3]] + 1.0 / Uspin[tUJindx[0][it][4]]);
          StdFace::GeneralJ(Jtmp, StdI::S2, StdI::S2, isite, jsite);
        }/*if (strcmp(StdI::model, "spin") == 0 )*/
        else {
          StdFace::Hopping(- Cphase * tUJ[0][it], jsite, isite, dR);
        }
      }/*Non-local term*/
    }/*for (it = 0; it < NtUJ[0]; it++)*/
    /*
    Coulomb integral (U)
    */
    for (it = 0; it < NtUJ[1]; it++) {
      /*
      Local term
      */
      if (tUJindx[1][it][0] == 0 && tUJindx[1][it][1] == 0 && tUJindx[1][it][2] == 0
        && tUJindx[1][it][3] == tUJindx[1][it][4])
      {
        StdI::Cintra.push_back(real(tUJ[1][it]));
        StdI::CintraIndx.push_back(StdI::NsiteUC*kCell + tUJindx[1][it][3]);
        /*
        Double-counting correction @f$0.5 U_{0ii} D_{0ii}@f$
        */
        if (StdI::double_counting == 1) {
          isite = StdI::NsiteUC*kCell + tUJindx[1][it][3];
          for (ispin = 0; ispin < 2; ispin++) {
            DenMat0 = DenMat[0][0][0][tUJindx[1][it][3]][tUJindx[1][it][3]];
            StdI::trans.push_back(0.5*real(tUJ[1][it])*DenMat0);
            StdI::transindx.push_back(std::vector<int>{isite, ispin, isite, ispin});
          }/*for (ispin = 0; ispin < 2; ispin++)*/
        }
      }/*Local term*/
      else {
        /*
        Non-local term
        */
        StdFace::FindSite(iW, iL, iH,
          tUJindx[1][it][0], tUJindx[1][it][1], tUJindx[1][it][2],
          tUJindx[1][it][3], tUJindx[1][it][4], &isite, &jsite, &Cphase, dR);
        StdFace::Coulomb(real(tUJ[1][it]), isite, jsite);
        /*
        Double-counting correction
        */
        if (StdI::double_counting == 1) {
          for (ispin = 0; ispin < 2; ispin++) {
            /*
            @f$sum_{(R,j)(>0,i)} U_{Rij} D_{0jj} (Local)@f$
            */
            DenMat0 = DenMat[0][0][0][tUJindx[1][it][4]][tUJindx[1][it][4]];
            StdI::trans.push_back(real(tUJ[1][it])*DenMat0);
            StdI::transindx.push_back(std::vector<int>{isite, ispin, isite, ispin});
            /*
            @f$sum_{(R,j)(>0,i)} U_{Rij} D_{0jj} (Local)@f$
            */
            DenMat0 = DenMat[0][0][0][tUJindx[1][it][3]][tUJindx[1][it][3]];
            StdI::trans.push_back(real(tUJ[1][it])*DenMat0);
            StdI::transindx.push_back(std::vector<int>{jsite, ispin, jsite, ispin});
          }/*for (ispin = 0; ispin < 2; ispin++)*/
          /*
          @f$-0.5U_{Rij} D_{Rjj}@f$
          */
          DenMat0 = DenMat[tUJindx[1][it][0]][tUJindx[1][it][1]][tUJindx[1][it][2]]
            [tUJindx[1][it][3]][tUJindx[1][it][4]];
          StdFace::Hopping(-0.5*Cphase * real(tUJ[1][it])*DenMat0, jsite, isite, dR);
        }/*if (StdI::double_counting == 1)*/
      }/*Non-local term*/
    }/*for (it = 0; it < NtUJ[0]; it++)*/
    /*
     Hund coupling (J)
    */
    for (it = 0; it < NtUJ[2]; it++) {
      /*
      Local term should not be computed
      */
      if (tUJindx[2][it][0] != 0 || tUJindx[2][it][1] != 0 || tUJindx[2][it][2] != 0
        || tUJindx[2][it][3] != tUJindx[2][it][4])
      {
        StdFace::FindSite(iW, iL, iH,
          tUJindx[2][it][0], tUJindx[2][it][1], tUJindx[2][it][2],
          tUJindx[2][it][3], tUJindx[2][it][4], &isite, &jsite, &Cphase, dR);

        StdI::Hund.push_back(real(tUJ[2][it]));
        StdI::HundIndx.push_back(std::vector<int>{isite, jsite});

        if (strcmp(StdI::model, "hubbard") == 0) {
          StdI::Ex.push_back(real(tUJ[2][it]));
          StdI::ExIndx.push_back(std::vector<int>{isite, jsite});

          StdI::PairHopp.push_back(real(tUJ[2][it]));
          StdI::PHIndx.push_back(std::vector<int>{isite, jsite});
          /*
          Double-counting correction
          */
          if (StdI::double_counting == 1) {
            for (ispin = 0; ispin < 2; ispin++) {
              /*
              @f$- \frac{1}{2}sum_{(R,j)(>0,i)} J_{Rij} D_{0jj}@f$
              */
              DenMat0 = DenMat[0][0][0][tUJindx[2][it][4]][tUJindx[2][it][4]];
              StdI::trans.push_back(-0.5*real(tUJ[2][it]) *DenMat0);
              StdI::transindx.push_back(std::vector<int>{isite, ispin, isite, ispin});
              /*
              @f$- \frac{1}{2}sum_{(R,j)(>0,i)} J_{Rij} D_{0jj}@f$
              */
              DenMat0 = DenMat[0][0][0][tUJindx[2][it][3]][tUJindx[2][it][3]];
              StdI::trans.push_back(-0.5*real(tUJ[2][it]) *DenMat0);
              StdI::transindx.push_back(std::vector<int>{jsite, ispin, jsite, ispin});
            }/*for (ispin = 0; ispin < 2; ispin++)*/
            /*
            @f$J_{Rij} (D_{Rjj}+2{\rm Re}[D_{Rjj])@f$
            */
            DenMat0 = DenMat[tUJindx[2][it][0]][tUJindx[2][it][1]][tUJindx[2][it][2]]
                            [tUJindx[2][it][3]][tUJindx[2][it][4]];
            StdFace::Hopping(
              0.5*Cphase * real(tUJ[2][it])*(DenMat0 + 2.0*real(DenMat0)), jsite, isite, dR);
          }/*if (StdI::double_counting == 1)*/
        }
        else {
#if defined(_mVMC)
          StdI::Ex.push_back(real(tUJ[2][it]));
#else
          StdI::Ex.push_back(-real(tUJ[2][it]));
#endif
          StdI::ExIndx.push_back(std::vector<int>{isite, jsite});
        }
      }/*Non-local term*/
    }/*for (it = 0; it < NtUJ[0]; it++)*/
  }/*for (kCell = 0; kCell < StdI::NCell; kCell++)*/

  fclose(fp);
  PrintUHFinitial(NtUJ, DenMat, tUJindx);
  StdFace::PrintXSF();
  StdFace::PrintGeometry();

  for (it = 0; it < NtUJ[0]; it++) free(tUJindx[0][it]);
  free(tUJindx[0]);
  free(tUJ[0]);
  for (it = 0; it < NtUJ[1]; it++) free(tUJindx[1][it]);
  free(tUJindx[1]);
  free(tUJ[1]); 
  for (it = 0; it < NtUJ[2]; it++) free(tUJindx[2][it]);
  free(tUJindx[2]);
  free(tUJ[2]); 
  if (strcmp(StdI::model, "spin") == 0) free(Uspin);

}/*void StdFace::Wannier90*/
