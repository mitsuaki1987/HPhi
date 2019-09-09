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

/*!
  @mainpage

  <H2>Introduction</H2>
  A numerical solver package for a wide range of quantum lattice models including Hubbard-type itinerant electron hamiltonians, quantum spin models, and Kondo-type hamiltonians for itinerant electrons coupled with quantum spins. The Lanczos algorithm for finding ground states and newly developed Lanczos-based algorithm for finite-temperature properties of these models are implemented for parallel computing. A broad spectrum of users including experimental researchers is cordially welcome.
  <HR>
  <H2>Developers</H2>
  Youhei Yamaji (Quantum-Phase Electronics Center, The University of Tokyo)\n
  Takahiro Misawa (Institute for Solid State Physics, The University of Tokyo)\n
  Synge Todo (Department of Physics, The University of Tokyo)\n
  Kazuyoshi Yoshimi (Institute for Solid State Physics, The University of Tokyo)\n
  Mitsuaki Kawamura (Institute for Solid State Physics, The University of Tokyo)\n
  Kota Ido (Department of Applied Physics, The University of Tokyo)\n
  Naoki Kawashima (Institute for Solid State Physics, The University of Tokyo)
  <HR>
  <H2>Methods</H2>
  - Lanczos algorithm
  - Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG) method : See LOBPCG_Main()
  - Thermal Pure Quantum (TPQ) state
  - Full Diagonalization
  - Shifted BiCG method : See CalcSpectrumByBiCG()
  - Lehmann's spectral representation : See CalcSpectrumByFullDiag()
  .
  <HR>
  <H2>Target models</H2>
  Hubbard model, Heisenberg model, Kondo lattice model, Kitaev model, Kitaev-Heisenberg model, multi-orbital Hubbard model
  <HR>
  <H2>Important functions and source files</H2>
  <ul>
    <li>mltply.cpp : Perform Hamiltonian-vector product</li>
    <ul>
      <li>mltplyHubbard.cpp : For Hubbard and Kondo system</li>
      <li>mltplySpin.cpp : For local spin system</li>
    </ul>
    <li>StdFace_main.cpp : Construct typical models</li>
    <li>global.h : Global variables</li>
  </ul>
  <HR>
  <H2>How to modify HPhi (Developer's note)</H2>
  - @ref page_codingrule
  - @ref page_addstandard
  - @ref page_addstandardval
  - @ref page_variable
  - @ref page_setmem
  - @ref page_cmake
  - @ref page_addcalcmod
  - @ref page_addmodpara
  - @ref page_addexpert
  - @ref page_time
  - @ref page_log
  - Some contrivances for HPhi (only in Jananese) http://qlms.github.io/HPhi/develop/tips.pdf
  .
  <HR>
  <H2>Link</H2>
  https://github.com/QLMS/HPhi
  <HR>
  <H2>Download</H2>
  https://github.com/QLMS/HPhi/releases
  <HR>
  <H2>Forum</H2>
  https://github.com/QLMS/HPhi/issues
  <HR>
  <H2>licence</H2>
  <B>GNU GPL version 3</B>\n
  This software is developed under the support of "Project for advancement of software usability in materials science" by The Institute for Solid State Physics, The University of Tokyo.\n

@page page_codingrule Coding rule

- Do not use TAB character. Use two spaces as an indent.
- Use @c default(none) for scoping of OpenMP-parallel region. E.g.
  @dontinclude CalcByLOBPCG.cpp
  @skip pragma
  @until 0.0
- Variable declared with @c const must not be included in @c firstprivate of OpenMP scoping.
  Use @c shared.
- For MPI parallelization, use the following functions for I/O and abortation:
  - wrapperMPI::Fgets() instead of @c fgets
  - @c fprintf(::MP::STDOUT,... instead of @c printf(...
  - wrapperMPI::Fopen() instead of @c fopen
  - wrapperMPI::Exit() instead of @c exit
- When you add new features into HPhi, please run <tt>make test</tt>,
  and check whether other features still work fine.
  Also, try <tt>make test MPIRUN="mpiexec -np 4"</tt> to check MPI feature.
.    

@page page_cmake Add new source-file, executable, scripts (handle CMake)

To build HPhi, CMake is required.
We have to modify the CMake configuration file when we add new sources, executables, scripts.

@section sec_newsource New source code

When we add a new source code, we have to add the file-name
into the following part of @c src/CMakeLists.txt.

\code{cmake}
set(SOURCES source1.cpp source2.cpp ...)
\endcode

@section sec_newexecutable New executable

When we add a new executable ("myprog" in this case),
we have to add following command in @c src/CMakeLists.txt.

\code{CMake}
set(SOURCES_MYPROG source1.cpp source2.cpp ...)
add_executable(myprog ${SOURCES_MYPROG})
target_link_libraries(myprog ${LAPACK_LIBRARIES} m)
if(MPI_FOUND)
target_link_libraries(myprog ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
install(TARGETS myprog RUNTIME DESTINATION bin)
\endcode

@section sec_newscript New script

When we add a new script written in python, sh, etc. ("myscript.sh" in this case)
into @c tool/, we have to add the following command in @c tool/CMakeLists.txt.

\code{CMake}
configure_file(myscript.sh myscript.sh COPYONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/myscript.sh DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
\endcode
*/

#include "sz.hpp"
#include "HPhiTrans.hpp"
#include "output_list.hpp"
#include "diagonalcalc.hpp"
#include "CalcByLOBPCG.hpp"
#include "CalcByFullDiag.hpp"
#include "CalcByTPQ.hpp"
#include "CalcSpectrum.hpp"
#include "check.hpp"
#include "CalcByTEM.hpp"
#include "readdef.hpp"
#include "StdFace/StdFace_main.hpp"
#include "wrapperMPI.hpp"
#include "splash.hpp"
#include "CalcTime.hpp"
#include "common/setmemory.hpp"
#include "xsetmem.hpp"
#include "log.hpp"
#include "DefCommon.hpp"
#include "global.hpp"
#include <sys/stat.h>
#include <cstdlib>

/** 
 * @brief Main program for HPhi
 * 
 * @param argc [in] argument count
 * @param argv [in] argument vector
 *
 * @version 2.1 Add Time evolution mode.
 * @version 1.2 Add calculation spectrum mode.
 * @version 1.0
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @retval -1 fail the calculation.
 * @retval 0 succeed the calculation.
 */
int main(int argc, char* argv[]){

  int mode=0;
  char cFileListName[D_FileNameMax];

  MP::STDOUT = stdout;
  JudgeDefType(argc, argv, &mode);

  if (mode == STANDARD_DRY_MODE) {
    MP::myrank = 0;
    MP::nproc = 1;
    MP::STDOUT = stdout;
    splash();
  }
  else wrapperMPI::Initialize(argc, argv);

  //Timer
  InitTimer();
  if (mode != STANDARD_DRY_MODE) StartTimer(0);
 
  //MakeDirectory for output
  struct stat tmpst;
  if (MP::myrank == 0) {
    if (stat("./output/", &tmpst) != 0) {
      if (mkdir("./output/", 0777) != 0) {
        fprintf(MP::STDOUT, "%s", "Error: Fail to make output folder in current directory. \n");
        wrapperMPI::Exit(-1);
      }
    }
  }/*if (MP::myrank == 0)*/

  strcpy(cFileListName, argv[2]);
  
  if(mode==STANDARD_MODE || mode == STANDARD_DRY_MODE){
    if (MP::myrank == 0) StdFace::main(argv[2]);
    strcpy(cFileListName, "namelist.def");
    if (mode == STANDARD_DRY_MODE){
      fprintf(stdout, "Dry run is Finished. \n\n");
      return 0;
    }
  }

  xsetmem::HEAD();
  if(ReadDefFileNInt(cFileListName)!=0){
    fprintf(MP::STDOUT, "%s", "Error: Definition files(*.def) are incomplete.\n");
    wrapperMPI::Exit(-1);
  }

  if (Def::nvec < Def::k_exct){
    fprintf(MP::STDOUT, "%s", "Error: nvec should be smaller than exct are incorrect.\n");
    fprintf(MP::STDOUT, "Error: nvec = %d, exct=%d.\n", Def::nvec, Def::k_exct);
    wrapperMPI::Exit(-1);
  }
  fprintf(MP::STDOUT, "%s", "\n######  Definition files are correct.  ######\n\n");
  
  /*ALLOCATE-------------------------------------------*/
  xsetmem::def();
  /*-----------------------------------------------------*/

  /*Read Def files.*/
  TimeKeeper("%s_TimeKeeper.dat", "Read File starts:   %s", "w");
  if (ReadDefFileIdxPara() != 0) {
    fprintf(MP::STDOUT,
            "Error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
    wrapperMPI::Exit(-1);
  }
  TimeKeeper("%s_TimeKeeper.dat", "Read File finishes: %s", "a");
  fprintf(MP::STDOUT, "%s", "\n######  Indices and Parameters of Definition files(*.def) are complete.  ######\n\n");

  /*Set convergence Factor*/
  SetConvergenceFactor();

  /*---------------------------*/
  if (HPhiTrans() != 0) {
    wrapperMPI::Exit(-1);
  }

  //Start Calculation
  if (Def::iFlgCalcSpec == DC::CALCSPEC_NOT || Def::iFlgCalcSpec == DC::CALCSPEC_SCRATCH) {
    
    if (check() == MPIFALSE) {
     wrapperMPI::Exit(-1);
    }
    
    /*LARGE VECTORS ARE ALLOCATED*/
    xsetmem::large();
    
    StartTimer(1000);
    if(sz(List::c1, List::c2_1, List::c2_2)!=0){
      wrapperMPI::Exit(-1);
    }

    StopTimer(1000);
    if(Def::WRITE==1){
      output_list();
      wrapperMPI::Exit(-2);
    }
    StartTimer(2000);
    diagonalcalc();
    StopTimer(2000);
      
    switch (Def::iCalcType) {
    case DC::CG:
      if (CalcByLOBPCG::main() != TRUE) {
          wrapperMPI::Exit(-3);
      }
      break;

    case DC::FullDiag:
      StartTimer(5000);
      if (Def::iFlgScaLAPACK == 0 && MP::nproc != 1) {
        fprintf(MP::STDOUT, "Error: Full Diagonalization by LAPACK is only allowed for one process.\n");
        wrapperMPI::Finalize();
      }
      if (CalcByFullDiag::main() != TRUE) {
        wrapperMPI::Finalize();
      }
      StopTimer(5000);
      break;

    case DC::TPQCalc:
      StartTimer(3000);
      if (CalcByTPQ(Step::NumAve, Param::ExpecInterval) != TRUE) {
        StopTimer(3000);
        wrapperMPI::Exit(-3);
      }
      StopTimer(3000);
      break;

    case DC::TimeEvolution:
      if (CalcByTEM(Param::ExpecInterval) != 0) {
        wrapperMPI::Exit(-3);
      }
      break;

    default:
      StopTimer(0);
      wrapperMPI::Exit(-3);
    }
  }

  if(Def::iFlgCalcSpec != DC::CALCSPEC_NOT){
    StartTimer(6000);
    CalcSpectrum();
    StopTimer(6000);
  }
  
  StopTimer(0);
  OutputTimer();
  wrapperMPI::Finalize();
  return 0;
}
