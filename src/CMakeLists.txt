# include guard
cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of HPhi.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
add_definitions(-D_HPhi)
add_definitions(-DDSFMT_MEXP=19937)

set(SOURCES FileIO.cpp HPhiMain.cpp HPhiTrans.cpp bitcalc.cpp check.cpp CheckMPI.cpp dSFMT.cpp diagonalcalc.cpp global.cpp log.cpp input.cpp output.cpp output_list.cpp readdef.cpp sz.cpp xsetmem.cpp wrapperMPI.cpp splash.cpp time.cpp)

set(SOURCES_MLTPLY Multiply.cpp mltply.cpp mltplySpin.cpp mltplyHubbard.cpp mltplyMPIHubbard.cpp mltplyMPISpin.cpp mltplyHubbardCore.cpp mltplySpinCore.cpp mltplyMPIHubbardCore.cpp mltplyMPISpinCore.cpp)
set(SOURCES_CG CalcByLOBPCG.cpp )
set(SOURCES_TPQ CalcByTPQ.cpp FirstMultiply.cpp)
set(SOURCES_DIAG CalcByFullDiag.cpp lapack_diag.cpp matrixlapack.cpp matrixlapack_magma.cpp matrixscalapack.cpp )
set(SOURCES_SPECTRUM CalcSpectrum.cpp CalcSpectrumByBiCG.cpp CalcSpectrumByFullDiag.cpp SingleEx.cpp SingleExHubbard.cpp PairEx.cpp PairExHubbard.cpp PairExSpin.cpp )
set(SOURCES_TEM CalcByTEM.cpp)
set(SOURCES_PHYS expec_cisajs.cpp expec_cisajscktaltdc.cpp expec_totalspin.cpp phys.cpp expec_energy_flct.cpp)

set(SOURCES_COMMON common/setmemory.cpp)

set(SOURCES_STDFACE StdFace/ChainLattice.cpp StdFace/FCOrtho.cpp StdFace/HoneycombLattice.cpp StdFace/Kagome.cpp StdFace/Ladder.cpp StdFace/Orthorhombic.cpp StdFace/Pyrochlore.cpp StdFace/SquareLattice.cpp StdFace/StdFace_main.cpp StdFace/StdFace_ModelUtil.cpp StdFace/TriangularLattice.cpp StdFace/Wannier90.cpp)

add_executable(HPhi++ ${SOURCES} ${SOURCES_CG} ${SOURCES_TPQ} ${SOURCES_SPECTRUM} ${SOURCES_TEM} ${SOURCES_MLTPLY} ${SOURCES_DIAG} ${SOURCES_PHYS} ${SOURCES_STDFACE} ${SOURCES_COMMON})
target_link_libraries(HPhi++ ${LAPACK_LIBRARIES} m)
if(MPI_FOUND)
  target_link_libraries(HPhi++ ${MPI_CXX_LIBRARIES})
endif(MPI_FOUND)

if(MAGMA_FOUND)
  target_link_libraries(HPhi++ ${MAGMA_CXX_LIBRARIES})
  add_definitions(-D_MAGMA)
endif(MAGMA_FOUND)
if(CUDA_FOUND)
  target_link_libraries(HPhi++ ${CUDA_LIBRARIES})
endif(CUDA_FOUND)

if(USE_SCALAPACK)
  string(REGEX REPLACE "-L[ ]+" "-L" sc_libs "${SCALAPACK_LIBRARIES}")
  string(REGEX REPLACE "[ ]+" ";" sc_libs "${sc_libs}")
  foreach(sc_lib IN LISTS sc_libs)
    target_link_libraries(HPhi++ ${sc_lib})
  endforeach(sc_lib)
  message(STATUS "SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")
endif(USE_SCALAPACK)


install(TARGETS HPhi++ RUNTIME DESTINATION bin)
