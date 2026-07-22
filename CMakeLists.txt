#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
cmake_minimum_required(VERSION 3.20)

project(hyperlite
  VERSION 1.99.0
  DESCRIPTION "SuperLite/hyperlite IMC-DDMC radiation transport for interacting transients"
  LANGUAGES Fortran)

#-- options
option(HYPERLITE_MPI "Build with MPI (cell-decomposed gas solve, replicated transport)" OFF)
option(HYPERLITE_OPENMP "Build with OpenMP threading" ON)

#-- default to a Release build
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

#-- all Fortran modules land in one directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

#-- compiler flags: the code base is legacy Fortran (real*8, argument
#-- mismatches in MPI calls); tightened per-phase during modernization
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-std=legacy -fallow-argument-mismatch)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -fcheck=all -fbacktrace")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback")
endif()

#-- dependencies
find_package(LAPACK REQUIRED)  # dgetrf/dgetri in GAS/physical_opacity_nlte.f

if(HYPERLITE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

if(HYPERLITE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  configure_file(mpimod_mpi.f ${PROJECT_BINARY_DIR}/mpimod.f COPYONLY)
else()
  configure_file(mpimod_ser.f ${PROJECT_BINARY_DIR}/mpimod.f COPYONLY)
endif()

#-- helper: propagate OpenMP/MPI usage to a target
function(hyperlite_target_deps tgt)
  if(HYPERLITE_OPENMP)
    target_link_libraries(${tgt} PUBLIC OpenMP::OpenMP_Fortran)
  endif()
  if(HYPERLITE_MPI)
    target_link_libraries(${tgt} PUBLIC MPI::MPI_Fortran)
  endif()
endfunction()

########################################################################
# module library (all .mod-producing sources; consumed by every component)
########################################################################
add_library(hyperlite_modules OBJECT
  elemdatamod.f
  physconstmod.f
  randommod.f
  ${PROJECT_BINARY_DIR}/mpimod.f
  miscmod.f
  ionsmod.f
  ffxsmod.f
  bfxsmod.f
  bbxsmod.f
  nltemod.f
  inputparmod.f
  groupmod.f
  gridmod.f
  gasmod.f
  inputstrmod.f
  particlemod.f90
  timingmod.f
  countersmod.f
  fluxmod.f
  totalsmod.f
  transportmod.f
  sourcemod.f
)
hyperlite_target_deps(hyperlite_modules)

########################################################################
# component libraries (former per-directory archive libs)
########################################################################
function(hyperlite_component name)
  add_library(${name} OBJECT ${ARGN})
  target_link_libraries(${name} PRIVATE hyperlite_modules)
  hyperlite_target_deps(${name})
endfunction()

hyperlite_component(hyperlite_misc
  MISC/binsrch.f90
  MISC/binsrch_decr.f90
  MISC/lcase.f
  MISC/locate.f
  MISC/memusg.f
  MISC/sorti.f
  MISC/specint.f
  MISC/warn.f
)

hyperlite_component(hyperlite_grid
  GRID/grid_setup.f
  GRID/grid_volume.f
  GRID/emitgroup.f
  GRID/emission_probability.f90
  GRID/leakage_opacity.f
  GRID/leakage_opacity11.f90
)

hyperlite_component(hyperlite_gas
  GAS/analytic_opacity.f90
  GAS/eos_update.f
  GAS/gas_setup.f
  GAS/gas_update.f
  GAS/opacity_planckmean.f
  GAS/opacity_rossmean.f
  GAS/physical_opacity.f
  GAS/physical_opacity_nlte.f
)

hyperlite_component(hyperlite_source
  SOURCE/boundary_source.f90
  SOURCE/interior_source.f90
  SOURCE/source_transformdirection.f90
  SOURCE/sourceenergy.f90
  SOURCE/sourcenumbers.f90
  SOURCE/vacancies.f90
)

hyperlite_component(hyperlite_transport
  TRANSPORT1/diffusion11.f90
  TRANSPORT1/transport11.f90
)

hyperlite_component(hyperlite_output
  OUTPUT/open_logfiles.f
  OUTPUT/output_flux.f
  OUTPUT/output_grid.f
  OUTPUT/output_grp.f
  OUTPUT/output_profile.f
  OUTPUT/output_source.f
)

########################################################################
# executable
########################################################################
add_executable(superlite
  superlite.f90
  banner.f
  particle_advance.f90
  fluxtally.f
  sourcetally.f
  dealloc_all.f
  read_bbxs_data.f
  read_nlte_data.f
)
target_link_libraries(superlite PRIVATE
  hyperlite_modules
  hyperlite_misc
  hyperlite_grid
  hyperlite_gas
  hyperlite_source
  hyperlite_transport
  hyperlite_output
  ${LAPACK_LIBRARIES}
)
hyperlite_target_deps(superlite)

install(TARGETS superlite RUNTIME DESTINATION bin)

########################################################################
# tests (CTest wrappers around the Phase-0 statistical harness)
########################################################################
enable_testing()
find_package(Python3 COMPONENTS Interpreter)

if(Python3_FOUND)
  set(SMOKE_RUNDIR ${PROJECT_BINARY_DIR}/smoke_run)

  #-- stage + run the W7 LTE smoke case (fixture for the checks below)
  add_test(NAME smoke-run
    COMMAND ${PROJECT_SOURCE_DIR}/tests/smoke_run.sh
            $<TARGET_FILE:superlite> ${PROJECT_SOURCE_DIR} ${SMOKE_RUNDIR})
  set_tests_properties(smoke-run PROPERTIES
    FIXTURES_SETUP smoke TIMEOUT 1800)

  #-- statistical regression vs the committed reference ensemble
  add_test(NAME regression-smoke
    COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/regression/compare.py
            ${SMOKE_RUNDIR} -v)
  set_tests_properties(regression-smoke PROPERTIES
    FIXTURES_REQUIRED smoke)

  #-- harness self-test: a perturbed run MUST fail the comparison
  add_test(NAME regression-selftest
    COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/regression/compare.py
            ${SMOKE_RUNDIR} --perturb 0.05)
  set_tests_properties(regression-selftest PROPERTIES
    FIXTURES_REQUIRED smoke)

  #-- MPI completion check (regression stays pinned to serial: see OVERHAUL.md
  #-- known defect on nmpi-dependent eout/evelo bias)
  if(HYPERLITE_MPI)
    add_test(NAME smoke-run-mpi
      COMMAND ${PROJECT_SOURCE_DIR}/tests/smoke_run.sh
              $<TARGET_FILE:superlite> ${PROJECT_SOURCE_DIR}
              ${PROJECT_BINARY_DIR}/smoke_run_mpi
              ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2)
    set_tests_properties(smoke-run-mpi PROPERTIES TIMEOUT 1800)
  endif()
else()
  message(STATUS "Python3 not found: regression tests disabled")
endif()

message(STATUS "hyperlite ${PROJECT_VERSION} configured "
  "(MPI=${HYPERLITE_MPI} OpenMP=${HYPERLITE_OPENMP} ${CMAKE_BUILD_TYPE})")
