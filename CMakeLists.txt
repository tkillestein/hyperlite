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
option(HYPERLITE_HDF5 "Build with HDF5 (output.h5; OFF compiles a no-op stub)" ON)

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
find_package(LAPACK REQUIRED)  # dgetrf/dgetri in src/gas/physical_opacity_nlte.f
if(HYPERLITE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS Fortran)  # output.h5 (+ atomic data, Phase 3)
  set(HYPERLITE_HDF5_IO_SRC src/io/hdf5_io.f90)
else()
  set(HYPERLITE_HDF5_IO_SRC src/io/hdf5_io_stub.f90)
endif()

if(HYPERLITE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

if(HYPERLITE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  configure_file(src/core/mpimod_mpi.f ${PROJECT_BINARY_DIR}/mpimod.f COPYONLY)
else()
  configure_file(src/core/mpimod_ser.f ${PROJECT_BINARY_DIR}/mpimod.f COPYONLY)
endif()

#-- helper: propagate OpenMP/MPI usage to a target
function(hyperlite_target_deps tgt)
  if(HYPERLITE_HDF5)
    target_link_libraries(${tgt} PUBLIC HDF5::HDF5)
  endif()
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
  src/core/elemdatamod.f
  src/core/physconstmod.f
  src/core/randommod.f
  ${PROJECT_BINARY_DIR}/mpimod.f
  src/core/miscmod.f
  src/gas/ionsmod.f
  src/gas/ffxsmod.f
  src/gas/bfxsmod.f
  src/gas/bbxsmod.f
  src/gas/nltemod.f
  src/io/inputparmod.f
  src/group/groupmod.f
  src/grid/gridmod.f
  src/gas/gasmod.f
  src/io/inputstrmod.f
  src/transport/particlemod.f90
  src/core/timingmod.f
  src/core/countersmod.f
  src/output/fluxmod.f
  src/output/totalsmod.f
  src/transport/transportmod.f
  src/source/sourcemod.f
  ${HYPERLITE_HDF5_IO_SRC}
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
  src/core/binsrch.f90
  src/core/binsrch_decr.f90
  src/core/lcase.f
  src/core/locate.f
  src/core/memusg.f
  src/core/sorti.f
  src/core/specint.f
  src/core/warn.f
)

hyperlite_component(hyperlite_grid
  src/grid/grid_setup.f
  src/grid/grid_volume.f
  src/grid/emitgroup.f
  src/grid/emission_probability.f90
  src/grid/leakage_opacity.f
  src/grid/leakage_opacity11.f90
)

hyperlite_component(hyperlite_gas
  src/gas/analytic_opacity.f90
  src/gas/eos_update.f
  src/gas/gas_setup.f
  src/gas/gas_update.f
  src/gas/opacity_planckmean.f
  src/gas/opacity_rossmean.f
  src/gas/physical_opacity.f
  src/gas/physical_opacity_nlte.f
)

hyperlite_component(hyperlite_source
  src/source/boundary_source.f90
  src/source/interior_source.f90
  src/source/source_transformdirection.f90
  src/source/sourceenergy.f90
  src/source/sourcenumbers.f90
  src/source/vacancies.f90
)

hyperlite_component(hyperlite_transport
  src/transport/diffusion11.f90
  src/transport/transport11.f90
)

hyperlite_component(hyperlite_output
  src/output/open_logfiles.f
  src/output/output_flux.f
  src/output/output_grid.f
  src/output/output_grp.f
  src/output/output_profile.f
  src/output/output_source.f
  src/output/output_h5.f90
)

########################################################################
# executable
########################################################################
add_executable(superlite
  src/superlite.f90
  src/banner.f
  src/transport/particle_advance.f90
  src/output/fluxtally.f
  src/source/sourcetally.f
  src/dealloc_all.f
  src/gas/read_bbxs_data.f
  src/gas/read_nlte_data.f
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

if(Python3_FOUND AND HYPERLITE_HDF5)
  #-- bundled atomic data, generated from the raw ASCII (needs python3 h5py)
  add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/data/atomic.h5
    COMMAND Python3::Interpreter
            ${PROJECT_SOURCE_DIR}/tools/hyperlite_tools/build_atomic_h5.py
            -d ${PROJECT_SOURCE_DIR}/data
    DEPENDS ${PROJECT_SOURCE_DIR}/tools/hyperlite_tools/build_atomic_h5.py
            ${PROJECT_SOURCE_DIR}/data/data.ion
            ${PROJECT_SOURCE_DIR}/data/data.zeta
            ${PROJECT_SOURCE_DIR}/data/data.bf_verner
            ${PROJECT_SOURCE_DIR}/data/data.ff_sutherland
    COMMENT "Generating data/atomic.h5")
  add_custom_target(atomic_h5 ALL DEPENDS ${PROJECT_SOURCE_DIR}/data/atomic.h5)
endif()

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

  #-- statistical regression reading output.h5 directly + h5<->ascii roundtrip
  if(HYPERLITE_HDF5)
    add_test(NAME regression-smoke-h5
      COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/regression/compare.py
              ${SMOKE_RUNDIR} --source h5 -v)
    set_tests_properties(regression-smoke-h5 PROPERTIES
      FIXTURES_REQUIRED smoke)

    add_test(NAME regression-roundtrip
      COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/regression/roundtrip.py
              ${SMOKE_RUNDIR})
    set_tests_properties(regression-roundtrip PROPERTIES
      FIXTURES_REQUIRED smoke)
  endif()

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
  "(MPI=${HYPERLITE_MPI} OpenMP=${HYPERLITE_OPENMP} HDF5=${HYPERLITE_HDF5} ${CMAKE_BUILD_TYPE})")
