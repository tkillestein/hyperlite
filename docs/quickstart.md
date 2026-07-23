# Quickstart

## Prerequisites

- CMake ≥ 3.23 and a Fortran compiler (gfortran ≥ 10; ifx also supported)
- LAPACK/BLAS
- HDF5 with Fortran bindings (`libhdf5-dev` / `hdf5-devel` / brew `hdf5`)
- Python 3 with `numpy` and `h5py` (atomic-data generation + regression tests)
- Optional: OpenMPI for the MPI build

On Debian/Ubuntu:

```bash
sudo apt-get install gfortran cmake liblapack-dev libblas-dev libhdf5-dev
pip install numpy h5py
```

## Build

```bash
git clone https://github.com/tkillestein/hyperlite
cd hyperlite
cmake --preset gfortran-openmp       # or gfortran-serial / -mpi / -debug
cmake --build --preset gfortran-openmp -j
```

The build also generates `data/atomic.h5` from the bundled raw atomic data
(needs `python3` + `h5py`).

## Test

```bash
ctest --preset gfortran-openmp
```

This runs the unit tests, a full W7 smoke case (a typical SN Ia deck,
~1-3 min threaded), and the statistical regression checks against the
committed reference ensemble — see the
[developer guide](developer.md#tests) for what each test does.

## Run the smoke case by hand

A run directory needs the input deck (`input.par`), the structure file
(`input.str`), and the atomic data:

```bash
mkdir run && cd run
ln -s ../data/* .
ln -s ../tests/cases/w7/* .
ln -sf input_w7.str input.str
cp input.par.lte input.par
../build/gfortran-openmp/superlite
```

Outputs land in the run directory as `output.h5` (and legacy ASCII
`output.*` tables while `in_io_ascii = t`). See
[Running a simulation](running.md) for deck parameters and parallel runs.
