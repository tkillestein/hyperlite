# Building

hyperlite uses modern CMake (≥ 3.23 for presets) with CTest. All builds go
through presets defined in `CMakePresets.json`.

## Presets

| Preset | Compiler | MPI | OpenMP | HDF5 | Notes |
|---|---|---|---|---|---|
| `gfortran-serial` | gfortran | – | – | ✓ | reference configuration for the regression ensemble |
| `gfortran-openmp` | gfortran | – | ✓ | ✓ | recommended for production single-node runs |
| `gfortran-mpi` | gfortran | ✓ | ✓ | ✓ | hybrid MPI+OpenMP |
| `gfortran-debug` | gfortran | – | ✓ | ✓ | `-g -O0 -fcheck=all -fbacktrace` |
| `ifx-release` | Intel ifx | – | ✓ | – (stub) | CI-only, non-blocking |

```bash
cmake --preset gfortran-openmp
cmake --build --preset gfortran-openmp -j
ctest  --preset gfortran-openmp
```

Build trees land in `build/<preset>/`; the executable is
`build/<preset>/superlite`.

## Options

| CMake option | Default | Effect |
|---|---|---|
| `HYPERLITE_MPI` | `OFF` | MPI build (cell-decomposed gas solve, replicated transport) |
| `HYPERLITE_OPENMP` | `ON` | OpenMP threading of the particle loop and opacity build |
| `HYPERLITE_HDF5` | `ON` | HDF5 I/O; `OFF` compiles a no-op stub (no `output.h5`, ASCII atomic data only) |

## Dependencies

- **LAPACK/BLAS** — required (NLTE opacity solve).
- **HDF5 (Fortran)** — required unless `HYPERLITE_HDF5=OFF`.
- **Python 3 + h5py** — build-time generation of `data/atomic.h5`
  (target `atomic_h5`); also the regression tests (plus numpy).
- **MPI** — only for `HYPERLITE_MPI=ON`.

## Compiler flags

The tree is held to `-std=f2018 -Wall` on gfortran with zero warnings.
The only exemption is the generated MPI `mpimod` (built with
`-std=legacy -fallow-argument-mismatch`, inherent to the pre-`mpi_f08`
MPI interface).

## Version stamping

`cmake/GenerateVersion.cmake` regenerates `versionmod.f90` at build time
from the project version and `git describe --tags --always --dirty`. The
stamp is printed in the run banner and written to `/meta` in `output.h5`.

## Install

```bash
cmake --install build/gfortran-openmp --prefix ~/.local
```

installs the `superlite` binary to `<prefix>/bin`.
