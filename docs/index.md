# hyperlite

**hyperlite** is a modernized fork of
[SuperLite](https://github.com/gururajw/superlite), a 1-D spherical,
multi-group **Implicit Monte Carlo / Discrete Diffusion Monte Carlo**
(IMC-DDMC) radiation-transport code for spectral synthesis of interacting
transients (supernovae), descended from the LANL/UW-Madison
[SuperNu](https://github.com/lanl/supernu) code.

SuperLite is described in
[Wagle et al. 2023, ApJ 953 132](https://doi.org/10.3847/1538-4357/acda23);
please cite it (and the Zenodo
record [10.5281/zenodo.8102236](https://doi.org/10.5281/zenodo.8102236)) if
you use this code.

## What the fork changes

The physics is SuperLite's; the engineering around it is rebuilt
(see the [changelog](changelog.md) and `OVERHAUL.md` in the repository for
the full plan and history):

- **Modern build**: CMake presets + CTest, gfortran and ifx, CI matrix
  (serial / OpenMP / MPI). The recursive-Make build is gone.
- **Fast**: SoA particle storage, per-cell DDMC tables, OpenMP on by
  default — **5.4× faster** on the reference case than the original
  serial code (1.5× serial, the rest from threading).
- **Reproducible**: a counter-based Philox-4x32-10 RNG keyed per particle
  makes results independent of thread count and scheduling; a statistical
  regression harness with a committed 20-seed reference ensemble guards
  every change.
- **Self-describing I/O**: a single HDF5 `output.h5` (plus optional legacy
  ASCII tables), and atomic data bundled as `data/atomic.h5` — no more
  `Atoms/` working-directory dependency.
- **Correctness fixes**: an LTE partition-function caching bug that biased
  opacities (and made results depend on the MPI rank count) is fixed —
  a genuine physics correction (+15% total luminosity on the reference
  case); a stale DDMC Doppler-shift cache is fixed.
- **Packaged tools**: the Python helpers ship as the `hyperlite-tools`
  package with rich-click CLIs.
- **Modern Fortran**: the whole tree builds under `-std=f2018 -Wall` with
  zero warnings.

## Where to go

- [Quickstart](quickstart.md) — build and run the smoke case in minutes.
- [Building](building.md) — presets, options, dependencies.
- [Running a simulation](running.md) — input decks, staging, parallelism.
- [I/O schema](io-schema.md) — the `output.h5` and `atomic.h5` layouts.
- [Developer guide](developer.md) — tests, benchmarks, RNG design.

## License

GPLv3. © 2023 Gururaj A. Wagle (SuperLite); portions © LANL (SuperNu
lineage; see `LANL_README`).
