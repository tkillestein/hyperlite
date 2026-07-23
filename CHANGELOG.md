All notable changes to hyperlite are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/); versions follow
[SemVer](https://semver.org/). hyperlite forked from
[SuperLite v1.0.1](https://github.com/gururajw/superlite).

## [2.0.0-rc] — 2026-07

The modernization overhaul (`OVERHAUL.md`, Phases 0–7). Physics
algorithms are SuperLite's; engineering is rebuilt end to end.

### Fixed (physics — results change beyond MC noise)

- **LTE partition-function caching**: the Saha EOS froze atomic partition
  functions at the first cell's temperature seen by each process, biasing
  ionization and every opacity downstream — and making results depend on
  the MPI rank count. Partition functions are now evaluated per cell.
  Total luminosity on the W7 reference case shifts **+15%**
  (2.72→3.13×10⁴² erg/s); opacities are byte-identical across rank counts.
- **Stale DDMC Doppler-shift cache** in `diffusion11`: out-of-lump
  wavelengths read a misaligned `specarr` cache. Fix zeroes two spurious
  far-UV edge bins and shifts 2100–2700 Å by 6–35% (L_tot conserved to
  2×10⁻⁵).
- **MPI vacancy ledger** in `sourcenumbers` (benign but wrong) and a stale
  buffer length in `allgather_leakage` (UB).
- LTE runs no longer dump uninitialized NLTE-only fields
  (`gas_capemit`/`capemitgrey`).

### Added

- **CMake + CTest build** with presets (`gfortran-serial/-openmp/-mpi/-debug`,
  `ifx-release`); GitHub Actions CI matrix. Legacy recursive Make removed.
- **Statistical regression harness** (`tests/regression/`): 20-seed
  committed reference ensemble, deterministic anchors + k·σ stochastic
  checks + conservation checks + perturbation self-test.
- **Unit tests** (`tests/unit/`): Philox known-answer vectors, `binsrch`,
  `specint` Planck integrals.
- **HDF5 I/O**: self-describing `output.h5` (schema v1); atomic data
  bundled as build-time-generated `data/atomic.h5` (`in_io_atomdata`),
  removing the `Atoms/` CWD dependency; `h5-to-ascii` exporter;
  `in_io_ascii` keeps the legacy tables.
- **Counter-based RNG**: Philox-4x32-10 keyed per (rank/seed, iteration,
  particle) — results independent of OpenMP thread count and scheduling.
  New `in_rnd_seed` parameter.
- **Performance**: SoA particle storage; per-cell DDMC tables replacing
  per-thread cache rebuilds; guided scheduling; OpenMP on by default
  (`in_nomp = 0` = all threads). **5.4×** total on the W7 case vs the
  1.x serial baseline (1.5× serial, 4.2× at 4 threads).
- **Benchmark harness** `hyperlite-bench`, reported in CI.
- **`hyperlite-tools` Python package** (uv/ruff/ty; rich-click CLIs):
  `stella2superlite`, `parse-hdf5-sndata`, `sgfilter`, `wlgenadd`,
  `h5-to-ascii`, `build-atomic-h5`, `hyperlite-bench`. `sgfilter` and
  `wlgenadd` ported Python 2 → 3.
- **Docs site** (mkdocs-material, built in CI): quickstart, build/run,
  physics overview, I/O schema, migration notes, developer guide.
- Build-time **version stamping** (git describe) in the banner and
  `output.h5` `/meta`.

### Changed

- Sources reorganized under `src/{core,io,grid,gas,group,source,transport,output}`;
  decks under `tests/cases/`; Python tools under `tools/`.
- Whole tree modernized to build under `-std=f2018 -Wall` with zero
  warnings (kind parameters, no statement functions/forall/extensions).
- Default OpenMP threading on (`in_nomp = 0`); input decks pinning
  `in_nomp = 1` behave as before.

### Removed

- Legacy Make machinery, dead `TRANSPORT2/3` references, `g95`/`xlf_bgp`
  toolchains, `Data/Tinit.dat` (uninitialized garbage), the mzran RNG,
  and the `< 1e-99` ASCII exponent-truncation hack (HDF5 path).

## [1.0.1] — 2023

SuperLite v1.0.1 (Wagle et al. 2023, ApJ 953 132), the fork point.
