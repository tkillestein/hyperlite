# SuperLite / hyperlite — Modernization Overhaul

> Architecture and execution plan for modernizing the SuperLite IMC–DDMC Monte
> Carlo radiation-transport code: toolchain, performance, I/O, packaging,
> documentation, and a regression-anchored test suite.

---

## 1. Context

SuperLite (this fork: **hyperlite**) is a 1‑D spherical, multi-group
**Implicit Monte Carlo / Discrete Diffusion Monte Carlo** (IMC–DDMC)
radiation-transport code for interacting transients (supernovae), descended from
the LANL/UW‑Madison **SuperNu** code. It is ~11.7k lines of Fortran (51
fixed-form `.f`, 16 free-form `.f90`), GPLv3, © 2023 Gururaj A. Wagle
(Wagle et al. 2023, ApJ 953 132; Zenodo 10.5281/zenodo.8102236).

The code works, but the surrounding engineering has aged badly and now blocks
both correctness confidence and performance:

- **Build rot.** The recursive Make build references directories that no longer
  exist (`TRANSPORT2/`, `TRANSPORT3/`) so a stock `make all` fails; the CMake
  build is a crude flat-glob stub; two of four compiler toolchains (`g95`,
  `xlf_bgp`) are dead; there is no CI and no packaging metadata.
- **No safety net.** `make check` points at a `Testsuite/first/test.sh` that
  does not exist; **no golden/reference outputs are committed anywhere**. There
  is nothing today that would catch a regression introduced by modernization.
- **Single-threaded by default.** OpenMP exists (77 directives) but ships
  disabled (`in_nomp = 1` in every input deck); the hot particle loop is an
  AoS layout that blocks SIMD; the RNG is not counter-based.
- **Fragile, non-standard I/O.** ~15 hand-rolled ASCII output files (with a
  precision-truncation hack for values `< 1e-99`), and ~27 MB of atomic data in
  brittle fixed-format ASCII read via a hard-coded `Atoms/` CWD dependency.

**Goal:** modernize the toolchain, dramatically speed up a *single* simulation
via vectorization + threading, move I/O to a self-describing standard format
(HDF5), reorganize/package/document the repo, and — first and underpinning
everything — build a regression test that proves the outputs stay physically
consistent through every change.

### Design decisions (locked)

| Area | Decision | Consequence |
|---|---|---|
| **Regression basis** | **Statistical** (within Monte Carlo noise) | Free to adopt a counter-based RNG immediately; reference is an *ensemble distribution*, not a single byte-exact file. |
| **Build system** | **Modern CMake + CTest** | Native Fortran module-dep tracking, presets, install rules, HDF5/MPI/OpenMP options; retire the bash dep-scanner and per-dir Makefiles. |
| **I/O format** | **HDF5 primary** + thin ASCII/CSV export | Single self-describing output + atomic-data bundle; matches astro convention (TARDIS, Sedona, ARTIS); fixes the precision hack and the `Atoms/` CWD dependency. |
| **Performance scope** | **SoA + threaded kernel** (stay in modern Fortran, keep hybrid MPI) | Particles AoS→SoA, counter-based RNG, OpenMP-by-default, batched/SIMD-friendly kernel, log-index math. A **GPU transport kernel is a stretch goal** (Phase 8) — the CPU refactor is deliberately built to enable it. |

---

## 2. Current-state assessment (facts that drive the plan)

**Hot path.** Essentially all compute is in one region:
`superlite.f90` main iteration loop → `particle_advance.f90` (loop over
`prt_npartmax` particles) → per-particle `do while` event loop dispatching, via a
**procedure pointer**, to `TRANSPORT1/transport11.f90` (IMC) or
`TRANSPORT1/diffusion11.f90` (DDMC). Three nested loop levels: particles (outer),
events-per-particle (middle, data-dependent/unbounded), frequency groups (inner,
inside DDMC). Secondary hotspot: the OpenMP opacity build in
`GAS/physical_opacity*.f`.

**Vectorization blockers** (`particle_advance.f90`, `particlemod.f90`):
1. **AoS particle storage** — `type(packet)` with 9 doubles;
   `ptcl = prt_particles(ipart)` gathers a strided struct. #1 blocker.
2. Data-dependent branchy control flow (`goto`, inner Newton–Raphson `do while`
   for the Doppler distance), the classic MC divergence problem.
3. **mzran RNG** (`randommod.f`) — scalar LCG-style recurrence, `pure`, per-thread
   stream via skip-ahead; **not counter-based**, so it cannot produce per-SIMD-lane
   draws.
4. Indirect addressing / gathers into `grd_cap(ig,ic)`, `grd_sig(ic)`,
   `grd_icell(...)`, and scatter-adds into `grd_tally`, `grd_jrad`.
5. Three data-dependent `binsrch` (`MISC/binsrch.f90`) calls at loop entry on
   *geometric* grids (`grp_wl` is log-spaced → replaceable by `log`-index math).

**Parallelism today.** MPI decomposes **cells** for the gas/opacity solve but
**replicates particles** for transport (each rank transports
`2**in_src_n2s / nmpi` packets with a disjoint RNG segment), combining tallies
with `mpi_reduce/allreduce/scatterv`. OpenMP parallelizes the particle loop with
`reduction(+:grd_tally,grd_jrad,...)` and `schedule(static,1)`, but ships
disabled. Grid/gas/group state is already **struct-of-arrays**; only the
particles are AoS. Note mixed precision: `grd_cap/capemit/emiss` are `real*4`.

**Reproducibility reality.** Output is deterministic only for fixed
`(nmpi, in_nomp)`; any reorder/rescheduling or an RNG swap changes the draw
order and breaks bit-identity. This is *why* we chose a statistical regression
basis — it lets the speedup work proceed.

**I/O seams** (all cleanly isolated, ideal wrap points):
- Inputs: `read_inputpars` (namelist, `inputparmod.f:288`), `read_inputstr`
  (`inputstrmod.f:85`, label→column table), `read_wlgrid` (`groupmod.f:55`),
  `read_fluxgrid` (`fluxmod.f:68`).
- Atomic data: `ions_read_data` (`ionsmod.f:652`), `bfxsmod.f`, `ffxsmod.f`,
  `bbxsmod.f` + `read_bbxs_data.f` / `read_nlte_data.f`. Expects `Atoms/` in CWD.
- Outputs: six self-contained `OUTPUT/output_*.f` writers, each already grouped
  by physical quantity.

**Known defects to sweep up early** (surfaced during exploration):
- Dead lib refs `TRANSPORT2/transport2.a`, `TRANSPORT3/transport3.a` break
  `make all`.
- `Testsuite/first/test.sh` missing → `make check` broken; **no reference
  outputs committed** despite `Testsuite/README` claiming otherwise.
- `OUTPUT/Makefile` typo `output_source.0` (should be `.o`).
- `GRID/grid_volume.f` on disk but orphaned from `GRID/Makefile`.
- Compiler-file usage message names a nonexistent
  `Makefile.compiler.intel-x86_64`.
- A4 test decks appear dimension-mismatched (`in_ndim=117` vs `.str` header
  `122`) — verify before using as regression inputs.
- `Data/Tinit.dat` is dumped uninitialized memory, referenced nowhere — delete.
- Missing `LANL_COPYING` / `LANL_README` referenced in `particle_advance.f90:3`
  — restore for license/provenance compliance.
- `Input/` ships `input.par.lte`/`.nlte` but the code opens literal `input.par`;
  the `make run` staging never creates `input.par` → as-shipped `make run` fails.

**Environment note.** The dev container here has **no Fortran compiler, no MPI,
no NumPy/h5py**. Building, reference-ensemble generation, and benchmarking must
run in CI (or a toolchain-equipped environment). CI is therefore a Phase‑0
deliverable, not an afterthought.

---

## 3. Target architecture

```
hyperlite/
├─ CMakeLists.txt              # top-level modern CMake
├─ CMakePresets.json           # gfortran / ifx / debug / mpi / openmp presets
├─ pyproject.toml              # packages the Python tools (hyperlite-tools)
├─ CITATION.cff, CHANGELOG.md, CONTRIBUTING.md, LICENSE, LANL_*
├─ src/
│  ├─ core/                    # kinds, physconst, random(Philox), timing, counters
│  ├─ io/                      # hdf5 read/write, namelist/yaml config, ascii export
│  ├─ grid/  gas/  group/  source/  transport/  output/   (former subdir libs)
│  └─ superlite.f90            # driver
├─ data/                       # atomic data -> generated data/atomic.h5 (+ raw kept out of tree or LFS)
├─ tools/                      # python package: stella2superlite, converters, plotting
├─ tests/
│  ├─ regression/              # harness + reference ensembles (HDF5) + tolerances
│  ├─ unit/                    # pFUnit unit tests (RNG, binsrch, EOS, specint...)
│  └─ cases/                   # small fast decks (smoke) + physics decks
├─ docs/                       # mkdocs/Sphinx: build, run, physics, I/O schema, dev
└─ .github/workflows/          # ci.yml (build+test matrix), release.yml
```

**Key module additions**
- `src/core/kinds.f90` — `dp = real64`, `sp = real32`, `i4 = int32`, etc.
  (`iso_fortran_env`), replacing `real*8`/`integer*2`/`logical*2` extensions.
- `src/core/random_philox.f90` — counter-based RNG (Philox‑4×32 or Threefry),
  keyed by `(stream_id, particle_global_id, event_counter)`; drop-in for the
  `rnd_r`/`rnd_i` interface so call sites are unchanged.
- `src/io/hdf5_io.f90` — thin wrappers over the HDF5 Fortran 2003 interface for
  the read/write seams listed above.
- `particlemod.f90` — particles become **SoA**: `prt_x(:), prt_mu(:), prt_om(:),
  prt_e(:), prt_e0(:), prt_wl(:), prt_isvacant(:), prt_ic(:), prt_ig(:), ...`
  (drop the trivially-zero `y,z` in 1‑D).

**HDF5 output schema** (single `output.h5`, one group per physical quantity):
```
/meta         : code version, git describe, build/run date, input echo, nmpi, nomp, seed
/grid         : igeom, dims, xarr/yarr/zarr, icell map, per-cell temp/radtemp/eraddens/nvol
/group        : wl boundaries, per-cell cap/capemit/emiss (opac dump)
/flux         : wl/mu/om axes, luminos[ng,nmu,nom], lumnum (counts), lumdev (MC error)
/profile      : ejecta structure columns + per-element mass fractions
/source       : src_number, src_luminos
```
A `tools/h5_to_ascii.py` exporter reproduces the legacy `output.*` tables for
humans/diffing. Atomic data collapses to a single `data/atomic.h5`
(ions/levels/lines/zeta/bf/ff as named datasets), read directly by the Fortran —
removing the fragile fixed-format reads and the `Atoms/` CWD dependency.

---

## 4. Performance strategy

The transport kernel is a textbook divergent Monte Carlo loop. The plan borrows
proven techniques from sibling radiative-transfer / MC-transport codes:

- **Get DDMC pulling its weight (biggest algorithmic lever).** DDMC replaces
  many small IMC steps in optically-thick cells with a single diffusion draw
  (Densmore, Cleveland; as in **SuperNu**). Audit `in_trn_tauddmc` / `in_taulump`
  and the per-cell `grp_t_cache` lump so thick cells actually take the cheap
  path; ensure `specintv` (flagged "this is slow" in `diffusion11.f90`) is
  cached and vectorized over groups.
- **AoS → SoA particles.** Enables SIMD across particles and cuts the
  gather/scatter cost of `ptcl = prt_particles(ipart)`. This is the enabling
  refactor; everything else builds on it. (Standard in GPU/Kokkos MC transport,
  e.g. OpenMC's event-based mode.)
- **Counter-based RNG (Philox/Threefry).** Per-lane, stateless draws keyed by
  particle id + event counter — removes the serial state recurrence and the
  per-thread state write-back, and makes results independent of scheduling
  (only the ensemble statistics matter). This is *the* decision that unblocks
  SIMD across particles.
- **Event-batched restructuring to defeat divergence.** Move from
  "one thread chases one particle through all its events" toward processing
  particles in batches by *event type* (transport vs diffusion vs boundary),
  the SIMT/event-based approach used to tame divergence in modern MC transport
  (Romano et al., OpenMC event mode). Keep the per-cell DDMC cache locality.
- **Kill the entry binary searches.** `grp_wl` is geometric and `grd_xarr`
  regular, so the three `binsrch` calls per particle become closed-form
  `log`-index arithmetic.
- **Threading by default with better load balance.** Enable OpenMP
  (`in_nomp` default > 1), switch the particle loop from `schedule(static,1)`
  to `guided`/`dynamic` to handle the wildly variable per-particle step counts
  (the load-imbalance issue ARTIS/Sedona address with nested MPI+OpenMP particle
  balancing).
- **Precision & memory.** Keep `grd_cap` in `real*4` where it's read-mostly to
  halve bandwidth on the dominant gather; align/pad SoA arrays for vector loads.

**GPU transport kernel (stretch goal, Phase 8).** Every choice in the CPU
refactor doubles as a GPU enabler, so the offload is an extension rather than a
rewrite:
- **SoA particles** map directly to coalesced device memory access.
- **Counter-based Philox** is stateless per lane — the standard RNG for GPU MC
  transport (no per-thread state to migrate host↔device).
- **Event-batched processing** is precisely the SIMT-divergence remedy GPUs need;
  the batches become kernel launches per event type (transport / diffusion /
  boundary / census), the design used by GPU MC transport codes (e.g. OpenMC's
  event-based GPU mode, Shift/Profugus, and Kokkos-based transport).

Later follow-on (beyond this overhaul): expansion-opacity line binning
(TARDIS/ARTIS-style) for the opacity build.

**Targets / how we'll know it worked:** on the reference smoke case, wall-clock
speedup measured at each phase, reported in CI. Concrete gate: ≥4× single-node
speedup (SoA + counter-RNG + OpenMP on a typical multicore CI runner) versus the
Phase‑0 serial baseline, with regression still passing within MC noise; plus an
OpenMP strong-scaling curve (efficiency reported, not gated to a hard number).

---

## 5. Phased implementation plan

Each phase is independently shippable and ends at a **success gate**. Phases 0–2
carry **zero intended change to physics output**; the regression harness (built
in Phase 0) guards every subsequent phase.

### Phase 0 — Stabilize, CI, and capture the baseline *(no behavior change)*
The safety net must exist before anything else moves.
1. Fix build rot so a stock build succeeds on **gfortran**:
   remove `TRANSPORT2/3` refs, fix the `output_source.0` typo, wire in
   `grid_volume.f`, fix the compiler-file usage message, add the missing
   `input.par` staging step. Compile with `-std=legacy -fallow-argument-mismatch`.
2. Delete `Data/Tinit.dat`; restore `LANL_COPYING`/`LANL_README` provenance.
3. Stand up **GitHub Actions CI**: a gfortran build job (matrix: serial + MPI +
   OpenMP) that compiles and runs the smallest case to completion.
4. Define the **smoke case** (root W7 deck, `in_ndim=64`, `in_grp_ng=250`,
   `in_src_n2s=20`, serial). Verify `in_ndim` vs `.str` header consistency.
5. **Generate the statistical reference:** run the smoke case *N≈20 times* with
   distinct seeds; store per-bin **mean and std** (spectrum bins, integrated
   luminosity, per-cell temperatures) plus the deterministic geometry/grid/group
   arrays as `tests/regression/reference/smoke.h5`. Record `(nmpi, nomp,
   code-version)` in metadata.
6. Write the **regression harness** (`tests/regression/compare.py` + a CTest
   wrapper): deterministic anchors (geometry, group grid, profile structure)
   compared at tight tolerance; stochastic fields required within `k·σ` of the
   reference ensemble (k≈4) **and** integrated/conserved quantities within a
   tight relative tolerance.
- **Gate 0:** clean gfortran build (serial/MPI/OpenMP) in CI; smoke case runs to
  completion; reference ensemble committed; harness runs green against the
  current code and *fails* on an injected perturbation (sanity check).

### Phase 1 — Modern CMake + CTest
1. Rewrite `CMakeLists.txt` as a proper multi-target project: `kinds`/core →
   per-component object libraries (`grid`, `gas`, `group`, `source`, `transport`,
   `output`, `io`, `misc`) → `superlite` executable. Use native CMake Fortran
   module-dependency tracking (retire `depend.sh`/`dependmodule.sh`).
2. `find_package` for MPI, OpenMP, HDF5; options `HYPERLITE_MPI`,
   `HYPERLITE_OPENMP` (mpimod selected via `configure_file` swap, as today).
3. Add `CMakePresets.json` (gfortran-release, gfortran-debug with
   `-fcheck=all`, ifx-release, +mpi/+openmp variants) and `install()` rules.
4. Register the Phase‑0 harness as a **CTest** test; delete the legacy
   Make/`check` machinery once parity is confirmed.
- **Gate 1:** `cmake --preset` builds all serial/MPI/OpenMP variants;
  `ctest` runs the regression green; ifx build attempted in CI (allowed to be
  non-blocking initially).

### Phase 2 — Repo restructure & Python packaging
1. Move sources into `src/**` per the target layout; update CMake source lists.
2. `pyproject.toml` packaging the tools as `hyperlite-tools`; **port
   `sgfilter.py` and `wlgenadd.py` from Python 2 → 3**; declare deps
   (numpy/scipy/matplotlib/h5py).
3. Move test decks under `tests/cases/`; fix/annotate the mismatched A4 decks.
- **Gate 2:** build green from the new layout; `pip install -e tools/` works and
  its console scripts run; regression still green (byte-for-byte identical binary
  → identical statistics).

### Phase 3 — HDF5 I/O *(behavior-preserving at the physics level)*
1. Add `src/io/hdf5_io.f90` wrappers; introduce the output schema (§3). Write
   `output.h5`; keep the legacy ASCII writers behind an `in_io_ascii` switch for
   one release, plus `tools/h5_to_ascii.py`.
2. Convert atomic data to `data/atomic.h5` via a one-time
   `tools/build_atomic_h5.py`; switch the Fortran readers to HDF5; remove the
   `Atoms/` CWD dependency. Remove the `< 1e-99` truncation hack (native binary).
3. Optional structure/grid inputs via HDF5 (with ASCII converters); config stays
   namelist (already standard-ish) — YAML front-end deferred.
4. Extend the harness to read HDF5 directly; verify HDF5 ↔ ASCII round-trip.
- **Gate 3:** HDF5 outputs reproduce the ASCII baseline within the same
  statistical tolerances; atomic-data HDF5 path produces identical opacities to
  the ASCII path (tight tol); regression green; runs no longer need `Atoms/`.

### Phase 4 — Fortran language modernization
1. Introduce `src/core/kinds.f90`; replace `real*8`/`integer*2`/`logical*2`
   with kind parameters module-by-module.
2. Remove statement functions (e.g. `dx(l)=...` in `particle_advance.f90`) and
   other obsolescent features; convert fixed-form `.f` → free-form where cheap,
   otherwise isolate legacy form.
3. Tighten flags toward `-std=f2018 -Wall`; triage warnings.
- **Gate 4:** builds under stricter flags on gfortran + ifx; warning count
  driven to an agreed floor; regression green throughout.

### Phase 5 — Performance I: data layout + counter-based RNG *(enabling)*
1. Convert `prt_particles` **AoS → SoA** in `particlemod.f90` and all call sites
   (`particle_advance.f90`, `SOURCE/*`, `transport11`/`diffusion11`).
2. Replace mzran with **Philox/Threefry** (`src/core/random_philox.f90`) keyed by
   particle global id + event counter; keep the `rnd_r`/`rnd_i` signatures.
3. **Re-baseline the reference ensemble** from the counter-RNG serial build
   (physics unchanged; draw order changes) and confirm the new ensemble is
   statistically consistent with the Phase‑0 ensemble (two-sample check on
   integrated luminosity + spectrum shape).
- **Gate 5:** SoA build correct; counter-RNG serial ensemble statistically
  consistent with the mzran baseline; regression green.

### Phase 6 — Performance II: threading + vectorization
1. Enable OpenMP by default (`in_nomp` default > 1); switch particle loop to
   `guided`/`dynamic` scheduling; per-lane RNG removes state contention.
2. Replace entry `binsrch` with log-index math; ensure `specintv`/DDMC group
   loops vectorize; align/pad SoA arrays; keep `grd_cap` in `real*4`.
3. Introduce event-batched processing of the particle loop to reduce divergence;
   preserve per-cell DDMC cache locality.
4. Add a **benchmark harness** (CI-reported wall-clock + OpenMP scaling on the
   smoke case and one mid-size case).
- **Gate 6:** ≥4× single-node speedup vs Phase‑0 serial baseline on the CI
  runner; regression within MC noise at 1, 2, and max threads and at
  `nmpi∈{1,2}`; scaling curve published.

### Phase 7 — Documentation & release
1. Rewrite `README`; author `docs/` (mkdocs or Sphinx): quickstart, build/run,
   physics overview, **I/O schema reference**, developer guide, migration notes
   (ASCII→HDF5).
2. `CHANGELOG.md`, `CONTRIBUTING.md`, `CITATION.cff`; version stamping via CMake
   + git; tagged release candidate; Zenodo metadata.
- **Gate 7:** docs build in CI; full matrix (gfortran+ifx × serial/MPI/OpenMP)
  green; regression + benchmark green; unit tests green; tag `v2.0.0-rc`.

### Phase 8 — GPU transport kernel *(stretch goal, post-release)*
Optional, undertaken only after Gate 6; the CPU refactor (SoA + counter-based
RNG + event batching) is its foundation, so this is an extension, not a rewrite.
1. **Choose the offload model.** Default: **OpenMP `target`** directives (stays
   in Fortran, single source, vendor-portable via gfortran/ifx/nvfortran).
   Alternative if more control is needed: extract the batched kernel to C++ and
   use **Kokkos** or vendor CUDA/HIP, called from Fortran via `iso_c_binding`.
   Decide with a spike benchmark before committing.
2. **Port the batched kernel.** Move the SoA particle arrays and the read-mostly
   grid/group/opacity state (`grd_cap` etc.) to the device; launch one kernel per
   event type over the particle batch; keep tallies in device memory with atomic
   or per-block-reduced scatter-adds into `grd_tally`/`grd_jrad`, reduced back to
   host per iteration. Philox draws are computed per-thread from
   `(particle_id, event_counter)` with no state transfer.
3. **Minimize host↔device traffic.** Opacities change once per outer iteration
   (`gas_update`), so upload them once per iteration and keep particles resident
   across the event loop; only tallies and census return to the host.
4. **Integrate with the existing parallelism.** GPU path selected at build time
   (`HYPERLITE_GPU=ON`) and/or runtime; one GPU per MPI rank (the existing cell
   decomposition + particle replication maps cleanly to multi-GPU); the CPU
   OpenMP path remains the default and fallback.
5. **CI.** Add a GPU runner (self-hosted or cloud) running the regression at
   loosened statistical tolerance plus the benchmark; keep it non-blocking for
   the main matrix so CPU CI never depends on GPU availability.
- **Gate 8:** GPU build produces regression results within MC noise of the CPU
  reference ensemble on the smoke + one physics case; documented device speedup
  vs the Phase‑6 CPU baseline on a representative case; GPU path is opt-in and the
  CPU path is unaffected when it is off.

---

## 6. Testing strategy (detail)

- **Regression (primary).** Statistical, per §Phase 0. Reference = HDF5 ensemble
  (mean/σ per bin + deterministic anchors). Pass = deterministic anchors within
  tight tol, stochastic fields within `k·σ`, conserved/integrated quantities
  within tight relative tol. Runs pinned to fixed `(nmpi, nomp)` per matrix cell.
  A deliberately-perturbed build must fail the harness (guards against a
  no-op/over-loose test).
- **Unit tests** (`tests/unit/`, **pFUnit**): RNG uniformity/independence and
  known-answer vectors for Philox; `binsrch`/log-index equivalence; `specint`
  Planck-integral against analytic values; EOS/opacity spot values; HDF5
  round-trip.
- **Smoke** (fast, every PR): the small W7 case to completion, few seconds.
- **Physics cases** (nightly): `s18`, `sn1999em-like`, `sn2017hcc-like`, A4
  LTE/NLTE — larger, looser statistical tolerances.
- **Benchmarks** (nightly): wall-clock + OpenMP scaling, tracked over time.

Because the environment lacks a Fortran toolchain, **all of the above run in CI**
(GitHub Actions); reference ensembles are generated once in CI and committed.

---

## 7. Risks & mitigations

| Risk | Mitigation |
|---|---|
| No golden outputs exist to anchor against | Phase 0 *generates* the reference ensemble from the current code before any change; harness must fail on an injected perturbation. |
| RNG swap changes results | Chosen statistical basis expects this; re-baseline in Phase 5 and prove two-sample consistency with the mzran ensemble. |
| Statistical test too loose (misses real regressions) | Tight tolerances on deterministic + conserved quantities; `k·σ` calibrated from the ensemble; perturbation self-test in CI. |
| AoS→SoA touches many call sites | Land behind Gate 5 with regression green at each step; SoA is mechanical and localized to particle access. |
| `real*4` opacity precision vs HDF5 native binary | Compare atomic-data HDF5 path against ASCII path at tight tol in Gate 3. |
| ifx/gfortran divergence | CI matrix from Phase 1; ifx non-blocking until Phase 4. |
| Atomic data size in git (~27 MB) | Store `data/atomic.h5` via Git LFS or a fetch script; keep raw ASCII out of the main tree. |
| GPU stretch (Phase 8): vendor portability / CI availability | Default to single-source OpenMP `target`; make the GPU path opt-in and non-blocking in CI; the release (`v2.0.0-rc`) ships at Gate 7 and does not depend on Phase 8. |

**Rollback:** every phase is a separate, revertible change gated by a green
regression + build matrix; no phase depends on an unshipped later phase.

---

## 8. Verification (end-to-end)

Because there is no Fortran toolchain in this environment, verification runs in
CI (or any box with gfortran/ifx + HDF5 + MPI). The canonical loop:

```bash
# Build (modern CMake, from Phase 1 on)
cmake --preset gfortran-openmp
cmake --build --preset gfortran-openmp -j

# Run the smoke case
ctest --preset gfortran-openmp -R regression-smoke --output-on-failure

# Full regression + benchmarks (nightly)
ctest --preset gfortran-openmp --output-on-failure
python tools/bench_report.py   # wall-clock + OpenMP scaling vs Phase-0 baseline
```

Pre-Phase-1 (legacy Make, Phase 0 only):
```bash
cp System/Makefile.compiler.gfortran Makefile.compiler
make all && make run           # after the Phase-0 build+staging fixes
python tests/regression/compare.py Run/output.* tests/regression/reference/smoke.h5
```

**Definition of done:** all seven gates passed; CI matrix
(gfortran+ifx × serial/MPI/OpenMP) green; regression within MC noise across the
thread/rank matrix; ≥4× single-node speedup on the smoke case vs the Phase‑0
serial baseline; HDF5 I/O with ASCII export; docs building; `v2.0.0-rc` tagged.
