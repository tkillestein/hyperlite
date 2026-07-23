# Developer guide

## Repository layout

```
├─ CMakeLists.txt / CMakePresets.json   # build (see Building)
├─ cmake/GenerateVersion.cmake          # build-time git version stamp
├─ src/
│  ├─ core/       # kinds, physical constants, Philox RNG, misc, timing, MPI shims
│  ├─ io/         # input namelist/structure readers, HDF5 wrappers (+stub)
│  ├─ grid/ gas/ group/ source/ transport/ output/
│  └─ superlite.f90                     # driver
├─ data/          # raw atomic data + generated atomic.h5 (gitignored)
├─ tools/         # hyperlite-tools Python package (uv + ruff + ty)
├─ tests/
│  ├─ unit/       # self-checking Fortran programs (CTest: unit-*)
│  ├─ regression/ # statistical harness + committed reference ensemble
│  └─ cases/      # input decks (w7 = smoke/regression case)
└─ docs/          # this site (mkdocs-material)
```

Conventions: modern Fortran 2018, `-Wall` clean; kinds from
`src/core/kindmod.f90` (`real(dp)` etc.); fixed-form `.f` files are
legacy-but-conforming; new code is free-form `.f90`.

## Tests

`ctest --preset <preset>` runs, in order:

| Test | What it checks |
|---|---|
| `unit-philox` | Random123 known-answer vectors for Philox-4x32-10; draw range; per-particle keying determinism; uniformity |
| `unit-binsrch` | interval search vs a linear-scan reference; edge/out-of-range behavior |
| `unit-specint` | Planck integrals against analytic values (π⁴/15, 2ζ(3)); Simpson additivity |
| `smoke-run` | the full W7 LTE case runs to completion (fixture for the rest) |
| `regression-smoke[-h5]` | statistical regression vs the reference ensemble (ASCII and HDF5 readers) |
| `regression-roundtrip` | HDF5 ↔ ASCII outputs of the same run agree |
| `regression-selftest` | a perturbed comparison **must fail** (guards test looseness) |
| `smoke-run-mpi` | (MPI preset) 2-rank completion check |

### Statistical regression

The reference is a **20-seed ensemble** (per-bin mean/σ + deterministic
anchors) in `tests/regression/reference/smoke.h5`, generated with
`make_reference.py`. `compare.py` requires deterministic anchors at tight
tolerance, stochastic fields within `k·σ` (k=4), integrated quantities
within tight relative tolerance, and energy conservation. Full details:
`tests/regression/README.md`.

Regenerate the reference **only** for an intended draw-sequence or
physics change, and demonstrate two-sample consistency
(`compare_ensembles.py`) or document the shift (see the re-baseline
history in the regression README).

## RNG design

`src/core/randommod.f90` implements counter-based **Philox-4x32-10**
(Salmon et al. 2011). Transport draws are a pure function of
`(key0, epoch, ipart)` — rank/seed, iteration, particle slot — so the
draw sequence is independent of OpenMP thread count and scheduling.
Service draws (source sampling) use a stateful stream per rank; scratch
streams are keyed per thread. The 32×32 multiply uses 16-bit limbs (no
unsigned/int128), verified against the Random123 known-answer vectors in
`unit-philox`.

## Benchmarks

```bash
uv run hyperlite-bench --exe build/gfortran-openmp/superlite --nomp 1,4,0
```

runs a reduced W7 case across thread counts (`--n2s` scales particle
count; `--baseline-exe` compares two builds). CI publishes a reduced
benchmark on every push. Phase-6 result on a 12-core laptop: 244 s
(Phase-0 serial) → 45 s (14 threads), **5.4×**.

## CI

`.github/workflows/ci.yml`: build+test matrix
(`gfortran-serial/-openmp/-mpi`), a threaded (`nomp=2`) regression run, a
benchmark step, the Python tools job (uv sync, ruff, ty, console-script
smoke tests, pip-install check), a docs build (`mkdocs build --strict`),
and a non-blocking ifx build.

## Releasing

1. Update `CHANGELOG.md`; bump `project(... VERSION ...)` in
   `CMakeLists.txt` (+ suffix), `tools/pyproject.toml`, `CITATION.cff`.
2. CI green on the release commit.
3. Tag: `git tag -a v2.0.0-rc -m "..." && git push --tags`. The banner
   and `output.h5` `/meta` pick the tag up via `git describe`.

See `CONTRIBUTING.md` for the contribution workflow, and `OVERHAUL.md`
for the full modernization plan and its phase history.
