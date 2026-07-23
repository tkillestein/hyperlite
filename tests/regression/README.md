# Statistical regression harness

SuperLite/hyperlite is a Monte Carlo code: its outputs are deterministic only
for a fixed `(nmpi, in_nomp)` and RNG. To allow modernization work (RNG swap,
threading, reordering) without losing regression coverage, the reference is a
**statistical ensemble**, not a single golden file.

## Files

- `sloutput.py` — parsers for the run outputs: the ASCII `output.*` tables
  and, via `read_run(dir, source='h5')`, the HDF5 `output.h5` (both yield
  identically-shaped fields).
- `roundtrip.py` — requires the HDF5 and ASCII outputs of the *same* run to
  agree to ASCII formatting precision (~5 significant digits).
- `make_reference.py` — runs the smoke case N times with distinct
  `in_rnd_seed` values and stores per-bin ensemble mean/std plus the
  deterministic anchors in `reference/smoke.h5`.
- `compare.py` — checks one run directory against the reference:
  - deterministic anchors (grid/group/flux axes) at tight relative tolerance;
  - stochastic fields (spectrum, temperatures, energy totals) within
    `k·σ` per bin (k=4), with a bounded outlier fraction and a hard
    `1.5k·σ` limit;
  - integrated quantities (total luminosity) within `k·σ` **and** a tight
    relative tolerance;
  - energy conservation (`eout + sflux ≈ 0`).
  - `--perturb FRAC` runs a self-test: the perturbed comparison **must
    fail**, guarding against an over-loose test.

## Regenerating the reference

Only regenerate when a change is *intended* to alter the draw sequence (e.g.
an RNG swap) and after demonstrating two-sample consistency with the previous
ensemble:

```bash
cmake --preset gfortran-serial && cmake --build --preset gfortran-serial -j
python tests/regression/make_reference.py \
    --exe build/gfortran-serial/superlite --n 20 --jobs 8
```

The reference records `(code version, date, n_ensemble, nmpi, nomp, case)` in
`/meta`. The committed reference was generated with `nmpi=1, nomp=1` from the
`tests/cases/w7/input.par.lte` W7 smoke case.

## Checking a run

```bash
ctest --preset gfortran-openmp   # smoke-run + regression-smoke[-h5] + roundtrip + selftest
```

or manually against any run directory:

```bash
python tests/regression/compare.py RUNDIR -v          # exit 0 = pass
python tests/regression/compare.py RUNDIR --source h5 -v  # read output.h5 instead
python tests/regression/compare.py RUNDIR --perturb 0.05  # self-test: must fail
python tests/regression/roundtrip.py RUNDIR           # h5 <-> ascii agreement
```
