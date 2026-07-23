# Contributing to hyperlite

Thanks for your interest! hyperlite is GPLv3; by contributing you agree
to license your work under the same terms.

## Getting set up

```bash
cmake --preset gfortran-openmp && cmake --build --preset gfortran-openmp -j
ctest --preset gfortran-openmp          # must be green before and after
cd tools && uv sync                     # Python tools + dev lint/type tools
```

See the [developer guide](https://tkillestein.github.io/hyperlite/developer/)
(or `docs/developer.md`) for layout, tests, and benchmarks.

## Ground rules

- **Fortran**: modern F2018, `-std=f2018 -Wall` clean on gfortran (the CI
  enforces zero warnings). Use the kinds from `src/core/kindmod.f90`.
  New files are free-form `.f90`.
- **Python** (`tools/`): `uv run ruff check .` and `uv run ty check` must
  pass; CLIs use rich-click (except the build/CI scripts
  `build_atomic_h5.py` and `bench.py`, which stay stdlib-argparse so they
  run without the package installed).
- **Every change runs the regression.** The statistical harness
  (`tests/regression/`) guards physics output. If your change is
  *intended* to alter results, say so explicitly, demonstrate consistency
  (`compare_ensembles.py`) or justify the shift, and re-baseline the
  reference in a dedicated commit documenting the physics.
- **Determinism**: transport must remain independent of OpenMP thread
  count (per-particle RNG keying). Don't introduce order-dependent
  tallies outside the established reductions.
- Add or extend a unit test (`tests/unit/`) when touching `src/core`
  numerics.

## Pull requests

1. Branch from `main`; keep PRs focused.
2. CI must be green (build matrix, ctest, tools lint/type, docs build).
3. Update `CHANGELOG.md` under an *Unreleased* heading, and the docs if
   behavior or interfaces change.
4. Reference issues where applicable; describe *why*, not just *what*.

## Reporting issues

Please include: the preset/compiler, `superlite` banner output (version +
git revision), the input deck, and for physics discrepancies the
`compare.py -v` output against the committed reference.
