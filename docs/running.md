# Running a simulation

## Anatomy of a run directory

`superlite` reads everything from the current working directory:

| File | Meaning |
|---|---|
| `input.par` | Fortran namelist (`&inputpars`) with all run parameters |
| `input.str` | ejecta structure: zones, velocities, mass, temperature, mass fractions |
| `atomic.h5` | atomic data (symlink `data/atomic.h5`; or set `in_io_atomdata='asci'` and provide `Atoms/` + `data.*` files) |
| `input.wlgrid` | optional custom wavelength grid (see `wlgenadd`) |

Staging a case from the repository:

```bash
mkdir run && cd run
ln -s ../data/* .                      # atomic.h5 (+ raw ASCII fallback)
ln -s ../tests/cases/w7/* .            # deck + structure
ln -sf input_w7.str input.str
cp input.par.lte input.par
../build/gfortran-openmp/superlite
```

Bundled cases under `tests/cases/`: `w7` (SN Ia, the smoke/regression
case), `s18`, `sn1999em-like`, `sn2017hcc-like`, `A4` — most with both
`input.par.lte` and `input.par.nlte` variants.

## Parameters of note

The full parameter list (with defaults) is the `inputpars` namelist in
`src/io/inputparmod.f`. Commonly touched:

| Parameter | Meaning |
|---|---|
| `in_name`, `in_comment` | run labels echoed to the outputs |
| `in_ndim` | number of zones (must match the `input.str` header) |
| `in_grp_ng` | number of wavelength groups |
| `in_src_n2s` | log2 of source particle count (`2**n2s` packets) |
| `in_L_bol` | bolometric luminosity of the inner boundary source |
| `in_nomp` | OpenMP threads; **0 = all available (default)**, 1 = serial |
| `in_rnd_seed` | RNG seed offset; 0 (default) reproduces the legacy stream family |
| `in_trn_tauddmc`, `in_taulump` | IMC→DDMC switching thresholds (mean free paths per cell) |
| `in_io_ascii` | also write the legacy ASCII `output.*` tables |
| `in_io_atomdata` | `'h5'` (default) or `'asci'` legacy atomic-data path |
| `in_nlte` | enable NLTE treatment |

## Parallelism

- **OpenMP** (default): threads share the particle loop; `in_nomp = 0`
  uses all available cores. Results are **independent of thread count**
  (per-particle RNG keying) — a 14-thread run is byte-identical to a
  serial run.
- **MPI** (`gfortran-mpi` build): ranks decompose *cells* for the
  gas/opacity solve and replicate *particles* for transport
  (`2**in_src_n2s / nmpi` packets per rank, disjoint RNG streams):

  ```bash
  mpirun -np 4 .../superlite
  ```

  MPI results agree with serial within Monte Carlo noise (opacities are
  byte-identical across rank counts since the partition-function fix).

- **Reproducibility**: for a fixed `(nmpi, in_rnd_seed)` the run is
  deterministic. Different seeds give statistically equivalent ensembles;
  the regression harness (`tests/regression/`) formalizes this.

## Outputs

A single self-describing `output.h5` ([schema](io-schema.md)), plus the
legacy ASCII tables while `in_io_ascii = t`. Convert between them with
`h5-to-ascii` ([tools](tools.md)); smooth flux spectra with `sgfilter`.
