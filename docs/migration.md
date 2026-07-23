# Migration notes: ASCII → HDF5

hyperlite 2.0 moves primary I/O to HDF5. Everything legacy still works
behind flags; this page maps old workflows to new ones.

## Outputs

| Before (SuperLite 1.x) | Now |
|---|---|
| ~15 `output.*` ASCII tables | single `output.h5` ([schema](io-schema.md)) |
| ASCII always written | `in_io_ascii = t` keeps writing the legacy tables (default **on** for the 2.0 release; will default off later) |
| values `< 1e-99` printed with a truncated exponent (`1.2345-100`) | full-precision native doubles in HDF5; the truncation hack is gone |

Regenerate the ASCII tables from any `output.h5`:

```bash
h5-to-ascii output.h5 -o some_dir/
```

Reading outputs in Python:

```python
import h5py
h = h5py.File('output.h5')
wl  = h['flux/wl'][:]            # group boundaries (cm)
lum = h['flux/luminos'][-1,0,0]  # last iteration, erg/s per group
```

## Atomic data

| Before | Now |
|---|---|
| `Atoms/` directory + `data.*` fixed-format files **in the CWD** | `data/atomic.h5`, generated at build time, symlinked into the run directory |
| hard failure if `Atoms/` missing | `in_io_atomdata = 'h5'` (default) reads the bundle; `'asci'` restores the legacy reader |

Both paths produce byte-identical opacities.

## Inputs

Unchanged: `input.par` (namelist) and `input.str` (ASCII structure)
remain the input formats. HDF5 structure inputs are a possible future
addition.

## Behavior changes to be aware of

These are deliberate corrections, not format changes — spectra differ
from SuperLite 1.x beyond Monte Carlo noise:

1. **LTE partition-function fix** (2.0): the Saha EOS no longer freezes
   partition functions at the first cell's temperature. Total luminosity
   shifts **+15%** on the W7 reference case, and results no longer depend
   on the MPI rank count.
2. **DDMC Doppler-shift cache fix** (2.0): removes two spurious far-UV
   edge bins and coherently shifts the 2100–2700 Å region by 6–35%
   (total luminosity conserved to 2×10⁻⁵).
3. **RNG swap** (mzran → Philox): draw sequences differ; ensembles are
   statistically consistent. For a fixed seed, results are now
   independent of OpenMP thread count.
