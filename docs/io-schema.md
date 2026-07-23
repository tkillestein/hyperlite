# I/O schema reference

## `output.h5`

One file per run, one HDF5 group per physical quantity. Static axes are
written once; per-iteration fields append along a trailing unlimited
dimension (dimension order below is Fortran/h5py-reported:
`[iter, ...]` when read with h5py).

### `/meta` (attributes)

| Attribute | Meaning |
|---|---|
| `code` | `"hyperlite"` |
| `version` | release version, e.g. `2.0.0-rc` |
| `git_revision` | `git describe --tags --always --dirty` at build time |
| `schema_version` | integer, currently `1` |
| `nmpi`, `nomp` | MPI ranks and OpenMP threads of the run |
| `rnd_seed` | `in_rnd_seed` |
| `niter` | number of iterations |

### `/grid`

Attributes: `igeom`, `nx`, `ny`, `nz`.

| Dataset | Shape | Meaning |
|---|---|---|
| `xarr`, `yarr`, `zarr` | `[n+1]` | cell boundary coordinates (cm; 1-D spherical: radii) |
| `icell` | `[nx,ny,nz]` | cell index map |
| `temp` | `[niter,ncell]` | gas temperature (K) |
| `radtemp` | `[niter,ncell]` | radiation temperature (K) |
| `eraddens` | `[niter,ncell]` | radiation energy density |
| `nvol` | `[niter,ncell]` | source particles per cell |
| `capgrey`, `capemitgrey` | `[niter,ncell]` | Planck-mean absorption/emission opacity (cm⁻¹) |
| `capross` | `[niter,ncell]` | Rosseland-mean opacity (cm⁻¹) |
| `sig` | `[niter,ncell]` | Thomson scattering coefficient (cm⁻¹) |
| `methodswap` | `[niter,ncell]` | IMC↔DDMC method flag |

(Opacity fields follow the `in_io_opacdump` setting.)

### `/group`

Attribute: `ng`. Dataset `wl [ng+1]`: group wavelength boundaries (cm).
Per-iteration `[niter,ncell,ng]` (real32): `cap` (absorption opacity),
`capemit` (emission opacity), `emiss` (emissivity).

### `/flux`

Attributes: `ng`, `nmu`, `nom`. Axes `wl [ng+1]`, `mu [nmu+1]`,
`om [nom+1]`.

| Dataset | Shape | Meaning |
|---|---|---|
| `luminos` | `[niter,nom,nmu,ng]` | emergent luminosity per group (erg/s) |
| `lumnum` | `[niter,nom,nmu,ng]` | escaping packet counts |
| `lumdev` | `[niter,nom,nmu,ng]` | Monte Carlo deviance (error estimate) |

### `/total`

`energy [niter,4]`: `[eout, evelo, sflux, sthermal]` — escaped energy,
velocity-field work term, surface flux, thermal source. Energy
conservation: `eout + sflux ≈ 0`.

### `/source`

`number [niter,ncell]` (int), `energy [niter,ncell]`: source packet
counts and energy per cell.

### `/profile` (with `in_io_profdump = t`)

Static ejecta structure: `x_left/x_right`, `vx_left/vx_right`, `mass`,
`rho`, `vol`, `avg_temp`, `rad_temp`, `ye`, `n_e`, `n_atom`, and
`massfr [nelem,ncell]` with the element symbol list as attribute
`elements`.

## Legacy ASCII tables

With `in_io_ascii = t` the historic `output.*` files are also written
(`output.grd_temp`, `output.flx_luminos`, …). `h5-to-ascii` regenerates
them from `output.h5`, numerically identical at ASCII precision — except
that the legacy writers *truncated* 3-digit exponents (`1.2345-100`);
the HDF5 path stores full native doubles.

## `data/atomic.h5`

Generated at build time by `build-atomic-h5` from the raw fixed-format
ASCII data (`data/data.ion`, `data.zeta`, `data.bf_verner`,
`data.ff_sutherland`, `data/Atoms/`), and read directly by the Fortran
(`in_io_atomdata = 'h5'`, the default):

| Group | Contents |
|---|---|
| `/ion` | ionization energies and partition-function level data per element/ion |
| `/zeta` | recombination ζ tables |
| `/bbxs/zNNiMM` | bound-bound line lists per species |
| `/nlte/zNNiMM` | NLTE level/line/collision data per species |
| `/bf` | Verner bound-free cross-section fits |
| `/ff` | Sutherland free-free Gaunt-factor table |

Both paths (`'h5'` and `'asci'`) produce byte-identical opacities; the
HDF5 bundle removes the `Atoms/` working-directory dependency.

!!! tip
    Inspect any of these files with `h5dump -n output.h5` or
    `h5py`: `h5py.File('output.h5')['flux/luminos'][...]`.
