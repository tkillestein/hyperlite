# Physics overview

!!! note
    This is an orientation for users of the code, not a derivation. The
    method is described in
    [Wagle et al. 2023, ApJ 953 132](https://doi.org/10.3847/1538-4357/acda23)
    (SuperLite) and its
    SuperNu ancestry in Wollaeger et al. 2013/2014.

## What the code computes

SuperLite/hyperlite performs **post-processing spectral synthesis**: given
a snapshot of supernova ejecta (radii, velocities, densities,
temperatures, composition — from e.g. a Stella or FLASH model), it
transports photon packets through the ejecta in 1-D spherical geometry
and tallies the emergent multi-group flux (spectrum) and radiation-field
diagnostics. It targets **interacting transients** (SNe IIn, superluminous
SNe), where circumstellar interaction powers the light curve.

## IMC-DDMC hybrid transport

- **Implicit Monte Carlo (IMC)** (Fleck & Cummings 1971): photon packets
  undergo effective scattering/absorption with an implicit treatment of
  emission that stabilizes the radiation-matter coupling.
- **Discrete Diffusion Monte Carlo (DDMC)** (Densmore et al. 2007;
  Abdikamalov et al. 2012): in optically thick cells a packet takes
  single diffusion steps between cells instead of many short IMC flights.
  Cells switch method where the cell optical depth exceeds
  `in_trn_tauddmc` mean free paths; frequency groups are *lumped* above
  `in_taulump` (Wollaeger & van Rossum 2014).
- Packets convert between IMC and DDMC at cell/regime interfaces with
  the standard asymptotic boundary conditions.

The transport loop (per iteration) is: source packets from the boundary
luminosity `in_L_bol` and interior emission → advance every packet
through IMC flights / DDMC hops until census, escape, or absorption →
tally flux, energy deposition, and radiation-field moments.

## Gas state and opacities

Per iteration, the gas solver computes for each cell:

- **EOS/ionization**: Saha LTE by default, or NLTE excitation/ionization
  for selected species (`in_nlte`), solving level populations against the
  Monte Carlo radiation field estimator (`grd_jrad`).
- **Opacities**: bound-bound (line list), bound-free (Verner), free-free
  (Sutherland), and Thomson scattering, binned into `in_grp_ng`
  multigroup opacities on a log wavelength grid, plus Planck and
  Rosseland means. Atomic data comes from `data/atomic.h5`.

Iterations repeat (`in_niter`) until the temperature/radiation field
estimates are converged for the snapshot.

## Monte Carlo statistics

Emergent spectra carry per-bin Monte Carlo noise (`/flux/lumdev` tallies
the deviance; `/flux/lumnum` the packet counts). Two runs with different
seeds — or thread/rank counts in earlier versions — agree only within
this noise, which is why the regression suite compares against a
[statistical ensemble](developer.md#statistical-regression) rather than
a golden file. With hyperlite's per-particle counter-based RNG, results
are exactly reproducible for a fixed `(nmpi, in_rnd_seed)` regardless of
thread count.

## Lineage

- **SuperNu** (Wollaeger & van Rossum): time-dependent IMC-DDMC for SN
  light curves and spectra (LANL, GPLv3).
- **SuperLite v1.0** (Wagle et al. 2023): SuperNu adapted to steady-state
  snapshot spectral synthesis for interacting/superluminous SNe, adding
  NLTE.
- **hyperlite** (this fork): engineering modernization — build system,
  performance (SoA + Philox + threading), HDF5 I/O, tests/CI, packaging —
  plus correctness fixes (LTE partition-function cache, DDMC Doppler
  cache, MPI vacancy ledger). Physics algorithms are otherwise unchanged.
