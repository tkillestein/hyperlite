"""Parsers for SuperLite/hyperlite ASCII output files.

This file is part of SuperLite. SuperLite is released under the terms of the
GNU GPLv3, see COPYING.

Each ``output.*`` file is a whitespace-separated table of floats, with
optional ``#``-prefixed header lines.  Rows may be ragged (e.g.
``output.flx_grid`` holds one row of wavelength edges and one row of mu
edges), so the generic reader returns a list of 1-D arrays; helpers stack
them when rectangular.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

# Fortran writes 3-digit exponents without the 'E' (e.g. 1.234-100).
_FORTRAN_FLOAT_FIX = re.compile(r'(\d)([+-]\d{3})\b')


def read_rows(path: str | Path) -> list[np.ndarray]:
    """Read all non-comment rows of a SuperLite output file as float arrays."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line = _FORTRAN_FLOAT_FIX.sub(r'\1E\2', line)
            rows.append(np.array([float(tok) for tok in line.split()]))
    return rows


def read_array(path: str | Path) -> np.ndarray:
    """Read a rectangular output file into a 2-D array (squeezed to 1-D)."""
    rows = read_rows(path)
    if not rows:
        raise ValueError(f'no data rows in {path}')
    n = len(rows[0])
    if any(len(r) != n for r in rows):
        raise ValueError(f'ragged rows in {path}; use read_rows()')
    return np.vstack(rows).squeeze()


def read_run(rundir: str | Path) -> dict[str, np.ndarray]:
    """Parse the fields of one run directory used by the regression test.

    Returns a dict with deterministic anchors (grid/group/flux axes) and
    stochastic fields (spectra, temperatures, energy totals).
    """
    rundir = Path(rundir)
    out: dict[str, np.ndarray] = {}

    # -- deterministic anchors
    flx = read_rows(rundir / 'output.flx_grid')
    out['flx_wl'] = flx[0]
    out['flx_mu'] = flx[1]
    out['grp_wl'] = read_array(rundir / 'output.grp_grid')
    grd = read_rows(rundir / 'output.grd_grid')
    out['grd_xarr'] = grd[0]  # cell boundary radii

    # -- stochastic fields
    for name in ('flx_luminos', 'flx_lumnum', 'grd_temp', 'grd_radtemp',
                 'grd_eraddens'):
        out[name] = read_array(rundir / f'output.{name}')

    # -- energy totals: single row [eout, evelo, sflux, sthermal]
    out['tot_energy'] = np.atleast_1d(read_array(rundir / 'output.tot_energy'))

    # -- derived integrated quantities
    out['L_tot'] = np.atleast_1d(out['flx_luminos'].sum())
    return out


ANCHOR_FIELDS = ('flx_wl', 'flx_mu', 'grp_wl', 'grd_xarr')
STOCHASTIC_FIELDS = ('flx_luminos', 'flx_lumnum', 'grd_temp', 'grd_radtemp',
                     'grd_eraddens', 'tot_energy')
INTEGRATED_FIELDS = ('L_tot',)
