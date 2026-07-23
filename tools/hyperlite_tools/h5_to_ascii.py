#!/usr/bin/env python3
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""Export a hyperlite output.h5 back to the legacy ASCII output.* tables.

Reproduces the layout of the legacy Fortran writers (one row block per
iteration for appendable fields).  Values are numerically identical to
the HDF5 data at full double precision; the text formatting matches the
legacy files except that 3-digit exponents keep their 'E' (the Fortran
writers drop it, e.g. 1.2345-100).
"""
import argparse
from pathlib import Path

import h5py
import numpy as np


def _wrow(f, arr, fmt='%12.4E'):
    f.write(''.join(fmt % v for v in np.atleast_1d(arr)) + '\n')


def _wirow(f, arr):
    f.write(''.join('%12d' % v for v in np.atleast_1d(arr)) + '\n')


def export(h5file, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with h5py.File(h5file, 'r') as h:
        _export_grid(h, outdir)
        _export_group(h, outdir)
        _export_flux(h, outdir)
        _export_total(h, outdir)
        _export_source(h, outdir)
        if 'profile' in h:
            _export_profile(h, outdir)


def _export_grid(h, outdir):
    g = h['grid']
    nx, ny, nz = (int(g.attrs[k]) for k in ('nx', 'ny', 'nz'))
    ncell = h['grid/icell'][()].size
    ncpr = nx
    nrow = -(-ncell // ncpr)
    with open(outdir / 'output.grd_grid', 'w') as f:
        f.write(' # %11d\n' % int(g.attrs['igeom']))
        f.write(' # %11d%12d%12d\n' % (nx, ny, nz))
        f.write(' # %11d%12d%12d\n' % (nrow * ncpr, nrow, ncpr))
        _wrow(f, g['xarr'][()])
        _wrow(f, g['yarr'][()])
        _wrow(f, g['zarr'][()])
        icell = g['icell'][()]  # h5py order (nz, ny, nx)
        for k in range(icell.shape[0]):
            for j in range(icell.shape[1]):
                _wirow(f, icell[k, j, :])
    # -- per-iteration cell fields, padded to nrow*ncpr like the writer
    fields = [('nvol', 'output.grd_nvol', True),
              ('temp', 'output.grd_temp', False),
              ('radtemp', 'output.grd_radtemp', False),
              ('capgrey', 'output.grd_capgrey', False),
              ('capemitgrey', 'output.grd_capemitgrey', False),
              ('capross', 'output.grd_capross', False),
              ('sig', 'output.grd_sig', False),
              ('eraddens', 'output.grd_eraddens', False),
              ('methodswap', 'output.grd_methodswap', True)]
    for name, fname, isint in fields:
        if name not in g:
            continue
        data = g[name][()]  # (niter, ncell)
        pad = nrow * ncpr - ncell
        with open(outdir / fname, 'w') as f:
            for it in range(data.shape[0]):
                row = np.pad(data[it], (0, pad))
                for i in range(nrow):
                    chunk = row[i * ncpr:(i + 1) * ncpr]
                    _wirow(f, chunk) if isint else _wrow(f, chunk)


def _export_group(h, outdir):
    g = h['group']
    ng = int(g.attrs['ng'])
    with open(outdir / 'output.grp_grid', 'w') as f:
        f.write('#%5d\n' % ng)
        _wrow(f, g['wl'][()])
    titles = {'cap': 'Multigroup opacity',
              'capemit': 'Multigroup "emission" opacity',
              'emiss': 'Multigroup emissivity'}
    for name, title in titles.items():
        if name not in g:
            continue
        data = g[name][()]  # (niter, ncell, ng)
        with open(outdir / f'output.grp_{name}', 'w') as f:
            for it in range(data.shape[0]):
                f.write('#%s\n' % title.rjust(20 if name == 'cap' else 30))
                f.write('#%5s%5s\n' % ('ng', 'nr'))
                f.write('#%5d%5d\n' % (ng, data.shape[1]))
                for i in range(data.shape[1]):
                    _wrow(f, data[it, i, :])


def _export_flux(h, outdir):
    g = h['flux']
    ng, nmu, nom = (int(g.attrs[k]) for k in ('ng', 'nmu', 'nom'))
    with open(outdir / 'output.flx_grid', 'w') as f:
        f.write(' # %11d%12d%12d\n' % (ng, nmu, nom))
        _wrow(f, g['wl'][()])
        _wrow(f, g['mu'][()])
        _wrow(f, g['om'][()])
    for name, isint in (('luminos', False), ('lumnum', True),
                        ('lumdev', False)):
        data = g[name][()]  # (niter, nmu*nom, ng)
        with open(outdir / f'output.flx_{name}', 'w') as f:
            for it in range(data.shape[0]):
                for jk in range(data.shape[1]):
                    row = data[it, jk, :]
                    if isint:
                        _wirow(f, row)
                    else:
                        # legacy truncation guard: zero values below 1e-99
                        _wrow(f, np.where(row > 1e-99, row, 0.0))


def _export_total(h, outdir):
    data = h['total/energy'][()]  # (niter, 4)
    with open(outdir / 'output.tot_energy', 'w') as f:
        f.write('#%5d\n' % 4)
        f.write('#' + ''.join(s.rjust(12) for s in
                              ('eout', 'evelo', 'sflux', 'sthermal')) + '\n')
        for it in range(data.shape[0]):
            f.write(' ' + ''.join('%12.4E' % v for v in data[it]) + '\n')


def _export_source(h, outdir):
    number = h['source/number'][()]   # (niter, ng, ncell)
    energy = h['source/energy'][()]
    with open(outdir / 'output.src_number', 'w') as f:
        for it in range(number.shape[0]):
            for i in range(number.shape[2]):
                _wirow(f, number[it, :, i])
    with open(outdir / 'output.src_luminos', 'w') as f:
        for it in range(energy.shape[0]):
            for i in range(energy.shape[2]):
                _wrow(f, energy[it, :, i], fmt='%14.4E')


def _export_profile(h, outdir):
    g = h['profile']
    cols = ['x_left', 'x_right', 'vx_left', 'vx_right', 'mass', 'rho',
            'vol', 'avg_temp', 'rad_temp', 'ye', 'n_e', 'n_atom']
    massfr = g['massfr'][()]  # (ncell, nelem)
    elements = g['massfr'].attrs['elements']
    if isinstance(elements, bytes):
        elements = elements.decode()
    symbols = elements.split()
    # legacy writer skips all-zero elements
    keep = [l for l in range(massfr.shape[1]) if massfr[:, l].any()]
    with open(outdir / 'output.profile', 'w') as f:
        f.write('#     x_left' + ''.join(c.rjust(12) for c in cols[1:]))
        f.write(''.join(symbols[l].rjust(12) for l in keep) + '\n')
        data = np.column_stack([g[c][()] for c in cols]
                               + [massfr[:, l] for l in keep])
        for row in data:
            _wrow(f, row)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('h5file', type=Path, help='hyperlite output.h5 file')
    ap.add_argument('-o', '--outdir', type=Path, default=Path('.'),
                    help='directory for the output.* files (default: cwd)')
    args = ap.parse_args()
    export(args.h5file, args.outdir)
    print(f'exported {args.h5file} -> {args.outdir}/output.*')


if __name__ == '__main__':
    main()
