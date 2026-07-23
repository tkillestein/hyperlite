#!/usr/bin/env python3
"""Two-sample consistency check between two reference ensembles.

This file is part of SuperLite. SuperLite is released under the terms of the
GNU GPLv3, see COPYING.

Used when re-baselining the regression reference after an intended change
to the draw sequence (e.g. the Phase-5 mzran -> Philox RNG swap): the new
ensemble must be statistically consistent with the old one, i.e. per-bin
mean differences within k standard errors (pooled), with anchors identical.

    python tests/regression/compare_ensembles.py OLD.h5 NEW.h5

Exit code 0 = consistent, 1 = ensembles differ, 2 = usage error.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sloutput  # noqa: E402

ANCHOR_RTOL = 1e-6


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('old', type=Path)
    ap.add_argument('new', type=Path)
    ap.add_argument('-k', type=float, default=5.0,
                    help='max allowed |z| per bin (default 5; ~450 bins)')
    ap.add_argument('--frac3', type=float, default=0.02,
                    help='allowed fraction of bins with |z| > 3')
    args = ap.parse_args()

    errors: list[str] = []
    with h5py.File(args.old) as fo, h5py.File(args.new) as fn:
        n_o = int(fo['meta'].attrs['n_ensemble'])
        n_n = int(fn['meta'].attrs['n_ensemble'])
        print(f'old: {args.old} (n={n_o})  {dict(fo["meta"].attrs)}')
        print(f'new: {args.new} (n={n_n})')

        for field in sloutput.ANCHOR_FIELDS:
            a, b = fo['anchors'][field][()], fn['anchors'][field][()]
            if a.shape != b.shape or not np.allclose(a, b, rtol=ANCHOR_RTOL):
                errors.append(f'anchor {field} differs')

        for field in (sloutput.STOCHASTIC_FIELDS + sloutput.INTEGRATED_FIELDS):
            mo = fo['stochastic'][field]['mean'][()]
            so = fo['stochastic'][field]['std'][()]
            mn = fn['stochastic'][field]['mean'][()]
            sn = fn['stochastic'][field]['std'][()]
            if mo.shape != mn.shape:
                errors.append(f'{field}: shape {mn.shape} != {mo.shape}')
                continue
            # -- standard error of the mean difference; floor degenerate bins
            scale = max(np.max(np.abs(mo)), 1e-300)
            se = np.sqrt(so**2/n_o + sn**2/n_n)
            se = np.maximum(se, 1e-3*scale/np.sqrt(min(n_o, n_n)))
            z = np.abs(mn - mo)/se
            n3 = int(np.sum(z > 3))
            print(f'  {field:14s} max|z| {z.max():5.2f}  '
                  f'{n3}/{z.size} bins > 3')
            if z.max() > args.k:
                errors.append(f'{field}: max|z| {z.max():.2f} > {args.k}')
            elif n3 > max(1, args.frac3*z.size):
                errors.append(f'{field}: {n3}/{z.size} bins with |z| > 3')

    for e in errors:
        print(f'FAIL: {e}')
    print('ENSEMBLES CONSISTENT' if not errors else 'ENSEMBLES DIFFER')
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
