#!/usr/bin/env python3
"""Statistical regression check: compare one run against the reference ensemble.

This file is part of SuperLite. SuperLite is released under the terms of the
GNU GPLv3, see COPYING.

Pass criteria
-------------
* **Deterministic anchors** (grid/group/flux axes): must match the reference
  to tight relative tolerance (they are seed-independent).
* **Stochastic fields** (spectra, temperatures, energy totals): each bin must
  lie within ``k·sigma`` of the ensemble mean (k=4 by default); at most
  ``--outlier-frac`` of bins may lie between k·sigma and the hard limit of
  1.5·k·sigma (expected tail for ~hundreds of Gaussian bins).
* **Integrated quantities** (total luminosity): within k·sigma AND a tight
  relative tolerance of the ensemble mean.
* **Conservation**: eout + sflux must cancel to a tight relative tolerance.

Exit code 0 = pass, 1 = regression detected, 2 = usage/parse error.

Usage:

    python tests/regression/compare.py RUNDIR [REFERENCE.h5]
    python tests/regression/compare.py RUNDIR --perturb 0.05   # self-test: must FAIL
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sloutput  # noqa: E402

ANCHOR_RTOL = 1e-6      # ASCII outputs carry ~5 significant digits
INTEGRATED_RTOL = 5e-3  # tight relative tolerance on conserved/integrated
CONSERVATION_RTOL = 1e-3


def check_anchors(run: dict, ref: h5py.File, errors: list[str]) -> None:
    for field in sloutput.ANCHOR_FIELDS:
        want = ref['anchors'][field][()]
        got = run[field]
        if got.shape != want.shape:
            errors.append(f'anchor {field}: shape {got.shape} != {want.shape}')
        elif not np.allclose(got, want, rtol=ANCHOR_RTOL, atol=0.0):
            worst = np.max(np.abs(got - want) / np.maximum(np.abs(want), 1e-300))
            errors.append(f'anchor {field}: max rel dev {worst:.2e} > {ANCHOR_RTOL}')


def check_stochastic(run: dict, ref: h5py.File, k: float, outlier_frac: float,
                     errors: list[str], verbose: bool) -> None:
    for field in sloutput.STOCHASTIC_FIELDS:
        mean = ref['stochastic'][field]['mean'][()]
        std = ref['stochastic'][field]['std'][()]
        got = run[field]
        if got.shape != mean.shape:
            errors.append(f'{field}: shape {got.shape} != {mean.shape}')
            continue
        # -- effective sigma: floor degenerate bins (zero variance in the
        #    ensemble, e.g. empty spectral bins) at a small fraction of the
        #    field scale so they cannot divide by zero or over-trigger.
        scale = np.max(np.abs(mean))
        sigma = np.maximum(std, 1e-3 * scale)
        ndev = np.abs(got - mean) / sigma
        nbad = int(np.sum(ndev > k))
        nhard = int(np.sum(ndev > 1.5 * k))
        if verbose:
            print(f'  {field:14s} max|dev| {ndev.max():5.2f} sigma, '
                  f'{nbad}/{ndev.size} bins > {k} sigma')
        if nhard:
            errors.append(f'{field}: {nhard} bins beyond {1.5*k:.0f} sigma '
                          f'(max {ndev.max():.1f})')
        elif nbad > max(1, outlier_frac * ndev.size):
            errors.append(f'{field}: {nbad}/{ndev.size} bins beyond {k} sigma')


def check_integrated(run: dict, ref: h5py.File, k: float,
                     errors: list[str], verbose: bool) -> None:
    for field in sloutput.INTEGRATED_FIELDS:
        mean = ref['stochastic'][field]['mean'][()]
        std = ref['stochastic'][field]['std'][()]
        got = run[field]
        dev = np.abs(got - mean)
        tol = np.maximum(k * std, INTEGRATED_RTOL * np.abs(mean))
        if verbose:
            print(f'  {field:14s} value {got[0]:.4e} ref {mean[0]:.4e} '
                  f'({dev[0]/max(std[0], 1e-300):.2f} sigma)')
        if np.any(dev > tol):
            errors.append(f'integrated {field}: |{got}-{mean}| > {tol}')


def check_conservation(run: dict, errors: list[str]) -> None:
    eout, _evelo, sflux, _sthermal = run['tot_energy']
    if abs(eout + sflux) > CONSERVATION_RTOL * abs(eout):
        errors.append(f'conservation: eout+sflux = {eout + sflux:.3e} '
                      f'(eout {eout:.3e})')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('rundir', type=Path, help='run directory with output.* files')
    ap.add_argument('reference', type=Path, nargs='?',
                    default=Path(__file__).parent / 'reference' / 'smoke.h5')
    ap.add_argument('-k', type=float, default=4.0, help='sigma multiplier')
    ap.add_argument('--outlier-frac', type=float, default=0.01,
                    help='allowed fraction of bins in (k, 1.5k] sigma')
    ap.add_argument('--perturb', type=float, default=0.0, metavar='FRAC',
                    help='self-test: scale stochastic fields by (1+FRAC); '
                         'the comparison is then EXPECTED to fail')
    ap.add_argument('-v', '--verbose', action='store_true')
    args = ap.parse_args()

    try:
        run = sloutput.read_run(args.rundir)
    except (OSError, ValueError) as e:
        print(f'ERROR: cannot parse run outputs: {e}', file=sys.stderr)
        return 2
    if not args.reference.exists():
        print(f'ERROR: reference not found: {args.reference}', file=sys.stderr)
        return 2

    if args.perturb:
        for field in sloutput.STOCHASTIC_FIELDS + sloutput.INTEGRATED_FIELDS:
            run[field] = run[field] * (1.0 + args.perturb)

    errors: list[str] = []
    with h5py.File(args.reference, 'r') as ref:
        if args.verbose:
            meta = dict(ref['meta'].attrs)
            print(f'reference: {args.reference} {meta}')
        check_anchors(run, ref, errors)
        check_stochastic(run, ref, args.k, args.outlier_frac, errors,
                         args.verbose)
        check_integrated(run, ref, args.k, errors, args.verbose)
    check_conservation(run, errors)

    failed = bool(errors)
    for e in errors:
        print(f'FAIL: {e}')
    if args.perturb:
        if failed:
            print(f'SELF-TEST PASS: perturbation {args.perturb} was detected')
            return 0
        print('SELF-TEST FAIL: perturbation was NOT detected '
              '(regression test is too loose)')
        return 1
    print('REGRESSION PASS' if not failed else 'REGRESSION FAIL')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
