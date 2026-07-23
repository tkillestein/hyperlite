#!/usr/bin/env python3
"""HDF5 <-> ASCII round-trip check for one run directory.

This file is part of SuperLite. SuperLite is released under the terms of the
GNU GPLv3, see COPYING.

Reads the regression fields from both output.h5 and the legacy output.*
tables of the SAME run and requires them to agree to ASCII formatting
precision (the e12.4 writers carry 5 significant digits; values below
1e-99 are zeroed in ASCII by the legacy truncation guard).

Exit code 0 = formats agree, 1 = mismatch, 2 = parse error.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sloutput  # noqa: E402

RTOL = 1.2e-4   # half-ulp of a 5-significant-digit mantissa, with margin
ATOL = 1e-90    # legacy writers zero values below 1e-99


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('rundir', type=Path,
                    help='run directory with output.h5 AND output.* files')
    args = ap.parse_args()

    try:
        ascii_run = sloutput.read_run(args.rundir, source='ascii')
        h5_run = sloutput.read_run(args.rundir, source='h5')
    except (OSError, ValueError) as e:
        print(f'ERROR: cannot parse run outputs: {e}', file=sys.stderr)
        return 2

    fields = (sloutput.ANCHOR_FIELDS + sloutput.STOCHASTIC_FIELDS
              + sloutput.INTEGRATED_FIELDS)
    errors = []
    for field in fields:
        a, h = ascii_run[field], h5_run[field]
        if a.shape != h.shape:
            errors.append(f'{field}: shape ascii {a.shape} != h5 {h.shape}')
            continue
        if not np.allclose(h, a, rtol=RTOL, atol=ATOL):
            worst = np.max(np.abs(h - a) / (np.abs(a) + 1e-300))
            errors.append(f'{field}: max rel dev {worst:.2e} > {RTOL}')
        else:
            print(f'  {field:14s} ok')

    for e in errors:
        print(f'FAIL: {e}')
    print('ROUNDTRIP PASS' if not errors else 'ROUNDTRIP FAIL')
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
