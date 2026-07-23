#!/usr/bin/env python3
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""One-time converter: raw ASCII atomic data -> data/atomic.h5.

Bundles the fixed-format files read by the legacy Fortran readers into a
single self-describing HDF5 file, preserving raw values and units so the
Fortran post-processing (unit conversion, level remapping, line assembly)
is unchanged:

  data.ion            -> /ion   (levels per ion; ionization potential in Ky)
  data.zeta           -> /zeta  (nebular-approximation zeta table)
  Atoms/data.atom.*   -> /bbxs/zZZiII (line-list levels+lines, real*4)
  Atoms/data.nlte.*   -> /nlte/zZZiII (NLTE levels/lines/collisions/PI-RR)
  data.bf_verner      -> /bf    (Verner 1995/1996 photoionization tables)
  data.ff_sutherland  -> /ff    (Sutherland 1998 free-free gaunt factors)

Column layouts mirror the Fortran read formats (data.atom/data.ion are
fixed-column; the NLTE file is token-separated).
"""
import argparse
import re
from pathlib import Path

import h5py
import numpy as np

# element symbols 1..30 as in elemdatamod (lowercased in file names)
SYMBOLS = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne',
           'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca',
           'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn']

KY_RE = re.compile(r'\(Ky\):\s*([0-9.Ee+-]+)')


def parse_ion(path):
    """data.ion: per-ion blocks (icod header, Ky potential, level table)."""
    ions = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        tok = lines[i].split()
        if not tok:
            i += 1
            continue
        icod = int(tok[0])
        m = KY_RE.search(lines[i + 1])
        chi_ion = float(m.group(1)) if m else 0.0
        nlev = int(lines[i + 3].split()[0])
        lev = lines[i + 4:i + 4 + nlev]
        # Fortran format (f11.3,1x,a12,i5,i5,l6):
        # chi cols 1-11, label 13-24, ilevel 25-29 (unused), g 30-34,
        # meta 35-40
        chi = np.array([float(s[0:11]) for s in lev])
        g = np.array([int(s[29:34]) for s in lev], dtype=np.int32)
        meta = np.array([1 if 'T' in s[34:40] else 0 for s in lev],
                        dtype=np.int32)
        ions.append((icod, chi_ion, chi, g, meta))
        i += 4 + nlev
    return ions


def parse_zeta(path):
    with open(path) as f:
        lines = f.readlines()
    ntemp = int(lines[1].split()[0])
    temp = np.array([float(s) for s in lines[2:2 + ntemp]])
    nzeta = int(lines[2 + ntemp].split()[0])
    rows = lines[4 + ntemp:4 + ntemp + nzeta]
    iz = np.array([int(r.split()[0]) for r in rows], dtype=np.int32)
    ii = np.array([int(r.split()[1]) for r in rows], dtype=np.int32)
    val = np.array([[float(v) for v in r.split()[2:2 + ntemp]] for r in rows])
    return temp, iz, ii, val


def parse_atom(path):
    """data.atom.*: header + fixed-column level and line tables (real*4)."""
    with open(path) as f:
        lines = f.readlines()
    nlevel, nline = (int(t) for t in lines[3].split()[:2])
    lev = lines[4:4 + nlevel]
    lin = lines[4 + nlevel:4 + nlevel + nline]
    # levels: (f11.3,13x,2i5)
    lev_chi = np.array([float(s[0:11]) for s in lev], dtype=np.float32)
    lev_id = np.array([int(s[24:29]) for s in lev], dtype=np.int32)
    lev_g = np.array([int(s[29:34]) for s in lev], dtype=np.int32)
    # lines: (2i5,f7.3)
    lin_lev1 = np.array([int(s[0:5]) for s in lin], dtype=np.int32)
    lin_lev2 = np.array([int(s[5:10]) for s in lin], dtype=np.int32)
    lin_f = np.array([float(s[10:17]) for s in lin], dtype=np.float32)
    if 4 + nlevel + nline != len([s for s in lines if s.strip()]):
        raise ValueError(f'{path}: trailing data')
    return (nlevel, nline, lev_chi, lev_id, lev_g,
            lin_lev1, lin_lev2, lin_f)


def parse_nlte(path):
    """data.nlte.*: token-separated tables (file contains tabs)."""
    with open(path) as f:
        lines = f.readlines()
    nlevel, nline, ncoll = (int(t) for t in lines[3].split()[:3])
    o = 4
    lev = [s.split() for s in lines[o:o + nlevel]]
    o += nlevel
    lin = [s.split() for s in lines[o:o + nline]]
    o += nline
    coll = [s.split() for s in lines[o:o + ncoll]]
    o += ncoll
    pirr = [s.split() for s in lines[o:o + nlevel]]
    if o + nlevel != len([s for s in lines if s.strip()]):
        raise ValueError(f'{path}: trailing data')
    d = {}
    # levels: id  label  g  chi  n
    d['lev_id'] = np.array([int(t[0]) for t in lev], dtype=np.int32)
    d['lev_g'] = np.array([int(t[2]) for t in lev], dtype=np.int32)
    d['lev_chi'] = np.array([float(t[3]) for t in lev])
    d['lev_n'] = np.array([int(t[4]) for t in lev], dtype=np.int32)
    # radiative lines: lev1 lev2 f wl0
    d['lin_lev1'] = np.array([int(t[0]) for t in lin], dtype=np.int32)
    d['lin_lev2'] = np.array([int(t[1]) for t in lin], dtype=np.int32)
    d['lin_f'] = np.array([float(t[2]) for t in lin])
    d['lin_wl0'] = np.array([float(t[3]) for t in lin])
    # EIE collisions: lev1 lev2 C0 C1 C2 C3 n delE
    d['coll_lev1'] = np.array([int(t[0]) for t in coll], dtype=np.int32)
    d['coll_lev2'] = np.array([int(t[1]) for t in coll], dtype=np.int32)
    for k, c in enumerate(('C0', 'C1', 'C2', 'C3')):
        d[f'coll_{c}'] = np.array([float(t[2 + k]) for t in coll])
    d['coll_n'] = np.array([float(t[6]) for t in coll], dtype=np.float32)
    d['coll_delE'] = np.array([float(t[7]) for t in coll])
    # PI/RR per level: lev (C0..C3 n)_RR (C0..C3 n)_PI delE
    d['pi_lev'] = np.array([int(t[0]) for t in pirr], dtype=np.int32)
    for k, c in enumerate(('C0', 'C1', 'C2', 'C3')):
        d[f'rr_{c}'] = np.array([float(t[1 + k]) for t in pirr])
    d['rr_n'] = np.array([float(t[5]) for t in pirr], dtype=np.float32)
    for k, c in enumerate(('C0', 'C1', 'C2', 'C3')):
        d[f'pi_{c}'] = np.array([float(t[6 + k]) for t in pirr])
    d['pi_n'] = np.array([float(t[10]) for t in pirr], dtype=np.float32)
    d['pi_delE'] = np.array([float(t[11]) for t in pirr])
    return nlevel, nline, ncoll, d


def build(datadir, outfile):
    datadir = Path(datadir)
    with h5py.File(outfile, 'w') as h:
        h.create_group('meta')
        h['meta'].attrs['schema_version'] = 1
        h['meta'].attrs['source'] = 'built from raw ASCII by build_atomic_h5'

        # -- /ion + /zeta
        ions = parse_ion(datadir / 'data.ion')
        g = h.create_group('ion')
        g.attrs['nion'] = len(ions)
        g['icod'] = np.array([e[0] for e in ions], dtype=np.int32)
        g['chi_ion'] = np.array([e[1] for e in ions])
        nlev = np.array([len(e[2]) for e in ions], dtype=np.int32)
        g['nlev'] = nlev
        g['offset'] = np.cumsum(np.concatenate(([1], nlev[:-1]))).astype(
            np.int32)  # 1-based start index per ion
        g['lev_chi'] = np.concatenate([e[2] for e in ions])
        g['lev_g'] = np.concatenate([e[3] for e in ions])
        g['lev_meta'] = np.concatenate([e[4] for e in ions])

        temp, iz, ii, val = parse_zeta(datadir / 'data.zeta')
        g = h.create_group('zeta')
        g.attrs['ntemp'] = len(temp)
        g.attrs['nzeta'] = len(iz)
        g['temp'], g['iz'], g['ii'], g['val'] = temp, iz, ii, val

        # -- /bbxs from Atoms/data.atom.*
        g = h.create_group('bbxs')
        natom = 0
        for izel, sym in enumerate(SYMBOLS, start=1):
            for ion in range(1, izel + 1):
                p = datadir / 'Atoms' / f'data.atom.{sym}{ion}'
                if not p.exists():
                    continue
                (nlevel, nline, lev_chi, lev_id, lev_g,
                 lin_lev1, lin_lev2, lin_f) = parse_atom(p)
                gg = g.create_group(f'z{izel:02d}i{ion:02d}')
                gg.attrs['nlevel'] = nlevel
                gg.attrs['nline'] = nline
                gg['lev_chi'], gg['lev_id'], gg['lev_g'] = \
                    lev_chi, lev_id, lev_g
                gg['lin_lev1'], gg['lin_lev2'], gg['lin_f'] = \
                    lin_lev1, lin_lev2, lin_f
                natom += 1

        # -- /nlte from Atoms/data.nlte.*
        g = h.create_group('nlte')
        for izel, sym in enumerate(SYMBOLS, start=1):
            for ion in range(1, izel + 1):
                p = datadir / 'Atoms' / f'data.nlte.{sym}{ion}'
                if not p.exists():
                    continue
                nlevel, nline, ncoll, d = parse_nlte(p)
                gg = g.create_group(f'z{izel:02d}i{ion:02d}')
                gg.attrs['nlevel'] = nlevel
                gg.attrs['nline'] = nline
                gg.attrs['ncoll'] = ncoll
                for k, v in d.items():
                    gg[k] = v

        # -- /bf: Verner tables (values kept real*4 like the Fortran reader)
        raw = np.loadtxt(datadir / 'data.bf_verner', comments='#',
                         dtype=np.float32)
        if raw.shape != (1699 + 465, 9):
            raise ValueError(f'data.bf_verner: unexpected shape {raw.shape}')
        g = h.create_group('bf')
        g['ph1d'] = raw[:1699]
        g['ph2d'] = raw[1699:]

        # -- /ff: Sutherland gaunt factors, 81 u values x 41 g^2 values
        raw = np.loadtxt(datadir / 'data.ff_sutherland', skiprows=5)
        if raw.shape != (81 * 41, 3):
            raise ValueError(f'data.ff_sutherland: unexpected shape '
                             f'{raw.shape}')
        h.create_group('ff')
        h['ff/gff'] = raw[:, 2].reshape(41, 81)

    return natom


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('-d', '--datadir', type=Path, default=Path('data'),
                    help='directory with the raw ASCII data (default: data)')
    ap.add_argument('-o', '--out', type=Path, default=None,
                    help='output file (default: DATADIR/atomic.h5)')
    args = ap.parse_args()
    out = args.out or args.datadir / 'atomic.h5'
    natom = build(args.datadir, out)
    print(f'wrote {out} ({natom} line-list ions)')


if __name__ == '__main__':
    main()
