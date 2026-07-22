#!/usr/bin/env python3
"""Generate the statistical regression reference ensemble.

This file is part of SuperLite. SuperLite is released under the terms of the
GNU GPLv3, see COPYING.

Runs the smoke case N times with distinct RNG seeds (via ``in_rnd_seed``),
collects the outputs, and stores per-bin ensemble mean/std of the stochastic
fields plus the deterministic anchors in an HDF5 reference file.

Usage (from the repo root, after building ``superlite``):

    python tests/regression/make_reference.py \
        --exe ./superlite --n 20 --jobs 8 \
        --out tests/regression/reference/smoke.h5
"""
from __future__ import annotations

import argparse
import datetime
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sloutput  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def stage_run(workdir: Path, exe: Path, input_par: Path, seed: int) -> Path:
    """Create a run directory for one ensemble member."""
    rundir = workdir / f'seed{seed:03d}'
    rundir.mkdir(parents=True)
    for src in (REPO / 'Data').iterdir():
        (rundir / src.name).symlink_to(src)
    for src in (REPO / 'Input').iterdir():
        (rundir / src.name).symlink_to(src)
    (rundir / 'input.str').symlink_to(rundir / 'input_w7.str')
    # -- inject the ensemble seed into the namelist
    text = input_par.read_text()
    assert 'in_rnd_seed' not in text, 'input deck already sets in_rnd_seed'
    text = text.replace('&inputpars', f'&inputpars\n in_rnd_seed = {seed}', 1)
    (rundir / 'input.par').write_text(text)
    return rundir


def run_member(rundir: Path, exe: Path) -> dict[str, np.ndarray]:
    log = rundir / 'stdout.log'
    with open(log, 'w') as f:
        subprocess.run([str(exe)], cwd=rundir, stdout=f, stderr=subprocess.STDOUT,
                       check=True)
    if 'SuperLite finished' not in log.read_text():
        raise RuntimeError(f'{rundir}: run did not finish cleanly')
    return sloutput.read_run(rundir)


def git_describe() -> str:
    try:
        return subprocess.run(
            ['git', '-C', str(REPO), 'describe', '--tags', '--long', '--dirty',
             '--always'], capture_output=True, text=True, check=True,
        ).stdout.strip()
    except subprocess.CalledProcessError:
        return 'unknown'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe', type=Path,
                    default=REPO / 'build' / 'gfortran-serial' / 'superlite')
    ap.add_argument('--input-par', type=Path,
                    default=REPO / 'Input' / 'input.par.lte')
    ap.add_argument('--n', type=int, default=20, help='ensemble size')
    ap.add_argument('--jobs', type=int, default=4, help='concurrent runs')
    ap.add_argument('--out', type=Path,
                    default=Path(__file__).parent / 'reference' / 'smoke.h5')
    ap.add_argument('--keep-runs', action='store_true',
                    help='keep the per-seed run directories')
    args = ap.parse_args()

    exe = args.exe.resolve()
    if not exe.exists():
        ap.error(f'executable not found: {exe}')

    workdir = Path(tempfile.mkdtemp(prefix='hyperlite-ensemble-'))
    print(f'ensemble workdir: {workdir}')

    rundirs = [stage_run(workdir, exe, args.input_par, seed)
               for seed in range(args.n)]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda d: run_member(d, exe), rundirs))
    print(f'{args.n} ensemble members finished')

    # -- deterministic anchors must be identical across the ensemble
    ref = results[0]
    for field in sloutput.ANCHOR_FIELDS:
        for i, res in enumerate(results[1:], 1):
            if not np.array_equal(ref[field], res[field]):
                raise RuntimeError(f'anchor {field} differs in member {i}')

    # -- write reference
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.out, 'w') as f:
        meta = f.create_group('meta')
        meta.attrs['code_version'] = git_describe()
        meta.attrs['date'] = datetime.datetime.now(datetime.UTC).isoformat()
        meta.attrs['n_ensemble'] = args.n
        meta.attrs['nmpi'] = 1
        meta.attrs['nomp'] = 1
        meta.attrs['case'] = str(args.input_par.name)

        g = f.create_group('anchors')
        for field in sloutput.ANCHOR_FIELDS:
            g.create_dataset(field, data=ref[field])

        g = f.create_group('stochastic')
        for field in sloutput.STOCHASTIC_FIELDS + sloutput.INTEGRATED_FIELDS:
            stack = np.stack([res[field] for res in results])
            grp = g.create_group(field)
            grp.create_dataset('mean', data=stack.mean(axis=0))
            grp.create_dataset('std', data=stack.std(axis=0, ddof=1))

    print(f'reference written: {args.out}')
    if not args.keep_runs:
        shutil.rmtree(workdir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
