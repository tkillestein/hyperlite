#!/usr/bin/env python3
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""Benchmark harness: wall-clock and OpenMP scaling for the W7 smoke case.

Stages the case in a temporary directory, runs the given superlite binary
at each requested thread count, and reports wall-clock time, speedup, and
parallel efficiency.  With --baseline-exe, also reports the speedup of
EXE over the baseline binary at each thread count.

    hyperlite-bench --exe build/gfortran-openmp/superlite --nomp 1,2,4,8
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def stage(repo, workdir, n2s, nomp):
    workdir.mkdir(parents=True, exist_ok=True)
    for sub in ('data',):
        for f in (repo/sub).iterdir():
            dst = workdir/f.name
            if not dst.exists():
                dst.symlink_to(f)
    for f in (repo/'tests'/'cases'/'w7').iterdir():
        dst = workdir/f.name
        if not dst.exists():
            dst.symlink_to(f)
    s = (workdir/'input.str')
    if s.is_symlink() or s.exists():
        s.unlink()
    s.symlink_to(workdir/'input_w7.str')
    text = (repo/'tests'/'cases'/'w7'/'input.par.lte').read_text()
    if n2s is not None:
        text = text.replace('in_src_n2s = 20', f'in_src_n2s = {n2s}')
    text = text.replace(' in_nomp = 1', f' in_nomp = {nomp}')
    (workdir/'input.par').write_text(text)


def run_once(exe, workdir):
    for f in workdir.glob('output.*'):
        f.unlink()
    t0 = time.monotonic()
    r = subprocess.run([str(exe)], cwd=workdir, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    t1 = time.monotonic()
    if r.returncode != 0 or 'SuperLite finished' not in r.stdout:
        print(r.stdout[-2000:], file=sys.stderr)
        raise RuntimeError(f'{exe}: run failed in {workdir}')
    return t1 - t0


def bench(exe, repo, tmp, nomps, n2s, repeats):
    times = {}
    for nomp in nomps:
        wd = tmp/f'nomp{nomp}'
        stage(repo, wd, n2s, nomp)
        best = min(run_once(exe, wd) for _ in range(repeats))
        times[nomp] = best
        print(f'  nomp={nomp:<3d} {best:8.2f} s', flush=True)
    return times


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe', type=Path, required=True)
    ap.add_argument('--baseline-exe', type=Path, default=None,
                    help='optional second binary to compare against')
    ap.add_argument('--repo', type=Path,
                    default=Path(__file__).resolve().parents[2])
    ap.add_argument('--nomp', type=str, default='1,2,4',
                    help='comma-separated thread counts')
    ap.add_argument('--n2s', type=int, default=None,
                    help='override in_src_n2s (default: deck value, 20)')
    ap.add_argument('--repeats', type=int, default=1,
                    help='runs per configuration (best-of)')
    args = ap.parse_args()
    nomps = [int(t) for t in args.nomp.split(',')]
    args.exe = args.exe.resolve()  # runs execute with cwd=workdir
    if args.baseline_exe:
        args.baseline_exe = args.baseline_exe.resolve()

    with tempfile.TemporaryDirectory(prefix='hyperlite-bench-') as td:
        tmp = Path(td)
        print(f'exe: {args.exe}')
        times = bench(args.exe, args.repo, tmp/'cur', nomps, args.n2s,
                      args.repeats)
        base = None
        if args.baseline_exe:
            print(f'baseline: {args.baseline_exe}')
            base = bench(args.baseline_exe, args.repo, tmp/'base', nomps,
                         args.n2s, args.repeats)

    t1 = times[nomps[0]]
    print()
    print('nomp   wall[s]   scaling  efficiency' +
          ('   vs-baseline' if base else ''))
    for nomp in nomps:
        s = t1/times[nomp]
        line = f'{nomp:4d} {times[nomp]:9.2f} {s:8.2f}x {s/nomp*100:9.0f}%'
        if base:
            line += f' {base[nomp]/times[nomp]:10.2f}x'
        print(line)
    if shutil.which('nproc'):
        ncpu = subprocess.run(['nproc'], capture_output=True,
                              text=True).stdout.strip()
        print(f'(host cores: {ncpu})')


if __name__ == '__main__':
    main()
