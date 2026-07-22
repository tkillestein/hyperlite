#!/bin/bash
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
# Stage and run the W7 LTE smoke case.
#   smoke_run.sh EXE SRCDIR RUNDIR [LAUNCHER...]
# LAUNCHER (e.g. "mpirun -n 2") is prepended to the executable if given.
set -e

EXE=${1:?usage: smoke_run.sh EXE SRCDIR RUNDIR [LAUNCHER...]}
SRC=${2:?missing SRCDIR}
RUNDIR=${3:?missing RUNDIR}
shift 3

rm -rf "$RUNDIR"
mkdir -p "$RUNDIR"
cd "$RUNDIR"
ln -s "$SRC"/Data/* .
ln -s "$SRC"/Input/* .
ln -sf input_w7.str input.str
cp "$SRC"/Input/input.par.lte input.par

"$@" "$EXE" | tee stdout.log
grep -q "SuperLite finished" stdout.log
