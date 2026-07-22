#!/bin/bash
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
# Smoke test: run the W7 LTE case and check it against the committed
# statistical reference ensemble (tests/regression/reference/smoke.h5).
# Called from the top-level Makefile as: test.sh $RUNDIR
set -e

RUNDIR=${1:?usage: test.sh RUNDIR}
TOP=$(cd "$(dirname "$0")"/../.. && pwd)

cd "$RUNDIR"
./superlite | tee stdout.log
grep -q "SuperLite finished" stdout.log

REF="$TOP/tests/regression/reference/smoke.h5"
if [[ -f $REF ]]; then
  python3 "$TOP/tests/regression/compare.py" "$RUNDIR" "$REF" -v
else
  echo "WARNING: no reference ensemble at $REF; smoke-only (run completed)"
fi
