#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""Console-script shims for the legacy top-to-bottom scripts.

stella2superlite.py and parse_hdf5_sndata.py execute at module level (no
main() function yet); runpy re-executes them as __main__ so they keep
working both as console scripts and via `python -m`. To be refactored
into proper entry points in a later phase.
"""
import runpy


def _run(module):
    runpy.run_module(f"hyperlite_tools.{module}", run_name="__main__")


def stella2superlite():
    _run("stella2superlite")


def parse_hdf5_sndata():
    _run("parse_hdf5_sndata")
