# hyperlite-tools

Python tools for the hyperlite (SuperLite) IMC-DDMC radiation-transport code.

## Install

Managed with [uv](https://docs.astral.sh/uv/) (lockfile committed):

```bash
cd tools/
uv sync                 # install package + dev tools into tools/.venv
uv run sgfilter --help  # run a console script
```

Plain pip works too:

```bash
pip install -e tools/   # from the repository root
```

## Development

```bash
uv run ruff check .     # lint
uv run ty check         # type-check
```

Both run in CI on every push.

## Console scripts

| Command | Purpose |
|---|---|
| `stella2superlite` | Convert a Stella profile file to a SuperLite `input.str` deck |
| `sgfilter` | Smooth `output.flx_luminos` flux data with a Savitzky-Golay filter |
| `wlgenadd` | Generate an `input.wlgrid` (or `input.fluxwl`) custom wavelength grid |
| `parse-hdf5-sndata` | Coarsen HDF5 (FLASH-style block AMR) supernova data and convert to ASCII |

Each accepts `--help`. `hyperlite_tools.myfuncts` provides shared helper
functions (Stella profile loading, safe logs, file editing).

`stella2superlite.py` and `parse_hdf5_sndata.py` still run top-to-bottom at
module level (wrapped via `hyperlite_tools._run`); they will be refactored
into proper entry points in a later phase.
