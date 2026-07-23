# Python tools

The helpers ship as the **`hyperlite-tools`** package (`tools/`),
managed with [uv](https://docs.astral.sh/uv/) (lockfile committed):

```bash
cd tools/
uv sync                  # package + dev tools into tools/.venv
uv run sgfilter --help
```

or plain pip, from the repository root:

```bash
pip install -e tools/
```

All user-facing commands use [rich-click](https://github.com/ewels/rich-click)
styled `--help`.

## Console scripts

| Command | Purpose |
|---|---|
| `stella2superlite` | convert a Stella profile snapshot (`mesa.dayNNN_post_Lbol_max.data`) into a SuperLite `input.str` + `input.par.lte/.nlte` deck, with diagnostic profile plots |
| `parse-hdf5-sndata` | coarsen a FLASH-style block-AMR HDF5 snapshot to a block-averaged Cartesian structure file + profile plots |
| `sgfilter` | smooth a flux spectrum (`output.flx_luminos`) with a Savitzky-Golay filter |
| `wlgenadd` | generate an `input.wlgrid` custom wavelength grid (log base grid + high-resolution insert) |
| `h5-to-ascii` | export an `output.h5` back to the legacy ASCII `output.*` tables |
| `build-atomic-h5` | (re)generate `data/atomic.h5` from the raw ASCII atomic data — run automatically by the build |
| `hyperlite-bench` | wall-clock benchmark of the W7 smoke case across OpenMP thread counts (used by CI) |

`build-atomic-h5` and `hyperlite-bench` deliberately stay
argparse/stdlib-only so CMake and CI can invoke them as bare scripts
(`python3 tools/hyperlite_tools/...py`) without the package installed.

## Development

```bash
uv run ruff check .   # lint
uv run ty check       # type-check
```

Both run in CI on every push, along with `--help` smoke tests of every
console script and a `pip install -e tools/` check.
