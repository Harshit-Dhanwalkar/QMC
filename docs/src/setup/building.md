# Building QMC

The top-level Makefile manages all build targets.

## Build Everything

To compile core libraries, physics modules, all examples, unit tests, and the demo driver

```bash
make all
```

## Build Targets

- `make demo`: Builds and launches the interactive terminal demo driver (`build/main`)
- `make examples`: Compiles all example binaries under `build/eg_*`.
- `make tests`: Compiles all unit tests under `build/test_*`.
- `make run-tests`: Builds and runs the full unit test suite with `pass/fail logging`.
- `make run-examples`: Executes all examples in sequence.
- `make clean`: Cleans the `build/` and `output/` directories.

## Makefile Build Flags

Pass configuration parameters directly to `make`

```make
# Force GNUplot backend instead of GR
make all PLOT_BACKEND=GNUPLOT

# Use custom GR installation prefix
make all GR_PREFIX=/usr/local/gr

# Use LuaLaTeX instead of pdflatex
make all LATEX_COMPILER=lualatex
```
