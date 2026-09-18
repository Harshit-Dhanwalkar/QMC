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

---

## Building and Installing as a Library

QMC can also be built as a static library and installed system-wide (or to
any prefix) for use from other projects, discoverable via `pkg-config`:

```bash
make lib                          # build/libqmc.a
make install PREFIX=/usr/local    # lib + headers under include/qmc/ + qmc.pc
make uninstall PREFIX=/usr/local  # removes build
```

`PREFIX` defaults to `/usr/local`. The generated `qmc.pc` records whatever
`USE_LAPACK`/`PLOT_BACKEND` the install build was actually made with, so a
downstream project's own build gets the right link flags automatically:

```bash
pkg-config --cflags --libs qmc
gcc myprogram.c $(pkg-config --cflags --libs qmc) -fopenmp -o myprogram
```

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
