# Plotting & Data Export

QMC separates _data export_ (writing raw results to disk) from _plotting_ (rendering a figure), and further separates plotting itself into swappable backends selected at compile time via `PLOT_BACKEND=`.

> **Note:** `export/plot.h` itself hasn't been reviewed yet - everything below about `plot_opts_t`'s exact fields and `plot_format_t`'s exact enum values is reconstructed from how the four backend `.c` files use them, not read directly from the header. The shared interface (`plot_line`, `plot_lines`) is solid - all four backends implement identical signatures - but treat field names as "consistent with observed usage" rather than "verified against the header" until `plot.h` is reviewed directly.

## Shared interface

Every backend implements the same two entry points:

```c
int plot_line(const char *filename, plot_format_t format, const double *x,
const double *y, int n, const plot_opts_t *opts);

int plot_lines(const char *filename, plot_format_t format, const double *x,
const double **ys, int n_series, int n_pts, const char **labels,
const plot_opts_t \*opts);
```

`plot_opts_t` (`export/plot.h`):

```c
typedef struct {
int width, height;
double xmin, xmax;
double ymin, ymax;
const char *title;
const char *xlabel;
const char *ylabel;
const char *color; // lowercase name ("red","blue",...); GR backend
                   // only, maps internally to a palette index
double line_width;
int tex_labels;    // 1 = render title/xlabel/ylabel/legend as TeX
                   // math (GR backend only; gnuplot/matplotlib
                   // always render plain text). Default 0.
} plot_opts_t;
```

`plot_format_t`: `PLOT_FORMAT_PNG`, `PDF`, `SVG`, `EPS`, `JPEG`, `WINDOW` (interactive display, not a file).

## Backends

Selected at compile time (`PLOT_BACKEND=GR|GNUPLOT|MATPLOTLIB|NONE`); only one is linked into a given build.

### GR (`plot_gr.c`, `-DUSE_GR`)

Thin translation layer onto a `gr/gr_plot.h` wrapper (`gr_plot_to_file`, `gr_plot_lines`) - the actual GR calls live in that wrapper, not here. This is the primary backend per project conventions (built locally with Qt disabled, `GKS_WSTYPE=nul` for headless file output).

### Gnuplot (`plot_gnuplot.c`, `-DUSE_GNUPLOT`)

Drives `gnuplot` as a subprocess via a pipe wrapper (`gnuplot/gnuplot_pipe.h`). Checks `command -v gnuplot` before doing anything and fails gracefully with an install hint (`sudo apt install gnuplot`) if it's missing, rather than crashing. Builds one `pngcairo`/`pdfcairo`/`svg`/`epscairo`/`jpeg`/`wxt` terminal string per format, writes data inline via repeated `gnuplot_cmd(gp, "%.10e %.10e", ...)` calls terminated with `e`, and supports multi-series plots with a `set key outside right` legend. This is what CI (`ci.yml`) actually exercises, since GR isn't installed there.

### Matplotlib (`plot_matplotlib.c`, `PLOT_BACKEND=MATPLOTLIB`)

> Drives a `python3` subprocess (`matplotlib/matplotlib_pipe.h`), building up a Python script line-by-line via `matplotlib_cmd`. Checks both that `python3` is on `PATH` and that the `matplotlib` module actually imports before proceeding, with distinct error messages for each failure mode (`pip install matplotlib` hint on the latter). Includes a `py_repr()` helper that escapes backslashes and single quotes to safely build Python string literals from C strings - worth knowing about if a title/label ever contains one of those characters.

> **Known gap:** `PLOT_FORMAT_WINDOW` (interactive display) is explicitly
> unimplemented for this backend - both `plot_line` and `plot_lines` return
> `-1` immediately with a stderr message for that format, rather than
> silently no-op'ing.

### None (`plot_none.c`, `PLOT_BACKEND=NONE`)

Silent no-op - every call returns `0` immediately without touching the filesystem. For headless CI, benchmarking, or environments with none of the above available.

## Data export

`export/` also holds raw data writers, independent of the plotting backend chosen:

```c
// CSV, with a header row
int csv_write_1d(const char *filename, const double *x, const double *y, int n,
                 const char *xlabel, const char *ylabel);
int csv_write_matrix(const char *filename, const double *data, int rows,
                     int cols, const char **col_headers);
```

Both write into `QMC_OUTPUT_DIR` (default `"output"`), 6-significant-figure scientific notation (`% .6e`). This is the format `.dat` files referenced throughout the physics pages (e.g. `harmonic_energies.dat`) presumably come from.

> **Not yet implemented:**
>
> - `hdf5_write_1d` (`hdf5_writer.c`) is a stub - the real HDF5 implementation is commented out in full, and the active function just prints `"HDF5 not is not supported yet."` to stderr and returns `-1`.
> - `netcdf_writer.c` is completely empty (just a `TODO` comment) - no NetCDF export exists at all yet, despite `CMakeLists.txt`/`Makefile` presumably listing it as a source file.
