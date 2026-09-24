/*
 * No-op plot backend
 *
 * PLOT_BACKEND=NONE
 *
 * All plot calls succeed silently without generating any output
 * For headless CI, benchmarking or not ploting backend is available
 * All other formats succeed silently without producing any file
 */

#include "plot.h"
#include "plot_internal.h"

const char *plot_backend_name(void) { return "none"; }

int plot_line(const char *filename, plot_format_t format, const double *x,
              const double *y, int n, const plot_opts_t *opts) {
  if (format == PLOT_FORMAT_TEXT) {
    char path[1024];
    plot_status_t st =
        plot_output_path(filename, format, opts, path, sizeof path);
    if (st != PLOT_OK) {
      return st;
    }

    return (int)plot_write_text_1d(path, x, y, n, opts);
  }

  (void)filename;
  (void)x;
  (void)y;
  (void)n;
  (void)opts;

  return PLOT_OK;
}

int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts) {
  if (format == PLOT_FORMAT_TEXT) {
    char path[1024];
    plot_status_t st =
        plot_output_path(filename, format, opts, path, sizeof path);
    if (st != PLOT_OK) {
      return st;
    }
    return (int)plot_write_text_nd(path, x, ys, n_series, n_pts, labels, opts);
  }

  (void)filename;
  (void)x;
  (void)ys;
  (void)n_series;
  (void)n_pts;
  (void)labels;
  (void)opts;

  return PLOT_OK;
}
