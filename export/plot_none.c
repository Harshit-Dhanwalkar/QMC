/*
 * No-op plot backend
 *
 * PLOT_BACKEND=NONE.
 *
 * All plot calls succeed silently without generating any output.
 * For headless CI, benchmarking, or when neither GR nor gnuplot
 * is available.
 */

#include "plot.h"
#include <stdio.h>

int plot_line(const char *filename, plot_format_t format, const double *x,
              const double *y, int n, const plot_opts_t *opts) {
  (void)filename;
  (void)format;
  (void)x;
  (void)y;
  (void)n;
  (void)opts;

  return 0; // silent no-op
}

int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts) {
  (void)filename;
  (void)format;
  (void)x;
  (void)ys;
  (void)n_series;
  (void)n_pts;
  (void)labels;
  (void)opts;

  return 0;
}
