/*
 * Interactive (on-screen window) output

 * NOTE: this implements a single draw-and-refresh, not a
 * persistent interactive event loop (mouse/keyboard-driven redraw,
 * staying open until closed, etc.).
 * TODO: Look at GR's API
 */

#include "../gr_plot.h"
#include "../gr_plot_internal.h"

int gr_plot_interactive(const double *x, const double **ys, int n_series,
                        int n_pts, const char **labels,
                        const gr_plot_opt_t *opts) {
  double xmin, xmax, ymin, ymax;
  gr_compute_ranges(x, ys, n_series, n_pts, opts, &xmin, &xmax, &ymin, &ymax);

  gr_clearws();
  gr_draw_series(x, ys, n_series, n_pts, labels, opts, xmin, xmax, ymin, ymax);
  gr_updatews();

  return 0;
}
