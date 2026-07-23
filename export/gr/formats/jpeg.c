/*
 * JPEG output. Same screen-pixel default sizing as png.c.
 */

#include "../gr_plot_internal.h"
#include "../gr_plot.h"
#include <stdlib.h>

int gr_plot_jpeg(const char *filename, const double *x, const double **ys,
                 int n_series, int n_pts, const char **labels,
                 const gr_plot_opt_t *opts) {
  char path[512];
  gr_build_sized_path(path, sizeof path, filename, ".jpg", opts, 800, 600);
  setenv("GKS_WSTYPE", "nul", 0);

  return gr_emit_single_page(path, x, ys, n_series, n_pts, labels, opts);
}
