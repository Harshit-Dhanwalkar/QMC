#ifndef QMC_GR_PLOT_INTERNAL_H
#define QMC_GR_PLOT_INTERNAL_H

#include "gr_plot.h"
#include <stddef.h>
#include "../../third_party/gr/install/include/gr.h"
#include "../../third_party/gr/install/include/gr3.h"
#include "../../third_party/gr/install/include/gks.h"


#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

extern const int gr_colors[];
extern const int gr_n_colors;

/* Auto-scale range with padding */
void gr_auto_range(const double *v, int n, double *lo, double *hi,
                   double pad_frac);

void gr_setup_viewport_and_window(double xmin, double xmax, double ymin,
                                  double ymax);
void gr_draw_labels(const gr_plot_opt_t *opts);

/*
 * Draws axes + grid + one-or-more polylines (+ legend if n_series>1
 * and labels given)
 */
void gr_draw_series(const double *x, const double **ys, int n_series, int n_pts,
                    const char **labels, const gr_plot_opt_t *opts, double xmin,
                    double xmax, double ymin, double ymax);

/*
 * Computes combined axis range for n_series data series
 */
void gr_compute_ranges(const double *x, const double **ys, int n_series,
                       int n_pts, const gr_plot_opt_t *opts, double *xmin,
                       double *xmax, double *ymin, double *ymax);

/*
 * Builds "<QMC_OUTPUT_DIR>/<filename><ext>|<W>x<H>" - QMC_OUTPUT_DIR
 * prefix plus GR's own "|WxH" canvas-size suffix convention on
 * gr_beginprint.
 */
void gr_build_sized_path(char *path, size_t path_size, const char *filename,
                         const char *ext, const gr_plot_opt_t *opts,
                         int default_w, int default_h);

/*
 * Generic single-page file emitter: opens path, draws via
 * gr_draw_series, closes.
 */
int gr_emit_single_page(const char *path, const double *x, const double **ys,
                        int n_series, int n_pts, const char **labels,
                        const gr_plot_opt_t *opts);

#endif
