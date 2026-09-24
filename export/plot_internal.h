#ifndef QMC_PLOT_INTERNAL_H
#define QMC_PLOT_INTERNAL_H

#include "plot.h"

/* Shared PLOT_FORMAT_TEXT writer used by every backend */
plot_status_t plot_write_text_1d(const char *path, const double *x,
                                 const double *y, int n,
                                 const plot_opts_t *opts);

plot_status_t plot_write_text_nd(const char *path, const double *x,
                                 const double **ys, int n_series, int n_pts,
                                 const char **labels, const plot_opts_t *opts);

#endif /* QMC_PLOT_INTERNAL_H */
