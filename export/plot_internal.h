#ifndef QMC_PLOT_INTERNAL_H
#define QMC_PLOT_INTERNAL_H

#include "plot.h"

/* Shared between plot backends' translation units but not part of public API:
 * keep out of libqmc.so's dynamic symbol table (tests/consumer/ symbols.sh
 * audits exported surface) */
#if defined(__GNUC__) || defined(__clang__)
#define QMC_INTERNAL __attribute__((visibility("hidden")))
#else
#define QMC_INTERNAL
#endif

/* Shared PLOT_FORMAT_TEXT writer used by every backend */
QMC_INTERNAL plot_status_t plot_write_text_1d(const char *path, const double *x,
                                              const double *y, int n,
                                              const plot_opts_t *opts);

QMC_INTERNAL plot_status_t plot_write_text_nd(const char *path, const double *x,
                                              const double **ys, int n_series,
                                              int n_pts, const char **labels,
                                              const plot_opts_t *opts);

#endif /* QMC_PLOT_INTERNAL_H */
