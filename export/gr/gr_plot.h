#ifndef QMC_GR_PLOT_H
#define QMC_GR_PLOT_H

typedef enum {
  GR_FORMAT_PNG,
  GR_FORMAT_PDF,
  GR_FORMAT_SVG,
  GR_FORMAT_EPS,
  GR_FORMAT_JPEG,
  GR_FORMAT_WINDOW,
} gr_format_t;

typedef struct {
  int width, height;
  double xmin, xmax;
  double ymin, ymax;
  const char *title;
  const char *xlabel;
  const char *ylabel;
  const char *color; // lowercase color name.
                     // NOTE: Only used for single-series plots; ignored/falls
                     // back to default if NULL or unrecognized.
  double line_width;
  int tex_text; // 1 = render labels as TeX math
} gr_plot_opt_t;

/*
 * Plot y(x) using GR framework.
 */
int gr_plot_to_file(const char *filename, gr_format_t format, const double *x,
                    const double *y, int n, const gr_plot_opt_t *opts);

/*
 * Plot multiple series.
 * ys[k] = k-th series data, labels[k] = legend entry or NULL.
 */
int gr_plot_lines(const char *filename, gr_format_t format, const double *x,
                  const double **ys, int n_series, int n_pts,
                  const char **labels, const gr_plot_opt_t *opts);

#endif /* QMC_GR_PLOT_H */
