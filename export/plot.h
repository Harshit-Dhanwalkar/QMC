#ifndef QMC_PLOT_H
#define QMC_PLOT_H

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

typedef enum {
  PLOT_FORMAT_PNG,
  PLOT_FORMAT_PDF,
  PLOT_FORMAT_SVG,
  PLOT_FORMAT_EPS,
  PLOT_FORMAT_JPEG,
  PLOT_FORMAT_WINDOW, // interactive display
} plot_format_t;

typedef struct {
  int width, height; // pixels (PNG) or points (PDF/SVG)
  double xmin, xmax; // 0,0 = auto-scale
  double ymin, ymax;
  const char *title;
  const char *xlabel;
  const char *ylabel;
  const char *color; // NULL = default
  double line_width; // 0 = default
  int tex_text; // 1 = render (GR backend only) for title/xlabel/ylabel/legend
                // as TeX math 0 = plain text (default).
} plot_opts_t;

/*
 * plot_line: plot y(x) to a file (or window).
 * Returns 0 on success, -1 on error.
 */
int plot_line(const char *filename, plot_format_t format, const double *x,
              const double *y, int n, const plot_opts_t *opts);

/*
 * plot_lines: plot multiple series on one axes.
 * ys[k] is the k-th series.
 * labels[k] is the legend label (NULL = no legend).
 */
int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts);

#endif
