#ifndef QMC_PLOT_H
#define QMC_PLOT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

typedef enum {
  PLOT_OK = 0,
  PLOT_ERR_INVALID_ARGUMENT = -1,
  PLOT_ERR_BACKEND_UNAVAILABLE = -2,
  PLOT_ERR_OPEN = -3,
  PLOT_ERR_WRITE = -4,
  PLOT_ERR_FORMAT_UNSUPPORTED = -5,
  PLOT_ERR_NOT_IMPLEMENTED = -6,
  PLOT_ERR_PATH_TOO_LONG = -7,
  PLOT_ERR_INTERNAL = -8
} plot_status_t;

typedef enum {
  PLOT_FORMAT_PNG,
  PLOT_FORMAT_PDF,
  PLOT_FORMAT_SVG,
  PLOT_FORMAT_EPS,
  PLOT_FORMAT_JPEG,
  PLOT_FORMAT_WINDOW, // interactive display
  PLOT_FORMAT_TEXT    // CSV-ish fallback
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
  int tex_text; // GR backend only: 1 = render for title/xlabel/ylabel/legend as
                // TeX math 0 = plain text (default)

  /*
   * Output directory prefix:
   *   NULL          -> use QMC_OUTPUT_DIR (default)
   *   ""            -> no prefix; use filename verbatim
   *   "/some/path"  -> prefix "/some/path/"
   */
  const char *output_dir;
} plot_opts_t;

/* Human-readable message for a status code */
const char *plot_strerror(plot_status_t status);

/* Name of compiled-in backend: "gr", "gnuplot", "matplotlib", "none" */
const char *plot_backend_name(void);

/*
 * Resolve full output path plot_line/plot_lines would write to
 *
 * Returns PLOT_ERR_FORMAT_UNSUPPORTED for PLOT_FORMAT_WINDOW (no file), or
 * PLOT_ERR_PATH_TOO_LONG if the result would not fit in `out`
 */
plot_status_t plot_output_path(const char *filename, plot_format_t format,
                               const plot_opts_t *opts, char *out,
                               size_t out_size);

/*
 * plot_line: plot y(x) to a file (or window)
 * 
 * Returns PLOT_OK on success, a negative plot_status_t on error
 */
int plot_line(const char *filename, plot_format_t format, const double *x,
              const double *y, int n, const plot_opts_t *opts);

/*
 * plot_lines: plot multiple series on one axes
 *  - ys[k]     : k-th series
 *  - labels[k] : its legend label (NULL = no legend)
 */
int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts);

#ifdef __cplusplus
}
#endif

#endif /* QMC_PLOT_H */
