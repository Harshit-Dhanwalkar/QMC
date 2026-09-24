/*
 * GR framework backend
 * NOTE: Compiled with -DUSE_GR flag
 *
 * PLOT_FORMAT_TEXT is handled locally so that a GR-built binary still produces
 * CSV output for headless consumers.
 */

#include "gr/gr_plot.h"
#include "plot.h"
#include "plot_internal.h"

#include <string.h>

static gr_format_t to_gr_fmt(plot_format_t f) {
  switch (f) {
  case PLOT_FORMAT_PDF:
    return GR_FORMAT_PDF;
  case PLOT_FORMAT_SVG:
    return GR_FORMAT_SVG;
  case PLOT_FORMAT_PNG:
    return GR_FORMAT_PNG;
  case PLOT_FORMAT_JPEG:
    return GR_FORMAT_JPEG;
  case PLOT_FORMAT_EPS:
    return GR_FORMAT_EPS;
  case PLOT_FORMAT_WINDOW:
    return GR_FORMAT_WINDOW;
  default:
    return GR_FORMAT_PNG;
  }
}

static gr_plot_opt_t to_gr_opts(const plot_opts_t *o) {
  gr_plot_opt_t g = {0};
  if (!o) {
    return g;
  }

  g.xmin = o->xmin;
  g.xmax = o->xmax;
  g.ymin = o->ymin;
  g.ymax = o->ymax;
  g.title = o->title;
  g.xlabel = o->xlabel;
  g.ylabel = o->ylabel;
  g.width = o->width;
  g.height = o->height;
  g.line_width = o->line_width;
  g.color = o->color;
  g.tex_text = o->tex_text;

  return g;
}

const char *plot_backend_name(void) { return "gr"; }

int plot_line(const char *filename, plot_format_t format, const double *x,
              const double *y, int n, const plot_opts_t *opts) {
  if (!filename || (n > 0 && (!x || !y))) {
    return PLOT_ERR_INVALID_ARGUMENT;
  }

  if (format == PLOT_FORMAT_TEXT) {
    char path[1024];
    plot_status_t st =
        plot_output_path(filename, format, opts, path, sizeof path);
    if (st != PLOT_OK) {
      return st;
    }

    return (int)plot_write_text_1d(path, x, y, n, opts);
  }

  // WINDOW writes to a display, no path needed
  if (format == PLOT_FORMAT_WINDOW) {
    gr_plot_opt_t gr_opts = to_gr_opts(opts);

    return (int)gr_plot_to_file(filename, GR_FORMAT_WINDOW, x, y, n, &gr_opts);
  }

  char path[1024];
  plot_status_t st =
      plot_output_path(filename, format, opts, path, sizeof path);
  if (st != PLOT_OK) {
    return st;
  }

  gr_plot_opt_t gr_opts = to_gr_opts(opts);
  int rc = gr_plot_to_file(path, to_gr_fmt(format), x, y, n, &gr_opts);

  return rc == 0 ? PLOT_OK : PLOT_ERR_WRITE;
}

int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts) {
  if (!filename || n_series <= 0 || n_pts < 0 || (n_pts > 0 && !ys)) {
    return PLOT_ERR_INVALID_ARGUMENT;
  }

  if (format == PLOT_FORMAT_TEXT) {
    char path[1024];
    plot_status_t st =
        plot_output_path(filename, format, opts, path, sizeof path);
    if (st != PLOT_OK) {
      return st;
    }

    return (int)plot_write_text_nd(path, x, ys, n_series, n_pts, labels, opts);
  }

  if (format == PLOT_FORMAT_WINDOW) {
    gr_plot_opt_t gr_opts = to_gr_opts(opts);

    return (int)gr_plot_lines(filename, GR_FORMAT_WINDOW, x, ys, n_series,
                              n_pts, labels, &gr_opts);
  }

  char path[1024];
  plot_status_t st =
      plot_output_path(filename, format, opts, path, sizeof path);
  if (st != PLOT_OK) {
    return st;
  }

  gr_plot_opt_t gr_opts = to_gr_opts(opts);
  int rc = gr_plot_lines(path, to_gr_fmt(format), x, ys, n_series, n_pts,

                         labels, &gr_opts);
  return rc == 0 ? PLOT_OK : PLOT_ERR_WRITE;
}
