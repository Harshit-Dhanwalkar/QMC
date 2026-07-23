/*
 * GR framework backend
 * Compiled when -DUSE_GR is set.
 */

#include "gr/gr_plot.h"
#include "plot.h"

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
  if (!o)
    return g;

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
  g.tex_labels = 1;

  return g;
}

int plot_line(const char *filename, plot_format_t format, const double *x,
              const double *y, int n, const plot_opts_t *opts) {
  gr_format_t gr_fmt = to_gr_fmt(format);
  gr_plot_opt_t gr_opts = to_gr_opts(opts);

  return gr_plot_to_file(filename, gr_fmt, x, y, n, &gr_opts);
}

int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts) {
  gr_format_t gr_fmt = to_gr_fmt(format);
  gr_plot_opt_t gr_opts = to_gr_opts(opts);

  return gr_plot_lines(filename, gr_fmt, x, ys, n_series, n_pts, labels,
                       &gr_opts);
}
