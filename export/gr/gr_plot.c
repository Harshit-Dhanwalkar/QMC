/*
 * GR framework backend: shared plotting core + format dispatch.
 */

#include "gr_plot.h"
#include "../../third_party/gr/install/include/gks.h"
#include "../../third_party/gr/install/include/gr.h"
#include "gr_plot_internal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int gr_colors[] = {4, 2, 3, 5, 6, 7, 8, 9};
const int gr_n_colors = (int)(sizeof(gr_colors) / sizeof(gr_colors[0]));

static int gr_color_from_name(const char *name) {
  if (!name) {
    return -1;
  }

  static const struct {
    const char *name;
    int idx;
  } table[] = {
      {"black", 1}, {"red", 2},    {"green", 3}, {"blue", 4}, {"magenta", 5},
      {"cyan", 6},  {"yellow", 7}, {"gray", 8},  {"grey", 8}, {"orange", 9},
  };

  for (size_t i = 0; i < sizeof(table) / sizeof(table[0]); i++) {
    if (strcmp(name, table[i].name) == 0) {
      return table[i].idx;
    }
  }

  return -1;
}

void gr_auto_range(const double *v, int n, double *lo, double *hi,
                   double pad_frac) {
  *lo = *hi = v[0];
  for (int i = 1; i < n; i++) {
    if (v[i] < *lo) {
      *lo = v[i];
    }

    if (v[i] > *hi) {
      *hi = v[i];
    }
  }

  double pad = (*hi - *lo) * pad_frac;
  if (pad < 1e-12) {
    pad = 1.0;
  }

  *lo -= pad;
  *hi += pad;
}

void gr_setup_viewport_and_window(double xmin, double xmax, double ymin,
                                  double ymax) {
  gr_setviewport(0.12, 0.93, 0.12, 0.90);
  gr_setwindow(xmin, xmax, ymin, ymax);
}

void gr_draw_labels(const gr_plot_opt_t *opts) {
  if (!opts) {
    return;
  }

  double vp[4], wn[4];
  gr_inqviewport(&vp[0], &vp[1], &vp[2], &vp[3]);
  gr_inqwindow(&wn[0], &wn[1], &wn[2], &wn[3]);

  gr_setviewport(0.0, 1.0, 0.0, 1.0);
  gr_setwindow(0.0, 1.0, 0.0, 1.0);

  gr_settextfontprec(232, 3);
  gr_settextcolorind(1);
  gr_setcharup(0.0, 1.0);

  int tex = opts->tex_text;

  if (opts->title && opts->title[0]) {
    gr_setcharheight(0.030);
    gr_settextalign(GKS_K_TEXT_HALIGN_CENTER, GKS_K_TEXT_VALIGN_TOP);
    if (tex) {
      char buf[160];
      snprintf(buf, sizeof buf, "$%s$", opts->title);

      gr_text(0.5, 0.98, buf);
    } else {
      gr_text(0.5, 0.98, (char *)opts->title);
    }
  }

  if (opts->xlabel && opts->xlabel[0]) {
    gr_setcharheight(0.022);
    gr_settextalign(GKS_K_TEXT_HALIGN_CENTER, GKS_K_TEXT_VALIGN_BOTTOM);
    if (tex) {
      char buf[160];
      snprintf(buf, sizeof buf, "$%s$", opts->xlabel);

      gr_text(0.5, 0.02, buf);
    } else {
      gr_text(0.5, 0.02, (char *)opts->xlabel);
    }
  }

  if (opts->ylabel && opts->ylabel[0]) {
    gr_setcharheight(0.022);
    gr_setcharup(-1.0, 0.0);
    gr_settextalign(GKS_K_TEXT_HALIGN_CENTER, GKS_K_TEXT_VALIGN_TOP);

    if (tex) {
      char buf[160];
      snprintf(buf, sizeof buf, "$%s$", opts->ylabel);

      gr_text(0.02, 0.5, buf);
    } else {
      gr_text(0.02, 0.5, (char *)opts->ylabel);
    }

    gr_setcharup(0.0, 1.0);
  }

  gr_setviewport(vp[0], vp[1], vp[2], vp[3]);
  gr_setwindow(wn[0], wn[1], wn[2], wn[3]);
}

void gr_compute_ranges(const double *x, const double **ys, int n_series,
                       int n_pts, const gr_plot_opt_t *opts, double *xmin,
                       double *xmax, double *ymin, double *ymax) {
  gr_auto_range(x, n_pts, xmin, xmax, 0.05);
  if (opts && opts->xmin < opts->xmax) {
    *xmin = opts->xmin;
    *xmax = opts->xmax;
  }

  *ymin = ys[0][0];
  *ymax = ys[0][0];
  for (int s = 0; s < n_series; s++) {
    for (int i = 0; i < n_pts; i++) {
      if (ys[s][i] < *ymin) {
        *ymin = ys[s][i];
      }

      if (ys[s][i] > *ymax) {
        *ymax = ys[s][i];
      }
    }
  }

  double pad = (*ymax - *ymin) * 0.08;
  if (pad < 1e-12) {
    pad = 1.0;
  }

  *ymin -= pad;
  *ymax += pad;
  if (opts && opts->ymin < opts->ymax) {
    *ymin = opts->ymin;
    *ymax = opts->ymax;
  }
}

void gr_draw_series(const double *x, const double **ys, int n_series, int n_pts,
                    const char **labels, const gr_plot_opt_t *opts, double xmin,
                    double xmax, double ymin, double ymax) {
  gr_setup_viewport_and_window(xmin, xmax, ymin, ymax);
  gr_setlinewidth(1.0);
  gr_setlinecolorind(1);

  double xtick = gr_tick(xmin, xmax);
  double ytick = gr_tick(ymin, ymax);

  double x_org = (xmin <= 0.0 && 0.0 <= xmax) ? 0.0 : xmin;
  double y_org = (ymin <= 0.0 && 0.0 <= ymax) ? 0.0 : ymin;

  gr_axes(xtick, ytick, x_org, y_org, 1, 1, -0.01);
  gr_grid(xtick, ytick, x_org, y_org, 1, 1);

  int base_color = 4;
  if (opts && opts->color) {
    int mapped = gr_color_from_name(opts->color);

    if (mapped >= 0) {
      base_color = mapped;
    }
  }

  double lw = (opts && opts->line_width > 0) ? opts->line_width : 2.0;
  for (int s = 0; s < n_series; s++) {
    gr_setlinewidth(lw);
    gr_setlinecolorind(n_series == 1 ? base_color : gr_colors[s % gr_n_colors]);
    gr_polyline(n_pts, (double *)x, (double *)ys[s]);
  }

  if (labels && n_series > 1) {
    gr_settextfontprec(232, 3);
    gr_setcharheight(0.018);
    int tex = opts && opts->tex_text;

    for (int s = 0; s < n_series && labels[s]; s++) {
      gr_setlinecolorind(gr_colors[s % gr_n_colors]);
      double lx = xmax - (xmax - xmin) * 0.28;
      double ly = ymax - (ymax - ymin) * (0.05 + s * 0.07);
      double leg_x[2] = {lx, lx + (xmax - xmin) * 0.06};
      double leg_y[2] = {ly, ly};
      gr_polyline(2, leg_x, leg_y);
      gr_settextcolorind(1);
      gr_settextalign(GKS_K_TEXT_HALIGN_LEFT, GKS_K_TEXT_VALIGN_HALF);

      char lbl[128];
      if (tex) {
        snprintf(lbl, sizeof lbl, "$%s$", labels[s]);
      } else {
        snprintf(lbl, sizeof lbl, "%s", labels[s]);
      }

      gr_text(lx + (xmax - xmin) * 0.07, ly, lbl);
    }
  }

  gr_draw_labels(opts);
}

void gr_build_sized_path(char *path, size_t path_size, const char *filename,
                         const char *ext, const gr_plot_opt_t *opts,
                         int default_w, int default_h) {
  int w = (opts && opts->width > 0) ? opts->width : default_w;
  int h = (opts && opts->height > 0) ? opts->height : default_h;

  (void)opts;
  (void)default_w;
  (void)default_h;
  snprintf(path, path_size, "%s/%s%s", QMC_OUTPUT_DIR, filename, ext);

  // TEST: For EPS, GR may not support the |WxH suffix; test it without suffix
  // for EPS.
  // snprintf(path, path_size, "%s/%s%s|x%d", QMC_OUTPUT_DIR, filename, ext, w,
  // h);
}

int gr_emit_single_page(const char *path, const double *x, const double **ys,
                        int n_series, int n_pts, const char **labels,
                        const gr_plot_opt_t *opts) {
  double xmin, xmax, ymin, ymax;
  gr_compute_ranges(x, ys, n_series, n_pts, opts, &xmin, &xmax, &ymin, &ymax);

  gr_beginprint((char *)path);
  gr_draw_series(x, ys, n_series, n_pts, labels, opts, xmin, xmax, ymin, ymax);
  gr_endprint();

  return 0;
}

int gr_plot_png(const char *filename, const double *x, const double **ys,
                int n_series, int n_pts, const char **labels,
                const gr_plot_opt_t *opts);

int gr_plot_jpeg(const char *filename, const double *x, const double **ys,
                 int n_series, int n_pts, const char **labels,
                 const gr_plot_opt_t *opts);

int gr_plot_svg(const char *filename, const double *x, const double **ys,
                int n_series, int n_pts, const char **labels,
                const gr_plot_opt_t *opts);

int gr_plot_pdf(const char *filename, const double *x, const double **ys,
                int n_series, int n_pts, const char **labels,
                const gr_plot_opt_t *opts);

int gr_plot_interactive(const double *x, const double **ys, int n_series,
                        int n_pts, const char **labels,
                        const gr_plot_opt_t *opts);

static int gr_plot_eps(const char *filename, const double *x, const double **ys,
                       int n_series, int n_pts, const char **labels,
                       const gr_plot_opt_t *opts) {
  char path[512];
  gr_build_sized_path(path, sizeof path, filename, ".eps", opts, 1000, 750);
  setenv("GKS_WSTYPE", "nul", 0);

  return gr_emit_single_page(path, x, ys, n_series, n_pts, labels, opts);
}

static int dispatch(const char *filename, gr_format_t format, const double *x,
                    const double **ys, int n_series, int n_pts,
                    const char **labels, const gr_plot_opt_t *opts) {
  switch (format) {
  case GR_FORMAT_PNG:
    return gr_plot_png(filename, x, ys, n_series, n_pts, labels, opts);
  case GR_FORMAT_JPEG:
    return gr_plot_jpeg(filename, x, ys, n_series, n_pts, labels, opts);
  case GR_FORMAT_SVG:
    return gr_plot_svg(filename, x, ys, n_series, n_pts, labels, opts);
  case GR_FORMAT_PDF:
    return gr_plot_pdf(filename, x, ys, n_series, n_pts, labels, opts);
  case GR_FORMAT_EPS:
    return gr_plot_eps(filename, x, ys, n_series, n_pts, labels, opts);
  case GR_FORMAT_WINDOW:
    return gr_plot_interactive(x, ys, n_series, n_pts, labels, opts);
  default:
    return gr_plot_png(filename, x, ys, n_series, n_pts, labels, opts);
  }
}

int gr_plot_to_file(const char *filename, gr_format_t format, const double *x,
                    const double *y, int n, const gr_plot_opt_t *opts) {
  if (!filename || !x || !y || n < 2) {
    return -1;
  }

  const double *ys[1] = {y};

  return dispatch(filename, format, x, ys, 1, n, NULL, opts);
}

int gr_plot_lines(const char *filename, gr_format_t format, const double *x,
                  const double **ys, int n_series, int n_pts,
                  const char **labels, const gr_plot_opt_t *opts) {
  if (!filename || !x || !ys || n_series < 1 || n_pts < 2) {
    return -1;
  }

  return dispatch(filename, format, x, ys, n_series, n_pts, labels, opts);
}
