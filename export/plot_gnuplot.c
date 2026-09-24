/*
 * GNUplot pipe backend
 * Compiled with -DUSE_GNUPLOT flag
 *
 * Requires gnuplot on PATH
 */

#include "plot.h"
#include "gnuplot/gnuplot_pipe.h"
#include "plot_internal.h"

#include <stdio.h>
#include <stdlib.h>

static int gnuplot_available(void) {
  return system("command -v gnuplot > /dev/null 2>&1") == 0;
}

static const char *gp_terminal(plot_format_t fmt) {
  switch (fmt) {
  case PLOT_FORMAT_PNG:
    return "pngcairo enhanced font 'Sans,11' size 800,600";
  case PLOT_FORMAT_PDF:
    return "pdfcairo enhanced font 'Sans,11'";
  case PLOT_FORMAT_SVG:
    return "svg enhanced font 'Sans,11'";
  case PLOT_FORMAT_EPS:
    return "epscairo enhanced font 'Sans,11'";
  case PLOT_FORMAT_JPEG:
    return "jpeg size 800,600";
  case PLOT_FORMAT_WINDOW:
    return "wxt persist";
  default:
    return "pngcairo enhanced";
  }
}

const char *plot_backend_name(void) { return "gnuplot"; }

static void apply_opts(gnuplot_t *gp, const plot_opts_t *opts) {
  if (!opts) {
    return;
  }
  if (opts->title) {
    gnuplot_cmd(gp, "set title  '%s'", opts->title);
  }
  if (opts->xlabel) {
    gnuplot_cmd(gp, "set xlabel '%s'", opts->xlabel);
  }
  if (opts->ylabel) {
    gnuplot_cmd(gp, "set ylabel '%s'", opts->ylabel);
  }
  if (opts->xmin != 0.0 || opts->xmax != 0.0) {
    gnuplot_cmd(gp, "set xrange [%g:%g]", opts->xmin, opts->xmax);
  }
  if (opts->ymin != 0.0 || opts->ymax != 0.0) {
    gnuplot_cmd(gp, "set yrange [%g:%g]", opts->ymin, opts->ymax);
  }
}

static const char *gp_ext(plot_format_t fmt) {
  switch (fmt) {
  case PLOT_FORMAT_PNG:
    return ".png";
  case PLOT_FORMAT_PDF:
    return ".pdf";
  case PLOT_FORMAT_SVG:
    return ".svg";
  case PLOT_FORMAT_EPS:
    return ".eps";
  case PLOT_FORMAT_JPEG:
    return ".jpg";
  default:
    return ".png";
  }
}

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

  if (!gnuplot_available()) {
    fprintf(stderr,
            "plot(gnuplot): gnuplot not found - skipping '%s'\n"
            "      Install: sudo apt install gnuplot\n",
            filename);

    return PLOT_ERR_BACKEND_UNAVAILABLE;
  }

  gnuplot_t *gp = gnuplot_open();
  if (!gp) {
    return PLOT_ERR_BACKEND_UNAVAILABLE;
  }

  gnuplot_cmd(gp, "set terminal %s", gp_terminal(format));

  if (format != PLOT_FORMAT_WINDOW) {
    char path[1024];
    plot_status_t st =
        plot_output_path(filename, format, opts, path, sizeof path);
    if (st != PLOT_OK) {
      gnuplot_close(gp);

      return st;
    }

    gnuplot_cmd(gp, "set output '%s'", path);
  }

  apply_opts(gp, opts);
  gnuplot_cmd(gp, "set grid");
  gnuplot_cmd(gp, "plot '-' with lines lw 2 notitle");
  for (int i = 0; i < n; i++) {
    gnuplot_cmd(gp, "%.10e %.10e", x[i], y[i]);
  }
  gnuplot_cmd(gp, "e");

  if (format != PLOT_FORMAT_WINDOW) {
    gnuplot_cmd(gp, "set output");
  }

  gnuplot_close(gp);

  return PLOT_OK;
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

  if (!gnuplot_available()) {
    fprintf(stderr, "plot(gnuplot): gnuplot not found - skipping '%s'\n",
            filename);

    return PLOT_ERR_BACKEND_UNAVAILABLE;
  }

  gnuplot_t *gp = gnuplot_open();
  if (!gp) {
    return PLOT_ERR_BACKEND_UNAVAILABLE;
  }

  gnuplot_cmd(gp, "set terminal %s", gp_terminal(format));

  if (format != PLOT_FORMAT_WINDOW) {
    char path[1024];
    plot_status_t st =
        plot_output_path(filename, format, opts, path, sizeof path);
    if (st != PLOT_OK) {
      gnuplot_close(gp);

      return st;
    }

    gnuplot_cmd(gp, "set output '%s'", path);
  }

  apply_opts(gp, opts);
  gnuplot_cmd(gp, "set grid");
  gnuplot_cmd(gp, "set key outside right");

  for (int s = 0; s < n_series; s++) {
    const char *lbl = (labels && labels[s]) ? labels[s] : "";
    if (s == 0) {
      gnuplot_cmd(gp, "plot '-' with lines lw 2 title '%s'%s", lbl,
                  (n_series > 1) ? " \\" : "");
    } else {
      gnuplot_cmd(gp, "  , '-' with lines lw 2 title '%s'%s", lbl,
                  (s < n_series - 1) ? " \\" : "");
    }
  }

  for (int s = 0; s < n_series; s++) {
    for (int i = 0; i < n_pts; i++) {
      gnuplot_cmd(gp, "%.10e %.10e", x[i], ys[s][i]);
    }

    gnuplot_cmd(gp, "e");
  }

  if (format != PLOT_FORMAT_WINDOW) {
    gnuplot_cmd(gp, "set output");
  }

  gnuplot_close(gp);

  return PLOT_OK;
}
