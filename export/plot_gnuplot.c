/*
 * gnuplot pipe backend
 * Compiled when -DUSE_GNUPLOT is set.
 *
 * Requires gnuplot on PATH: sudo apt install gnuplot
 */

#include "plot.h"
#include "gnuplot/gnuplot_pipe.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  if (!gnuplot_available()) {
    fprintf(stderr,
            "plot: gnuplot not found - skipping '%s'\n"
            "      Install: sudo apt install gnuplot\n",
            filename);
    return -1;
  }

  gnuplot_t *gp = gnuplot_open();
  if (!gp)
    return -1;

  gnuplot_cmd(gp, "set terminal %s", gp_terminal(format));

  if (format != PLOT_FORMAT_WINDOW) {
    char path[512];
    snprintf(path, sizeof path, "%s/%s%s", QMC_OUTPUT_DIR, filename,
             gp_ext(format));
    gnuplot_cmd(gp, "set output '%s'", path);
  }

  if (opts) {
    if (opts->title)
      gnuplot_cmd(gp, "set title  '%s'", opts->title);
    if (opts->xlabel)
      gnuplot_cmd(gp, "set xlabel '%s'", opts->xlabel);
    if (opts->ylabel)
      gnuplot_cmd(gp, "set ylabel '%s'", opts->ylabel);
  }
  gnuplot_cmd(gp, "set grid");
  gnuplot_cmd(gp, "plot '-' with lines lw 2 notitle");
  for (int i = 0; i < n; i++)
    gnuplot_cmd(gp, "%.10e %.10e", x[i], y[i]);
  gnuplot_cmd(gp, "e");

  if (format != PLOT_FORMAT_WINDOW)
    gnuplot_cmd(gp, "set output");

  gnuplot_close(gp);
  return 0;
}

int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts) {
  if (!gnuplot_available()) {
    fprintf(stderr, "plot: gnuplot not found - skipping '%s'\n", filename);
    return -1;
  }

  gnuplot_t *gp = gnuplot_open();
  if (!gp)
    return -1;

  gnuplot_cmd(gp, "set terminal %s", gp_terminal(format));

  if (format != PLOT_FORMAT_WINDOW) {
    char path[512];
    snprintf(path, sizeof path, "%s/%s%s", QMC_OUTPUT_DIR, filename,
             gp_ext(format));
    gnuplot_cmd(gp, "set output '%s'", path);
  }

  if (opts) {
    if (opts->title)
      gnuplot_cmd(gp, "set title  '%s'", opts->title);
    if (opts->xlabel)
      gnuplot_cmd(gp, "set xlabel '%s'", opts->xlabel);
    if (opts->ylabel)
      gnuplot_cmd(gp, "set ylabel '%s'", opts->ylabel);
  }
  gnuplot_cmd(gp, "set grid");
  gnuplot_cmd(gp, "set key outside right");

  /* Build inline plot command */
  for (int s = 0; s < n_series; s++) {
    const char *lbl = (labels && labels[s]) ? labels[s] : "";
    if (s == 0)
      gnuplot_cmd(gp, "plot '-' with lines lw 2 title '%s'%s", lbl,
                  (n_series > 1) ? " \\" : "");
    else
      gnuplot_cmd(gp, "  , '-' with lines lw 2 title '%s'%s", lbl,
                  (s < n_series - 1) ? " \\" : "");
  }
  for (int s = 0; s < n_series; s++) {
    for (int i = 0; i < n_pts; i++)
      gnuplot_cmd(gp, "%.10e %.10e", x[i], ys[s][i]);
    gnuplot_cmd(gp, "e");
  }

  if (format != PLOT_FORMAT_WINDOW)
    gnuplot_cmd(gp, "set output");

  gnuplot_close(gp);
  return 0;
}
