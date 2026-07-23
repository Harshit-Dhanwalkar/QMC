/*
 * GNU plot backend
 */

#include "gnuplot_pipe.h"
#include "../plot.h"
#include "/home/harshitpd/Documents/GITHUB/QMC/core/vector.h"
// #include <math.h>
#include <stdarg.h>
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

gnuplot_t *gnuplot_open(void) {
  gnuplot_t *gp = malloc(sizeof *gp);
  if (!gp)
    return NULL;

  gp->pipe = popen("gnuplot", "w");
  if (!gp->pipe) {
    free(gp);
    return NULL;
  }

  return gp;
}

void gnuplot_close(gnuplot_t *gp) {
  if (!gp)
    return;

  if (gp->pipe)
    pclose(gp->pipe);

  free(gp);
}

void gnuplot_cmd(gnuplot_t *gp, const char *cmd, ...) {
  if (!gp || !gp->pipe)
    return;

  va_list args;
  va_start(args, cmd);
  vfprintf(gp->pipe, cmd, args);
  va_end(args);
  fprintf(gp->pipe, "\n");
  fflush(gp->pipe);
}

// TODO:
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

  char path[512];
  if (format != PLOT_FORMAT_WINDOW)
    snprintf(path, sizeof path, "%s/%s%s", QMC_OUTPUT_DIR, filename,
             gp_ext(format));

  gnuplot_cmd(gp, "set terminal %s", gp_terminal(format));
  if (format != PLOT_FORMAT_WINDOW)
    gnuplot_cmd(gp, "set output '%s'", path);

  if (opts && opts->title)
    gnuplot_cmd(gp, "set title  '%s'", opts->title);
  if (opts && opts->xlabel)
    gnuplot_cmd(gp, "set xlabel '%s'", opts->xlabel);
  if (opts && opts->ylabel)
    gnuplot_cmd(gp, "set ylabel '%s'", opts->ylabel);
  gnuplot_cmd(gp, "set grid");

  double lw = (opts && opts->line_width > 0) ? opts->line_width : 2.0;
  if (opts && opts->color)
    gnuplot_cmd(gp, "plot '-' with lines lw %g lc rgb '%s' notitle", lw,
                opts->color);
  else
    gnuplot_cmd(gp, "plot '-' with lines lw %g notitle", lw);

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

  char path[512];
  if (format != PLOT_FORMAT_WINDOW)
    snprintf(path, sizeof path, "%s/%s%s", QMC_OUTPUT_DIR, filename,
             gp_ext(format));

  gnuplot_cmd(gp, "set terminal %s", gp_terminal(format));
  if (format != PLOT_FORMAT_WINDOW)
    gnuplot_cmd(gp, "set output '%s'", path);

  if (opts && opts->title)
    gnuplot_cmd(gp, "set title  '%s'", opts->title);
  if (opts && opts->xlabel)
    gnuplot_cmd(gp, "set xlabel '%s'", opts->xlabel);
  if (opts && opts->ylabel)
    gnuplot_cmd(gp, "set ylabel '%s'", opts->ylabel);
  gnuplot_cmd(gp, "set grid");
  gnuplot_cmd(gp, "set key outside right");

  double lw = (opts && opts->line_width > 0) ? opts->line_width : 2.0;
  for (int s = 0; s < n_series; s++) {
    char lbl[128];
    if (labels && labels[s])
      snprintf(lbl, sizeof lbl, "'%s'", labels[s]);
    else
      snprintf(lbl, sizeof lbl, "notitle");
    const char *cont = (s < n_series - 1) ? " \\" : "";
    if (s == 0)
      gnuplot_cmd(gp, "plot '-' with lines lw %g title %s%s", lw, lbl, cont);
    else
      gnuplot_cmd(gp, "  , '-' with lines lw %g title %s%s", lw, lbl, cont);
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

// Data savers
static FILE *open_output(const char *filename) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  FILE *f = fopen(path, "w");
  if (!f)
    fprintf(stderr, "plot: cannot open '%s'\n", path);
  return f;
}

void save_wavefunction(const char *filename, const double *x,
                       const cvector_t *psi, int n) {
  FILE *f = open_output(filename);
  if (!f)
    return;
  fprintf(f, "# x  Re(\\psi)  Im(\\psi)  |\\psi|^2\n");
  for (int i = 0; i < n; i++) {
    double re = psi->data[i].re, im = psi->data[i].im;
    fprintf(f, "%.8e  %.8e  %.8e  %.8e\n", x[i], re, im, re * re + im * im);
  }
  fclose(f);
}

void save_potential(const char *filename, const double *x, const double *V,
                    int n) {
  FILE *f = open_output(filename);
  if (!f)
    return;
  fprintf(f, "# x  V(x)\n");
  for (int i = 0; i < n; i++)
    fprintf(f, "%.8e  %.8e\n", x[i], V[i]);
  fclose(f);
}

void save_eigenvalues(const char *filename, const double *eigenvals, int n) {
  FILE *f = open_output(filename);
  if (!f)
    return;
  fprintf(f, "# n  E_n\n");
  for (int i = 0; i < n; i++)
    fprintf(f, "%d  %.10e\n", i, eigenvals[i]);
  fclose(f);
}
