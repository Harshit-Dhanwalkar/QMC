/*
 * Matplotlib (Python subprocess pipe) backend for PLOT_BACKEND=MATPLOTLIB.
 *
 * TODO: Implement PLOT_FORMAT_WINDOW for matplotlib's interactive
 * plt.show() GUI backend (TkAgg/Qt5Agg)
 *
 * Returns -1 for WINDOW mode
 */

#include "matplotlib/matplotlib_pipe.h"
#include "plot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int python3_available(void) {
  return system("command -v python3 > /dev/null 2>&1") == 0;
}

static int matplotlib_module_available(void) {
  return system("python3 -c 'import matplotlib' > /dev/null 2>&1") == 0;
}

static const char *mpl_ext(plot_format_t fmt) {
  switch (fmt) {
  case PLOT_FORMAT_PNG:
    return "png";
  case PLOT_FORMAT_PDF:
    return "pdf";
  case PLOT_FORMAT_SVG:
    return "svg";
  case PLOT_FORMAT_EPS:
    return "eps";
  case PLOT_FORMAT_JPEG:
    return "jpg";
  default:
    return "png";
  }
}

// Builds Python single-quoted string literal
static void py_repr(char *out, size_t out_size, const char *s) {
  size_t o = 0;
  if (o + 1 < out_size)
    out[o++] = '\'';
  for (const char *p = s; *p && o + 2 < out_size; p++) {
    if (*p == '\\' || *p == '\'')
      out[o++] = '\\';
    out[o++] = *p;
  }
  if (o + 1 < out_size)
    out[o++] = '\'';
  out[o] = '\0';
}

static void apply_common_opts(matplotlib_t *mp, const plot_opts_t *opts) {
  int w = (opts && opts->width > 0) ? opts->width : 800;
  int h = (opts && opts->height > 0) ? opts->height : 600;
  matplotlib_cmd(mp, "fig = plt.figure(figsize=(%g, %g), dpi=100)", w / 100.0,
                 h / 100.0);

  char buf[256];
  if (opts && opts->title) {
    py_repr(buf, sizeof buf, opts->title);
    matplotlib_cmd(mp, "plt.title(%s)", buf);
  }
  if (opts && opts->xlabel) {
    py_repr(buf, sizeof buf, opts->xlabel);
    matplotlib_cmd(mp, "plt.xlabel(%s)", buf);
  }
  if (opts && opts->ylabel) {
    py_repr(buf, sizeof buf, opts->ylabel);
    matplotlib_cmd(mp, "plt.ylabel(%s)", buf);
  }
  matplotlib_cmd(mp, "plt.grid(True)");
}

int plot_line(const char *filename, plot_format_t format, const double *x,
              const double *y, int n, const plot_opts_t *opts) {
  if (format == PLOT_FORMAT_WINDOW) {
    // TODO: Implement GUI backend
    fprintf(stderr,
            "plot: matplotlib backend does not implement PLOT_FORMAT_WINDOW\n");
    return -1;
  }

  if (!python3_available()) {
    fprintf(stderr, "plot: python3 not found - skipping '%s'\n", filename);
    return -1;
  }

  if (!matplotlib_module_available()) {
    fprintf(stderr,
            "plot: matplotlib module not importable - skipping '%s'\n"
            "      Install: pip install matplotlib\n",
            filename);
    return -1;
  }

  matplotlib_t *mp = matplotlib_open();
  if (!mp)
    return -1;

  apply_common_opts(mp, opts);

  matplotlib_cmd(mp, "x = []");
  matplotlib_cmd(mp, "y = []");
  for (int i = 0; i < n; i++)
    matplotlib_cmd(mp, "x.append(%.10e); y.append(%.10e)", x[i], y[i]);

  double lw = (opts && opts->line_width > 0) ? opts->line_width : 2.0;
  if (opts && opts->color) {
    char cbuf[64];
    py_repr(cbuf, sizeof cbuf, opts->color);
    matplotlib_cmd(mp, "plt.plot(x, y, linewidth=%g, color=%s)", lw, cbuf);
  } else {
    matplotlib_cmd(mp, "plt.plot(x, y, linewidth=%g)", lw);
  }

  char path[512], pbuf[600];
  snprintf(path, sizeof path, "%s/%s.%s", QMC_OUTPUT_DIR, filename,
           mpl_ext(format));
  py_repr(pbuf, sizeof pbuf, path);
  matplotlib_cmd(mp, "plt.savefig(%s, dpi=100)", pbuf);
  matplotlib_cmd(mp, "plt.close(fig)");

  matplotlib_close(mp);
  return 0;
}

int plot_lines(const char *filename, plot_format_t format, const double *x,
               const double **ys, int n_series, int n_pts, const char **labels,
               const plot_opts_t *opts) {
  if (format == PLOT_FORMAT_WINDOW) {
    // TODO: Implement GUI backend
    fprintf(stderr,
            "plot: matplotlib backend does not implement PLOT_FORMAT_WINDOW\n");
    return -1;
  }
  if (!python3_available() || !matplotlib_module_available()) {
    fprintf(stderr, "plot: python3/matplotlib not available - skipping '%s'\n",
            filename);
    return -1;
  }

  matplotlib_t *mp = matplotlib_open();
  if (!mp)
    return -1;

  apply_common_opts(mp, opts);

  double lw = (opts && opts->line_width > 0) ? opts->line_width : 2.0;
  matplotlib_cmd(mp, "xs = []");
  for (int i = 0; i < n_pts; i++)
    matplotlib_cmd(mp, "xs.append(%.10e)", x[i]);

  int any_labels = 0;
  for (int s = 0; s < n_series; s++) {
    matplotlib_cmd(mp, "ys%d = []", s);
    for (int i = 0; i < n_pts; i++)
      matplotlib_cmd(mp, "ys%d.append(%.10e)", s, ys[s][i]);

    if (labels && labels[s]) {
      char lbuf[128];
      py_repr(lbuf, sizeof lbuf, labels[s]);
      matplotlib_cmd(mp, "plt.plot(xs, ys%d, linewidth=%g, label=%s)", s, lw,
                     lbuf);
      any_labels = 1;
    } else {
      matplotlib_cmd(mp, "plt.plot(xs, ys%d, linewidth=%g)", s, lw);
    }
  }

  if (any_labels)
    matplotlib_cmd(mp, "plt.legend()");

  char path[512], pbuf[600];
  snprintf(path, sizeof path, "%s/%s.%s", QMC_OUTPUT_DIR, filename,
           mpl_ext(format));
  py_repr(pbuf, sizeof pbuf, path);
  matplotlib_cmd(mp, "plt.savefig(%s, dpi=100)", pbuf);
  matplotlib_cmd(mp, "plt.close('all')");

  matplotlib_close(mp);
  return 0;
}
