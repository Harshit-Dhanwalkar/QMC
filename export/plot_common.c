/*
 * Shared helpers for all plot backends: strerror, output-path resolution,
 * and PLOT_FORMAT_TEXT fallback writer
 */

#include "plot.h"
#include "plot_internal.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Extension per format
 *
 * NULL for WINDOW and unrecognized 
 */
static const char *format_extension(plot_format_t format) {
  switch (format) {
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
  case PLOT_FORMAT_TEXT:
    return "csv";
  case PLOT_FORMAT_WINDOW:
  default:
    return NULL;
  }
}

const char *plot_strerror(plot_status_t status) {
  switch (status) {
  case PLOT_OK:
    return "ok";
  case PLOT_ERR_INVALID_ARGUMENT:
    return "invalid argument";
  case PLOT_ERR_BACKEND_UNAVAILABLE:
    return "plot backend not available";
  case PLOT_ERR_OPEN:
    return "could not open output file";
  case PLOT_ERR_WRITE:
    return "write failed";
  case PLOT_ERR_FORMAT_UNSUPPORTED:
    return "format not supported";
  case PLOT_ERR_NOT_IMPLEMENTED:
    return "feature not implemented";
  case PLOT_ERR_PATH_TOO_LONG:
    return "output path too long";
  case PLOT_ERR_INTERNAL:
    return "internal error";
  default:
    return "unknown plot error";
  }
}

plot_status_t plot_output_path(const char *filename, plot_format_t format,
                               const plot_opts_t *opts, char *out,
                               size_t out_size) {
  if (!filename || !out || out_size == 0) {
    return PLOT_ERR_INVALID_ARGUMENT;
  }

  const char *ext = format_extension(format);
  if (!ext) {
    return PLOT_ERR_FORMAT_UNSUPPORTED;
  }

  const char *dir =
      (opts && opts->output_dir) ? opts->output_dir : QMC_OUTPUT_DIR;

  int n;
  if (dir && dir[0] != '\0') {
    n = snprintf(out, out_size, "%s/%s.%s", dir, filename, ext);
  } else {
    n = snprintf(out, out_size, "%s.%s", filename, ext);
  }

  if (n < 0 || (size_t)n >= out_size) {
    return PLOT_ERR_PATH_TOO_LONG;
  }

  return PLOT_OK;
}

// PLOT_FORMAT_TEXT writers

/*
 * Full-precision decimal formatting: try 1..17 significant digits and stop
 * at the first whose strtod round-trips bit-exactly
 */
static void write_double(FILE *f, double v) {
  char buf[64];
  for (int precision = 1; precision <= 17; precision++) {
    snprintf(buf, sizeof buf, "%.*g", precision, v);
    if (strtod(buf, NULL) == v) {
      break;
    }
  }

  fputs(buf, f);
}

plot_status_t plot_write_text_1d(const char *path, const double *x,
                                 const double *y, int n,
                                 const plot_opts_t *opts) {
  if (!path || (n > 0 && (!x || !y))) {
    return PLOT_ERR_INVALID_ARGUMENT;
  }

  FILE *f = fopen(path, "w");
  if (!f) {
    return PLOT_ERR_OPEN;
  }

  fputs(opts && opts->xlabel ? opts->xlabel : "x", f);
  fputc(',', f);
  fputs(opts && opts->ylabel ? opts->ylabel : "y", f);
  fputc('\n', f);

  for (int i = 0; i < n; i++) {
    write_double(f, x[i]);
    fputc(',', f);
    write_double(f, y[i]);
    fputc('\n', f);
  }

  if (ferror(f)) {
    fclose(f);
    return PLOT_ERR_WRITE;
  }
  if (fclose(f) != 0) {
    return PLOT_ERR_WRITE;
  }

  return PLOT_OK;
}

plot_status_t plot_write_text_nd(const char *path, const double *x,
                                 const double **ys, int n_series, int n_pts,
                                 const char **labels, const plot_opts_t *opts) {
  if (!path || n_series <= 0 || n_pts < 0) {
    return PLOT_ERR_INVALID_ARGUMENT;
  }
  if (n_pts > 0 && (!x || !ys)) {
    return PLOT_ERR_INVALID_ARGUMENT;
  }

  FILE *f = fopen(path, "w");
  if (!f) {
    return PLOT_ERR_OPEN;
  }

  fputs(opts && opts->xlabel ? opts->xlabel : "x", f);
  for (int s = 0; s < n_series; s++) {
    fputc(',', f);
    if (labels && labels[s]) {
      fputs(labels[s], f);
    } else {
      char buf[32];
      snprintf(buf, sizeof buf, "series%d", s);

      fputs(buf, f);
    }
  }

  fputc('\n', f);

  for (int i = 0; i < n_pts; i++) {
    write_double(f, x[i]);
    for (int s = 0; s < n_series; s++) {
      fputc(',', f);
      write_double(f, ys[s][i]);
    }

    fputc('\n', f);
  }

  if (ferror(f)) {
    fclose(f);
    return PLOT_ERR_WRITE;
  }

  if (fclose(f) != 0) {
    return PLOT_ERR_WRITE;
  }

  return PLOT_OK;
}
