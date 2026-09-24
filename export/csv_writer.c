#include "csv_writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

#define CSV_MAX_PATH 4096
#define CSV_MAX_COLS_PER_ROW 4096

struct csv_writer {
  FILE *file;
  csv_options_t opts;
  int cols_written_this_row;
  int wrote_anything_this_row;
  csv_status_t last_error;
};

csv_options_t csv_options_default(void) {
  csv_options_t o = {
      .delimiter = ',',
      .precision = 15,
      .scientific = 1,
      .write_header = 1,
      .nonfinite = CSV_NONFINITE_LITERAL,
  };

  return o;
}

csv_writer_t *csv_open(const char *path, const csv_options_t *options) {
  if (!path || path[0] == '\0') {
    return NULL;
  }

  csv_writer_t *w = calloc(1, sizeof *w);
  if (!w) {
    return NULL;
  }

  w->opts = options ? *options : csv_options_default();

  w->file = fopen(path, "w");
  if (!w->file) {
    free(w);
    return NULL;
  }

  return w;
}

csv_writer_t *csv_open_append(const char *path, const csv_options_t *options) {
  if (!path || path[0] == '\0') {
    return NULL;
  }

  csv_writer_t *w = calloc(1, sizeof *w);
  if (!w) {
    return NULL;
  }

  w->opts = options ? *options : csv_options_default();
  w->opts.write_header = 0;

  w->file = fopen(path, "a");
  if (!w->file) {
    free(w);

    return NULL;
  }

  return w;
}

csv_status_t csv_close(csv_writer_t *writer) {
  if (!writer) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  csv_status_t rc = writer->last_error;

  if (writer->file) {
    if (fflush(writer->file) != 0 && rc == CSV_OK) {
      rc = CSV_ERR_IO;
    }
    if (fclose(writer->file) != 0 && rc == CSV_OK) {
      rc = CSV_ERR_IO;
    }
  }

  free(writer);

  return rc;
}

csv_status_t csv_status(csv_writer_t *writer) {
  return writer ? writer->last_error : CSV_ERR_INVALID_ARGUMENT;
}

/* Quote a field if it contains the delimiter, a quote, CR, or LF.
 * RFC 4180: fields containing these must be double-quoted, and inner
 * quotes doubled. */
static csv_status_t csv_write_escaped(csv_writer_t *writer, const char *s) {
  if (!writer || !writer->file) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  if (!s) {
    s = "";
  }

  int needs_quotes = 0;
  for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
    if (*p == (unsigned char)writer->opts.delimiter || *p == '"' ||
        *p == '\n' || *p == '\r') {
      needs_quotes = 1;

      break;
    }
  }

  int rc;
  if (needs_quotes) {
    rc = fputc('"', writer->file);
    if (rc == EOF) {
      writer->last_error = CSV_ERR_IO;

      return CSV_ERR_IO;
    }

    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
      if (*p == '"') {
        if (fputc('"', writer->file) == EOF ||
            fputc('"', writer->file) == EOF) {
          writer->last_error = CSV_ERR_IO;

          return CSV_ERR_IO;
        }
      } else {
        if (fputc((int)*p, writer->file) == EOF) {
          writer->last_error = CSV_ERR_IO;

          return CSV_ERR_IO;
        }
      }
    }
    if (fputc('"', writer->file) == EOF) {
      writer->last_error = CSV_ERR_IO;

      return CSV_ERR_IO;
    }
  } else {
    if (fputs(s, writer->file) == EOF) {
      writer->last_error = CSV_ERR_IO;

      return CSV_ERR_IO;
    }
  }

  return CSV_OK;
}

static csv_status_t csv_write_separator(csv_writer_t *writer) {
  if (!writer || !writer->file) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  if (writer->wrote_anything_this_row) {
    if (fputc((int)writer->opts.delimiter, writer->file) == EOF) {
      writer->last_error = CSV_ERR_IO;

      return CSV_ERR_IO;
    }
  } else {
    writer->wrote_anything_this_row = 1;
  }

  return CSV_OK;
}

csv_status_t csv_write_string(csv_writer_t *writer, const char *value) {
  if (!writer) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  csv_status_t rc = csv_write_separator(writer);
  if (rc != CSV_OK) {
    return rc;
  }

  return csv_write_escaped(writer, value);
}

csv_status_t csv_write_double(csv_writer_t *writer, double value) {
  if (!writer || !writer->file) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  csv_status_t rc = csv_write_separator(writer);
  if (rc != CSV_OK) {
    return rc;
  }

  if (!isfinite(value)) {
    switch (writer->opts.nonfinite) {
    case CSV_NONFINITE_EMPTY:
      return CSV_OK; /* empty field */
    case CSV_NONFINITE_ERROR:
      writer->last_error = CSV_ERR_NONFINITE;
      return CSV_ERR_NONFINITE;
    case CSV_NONFINITE_LITERAL:
    default:
      if (isnan(value)) {
        if (fputs("nan", writer->file) == EOF) {
          writer->last_error = CSV_ERR_IO;
          return CSV_ERR_IO;
        }
      } else if (value > 0) {
        if (fputs("inf", writer->file) == EOF) {
          writer->last_error = CSV_ERR_IO;
          return CSV_ERR_IO;
        }
      } else {
        if (fputs("-inf", writer->file) == EOF) {
          writer->last_error = CSV_ERR_IO;
          return CSV_ERR_IO;
        }
      }
      return CSV_OK;
    }
  }

  int prec = writer->opts.precision > 0 ? writer->opts.precision : 15;
  char fmt[16];
  snprintf(fmt, sizeof fmt, writer->opts.scientific ? "%%.%de" : "%%.%dg",
           prec);

  if (fprintf(writer->file, fmt, value) < 0) {
    writer->last_error = CSV_ERR_IO;

    return CSV_ERR_IO;
  }

  return CSV_OK;
}

csv_status_t csv_write_int(csv_writer_t *writer, long value) {
  if (!writer || !writer->file) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  csv_status_t rc = csv_write_separator(writer);
  if (rc != CSV_OK) {
    return rc;
  }

  if (fprintf(writer->file, "%ld", value) < 0) {
    writer->last_error = CSV_ERR_IO;

    return CSV_ERR_IO;
  }

  return CSV_OK;
}

csv_status_t csv_end_row(csv_writer_t *writer) {
  if (!writer || !writer->file) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  if (fputc('\n', writer->file) == EOF) {
    writer->last_error = CSV_ERR_IO;

    return CSV_ERR_IO;
  }

  writer->wrote_anything_this_row = 0;

  return CSV_OK;
}

csv_status_t csv_write_row(csv_writer_t *writer, const double *values,
                           size_t n) {
  if (!writer || (n > 0 && !values)) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  for (size_t i = 0; i < n; i++) {
    csv_status_t rc = csv_write_double(writer, values[i]);
    if (rc != CSV_OK) {
      return rc;
    }
  }

  return csv_end_row(writer);
}

const char *csv_strerror(csv_status_t status) {
  switch (status) {
  case CSV_OK:
    return "ok";
  case CSV_ERR_INVALID_ARGUMENT:
    return "invalid argument";
  case CSV_ERR_OPEN:
    return "could not open file";
  case CSV_ERR_IO:
    return "I/O error";
  case CSV_ERR_PATH_TOO_LONG:
    return "path too long";
  case CSV_ERR_INVALID_STATE:
    return "invalid state";
  case CSV_ERR_NONFINITE:
    return "non-finite value";
  default:
    return "unknown";
  }
}

static int build_output_path(const char *filename, char *out, size_t out_size) {
  (void)mkdir(QMC_OUTPUT_DIR, 0755);
  int n = snprintf(out, out_size, "%s/%s", QMC_OUTPUT_DIR, filename);
  if (n < 0 || (size_t)n >= out_size) {
    return -1;
  }

  return 0;
}

int csv_write_1d(const char *filename, const double *x, const double *y, int n,
                 const char *xlabel, const char *ylabel) {
  if (!filename || (n > 0 && (!x || !y))) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  char path[CSV_MAX_PATH];
  if (build_output_path(filename, path, sizeof path) != 0) {
    return CSV_ERR_PATH_TOO_LONG;
  }

  csv_options_t opts = csv_options_default();
  csv_writer_t *w = csv_open(path, &opts);
  if (!w) {
    return CSV_ERR_OPEN;
  }

  csv_write_string(w, xlabel ? xlabel : "x");
  csv_write_string(w, ylabel ? ylabel : "y");
  csv_end_row(w);

  for (int i = 0; i < n; i++) {
    csv_write_double(w, x[i]);
    csv_write_double(w, y[i]);
    csv_end_row(w);
  }

  return (int)csv_close(w);
}

int csv_write_matrix(const char *filename, const double *data, int rows,
                     int cols, const char **col_headers) {
  if (!filename || rows < 0 || cols < 0 || (rows > 0 && cols > 0 && !data)) {
    return CSV_ERR_INVALID_ARGUMENT;
  }

  char path[CSV_MAX_PATH];
  if (build_output_path(filename, path, sizeof path) != 0) {
    return CSV_ERR_PATH_TOO_LONG;
  }

  csv_options_t opts = csv_options_default();
  csv_writer_t *w = csv_open(path, &opts);
  if (!w) {
    return CSV_ERR_OPEN;
  }

  if (col_headers) {
    for (int j = 0; j < cols; j++) {
      csv_write_string(w, col_headers[j]);
    }
  } else {
    for (int j = 0; j < cols; j++) {
      char tmp[32];
      snprintf(tmp, sizeof tmp, "col%d", j);
      csv_write_string(w, tmp);
    }
  }

  csv_end_row(w);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      csv_write_double(w, data[i * cols + j]);
    }

    csv_end_row(w);
  }

  return (int)csv_close(w);
}
