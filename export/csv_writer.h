#ifndef QMC_CSV_WRITER_H
#define QMC_CSV_WRITER_H

#include <stddef.h>

#include "export_schema.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  CSV_OK = 0,
  CSV_ERR_INVALID_ARGUMENT = -1,
  CSV_ERR_OPEN = -2,
  CSV_ERR_IO = -3,
  CSV_ERR_PATH_TOO_LONG = -4,
  CSV_ERR_INVALID_STATE = -5,
  CSV_ERR_NONFINITE = -6
} csv_status_t;

typedef enum {
  CSV_NONFINITE_LITERAL, /* nan / inf / -inf */
  CSV_NONFINITE_EMPTY,   /* empty field */
  CSV_NONFINITE_ERROR    /* fail the write */
} csv_nonfinite_policy_t;

typedef struct {
  char delimiter;                   /* default ',' */
  int precision;                    /* default 15 */
  int scientific;                   /* 1 -> %.*e, 0 -> %.*g */
  int write_header;                 /* default 1 */
  csv_nonfinite_policy_t nonfinite; /* default CSV_NONFINITE_LITERAL */
} csv_options_t;

/* Returns default options (delimiter=',', precision=15, scientific=1,
 * write_header=1, nonfinite=LITERAL) */
csv_options_t csv_options_default(void);

/* Opaque handle for the streaming writer */
typedef struct csv_writer csv_writer_t;

/* Open a file for writing. Path is used verbatim */
csv_writer_t *csv_open(const char *path, const csv_options_t *options);

/* Appends to an existing file */
csv_writer_t *csv_open_append(const char *path, const csv_options_t *options);

/* Flush, close, free
 *
 * Returns CSV_OK only if flush+close+checks all succeeded
 */
csv_status_t csv_close(csv_writer_t *writer);

/* Streaming row API */
csv_status_t csv_write_string(csv_writer_t *writer, const char *value);
csv_status_t csv_write_double(csv_writer_t *writer, double value);
csv_status_t csv_status(csv_writer_t *writer); /* last error, or CSV_OK */
csv_status_t csv_write_int(csv_writer_t *writer, long value);
csv_status_t csv_end_row(csv_writer_t *writer);
csv_status_t csv_write_row(csv_writer_t *writer, const double *values,
                           size_t n);

/* Column descriptor (name / unit / description), shared with the JSON and HDF5
 * exporters. Only `name` ever reaches the CSV itself - units and descriptions
 * belong in the JSON sidecar (see json_write_columns) so the CSV header stays
 * plain and parseable */
typedef export_column_t csv_column_t;

/* Writes one header row from `columns[i].name` and ends the row. Names go
 * through the normal RFC 4180 escaping. Use with streaming csv_write_row() */
csv_status_t csv_write_header_columns(csv_writer_t *writer,
                                      const csv_column_t *columns, size_t n);

/* One-shot table export (row-major `data`, rows x ncols) with a descriptor
 * header. `path` is used verbatim (no implicit QMC_OUTPUT_DIR) and `options`
 * may be NULL for defaults. rows == 0 writes a header-only file. Returns CSV_OK
 * or a negative csv_status_t */
csv_status_t csv_write_table(const char *path, const csv_column_t *columns,
                             size_t ncols, const double *data, size_t rows,
                             const csv_options_t *options);

int csv_write_1d(const char *filename, const double *x, const double *y,
                 size_t n, const char *xlabel, const char *ylabel);
int csv_write_matrix(const char *filename, const double *data, size_t rows,
                     size_t cols, const char **col_headers);

const char *csv_strerror(csv_status_t status);

#ifdef __cplusplus
}
#endif

#endif /* QMC_CSV_WRITER_H */
