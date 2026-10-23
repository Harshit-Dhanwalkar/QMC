# `csv_writer`

> Generated from `export/csv_writer.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Streaming CSV writer: configurable delimiter and precision, RFC 4180 escaping, append mode, column descriptors.

**Category:** Export, plotting & LaTeX

```c
#include <qmc/export/csv_writer.h>
```

**Start with:** [`csv_options_default`](#csv_options_default), [`csv_open`](#csv_open), [`csv_write_row`](#csv_write_row), [`csv_close`](#csv_close)

- Test (every call with expected results): `tests/test_csv_writer.c`
- Runtime help: `qmc_help(stdout, "csv_writer");` or `qmc-help csv_writer`

## Types

### `csv_status_t`

```c
typedef enum {
  CSV_OK = 0,
  CSV_ERR_INVALID_ARGUMENT = -1,
  CSV_ERR_OPEN = -2,
  CSV_ERR_IO = -3,
  CSV_ERR_PATH_TOO_LONG = -4,
  CSV_ERR_INVALID_STATE = -5,
  CSV_ERR_NONFINITE = -6
} csv_status_t;
```

### `csv_nonfinite_policy_t`

```c
typedef enum {
  CSV_NONFINITE_LITERAL, /* nan / inf / -inf */
  CSV_NONFINITE_EMPTY,   /* empty field */
  CSV_NONFINITE_ERROR    /* fail the write */
} csv_nonfinite_policy_t;
```

### `csv_options_t`

```c
typedef struct {
  char delimiter;                   /* default ',' */
  int precision;                    /* default 15 */
  int scientific;                   /* 1 -> %.*e, 0 -> %.*g */
  int write_header;                 /* default 1 */
  csv_nonfinite_policy_t nonfinite; /* default CSV_NONFINITE_LITERAL */
} csv_options_t;
```

### `csv_writer_t`

```text
Opaque handle for the streaming writer
```

```c
typedef struct csv_writer csv_writer_t;
```

### `csv_column_t`

```text
Column descriptor (name / unit / description), shared with the JSON and HDF5
exporters. Only `name` ever reaches the CSV itself - units and descriptions
belong in the JSON sidecar (see json_write_columns) so the CSV header stays
plain and parseable
```

```c
typedef export_column_t csv_column_t;
```

## Functions (15)

[`csv_options_default`](#csv_options_default), [`csv_open`](#csv_open), [`csv_open_append`](#csv_open_append), [`csv_close`](#csv_close), [`csv_write_string`](#csv_write_string), [`csv_write_double`](#csv_write_double), [`csv_status`](#csv_status), [`csv_write_int`](#csv_write_int), [`csv_end_row`](#csv_end_row), [`csv_write_row`](#csv_write_row), [`csv_write_header_columns`](#csv_write_header_columns), [`csv_write_table`](#csv_write_table), [`csv_write_1d`](#csv_write_1d), [`csv_write_matrix`](#csv_write_matrix), [`csv_strerror`](#csv_strerror)

### `csv_options_default`

```c
csv_options_t csv_options_default(void);
```

```text
Returns default options (delimiter=',', precision=15, scientific=1,
write_header=1, nonfinite=LITERAL)
```

### `csv_open`

```c
csv_writer_t *csv_open(const char *path, const csv_options_t *options);
```

```text
Open a file for writing. Path is used verbatim
```

### `csv_open_append`

```c
csv_writer_t *csv_open_append(const char *path, const csv_options_t *options);
```

```text
Appends to an existing file
```

### `csv_close`

```c
csv_status_t csv_close(csv_writer_t *writer);
```

```text
Flush, close, free

Returns CSV_OK only if flush+close+checks all succeeded
```

### `csv_write_string`

```c
csv_status_t csv_write_string(csv_writer_t *writer, const char *value);
```

```text
Streaming row API
```

### `csv_write_double`

```c
csv_status_t csv_write_double(csv_writer_t *writer, double value);
```

_No description in header yet._

### `csv_status`

```c
csv_status_t csv_status(csv_writer_t *writer);
```

_No description in header yet._

### `csv_write_int`

```c
csv_status_t csv_write_int(csv_writer_t *writer, long value);
```

```text
last error, or CSV_OK
```

### `csv_end_row`

```c
csv_status_t csv_end_row(csv_writer_t *writer);
```

_No description in header yet._

### `csv_write_row`

```c
csv_status_t csv_write_row(csv_writer_t *writer, const double *values,
                           size_t n);
```

_No description in header yet._

### `csv_write_header_columns`

```c
csv_status_t csv_write_header_columns(csv_writer_t *writer,
                                      const csv_column_t *columns, size_t n);
```

```text
Writes one header row from `columns[i].name` and ends the row. Names go
through the normal RFC 4180 escaping. Use with streaming csv_write_row()
```

### `csv_write_table`

```c
csv_status_t csv_write_table(const char *path, const csv_column_t *columns,
                             size_t ncols, const double *data, size_t rows,
                             const csv_options_t *options);
```

```text
One-shot table export (row-major `data`, rows x ncols) with a descriptor
header. `path` is used verbatim (no implicit QMC_OUTPUT_DIR) and `options`
may be NULL for defaults. rows == 0 writes a header-only file. Returns CSV_OK
or a negative csv_status_t
```

### `csv_write_1d`

```c
int csv_write_1d(const char *filename, const double *x, const double *y,
                 size_t n, const char *xlabel, const char *ylabel);
```

_No description in header yet._

### `csv_write_matrix`

```c
int csv_write_matrix(const char *filename, const double *data, size_t rows,
                     size_t cols, const char **col_headers);
```

_No description in header yet._

### `csv_strerror`

```c
const char *csv_strerror(csv_status_t status);
```

_No description in header yet._
