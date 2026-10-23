# `json_writer`

> Generated from `export/json_writer.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Streaming JSON writer for results and metadata: objects, arrays, quantities with units, schema header.

**Category:** Export, plotting & LaTeX

```c
#include <qmc/export/json_writer.h>
```

**Start with:** [`json_open`](#json_open), [`json_begin_object`](#json_begin_object), [`json_write_field`](#json_write_field), [`json_close`](#json_close)

- Test (every call with expected results): `tests/test_json_writer.c`
- Runtime help: `qmc_help(stdout, "json_writer");` or `qmc-help json_writer`

## Overview

```text
Minimal streaming JSON metadata writer

Conventions:
  - Doubles are written with shortest %g precision (1..17 digits)
    whose decimal text round-trips bit-exactly. NaN/Inf are written as
    `null`, since neither is a valid JSON number (RFC 8259)
  - Strings are escaped per RFC 8259 §7 (quote, backslash, control chars)
    A NULL string is written as `""`
  - Object keys must be non-NULL
```

## Types

### `json_status_t`

```c
typedef enum {
  JSON_OK = 0,
  JSON_ERR_INVALID_ARGUMENT = -1,
  JSON_ERR_OPEN = -2,
  JSON_ERR_IO = -3,
  JSON_ERR_DEPTH = -4,
  JSON_ERR_INVALID_STATE = -5,
  JSON_ERR_NONFINITE = -6,
  JSON_ERR_MEMORY = -7
} json_status_t;
```

### `json_nonfinite_policy_t`

```text
NaN/Inf handling for json_write_double()/json_write_double_array():
  JSON_NONFINITE_NULL  - write `null` (default, matches historical
                         behavior; RFC 8259 has no NaN/Inf literal)
  JSON_NONFINITE_ERROR - fail write with JSON_ERR_NONFINITE instead,
                         for callers where silently mapping a NaN/Inf
                         result to `null` could hide a numerical failure
```

```c
typedef enum {
  JSON_NONFINITE_NULL,
  JSON_NONFINITE_ERROR
} json_nonfinite_policy_t;
```

### `json_options_t`

```c
typedef struct {
  int pretty;      /* 1 = multi-line with indent; 0 = compact (default 1) */
  unsigned indent; /* spaces per level; ignored if !pretty (default 2) */
  json_nonfinite_policy_t nonfinite; /* default JSON_NONFINITE_NULL */
} json_options_t;
```

### `json_writer_t`

```c
typedef struct json_writer json_writer_t;
```

### `json_field_type_t`

```c
typedef enum {
  JSON_FIELD_STRING,      /* value.s: NUL-terminated C string */
  JSON_FIELD_INT,         /* value.i: a signed integer */
  JSON_FIELD_DOUBLE,      /* value.d: a double (NaN/Inf are written as null) */
  JSON_FIELD_DOUBLE_ARRAY /* value.arr: a flat array of doubles */
} json_field_type_t;
```

### `json_field_t`

```c
typedef struct {
  const char *key;
  json_field_type_t type;
  union {
    const char *s;
    long i;
    double d;
    struct {
      const double *data;
      size_t len;
    } arr;
  } value;
} json_field_t;
```

## Functions (27)

[`json_options_default`](#json_options_default), [`json_open`](#json_open), [`json_open_with_options`](#json_open_with_options), [`json_close`](#json_close), [`json_status`](#json_status), [`json_strerror`](#json_strerror), [`json_begin_object`](#json_begin_object), [`json_end_object`](#json_end_object), [`json_begin_array`](#json_begin_array), [`json_end_array`](#json_end_array), [`json_write_string`](#json_write_string), [`json_write_int`](#json_write_int), [`json_write_double`](#json_write_double), [`json_write_bool`](#json_write_bool), [`json_write_null`](#json_write_null), [`json_write_quantity`](#json_write_quantity), [`json_write_double_array`](#json_write_double_array), [`json_write_int_array`](#json_write_int_array), [`json_write_string_array`](#json_write_string_array), [`json_field_string`](#json_field_string), [`json_field_int`](#json_field_int), [`json_field_double`](#json_field_double), [`json_field_double_array`](#json_field_double_array), [`json_write_field`](#json_write_field), [`json_write_metadata`](#json_write_metadata), [`json_write_schema_header`](#json_write_schema_header), [`json_write_columns`](#json_write_columns)

### `json_options_default`

```c
json_options_t json_options_default(void);
```

_No description in the header yet._

### `json_open`

```c
json_writer_t *json_open(const char *path);
```

```text
Open `path` for writing. `path` is used verbatim - no QMC_OUTPUT_DIR prefix

Returns NULL on failure (e.g. fopen failed)
```

### `json_open_with_options`

```c
json_writer_t *json_open_with_options(const char *path,
                                      const json_options_t *options);
```

_No description in the header yet._

### `json_close`

```c
json_status_t json_close(json_writer_t *writer);
```

```text
Flush, close, free. Safe to call exactly once

If any prior write returned an error, json_close still closes file and
returns that last error status. Returns JSON_OK only if file flushed and
closed cleanly and no prior error occurred
```

### `json_status`

```c
json_status_t json_status(const json_writer_t *writer);
```

```text
Last error recorded on writer (JSON_OK if none)
```

### `json_strerror`

```c
const char *json_strerror(json_status_t status);
```

```text
Human-readable status string
```

### `json_begin_object`

```c
json_status_t json_begin_object(json_writer_t *writer, const char *key);
```

```text
Nested objects
`json_begin_object` writes `"key": {`; writes inside go into that nested
object until matching `json_end_object`

Pass key=NULL for an object opened as an *array element* (i.e. when
enclosing container is an array opened with json_begin_array, not an
object) - array elements have no key. A non-NULL key inside an array, or
a NULL key inside an object, is JSON_ERR_INVALID_ARGUMENT
```

### `json_end_object`

```c
json_status_t json_end_object(json_writer_t *writer);
```

_No description in the header yet._

### `json_begin_array`

```c
json_status_t json_begin_array(json_writer_t *writer, const char *key);
```

```text
Nested arrays: `json_begin_array` writes `"key": [`; writes inside go into
that array (as unkeyed elements - pass key=NULL to every write call made
while inside it) until matching `json_end_array`. Like json_begin_object,
key=NULL is required/allowed exactly when enclosing container is itself an
array. Unlike json_write_double_array / json_write_int_array /
json_write_string_array (flat arrays of one primitive type written inline in
a single call), this streams arbitrary mixed or nested elements one at a time
- e.g. an array of objects:

  json_begin_array(w, "results");
  json_begin_object(w, NULL);
  json_write_double(w, "energy", -2.9037);
  json_end_object(w);
  json_begin_object(w, NULL);
  json_write_double(w, "energy", -2.9012);
  json_end_object(w);
  json_end_array(w);
```

### `json_end_array`

```c
json_status_t json_end_array(json_writer_t *writer);
```

_No description in the header yet._

### `json_write_string`

```c
json_status_t json_write_string(json_writer_t *writer, const char *key,
                                const char *value);
```

```text
Scalar writes at current depth. Pass key=NULL when writing an
unkeyed element directly inside an array opened with json_begin_array
```

### `json_write_int`

```c
json_status_t json_write_int(json_writer_t *writer, const char *key,
                             long value);
```

_No description in the header yet._

### `json_write_double`

```c
json_status_t json_write_double(json_writer_t *writer, const char *key,
                                double value);
```

_No description in the header yet._

### `json_write_bool`

```c
json_status_t json_write_bool(json_writer_t *writer, const char *key,
                              int value);
```

_No description in the header yet._

### `json_write_null`

```c
json_status_t json_write_null(json_writer_t *writer, const char *key);
```

_No description in the header yet._

### `json_write_quantity`

```c
json_status_t json_write_quantity(json_writer_t *writer, const char *key,
                                  double value, const char *unit,
                                  const char *description);
```

```text
Write `{"value": value, "unit": unit, "description": description}` as a
single nested object at `key` - a lightweight way to attach units/
descriptions to individual scientific quantities without wrapping every
scalar this way by default. `unit`/`description` may be NULL to omit
that member entirely (not written as "unit": null - simply absent)
```

### `json_write_double_array`

```c
json_status_t json_write_double_array(json_writer_t *writer, const char *key,
                                      const double *data, size_t len);
```

```text
Flat arrays written inline as a single value
```

### `json_write_int_array`

```c
json_status_t json_write_int_array(json_writer_t *writer, const char *key,
                                   const long *data, size_t len);
```

_No description in the header yet._

### `json_write_string_array`

```c
json_status_t json_write_string_array(json_writer_t *writer, const char *key,
                                      const char *const *data, size_t len);
```

_No description in the header yet._

### `json_field_string`

```c
json_field_t json_field_string(const char *key, const char *value);
```

_No description in the header yet._

### `json_field_int`

```c
json_field_t json_field_int(const char *key, long value);
```

_No description in the header yet._

### `json_field_double`

```c
json_field_t json_field_double(const char *key, double value);
```

_No description in the header yet._

### `json_field_double_array`

```c
json_field_t json_field_double_array(const char *key, const double *data,
                                     size_t len);
```

_No description in the header yet._

### `json_write_field`

```c
json_status_t json_write_field(json_writer_t *writer,
                               const json_field_t *field);
```

```text
Write a pre-built json_field_t at current depth of an open writer
```

### `json_write_metadata`

```c
int json_write_metadata(const char *filename, const json_field_t *fields,
                        size_t n_fields);
```

```text
Write `fields` as a single flat JSON object to QMC_OUTPUT_DIR/filename

Returns 0 on success, -1 on failure
```

### `json_write_schema_header`

```c
json_status_t json_write_schema_header(json_writer_t *writer,
                                       int include_timestamp);
```

```text
Writes standard QMC provenance block into current object:

  "schema_version": 1,
  "library": { "name": "QMC", "version": "..." },
  "created": "YYYY-MM-DDTHH:MM:00Z"        (only if include_timestamp)

Call it first, right after json_open() (or inside any begin_object), so
every run.json is self-identifying. Pass include_timestamp = 0 for
byte-reproducible output (e.g. golden-file tests)
```

### `json_write_columns`

```c
json_status_t json_write_columns(json_writer_t *writer, const char *key,
                                 const export_column_t *columns, size_t n);
```

```text
Writes an array of column descriptors under `key`:
  "key": [ {"name": "energy", "unit": "Ha", "description": "..."}, ... ]
`unit` / `description` members are omitted when NULL. Pairs with
csv_write_header_columns(): CSV carries plain names, this array carries
units and descriptions
```
