#ifndef QMC_JSON_WRITER_H
#define QMC_JSON_WRITER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Minimal streaming JSON metadata writer
 *
 * Conventions:
 *   - Doubles are written with the shortest %g precision (1..17 digits)
 *     whose decimal text round-trips bit-exactly. NaN/Inf are written as
 *     `null`, since neither is a valid JSON number (RFC 8259)
 *   - Strings are escaped per RFC 8259 §7 (quote, backslash, control chars)
 *     A NULL string is written as `""`.
 *   - Object keys must be non-NULL
 */

typedef enum {
  JSON_OK = 0,
  JSON_ERR_INVALID_ARGUMENT = -1,
  JSON_ERR_OPEN = -2,
  JSON_ERR_IO = -3,
  JSON_ERR_DEPTH = -4,
  JSON_ERR_INVALID_STATE = -5
} json_status_t;

typedef struct {
  int pretty;      /* 1 = multi-line with indent; 0 = compact (default 1) */
  unsigned indent; /* spaces per level; ignored if !pretty (default 2) */
} json_options_t;

json_options_t json_options_default(void);

typedef struct json_writer json_writer_t;

/* Open `path` for writing. `path` is used verbatim - no QMC_OUTPUT_DIR prefix
 *
 * Returns NULL on failure (e.g. fopen failed) */
json_writer_t *json_open(const char *path);
json_writer_t *json_open_with_options(const char *path,
                                      const json_options_t *options);

/* Flush, close, free. Safe to call exactly once
 *
 * If any prior write returned an error, json_close still closes the file  and
 * returns that last error status. Returns JSON_OK only if the file flushed and
 * closed cleanly and no prior error occurred */
json_status_t json_close(json_writer_t *writer);

/* Last error recorded on the writer (JSON_OK if none) */
json_status_t json_status(const json_writer_t *writer);

/* Human-readable status string */
const char *json_strerror(json_status_t status);

/* Nested objects
 * `json_begin_object` writes `"key": {`; writes inside go into that nested
 * object until the matching `json_end_object. */
json_status_t json_begin_object(json_writer_t *writer, const char *key);
json_status_t json_end_object(json_writer_t *writer);

/* Scalar writes at the current depth */
json_status_t json_write_string(json_writer_t *writer, const char *key,
                                const char *value);
json_status_t json_write_int(json_writer_t *writer, const char *key,
                             long value);
json_status_t json_write_double(json_writer_t *writer, const char *key,
                                double value);
json_status_t json_write_bool(json_writer_t *writer, const char *key,
                              int value);
json_status_t json_write_null(json_writer_t *writer, const char *key);

/* Flat arrays written inline as a single value */
json_status_t json_write_double_array(json_writer_t *writer, const char *key,
                                      const double *data, size_t len);
json_status_t json_write_int_array(json_writer_t *writer, const char *key,
                                   const long *data, size_t len);

typedef enum {
  JSON_FIELD_STRING,      /* value.s: NUL-terminated C string */
  JSON_FIELD_INT,         /* value.i: a signed integer */
  JSON_FIELD_DOUBLE,      /* value.d: a double (NaN/Inf are written as null) */
  JSON_FIELD_DOUBLE_ARRAY /* value.arr: a flat array of doubles */
} json_field_type_t;

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

json_field_t json_field_string(const char *key, const char *value);
json_field_t json_field_int(const char *key, long value);
json_field_t json_field_double(const char *key, double value);
json_field_t json_field_double_array(const char *key, const double *data,
                                     size_t len);

/* Write a pre-built json_field_t at the current depth of an open writer */
json_status_t json_write_field(json_writer_t *writer,
                               const json_field_t *field);

/*
 * Write `fields` as a single flat JSON object to QMC_OUTPUT_DIR/filename
 *
 * Returns 0 on success, -1 on failure
 */
int json_write_metadata(const char *filename, const json_field_t *fields,
                        size_t n_fields);

#ifdef __cplusplus
}
#endif

#endif /* QMC_JSON_WRITER_H */
