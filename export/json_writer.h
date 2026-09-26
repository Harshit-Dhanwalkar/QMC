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
 *   - Doubles are written with shortest %g precision (1..17 digits)
 *     whose decimal text round-trips bit-exactly. NaN/Inf are written as
 *     `null`, since neither is a valid JSON number (RFC 8259)
 *   - Strings are escaped per RFC 8259 §7 (quote, backslash, control chars)
 *     A NULL string is written as `""`
 *   - Object keys must be non-NULL
 */

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

/* NaN/Inf handling for json_write_double()/json_write_double_array():
 *   JSON_NONFINITE_NULL  - write `null` (default, matches historical
 *                          behavior; RFC 8259 has no NaN/Inf literal)
 *   JSON_NONFINITE_ERROR - fail write with JSON_ERR_NONFINITE instead,
 *                          for callers where silently mapping a NaN/Inf
 *                          result to `null` could hide a numerical failure
 */
typedef enum {
  JSON_NONFINITE_NULL,
  JSON_NONFINITE_ERROR
} json_nonfinite_policy_t;

typedef struct {
  int pretty;      /* 1 = multi-line with indent; 0 = compact (default 1) */
  unsigned indent; /* spaces per level; ignored if !pretty (default 2) */
  json_nonfinite_policy_t nonfinite; /* default JSON_NONFINITE_NULL */
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
 * If any prior write returned an error, json_close still closes file and
 * returns that last error status. Returns JSON_OK only if file flushed and
 * closed cleanly and no prior error occurred */
json_status_t json_close(json_writer_t *writer);

/* Last error recorded on writer (JSON_OK if none) */
json_status_t json_status(const json_writer_t *writer);

/* Human-readable status string */
const char *json_strerror(json_status_t status);

/* Nested objects
 * `json_begin_object` writes `"key": {`; writes inside go into that nested
 * object until matching `json_end_object`
 *
 * Pass key=NULL for an object opened as an *array element* (i.e. when
 * enclosing container is an array opened with json_begin_array, not an
 * object) - array elements have no key. A non-NULL key inside an array, or
 * a NULL key inside an object, is JSON_ERR_INVALID_ARGUMENT
 */
json_status_t json_begin_object(json_writer_t *writer, const char *key);
json_status_t json_end_object(json_writer_t *writer);

/* Nested arrays: `json_begin_array` writes `"key": [`; writes inside go into
 * that array (as unkeyed elements - pass key=NULL to every write call made
 * while inside it) until matching `json_end_array`. Like json_begin_object,
 * key=NULL is required/allowed exactly when enclosing container is itself an
 * array. Unlike json_write_double_array / json_write_int_array /
 * json_write_string_array (flat arrays of one primitive type written inline in
 * a single call), this streams arbitrary mixed or nested elements one at a time
 * - e.g. an array of objects:
 *
 *   json_begin_array(w, "results");
 *   json_begin_object(w, NULL);
 *   json_write_double(w, "energy", -2.9037);
 *   json_end_object(w);
 *   json_begin_object(w, NULL);
 *   json_write_double(w, "energy", -2.9012);
 *   json_end_object(w);
 *   json_end_array(w);
 */
json_status_t json_begin_array(json_writer_t *writer, const char *key);
json_status_t json_end_array(json_writer_t *writer);

/* Scalar writes at current depth. Pass key=NULL when writing an
 * unkeyed element directly inside an array opened with json_begin_array */
json_status_t json_write_string(json_writer_t *writer, const char *key,
                                const char *value);
json_status_t json_write_int(json_writer_t *writer, const char *key,
                             long value);
json_status_t json_write_double(json_writer_t *writer, const char *key,
                                double value);
json_status_t json_write_bool(json_writer_t *writer, const char *key,
                              int value);
json_status_t json_write_null(json_writer_t *writer, const char *key);

/* Write `{"value": value, "unit": unit, "description": description}` as a
 * single nested object at `key` - a lightweight way to attach units/
 * descriptions to individual scientific quantities without wrapping every
 * scalar this way by default. `unit`/`description` may be NULL to omit
 * that member entirely (not written as "unit": null - simply absent) */
json_status_t json_write_quantity(json_writer_t *writer, const char *key,
                                  double value, const char *unit,
                                  const char *description);

/* Flat arrays written inline as a single value */
json_status_t json_write_double_array(json_writer_t *writer, const char *key,
                                      const double *data, size_t len);
json_status_t json_write_int_array(json_writer_t *writer, const char *key,
                                   const long *data, size_t len);
json_status_t json_write_string_array(json_writer_t *writer, const char *key,
                                      const char *const *data, size_t len);

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

/* Write a pre-built json_field_t at current depth of an open writer */
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
