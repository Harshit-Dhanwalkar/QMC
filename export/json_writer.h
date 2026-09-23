#ifndef QMC_JSON_WRITER_H
#define QMC_JSON_WRITER_H

#include <stddef.h>

// Minimal JSON metadata writer
//
// WARN: Not a general JSON writer: no nested objects/arrays, no
// arbitrary value types. If a use case needs more than this, it
// has outgrown "run metadata"

typedef enum {
  JSON_FIELD_STRING,      // value.s: NUL-terminated C string
  JSON_FIELD_INT,         // value.i: a signed integer
  JSON_FIELD_DOUBLE,      // value.d: a double (NaN/Inf are
                          // written as null, since neither is
                          // valid JSON)
  JSON_FIELD_DOUBLE_ARRAY // value.arr: a flat array of doubles
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

// Convenience constructors, so callers don't have to hand-fill union
json_field_t json_field_string(const char *key, const char *value);
json_field_t json_field_int(const char *key, long value);
json_field_t json_field_double(const char *key, double value);
json_field_t json_field_double_array(const char *key, const double *data,
                                     size_t len);

// Write `fields` as a single flat JSON object to QMC_OUTPUT_DIR/filename,
// e.g. {"grid_n": 160, "seed": 42, "energies": [1.0, 2.0]}
// Returns 0 on success, -1 if the output file could not be opened
int json_write_metadata(const char *filename, const json_field_t *fields,
                        size_t n_fields);

#endif // QMC_JSON_WRITER_H
