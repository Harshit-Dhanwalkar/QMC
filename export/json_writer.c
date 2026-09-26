#include "json_writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

#define JSON_MAX_DEPTH 32
#define JSON_MAX_PATH 512

/* Write `s` as a JSON string literal (with surrounding quotes), escaping
 * characters per RFC 8259 §7
 * NULL writes as ""
 */
static void write_json_string(FILE *f, const char *s) {
  fputc('"', f);
  if (s) {
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
      switch (*p) {
      case '"':
        fputs("\\\"", f);
        break;
      case '\\':
        fputs("\\\\", f);
        break;
      case '\n':
        fputs("\\n", f);
        break;
      case '\r':
        fputs("\\r", f);
        break;
      case '\t':
        fputs("\\t", f);
        break;
      default:
        if (*p < 0x20) {
          fprintf(f, "\\u%04x", *p);
        } else {
          fputc((int)*p, f);
        }
      }
    }
  }

  fputc('"', f);
}

/* Write a double as a JSON number, or `null` for NaN/Inf
 *
 * Uses shortest decimal precision (1..17 significant digits) whose value reads
 * back bit-for-bit identical via strtod (Reference: Steele & White 1990), so
 * e.g. -2.9037 is written as "-2.9037" rather than "-2.903700..."
 */
static void write_json_double(FILE *f, double d) {
  if (isnan(d) || isinf(d)) {
    fputs("null", f);

    return;
  }

  char buf[64];
  for (int precision = 1; precision <= 17; precision++) {
    snprintf(buf, sizeof buf, "%.*g", precision, d);

    if (strtod(buf, NULL) == d) {
      break;
    }
  }

  fputs(buf, f);
}

static void write_json_bool(FILE *f, int value) {
  fputs(value ? "true" : "false", f);
}

/* Streaming writer */
typedef enum { JSON_CTX_OBJECT, JSON_CTX_ARRAY } json_ctx_t;

struct json_writer {
  FILE *f;
  json_options_t opts;
  int depth;                      /* 0 = inside top-level object */
  int need_comma[JSON_MAX_DEPTH]; /* per-level: has a value been written? */
  json_ctx_t ctx[JSON_MAX_DEPTH]; /* per-level: object or array container */
  json_status_t last_error;
  int closed;
};

json_options_t json_options_default(void) {
  json_options_t o = {
      .pretty = 1, .indent = 2, .nonfinite = JSON_NONFINITE_NULL};

  return o;
}

static json_status_t set_error(json_writer_t *w, json_status_t err) {
  if (w->last_error == JSON_OK) {
    w->last_error = err;
  }

  return w->last_error;
}

json_writer_t *json_open(const char *path) {
  return json_open_with_options(path, NULL);
}

json_writer_t *json_open_with_options(const char *path,
                                      const json_options_t *options) {
  if (!path || path[0] == '\0') {
    return NULL;
  }

  FILE *f = fopen(path, "w");
  if (!f) {
    return NULL;
  }

  json_writer_t *w = calloc(1, sizeof *w);
  if (!w) {
    fclose(f);

    return NULL;
  }

  w->f = f;
  w->opts = options ? *options : json_options_default();
  if (w->opts.indent == 0) {
    w->opts.indent = 2;
  }
  w->depth = 0;
  w->need_comma[0] = 0;
  w->ctx[0] = JSON_CTX_OBJECT; // top-level container is always an object
  w->last_error = JSON_OK;
  w->closed = 0;

  if (fputc('{', f) == EOF) {
    set_error(w, JSON_ERR_IO);
  }

  return w;
}

json_status_t json_status(const json_writer_t *writer) {
  return writer ? writer->last_error : JSON_ERR_INVALID_ARGUMENT;
}

const char *json_strerror(json_status_t status) {
  switch (status) {
  case JSON_OK:
    return "ok";
  case JSON_ERR_INVALID_ARGUMENT:
    return "invalid argument";
  case JSON_ERR_OPEN:
    return "could not open file";
  case JSON_ERR_IO:
    return "I/O error";
  case JSON_ERR_DEPTH:
    return "maximum nesting depth exceeded";
  case JSON_ERR_INVALID_STATE:
    return "invalid writer state (mismatched begin/end?)";
  case JSON_ERR_NONFINITE:
    return "NaN/Inf value rejected under JSON_NONFINITE_ERROR policy";
  case JSON_ERR_MEMORY:
    return "out of memory";
  default:
    return "unknown error";
  }
}

json_status_t json_close(json_writer_t *writer) {
  if (!writer) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (writer->closed) {
    json_status_t rc = writer->last_error;

    free(writer);

    return rc;
  }

  json_status_t rc = writer->last_error;

  if (writer->depth != 0 && rc == JSON_OK) {
    rc = JSON_ERR_INVALID_STATE;
  }

  if (writer->opts.pretty) {
    if (fputc('\n', writer->f) == EOF && rc == JSON_OK) {
      rc = JSON_ERR_IO;
    }
  }
  if (fputc('}', writer->f) == EOF && rc == JSON_OK) {
    rc = JSON_ERR_IO;
  }
  if (fputc('\n', writer->f) == EOF && rc == JSON_OK) {
    rc = JSON_ERR_IO;
  }
  if (fflush(writer->f) != 0 && rc == JSON_OK) {
    rc = JSON_ERR_IO;
  }
  if (fclose(writer->f) != 0 && rc == JSON_OK) {
    rc = JSON_ERR_IO;
  }

  writer->closed = 1;

  free(writer);

  return rc;
}

/* Emit a comma if needed at current depth, then a newline + indent
 * (when pretty). Also flips need_comma bit for this depth
 */
static json_status_t element_prefix(json_writer_t *w) {
  if (w->need_comma[w->depth]) {
    if (fputc(',', w->f) == EOF) {
      return set_error(w, JSON_ERR_IO);
    }
  }
  if (w->opts.pretty) {
    if (fputc('\n', w->f) == EOF) {
      return set_error(w, JSON_ERR_IO);
    }

    unsigned n = (unsigned)(w->depth + 1) * w->opts.indent;
    for (unsigned i = 0; i < n; i++) {
      if (fputc(' ', w->f) == EOF) {
        return set_error(w, JSON_ERR_IO);
      }
    }
  }

  w->need_comma[w->depth] = 1;

  return JSON_OK;
}

static json_status_t write_key(json_writer_t *w, const char *key) {
  if (!key) {
    return JSON_OK; // unkeyed array element - nothing to write
  }

  write_json_string(w->f, key);

  if (fputs(": ", w->f) == EOF) {
    return set_error(w, JSON_ERR_IO);
  }

  return JSON_OK;
}

/* Every element-writing function (scalars, begin_object, begin_array) must
 * agree with enclosing container on whether key is given: object members
 * require one, array elements must not have one
 */
static int key_mismatches_context(json_writer_t *w, const char *key) {
  int in_array = w->ctx[w->depth] == JSON_CTX_ARRAY;

  return (in_array && key != NULL) || (!in_array && key == NULL);
}

json_status_t json_begin_object(json_writer_t *writer, const char *key) {
  if (!writer || writer->closed) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (writer->depth + 1 >= JSON_MAX_DEPTH) {
    return set_error(writer, JSON_ERR_DEPTH);
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  if (fputc('{', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  writer->depth++;
  writer->need_comma[writer->depth] = 0;
  writer->ctx[writer->depth] = JSON_CTX_OBJECT;

  return JSON_OK;
}

json_status_t json_end_object(json_writer_t *writer) {
  if (!writer || writer->closed) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (writer->depth <= 0 || writer->ctx[writer->depth] != JSON_CTX_OBJECT) {
    return set_error(writer, JSON_ERR_INVALID_STATE);
  }

  if (writer->opts.pretty) {
    if (fputc('\n', writer->f) == EOF) {
      return set_error(writer, JSON_ERR_IO);
    }

    unsigned n = (unsigned)writer->depth * writer->opts.indent;
    for (unsigned i = 0; i < n; i++) {
      if (fputc(' ', writer->f) == EOF) {
        return set_error(writer, JSON_ERR_IO);
      }
    }
  }

  if (fputc('}', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  writer->depth--;

  return JSON_OK;
}

json_status_t json_begin_array(json_writer_t *writer, const char *key) {
  if (!writer || writer->closed) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (writer->depth + 1 >= JSON_MAX_DEPTH) {
    return set_error(writer, JSON_ERR_DEPTH);
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  if (fputc('[', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  writer->depth++;
  writer->need_comma[writer->depth] = 0;
  writer->ctx[writer->depth] = JSON_CTX_ARRAY;

  return JSON_OK;
}

json_status_t json_end_array(json_writer_t *writer) {
  if (!writer || writer->closed) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (writer->depth <= 0 || writer->ctx[writer->depth] != JSON_CTX_ARRAY) {
    return set_error(writer, JSON_ERR_INVALID_STATE);
  }

  if (writer->opts.pretty) {
    if (fputc('\n', writer->f) == EOF) {
      return set_error(writer, JSON_ERR_IO);
    }

    unsigned n = (unsigned)writer->depth * writer->opts.indent;
    for (unsigned i = 0; i < n; i++) {
      if (fputc(' ', writer->f) == EOF) {
        return set_error(writer, JSON_ERR_IO);
      }
    }
  }

  if (fputc(']', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  writer->depth--;

  return JSON_OK;
}

json_status_t json_write_string(json_writer_t *writer, const char *key,
                                const char *value) {
  if (!writer || writer->closed || key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  write_json_string(writer->f, value);

  return JSON_OK;
}

json_status_t json_write_int(json_writer_t *writer, const char *key,
                             long value) {
  if (!writer || !key || writer->closed) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  if (fprintf(writer->f, "%ld", value) < 0) {
    return set_error(writer, JSON_ERR_IO);
  }

  return JSON_OK;
}

json_status_t json_write_double(json_writer_t *writer, const char *key,
                                double value) {
  if (!writer || writer->closed || key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if ((isnan(value) || isinf(value)) &&
      writer->opts.nonfinite == JSON_NONFINITE_ERROR) {
    return set_error(writer, JSON_ERR_NONFINITE);
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  write_json_double(writer->f, value);

  return JSON_OK;
}

json_status_t json_write_bool(json_writer_t *writer, const char *key,
                              int value) {
  if (!writer || writer->closed || key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  write_json_bool(writer->f, value);

  return JSON_OK;
}

json_status_t json_write_null(json_writer_t *writer, const char *key) {
  if (!writer || (!data && len > 0) || writer->closed ||
      key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if (writer->opts.nonfinite == JSON_NONFINITE_ERROR) {
    for (size_t i = 0; i < len; i++) {
      if (isnan(data[i]) || isinf(data[i])) {
        return set_error(writer, JSON_ERR_NONFINITE);
      }
    }
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  if (fputs("null", writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  return JSON_OK;
}

json_status_t json_write_double_array(json_writer_t *writer, const char *key,
                                      const double *data, size_t len) {
  if (!writer || (!data && len > 0) || writer->closed ||
      key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  if (fputc('[', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  for (size_t i = 0; i < len; i++) {
    if (i > 0 && fputs(", ", writer->f) == EOF) {
      return set_error(writer, JSON_ERR_IO);
    }

    write_json_double(writer->f, data[i]);
  }

  if (fputc(']', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  return JSON_OK;
}

json_status_t json_write_int_array(json_writer_t *writer, const char *key,
                                   const long *data, size_t len) {
  if (!writer || (!data && len > 0) || writer->closed ||
      key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  if (fputc('[', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  for (size_t i = 0; i < len; i++) {
    if (i > 0 && fputs(", ", writer->f) == EOF) {
      return set_error(writer, JSON_ERR_IO);
    }

    if (fprintf(writer->f, "%ld", data[i]) < 0) {
      return set_error(writer, JSON_ERR_IO);
    }
  }

  if (fputc(']', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  return JSON_OK;
}

json_status_t json_write_string_array(json_writer_t *writer, const char *key,
                                      const char *const *data, size_t len) {
  if (!writer || (!data && len > 0) || writer->closed ||
      key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  json_status_t rc = element_prefix(writer);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = write_key(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  if (fputc('[', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  for (size_t i = 0; i < len; i++) {
    if (i > 0 && fputs(", ", writer->f) == EOF) {
      return set_error(writer, JSON_ERR_IO);
    }

    write_json_string(writer->f, data[i]);
  }

  if (fputc(']', writer->f) == EOF) {
    return set_error(writer, JSON_ERR_IO);
  }

  return JSON_OK;
}

json_status_t json_write_quantity(json_writer_t *writer, const char *key,
                                  double value, const char *unit,
                                  const char *description) {
  if (!writer || writer->closed || key_mismatches_context(writer, key)) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  if ((isnan(value) || isinf(value)) &&
      writer->opts.nonfinite == JSON_NONFINITE_ERROR) {
    return set_error(writer, JSON_ERR_NONFINITE);
  }

  json_status_t rc = json_begin_object(writer, key);
  if (rc != JSON_OK) {
    return rc;
  }

  rc = json_write_double(writer, "value", value);
  if (rc != JSON_OK) {
    return rc;
  }

  if (unit) {
    rc = json_write_string(writer, "unit", unit);
    if (rc != JSON_OK) {
      return rc;
    }
  }

  if (description) {
    rc = json_write_string(writer, "description", description);
    if (rc != JSON_OK) {
      return rc;
    }
  }

  return json_end_object(writer);
}

json_field_t json_field_string(const char *key, const char *value) {
  json_field_t f;
  f.key = key;
  f.type = JSON_FIELD_STRING;
  f.value.s = value;

  return f;
}

json_field_t json_field_int(const char *key, long value) {
  json_field_t f;
  f.key = key;
  f.type = JSON_FIELD_INT;
  f.value.i = value;

  return f;
}

json_field_t json_field_double(const char *key, double value) {
  json_field_t f;
  f.key = key;
  f.type = JSON_FIELD_DOUBLE;
  f.value.d = value;

  return f;
}

json_field_t json_field_double_array(const char *key, const double *data,
                                     size_t len) {
  json_field_t f;
  f.key = key;
  f.type = JSON_FIELD_DOUBLE_ARRAY;
  f.value.arr.data = data;
  f.value.arr.len = len;

  return f;
}

json_status_t json_write_field(json_writer_t *writer,
                               const json_field_t *field) {
  if (!writer || !field) {
    return JSON_ERR_INVALID_ARGUMENT;
  }

  switch (field->type) {
  case JSON_FIELD_STRING:
    return json_write_string(writer, field->key, field->value.s);
  case JSON_FIELD_INT:
    return json_write_int(writer, field->key, field->value.i);
  case JSON_FIELD_DOUBLE:
    return json_write_double(writer, field->key, field->value.d);
  case JSON_FIELD_DOUBLE_ARRAY:
    return json_write_double_array(writer, field->key, field->value.arr.data,
                                   field->value.arr.len);
  default:
    return JSON_ERR_INVALID_ARGUMENT;
  }
}

int json_write_metadata(const char *filename, const json_field_t *fields,
                        size_t n_fields) {
  if (!filename || filename[0] == '\0' || (n_fields > 0 && !fields)) {
    return -1;
  }

  (void)mkdir(QMC_OUTPUT_DIR, 0755);

  char path[JSON_MAX_PATH];
  int n = snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  if (n < 0 || (size_t)n >= sizeof path) {
    return -1;
  }

  json_writer_t *w = json_open(path);
  if (!w) {
    return -1;
  }

  for (size_t i = 0; i < n_fields; i++) {
    if (json_write_field(w, &fields[i]) != JSON_OK) {
      (void)json_close(w);

      return -1;
    }
  }

  json_status_t rc = json_close(w);

  return (rc == JSON_OK) ? 0 : -1;
}
