#include "json_writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

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

// Write `s` as a JSON string literal (with surrounding quotes), escaping
// the characters JSON requires (RFC 8259 section 7): backslash, double
// quote, and the C0 control characters. A NULL string writes as "".
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

// Write a double as a JSON number, or `null` for NaN/Inf since neither is
// valid JSON (RFC 8259 has no representation for them) - silently emitting
// a fabricated finite value in their place would be worse than admitting
// the value isn't representable.
//
// Uses the shortest decimal precision (1 to 17 significant digits) whose
// value reads back bit-for-bit identical via strtod, rather than always
// emitting the maximum 17 digits - so e.g. -2.9037 is written as
// "-2.9037" instead of "-2.9037000000000002". 17 significant digits is
// the well-known bound (Steele & White 1990) below which every double is
// guaranteed representable and round-trippable, so trying 1..17 in order
// and stopping at the first exact match always terminates and is always
// exact - this is purely a formatting choice, never a precision loss.
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

static void write_json_field(FILE *f, const json_field_t *field) {
  write_json_string(f, field->key);
  fputs(": ", f);

  switch (field->type) {
  case JSON_FIELD_STRING:
    write_json_string(f, field->value.s);
    break;

  case JSON_FIELD_INT:
    fprintf(f, "%ld", field->value.i);
    break;

  case JSON_FIELD_DOUBLE:
    write_json_double(f, field->value.d);
    break;

  case JSON_FIELD_DOUBLE_ARRAY:
    fputc('[', f);
    for (size_t i = 0; i < field->value.arr.len; i++) {
      if (i > 0) {
        fputs(", ", f);
      }
      write_json_double(f, field->value.arr.data[i]);
    }
    fputc(']', f);
    break;
  }
}

int json_write_metadata(const char *filename, const json_field_t *fields,
                        size_t n_fields) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);

  FILE *f = fopen(path, "w");
  if (!f) {
    return -1;
  }

  fputs("{\n", f);
  for (size_t i = 0; i < n_fields; i++) {
    fputs("  ", f);
    write_json_field(f, &fields[i]);
    fputs((i + 1 < n_fields) ? ",\n" : "\n", f);
  }
  fputs("}\n", f);

  fclose(f);

  return 0;
}
