/*
 * Test: JSON metadata writer (export/json_writer.c)
 *
 * 1. Exact-output checks for each field type (string, int, double, double
 *    array), including string-escaping edge cases and NaN/Inf -> null
 * 2. Invalid-input handling (NULL filename must not crash; 0 fields must still
 *    produce a valid, empty-bodied JSON object)
 */

#include "../export/json_writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

static int failures = 0;

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

// Read a whole file into a NUL-terminated heap buffer
static char *slurp(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) {
    return NULL;
  }

  fseek(f, 0, SEEK_END);
  long len = ftell(f);
  fseek(f, 0, SEEK_SET);

  char *buf = malloc((size_t)len + 1);
  if (!buf) {
    fclose(f);

    return NULL;
  }

  size_t got = fread(buf, 1, (size_t)len, f);
  buf[got] = '\0';
  fclose(f);

  return buf;
}

static void test_field_types_and_escaping(void) {
  printf("  === Test field_types_and_escaping ===\n");

  const double energies[3] = {-2.9037, 0.0, 1.5};

  json_field_t fields[6] = {
      json_field_string("method", "DMC"),
      json_field_string("note", "quote:\" backslash:\\ tab:\tend"),
      json_field_int("n_walkers", 200),
      json_field_double("target_energy", -2.9037),
      json_field_double_array("energies", energies, 3),
      json_field_double("bad_value", NAN),
  };

  int rc = json_write_metadata("test_json_writer_types.json", fields, 6);
  check_true(rc == 0, "json_write_metadata returns 0 on success");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
           "test_json_writer_types.json");
  char *content = slurp(path);
  check_true(content != NULL, "output file was created and is readable");
  if (!content) {
    return;
  }

  check_true(strstr(content, "\"method\": \"DMC\"") != NULL,
             "string field round-trips exactly");
  check_true(
      strstr(content, "\"note\": \"quote:\\\" backslash:\\\\ tab:\\tend\"") !=
          NULL,
      "quotes/backslash/tab are escaped correctly");
  check_true(strstr(content, "\"n_walkers\": 200") != NULL,
             "int field round-trips exactly");
  check_true(strstr(content, "\"target_energy\": -2.9037") != NULL,
             "double field round-trips exactly");
  check_true(strstr(content, "\"energies\": [-2.9037, 0, 1.5]") != NULL,
             "double array field round-trips exactly");
  check_true(strstr(content, "\"bad_value\": null") != NULL,
             "NaN is written as JSON null, not a fabricated number");

  // Braces present and the object is exactly one level deep
  check_true(content[0] == '{', "output starts with an opening brace");
  check_true(strchr(content, '}') != NULL, "output has a closing brace");
  check_true(strchr(content, '{') == content,
             "no nested '{' - flat object only, as documented");

  free(content);
}

static void test_empty_and_invalid_input(void) {
  printf("  === Test empty_and_invalid_input ===\n");

  int rc = json_write_metadata("test_json_writer_empty.json", NULL, 0);
  check_true(rc == 0, "0 fields still succeeds");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
           "test_json_writer_empty.json");
  char *content = slurp(path);
  check_true(content != NULL, "empty-metadata file was created");
  if (content) {
    check_true(strchr(content, '{') != NULL && strchr(content, '}') != NULL,
               "empty metadata is still a valid (empty) JSON object");

    free(content);
  }

  // A NULL string value for a string field must not crash - written as ""
  json_field_t f = json_field_string("label", NULL);
  rc = json_write_metadata("test_json_writer_null_string.json", &f, 1);
  check_true(rc == 0, "NULL string value does not crash, writes as \"\"");

  rc = json_write_metadata("/nonexistent_dir_xyz/cant_write.json", &f, 1);
  check_true(rc == -1, "unwritable path returns -1 rather than crashing");
}

static void test_nested_objects(void) {
  printf("  === Test nested_objects ===\n");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
           "test_json_writer_nested.json");

  json_writer_t *w = json_open(path);
  check_true(w != NULL, "json_open succeeds");
  if (!w) {
    return;
  }

  json_write_int(w, "schema_version", 1);

  check_true(json_begin_object(w, "run") == JSON_OK,
             "json_begin_object('run') succeeds");
  json_write_string(w, "method", "DMC");
  json_write_int(w, "seed", 12345);
  check_true(json_end_object(w) == JSON_OK, "json_end_object closes 'run'");

  check_true(json_begin_object(w, "parameters") == JSON_OK,
             "json_begin_object('parameters') succeeds");
  json_write_double(w, "time_step", 0.01);
  check_true(json_end_object(w) == JSON_OK,
             "json_end_object closes 'parameters'");

  // Mismatched end: closing an object that isn't open at this depth. This
  // sets a sticky error on the writer, which json_close() must propagate.
  check_true(json_end_object(w) == JSON_ERR_INVALID_STATE,
             "extra json_end_object at top level is rejected");

  check_true(json_close(w) == JSON_ERR_INVALID_STATE,
             "json_close propagates the sticky error from the rejected "
             "extra end_object");

  char *content = slurp(path);
  check_true(content != NULL, "nested output file readable");
  if (content) {
    check_true(strstr(content, "\"run\": {") != NULL,
               "'run' nested object present");
    check_true(strstr(content, "\"method\": \"DMC\"") != NULL,
               "nested string field present");
    check_true(strstr(content, "\"parameters\": {") != NULL,
               "'parameters' nested object present");
    free(content);
  }
}

static void test_arrays_of_objects(void) {
  printf("  === Test arrays_of_objects ===\n");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
           "test_json_writer_array.json");

  json_writer_t *w = json_open(path);
  check_true(w != NULL, "json_open succeeds");
  if (!w) {
    return;
  }

  check_true(json_begin_array(w, "results") == JSON_OK,
             "json_begin_array('results') succeeds");

  check_true(json_begin_object(w, NULL) == JSON_OK,
             "unkeyed json_begin_object inside an array succeeds");
  check_true(json_write_double(w, "energy", -2.9037) == JSON_OK,
             "keyed field inside an array-element object succeeds");
  check_true(json_end_object(w) == JSON_OK,
             "first array-element object closes");

  check_true(json_begin_object(w, NULL) == JSON_OK,
             "second unkeyed json_begin_object succeeds");
  json_write_double(w, "energy", -2.9012);
  check_true(json_end_object(w) == JSON_OK,
             "second array-element object closes");

  // A keyed value directly inside an array is a usage error
  check_true(json_write_int(w, "bad_key", 1) == JSON_ERR_INVALID_ARGUMENT,
             "keyed value directly inside an array is rejected");

  // An unkeyed scalar directly inside an array IS legitimate - it's simply
  // a bare element, same as what json_write_int_array would produce
  check_true(json_write_int(w, NULL, 7) == JSON_OK,
             "unkeyed json_write_int inside an array succeeds");

  check_true(json_end_array(w) == JSON_OK, "json_end_array closes 'results'");
  check_true(json_close(w) == JSON_OK, "json_close succeeds");

  char *content = slurp(path);
  check_true(content != NULL, "array-of-objects output file readable");
  if (content) {
    check_true(strstr(content, "\"results\": [") != NULL,
               "'results' array present");
    check_true(strstr(content, "-2.9037") != NULL &&
                   strstr(content, "-2.9012") != NULL,
               "both array-element objects' values are present");
    check_true(strstr(content, "7") != NULL,
               "trailing unkeyed scalar element is present in the array");
    free(content);
  }
}

static void test_string_array(void) {
  printf("  === Test string_array ===\n");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
           "test_json_writer_string_array.json");

  json_writer_t *w = json_open(path);
  const char *tags[3] = {"alpha", "quote:\"in\\here", "beta"};
  check_true(json_write_string_array(w, "tags", tags, 3) == JSON_OK,
             "json_write_string_array succeeds");
  check_true(json_close(w) == JSON_OK, "json_close succeeds");

  char *content = slurp(path);
  check_true(content != NULL, "string array output file readable");
  if (content) {
    check_true(strstr(content, "\"tags\": [\"alpha\", "
                               "\"quote:\\\"in\\\\here\", \"beta\"]") != NULL,
               "string array elements are individually escaped");
    free(content);
  }
}

static void test_quantity(void) {
  printf("  === Test quantity (value/unit/description) ===\n");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
           "test_json_writer_quantity.json");

  json_writer_t *w = json_open(path);
  check_true(json_write_quantity(w, "time_step", 0.01, "a.u.",
                                 "Imaginary-time step") == JSON_OK,
             "json_write_quantity with unit+description succeeds");
  check_true(json_write_quantity(w, "bare_value", 42.0, NULL, NULL) == JSON_OK,
             "json_write_quantity with NULL unit/description omits them");
  check_true(json_close(w) == JSON_OK, "json_close succeeds");

  char *content = slurp(path);
  check_true(content != NULL, "quantity output file readable");
  if (content) {
    check_true(strstr(content, "\"time_step\": {") != NULL,
               "'time_step' is written as a nested object");
    check_true(strstr(content, "\"value\": 0.01") != NULL,
               "quantity 'value' member present");
    check_true(strstr(content, "\"unit\": \"a.u.\"") != NULL,
               "quantity 'unit' member present");
    check_true(strstr(content, "\"description\": \"Imaginary-time step\"") !=
                   NULL,
               "quantity 'description' member present");
    check_true(strstr(content, "\"bare_value\": {\n    \"value\": 42") !=
                       NULL ||
                   strstr(content, "\"bare_value\": {") != NULL,
               "'bare_value' quantity written even with no unit/description");
    free(content);
  }
}

static void test_nonfinite_policy(void) {
  printf("  === Test nonfinite_policy ===\n");

  // Default (JSON_NONFINITE_NULL): unchanged historical behavior
  {
    char path[512];
    snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
             "test_json_writer_nonfinite_null.json");
    json_writer_t *w = json_open(path);
    check_true(json_write_double(w, "x", INFINITY) == JSON_OK,
               "default policy: +Inf write itself succeeds");
    check_true(json_close(w) == JSON_OK, "default policy: json_close succeeds");

    char *content = slurp(path);
    check_true(content && strstr(content, "\"x\": null") != NULL,
               "default policy still writes +Inf as null");
    free(content);
  }

  // JSON_NONFINITE_ERROR: the write itself fails, and it's sticky through
  // json_close()
  {
    char path[512];
    snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
             "test_json_writer_nonfinite_error.json");
    json_options_t opts = json_options_default();
    opts.nonfinite = JSON_NONFINITE_ERROR;
    json_writer_t *w = json_open_with_options(path, &opts);

    json_status_t rc = json_write_double(w, "x", NAN);
    check_true(rc == JSON_ERR_NONFINITE,
               "ERROR policy: NaN write fails with JSON_ERR_NONFINITE");
    check_true(json_status(w) == JSON_ERR_NONFINITE,
               "ERROR policy: json_status reflects the sticky error");
    check_true(json_close(w) == JSON_ERR_NONFINITE,
               "ERROR policy: json_close propagates the sticky error");
  }

  // JSON_NONFINITE_ERROR also rejects a nonfinite element inside a double
  // array, atomically (nothing from that array call is written)
  {
    char path[512];
    snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR,
             "test_json_writer_nonfinite_error_array.json");
    json_options_t opts = json_options_default();
    opts.nonfinite = JSON_NONFINITE_ERROR;
    json_writer_t *w = json_open_with_options(path, &opts);

    const double vals[3] = {1.0, NAN, 2.0};
    json_status_t rc = json_write_double_array(w, "vals", vals, 3);
    check_true(rc == JSON_ERR_NONFINITE,
               "ERROR policy: array containing NaN fails atomically");
    json_close(w);
  }
}

int main(void) {
  test_field_types_and_escaping();
  test_empty_and_invalid_input();
  test_nested_objects();
  test_arrays_of_objects();
  test_string_array();
  test_quantity();
  test_nonfinite_policy();

  if (failures == 0) {
    printf("\nAll test_json_writer checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
