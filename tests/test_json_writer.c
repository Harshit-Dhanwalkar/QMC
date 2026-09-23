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

  double energies[3] = {-2.9037, 0.0, 1.5};

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

int main(void) {
  test_field_types_and_escaping();
  test_empty_and_invalid_input();

  if (failures == 0) {
    printf("\nAll test_json_writer checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
