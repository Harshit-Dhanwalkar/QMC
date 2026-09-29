/*
 * Test: CSV writer (export/csv_writer.c)
 *
 * Covers the streaming writer (csv_open/csv_write_string/csv_write_double/
 * csv_end_row), RFC 4180 escaping, configurable delimiter/precision/ nonfinite
 * policy, append mode, and the legacy csv_write_1d/csv_write_matrix wrappers
 * (now size_t-dimensioned)
 */

#include "../export/csv_writer.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

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

/* Read an entire small file into a caller-provided buffer, NUL-terminated
 *
 * Returns 1 on success, 0 on failure (including truncation)
 */
static int slurp_file(const char *path, char *buf, size_t buf_size) {
  FILE *f = fopen(path, "r");
  if (!f) {
    return 0;
  }

  size_t n = fread(buf, 1, buf_size - 1, f);
  int truncated = (fgetc(f) != EOF);

  fclose(f);

  if (truncated) {
    return 0;
  }

  buf[n] = '\0';

  return 1;
}

static void test_streaming_basic(void) {
  printf("  === Test streaming_basic ===\n");

  const char *path = QMC_OUTPUT_DIR "/test_csv_streaming.csv";
  csv_writer_t *w = csv_open(path, NULL);
  check_true(w != NULL, "csv_open succeeds with default options");
  if (!w) {
    return;
  }

  check_true(csv_write_string(w, "step") == CSV_OK, "write header 'step'");
  check_true(csv_write_string(w, "energy") == CSV_OK, "write header 'energy'");
  check_true(csv_end_row(w) == CSV_OK, "end header row");

  check_true(csv_write_int(w, 0) == CSV_OK, "write int step 0");
  check_true(csv_write_double(w, -2.5) == CSV_OK, "write double -2.5");
  check_true(csv_end_row(w) == CSV_OK, "end data row");

  const double row[2] = {1.0, -2.9};
  check_true(csv_write_row(w, row, 2) == CSV_OK, "csv_write_row succeeds");

  check_true(csv_close(w) == CSV_OK, "csv_close succeeds cleanly");

  char content[512];
  check_true(slurp_file(path, content, sizeof content), "file readable");
  check_true(strstr(content, "step,energy\n") != NULL,
             "header row written correctly");
  check_true(strstr(content, "0,") != NULL, "int field present");
}

static void test_escaping(void) {
  printf("  === Test escaping (RFC 4180) ===\n");

  const char *path = QMC_OUTPUT_DIR "/test_csv_escaping.csv";
  csv_writer_t *w = csv_open(path, NULL);
  check_true(w != NULL, "csv_open succeeds");
  if (!w) {
    return;
  }

  csv_write_string(w, "Energy, total");
  csv_end_row(w);
  csv_write_string(w, "He said \"hello\"");
  csv_end_row(w);
  csv_write_string(w, "plain");
  csv_end_row(w);
  csv_close(w);

  char content[512];
  check_true(slurp_file(path, content, sizeof content), "file readable");
  check_true(strstr(content, "\"Energy, total\"") != NULL,
             "comma-containing field is quoted");
  check_true(strstr(content, "\"He said \"\"hello\"\"\"") != NULL,
             "embedded quotes are doubled per RFC 4180");
  check_true(strstr(content, "\nplain\n") != NULL,
             "plain field is left unquoted");
}

static void test_delimiter_and_precision(void) {
  printf("  === Test delimiter_and_precision ===\n");

  const char *path = QMC_OUTPUT_DIR "/test_csv_tsv.tsv";
  csv_options_t opts = csv_options_default();
  opts.delimiter = '\t';
  opts.precision = 3;
  opts.scientific = 0;

  csv_writer_t *w = csv_open(path, &opts);
  check_true(w != NULL, "csv_open succeeds with custom options");
  if (!w) {
    return;
  }

  csv_write_string(w, "a");
  csv_write_string(w, "b");
  csv_end_row(w);
  csv_write_double(w, 3.14159265);
  csv_write_double(w, 2.71828);
  csv_end_row(w);
  csv_close(w);

  char content[256];
  check_true(slurp_file(path, content, sizeof content), "file readable");
  check_true(strstr(content, "a\tb\n") != NULL,
             "tab delimiter used instead of comma");
  check_true(strstr(content, "3.14") != NULL,
             "precision=3, scientific=0 formats as %.3g-style");
}

static void test_nonfinite_policy(void) {
  printf("  === Test nonfinite_policy ===\n");

  // LITERAL (default): nan/inf/-inf as text
  {
    const char *path = QMC_OUTPUT_DIR "/test_csv_nonfinite_literal.csv";
    csv_writer_t *w = csv_open(path, NULL);
    csv_write_double(w, NAN);
    csv_write_double(w, INFINITY);
    csv_write_double(w, -INFINITY);
    csv_end_row(w);
    csv_close(w);

    char content[256];
    check_true(slurp_file(path, content, sizeof content), "file readable");
    check_true(strstr(content, "nan") != NULL, "NaN written as 'nan'");
    check_true(strstr(content, ",inf,") != NULL, "+Inf written as 'inf'");
    check_true(strstr(content, "-inf") != NULL, "-Inf written as '-inf'");
  }

  // EMPTY: nonfinite values become empty fields
  {
    const char *path = QMC_OUTPUT_DIR "/test_csv_nonfinite_empty.csv";
    csv_options_t opts = csv_options_default();
    opts.nonfinite = CSV_NONFINITE_EMPTY;
    csv_writer_t *w = csv_open(path, &opts);
    csv_write_double(w, 1.0);
    csv_write_double(w, NAN);
    csv_write_double(w, 2.0);
    csv_end_row(w);
    csv_close(w);

    char content[256];
    check_true(slurp_file(path, content, sizeof content), "file readable");
    check_true(strstr(content, "1.000000000000000e+00,,"
                               "2.000000000000000e+00") != NULL,
               "NaN becomes an empty field under CSV_NONFINITE_EMPTY");
  }

  // ERROR: nonfinite values fail write
  {
    csv_options_t opts = csv_options_default();
    opts.nonfinite = CSV_NONFINITE_ERROR;
    csv_writer_t *w =
        csv_open(QMC_OUTPUT_DIR "/test_csv_nonfinite_error.csv", &opts);
    csv_status_t rc = csv_write_double(w, NAN);
    check_true(rc == CSV_ERR_NONFINITE,
               "NaN write fails with CSV_ERR_NONFINITE under ERROR policy");
    check_true(csv_status(w) == CSV_ERR_NONFINITE,
               "csv_status reflects the sticky error");
    check_true(csv_close(w) == CSV_ERR_NONFINITE,
               "csv_close propagates the sticky error");
  }
}

static void test_append_mode(void) {
  printf("  === Test append_mode ===\n");

  const char *path = QMC_OUTPUT_DIR "/test_csv_append.csv";

  csv_writer_t *w = csv_open(path, NULL);
  csv_write_string(w, "step");
  csv_write_string(w, "energy");
  csv_end_row(w);
  csv_write_int(w, 0);
  csv_write_double(w, -1.0);
  csv_end_row(w);
  csv_close(w);

  csv_writer_t *a = csv_open_append(path, NULL);
  check_true(a != NULL, "csv_open_append succeeds");
  if (a) {
    csv_write_int(a, 1);
    csv_write_double(a, -1.5);
    csv_end_row(a);
    csv_close(a);
  }

  char content[512];
  check_true(slurp_file(path, content, sizeof content), "file readable");

  int header_count = 0;
  for (const char *p = content; (p = strstr(p, "step,energy")); p++) {
    header_count++;
  }

  check_true(header_count == 1, "append mode did not re-write the header row");
  check_true(strstr(content, "1,-1.5") != NULL,
             "appended row is present after the original rows");
}

static void test_null_and_error_handling(void) {
  printf("  === Test null_and_error_handling ===\n");

  check_true(csv_open(NULL, NULL) == NULL, "csv_open(NULL path) returns NULL");
  check_true(csv_open("", NULL) == NULL, "csv_open(empty path) returns NULL");

  csv_status_t rc = csv_write_row(NULL, NULL, 0);
  check_true(rc == CSV_ERR_INVALID_ARGUMENT,
             "csv_write_row(NULL writer) returns CSV_ERR_INVALID_ARGUMENT");

  // n > 0 with a NULL values pointer must not crash
  csv_writer_t *w = csv_open(QMC_OUTPUT_DIR "/test_csv_null_row.csv", NULL);
  if (w) {
    rc = csv_write_row(w, NULL, 3);
    check_true(rc == CSV_ERR_INVALID_ARGUMENT,
               "csv_write_row(n>0, NULL values) returns invalid-argument, "
               "not a crash");
    csv_close(w);
  }
}

static void test_legacy_write_1d_and_matrix(void) {
  printf("  === Test legacy_write_1d_and_matrix (size_t API) ===\n");

  const double x[4] = {0.0, 1.0, 2.0, 3.0};
  const double y[4] = {0.0, 1.0, 4.0, 9.0};

  // size_t dimensions (item 9 of the CSV update note)
  int rc = csv_write_1d("test_csv_legacy_1d.csv", x, y, (size_t)4, "x", "y");
  check_true(rc == CSV_OK, "csv_write_1d succeeds with size_t n");

  char path1[512];
  snprintf(path1, sizeof path1, "%s/%s", QMC_OUTPUT_DIR,
           "test_csv_legacy_1d.csv");
  char content[512];
  check_true(slurp_file(path1, content, sizeof content),
             "csv_write_1d output file readable");
  check_true(strstr(content, "x,y\n") != NULL,
             "csv_write_1d header uses given labels");

  // Zero-length is a valid empty table
  rc = csv_write_1d("test_csv_legacy_1d_empty.csv", NULL, NULL, (size_t)0, NULL,
                    NULL);
  check_true(rc == CSV_OK, "csv_write_1d with n=0 and NULL x/y succeeds");

  const double matrix[6] = {1, 2, 3, 4, 5, 6};
  const char *headers[3] = {"a", "b", "c"};
  rc = csv_write_matrix("test_csv_legacy_matrix.csv", matrix, (size_t)2,
                        (size_t)3, headers);
  check_true(rc == CSV_OK, "csv_write_matrix succeeds with size_t rows/cols");

  char path2[512];
  snprintf(path2, sizeof path2, "%s/%s", QMC_OUTPUT_DIR,
           "test_csv_legacy_matrix.csv");
  check_true(slurp_file(path2, content, sizeof content),
             "csv_write_matrix output file readable");
  check_true(strstr(content, "a,b,c\n") != NULL,
             "csv_write_matrix header uses given column headers");
  check_true(strstr(content, "1.000000000000000e+00,"
                             "2.000000000000000e+00,"
                             "3.000000000000000e+00\n") != NULL,
             "csv_write_matrix first data row correct (default precision/"
             "scientific formatting)");
  check_true(strstr(content, "4.000000000000000e+00,"
                             "5.000000000000000e+00,"
                             "6.000000000000000e+00\n") != NULL,
             "csv_write_matrix second data row correct");

  // NULL data with rows>0 and cols>0 must be rejected, not crash
  rc = csv_write_matrix("test_csv_legacy_matrix_null.csv", NULL, (size_t)2,
                        (size_t)3, NULL);
  check_true(rc == CSV_ERR_INVALID_ARGUMENT,
             "csv_write_matrix(NULL data, rows>0, cols>0) rejected cleanly");
}

static void test_column_descriptors_and_table(void) {
  printf("  === Test column descriptors + csv_write_table ===\n");

  const csv_column_t cols[3] = {
      {"step", NULL, NULL},
      {"Energy, total", "Ha", "Total energy"},
      {"say \"hi\"", "Ha^2", "Variance"},
  };
  const double data[6] = {0.0, -2.5, 0.25, 1.0, -2.75, 0.125};

  const char *path = QMC_OUTPUT_DIR "/test_csv_table.csv";
  csv_options_t opts = csv_options_default();
  opts.precision = 3;
  opts.scientific = 0;
  csv_status_t rc = csv_write_table(path, cols, 3, data, 2, &opts);
  check_true(rc == CSV_OK, "csv_write_table succeeds");

  char buf[512];
  check_true(slurp_file(path, buf, sizeof buf), "table file readable");
  check_true(strcmp(buf, "step,\"Energy, total\",\"say \"\"hi\"\"\"\n"
                         "0,-2.5,0.25\n1,-2.75,0.125\n") == 0,
             "descriptor names escaped per RFC 4180, %g values written");
  check_true(strstr(buf, "Ha") == NULL,
             "units/descriptions never leak into the CSV itself");

  // rows == 0 -> header-only file, data may be NULL
  rc = csv_write_table(path, cols, 3, NULL, 0, NULL);
  check_true(rc == CSV_OK, "rows == 0 with NULL data is valid");
  check_true(slurp_file(path, buf, sizeof buf) &&
                 strcmp(buf, "step,\"Energy, total\",\"say \"\"hi\"\"\"\n") ==
                     0,
             "rows == 0 produces a header-only file");

  check_true(csv_write_table(NULL, cols, 3, data, 2, NULL) ==
                 CSV_ERR_INVALID_ARGUMENT,
             "NULL path rejected");
  check_true(csv_write_table(path, NULL, 3, data, 2, NULL) ==
                 CSV_ERR_INVALID_ARGUMENT,
             "NULL columns rejected");
  check_true(csv_write_table(path, cols, 0, data, 2, NULL) ==
                 CSV_ERR_INVALID_ARGUMENT,
             "zero columns rejected");
  check_true(csv_write_table(path, cols, 3, NULL, 2, NULL) ==
                 CSV_ERR_INVALID_ARGUMENT,
             "NULL data with rows > 0 rejected");

  const csv_column_t bad[1] = {{NULL, "Ha", NULL}};
  check_true(csv_write_table(path, bad, 1, data, 1, NULL) ==
                 CSV_ERR_INVALID_ARGUMENT,
             "column with NULL name rejected");

  check_true(csv_write_header_columns(NULL, cols, 3) ==
                 CSV_ERR_INVALID_ARGUMENT,
             "csv_write_header_columns(NULL writer) rejected");
}

static void test_bit_exact_round_trip(void) {
  printf("  === Test CSV bit-exact round trip (precision 17) ===\n");

  const double vals[10] = {0.1,
                           1.0 / 3.0,
                           -2.9037123456789012,
                           1e-300,
                           DBL_MAX,
                           DBL_MIN,
                           DBL_MIN / 4.0, /* subnormal */
                           3.14159265358979323846,
                           -0.0,
                           123456789.987654321};
  const csv_column_t col[1] = {{"v", NULL, NULL}};

  const char *path = QMC_OUTPUT_DIR "/test_csv_roundtrip.csv";
  csv_options_t opts = csv_options_default();
  opts.precision = 17;
  opts.scientific = 1;
  check_true(csv_write_table(path, col, 1, vals, 10, &opts) == CSV_OK,
             "write succeeds");

  FILE *f = fopen(path, "r");
  check_true(f != NULL, "file readable");
  if (!f) {
    return;
  }

  char line[128];
  int header_ok = (fgets(line, sizeof line, f) != NULL);
  int all_match = header_ok;
  for (size_t i = 0; i < 10 && all_match; i++) {
    if (!fgets(line, sizeof line, f)) {
      all_match = 0;

      break;
    }

    double back = strtod(line, NULL);
    // Intentional bit-exact comparison, incl. -0.0 vs 0.0 which == treats as
    // equal NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
    if (memcmp(&back, &vals[i], sizeof back) != 0) {
      all_match = 0;
    }
  }

  fclose(f);
  check_true(all_match,
             "every value (incl. -0.0, subnormal, DBL_MAX) reads back "
             "bit-for-bit identical");
}

int main(void) {
  (void)mkdir(QMC_OUTPUT_DIR, 0755);

  test_streaming_basic();
  test_escaping();
  test_delimiter_and_precision();
  test_nonfinite_policy();
  test_append_mode();
  test_null_and_error_handling();
  test_legacy_write_1d_and_matrix();
  test_column_descriptors_and_table();
  test_bit_exact_round_trip();


  if (failures == 0) {
    printf("\nAll test_csv_writer checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
