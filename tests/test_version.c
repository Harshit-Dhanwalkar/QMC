/*
 * Test: library version mechanics
 *
 * 1. QMC_VERSION_STRING, MAJOR/MINOR/PATCH macros and QMC_VERSION_NUMBER agree
 *    with each other
 * 2. qmc_version()/qmc_version_number() report same values headers were
 *    compiled with
 * 3. QMC_LIBRARY_VERSION (written into JSON/HDF5 exports) is same string
 */
#include <stdio.h>
#include <string.h>

#include "../core/qmc_version.h"
#include "../export/export_schema.h"

static int failures = 0;

static void check(int cond, const char *what) {
  if (!cond) {
    fprintf(stderr, "FAIL: %s\n", what);
    failures++;
  }
}

int main(void) {
  int major = -1;
  int minor = -1;
  int patch = -1;
  char tail = 0;
  int n =
      sscanf(QMC_VERSION_STRING, "%d.%d.%d%c", &major, &minor, &patch, &tail);

  check(n == 3, "QMC_VERSION_STRING is exactly MAJOR.MINOR.PATCH");
  check(major == QMC_VERSION_MAJOR, "string major == QMC_VERSION_MAJOR");
  check(minor == QMC_VERSION_MINOR, "string minor == QMC_VERSION_MINOR");
  check(patch == QMC_VERSION_PATCH, "string patch == QMC_VERSION_PATCH");
  check(QMC_VERSION_NUMBER == major * 10000 + minor * 100 + patch,
        "QMC_VERSION_NUMBER == MAJOR*10000 + MINOR*100 + PATCH");

  check(strcmp(qmc_version(), QMC_VERSION_STRING) == 0,
        "qmc_version() == QMC_VERSION_STRING");
  check(qmc_version_number() == QMC_VERSION_NUMBER,
        "qmc_version_number() == QMC_VERSION_NUMBER");
  check(strcmp(QMC_LIBRARY_VERSION, QMC_VERSION_STRING) == 0,
        "QMC_LIBRARY_VERSION == QMC_VERSION_STRING");

  if (failures != 0) {
    fprintf(stderr, "test_version: %d check(s) failed\n", failures);
    return 1;
  }
  printf("All test_version checks passed (QMC %s).\n", qmc_version());

  return 0;
}
