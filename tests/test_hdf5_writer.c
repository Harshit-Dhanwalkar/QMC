/*
 * Test: HDF5 writer (export/hdf5_writer.c)
 *
 * Two build modes, both exercised same way from this one file:
 *  - Default (no USE_HDF5): hdf5_write_1d/hdf5_write_matrix must return -1
 *    without crashing
 *  - make USE_HDF5=1: writes must succeed and data must read back
 *    bit-for-bit identical via real HDF5 C API (H5Fopen/H5Dread)
 */

#include "../export/hdf5_writer.h"

#include <math.h>
#include <stdio.h>

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

#ifdef USE_HDF5

#include <hdf5.h>

static void test_1d_round_trip(void) {
  printf("  === Test 1d_round_trip (USE_HDF5=1) ===\n");

  const double data[5] = {1.0, 2.5, -3.0, 0.0, 42.125};
  int rc = hdf5_write_1d("test_hdf5_1d.h5", "values", data, 5);
  check_true(rc == 0, "hdf5_write_1d returns 0 on success");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, "test_hdf5_1d.h5");

  hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
  check_true(file >= 0, "written file can be reopened with H5Fopen");
  if (file < 0) {
    return;
  }

  hid_t dset = H5Dopen2(file, "values", H5P_DEFAULT);
  check_true(dset >= 0, "dataset 'values' exists in the file");

  double readback[5] = {0};
  herr_t status =
      H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, readback);
  check_true(status >= 0, "H5Dread succeeds");

  int all_match = 1;
  for (int i = 0; i < 5; i++) {
    if (readback[i] != data[i]) {
      all_match = 0;
    }
  }
  check_true(all_match, "read-back values are bit-for-bit identical");

  H5Dclose(dset);
  H5Fclose(file);
}

static void test_matrix_round_trip(void) {
  printf("  === Test matrix_round_trip (USE_HDF5=1) ===\n");

  // 2x3 row-major matrix
  const double data[6] = {1, 2, 3, 4, 5, 6};
  int rc = hdf5_write_matrix("test_hdf5_matrix.h5", "grid", data, 2, 3);
  check_true(rc == 0, "hdf5_write_matrix returns 0 on success");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, "test_hdf5_matrix.h5");

  hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dset = H5Dopen2(file, "grid", H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);

  hsize_t dims[2] = {0, 0};
  H5Sget_simple_extent_dims(space, dims, NULL);
  check_true(dims[0] == 2 && dims[1] == 3,
             "dataset shape is (2, 3) as written");

  double readback[6] = {0};
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, readback);

  int all_match = 1;
  for (int i = 0; i < 6; i++) {
    if (readback[i] != data[i]) {
      all_match = 0;
    }
  }
  check_true(all_match, "row-major matrix read back identical");

  H5Sclose(space);
  H5Dclose(dset);
  H5Fclose(file);
}

#else // !USE_HDF5

static void test_stub_returns_failure_without_crashing(void) {
  printf("  === Test stub_returns_failure_without_crashing (default build, no "
         "USE_HDF5) ===\n");

  const double data[3] = {1.0, 2.0, 3.0};
  int rc = hdf5_write_1d("test_hdf5_stub.h5", "values", data, 3);
  check_true(rc == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_1d stub returns HDF5_ERR_NOT_SUPPORTED (not a "
             "crash, not a silently-empty file)");

  rc = hdf5_write_matrix("test_hdf5_stub.h5", "grid", data, 1, 3);
  check_true(rc == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_matrix stub returns HDF5_ERR_NOT_SUPPORTED");

  // NULL/zero-length inputs must not crash stub either
  rc = hdf5_write_1d(NULL, NULL, NULL, 0);
  check_true(rc == HDF5_ERR_NOT_SUPPORTED,
             "stub handles NULL/zero-length input without crashing");
}

#endif // USE_HDF5

int main(void) {
#ifdef USE_HDF5
  test_1d_round_trip();
  test_matrix_round_trip();
#else
  test_stub_returns_failure_without_crashing();
  printf("\n(built without USE_HDF5 - rebuild with 'make USE_HDF5=1 "
         "build/test_hdf5_writer' to exercise read-back tests)\n");
#endif

  if (failures == 0) {
    printf("\nAll test_hdf5_writer checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
