/*
 * Test: HDF5 writer (export/hdf5_writer.c)
 *
 * Two build modes, both exercised same way from this one file:
 *  - Default (no USE_HDF5): hdf5_write_1d/hdf5_write_matrix must return -1
 *    without crashing
 *  - make USE_HDF5=1: writes must succeed and data must read back bit-for-bit
 *    identical via HDF5 C API (H5Fopen/H5Dread)
 */

#include "../export/hdf5_writer.h"

#include <math.h>
#include <stdio.h>
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

#ifdef USE_HDF5

#include <hdf5.h>

static void test_1d_round_trip(void) {
  printf("  === Test 1D round trip (USE_HDF5=1) ===\n");

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
  printf("  === Test matrix round trip (USE_HDF5=1) ===\n");

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

static void test_groups_and_generic_types(void) {
  printf("  === Test groups and generic_types (USE_HDF5=1) ===\n");

  hdf5_writer_t *w = hdf5_open(QMC_OUTPUT_DIR "/test_hdf5_groups.h5");
  check_true(w != NULL, "hdf5_open succeeds");
  if (!w) {
    return;
  }

  // Dataset paths with '/' should auto-create missing intermediate groups, no
  // explicit hdf5_create_group() call needed */
  const double energy[3] = {-2.90, -2.91, -2.93};
  hdf5_status_t rc = hdf5_writer_write_1d(w, "simulation/energy", energy, 3);
  check_true(rc == HDF5_OK,
             "hdf5_writer_write_1d into a nested group path succeeds");

  const float positions_f32[4] = {0.0f, 1.5f, 3.0f, -2.25f};
  rc = hdf5_writer_write_f32_1d(w, "particles/positions", positions_f32, 4);
  check_true(rc == HDF5_OK, "hdf5_writer_write_f32_1d succeeds");

  const int32_t counts_i32[3] = {10, 20, 30};
  rc = hdf5_writer_write_i32_1d(w, "counts", counts_i32, 3);
  check_true(rc == HDF5_OK, "hdf5_writer_write_i32_1d succeeds");

  const int64_t seeds_i64[2] = {123456789012LL, -42LL};
  rc = hdf5_writer_write_i64_1d(w, "metadata/seeds", seeds_i64, 2);
  check_true(rc == HDF5_OK, "hdf5_writer_write_i64_1d succeeds");

  // Explicit group creation, including a no-op re-creation
  rc = hdf5_create_group(w, "explicit_group");
  check_true(rc == HDF5_OK, "hdf5_create_group creates a new group");
  rc = hdf5_create_group(w, "explicit_group");
  check_true(rc == HDF5_OK,
             "hdf5_create_group on an existing group is a no-op success");

  check_true(hdf5_close(w) == HDF5_OK, "hdf5_close succeeds");

  hid_t file = H5Fopen(QMC_OUTPUT_DIR "/test_hdf5_groups.h5", H5F_ACC_RDONLY,
                       H5P_DEFAULT);
  check_true(file >= 0, "written file can be reopened");
  if (file < 0) {
    return;
  }

  check_true(H5Lexists(file, "simulation", H5P_DEFAULT) > 0,
             "group 'simulation' was auto-created");
  check_true(H5Lexists(file, "simulation/energy", H5P_DEFAULT) > 0,
             "dataset 'simulation/energy' exists");
  check_true(H5Lexists(file, "particles/positions", H5P_DEFAULT) > 0,
             "dataset 'particles/positions' exists");
  check_true(H5Lexists(file, "explicit_group", H5P_DEFAULT) > 0,
             "explicitly-created group 'explicit_group' exists");

  hid_t dset = H5Dopen2(file, "particles/positions", H5P_DEFAULT);
  check_true(dset >= 0, "'particles/positions' dataset opens");
  if (dset >= 0) {
    float readback[4] = {0};
    H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, readback);

    int all_match = 1;
    for (int i = 0; i < 4; i++) {
      if (readback[i] != positions_f32[i]) {
        all_match = 0;
      }
    }
    check_true(all_match, "f32 dataset reads back bit-for-bit identical");

    H5Dclose(dset);
  }

  dset = H5Dopen2(file, "metadata/seeds", H5P_DEFAULT);
  check_true(dset >= 0, "'metadata/seeds' dataset opens");
  if (dset >= 0) {
    int64_t readback[2] = {0};
    H5Dread(dset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, readback);
    check_true(readback[0] == seeds_i64[0] && readback[1] == seeds_i64[1],
               "i64 dataset reads back bit-for-bit identical");
    H5Dclose(dset);
  }

  H5Fclose(file);
}

static void test_attributes(void) {
  printf("  === Test attributes (USE_HDF5=1) ===\n");

  hdf5_writer_t *w = hdf5_open(QMC_OUTPUT_DIR "/test_hdf5_attrs.h5");
  check_true(w != NULL, "hdf5_open succeeds");
  if (!w) {
    return;
  }

  const double energy[3] = {-2.90, -2.91, -2.93};
  hdf5_writer_write_1d(w, "simulation/energy", energy, 3);

  hdf5_status_t rc =
      hdf5_write_attribute_string(w, "simulation/energy", "units", "Ha");
  check_true(rc == HDF5_OK, "hdf5_write_attribute_string on a dataset");

  rc = hdf5_write_attribute_double(w, "simulation/energy", "time_step", 0.01);
  check_true(rc == HDF5_OK, "hdf5_write_attribute_double on a dataset");

  // Root-group (file-level) metadata via "/"
  rc = hdf5_write_attribute_string(w, "/", "library_version", "0.9.0");
  check_true(rc == HDF5_OK, "hdf5_write_attribute_string on root group");

  // Overwriting an existing attribute should replace, not fail
  rc = hdf5_write_attribute_string(w, "simulation/energy", "units", "eV");
  check_true(rc == HDF5_OK, "re-writing an existing attribute succeeds");

  // Attribute on a nonexistent object should fail cleanly
  rc = hdf5_write_attribute_string(w, "does/not/exist", "units", "Ha");
  check_true(rc != HDF5_OK,
             "attribute write on a nonexistent object path fails cleanly");

  check_true(hdf5_close(w) == HDF5_OK, "hdf5_close succeeds");

  hid_t file = H5Fopen(QMC_OUTPUT_DIR "/test_hdf5_attrs.h5", H5F_ACC_RDONLY,
                       H5P_DEFAULT);
  check_true(file >= 0, "written file can be reopened");
  if (file < 0) {
    return;
  }

  hid_t dset = H5Dopen2(file, "simulation/energy", H5P_DEFAULT);
  check_true(dset >= 0, "dataset opens for attribute verification");
  if (dset >= 0) {
    hid_t attr = H5Aopen(dset, "units", H5P_DEFAULT);
    check_true(attr >= 0, "'units' attribute exists");
    if (attr >= 0) {
      hid_t atype = H5Aget_type(attr);
      size_t sz = H5Tget_size(atype);
      char buf[64] = {0};

      if (sz < sizeof buf) {
        H5Aread(attr, atype, buf);
      }

      check_true(strcmp(buf, "eV") == 0,
                 "'units' attribute reflects the overwritten value 'eV'");

      H5Tclose(atype);
      H5Aclose(attr);
    }

    attr = H5Aopen(dset, "time_step", H5P_DEFAULT);
    check_true(attr >= 0, "'time_step' attribute exists");
    if (attr >= 0) {
      double val = 0;
      H5Aread(attr, H5T_NATIVE_DOUBLE, &val);

      check_true(val == 0.01, "'time_step' attribute value matches");

      H5Aclose(attr);
    }

    H5Dclose(dset);
  }

  hid_t root = H5Gopen2(file, "/", H5P_DEFAULT);
  if (root >= 0) {
    hid_t attr = H5Aopen(root, "library_version", H5P_DEFAULT);
    check_true(attr >= 0, "root-group 'library_version' attribute exists");
    if (attr >= 0) {
      H5Aclose(attr);
    }

    H5Gclose(root);
  }

  H5Fclose(file);
}

static void test_chunking_and_compression(void) {
  printf("  === Test chunking and_ ompression (USE_HDF5=1) ===\n");

  hdf5_writer_t *w = hdf5_open(QMC_OUTPUT_DIR "/test_hdf5_chunked.h5");
  check_true(w != NULL, "hdf5_open succeeds");
  if (!w) {
    return;
  }

  enum { N = 200 };
  double series[N];
  for (int i = 0; i < N; i++) {
    series[i] = (double)i * 0.5;
  }

  hdf5_dataset_options_t opts = hdf5_dataset_options_default();
  opts.chunk_size = 32;
  opts.compression_level = 6;
  opts.shuffle = 1;

  const size_t dims[1] = {(size_t)N};
  hdf5_status_t rc =
      hdf5_write_dataset(w, "series", HDF5_TYPE_F64, 1, dims, series, &opts);
  check_true(rc == HDF5_OK, "hdf5_write_dataset with chunking+compression");

  // NULL options must fall back to defaults (no chunking) without error
  rc = hdf5_write_dataset(w, "series_plain", HDF5_TYPE_F64, 1, dims, series,
                          NULL);
  check_true(rc == HDF5_OK, "hdf5_write_dataset with NULL options (defaults)");

  check_true(hdf5_close(w) == HDF5_OK, "hdf5_close succeeds");

  hid_t file = H5Fopen(QMC_OUTPUT_DIR "/test_hdf5_chunked.h5", H5F_ACC_RDONLY,
                       H5P_DEFAULT);
  check_true(file >= 0, "written file can be reopened");
  if (file < 0) {
    return;
  }

  hid_t dset = H5Dopen2(file, "series", H5P_DEFAULT);
  check_true(dset >= 0, "chunked dataset opens");
  if (dset >= 0) {
    hid_t dcpl = H5Dget_create_plist(dset);
    check_true(H5Pget_layout(dcpl) == H5D_CHUNKED,
               "chunked dataset actually uses H5D_CHUNKED layout");

    double readback[N] = {0};
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, readback);
    int all_match = 1;
    for (int i = 0; i < N; i++) {
      if (readback[i] != series[i]) {
        all_match = 0;
      }
    }
    check_true(all_match,
               "chunked+compressed dataset reads back bit-for-bit identical");

    H5Pclose(dcpl);
    H5Dclose(dset);
  }

  H5Fclose(file);
}

#else // !USE_HDF5

static void test_apis_stub_without_crashing(void) {
  printf("  === Test apis stub without crashing (default build, no "
         "USE_HDF5) ===\n");

  hdf5_writer_t *w = hdf5_open("unused.h5");
  check_true(w == NULL, "hdf5_open stub returns NULL");

  const double data[3] = {1.0, 2.0, 3.0};
  const size_t dims[1] = {3};
  hdf5_status_t rc =
      hdf5_write_dataset(w, "x", HDF5_TYPE_F64, 1, dims, data, NULL);
  check_true(rc == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_dataset stub returns HDF5_ERR_NOT_SUPPORTED");

  rc = hdf5_create_group(w, "g");
  check_true(rc == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_create_group stub returns HDF5_ERR_NOT_SUPPORTED");

  rc = hdf5_write_attribute_string(w, "/", "units", "Ha");
  check_true(rc == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_attribute_string stub returns HDF5_ERR_NOT_SUPPORTED");

  rc = hdf5_write_attribute_double(w, "/", "time_step", 0.01);
  check_true(rc == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_attribute_double stub returns HDF5_ERR_NOT_SUPPORTED");
}

static void test_stub_returns_failure_without_crashing(void) {
  printf("  === Test stub returns failure without_crashing (default build, no "
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
  test_groups_and_generic_types();
  test_attributes();
  test_chunking_and_compression();
#else
  test_stub_returns_failure_without_crashing();
  test_apis_stub_without_crashing();
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
