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
#include <stdint.h>
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

static hid_t open_ro(const char *path) {
  return H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
}

static void test_appendable_datasets(void) {
  printf("  === Test appendable datasets (USE_HDF5=1) ===\n");

  const char *path = QMC_OUTPUT_DIR "/test_hdf5_append.h5";
  hdf5_writer_t *w = hdf5_open(path);
  check_true(w != NULL, "hdf5_open succeeds");
  if (!w) {
    return;
  }

  // 1D f64, appended in uneven pieces incl. a zero-row no-op
  hdf5_status_t st = HDF5_ERR_INVALID_ARGUMENT;
  hdf5_dataset_options_t small_chunks = hdf5_dataset_options_default();
  small_chunks.chunk_size = 4;
  hdf5_dataset_t *e = hdf5_dataset_create(w, "observables/energy",
                                          HDF5_TYPE_F64, 1, &small_chunks, &st);
  check_true(e != NULL && st == HDF5_OK, "hdf5_dataset_create (1D) succeeds");
  if (!e) {
    hdf5_close(w);

    return;
  }

  double expect[13];
  for (int i = 0; i < 13; i++) {
    expect[i] = -2.9 + 0.001 * i;
  }
  check_true(hdf5_dataset_append(e, expect, 5) == HDF5_OK, "append 5 rows");
  check_true(hdf5_dataset_append(e, expect + 5, 1) == HDF5_OK, "append 1 row");
  check_true(hdf5_dataset_append(e, NULL, 0) == HDF5_OK,
             "zero-row append with NULL data is a no-op");
  check_true(hdf5_dataset_append(e, expect + 6, 7) == HDF5_OK, "append 7 rows");
  check_true(hdf5_dataset_rows(e) == 13, "row counter tracks appends");
  check_true(hdf5_dataset_append(e, NULL, 3) == HDF5_ERR_INVALID_ARGUMENT,
             "NULL data with nrows > 0 rejected");
  check_true(hdf5_dataset_rows(e) == 13, "rejected append leaves count alone");
  check_true(hdf5_flush(w) == HDF5_OK, "hdf5_flush succeeds mid-run");

  // 2D i32, row_width 3, compressed + shuffled
  hdf5_dataset_options_t comp = hdf5_dataset_options_default();
  comp.chunk_size = 2;
  comp.compression_level = 5;
  comp.shuffle = 1;
  hdf5_dataset_t *m =
      hdf5_dataset_create(w, "configs/counts", HDF5_TYPE_I32, 3, &comp, &st);
  check_true(m != NULL && st == HDF5_OK, "hdf5_dataset_create (2D) succeeds");

  int32_t rows_i32[24];
  for (int i = 0; i < 24; i++) {
    rows_i32[i] = i * 7 - 40;
  }
  if (m) {
    check_true(hdf5_dataset_append(m, rows_i32, 4) == HDF5_OK &&
                   hdf5_dataset_append(m, rows_i32 + 12, 4) == HDF5_OK,
               "append 2 x 4 rows of width 3");
    check_true(hdf5_dataset_rows(m) == 8, "2D row counter correct");
    check_true(hdf5_dataset_close(m) == HDF5_OK, "close 2D dataset");
  }

  // Bad arguments
  check_true(hdf5_dataset_create(w, "z", HDF5_TYPE_F64, 0, NULL, &st) == NULL &&
                 st == HDF5_ERR_INVALID_ARGUMENT,
             "row_width == 0 rejected, reason reported");
  check_true(hdf5_dataset_create(NULL, "z", HDF5_TYPE_F64, 1, NULL, NULL) ==
                 NULL,
             "NULL writer rejected (status_out may be NULL)");
  check_true(hdf5_dataset_rows(NULL) == 0, "rows(NULL) == 0");

  // Default options: 1024-row chunks
  hdf5_dataset_t *d =
      hdf5_dataset_create(w, "defaults", HDF5_TYPE_I64, 1, NULL, &st);
  check_true(d != NULL && st == HDF5_OK, "create with NULL options succeeds");

  if (d) {
    const int64_t big[3] = {INT64_MIN, 0, INT64_MAX};
    check_true(hdf5_dataset_append(d, big, 3) == HDF5_OK, "append i64 rows");
    check_true(hdf5_dataset_close(d) == HDF5_OK, "close i64 dataset");
  }
  check_true(hdf5_dataset_close(e) == HDF5_OK, "close 1D dataset");
  check_true(hdf5_close(w) == HDF5_OK, "hdf5_close succeeds");

  // read everything back independently through the HDF5 C API
  hid_t file = open_ro(path);
  check_true(file >= 0, "file reopens read-only");
  if (file < 0) {
    return;
  }

  hid_t ds = H5Dopen2(file, "observables/energy", H5P_DEFAULT);
  hid_t sp = H5Dget_space(ds);
  hsize_t dims[2] = {0, 0}, maxd[2] = {0, 0};
  H5Sget_simple_extent_dims(sp, dims, maxd);
  check_true(H5Sget_simple_extent_ndims(sp) == 1 && dims[0] == 13,
             "1D dataset has shape (13,)");
  check_true(maxd[0] == H5S_UNLIMITED, "1D dataset is extendible (unlimited)");

  hid_t dcpl = H5Dget_create_plist(ds);
  check_true(H5Pget_layout(dcpl) == H5D_CHUNKED, "1D dataset is chunked");

  hsize_t ck[2] = {0, 0};
  H5Pget_chunk(dcpl, 1, ck);
  check_true(ck[0] == 4, "options->chunk_size is rows-per-chunk");

  double got[13];
  H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, got);

  // NOTE: Intentional bit-exact comparison of the whole array as written/read
  // back NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
  check_true(memcmp(got, expect, sizeof got) == 0,
             "1D data reads back bit-for-bit, in append order");
  H5Pclose(dcpl);
  H5Sclose(sp);
  H5Dclose(ds);

  ds = H5Dopen2(file, "configs/counts", H5P_DEFAULT);
  sp = H5Dget_space(ds);
  H5Sget_simple_extent_dims(sp, dims, maxd);
  check_true(H5Sget_simple_extent_ndims(sp) == 2 && dims[0] == 8 &&
                 dims[1] == 3,
             "2D dataset has shape (8, 3)");
  check_true(maxd[0] == H5S_UNLIMITED && maxd[1] == 3,
             "only the row axis is extendible");

  int32_t got2[24];
  H5Dread(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, got2);
  int ok2 = (memcmp(got2, rows_i32, 12 * sizeof(int32_t)) == 0) &&
            (memcmp(got2 + 12, rows_i32 + 12, 12 * sizeof(int32_t)) == 0);
  check_true(ok2, "2D compressed data reads back identical, row-major");

  dcpl = H5Dget_create_plist(ds);
  check_true(H5Pget_nfilters(dcpl) == 2, "shuffle + deflate filters attached");

  H5Pclose(dcpl);
  H5Sclose(sp);
  H5Dclose(ds);

  ds = H5Dopen2(file, "defaults", H5P_DEFAULT);
  dcpl = H5Dget_create_plist(ds);
  H5Pget_chunk(dcpl, 1, ck);
  check_true(ck[0] == 1024, "default appendable chunk is 1024 rows");

  int64_t g64[3];
  H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, g64);
  check_true(g64[0] == INT64_MIN && g64[1] == 0 && g64[2] == INT64_MAX,
             "i64 extremes survive");

  H5Pclose(dcpl);
  H5Dclose(ds);
  H5Fclose(file);
}

static void test_strings_int_attrs_and_schema_header(void) {
  printf("  === Test string dataset, int attribute, schema header ===\n");

  const char *path = QMC_OUTPUT_DIR "/test_hdf5_schema.h5";
  hdf5_writer_t *w = hdf5_open(path);
  check_true(w != NULL, "hdf5_open succeeds");
  if (!w) {
    return;
  }

  check_true(hdf5_write_schema_header(w, 0) == HDF5_OK,
             "schema header (no timestamp) succeeds");
  check_true(hdf5_write_string(w, "metadata/method", "DMC") == HDF5_OK,
             "hdf5_write_string into nested group succeeds");
  check_true(hdf5_write_string(w, "metadata/method", "VMC longer name") ==
                 HDF5_OK,
             "hdf5_write_string replaces an existing dataset");
  check_true(hdf5_create_group(w, "g") == HDF5_OK, "group created");
  check_true(hdf5_write_attribute_int(w, "g", "n_walkers", 200) == HDF5_OK,
             "hdf5_write_attribute_int succeeds");
  check_true(hdf5_write_attribute_int(w, "g", "n_walkers", 400) == HDF5_OK,
             "hdf5_write_attribute_int replaces existing attribute");
  check_true(hdf5_write_string(NULL, "s", "x") == HDF5_ERR_INVALID_ARGUMENT &&
                 hdf5_write_string(w, "s", NULL) == HDF5_ERR_INVALID_ARGUMENT,
             "hdf5_write_string rejects NULL args");
  check_true(hdf5_write_schema_header(NULL, 0) == HDF5_ERR_INVALID_ARGUMENT,
             "schema header rejects NULL writer");
  check_true(hdf5_close(w) == HDF5_OK, "hdf5_close succeeds");

  hid_t file = open_ro(path);
  check_true(file >= 0, "file reopens read-only");
  if (file < 0) {
    return;
  }

  hid_t root = H5Gopen2(file, "/", H5P_DEFAULT);
  long ver = -1;
  hid_t a = H5Aopen(root, "schema_version", H5P_DEFAULT);
  if (a >= 0) {
    H5Aread(a, H5T_NATIVE_LONG, &ver);

    H5Aclose(a);
  }
  check_true(ver == QMC_EXPORT_SCHEMA_VERSION, "root schema_version == 1");

  char buf[64] = {0};
  a = H5Aopen(root, "library", H5P_DEFAULT);
  hid_t at = (a >= 0) ? H5Aget_type(a) : -1;
  if (a >= 0) {
    H5Aread(a, at, buf);

    H5Tclose(at);
    H5Aclose(a);
  }
  check_true(strcmp(buf, QMC_LIBRARY_NAME) == 0, "root library attribute");

  memset(buf, 0, sizeof buf);
  a = H5Aopen(root, "library_version", H5P_DEFAULT);
  at = (a >= 0) ? H5Aget_type(a) : -1;
  if (a >= 0) {
    H5Aread(a, at, buf);

    H5Tclose(at);
    H5Aclose(a);
  }
  check_true(strcmp(buf, QMC_LIBRARY_VERSION) == 0,
             "root library_version attribute");
  check_true(H5Aexists(root, "created") == 0,
             "no created attribute when timestamp disabled");
  H5Gclose(root);

  hid_t ds = H5Dopen2(file, "metadata/method", H5P_DEFAULT);
  memset(buf, 0, sizeof buf);
  if (ds >= 0) {
    hid_t t = H5Dget_type(ds);
    H5Dread(ds, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

    H5Tclose(t);
    H5Dclose(ds);
  }
  check_true(strcmp(buf, "VMC longer name") == 0,
             "string dataset holds the replacement value");

  hid_t g = H5Gopen2(file, "g", H5P_DEFAULT);
  long nw = -1;
  a = H5Aopen(g, "n_walkers", H5P_DEFAULT);
  if (a >= 0) {
    H5Aread(a, H5T_NATIVE_LONG, &nw);

    H5Aclose(a);
  }
  check_true(nw == 400, "int attribute holds the replacement value");

  H5Gclose(g);
  H5Fclose(file);

  // Timestamp variant
  w = hdf5_open(path);
  check_true(w && hdf5_write_schema_header(w, 1) == HDF5_OK,
             "schema header (with timestamp) succeeds");
  if (w) {
    hdf5_close(w);
  }

  file = open_ro(path);
  root = (file >= 0) ? H5Gopen2(file, "/", H5P_DEFAULT) : -1;
  memset(buf, 0, sizeof buf);
  if (root >= 0 && H5Aexists(root, "created") > 0) {
    a = H5Aopen(root, "created", H5P_DEFAULT);
    at = H5Aget_type(a);
    H5Aread(a, at, buf);

    H5Tclose(at);
    H5Aclose(a);
  }
  check_true(strlen(buf) == 20 && buf[4] == '-' && buf[10] == 'T' &&
                 buf[19] == 'Z',
             "created is ISO 8601 UTC");

  if (root >= 0) {
    H5Gclose(root);
  }
  if (file >= 0) {
    H5Fclose(file);
  }
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

  const hdf5_dataset_t *ds =
      hdf5_dataset_create(w, "x", HDF5_TYPE_F64, 1, NULL, &rc);
  check_true(ds == NULL && rc == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_dataset_create stub returns NULL + NOT_SUPPORTED");
  check_true(hdf5_dataset_append(NULL, data, 1) == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_dataset_append stub returns HDF5_ERR_NOT_SUPPORTED");
  check_true(hdf5_dataset_rows(NULL) == 0, "hdf5_dataset_rows stub returns 0");
  check_true(hdf5_dataset_close(NULL) == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_dataset_close stub returns HDF5_ERR_NOT_SUPPORTED");
  check_true(hdf5_write_string(w, "s", "v") == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_string stub returns HDF5_ERR_NOT_SUPPORTED");
  check_true(hdf5_write_attribute_int(w, "/", "n", 1) == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_attribute_int stub returns HDF5_ERR_NOT_SUPPORTED");
  check_true(hdf5_write_schema_header(w, 0) == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_write_schema_header stub returns HDF5_ERR_NOT_SUPPORTED");
  check_true(hdf5_flush(w) == HDF5_ERR_NOT_SUPPORTED,
             "hdf5_flush stub returns HDF5_ERR_NOT_SUPPORTED");
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
  test_appendable_datasets();
  test_strings_int_attrs_and_schema_header();
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
