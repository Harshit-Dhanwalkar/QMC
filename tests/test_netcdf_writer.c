/*
 * Test: NetCDF writer (export/netcdf_writer.c)
 *
 * Two build modes, both exercised same way from this one file:
 *  - Default (no USE_NETCDF): netcdf_write_1d/netcdf_write_matrix must
 *    return -1 without crashing
 *  - make USE_NETCDF=1: writes must succeed and data must read back
 *    bit-for-bit identical via NetCDF C API (nc_open/nc_get_var)
 */

#include "../export/netcdf_writer.h"

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

#ifdef USE_NETCDF

#include <netcdf.h>

static void test_1d_round_trip(void) {
  printf("test_1d_round_trip (USE_NETCDF=1):\n");

  double data[5] = {1.0, 2.5, -3.0, 0.0, 42.125};
  int rc = netcdf_write_1d("test_netcdf_1d.nc", "values", data, 5);
  check_true(rc == 0, "netcdf_write_1d returns 0 on success");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, "test_netcdf_1d.nc");

  int ncid;
  check_true(nc_open(path, NC_NOWRITE, &ncid) == NC_NOERR,
             "written file can be reopened with nc_open");

  int varid;
  check_true(nc_inq_varid(ncid, "values", &varid) == NC_NOERR,
             "variable 'values' exists in the file");

  double readback[5] = {0};
  check_true(nc_get_var_double(ncid, varid, readback) == NC_NOERR,
             "nc_get_var_double succeeds");

  int all_match = 1;
  for (int i = 0; i < 5; i++) {
    if (readback[i] != data[i]) {
      all_match = 0;
    }
  }
  check_true(all_match, "read-back values are bit-for-bit identical");

  nc_close(ncid);
}

static void test_matrix_round_trip(void) {
  printf("test_matrix_round_trip (USE_NETCDF=1):\n");

  // 2x3 row-major matrix
  double data[6] = {1, 2, 3, 4, 5, 6};
  int rc = netcdf_write_matrix("test_netcdf_matrix.nc", "grid", data, 2, 3);
  check_true(rc == 0, "netcdf_write_matrix returns 0 on success");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, "test_netcdf_matrix.nc");

  int ncid;
  nc_open(path, NC_NOWRITE, &ncid);

  int varid;
  nc_inq_varid(ncid, "grid", &varid);

  int dimids[2];
  nc_inq_vardimid(ncid, varid, dimids);
  size_t rows = 0, cols = 0;
  nc_inq_dimlen(ncid, dimids[0], &rows);
  nc_inq_dimlen(ncid, dimids[1], &cols);
  check_true(rows == 2 && cols == 3, "variable shape is (2, 3) as written");

  double readback[6] = {0};
  nc_get_var_double(ncid, varid, readback);

  int all_match = 1;
  for (int i = 0; i < 6; i++) {
    if (readback[i] != data[i]) {
      all_match = 0;
    }
  }
  check_true(all_match, "row-major matrix read back identical");

  nc_close(ncid);
}

#else // !USE_NETCDF

static void test_stub_returns_failure_without_crashing(void) {
  printf("test_stub_returns_failure_without_crashing (default build, no "
         "USE_NETCDF):\n");

  double data[3] = {1.0, 2.0, 3.0};
  int rc = netcdf_write_1d("test_netcdf_stub.nc", "values", data, 3);
  check_true(rc == -1, "netcdf_write_1d stub returns -1 (not a crash, not a "
                       "silently-empty file)");

  rc = netcdf_write_matrix("test_netcdf_stub.nc", "grid", data, 1, 3);
  check_true(rc == -1, "netcdf_write_matrix stub returns -1");

  rc = netcdf_write_1d(NULL, NULL, NULL, 0);
  check_true(rc == -1, "stub handles NULL/zero-length input without crashing");
}

#endif // USE_NETCDF

int main(void) {
#ifdef USE_NETCDF
  test_1d_round_trip();
  test_matrix_round_trip();
#else
  test_stub_returns_failure_without_crashing();
  printf("\n(built without USE_NETCDF - rebuild with 'make USE_NETCDF=1 "
         "build/test_netcdf_writer' to exercise the real read-back tests)\n");
#endif

  if (failures == 0) {
    printf("\nAll test_netcdf_writer checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
