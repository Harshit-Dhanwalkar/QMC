#include "netcdf_writer.h"

#include <stdio.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

#ifdef USE_NETCDF

#include <netcdf.h>

// NOTE: nc_* calls return NC_NOERR (0) on success and a negative error code
// otherwise; check() centralizes the "print and bail" handling so each
// call site stays a single line
static int check(int status, const char *what) {
  if (status != NC_NOERR) {
    fprintf(stderr, "netcdf_writer: %s: %s\n", what, nc_strerror(status));
    return -1;
  }

  return 0;
}

static int write_variable(const char *filename, const char *varname,
                          const double *data, int ndims, const size_t *shape) {
  if (!filename || !varname || !data) {
    fprintf(stderr, "netcdf_writer: NULL filename/varname/data\n");
    return -1;
  }

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);

  int ncid;
  if (check(nc_create(path, NC_CLOBBER, &ncid), "nc_create") != 0) {
    return -1;
  }

  int dimids[2];
  static const char *dim_names[2] = {"dim0", "dim1"};
  for (int i = 0; i < ndims; i++) {
    if (check(nc_def_dim(ncid, dim_names[i], shape[i], &dimids[i]),
              "nc_def_dim") != 0) {
      nc_close(ncid);

      return -1;
    }
  }

  int varid;
  if (check(nc_def_var(ncid, varname, NC_DOUBLE, ndims, dimids, &varid),
            "nc_def_var") != 0) {
    nc_close(ncid);

    return -1;
  }

  if (check(nc_enddef(ncid), "nc_enddef") != 0) {
    nc_close(ncid);

    return -1;
  }

  if (check(nc_put_var_double(ncid, varid, data), "nc_put_var_double") != 0) {
    nc_close(ncid);

    return -1;
  }

  return check(nc_close(ncid), "nc_close");
}

int netcdf_write_1d(const char *filename, const char *varname,
                    const double *data, size_t len) {
  const size_t shape[1] = {len};

  return write_variable(filename, varname, data, 1, shape);
}

int netcdf_write_matrix(const char *filename, const char *varname,
                        const double *data, size_t rows, size_t cols) {
  const size_t shape[2] = {rows, cols};

  return write_variable(filename, varname, data, 2, shape);
}

#else // !USE_NETCDF

// Built without USE_NETCDF
int netcdf_write_1d(const char *filename, const char *varname,
                    const double *data, size_t len) {
  (void)filename;
  (void)varname;
  (void)data;
  (void)len;

  fprintf(stderr, "netcdf_writer: not built with NetCDF support (rebuild with "
                  "'make USE_NETCDF=1', which requires libnetcdf-dev)\n");

  return -1;
}

int netcdf_write_matrix(const char *filename, const char *varname,
                        const double *data, size_t rows, size_t cols) {
  (void)filename;
  (void)varname;
  (void)data;
  (void)rows;
  (void)cols;

  fprintf(stderr, "netcdf_writer: not built with NetCDF support (rebuild with "
                  "'make USE_NETCDF=1', which requires libnetcdf-dev)\n");

  return -1;
}

#endif // USE_NETCDF
