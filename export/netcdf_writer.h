#ifndef QMC_NETCDF_WRITER_H
#define QMC_NETCDF_WRITER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  NETCDF_OK = 0,
  NETCDF_ERR_INVALID_ARGUMENT = -1,
  NETCDF_ERR_OPEN = -2,
  NETCDF_ERR_CREATE_VARIABLE = -3,
  NETCDF_ERR_WRITE = -4,
  NETCDF_ERR_CLOSE = -5,
  NETCDF_ERR_NOT_SUPPORTED = -6,
  NETCDF_ERR_PATH_TOO_LONG = -7
} netcdf_status_t;

/*
 * Persistent NetCDF writer
 *
 * Preferred over the legacy one-shot API: netcdf_open() creates file once
 * (truncating if it already exists) and subsequent () calls add variables to
 * that same file without touching previously-written variables
 *
 * NOTE: Files are written in NETCDF4 format (HDF5-backed). This is required so
 * a variable can be replaced if a same-named variable is written again (classic
 * NetCDF3 does not permit variable deletion). All modern readers (netCDF4
 * Python, xarray, NCO, panoply, ...) read NETCDF4 files
 */
typedef struct netcdf_writer netcdf_writer_t;

netcdf_writer_t *netcdf_open(const char *path);
netcdf_status_t netcdf_close(netcdf_writer_t *writer);

netcdf_status_t netcdf_writer_write_1d(netcdf_writer_t *writer,
                                       const char *varname, const double *data,
                                       size_t len);
netcdf_status_t netcdf_writer_write_matrix(netcdf_writer_t *writer,
                                           const char *varname,
                                           const double *data, size_t rows,
                                           size_t cols);

/*
 * Writes into QMC_OUTPUT_DIR/filename. On the first call for a given
 * file, creates the file; on subsequent calls to the same file, opens it
 * in write mode and adds/replaces only the named variable. This means:
 *
 *     netcdf_write_1d("results.nc", "energy",   e, n);   // creates file
 *     netcdf_write_1d("results.nc", "variance", v, n);   // adds variance
 *
 * NOTE: leaves a file containing both "energy" and "variance"
 *
 * Returns NETCDF_OK (0) on success, a negative netcdf_status_t on failure
 */
int netcdf_write_1d(const char *filename, const char *varname,
                    const double *data, size_t len);
int netcdf_write_matrix(const char *filename, const char *varname,
                        const double *data, size_t rows, size_t cols);

const char *netcdf_strerror(netcdf_status_t status);

#ifdef __cplusplus
}
#endif

#endif /* QMC_NETCDF_WRITER_H */
