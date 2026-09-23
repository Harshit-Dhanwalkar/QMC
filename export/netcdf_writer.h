#ifndef QMC_NETCDF_WRITER_H
#define QMC_NETCDF_WRITER_H

#include <stddef.h>

// Both functions write to QMC_OUTPUT_DIR/filename
// Return 0 on success, -1 on failure

// Write a single 1D variable named `varname` inside a new NetCDF file
int netcdf_write_1d(const char *filename, const char *varname,
                    const double *data, size_t len);

// Write a single row-major 2D variable named `varname` inside a new NetCDF
// file, shaped (rows, cols)
int netcdf_write_matrix(const char *filename, const char *varname,
                        const double *data, size_t rows, size_t cols);

#endif // QMC_NETCDF_WRITER_H
