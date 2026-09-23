#ifndef QMC_HDF5_WRITER_H
#define QMC_HDF5_WRITER_H

#include <stddef.h>

// Both functions write to QMC_OUTPUT_DIR/filename
// Return 0 on success, -1 on failure

// Write a single 1D dataset named `dataset` inside a HDF5 file
int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len);

// Write a single row-major 2D dataset named `dataset` inside HDF5 file, shaped
// (rows, cols)
int hdf5_write_matrix(const char *filename, const char *dataset,
                      const double *data, size_t rows, size_t cols);

#endif // QMC_HDF5_WRITER_H
