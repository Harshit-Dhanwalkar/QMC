#ifndef QMC_HDF5_WRITER_H
#define QMC_HDF5_WRITER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  HDF5_OK = 0,
  HDF5_ERR_INVALID_ARGUMENT = -1,
  HDF5_ERR_OPEN = -2,
  HDF5_ERR_CREATE_DATASET = -3,
  HDF5_ERR_WRITE = -4,
  HDF5_ERR_CLOSE = -5,
  HDF5_ERR_NOT_SUPPORTED = -6,
  HDF5_ERR_PATH_TOO_LONG = -7
} hdf5_status_t;

/*
 * Persistent HDF5 writer
 *
 * hdf5_open() creates the file once (truncating if it already exists) and
 * subsequent hdf5_writer_write_*() calls add datasets to that same file without
 * touching previously-written data
 */
typedef struct hdf5_writer hdf5_writer_t;

/* Open (create, or truncate-and-recreate) `path` for writing
 *
 * Returns NULL on failure
 */
hdf5_writer_t *hdf5_open(const char *path);

/* Flush and close
 *
 * Returns HDF5_OK only if the file closed cleanly
 */
hdf5_status_t hdf5_close(hdf5_writer_t *writer);

/* Write native-double datasets into an open writer
 *
 * Dataset whose name already exists in the file is replaced; all other datasets
 * in the file are left untouched
 */
hdf5_status_t hdf5_writer_write_1d(hdf5_writer_t *writer, const char *dataset,
                                   const double *data, size_t len);
hdf5_status_t hdf5_writer_write_matrix(hdf5_writer_t *writer,
                                       const char *dataset, const double *data,
                                       size_t rows, size_t cols);

/*
 * Writes into QMC_OUTPUT_DIR/filename. On first call, creates file; on
 * subsequent calls to the same file, opens it in read-write mode and
 * adds/replaces only the named dataset. This means:
 *
 *    hdf5_write_1d("results.h5", "energy",   e, n);   // creates file
 *    hdf5_write_1d("results.h5", "variance", v, n);   // adds variance
 *
 * leaves a file containing both "energy" and "variance"
 *
 * Returns HDF5_OK (0) on success, a negative hdf5_status_t on failure
 */
int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len);
int hdf5_write_matrix(const char *filename, const char *dataset,
                      const double *data, size_t rows, size_t cols);

/* Description of a status code */
const char *hdf5_strerror(hdf5_status_t status);

#ifdef __cplusplus
}
#endif

#endif /* QMC_HDF5_WRITER_H */
