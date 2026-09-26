#ifndef QMC_HDF5_WRITER_H
#define QMC_HDF5_WRITER_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

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

/* Element type for the generic hdf5_write_dataset() entry point */
typedef enum {
  HDF5_TYPE_F32,
  HDF5_TYPE_F64,
  HDF5_TYPE_I32,
  HDF5_TYPE_I64
} hdf5_type_t;

/* Per-dataset chunking/compression knobs
 *
 * chunk_size == 0 means "no chunking" (dataset is written contiguously,
 * hdf5_writer_write_1d/matrix), in which case compression_level/shuffle are
 * ignored. Otherwise chunk_size is used as target chunk extent along every
 * dimension, clamped to that dimension's actual size
 */
typedef struct {
  int compression_level; /* 0 = none, 1-9 = gzip (deflate) level */
  size_t chunk_size;     /* elements per chunk dimension; 0 = no chunking */
  int shuffle;           /* 1 = enable the shuffle filter before deflate */
} hdf5_dataset_options_t;

/* Returns default options: no chunking, no compression, no shuffle */
hdf5_dataset_options_t hdf5_dataset_options_default(void);

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
 * Generic dataset writer: any rank, any of hdf5_type_t element types, optional
 * chunking/compression via `options` (NULL -> defaults, i.e. no
 * chunking/compression, matching hdf5_writer_write_1d/matrix)
 *
 * `dataset` may contain '/' to place dataset inside a group path
 * ("simulation/energy"); any missing intermediate groups are created
 * automatically for that prefix
 */
hdf5_status_t hdf5_write_dataset(hdf5_writer_t *writer, const char *dataset,
                                 hdf5_type_t type, int rank, const size_t *dims,
                                 const void *data,
                                 const hdf5_dataset_options_t *options);

/* Typed convenience wrappers around hdf5_write_dataset(), 1D case, default
 * (no chunking/compression) options */
hdf5_status_t hdf5_writer_write_f32_1d(hdf5_writer_t *writer,
                                       const char *dataset, const float *data,
                                       size_t len);
hdf5_status_t hdf5_writer_write_i32_1d(hdf5_writer_t *writer,
                                       const char *dataset, const int32_t *data,
                                       size_t len);
hdf5_status_t hdf5_writer_write_i64_1d(hdf5_writer_t *writer,
                                       const char *dataset, const int64_t *data,
                                       size_t len);

/* Create a group (and any missing intermediate groups) at `path`, e.g.
 * "simulation/observables". A no-op (returns HDF5_OK) if it already exists */
hdf5_status_t hdf5_create_group(hdf5_writer_t *writer, const char *path);

/* Attach a string/double attribute to a dataset or group. Use "/" as
 * object_path for file-level (root group) metadata, e.g. schema/library
 * version. Replaces an existing attribute of same name on that object */
hdf5_status_t hdf5_write_attribute_string(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name,
                                          const char *value);
hdf5_status_t hdf5_write_attribute_double(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name, double value);

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
