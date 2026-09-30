#ifndef QMC_HDF5_WRITER_H
#define QMC_HDF5_WRITER_H

#include <stddef.h>
#include <stdint.h>

#include "export_schema.h"

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

/* Element type for generic hdf5_write_dataset() entry point */
typedef enum {
  HDF5_TYPE_F32,
  HDF5_TYPE_F64,
  HDF5_TYPE_I32,
  HDF5_TYPE_I64
} hdf5_type_t;

/* Per-dataset chunking/compression knobs
 *
 * chunk_size == 0 means "no explicit chunking". If compression_level > 0 or
 * shuffle is requested with chunk_size == 0, an implementation-chosen default
 * chunk extent (1024 elements per axis, clamped to dataset dimension) is used
 * instead - HDF5 requires a chunked layout for filters. Otherwise dataset is
 * written contiguously and compression/shuffle are ignored
 * When chunk_size > 0, it is used as target chunk extent along every dimension,
 * clamped to that dimension's actual size
 */
typedef struct {
  int compression_level; /* 0 = none, 1-9 = gzip (deflate) level */
  size_t chunk_size;     /* elements per chunk dimension; 0 = no chunking */
  int shuffle;           /* 1 = enable shuffle filter before deflate */
} hdf5_dataset_options_t;

/* Returns default options: no chunking, no compression, no shuffle */
hdf5_dataset_options_t hdf5_dataset_options_default(void);

/*
 * Persistent HDF5 writer
 *
 * hdf5_open() creates file once (truncating if it already exists) and
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
 * Returns HDF5_OK only if file closed cleanly
 */
hdf5_status_t hdf5_close(hdf5_writer_t *writer);

/* Write native-double datasets into an open writer
 *
 * Dataset whose name already exists in file is replaced; all odatasets
 * in file are left untouched
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

/* Integer attribute (same object_path rules as string/double variants) */
hdf5_status_t hdf5_write_attribute_int(hdf5_writer_t *writer,
                                       const char *object_path,
                                       const char *attr_name, long value);

/* Scalar (single-value) string dataset, e.g. hdf5_write_string(w,
 * "metadata/method", "DMC"). Replaces an existing dataset of same name;
 *  missing intermediate groups are created */
hdf5_status_t hdf5_write_string(hdf5_writer_t *writer, const char *dataset,
                                const char *value);

/*
 * Stamps file-level provenance as root attributes, mirroring JSON header:
 *   schema_version (int), library, library_version, [created (UTC ISO 8601)]
 *  `schema_version` is *file layout* version and is independent of library
 *   version. Pass include_timestamp = 0 for reproducible files
 */
hdf5_status_t hdf5_write_schema_header(hdf5_writer_t *writer,
                                       int include_timestamp);

/* Push buffered data to disk without closing - worth calling periodically in
 * long simulations so a crash keeps everything written so far */
hdf5_status_t hdf5_flush(hdf5_writer_t *writer);

/*
 * Appendable ("extendible") dataset: grows along its first axis as rows arrive,
 * so a long run can stream observables instead of buffering them all
 *
 *   hdf5_dataset_t *e = hdf5_dataset_create(w, "observables/energy",
 *                                           HDF5_TYPE_F64, 1, NULL, NULL);
 *   for (...) hdf5_dataset_append(e, &energy_i, 1);
 *   hdf5_dataset_close(e);   // close every dataset BEFORE hdf5_close(w)
 *
 * row_width == 1 gives a 1D dataset of shape (rows,); row_width > 1 gives a 2D
 * dataset of shape (rows, row_width) and each appended row supplies row_width
 * consecutive elements. Appendable datasets are always chunked: for them
 * options->chunk_size means *rows per chunk* (default 1024 when 0/NULL);
 * compression_level and shuffle apply as usual. An existing dataset of same
 * name is replaced. On failure returns NULL and, if status_out != NULL, stores
 * reason there
 */
typedef struct hdf5_dataset hdf5_dataset_t;

hdf5_dataset_t *hdf5_dataset_create(hdf5_writer_t *writer, const char *dataset,
                                    hdf5_type_t type, size_t row_width,
                                    const hdf5_dataset_options_t *options,
                                    hdf5_status_t *status_out);

/* Appends `nrows` rows (nrows * row_width elements of dataset's type).
 * nrows == 0 is a valid no-op (data may then be NULL) */
hdf5_status_t hdf5_dataset_append(hdf5_dataset_t *ds, const void *data,
                                  size_t nrows);

/* Rows written so far (0 for a NULL handle) */
size_t hdf5_dataset_rows(const hdf5_dataset_t *ds);

/* Release handle; data stays in file */
hdf5_status_t hdf5_dataset_close(hdf5_dataset_t *ds);

/*
 * Writes into QMC_OUTPUT_DIR/filename. On first call, creates file; on
 * subsequent calls to same file, opens it in read-write mode and
 * adds/replaces only named dataset. This means:
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
