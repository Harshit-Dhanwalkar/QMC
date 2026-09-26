#include "hdf5_writer.h"

#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

#define HDF5_MAX_PATH 4096

const char *hdf5_strerror(hdf5_status_t status) {
  switch (status) {
  case HDF5_OK:
    return "ok";
  case HDF5_ERR_INVALID_ARGUMENT:
    return "invalid argument";
  case HDF5_ERR_OPEN:
    return "could not open/create file";
  case HDF5_ERR_CREATE_DATASET:
    return "could not create dataset";
  case HDF5_ERR_WRITE:
    return "dataset write failed";
  case HDF5_ERR_CLOSE:
    return "file close failed";
  case HDF5_ERR_NOT_SUPPORTED:
    return "not built with HDF5 support";
  case HDF5_ERR_PATH_TOO_LONG:
    return "path too long";
  default:
    return "unknown error";
  }
}

hdf5_dataset_options_t hdf5_dataset_options_default(void) {
  hdf5_dataset_options_t o = {
      .compression_level = 0,
      .chunk_size = 0,
      .shuffle = 0,
  };

  return o;
}

#ifdef USE_HDF5

#include <hdf5.h>

#define HDF5_MAX_RANK 8

struct hdf5_writer {
  hid_t file;
};

static hid_t hdf5_native_type(hdf5_type_t type) {
  switch (type) {
  case HDF5_TYPE_F32:
    return H5T_NATIVE_FLOAT;
  case HDF5_TYPE_I32:
    return H5T_NATIVE_INT32;
  case HDF5_TYPE_I64:
    return H5T_NATIVE_INT64;
  case HDF5_TYPE_F64:
  default:
    return H5T_NATIVE_DOUBLE;
  }
}

/* Create any missing intermediate groups implied by `dataset_path` ("a/b/c"
 * needs groups "a" and "a/b" to exist before dataset "a/b/c" can be created) */
static hdf5_status_t ensure_parent_groups(hid_t file,
                                          const char *dataset_path) {
  const char *last_slash = strrchr(dataset_path, '/');
  if (!last_slash) {
    return HDF5_OK;
  }

  size_t parent_len = (size_t)(last_slash - dataset_path);
  if (parent_len == 0) {
    return HDF5_OK; /* leading "/name" - root group already exists */
  }

  char work[HDF5_MAX_PATH];
  if (parent_len >= sizeof work) {
    return HDF5_ERR_PATH_TOO_LONG;
  }

  memcpy(work, dataset_path, parent_len);
  work[parent_len] = '\0';

  char cumulative[HDF5_MAX_PATH];
  cumulative[0] = '\0';

  char *saveptr = NULL;
  const char *tok = strtok_r(work, "/", &saveptr);
  while (tok) {
    size_t cur_len = strlen(cumulative);
    size_t tok_len = strlen(tok);
    if (cur_len + 1 + tok_len + 1 > sizeof cumulative) {
      return HDF5_ERR_PATH_TOO_LONG;
    }

    if (cur_len > 0) {
      cumulative[cur_len] = '/';
      cur_len++;
    }

    memcpy(cumulative + cur_len, tok, tok_len + 1); /* +1 copies NUL */

    if (H5Lexists(file, cumulative, H5P_DEFAULT) <= 0) {
      hid_t g =
          H5Gcreate2(file, cumulative, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (g < 0) {
        return HDF5_ERR_CREATE_DATASET;
      }
      H5Gclose(g);
    }

    tok = strtok_r(NULL, "/", &saveptr);
  }

  return HDF5_OK;
}

/* Write (or replace) a single dataset in an already-open file */
static hdf5_status_t write_dataset_impl(hdf5_writer_t *w, const char *dataset,
                                        const double *data, int rank,
                                        const hsize_t *dims) {
  if (!w || !dataset || !data) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  hdf5_status_t grc = ensure_parent_groups(w->file, dataset);
  if (grc != HDF5_OK) {
    return grc;
  }

  // Replace if a dataset of this name already exists
  if (H5Lexists(w->file, dataset, H5P_DEFAULT) > 0) {
    if (H5Ldelete(w->file, dataset, H5P_DEFAULT) < 0) {
      return HDF5_ERR_CREATE_DATASET;
    }
  }

  hid_t space = H5Screate_simple(rank, dims, NULL);
  if (space < 0) {
    return HDF5_ERR_CREATE_DATASET;
  }

  hid_t dset = H5Dcreate2(w->file, dataset, H5T_NATIVE_DOUBLE, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) {
    H5Sclose(space);

    return HDF5_ERR_CREATE_DATASET;
  }

  herr_t status =
      H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dset);
  H5Sclose(space);

  return (status < 0) ? HDF5_ERR_WRITE : HDF5_OK;
}

hdf5_writer_t *hdf5_open(const char *path) {
  if (!path || path[0] == '\0') {
    return NULL;
  }

  hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    return NULL;
  }

  hdf5_writer_t *w = malloc(sizeof *w);
  if (!w) {
    H5Fclose(file);

    return NULL;
  }

  w->file = file;

  return w;
}

hdf5_status_t hdf5_close(hdf5_writer_t *writer) {
  if (!writer) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  hdf5_status_t rc = HDF5_OK;
  if (writer->file >= 0) {
    if (H5Fclose(writer->file) < 0) {
      rc = HDF5_ERR_CLOSE;
    }
  }

  free(writer);

  return rc;
}

hdf5_status_t hdf5_writer_write_1d(hdf5_writer_t *writer, const char *dataset,
                                   const double *data, size_t len) {
  const hsize_t dims[1] = {(hsize_t)len};

  return write_dataset_impl(writer, dataset, data, 1, dims);
}

hdf5_status_t hdf5_writer_write_matrix(hdf5_writer_t *writer,
                                       const char *dataset, const double *data,
                                       size_t rows, size_t cols) {
  const hsize_t dims[2] = {(hsize_t)rows, (hsize_t)cols};

  return write_dataset_impl(writer, dataset, data, 2, dims);
}

hdf5_status_t hdf5_write_dataset(hdf5_writer_t *w, const char *dataset,
                                 hdf5_type_t type, int rank, const size_t *dims,
                                 const void *data,
                                 const hdf5_dataset_options_t *options) {
  if (!w || !dataset || !dims || !data || rank <= 0 || rank > HDF5_MAX_RANK) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  hdf5_status_t grc = ensure_parent_groups(w->file, dataset);
  if (grc != HDF5_OK) {
    return grc;
  }

  if (H5Lexists(w->file, dataset, H5P_DEFAULT) > 0) {
    if (H5Ldelete(w->file, dataset, H5P_DEFAULT) < 0) {
      return HDF5_ERR_CREATE_DATASET;
    }
  }

  hsize_t h5dims[HDF5_MAX_RANK];
  for (int i = 0; i < rank; i++) {
    h5dims[i] = (hsize_t)dims[i];
  }

  hid_t space = H5Screate_simple(rank, h5dims, NULL);
  if (space < 0) {
    return HDF5_ERR_CREATE_DATASET;
  }

  hdf5_dataset_options_t opts =
      options ? *options : hdf5_dataset_options_default();

  hid_t dcpl = H5P_DEFAULT;
  int owns_dcpl = 0;
  if (opts.chunk_size > 0) {
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0) {
      H5Sclose(space);

      return HDF5_ERR_CREATE_DATASET;
    }

    owns_dcpl = 1;

    hsize_t chunk_dims[HDF5_MAX_RANK];
    for (int i = 0; i < rank; i++) {
      hsize_t c = (hsize_t)opts.chunk_size;
      chunk_dims[i] = (c < h5dims[i]) ? c : h5dims[i];
      if (chunk_dims[i] == 0) {
        chunk_dims[i] = 1; /* HDF5 rejects a zero-sized chunk dimension */
      }
    }

    if (H5Pset_chunk(dcpl, rank, chunk_dims) < 0) {
      H5Pclose(dcpl);
      H5Sclose(space);

      return HDF5_ERR_CREATE_DATASET;
    }

    if (opts.shuffle && H5Pset_shuffle(dcpl) < 0) {
      H5Pclose(dcpl);
      H5Sclose(space);

      return HDF5_ERR_CREATE_DATASET;
    }

    if (opts.compression_level > 0 &&
        H5Pset_deflate(dcpl, (unsigned)opts.compression_level) < 0) {
      H5Pclose(dcpl);
      H5Sclose(space);

      return HDF5_ERR_CREATE_DATASET;
    }
  }

  hid_t native_type = hdf5_native_type(type);
  hid_t dset = H5Dcreate2(w->file, dataset, native_type, space, H5P_DEFAULT,
                          dcpl, H5P_DEFAULT);

  if (owns_dcpl) {
    H5Pclose(dcpl);
  }

  if (dset < 0) {
    H5Sclose(space);

    return HDF5_ERR_CREATE_DATASET;
  }

  herr_t status =
      H5Dwrite(dset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dset);
  H5Sclose(space);

  return (status < 0) ? HDF5_ERR_WRITE : HDF5_OK;
}

hdf5_status_t hdf5_writer_write_f32_1d(hdf5_writer_t *writer,
                                       const char *dataset, const float *data,
                                       size_t len) {
  const size_t dims[1] = {len};

  return hdf5_write_dataset(writer, dataset, HDF5_TYPE_F32, 1, dims, data,
                            NULL);
}

hdf5_status_t hdf5_writer_write_i32_1d(hdf5_writer_t *writer,
                                       const char *dataset, const int32_t *data,
                                       size_t len) {
  const size_t dims[1] = {len};

  return hdf5_write_dataset(writer, dataset, HDF5_TYPE_I32, 1, dims, data,
                            NULL);
}

hdf5_status_t hdf5_writer_write_i64_1d(hdf5_writer_t *writer,
                                       const char *dataset, const int64_t *data,
                                       size_t len) {
  const size_t dims[1] = {len};

  return hdf5_write_dataset(writer, dataset, HDF5_TYPE_I64, 1, dims, data,
                            NULL);
}

hdf5_status_t hdf5_create_group(hdf5_writer_t *writer, const char *path) {
  if (!writer || !path || path[0] == '\0') {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  if (H5Lexists(writer->file, path, H5P_DEFAULT) > 0) {
    return HDF5_OK; /* already exists */
  }

  /* Reuse ensure_parent_groups to create `path` itself (plus any missing groups
   * above it), by probing as if `path` were parent of a one-level-deeper name -
   * ensure_parent_groups creates every path component strictly above  probe's
   * last '/', which is exactly path` itself */
  char probe[HDF5_MAX_PATH];

  int n = snprintf(probe, sizeof probe, "%s/.", path);
  if (n < 0 || (size_t)n >= sizeof probe) {
    return HDF5_ERR_PATH_TOO_LONG;
  }

  return ensure_parent_groups(writer->file, probe);
}

static hdf5_status_t
hdf5_write_attribute_common(hdf5_writer_t *writer, const char *object_path,
                            const char *attr_name, hid_t attr_type,
                            const void *value, size_t value_size) {
  hid_t obj = H5Oopen(writer->file, object_path, H5P_DEFAULT);
  if (obj < 0) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  hid_t space = H5Screate(H5S_SCALAR);
  if (space < 0) {
    H5Oclose(obj);

    return HDF5_ERR_WRITE;
  }

  if (H5Aexists(obj, attr_name) > 0) {
    H5Adelete(obj, attr_name);
  }

  hdf5_status_t rc = HDF5_OK;
  hid_t attr =
      H5Acreate2(obj, attr_name, attr_type, space, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0 || H5Awrite(attr, attr_type, value) < 0) {
    rc = HDF5_ERR_WRITE;
  }

  (void)value_size;

  if (attr >= 0) {
    H5Aclose(attr);
  }
  H5Sclose(space);
  H5Oclose(obj);

  return rc;
}

hdf5_status_t hdf5_write_attribute_string(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name,
                                          const char *value) {
  if (!writer || !object_path || !attr_name || !value) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  hid_t type = H5Tcopy(H5T_C_S1);
  if (type < 0 || H5Tset_size(type, strlen(value) + 1) < 0) {
    if (type >= 0) {
      H5Tclose(type);
    }

    return HDF5_ERR_WRITE;
  }

  hdf5_status_t rc = hdf5_write_attribute_common(writer, object_path, attr_name,
                                                 type, value, 0);

  H5Tclose(type);

  return rc;
}

hdf5_status_t hdf5_write_attribute_double(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name, double value) {
  if (!writer || !object_path || !attr_name) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  return hdf5_write_attribute_common(writer, object_path, attr_name,
                                     H5T_NATIVE_DOUBLE, &value, sizeof value);
}

static int legacy_write_dataset(const char *filename, const char *dataset,
                                const double *data, int rank,
                                const hsize_t *dims) {
  if (!filename || !dataset || !data) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  (void)mkdir(QMC_OUTPUT_DIR, 0755);

  char path[HDF5_MAX_PATH];
  int n = snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  if (n < 0 || (size_t)n >= sizeof path) {
    return HDF5_ERR_PATH_TOO_LONG;
  }

  /* Open an existing file if present (so prior datasets survive), else
   * create a new one. Subsequent dataset write replaces any dataset
   * of same name, but does not touch other datasets */
  hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0) {
    file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  if (file < 0) {
    return HDF5_ERR_OPEN;
  }

  hdf5_writer_t stack_writer = {.file = file};
  hdf5_status_t rc =
      write_dataset_impl(&stack_writer, dataset, data, rank, dims);

  if (H5Fclose(file) < 0 && rc == HDF5_OK) {
    rc = HDF5_ERR_CLOSE;
  }

  return (int)rc;
}

int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len) {
  const hsize_t dims[1] = {(hsize_t)len};

  return legacy_write_dataset(filename, dataset, data, 1, dims);
}

int hdf5_write_matrix(const char *filename, const char *dataset,
                      const double *data, size_t rows, size_t cols) {
  const hsize_t dims[2] = {(hsize_t)rows, (hsize_t)cols};

  return legacy_write_dataset(filename, dataset, data, 2, dims);
}

#else /* !USE_HDF5 */

struct hdf5_writer {
  // cppcheck-suppress unusedStructMember
  int unused;
};

static const char *kNoHdf5Msg =
    "hdf5_writer: not built with HDF5 support "
    "(rebuild with 'make USE_HDF5=1', which requires libhdf5-dev)\n";

hdf5_writer_t *hdf5_open(const char *path) {
  (void)path;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return NULL;
}

hdf5_status_t hdf5_close(hdf5_writer_t *writer) {
  (void)writer;

  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_writer_write_1d(hdf5_writer_t *writer, const char *dataset,
                                   const double *data, size_t len) {
  (void)writer;
  (void)dataset;
  (void)data;
  (void)len;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_writer_write_matrix(hdf5_writer_t *writer,
                                       const char *dataset, const double *data,
                                       size_t rows, size_t cols) {
  (void)writer;
  (void)dataset;
  (void)data;
  (void)rows;
  (void)cols;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len) {
  (void)filename;
  (void)dataset;
  (void)data;
  (void)len;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

int hdf5_write_matrix(const char *filename, const char *dataset,
                      const double *data, size_t rows, size_t cols) {
  (void)filename;
  (void)dataset;
  (void)data;
  (void)rows;
  (void)cols;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_write_dataset(hdf5_writer_t *writer, const char *dataset,
                                 hdf5_type_t type, int rank,
                                 const size_t *dims, const void *data,
                                 const hdf5_dataset_options_t *options) {
  (void)writer;
  (void)dataset;
  (void)type;
  (void)rank;
  (void)dims;
  (void)data;
  (void)options;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_writer_write_f32_1d(hdf5_writer_t *writer,
                                       const char *dataset, const float *data,
                                       size_t len) {
  (void)writer;
  (void)dataset;
  (void)data;
  (void)len;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_writer_write_i32_1d(hdf5_writer_t *writer,
                                       const char *dataset,
                                       const int32_t *data, size_t len) {
  (void)writer;
  (void)dataset;
  (void)data;
  (void)len;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_writer_write_i64_1d(hdf5_writer_t *writer,
                                       const char *dataset,
                                       const int64_t *data, size_t len) {
  (void)writer;
  (void)dataset;
  (void)data;
  (void)len;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_create_group(hdf5_writer_t *writer, const char *path) {
  (void)writer;
  (void)path;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_write_attribute_string(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name,
                                          const char *value) {
  (void)writer;
  (void)object_path;
  (void)attr_name;
  (void)value;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

hdf5_status_t hdf5_write_attribute_double(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name,
                                          double value) {
  (void)writer;
  (void)object_path;
  (void)attr_name;
  (void)value;

  fprintf(stderr, "%s", kNoHdf5Msg);
  return HDF5_ERR_NOT_SUPPORTED;
}

#endif /* USE_HDF5 */
