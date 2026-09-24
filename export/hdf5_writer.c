#include "hdf5_writer.h"

#include <stdio.h>
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

#ifdef USE_HDF5

#include <hdf5.h>

struct hdf5_writer {
  hid_t file;
};

/* Write (or replace) a single dataset in an already-open file. */
static hdf5_status_t write_dataset_impl(hdf5_writer_t *w, const char *dataset,
                                        const double *data, int rank,
                                        const hsize_t *dims) {
  if (!w || !dataset || !data) {
    return HDF5_ERR_INVALID_ARGUMENT;
  }

  /* Replace if a dataset of this name already exists. */
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
   * create a new one. The subsequent dataset write replaces any dataset
   * of the same name, but does not touch other datasets */
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

#endif /* USE_HDF5 */
