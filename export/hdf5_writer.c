#include "hdf5_writer.h"

#include <stdio.h>
#include <stdlib.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

#ifdef USE_HDF5

#include <hdf5.h>

static int write_dataset(const char *filename, const char *dataset,
                         const double *data, int rank, const hsize_t *dims) {
  if (!filename || !dataset || !data) {
    fprintf(stderr, "hdf5_writer: NULL filename/dataset/data\n");
    return -1;
  }

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);

  hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr, "hdf5_writer: could not create '%s'\n", path);
    return -1;
  }

  hid_t space = H5Screate_simple(rank, dims, NULL);
  if (space < 0) {
    H5Fclose(file);

    return -1;
  }

  hid_t dset = H5Dcreate2(file, dataset, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) {
    H5Sclose(space);
    H5Fclose(file);

    return -1;
  }

  herr_t status =
      H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dset);
  H5Sclose(space);
  H5Fclose(file);

  return (status < 0) ? -1 : 0;
}

int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len) {
  const hsize_t dims[1] = {(hsize_t)len};

  return write_dataset(filename, dataset, data, 1, dims);
}

int hdf5_write_matrix(const char *filename, const char *dataset,
                      const double *data, size_t rows, size_t cols) {
  const hsize_t dims[2] = {(hsize_t)rows, (hsize_t)cols};

  return write_dataset(filename, dataset, data, 2, dims);
}

#else // !USE_HDF5

// Built without USE_HDF5
int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len) {
  (void)filename;
  (void)dataset;
  (void)data;
  (void)len;
  fprintf(stderr, "hdf5_writer: not built with HDF5 support (rebuild with "
                  "'make USE_HDF5=1', which requires libhdf5-dev)\n");

  return -1;
}

int hdf5_write_matrix(const char *filename, const char *dataset,
                      const double *data, size_t rows, size_t cols) {
  (void)filename;
  (void)dataset;
  (void)data;
  (void)rows;
  (void)cols;
  fprintf(stderr, "hdf5_writer: not built with HDF5 support (rebuild with "
                  "'make USE_HDF5=1', which requires libhdf5-dev)\n");

  return -1;
}

#endif // USE_HDF5
