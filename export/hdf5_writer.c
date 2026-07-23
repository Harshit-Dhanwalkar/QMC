// TODO: Implement data export to hdf5
#include <stdio.h>
#include <stdlib.h>

// #include <hdf5.h>
// int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
//                   size_t len) {
//   hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//   if (file < 0)
//     return -1;
//   hsize_t dims[1] = {len};
//   hid_t space = H5Screate_simple(1, dims, NULL);
//   hid_t dset = H5Dcreate(file, dataset, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT,
//                          H5P_DEFAULT, H5P_DEFAULT);
//   H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
//   H5Dclose(dset);
//   H5Sclose(space);
//   H5Fclose(file);
//   return 0;

int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len) {
  fprintf(stderr, "HDF5 not is not supported yet.\n");
  return -1;
}
