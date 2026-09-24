#include "netcdf_writer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

#define NETCDF_MAX_PATH 4096
#define NETCDF_DIM_NAME_MAX 256

const char *netcdf_strerror(netcdf_status_t status) {
  switch (status) {
  case NETCDF_OK:
    return "ok";
  case NETCDF_ERR_INVALID_ARGUMENT:
    return "invalid argument";
  case NETCDF_ERR_OPEN:
    return "could not open/create file";
  case NETCDF_ERR_CREATE_VARIABLE:
    return "could not create variable";
  case NETCDF_ERR_WRITE:
    return "variable write failed";
  case NETCDF_ERR_CLOSE:
    return "file close failed";
  case NETCDF_ERR_NOT_SUPPORTED:
    return "not built with NetCDF support";
  case NETCDF_ERR_PATH_TOO_LONG:
    return "path too long";
  default:
    return "unknown error";
  }
}

#ifdef USE_NETCDF

#include <netcdf.h>

/* NETCDF4 (HDF5-backed) is required so that a variable can be deleted and
 * re-created when its name is reused. Classic NetCDF3 files do not support
 * nc_del_var / nc_del_dim, which would make the "replace an existing variable"
 * semantics below impossible */
#define NETCDF_FILE_FORMAT (NC_NETCDF4)

struct netcdf_writer {
  int ncid;
};

// HACK: "replace a variable with a new shape, delete old dims" was intention of
// mine to implement but due to underlying NetCDF C API only grew tool for it
// (nc_del_var, nc_del_dim) in version 4.9.0. So version of my intent cannot
// link on the CI runner either
// /*
//  * Define (or redefine) a variable and write its data
//  *
//  * Each variable owns a private set of dimensions named "<varname>_dim0",
//  * "<varname>_dim1", ... so that two variables with same shape do not collide
//  * on dimension names, and so that replacing a variable is a well-defined
//  * delete-then-redefine operation
//  *
//  * Caller must not have the file in define mode on entry
//  */
// static netcdf_status_t define_and_write_variable(netcdf_writer_t *w,
//                                                  const char *varname,
//                                                  const double *data, int
//                                                  ndims, const size_t *shape)
//                                                  {
//   int rc = nc_redef(w->ncid);
//   if (rc != NC_NOERR && rc != NC_EINDEFINE) {
//     return NETCDF_ERR_CREATE_VARIABLE;
//   }
//
//   /* If a variable of this name already exists, remove it
//    * Delete variable first (that drops its reference to the dimensions), then
//    * try to remove the dimensions we previously defined for it. */
//   int old_varid;
//   if (nc_inq_varid(w->ncid, varname, &old_varid) == NC_NOERR) {
//     if (nc_del_var(w->ncid, old_varid) != NC_NOERR) {
//       nc_enddef(w->ncid);
//
//       return NETCDF_ERR_CREATE_VARIABLE;
//     }
//
//     char dimname[NETCDF_DIM_NAME_MAX];
//     for (int i = 0; i < ndims; i++) {
//       snprintf(dimname, sizeof dimname, "%s_dim%d", varname, i);
//       int dimid;
//       if (nc_inq_dimid(w->ncid, dimname, &dimid) == NC_NOERR) {
//         // best-effort; will fail if some other variable references it
//         (void)nc_del_dim(w->ncid, dimid);
//       }
//     }
//   }
//
//   int dimids[2];
//   char dimname[NETCDF_DIM_NAME_MAX];
//   for (int i = 0; i < ndims; i++) {
//     snprintf(dimname, sizeof dimname, "%s_dim%d", varname, i);
//     if (nc_def_dim(w->ncid, dimname, shape[i], &dimids[i]) != NC_NOERR) {
//       nc_enddef(w->ncid);
//
//       return NETCDF_ERR_CREATE_VARIABLE;
//     }
//   }
//
//   int varid;
//   if (nc_def_var(w->ncid, varname, NC_DOUBLE, ndims, dimids, &varid) !=
//       NC_NOERR) {
//     nc_enddef(w->ncid);
//
//     return NETCDF_ERR_CREATE_VARIABLE;
//   }
//
//   if (nc_enddef(w->ncid) != NC_NOERR) {
//     return NETCDF_ERR_CREATE_VARIABLE;
//   }
//
//   if (nc_put_var_double(w->ncid, varid, data) != NC_NOERR) {
//     return NETCDF_ERR_WRITE;
//   }
//
//   return NETCDF_OK;
// }

/*
 * Define (or redefine) a variable and write its data
 *
 * Each variable owns a private set of dimensions named "<varname>_dim0",
 * "<varname>_dim1", ... so that two variables with same shape do not collide on
 * dimension names, and so that replacing a variable is a well-defined
 * delete-then-redefine operation
 *
 * Caller must not have file in define mode on entry
 */
/*
 * Define (or update) a variable and write its data
 *
 * Each variable owns a private set of dimensions named "<varname>_dim0",
 * "<varname>_dim1", ... so two variables with same shape do not
 * collide on dimension names
 *
 * NOTE:: Strategy (works on all NetCDF 4.x; nc_del_var/nc_del_dim pair is
 * only available from NetCDF 4.9.0 onward, so it cannot be used here):
 *   - Variable does not exist yet:
 *       define it (dimensions + variable), enddef, write
 *   - Variable exists, shape matches what caller requested:
 *       write into it in place. No redefinition needed
 *   - Variable exists, shape differs:
 *       return NETCDF_ERR_CREATE_VARIABLE. Changing a variable's shape
 *       would require deleting and recreating it, which pre-4.9 NetCDF
 *       does not support
 *
 * Caller must not have file in define mode on entry
 */
static netcdf_status_t define_and_write_variable(netcdf_writer_t *w,
                                                 const char *varname,
                                                 const double *data, int ndims,
                                                 const size_t *shape) {
  int varid;

  // Case 1: variable already exists
  if (nc_inq_varid(w->ncid, varname, &varid) == NC_NOERR) {
    /* Shape must match to allow in-place overwrite. */
    int vndims = 0;
    if (nc_inq_varndims(w->ncid, varid, &vndims) != NC_NOERR ||
        vndims != ndims) {
      return NETCDF_ERR_CREATE_VARIABLE;
    }

    int dimids[2];
    if (nc_inq_vardimid(w->ncid, varid, dimids) != NC_NOERR) {
      return NETCDF_ERR_CREATE_VARIABLE;
    }

    for (int i = 0; i < ndims; i++) {
      size_t len = 0;
      if (nc_inq_dimlen(w->ncid, dimids[i], &len) != NC_NOERR ||
          len != shape[i]) {
        return NETCDF_ERR_CREATE_VARIABLE;
      }
    }

    // Shape matches: overwrite data in place
    if (nc_put_var_double(w->ncid, varid, data) != NC_NOERR) {
      return NETCDF_ERR_WRITE;
    }

    return NETCDF_OK;
  }

  // Case 2: variable does not exist yet: define it
  int rc = nc_redef(w->ncid);
  if (rc != NC_NOERR && rc != NC_EINDEFINE) {
    return NETCDF_ERR_CREATE_VARIABLE;
  }

  int dimids[2];
  char dimname[NETCDF_DIM_NAME_MAX];
  for (int i = 0; i < ndims; i++) {
    snprintf(dimname, sizeof dimname, "%s_dim%d", varname, i);
    if (nc_def_dim(w->ncid, dimname, shape[i], &dimids[i]) != NC_NOERR) {
      nc_enddef(w->ncid);

      return NETCDF_ERR_CREATE_VARIABLE;
    }
  }

  if (nc_def_var(w->ncid, varname, NC_DOUBLE, ndims, dimids, &varid) !=
      NC_NOERR) {
    nc_enddef(w->ncid);

    return NETCDF_ERR_CREATE_VARIABLE;
  }

  if (nc_enddef(w->ncid) != NC_NOERR) {
    return NETCDF_ERR_CREATE_VARIABLE;
  }

  if (nc_put_var_double(w->ncid, varid, data) != NC_NOERR) {
    return NETCDF_ERR_WRITE;
  }

  return NETCDF_OK;
}

netcdf_writer_t *netcdf_open(const char *path) {
  if (!path || path[0] == '\0') {
    return NULL;
  }

  int ncid;
  if (nc_create(path, NC_CLOBBER | NETCDF_FILE_FORMAT, &ncid) != NC_NOERR) {
    return NULL;
  }

  netcdf_writer_t *w = malloc(sizeof *w);
  if (!w) {
    nc_close(ncid);

    return NULL;
  }

  w->ncid = ncid;

  return w;
}

netcdf_status_t netcdf_close(netcdf_writer_t *writer) {
  if (!writer) {
    return NETCDF_ERR_INVALID_ARGUMENT;
  }

  netcdf_status_t rc = NETCDF_OK;
  if (writer->ncid >= 0) {
    if (nc_close(writer->ncid) != NC_NOERR) {
      rc = NETCDF_ERR_CLOSE;
    }
  }

  free(writer);
  return rc;
}

netcdf_status_t netcdf_writer_write_1d(netcdf_writer_t *writer,
                                       const char *varname, const double *data,
                                       size_t len) {
  if (!writer || !varname || (!data && len > 0)) {
    return NETCDF_ERR_INVALID_ARGUMENT;
  }

  const size_t shape[1] = {len};

  return define_and_write_variable(writer, varname, data, 1, shape);
}

netcdf_status_t netcdf_writer_write_matrix(netcdf_writer_t *writer,
                                           const char *varname,
                                           const double *data, size_t rows,
                                           size_t cols) {
  if (!writer || !varname || (!data && rows > 0 && cols > 0)) {
    return NETCDF_ERR_INVALID_ARGUMENT;
  }

  const size_t shape[2] = {rows, cols};

  return define_and_write_variable(writer, varname, data, 2, shape);
}

static netcdf_status_t legacy_write_variable(const char *filename,
                                             const char *varname,
                                             const double *data, int ndims,
                                             const size_t *shape) {
  if (!filename || !varname || !data) {
    return NETCDF_ERR_INVALID_ARGUMENT;
  }

  (void)mkdir(QMC_OUTPUT_DIR, 0755);

  char path[NETCDF_MAX_PATH];
  int n = snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  if (n < 0 || (size_t)n >= sizeof path) {
    return NETCDF_ERR_PATH_TOO_LONG;
  }

  // Open existing file if present else create a new one
  int ncid;
  if (nc_open(path, NC_WRITE, &ncid) != NC_NOERR) {
    if (nc_create(path, NC_CLOBBER | NETCDF_FILE_FORMAT, &ncid) != NC_NOERR) {
      return NETCDF_ERR_OPEN;
    }
  }

  netcdf_writer_t stack_writer = {.ncid = ncid};
  netcdf_status_t rc =
      define_and_write_variable(&stack_writer, varname, data, ndims, shape);

  if (nc_close(ncid) != NC_NOERR && rc == NETCDF_OK) {
    rc = NETCDF_ERR_CLOSE;
  }

  return rc;
}

int netcdf_write_1d(const char *filename, const char *varname,
                    const double *data, size_t len) {
  const size_t shape[1] = {len};

  return (int)legacy_write_variable(filename, varname, data, 1, shape);
}

int netcdf_write_matrix(const char *filename, const char *varname,
                        const double *data, size_t rows, size_t cols) {
  const size_t shape[2] = {rows, cols};

  return (int)legacy_write_variable(filename, varname, data, 2, shape);
}

#else /* !USE_NETCDF */

struct netcdf_writer {
  int unused;
};

static const char *kNoNetcdfMsg =
    "netcdf_writer: not built with NetCDF support "
    "(rebuild with 'make USE_NETCDF=1', which requires libnetcdf-dev)\n";

netcdf_writer_t *netcdf_open(const char *path) {
  (void)path;

  fprintf(stderr, "%s", kNoNetcdfMsg);
  return NULL;
}

netcdf_status_t netcdf_close(netcdf_writer_t *writer) {
  (void)writer;

  return NETCDF_ERR_NOT_SUPPORTED;
}

netcdf_status_t netcdf_writer_write_1d(netcdf_writer_t *writer,
                                       const char *varname, const double *data,
                                       size_t len) {
  (void)writer;
  (void)varname;
  (void)data;
  (void)len;

  fprintf(stderr, "%s", kNoNetcdfMsg);
  return NETCDF_ERR_NOT_SUPPORTED;
}

netcdf_status_t netcdf_writer_write_matrix(netcdf_writer_t *writer,
                                           const char *varname,
                                           const double *data, size_t rows,
                                           size_t cols) {
  (void)writer;
  (void)varname;
  (void)data;
  (void)rows;
  (void)cols;

  fprintf(stderr, "%s", kNoNetcdfMsg);
  return NETCDF_ERR_NOT_SUPPORTED;
}

int netcdf_write_1d(const char *filename, const char *varname,
                    const double *data, size_t len) {
  (void)filename;
  (void)varname;
  (void)data;
  (void)len;

  fprintf(stderr, "%s", kNoNetcdfMsg);
  return NETCDF_ERR_NOT_SUPPORTED;
}

int netcdf_write_matrix(const char *filename, const char *varname,
                        const double *data, size_t rows, size_t cols) {
  (void)filename;
  (void)varname;
  (void)data;
  (void)rows;
  (void)cols;

  fprintf(stderr, "%s", kNoNetcdfMsg);
  return NETCDF_ERR_NOT_SUPPORTED;
}

#endif /* USE_NETCDF */
