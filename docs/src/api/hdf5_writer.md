# `hdf5_writer`

> Generated from `export/hdf5_writer.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

HDF5 export: typed datasets, groups, attributes, chunking and compression, appendable datasets (needs USE_HDF5=1).

**Category:** Export, plotting & LaTeX

```c
#include <qmc/export/hdf5_writer.h>   /* not part of <qmc/qmc.h>: include it directly */
```

**Start with:** [`hdf5_open`](#hdf5_open), [`hdf5_write_dataset`](#hdf5_write_dataset), [`hdf5_close`](#hdf5_close)

- Test (every call with expected results): `tests/test_hdf5_writer.c`
- Runtime help: `qmc_help(stdout, "hdf5_writer");` or `qmc-help hdf5_writer`

## Types

### `hdf5_status_t`

```c
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
```

### `hdf5_type_t`

```text
Element type for generic hdf5_write_dataset() entry point
```

```c
typedef enum {
  HDF5_TYPE_F32,
  HDF5_TYPE_F64,
  HDF5_TYPE_I32,
  HDF5_TYPE_I64
} hdf5_type_t;
```

### `hdf5_dataset_options_t`

```text
Per-dataset chunking/compression knobs

chunk_size == 0 means "no explicit chunking". If compression_level > 0 or
shuffle is requested with chunk_size == 0, an implementation-chosen default
chunk extent (1024 elements per axis, clamped to dataset dimension) is used
instead - HDF5 requires a chunked layout for filters. Otherwise dataset is
written contiguously and compression/shuffle are ignored
When chunk_size > 0, it is used as target chunk extent along every dimension,
clamped to that dimension's actual size
```

```c
typedef struct {
  int compression_level; /* 0 = none, 1-9 = gzip (deflate) level */
  size_t chunk_size;     /* elements per chunk dimension; 0 = no chunking */
  int shuffle;           /* 1 = enable shuffle filter before deflate */
} hdf5_dataset_options_t;
```

### `hdf5_writer_t`

```text
Persistent HDF5 writer

hdf5_open() creates file once (truncating if it already exists) and
subsequent hdf5_writer_write_*() calls add datasets to that same file without
touching previously-written data
```

```c
typedef struct hdf5_writer hdf5_writer_t;
```

### `hdf5_dataset_t`

```text
Appendable ("extendible") dataset: grows along its first axis as rows arrive,
so a long run can stream observables instead of buffering them all

  hdf5_dataset_t *e = hdf5_dataset_create(w, "observables/energy",
                                          HDF5_TYPE_F64, 1, NULL, NULL);
  for (...) hdf5_dataset_append(e, &energy_i, 1);
  hdf5_dataset_close(e);   // close every dataset BEFORE hdf5_close(w)

row_width == 1 gives a 1D dataset of shape (rows,); row_width > 1 gives a 2D
dataset of shape (rows, row_width) and each appended row supplies row_width
consecutive elements. Appendable datasets are always chunked: for them
options->chunk_size means *rows per chunk* (default 1024 when 0/NULL);
compression_level and shuffle apply as usual. An existing dataset of same
name is replaced. On failure returns NULL and, if status_out != NULL, stores
reason there
```

```c
typedef struct hdf5_dataset hdf5_dataset_t;
```

## Functions (23)

[`hdf5_dataset_options_default`](#hdf5_dataset_options_default), [`hdf5_open`](#hdf5_open), [`hdf5_close`](#hdf5_close), [`hdf5_writer_write_1d`](#hdf5_writer_write_1d), [`hdf5_writer_write_matrix`](#hdf5_writer_write_matrix), [`hdf5_write_dataset`](#hdf5_write_dataset), [`hdf5_writer_write_f32_1d`](#hdf5_writer_write_f32_1d), [`hdf5_writer_write_i32_1d`](#hdf5_writer_write_i32_1d), [`hdf5_writer_write_i64_1d`](#hdf5_writer_write_i64_1d), [`hdf5_create_group`](#hdf5_create_group), [`hdf5_write_attribute_string`](#hdf5_write_attribute_string), [`hdf5_write_attribute_double`](#hdf5_write_attribute_double), [`hdf5_write_attribute_int`](#hdf5_write_attribute_int), [`hdf5_write_string`](#hdf5_write_string), [`hdf5_write_schema_header`](#hdf5_write_schema_header), [`hdf5_flush`](#hdf5_flush), [`hdf5_dataset_create`](#hdf5_dataset_create), [`hdf5_dataset_append`](#hdf5_dataset_append), [`hdf5_dataset_rows`](#hdf5_dataset_rows), [`hdf5_dataset_close`](#hdf5_dataset_close), [`hdf5_write_1d`](#hdf5_write_1d), [`hdf5_write_matrix`](#hdf5_write_matrix), [`hdf5_strerror`](#hdf5_strerror)

### `hdf5_dataset_options_default`

```c
hdf5_dataset_options_t hdf5_dataset_options_default(void);
```

```text
Returns default options: no chunking, no compression, no shuffle
```

### `hdf5_open`

```c
hdf5_writer_t *hdf5_open(const char *path);
```

```text
Open (create, or truncate-and-recreate) `path` for writing

Returns NULL on failure
```

### `hdf5_close`

```c
hdf5_status_t hdf5_close(hdf5_writer_t *writer);
```

```text
Flush and close

Returns HDF5_OK only if file closed cleanly
```

### `hdf5_writer_write_1d`

```c
hdf5_status_t hdf5_writer_write_1d(hdf5_writer_t *writer, const char *dataset,
                                   const double *data, size_t len);
```

```text
Write native-double datasets into an open writer

Dataset whose name already exists in file is replaced; all odatasets
in file are left untouched
```

### `hdf5_writer_write_matrix`

```c
hdf5_status_t hdf5_writer_write_matrix(hdf5_writer_t *writer,
                                       const char *dataset, const double *data,
                                       size_t rows, size_t cols);
```

_No description in header yet._

### `hdf5_write_dataset`

```c
hdf5_status_t hdf5_write_dataset(hdf5_writer_t *writer, const char *dataset,
                                 hdf5_type_t type, int rank, const size_t *dims,
                                 const void *data,
                                 const hdf5_dataset_options_t *options);
```

```text
Generic dataset writer: any rank, any of hdf5_type_t element types, optional
chunking/compression via `options` (NULL -> defaults, i.e. no
chunking/compression, matching hdf5_writer_write_1d/matrix)

`dataset` may contain '/' to place dataset inside a group path
("simulation/energy"); any missing intermediate groups are created
automatically for that prefix
```

### `hdf5_writer_write_f32_1d`

```c
hdf5_status_t hdf5_writer_write_f32_1d(hdf5_writer_t *writer,
                                       const char *dataset, const float *data,
                                       size_t len);
```

```text
Typed convenience wrappers around hdf5_write_dataset(), 1D case, default
(no chunking/compression) options
```

### `hdf5_writer_write_i32_1d`

```c
hdf5_status_t hdf5_writer_write_i32_1d(hdf5_writer_t *writer,
                                       const char *dataset, const int32_t *data,
                                       size_t len);
```

_No description in header yet._

### `hdf5_writer_write_i64_1d`

```c
hdf5_status_t hdf5_writer_write_i64_1d(hdf5_writer_t *writer,
                                       const char *dataset, const int64_t *data,
                                       size_t len);
```

_No description in header yet._

### `hdf5_create_group`

```c
hdf5_status_t hdf5_create_group(hdf5_writer_t *writer, const char *path);
```

```text
Create a group (and any missing intermediate groups) at `path`, e.g.
"simulation/observables". A no-op (returns HDF5_OK) if it already exists
```

### `hdf5_write_attribute_string`

```c
hdf5_status_t hdf5_write_attribute_string(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name,
                                          const char *value);
```

```text
Attach a string/double attribute to a dataset or group. Use "/" as
object_path for file-level (root group) metadata, e.g. schema/library
version. Replaces an existing attribute of same name on that object
```

### `hdf5_write_attribute_double`

```c
hdf5_status_t hdf5_write_attribute_double(hdf5_writer_t *writer,
                                          const char *object_path,
                                          const char *attr_name, double value);
```

_No description in header yet._

### `hdf5_write_attribute_int`

```c
hdf5_status_t hdf5_write_attribute_int(hdf5_writer_t *writer,
                                       const char *object_path,
                                       const char *attr_name, long value);
```

```text
Integer attribute (same object_path rules as string/double variants)
```

### `hdf5_write_string`

```c
hdf5_status_t hdf5_write_string(hdf5_writer_t *writer, const char *dataset,
                                const char *value);
```

```text
Scalar (single-value) string dataset, e.g. hdf5_write_string(w,
"metadata/method", "DMC"). Replaces an existing dataset of same name;
 missing intermediate groups are created
```

### `hdf5_write_schema_header`

```c
hdf5_status_t hdf5_write_schema_header(hdf5_writer_t *writer,
                                       int include_timestamp);
```

```text
Stamps file-level provenance as root attributes, mirroring JSON header:
  schema_version (int), library, library_version, [created (UTC ISO 8601)]
 `schema_version` is *file layout* version and is independent of library
  version. Pass include_timestamp = 0 for reproducible files
```

### `hdf5_flush`

```c
hdf5_status_t hdf5_flush(hdf5_writer_t *writer);
```

```text
Push buffered data to disk without closing - worth calling periodically in
long simulations so a crash keeps everything written so far
```

### `hdf5_dataset_create`

```c
hdf5_dataset_t *hdf5_dataset_create(hdf5_writer_t *writer, const char *dataset,
                                    hdf5_type_t type, size_t row_width,
                                    const hdf5_dataset_options_t *options,
                                    hdf5_status_t *status_out);
```

_No description in header yet._

### `hdf5_dataset_append`

```c
hdf5_status_t hdf5_dataset_append(hdf5_dataset_t *ds, const void *data,
                                  size_t nrows);
```

```text
Appends `nrows` rows (nrows * row_width elements of dataset's type).
nrows == 0 is a valid no-op (data may then be NULL)
```

### `hdf5_dataset_rows`

```c
size_t hdf5_dataset_rows(const hdf5_dataset_t *ds);
```

```text
Rows written so far (0 for a NULL handle)
```

### `hdf5_dataset_close`

```c
hdf5_status_t hdf5_dataset_close(hdf5_dataset_t *ds);
```

```text
Release handle; data stays in file
```

### `hdf5_write_1d`

```c
int hdf5_write_1d(const char *filename, const char *dataset, const double *data,
                  size_t len);
```

```text
Writes into QMC_OUTPUT_DIR/filename. On first call, creates file; on
subsequent calls to same file, opens it in read-write mode and
adds/replaces only named dataset. This means:

   hdf5_write_1d("results.h5", "energy",   e, n);   // creates file
   hdf5_write_1d("results.h5", "variance", v, n);   // adds variance

leaves a file containing both "energy" and "variance"

Returns HDF5_OK (0) on success, a negative hdf5_status_t on failure
```

### `hdf5_write_matrix`

```c
int hdf5_write_matrix(const char *filename, const char *dataset,
                      const double *data, size_t rows, size_t cols);
```

_No description in header yet._

### `hdf5_strerror`

```c
const char *hdf5_strerror(hdf5_status_t status);
```

```text
Description of a status code
```
