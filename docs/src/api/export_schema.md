# `export_schema`

> Generated from `export/export_schema.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Conventions shared by all exporters: schema version, library name and version, UTC timestamps, column descriptors.

**Category:** Export, plotting & LaTeX

```c
#include <qmc/export/export_schema.h>
```

**Start with:** [`export_timestamp_utc`](#export_timestamp_utc)

- Test (every call with expected results): `tests/test_json_writer.c`
- Test (every call with expected results): `tests/test_version.c`
- Runtime help: `qmc_help(stdout, "export_schema");` or `qmc-help export_schema`

## Overview

```text
Conventions shared by every exporter (JSON, HDF5, CSV sidecars)

schema_version describes layout of QMC's exported files and only changes on a
breaking layout change
```

## Types

### `export_column_t`

```text
Description of one tabular column. Shared by CSV writer (header row) and
JSON / HDF5 writers (unit + description metadata), so same descriptor array
can drive all three outputs. `unit` and `description` may be NULL
```

```c
typedef struct {
  const char *name;
  const char *unit;
  const char *description;
} export_column_t;
```

## Macros and constants

```c
#define QMC_EXPORT_SCHEMA_VERSION 1
```

```c
#define QMC_LIBRARY_NAME "QMC"
```

```c
#define QMC_LIBRARY_VERSION QMC_VERSION_STRING
```

```c
#define EXPORT_TIMESTAMP_SIZE 32
```

```text
Buffer size that always fits an ISO 8601 UTC timestamp + NUL
```

## Functions (1)

[`export_timestamp_utc`](#export_timestamp_utc)

### `export_timestamp_utc`

```c
int export_timestamp_utc(char *buf, size_t size);
```

```text
Writes current time as "YYYY-MM-DDTHH:MM:SSZ" (UTC) into buf

Returns 0 on success, -1 if buf is NULL/too small or clock fails
```
