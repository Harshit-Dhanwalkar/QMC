# `qmc_version`

> Generated from `core/qmc_version.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Library version macros (compile time) and qmc_version()/qmc_version_number() (run time).

**Category:** Core: data types & utilities

```c
#include <qmc/core/qmc_version.h>
```

**Start with:** [`qmc_version`](#qmc_version), [`qmc_version_number`](#qmc_version_number)

- Runtime help: `qmc_help(stdout, "qmc_version");` or `qmc-help qmc_version`

## Overview

```text
QMC library version: single source of truth

Makefile reads QMC_VERSION_STRING from this file for qmc.pc and
shared-library name (libqmc.so.<MAJOR>), and export/export_schema.h derives
QMC_LIBRARY_VERSION from it, so a release is a one-place edit: change
four numbers/strings below together (tests/test_version.c checks that they
agree)
```

## Macros and constants

```c
#define QMC_VERSION_MAJOR 0
```

```c
#define QMC_VERSION_MINOR 1
```

```c
#define QMC_VERSION_PATCH 0
```

```c
#define QMC_VERSION_STRING "0.1.0"
```

```c
#define QMC_VERSION_NUMBER                                                     \
  (QMC_VERSION_MAJOR * 10000 + QMC_VERSION_MINOR * 100 + QMC_VERSION_PATCH)
```

```text
MAJOR * 10000 + MINOR * 100 + PATCH
```

## Functions (2)

[`qmc_version`](#qmc_version), [`qmc_version_number`](#qmc_version_number)

### `qmc_version`

```c
const char *qmc_version(void);
```

```text
Version of the compiled library (not of the headers): "MAJOR.MINOR.PATCH"
```

### `qmc_version_number`

```c
int qmc_version_number(void);
```

```text
Same as QMC_VERSION_NUMBER, but of the compiled library
```
