# `catalog`

> Generated from `core/catalog.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Module catalog and help: discover modules and print where and how to use them (qmc_help).

**Category:** Core: data types & utilities

```c
#include <qmc/core/catalog.h>
```

**Start with:** [`qmc_module_at`](#qmc_module_at), [`qmc_module_find`](#qmc_module_find), [`qmc_module_search`](#qmc_module_search), [`qmc_category_at`](#qmc_category_at), [`qmc_help`](#qmc_help)

- Test (every call with expected results): `tests/test_catalog.c`
- Runtime help: `qmc_help(stdout, "catalog");` or `qmc-help catalog`

## Overview

```text
Module catalog and built-in help.

QMC has over a hundred modules. This catalog answers "what exists, what does
it do, where do I start, and where is a working example?" from inside your
program or from the shell, without opening the documentation:

  qmc_help(stdout, NULL);             overview of every category and module
  qmc_help(stdout, "tight_binding");  how to use one module
  qmc_help(stdout, "chemistry");      every module in a category
  qmc_help(stdout, "chern");          keyword search

The same data is available as a table (qmc_module_at, qmc_module_find,
qmc_module_search). The table is generated from the public headers by
tools/gen_api.py, so it cannot list a header that does not exist, and CI fails
when it is out of date.
```

## Types

### `qmc_module_t`

```c
typedef struct {
  const char *name;          /* module name, e.g. "tight_binding" */
  const char *header;        /* path under include/qmc/, e.g. "physics/tight_binding.h" */
  const char *category;      /* e.g. "Condensed matter & many-body" */
  const char *summary;       /* one line: what the module does */
  const char *key_functions; /* comma-separated entry points to start with */
  const char *guide;         /* hand-written guide page under docs/src/, or "" */
  const char *example;       /* example program in the repository, or "" */
  const char *test;          /* test showing every call with expected results, or "" */
  int in_qmc_h;              /* 1 if <qmc/qmc.h> includes it, 0 if you must include it directly */
} qmc_module_t;
```

## Macros and constants

```c
#define QMC_REPO_URL "https://github.com/Harshit-Dhanwalkar/QMC"
```

```c
#define QMC_DOCS_URL "https://harshit-dhanwalkar.github.io/QMC"
```

## Functions (7)

[`qmc_module_count`](#qmc_module_count), [`qmc_module_at`](#qmc_module_at), [`qmc_module_find`](#qmc_module_find), [`qmc_module_search`](#qmc_module_search), [`qmc_category_count`](#qmc_category_count), [`qmc_category_at`](#qmc_category_at), [`qmc_help`](#qmc_help)

### `qmc_module_count`

```c
size_t qmc_module_count(void);
```

```text
Number of modules in the catalog.
```

### `qmc_module_at`

```c
const qmc_module_t *qmc_module_at(size_t index);
```

```text
The index-th module (0 <= index < qmc_module_count()), grouped by category;
NULL if out of range.
```

### `qmc_module_find`

```c
const qmc_module_t *qmc_module_find(const char *name);
```

```text
Looks a module up by name. Case-insensitive; also accepts the header path
("physics/tight_binding.h", "qmc/physics/tight_binding.h" or
"tight_binding.h"). Returns NULL if there is no such module.
```

### `qmc_module_search`

```c
size_t qmc_module_search(const char *query, const qmc_module_t **out, size_t max_out);
```

```text
Case-insensitive substring search over module name, summary, category and
key functions. Writes up to max_out matches to out (may be NULL when
max_out is 0) and returns the TOTAL number of matches, which can exceed
max_out. A NULL or empty query matches nothing.
```

### `qmc_category_count`

```c
size_t qmc_category_count(void);
```

```text
Categories, in display order.
```

### `qmc_category_at`

```c
const char *qmc_category_at(size_t index);
```

_No description in the header yet._

### `qmc_help`

```c
size_t qmc_help(FILE *out, const char *topic);
```

```text
Prints help to out (stdout if NULL):
  topic NULL or ""  -> overview of all categories and modules
  a module name     -> what it does, #include line, where to start, guide,
                       example, test and API-reference location
  a category (or a part of its name, at least 4 letters) -> its modules
  anything else     -> keyword search; suggests qmc_help(out, NULL) if empty.
Returns the number of modules shown.
```
