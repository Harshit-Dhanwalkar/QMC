#ifndef QMC_EXPORT_SCHEMA_H
#define QMC_EXPORT_SCHEMA_H

#include "../core/qmc_version.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Conventions shared by every exporter (JSON, HDF5, CSV sidecars)
 *
 * schema_version describes layout of QMC's exported files and only changes on a
 * breaking layout change */
#define QMC_EXPORT_SCHEMA_VERSION 1

#define QMC_LIBRARY_NAME "QMC"

/* Derived from core/qmc_version.h, single source of truth for version */
#define QMC_LIBRARY_VERSION QMC_VERSION_STRING

/* Buffer size that always fits an ISO 8601 UTC timestamp + NUL */
#define EXPORT_TIMESTAMP_SIZE 32

/* Writes current time as "YYYY-MM-DDTHH:MM:SSZ" (UTC) into buf
 *
 * Returns 0 on success, -1 if buf is NULL/too small or clock fails */
int export_timestamp_utc(char *buf, size_t size);

/* Description of one tabular column. Shared by CSV writer (header row) and
 * JSON / HDF5 writers (unit + description metadata), so same descriptor array
 * can drive all three outputs. `unit` and `description` may be NULL */
typedef struct {
  const char *name;
  const char *unit;
  const char *description;
} export_column_t;

#ifdef __cplusplus
}
#endif

#endif /* QMC_EXPORT_SCHEMA_H */
