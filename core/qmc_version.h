#ifndef QMC_VERSION_H
#define QMC_VERSION_H

/*
 * QMC library version: single source of truth
 *
 * Makefile reads QMC_VERSION_STRING from this file for qmc.pc and
 * shared-library name (libqmc.so.<MAJOR>), and export/export_schema.h derives
 * QMC_LIBRARY_VERSION from it, so a release is a one-place edit: change
 * four numbers/strings below together (tests/test_version.c checks that they
 * agree)
 */
#define QMC_VERSION_MAJOR 0
#define QMC_VERSION_MINOR 1
#define QMC_VERSION_PATCH 0
#define QMC_VERSION_STRING "0.1.0"

// MAJOR * 10000 + MINOR * 100 + PATCH
#define QMC_VERSION_NUMBER                                                     \
  (QMC_VERSION_MAJOR * 10000 + QMC_VERSION_MINOR * 100 + QMC_VERSION_PATCH)

#ifdef __cplusplus
extern "C" {
#endif

// Version of the compiled library (not of the headers): "MAJOR.MINOR.PATCH"
const char *qmc_version(void);

// Same as QMC_VERSION_NUMBER, but of the compiled library
int qmc_version_number(void);

#ifdef __cplusplus
}
#endif

#endif /* QMC_VERSION_H */
