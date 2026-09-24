#include "vtk_writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

// Shortest decimal precision (1 to 17 significant digits) that round-trips
// exactly via strtod
static void write_double(FILE *f, double d) {
  char buf[64];
  for (int precision = 1; precision <= 17; precision++) {
    snprintf(buf, sizeof buf, "%.*g", precision, d);
    if (strtod(buf, NULL) == d) {
      break;
    }
  }

  fputs(buf, f);
}

int vtk_write_structured_points(const char *filename, const char *varname,
                                const double *data, size_t nx, size_t ny,
                                size_t nz, double spacing_x, double spacing_y,
                                double spacing_z, double origin_x,
                                double origin_y, double origin_z) {
  if (!filename || !varname || !data || nx == 0 || ny == 0 || nz == 0) {
    fprintf(stderr, "vtk_writer: NULL filename/varname/data, or a zero grid "
                    "dimension\n");
    return -1;
  }

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);

  FILE *f = fopen(path, "w");
  if (!f) {
    return -1;
  }

  fputs("# vtk DataFile Version 3.0\n", f);
  fprintf(f, "QMC scalar field: %s\n", varname);
  fputs("ASCII\n", f);
  fputs("DATASET STRUCTURED_POINTS\n", f);
  fprintf(f, "DIMENSIONS %zu %zu %zu\n", nx, ny, nz);

  fputs("ORIGIN ", f);
  write_double(f, origin_x);
  fputc(' ', f);
  write_double(f, origin_y);
  fputc(' ', f);
  write_double(f, origin_z);
  fputc('\n', f);

  fputs("SPACING ", f);
  write_double(f, spacing_x);
  fputc(' ', f);
  write_double(f, spacing_y);
  fputc(' ', f);
  write_double(f, spacing_z);
  fputc('\n', f);

  size_t n_points = nx * ny * nz;
  fprintf(f, "POINT_DATA %zu\n", n_points);
  fprintf(f, "SCALARS %s double 1\n", varname);
  fputs("LOOKUP_TABLE default\n", f);

  // VTK's structured-points point ordering is x fastest, then y, then z
  for (size_t i = 0; i < n_points; i++) {
    write_double(f, data[i]);
    fputc('\n', f);
  }

  fclose(f);

  return 0;
}

