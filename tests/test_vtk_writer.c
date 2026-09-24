/*
 * Test: VTK writer (export/vtk_writer.c).
 *
 * 1. A hand-rolled parser reads legacy VTK ASCII file back and checks
 *    DIMENSIONS/ORIGIN/SPACING and every data value - a round-trip check
 * 2. A 2D field (nz = 1) and a 3D field, since writer treats them uniformly
 *    and both need checking
 * 3. Invalid-input handling (zero dimension, NULL args) returns -1  without
 *    crashing
 */

#include "../export/vtk_writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

static int failures = 0;

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

static char *slurp(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) {
    return NULL;
  }

  fseek(f, 0, SEEK_END);
  long len = ftell(f);
  fseek(f, 0, SEEK_SET);

  char *buf = malloc((size_t)len + 1);
  if (!buf) {
    fclose(f);
    return NULL;
  }

  size_t got = fread(buf, 1, (size_t)len, f);
  buf[got] = '\0';
  fclose(f);

  return buf;
}

// Minimal independent parser for subset of legacy VTK ASCII format this writer
// produces
typedef struct {
  size_t nx, ny, nz;
  double origin[3];
  double spacing[3];
  double *values;
  size_t n_values;
} parsed_vtk_t;

static int parse_vtk(const char *content, parsed_vtk_t *out) {
  char line[256];
  const char *p = content;

  // Line 1: "# vtk DataFile Version 3.0"
  if (sscanf(p, "%255[^\n]", line) != 1 ||
      strncmp(line, "# vtk DataFile Version", 22) != 0) {
    return 0;
  }

  p = strchr(p, '\n') + 1;
  p = strchr(p, '\n') + 1; // skip title line

  if (sscanf(p, "%255[^\n]", line) != 1 || strcmp(line, "ASCII") != 0) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  if (sscanf(p, "%255[^\n]", line) != 1 ||
      strcmp(line, "DATASET STRUCTURED_POINTS") != 0) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  if (sscanf(p, "DIMENSIONS %zu %zu %zu", &out->nx, &out->ny, &out->nz) != 3) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  if (sscanf(p, "ORIGIN %lf %lf %lf", &out->origin[0], &out->origin[1],
             &out->origin[2]) != 3) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  if (sscanf(p, "SPACING %lf %lf %lf", &out->spacing[0], &out->spacing[1],
             &out->spacing[2]) != 3) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  size_t n_points = 0;
  if (sscanf(p, "POINT_DATA %zu", &n_points) != 1) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  char varname[128];
  char typname[32];
  int ncomp = 0;
  if (sscanf(p, "SCALARS %127s %31s %d", varname, typname, &ncomp) != 3) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  if (sscanf(p, "%255[^\n]", line) != 1 ||
      strcmp(line, "LOOKUP_TABLE default") != 0) {
    return 0;
  }

  p = strchr(p, '\n') + 1;

  out->values = malloc(n_points * sizeof(double));
  out->n_values = n_points;
  for (size_t i = 0; i < n_points; i++) {
    if (sscanf(p, "%lf", &out->values[i]) != 1) {
      return 0;
    }

    p = strchr(p, '\n');
    if (!p) {
      return (i == n_points - 1) ? 1 : 0;
    }

    p++;
  }

  return 1;
}

static void test_2d_round_trip(void) {
  printf("  === Test 2D round trip ===\n");

  // 3x2 grid (nx=3, ny=2), nz=1
  double data[6] = {1.0, 2.0, 3.0, 4.0, 5.0, -6.5};
  int rc = vtk_write_structured_points("test_vtk_2d.vtk", "psi", data, 3, 2, 1,
                                       0.1, 0.2, 1.0, -1.0, -2.0, 0.0);
  check_true(rc == 0, "vtk_write_structured_points returns 0 on success");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, "test_vtk_2d.vtk");
  char *content = slurp(path);
  check_true(content != NULL, "output file was created and is readable");
  if (!content) {
    return;
  }

  parsed_vtk_t parsed = {0};
  int ok = parse_vtk(content, &parsed);
  check_true(ok, "independent parser accepts the written file as valid "
                 "legacy VTK ASCII");

  if (ok) {
    check_true(parsed.nx == 3 && parsed.ny == 2 && parsed.nz == 1,
               "DIMENSIONS round-trips exactly (3 2 1)");
    check_true(parsed.origin[0] == -1.0 && parsed.origin[1] == -2.0 &&
                   parsed.origin[2] == 0.0,
               "ORIGIN round-trips exactly");
    check_true(parsed.spacing[0] == 0.1 && parsed.spacing[1] == 0.2 &&
                   parsed.spacing[2] == 1.0,
               "SPACING round-trips exactly");
    check_true(parsed.n_values == 6, "POINT_DATA count matches nx*ny*nz");

    int all_match = 1;
    for (int i = 0; i < 6; i++) {
      if (parsed.values[i] != data[i]) {
        all_match = 0;
      }
    }

    check_true(all_match,
               "every scalar value round-trips bit-for-bit identical");

    free(parsed.values);
  }

  free(content);
}

static void test_3d_field(void) {
  printf("  === Test 3D field ===\n");

  // 2x2x2 grid
  double data[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int rc = vtk_write_structured_points("test_vtk_3d.vtk", "density", data, 2, 2,
                                       2, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
  check_true(rc == 0, "3D write returns 0 on success");

  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, "test_vtk_3d.vtk");
  char *content = slurp(path);

  parsed_vtk_t parsed = {0};
  int ok = content && parse_vtk(content, &parsed);
  check_true(ok, "independent parser accepts the 3D file");
  if (ok) {
    check_true(parsed.nx == 2 && parsed.ny == 2 && parsed.nz == 2,
               "3D DIMENSIONS round-trips exactly (2 2 2)");
    check_true(parsed.n_values == 8, "POINT_DATA count is nx*ny*nz = 8");

    free(parsed.values);
  }

  free(content);
}

static void test_invalid_input(void) {
  printf("  === Test invalid input ===\n");

  double data[1] = {1.0};
  int rc = vtk_write_structured_points("x.vtk", "v", data, 0, 1, 1, 1, 1, 1, 0,
                                       0, 0);
  check_true(rc == -1, "zero grid dimension is rejected, not silently "
                       "written as an empty/garbage file");

  rc = vtk_write_structured_points(NULL, NULL, NULL, 1, 1, 1, 1, 1, 1, 0, 0, 0);
  check_true(rc == -1, "NULL filename/varname/data does not crash");

  rc = vtk_write_structured_points("/nonexistent_dir_xyz/cant_write.vtk", "v",
                                   data, 1, 1, 1, 1, 1, 1, 0, 0, 0);
  check_true(rc == -1, "unwritable path returns -1 rather than crashing");
}

int main(void) {
  printf(" > Test VTK writer ===\n");

  test_2d_round_trip();
  test_3d_field();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_vtk_writer checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
