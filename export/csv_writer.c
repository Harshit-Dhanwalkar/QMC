#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

// Write 1D data to CSV with header
int csv_write_1d(const char *filename, const double *x, const double *y, int n,
                 const char *xlabel, const char *ylabel) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  FILE *f = fopen(path, "w");
  if (!f) {
    return -1;
  }

  fprintf(f, "%s,%s\n", xlabel ? xlabel : "x", ylabel ? ylabel : "y");
  for (int i = 0; i < n; i++) {
    fprintf(f, "% .6e,% .6e\n", x[i], y[i]);
  }

  fclose(f);

  return 0;
}

// Write 2D matrix (row-major) to CSV
int csv_write_matrix(const char *filename, const double *data, int rows,
                     int cols, const char **col_headers) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  FILE *f = fopen(path, "w");
  if (!f) {
    return -1;
  }

  if (col_headers) {
    for (int j = 0; j < cols; j++) {
      fprintf(f, "%s%s", col_headers[j], (j == cols - 1) ? "\n" : ",");
    }
  } else {
    for (int j = 0; j < cols; j++) {
      fprintf(f, "col%d%s", j, (j == cols - 1) ? "\n" : ",");
    }
  }

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      fprintf(f, "% .6e%s", data[i * cols + j], (j == cols - 1) ? "\n" : ",");
    }
  }

  fclose(f);

  return 0;
}
