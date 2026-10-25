/*
 * ASCII terminal plotter: render 1D data as a text-mode line plot
 *
 * Purpose: quick visual inspection of numeric data without requiring an
 * external graphics library. Computes y-axis min/max, lays out a character grid
 * with a '+' origin at (x0, y0) = (2, height - 2), draws x and y axes as '-'
 * and '|', and marks each data point with '*'
 *
 * TEST: Implement ascii plotter for terminal
 * TODO: Implement 3d ploting (physic engine, example github origo)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Simple line plot for 1D data
void ascii_plot(const double *x, const double *y, int n, int height,
                int width) {
  if (!x || !y || n < 2 || height < 3 || width < 10) {
    return;
  }

  // Find min/max
  double ymin = y[0];
  double ymax = y[0];
  for (int i = 1; i < n; i++) {
    if (y[i] < ymin) {
      ymin = y[i];
    }
    if (y[i] > ymax) {
      ymax = y[i];
    }
  }

  double yrange = ymax - ymin;
  if (yrange < 1e-15) {
    yrange = 1.0;
  }

  // Allocate plot buffer
  char **plot = malloc((size_t)height * sizeof *plot);
  if (!plot) {
    fprintf(stderr, "Memory allocation failed\n");

    return;
  }

  for (int i = 0; i < height; i++) {
    plot[i] = malloc((size_t)(width + 1) * sizeof *plot[i]);
    if (!plot[i]) {
      fprintf(stderr, "Memory allocation failed\n");

      for (int j = 0; j < i; j++) {
        free(plot[j]);
      }
      free(plot);

      return;
    }

    memset(plot[i], ' ', (size_t)width);
    plot[i][width] = '\0';
  }

  // Draw axes
  int x0 = 2;
  int y0 = height - 2;
  for (int i = 0; i < width; i++) {
    plot[y0][i] = '-';
  }
  for (int i = 0; i < height; i++) {
    plot[i][x0] = '|';
  }
  plot[y0][x0] = '+';

  // Plot points
  for (int i = 0; i < n; i++) {
    int col = x0 + 1 + (int)((width - x0 - 2) * (i / (double)(n - 1)));
    int row = y0 - 1 - (int)((height - 3) * (y[i] - ymin) / yrange);
    if (col >= 0 && col < width && row >= 0 && row < height) {
      plot[row][col] = '*';
    }
  }

  // Print plot
  for (int i = 0; i < height; i++) {
    printf("%s\n", plot[i]);

    free(plot[i]);
  }
  free(plot);

  printf("Min=%.3g, Max=%.3g\n", ymin, ymax);
}
