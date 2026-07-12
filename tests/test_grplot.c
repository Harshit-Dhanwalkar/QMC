#include "../export/gr/gr_plot.h"
#include <gr.h>
#include <math.h>
#include <stdio.h>

int main() {
  int N = 100;
  double x[100], y[100];
  for (int i = 0; i < N; i++) {
    x[i] = (double)i / N * 10.0 - 5.0;
    y[i] = exp(-x[i] * x[i]);
  }
  gr_plot_opt_t opts = {0};
  opts.title = "Gaussian";
  opts.xlabel = "x";
  opts.ylabel = "f(x)";
  opts.width = 800;
  opts.height = 600;
  gr_plot_to_file("gaussian.pdf", GR_FORMAT_PDF, x, y, N, &opts);
  gr_plot_to_file("gaussian.svg", GR_FORMAT_SVG, x, y, N, &opts);
  gr_plot_to_file("gaussian.png", GR_FORMAT_PNG, x, y, N, &opts);
  gr_plot_to_file("gaussian.jpg", GR_FORMAT_JPEG, x, y, N, &opts);

  gr_closegks();

  return 0;
}
