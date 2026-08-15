#include "../export/gr/gr_plot.h"
#include "../third_party/gr/install/include/gr.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef SANITIZE_ENABLED
#include <sanitizer/lsan_interface.h>
#define LSAN_DISABLE() __lsan_disable()
#define LSAN_ENABLE() __lsan_enable()
#else
#define LSAN_DISABLE() ((void)0)
#define LSAN_ENABLE() ((void)0)
#endif

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

  LSAN_DISABLE();
  gr_plot_to_file("gaussian", GR_FORMAT_PDF, x, y, N, &opts);
  LSAN_ENABLE();

  gr_plot_to_file("gaussian", GR_FORMAT_SVG, x, y, N, &opts);
  gr_plot_to_file("gaussian", GR_FORMAT_PNG, x, y, N, &opts);
  gr_plot_to_file("gaussian", GR_FORMAT_JPEG, x, y, N, &opts);

  return 0;
}
