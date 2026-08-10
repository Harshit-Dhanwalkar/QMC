// FIX: sanitizer
#include "../export/gr/gr_plot.h"
// #include <gr.h>
#include "../third_party/gr/install/include/gr.h"
#include <math.h>
#include <sanitizer/lsan_interface.h>
#include <stdio.h>
#include <stdlib.h>

// #if __has_include(<sanitizer/lsan_interface.h>)
// #include <sanitizer/lsan_interface.h>
// #else
// void __lsan_disable(void);
// void __lsan_enable(void);
// #endif

int main() {
  // setenv("GKS_WSTYPE", "nul", 0);
  // GKS_WSTYPE=nul is set centrally inside gr_plot_to_file itself

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

// #if defined(__SANITIZER_ADDRESS__) || defined(__LEAK_SANITIZER__)
//   __lsan_disable();
// #endif

  __lsan_disable();
  gr_plot_to_file("gaussian", GR_FORMAT_PDF, x, y, N, &opts);
  __lsan_enable();

// #if defined(__SANITIZER_ADDRESS__) || defined(__LEAK_SANITIZER__)
//   __lsan_enable();
// #endif

  gr_plot_to_file("gaussian", GR_FORMAT_SVG, x, y, N, &opts);
  gr_plot_to_file("gaussian", GR_FORMAT_PNG, x, y, N, &opts);
  gr_plot_to_file("gaussian", GR_FORMAT_JPEG, x, y, N, &opts);

  // gr_closegks();

  return 0;
}
