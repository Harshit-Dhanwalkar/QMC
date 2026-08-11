// Normalisation
#include "../core/complex.h"
#include "../core/fft/fft.h"
#include "../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  printf(" > Testing FFT...\n");

  int N = 16;
  cvector_t *x = cvector_alloc(N);
  if (!x) {
    return 1;
  }

  for (int i = 0; i < N; i++) {
    double t = 2.0 * M_PI * i / N;

    x->data[i] = c_new(cos(t), sin(t)); // e^{it}
  }

  fft(x);

  // After FFT, only first bin should have peak.
  double max_abs = 0.0;
  int max_idx = -1;
  for (int i = 0; i < N; i++) {
    double a = c_abs(x->data[i]);

    if (a > max_abs) {
      max_abs = a;
      max_idx = i;
    }
  }

  printf("   Max peak at index %d with abs %f\n", max_idx, max_abs);
  if (max_idx == 1) {
    printf("   FFT test passed.\n");
  } else {
    printf("   FFT test failed.\n");
  }
  cvector_free(x);

  return 0;
}
