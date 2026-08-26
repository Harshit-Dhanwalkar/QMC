// Normalisation
#include "../core/complex.h"
#include "../core/fft/fft.h"
#include "../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static int test_fft_power_of_two_peak(void) {
  printf("Test: fft() of a pure exp(i*t) tone peaks at bin 1 (power-of-two "
         "radix-2 path)\n");

  int N = 16;
  cvector_t *x = cvector_alloc(N);
  for (int i = 0; i < N; i++) {
    double t = 2.0 * M_PI * i / N;

    x->data[i] = c_new(cos(t), sin(t)); // \exp^{it}
    x->data[i] = c_new(cos(t), sin(t));
  }

  fft(x);

  double max_abs = 0.0;
  int max_idx = -1;
  for (int i = 0; i < N; i++) {
    double a = c_abs(x->data[i]);

    if (a > max_abs) {
      max_abs = a;
      max_idx = i;
    }
  }

  check(max_idx == 1, "FFT peak lands at bin 1");
  cvector_free(x);

  return max_idx == 1 ? 0 : 1;
}

static void test_naive_dft_matches_reference(void) {
  printf("Test: fft() on a non-power-of-two length matches an independent "
         "reference DFT computed directly");

  int n = 17; // not a power of two
  cvector_t *x = cvector_alloc(n);
  for (int i = 0; i < n; i++) {
    x->data[i] = c_new(sin(0.3 * i) + 0.5 * cos(0.7 * i), 0.0);
  }

  cvector_t *ref = cvector_copy(x);

  fft(x);

  for (int k = 0; k < n; k++) {
    double re = 0.0, im = 0.0;

    for (int j = 0; j < n; j++) {
      double angle = -2.0 * M_PI * k * j / n;

      re += ref->data[j].re * cos(angle);
      im += ref->data[j].re * sin(angle);
    }

    double err = sqrt((re - x->data[k].re) * (re - x->data[k].re) +
                      (im - x->data[k].im) * (im - x->data[k].im));
    check(err < 1e-9, "non-power-of-two fft() output matches reference DFT");
  }

  cvector_free(x);
  cvector_free(ref);
}

static void test_naive_dft_not_a_noop(void) {
  printf("Test: fft() on a non-power-of-two length actually changes the "
         "data");

  int n = 11;
  cvector_t *x = cvector_alloc(n);
  for (int i = 0; i < n; i++) {
    x->data[i] = c_new(1.0 + i, 0.0); /* not already its own DFT */
  }

  cvector_t *before = cvector_copy(x);

  fft(x);

  int identical = 1;
  for (int i = 0; i < n; i++) {
    if (fabs(x->data[i].re - before->data[i].re) > 1e-12 ||
        fabs(x->data[i].im - before->data[i].im) > 1e-12) {

      identical = 0;
      break;
    }
  }

  check(!identical, "fft() transforms non-power-of-two input");

  cvector_free(x);
  cvector_free(before);
}

static void test_fft_shift_matches_numpy_convention(void) {
  printf("Test: fft_shift matches numpy.fft.fftshift's roll(x, n/2) "
         "convention for both even and odd n");

  // Even n=6: fftshift([0,1,2,3,4,5]) = [3,4,5,0,1,2]
  {
    int n = 6;
    cvector_t *x = cvector_alloc(n);
    for (int i = 0; i < n; i++) {
      x->data[i] = c_real((double)i);
    }

    fft_shift(x);
    const double expected[6] = {3, 4, 5, 0, 1, 2};
    for (int i = 0; i < n; i++) {
      check(fabs(x->data[i].re - expected[i]) < 1e-12,
            "even-n fft_shift matches numpy roll(x, n/2)");
    }

    cvector_free(x);
  }

  // Odd n=5: fftshift([0,1,2,3,4]) = [3,4,0,1,2]
  {
    int n = 5;
    cvector_t *x = cvector_alloc(n);
    for (int i = 0; i < n; i++) {
      x->data[i] = c_real((double)i);
    }

    fft_shift(x);
    const double expected[5] = {3, 4, 0, 1, 2};
    for (int i = 0; i < n; i++) {
      check(fabs(x->data[i].re - expected[i]) < 1e-12,
            "odd-n fft_shift matches numpy roll(x, n/2)");
    }

    cvector_free(x);
  }
}

int main(void) {
  printf(" > Testing FFT...\n\n");

  int peak_failed = test_fft_power_of_two_peak();

  test_naive_dft_matches_reference();
  test_naive_dft_not_a_noop();
  test_fft_shift_matches_numpy_convention();

  int total_failed = failures + peak_failed;
  if (total_failed == 0) {
    printf("\nAll test_dft checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_dft check(s) FAILED.\n", total_failed);
    return 1;
  }
}
