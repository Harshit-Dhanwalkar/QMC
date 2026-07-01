/*
Link to FFTW LIB
*/
#include "fft.h" // built‑in FFT
#include "fft_wrapper.h"

#ifdef USE_FFTW
#include <fftw3.h>
#endif

void fft_wrapper(cvector_t *x) {
#ifdef USE_FFTW
  int n = x->n;
  fftw_complex *in = fftw_malloc(n * sizeof(fftw_complex));
  fftw_complex *out = fftw_malloc(n * sizeof(fftw_complex));
  for (int i = 0; i < n; i++) {
    in[i][0] = x->data[i].re;
    in[i][1] = x->data[i].im;
  }
  fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  for (int i = 0; i < n; i++) {
    x->data[i].re = out[i][0];
    x->data[i].im = out[i][1];
  }
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
#else
  fft(x);
#endif
}

void ifft_wrapper(cvector_t *x) {
#ifdef USE_FFTW
  int n = x->n;
  fftw_complex *in = fftw_malloc(n * sizeof(fftw_complex));
  fftw_complex *out = fftw_malloc(n * sizeof(fftw_complex));
  for (int i = 0; i < n; i++) {
    in[i][0] = x->data[i].re;
    in[i][1] = x->data[i].im;
  }
  fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  for (int i = 0; i < n; i++) {
    x->data[i].re = out[i][0] / n;
    x->data[i].im = out[i][1] / n;
  }
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
#else
  ifft(x);
#endif
}

void fft_wrapper_normalized(cvector_t *x) {
#ifdef USE_FFTW
  fft_wrapper(x);
  double scale = 1.0 / sqrt(x->n);
  for (int i = 0; i < x->n; i++)
    x->data[i] = c_scale(x->data[i], scale);
#else
  fft_normalized(x);
#endif
}

void ifft_wrapper_normalized(cvector_t *x) {
#ifdef USE_FFTW
  ifft_wrapper(x);
  double scale = sqrt(x->n);
  for (int i = 0; i < x->n; i++)
    x->data[i] = c_scale(x->data[i], scale);
#else
  ifft_normalized(x);
#endif
}
