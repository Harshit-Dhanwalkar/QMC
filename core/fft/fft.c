/*
Fast Fourier transform
*/
#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Helpers
static unsigned int bit_reverse(unsigned int x, int log2n) {
  unsigned int y = 0;
  for (int i = 0; i < log2n; i++) {
    y <<= 1;
    y |= (x & 1);
    x >>= 1;
  }
  return y;
}

static void bit_reverse_permute(cvector_t *x) {
  int n = x->n;
  int log2n = (int)(log2(n) + 0.5);
  for (int i = 0; i < n; i++) {
    int j = bit_reverse(i, log2n);
    if (j > i) {
      complex_t tmp = x->data[i];
      x->data[i] = x->data[j];
      x->data[j] = tmp;
    }
  }
}

// public FFT routines
void fft(cvector_t *x) {
  int n = x->n;
  if (n < 2)
    return;
  // Check power of two
  if ((n & (n - 1)) != 0) {
    // HACK:
    // Not power of two: fallback to naive DFT?
    // For simplicity return; caller must ensure length is power of two.
    return;
  }

  bit_reverse_permute(x);

  for (int len = 2; len <= n; len <<= 1) {
    double angle = -2.0 * M_PI / len;
    complex_t wlen = {cos(angle), sin(angle)};
    for (int i = 0; i < n; i += len) {
      complex_t w = {1.0, 0.0};
      for (int j = 0; j < len / 2; j++) {
        complex_t u = x->data[i + j];
        complex_t v = c_mul(x->data[i + j + len / 2], w);
        x->data[i + j] = c_add(u, v);
        x->data[i + j + len / 2] = c_sub(u, v);
        w = c_mul(w, wlen);
      }
    }
  }
}

void ifft(cvector_t *x) {
  // Conjugate, forward FFT, conjugate, scale by 1/N
  int n = x->n;
  for (int i = 0; i < n; i++) {
    x->data[i].im = -x->data[i].im;
  }
  fft(x);
  for (int i = 0; i < n; i++) {
    x->data[i].im = -x->data[i].im;
    x->data[i].re /= n;
    x->data[i].im /= n;
  }
}

void fft_normalized(cvector_t *x) {
  fft(x);
  double scale = 1.0 / sqrt(x->n);
  for (int i = 0; i < x->n; i++) {
    x->data[i] = c_scale(x->data[i], scale);
  }
}

void ifft_normalized(cvector_t *x) {
  ifft(x);
  double scale = sqrt(x->n);
  for (int i = 0; i < x->n; i++) {
    x->data[i] = c_scale(x->data[i], scale);
  }
}

void fft_shift(cvector_t *x) {
  int n = x->n;
  int half = n / 2;
  cvector_t *tmp = cvector_alloc(n);
  if (!tmp)
    return;
  for (int i = 0; i < half; i++) {
    tmp->data[i] = x->data[i + half];
    tmp->data[i + half] = x->data[i];
  }
  memcpy(x->data, tmp->data, n * sizeof(complex_t));
  cvector_free(tmp);
}

// NOTE: Implemented in utils.c
// cvector_t *position_to_momentum(const cvector_t *psi_x, double dx) {
//   int n = psi_x->n;
//   cvector_t *psi_k = cvector_copy(psi_x); // copy to work on
//   if (!psi_k)
//     return NULL;
//
//   // Apply forward FFT (normalised)
//   fft_normalized(psi_k);
//
//   // Shift zero frequency to centre (so that k=0 is in the middle)
//   fft_shift(psi_k);
//
//   // NOTE:
//   // Compute k array and multiply by phase factor?
//   // \psi(k) = \int \phi(x) e^{-ikx} dx  with k = 2\pi n/(N dx) for n=0..N-1
//   // After shift, k runs from -\pi/dx to +\pi/dx.
//   // The normalisation is already correct such that \int|\phi|^2 dx =
//   // \int|psi|^2 dk because fft_normalized (1/sqrt(N)) and dk = 2\pi/(N dx)
//   return psi_k;
// }
