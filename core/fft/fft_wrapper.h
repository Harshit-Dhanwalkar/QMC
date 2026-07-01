#ifndef QMC_FFT_WRAPPER_H
#define QMC_FFT_WRAPPER_H

#include "../vector.h"

/* FFT wrapper that uses either internal or external (FFTW) library */
void fft_wrapper(cvector_t *x);
void ifft_wrapper(cvector_t *x);
void fft_wrapper_normalized(cvector_t *x);
void ifft_wrapper_normalized(cvector_t *x);

#endif
