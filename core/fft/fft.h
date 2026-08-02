#ifndef QMC_FFT_H
#define QMC_FFT_H

#include "../vector.h"

/*
 * Fast Fourier Transform (Cooley‑Tukey, radix‑2, iterative)
 * Works for power‑of‑two lengths.
 */

/* In‑place forward FFT (no normalization) */
void fft(cvector_t *x);

/* In‑place inverse FFT (no normalization) */
void ifft(cvector_t *x);

/* Normalized forward FFT (1 / \sqrt(N) factor) */
void fft_normalized(cvector_t *x);

/* Normalized inverse FFT (1 / \sqrt(N) factor) */
void ifft_normalized(cvector_t *x);

/* Utility: shift zero‑frequency component to centre of array */
void fft_shift(cvector_t *x);

/* Convert position‑space wavefunction to momentum space.
 * dx : grid spacing in position space.
 * The output array is ordered as k = 2 * \pi ** n / (N * dx) with n=0..N-1.
 * Normalisation is such that \int |\phi(x)|^2 dx = \int|\phi(k)|^2 dk.

 * Returns a newly allocated cvector_t in momentum space.
*/
cvector_t *position_to_momentum(const cvector_t *psi_x, double dx);

#endif
