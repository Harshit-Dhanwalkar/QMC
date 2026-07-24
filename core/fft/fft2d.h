#ifndef QMC_FFT2D_H
#define QMC_FFT2D_H

#include "../vector.h"

/*
 * 2D FFT built from 1D radix-2 fft()/ifft() via standard row-column (separable)
 * algorithm: FFT every row, then FFT every column
 */

// Forward 2D FFT, unnormalized
void fft2d(cvector_t *psi, int Nx, int Ny);

// Inverse 2D FFT, normalized by 1/(Nx*Ny)
void ifft2d(cvector_t *psi, int Nx, int Ny);

#endif // QMC_FFT2D_H
