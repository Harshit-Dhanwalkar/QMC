#ifndef QMC_FFT3D_H
#define QMC_FFT3D_H

#include "../vector.h"

/*
 * 3D FFT via row-column-depth (separable) decomposition; transform every line
 * along z, then every line along y, then every line along x.
 */

void fft3d(cvector_t *psi, int Nx, int Ny, int Nz);

void ifft3d(cvector_t *psi, int Nx, int Ny, int Nz);

#endif // QMC_FFT3D_H
