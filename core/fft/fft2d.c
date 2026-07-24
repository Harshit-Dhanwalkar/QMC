/*
2D FFT via row-column decomposition of 1D fft()/ifft().
*/

#include "fft2d.h"
#include "../vector.h"
#include "fft.h"
#include <stdlib.h>

static void transform_rows(cvector_t *psi, int Nx, int Ny, int inverse) {
  cvector_t *row = cvector_alloc(Ny);
  if (!row) {
    return;
  }

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      row->data[iy] = psi->data[ix * Ny + iy];
    }

    if (inverse) {
      ifft(row);
    } else {
      fft(row);
    }

    for (int iy = 0; iy < Ny; iy++) {
      psi->data[ix * Ny + iy] = row->data[iy];
    }
  }

  cvector_free(row);
}

static void transform_cols(cvector_t *psi, int Nx, int Ny, int inverse) {
  cvector_t *col = cvector_alloc(Nx);
  if (!col) {
    return;
  }

  for (int iy = 0; iy < Ny; iy++) {
    for (int ix = 0; ix < Nx; ix++) {
      col->data[ix] = psi->data[ix * Ny + iy];
    }

    if (inverse) {
      ifft(col);
    } else {
      fft(col);
    }

    for (int ix = 0; ix < Nx; ix++) {
      psi->data[ix * Ny + iy] = col->data[ix];
    }
  }

  cvector_free(col);
}

void fft2d(cvector_t *psi, int Nx, int Ny) {
  if (!psi || Nx < 1 || Ny < 1 || psi->n != Nx * Ny) {
    return;
  }

  transform_rows(psi, Nx, Ny, 0);
  transform_cols(psi, Nx, Ny, 0);
}

void ifft2d(cvector_t *psi, int Nx, int Ny) {
  if (!psi || Nx < 1 || Ny < 1 || psi->n != Nx * Ny) {
    return;
  }

  transform_rows(psi, Nx, Ny, 1);
  transform_cols(psi, Nx, Ny, 1);
}
