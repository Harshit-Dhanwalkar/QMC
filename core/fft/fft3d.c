/*
3D FFT via row-column-depth decomposition
*/

#include "fft3d.h"
#include "../vector.h"
#include "fft.h"
#include <stdlib.h>

// Transform every line along z (stride 1)
static void transform_z(cvector_t *psi, int Nx, int Ny, int Nz, int inverse) {
  cvector_t *line = cvector_alloc(Nz);
  if (!line) {
    return;
  }

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      int base = (ix * Ny + iy) * Nz;
      for (int iz = 0; iz < Nz; iz++) {
        line->data[iz] = psi->data[base + iz];
      }

      if (inverse) {
        ifft(line);
      } else {
        fft(line);
      }

      for (int iz = 0; iz < Nz; iz++) {
        psi->data[base + iz] = line->data[iz];
      }
    }
  }

  cvector_free(line);
}

// Transform every line along y (stride Nz)
static void transform_y(cvector_t *psi, int Nx, int Ny, int Nz, int inverse) {
  cvector_t *line = cvector_alloc(Ny);
  if (!line) {
    return;
  }

  for (int ix = 0; ix < Nx; ix++) {
    for (int iz = 0; iz < Nz; iz++) {
      int base = ix * Ny * Nz + iz;
      for (int iy = 0; iy < Ny; iy++) {
        line->data[iy] = psi->data[base + iy * Nz];
      }

      if (inverse) {
        ifft(line);
      } else {
        fft(line);
      }

      for (int iy = 0; iy < Ny; iy++) {
        psi->data[base + iy * Nz] = line->data[iy];
      }
    }
  }

  cvector_free(line);
}

// Transform every line along x (stride Ny*Nz)
static void transform_x(cvector_t *psi, int Nx, int Ny, int Nz, int inverse) {
  cvector_t *line = cvector_alloc(Nx);
  if (!line) {
    return;
  }

  for (int iy = 0; iy < Ny; iy++) {
    for (int iz = 0; iz < Nz; iz++) {
      int base = iy * Nz + iz;
      int stride = Ny * Nz;
      for (int ix = 0; ix < Nx; ix++) {
        line->data[ix] = psi->data[base + ix * stride];
      }

      if (inverse) {
        ifft(line);
      } else {
        fft(line);
      }

      for (int ix = 0; ix < Nx; ix++) {
        psi->data[base + ix * stride] = line->data[ix];
      }
    }
  }

  cvector_free(line);
}

void fft3d(cvector_t *psi, int Nx, int Ny, int Nz) {
  if (!psi || Nx < 1 || Ny < 1 || Nz < 1 || psi->n != Nx * Ny * Nz) {
    return;
  }

  transform_z(psi, Nx, Ny, Nz, 0);
  transform_y(psi, Nx, Ny, Nz, 0);
  transform_x(psi, Nx, Ny, Nz, 0);
}

void ifft3d(cvector_t *psi, int Nx, int Ny, int Nz) {
  if (!psi || Nx < 1 || Ny < 1 || Nz < 1 || psi->n != Nx * Ny * Nz) {
    return;
  }

  transform_z(psi, Nx, Ny, Nz, 1);
  transform_y(psi, Nx, Ny, Nz, 1);
  transform_x(psi, Nx, Ny, Nz, 1);
}
