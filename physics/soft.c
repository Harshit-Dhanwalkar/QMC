/*
Split-Operator Fourier Transform (SOFT) method, 2D and 3D.
*/

#include "soft.h"
#include "../core/complex.h"
#include "../core/fft/fft2d.h"
#include "../core/fft/fft3d.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>

int soft_evolve_2d(cvector_t *psi, const double *V, int Nx, int Ny, double dx,
                   double dy, double dt, int steps, double hbar, double mass) {
  if (!psi || !V || Nx < 2 || Ny < 2 || steps < 1 || mass <= 0.0 ||
      psi->n != Nx * Ny) {
    return -1;
  }

  double dkx = 2.0 * M_PI / (Nx * dx);
  double dky = 2.0 * M_PI / (Ny * dy);
  double hbar_over_2m = hbar * hbar / (2.0 * mass);

  double *kx2 = malloc(Nx * sizeof *kx2);
  double *ky2 = malloc(Ny * sizeof *ky2);
  if (!kx2 || !ky2) {
    free(kx2);
    free(ky2);

    return -1;
  }

  for (int ix = 0; ix < Nx; ix++) {
    double kx = (ix < Nx / 2) ? ix * dkx : (ix - Nx) * dkx;
    kx2[ix] = kx * kx;
  }

  for (int iy = 0; iy < Ny; iy++) {
    double ky = (iy < Ny / 2) ? iy * dky : (iy - Ny) * dky;
    ky2[iy] = ky * ky;
  }

  int n = Nx * Ny;
  for (int s = 0; s < steps; s++) {
    for (int i = 0; i < n; i++) {
      double phase = -V[i] * dt / 2.0;
      psi->data[i] = c_mul(psi->data[i], c_exp(c_imag(phase)));
    }

    fft2d(psi, Nx, Ny);
    for (int ix = 0; ix < Nx; ix++) {
      for (int iy = 0; iy < Ny; iy++) {
        int idx = ix * Ny + iy;
        double phase = -hbar_over_2m * (kx2[ix] + ky2[iy]) * dt;
        psi->data[idx] = c_mul(psi->data[idx], c_exp(c_imag(phase)));
      }
    }
    ifft2d(psi, Nx, Ny);

    for (int i = 0; i < n; i++) {
      double phase = -V[i] * dt / 2.0;
      psi->data[i] = c_mul(psi->data[i], c_exp(c_imag(phase)));
    }
  }

  free(kx2);
  free(ky2);

  return 0;
}

int soft_evolve_3d(cvector_t *psi, const double *V, int Nx, int Ny, int Nz,
                   double dx, double dy, double dz, double dt, int steps,
                   double hbar, double mass) {
  if (!psi || !V || Nx < 2 || Ny < 2 || Nz < 2 || steps < 1 || mass <= 0.0 ||
      psi->n != Nx * Ny * Nz) {
    return -1;
  }

  double dkx = 2.0 * M_PI / (Nx * dx);
  double dky = 2.0 * M_PI / (Ny * dy);
  double dkz = 2.0 * M_PI / (Nz * dz);
  double hbar_over_2m = hbar * hbar / (2.0 * mass);

  double *kx2 = malloc(Nx * sizeof *kx2);
  double *ky2 = malloc(Ny * sizeof *ky2);
  double *kz2 = malloc(Nz * sizeof *kz2);
  if (!kx2 || !ky2 || !kz2) {
    free(kx2);
    free(ky2);
    free(kz2);

    return -1;
  }

  for (int ix = 0; ix < Nx; ix++) {
    double kx = (ix < Nx / 2) ? ix * dkx : (ix - Nx) * dkx;
    kx2[ix] = kx * kx;
  }

  for (int iy = 0; iy < Ny; iy++) {
    double ky = (iy < Ny / 2) ? iy * dky : (iy - Ny) * dky;
    ky2[iy] = ky * ky;
  }

  for (int iz = 0; iz < Nz; iz++) {
    double kz = (iz < Nz / 2) ? iz * dkz : (iz - Nz) * dkz;
    kz2[iz] = kz * kz;
  }

  int n = Nx * Ny * Nz;
  for (int s = 0; s < steps; s++) {
    for (int i = 0; i < n; i++) {
      double phase = -V[i] * dt / 2.0;
      psi->data[i] = c_mul(psi->data[i], c_exp(c_imag(phase)));
    }

    fft3d(psi, Nx, Ny, Nz);
    for (int ix = 0; ix < Nx; ix++) {
      for (int iy = 0; iy < Ny; iy++) {
        for (int iz = 0; iz < Nz; iz++) {
          int idx = (ix * Ny + iy) * Nz + iz;
          double phase = -hbar_over_2m * (kx2[ix] + ky2[iy] + kz2[iz]) * dt;
          psi->data[idx] = c_mul(psi->data[idx], c_exp(c_imag(phase)));
        }
      }
    }
    ifft3d(psi, Nx, Ny, Nz);

    for (int i = 0; i < n; i++) {
      double phase = -V[i] * dt / 2.0;
      psi->data[i] = c_mul(psi->data[i], c_exp(c_imag(phase)));
    }
  }

  free(kx2);
  free(ky2);
  free(kz2);

  return 0;
}

double grid_norm(const cvector_t *psi, double cell_volume) {
  if (!psi) {
    return 0.0;
  }

  double sum = 0.0;
  for (int i = 0; i < psi->n; i++) {
    sum += c_abs2(psi->data[i]);
  }

  return sum * cell_volume;
}
