/*
Test : Split-Operator Fourier Transform (SOFT) method, 2D and 3D

Validation strategy:
1. Free-particle Gaussian wavepacket (V=0), 2D and 3D: a minimum-uncertainty
   Gaussian wavepacket spreads with a well-known closed form (e.g. Griffiths
   "Intro to QM" ch. 2; Cohen-Tannoudji complement G_I):
     <x>(t)  = x0 + (\hbar * k0/ m ) * t
     Var(x)(t) = \sigma_0^2 * (1 + (\hbar * t / (2 * m * \sigma_0^2))^2)
   independently in every Cartesian direction. This isolates and validates
   momentum-space kinetic propagator (fft2d/ fft3d + phase factors) with
   potential term switched off entirely.
2. 2D harmonic-oscillator coherent state: a Gaussian of width \sigma_0 =
   \sqrt(\hbar / (m * \omega)), displaced from origin, is an exact coherent
   state of the HO. It does NOT spread, and its centroid follows the trajectory
   exactly (Ehrenfest's theorem is exact for a harmonic potential):
   <x>(t) = x_0 * \cos(\omega *t),   <y>(t) = y_0 * \cos(\omega * t)
   The potential-operator half-steps together with kinetic propagator in same
   run.
3. Norm conservation (unitarity): Tr(|\psi|^2) stays 1 throughout every run.
*/

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/soft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

// Test 1: free-particle Gaussian spreading, 2D
static int test_free_particle_2d(void) {
  int Nx = 64, Ny = 64;
  double Lx = 24.0, Ly = 24.0;
  double dx = Lx / Nx, dy = Ly / Ny;
  double hbar = 1.0, mass = 1.0;
  double sigma0 = 1.2, x0 = -2.0, y0 = 0.5, kx0 = 0.8, ky0 = -0.3;

  double *x = malloc(Nx * sizeof(double));
  double *y = malloc(Ny * sizeof(double));
  for (int i = 0; i < Nx; i++) {
    x[i] = -Lx / 2 + i * dx;
  }
  for (int i = 0; i < Ny; i++) {
    y[i] = -Ly / 2 + i * dy;
  }

  cvector_t *psi = cvector_alloc(Nx * Ny);
  double *V = calloc(Nx * Ny, sizeof(double)); // free particle: V=0

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      double gx = exp(-(x[ix] - x0) * (x[ix] - x0) / (4 * sigma0 * sigma0));
      double gy = exp(-(y[iy] - y0) * (y[iy] - y0) / (4 * sigma0 * sigma0));
      double phase = kx0 * x[ix] + ky0 * y[iy];
      psi->data[ix * Ny + iy] = c_scale(c_new(cos(phase), sin(phase)), gx * gy);
    }
  }

  double scale = 1.0 / sqrt(grid_norm(psi, dx * dy));
  for (int i = 0; i < Nx * Ny; i++) {
    psi->data[i] = c_scale(psi->data[i], scale);
  }

  double T = 3.0, dt = 0.01;
  int fail = soft_evolve_2d(psi, V, Nx, Ny, dx, dy, dt, (int)(T / dt), hbar,
                            mass) != 0;

  double mean_x = 0.0, mean_y = 0.0;
  for (int ix = 0; ix < Nx; ix++)
    for (int iy = 0; iy < Ny; iy++) {
      double p = c_abs2(psi->data[ix * Ny + iy]) * dx * dy;
      mean_x += x[ix] * p;
      mean_y += y[iy] * p;
    }

  double var_x = 0.0, var_y = 0.0;
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      double p = c_abs2(psi->data[ix * Ny + iy]) * dx * dy;
      var_x += (x[ix] - mean_x) * (x[ix] - mean_x) * p;
      var_y += (y[iy] - mean_y) * (y[iy] - mean_y) * p;
    }
  }

  fail |= check_close(mean_x, x0 + (hbar * kx0 / mass) * T, 1e-3,
                      "2D free particle <x>(T)");
  fail |= check_close(mean_y, y0 + (hbar * ky0 / mass) * T, 1e-3,
                      "2D free particle <y>(T)");
  double var_exact =
      sigma0 * sigma0 * (1.0 + pow(hbar * T / (2 * mass * sigma0 * sigma0), 2));
  fail |= check_close(var_x, var_exact, 1e-2, "2D free particle Var(x)(T)");
  fail |= check_close(var_y, var_exact, 1e-2, "2D free particle Var(y)(T)");
  fail |= check_close(grid_norm(psi, dx * dy), 1.0, 1e-6,
                      "2D free particle norm conservation");

  cvector_free(psi);
  free(V);
  free(x);
  free(y);

  return fail;
}

// Test 2: free-particle Gaussian spreading, 3D
static int test_free_particle_3d(void) {
  int Nx = 32, Ny = 32, Nz = 32;
  double Lx = 20.0, Ly = 20.0, Lz = 20.0;
  double dx = Lx / Nx, dy = Ly / Ny, dz = Lz / Nz;
  double hbar = 1.0, mass = 1.0;
  double sigma0 = 1.3, x0 = -1.5, y0 = 0.5, z0 = 1.0;
  double kx0 = 0.5, ky0 = -0.4, kz0 = 0.3;

  double *x = malloc(Nx * sizeof(double));
  double *y = malloc(Ny * sizeof(double));
  double *z = malloc(Nz * sizeof(double));
  for (int i = 0; i < Nx; i++) {
    x[i] = -Lx / 2 + i * dx;
  }
  for (int i = 0; i < Ny; i++) {
    y[i] = -Ly / 2 + i * dy;
  }
  for (int i = 0; i < Nz; i++) {
    z[i] = -Lz / 2 + i * dz;
  }

  int N = Nx * Ny * Nz;
  cvector_t *psi = cvector_alloc(N);
  double *V = calloc(N, sizeof(double));

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int iz = 0; iz < Nz; iz++) {
        double gx = exp(-(x[ix] - x0) * (x[ix] - x0) / (4 * sigma0 * sigma0));
        double gy = exp(-(y[iy] - y0) * (y[iy] - y0) / (4 * sigma0 * sigma0));
        double gz = exp(-(z[iz] - z0) * (z[iz] - z0) / (4 * sigma0 * sigma0));
        double phase = kx0 * x[ix] + ky0 * y[iy] + kz0 * z[iz];
        psi->data[(ix * Ny + iy) * Nz + iz] =
            c_scale(c_new(cos(phase), sin(phase)), gx * gy * gz);
      }
    }
  }

  double scale = 1.0 / sqrt(grid_norm(psi, dx * dy * dz));
  for (int i = 0; i < N; i++) {
    psi->data[i] = c_scale(psi->data[i], scale);
  }

  double T = 2.0, dt = 0.01;
  int fail = soft_evolve_3d(psi, V, Nx, Ny, Nz, dx, dy, dz, dt, (int)(T / dt),
                            hbar, mass) != 0;

  double mean_x = 0.0;
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int iz = 0; iz < Nz; iz++) {
        mean_x +=
            x[ix] * c_abs2(psi->data[(ix * Ny + iy) * Nz + iz]) * dx * dy * dz;
      }
    }
  }

  fail |= check_close(mean_x, x0 + (hbar * kx0 / mass) * T, 1e-3,
                      "3D free particle <x>(T)");
  fail |= check_close(grid_norm(psi, dx * dy * dz), 1.0, 1e-6,
                      "3D free particle norm conservation");

  cvector_free(psi);
  free(V);
  free(x);
  free(y);
  free(z);

  return fail;
}

// Test 3: 2D harmonic-oscillator coherent state (no spreading, exact
// classical centroid motion via Ehrenfest's theorem).
static int test_ho_coherent_state_2d(void) {
  int Nx = 64, Ny = 64;
  double Lx = 20.0, Ly = 20.0;
  double dx = Lx / Nx, dy = Ly / Ny;
  double hbar = 1.0, mass = 1.0, omega = 1.0;
  double sigma0 = sqrt(hbar / (mass * omega));
  double x0 = 2.0, y0 = -1.5;

  double *x = malloc(Nx * sizeof(double));
  double *y = malloc(Ny * sizeof(double));
  for (int i = 0; i < Nx; i++) {
    x[i] = -Lx / 2 + i * dx;
  }
  for (int i = 0; i < Ny; i++) {
    y[i] = -Ly / 2 + i * dy;
  }

  cvector_t *psi = cvector_alloc(Nx * Ny);
  double *V = malloc(Nx * Ny * sizeof(double));
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      V[ix * Ny + iy] =
          0.5 * mass * omega * omega * (x[ix] * x[ix] + y[iy] * y[iy]);
    }
  }

  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      double gx = exp(-(x[ix] - x0) * (x[ix] - x0) / (4 * sigma0 * sigma0));
      double gy = exp(-(y[iy] - y0) * (y[iy] - y0) / (4 * sigma0 * sigma0));
      psi->data[ix * Ny + iy] = c_real(gx * gy);
    }
  }

  double scale = 1.0 / sqrt(grid_norm(psi, dx * dy));
  for (int i = 0; i < Nx * Ny; i++) {
    psi->data[i] = c_scale(psi->data[i], scale);
  }

  double T = 4.0, dt = 0.005;
  int fail = soft_evolve_2d(psi, V, Nx, Ny, dx, dy, dt, (int)(T / dt), hbar,
                            mass) != 0;

  double mean_x = 0.0, mean_y = 0.0;
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      double p = c_abs2(psi->data[ix * Ny + iy]) * dx * dy;
      mean_x += x[ix] * p;
      mean_y += y[iy] * p;
    }
  }

  fail |= check_close(mean_x, x0 * cos(omega * T), 1e-3,
                      "2D HO coherent state <x>(T) (Ehrenfest)");
  fail |= check_close(mean_y, y0 * cos(omega * T), 1e-3,
                      "2D HO coherent state <y>(T) (Ehrenfest)");
  fail |= check_close(grid_norm(psi, dx * dy), 1.0, 1e-6,
                      "2D HO coherent state norm conservation");

  cvector_free(psi);
  free(V);
  free(x);
  free(y);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Free-particle Gaussian wavepacket spreading (2D, V=0):\n");
  failed += test_free_particle_2d();

  printf("Free-particle Gaussian wavepacket spreading (3D, V=0):\n");
  failed += test_free_particle_3d();

  printf("2D harmonic-oscillator coherent state (Ehrenfest, no spreading):\n");
  failed += test_ho_coherent_state_2d();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
