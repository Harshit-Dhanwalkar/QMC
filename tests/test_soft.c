/*
 * Test: Split-Operator Fourier Transform (SOFT) method, 2D and 3D
 *
 * 1. Free-particle Gaussian wavepacket (V=0), 2D and 3D: a minimum-uncertainty
 *    Gaussian spreads with closed form
 *      <x>(t)    = x0 + (hbar*k0/m)*t
 *      Var(x)(t) = sigma0^2 * (1 + (hbar*t/(2*m*sigma0^2))^2)
 *    independently in every Cartesian direction.  Isolates momentum-space
 *    kinetic propagator (fft2d/fft3d + phase factors) with potential term
 *    switched off entirely
 * 2. 2D harmonic-oscillator coherent state: a Gaussian of width
 *    sigma0 = sqrt(hbar/(m*omega)), displaced from origin, is an exact
 *    coherent state.  It does not spread, and its centroid follows exact
 *    classical trajectory (Ehrenfest is exact for a harmonic potential):
 *      <x>(t) = x0*cos(omega*t),  <y>(t) = y0*cos(omega*t)
 * 3. Norm conservation (unitarity): grid_norm stays 1 throughout every run
 */

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/soft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

#define CHECK_ALLOC(ptr)                                                       \
  do {                                                                         \
    if (!(ptr)) {                                                              \
      fprintf(stderr, "OOM at %s:%d\n", __FILE__, __LINE__);                   \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);
  if (!(err <= tol)) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

static void test_free_particle_2d(void) {
  printf("  === Test: Free-particle Gaussian spreading (2D, V=0) ===\n");

  int Nx = RUNNING_ON_VALGRIND ? 32 : 64;
  int Ny = RUNNING_ON_VALGRIND ? 32 : 64;
  double Lx = RUNNING_ON_VALGRIND ? 12.0 : 24.0;
  double Ly = RUNNING_ON_VALGRIND ? 12.0 : 24.0;
  double dx = Lx / Nx, dy = Ly / Ny;
  double hbar = 1.0;
  double mass = 1.0;
  double sigma0 = 1.2;
  double x0 = -2.0;
  double y0 = 0.5;
  double kx0 = 0.8;
  double ky0 = -0.3;

  double *x = malloc(Nx * sizeof(double));
  double *y = malloc(Ny * sizeof(double));
  CHECK_ALLOC(x);
  CHECK_ALLOC(y);

  for (int i = 0; i < Nx; i++) {
    x[i] = -Lx / 2 + i * dx;
  }
  for (int i = 0; i < Ny; i++) {
    y[i] = -Ly / 2 + i * dy;
  }

  cvector_t *psi = cvector_alloc(Nx * Ny);
  double *V = calloc(Nx * Ny, sizeof(double));
  CHECK_ALLOC(V);

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

  double T = RUNNING_ON_VALGRIND ? 1.0 : 3.0;
  double dt = 0.01;
  int rc =
      soft_evolve_2d(psi, V, Nx, Ny, dx, dy, dt, (int)(T / dt), hbar, mass);
  check_true(rc == 0, "soft_evolve_2d ran");

  double mean_x = 0.0;
  double mean_y = 0.0;
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      double p = c_abs2(psi->data[ix * Ny + iy]) * dx * dy;

      mean_x += x[ix] * p;
      mean_y += y[iy] * p;
    }
  }

  double var_x = 0.0;
  double var_y = 0.0;
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      double p = c_abs2(psi->data[ix * Ny + iy]) * dx * dy;

      var_x += (x[ix] - mean_x) * (x[ix] - mean_x) * p;
      var_y += (y[iy] - mean_y) * (y[iy] - mean_y) * p;
    }
  }

  double tol = RUNNING_ON_VALGRIND ? 1e-1 : 1e-3;
  check_close(mean_x, x0 + (hbar * kx0 / mass) * T, tol,
              "2D free particle <x>(T)");
  check_close(mean_y, y0 + (hbar * ky0 / mass) * T, tol,
              "2D free particle <y>(T)");
  double var_exact =
      sigma0 * sigma0 * (1.0 + pow(hbar * T / (2 * mass * sigma0 * sigma0), 2));
  check_close(var_x, var_exact, tol * 10, "2D free particle Var(x)(T)");
  check_close(var_y, var_exact, tol * 10, "2D free particle Var(y)(T)");
  check_close(grid_norm(psi, dx * dy), 1.0, 1e-6,
              "2D free particle norm conservation");

  cvector_free(psi);
  free(V);
  free(x);
  free(y);
}

static void test_free_particle_3d(void) {
  printf("  === Test: Free-particle Gaussian spreading (3D, V=0) ===\n");

  int Nx = RUNNING_ON_VALGRIND ? 32 : 64;
  int Ny = RUNNING_ON_VALGRIND ? 32 : 64;
  int Nz = RUNNING_ON_VALGRIND ? 16 : 32;
  double Lx = RUNNING_ON_VALGRIND ? 10.0 : 20.0;
  double Ly = RUNNING_ON_VALGRIND ? 10.0 : 20.0;
  double Lz = RUNNING_ON_VALGRIND ? 10.0 : 20.0;
  double dx = Lx / Nx;
  double dy = Ly / Ny;
  double dz = Lz / Nz;
  double hbar = 1.0;
  double mass = 1.0;
  double sigma0 = 1.3;
  double x0 = -1.5;
  double y0 = 0.5;
  double z0 = 1.0;
  double kx0 = 0.5;
  double ky0 = -0.4;
  double kz0 = 0.3;

  double *x = malloc(Nx * sizeof(double));
  double *y = malloc(Ny * sizeof(double));
  double *z = malloc(Nz * sizeof(double));
  CHECK_ALLOC(x);
  CHECK_ALLOC(y);
  CHECK_ALLOC(z);

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
  CHECK_ALLOC(V);

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

  double T = RUNNING_ON_VALGRIND ? 1.0 : 2.0;
  double dt = 0.01;

  int rc = soft_evolve_3d(psi, V, Nx, Ny, Nz, dx, dy, dz, dt, (int)(T / dt),
                          hbar, mass);
  check_true(rc == 0, "soft_evolve_3d ran");

  double mean_x = 0.0;
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int iz = 0; iz < Nz; iz++) {
        mean_x +=
            x[ix] * c_abs2(psi->data[(ix * Ny + iy) * Nz + iz]) * dx * dy * dz;
      }
    }
  }

  double tol = RUNNING_ON_VALGRIND ? 1e-1 : 1e-3;
  check_close(mean_x, x0 + (hbar * kx0 / mass) * T, tol,
              "3D free particle <x>(T)");
  check_close(grid_norm(psi, dx * dy * dz), 1.0, 1e-6,
              "3D free particle norm conservation");

  cvector_free(psi);
  free(V);
  free(x);
  free(y);
  free(z);
}

static void test_ho_coherent_state_2d(void) {
  printf("  === Test: 2D Harmonic-oscillator coherent state (Ehrenfest, no "
         "spreading) ===\n");

  int Nx = RUNNING_ON_VALGRIND ? 32 : 64;
  int Ny = RUNNING_ON_VALGRIND ? 32 : 64;
  double Lx = RUNNING_ON_VALGRIND ? 10.0 : 20.0;
  double Ly = RUNNING_ON_VALGRIND ? 10.0 : 20.0;
  double dx = Lx / Nx, dy = Ly / Ny;
  double hbar = 1.0;
  double mass = 1.0;
  double omega = 1.0;
  double sigma0 = sqrt(hbar / (mass * omega));
  double x0 = 2.0;
  double y0 = -1.5;

  double *x = malloc(Nx * sizeof(double));
  double *y = malloc(Ny * sizeof(double));
  CHECK_ALLOC(x);
  CHECK_ALLOC(y);

  for (int i = 0; i < Nx; i++) {
    x[i] = -Lx / 2 + i * dx;
  }
  for (int i = 0; i < Ny; i++) {
    y[i] = -Ly / 2 + i * dy;
  }

  cvector_t *psi = cvector_alloc(Nx * Ny);
  double *V = malloc(Nx * Ny * sizeof(double));
  CHECK_ALLOC(V);

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

  double T = RUNNING_ON_VALGRIND ? 1.0 : 4.0;
  double dt = 0.005;
  int rc =
      soft_evolve_2d(psi, V, Nx, Ny, dx, dy, dt, (int)(T / dt), hbar, mass);
  check_true(rc == 0, "soft_evolve_2d ran");

  double mean_x = 0.0;
  double mean_y = 0.0;
  for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
      double p = c_abs2(psi->data[ix * Ny + iy]) * dx * dy;

      mean_x += x[ix] * p;
      mean_y += y[iy] * p;
    }
  }

  double tol = RUNNING_ON_VALGRIND ? 5e-2 : 1e-3;
  check_close(mean_x, x0 * cos(omega * T), tol,
              "2D HO coherent state <x>(T) (Ehrenfest)");
  check_close(mean_y, y0 * cos(omega * T), tol,
              "2D HO coherent state <y>(T) (Ehrenfest)");
  check_close(grid_norm(psi, dx * dy), 1.0, 1e-6,
              "2D HO coherent state norm conservation");

  cvector_free(psi);
  free(V);
  free(x);
  free(y);
}

int main(void) {
  printf(" > Split-operator Fourier method (SOFT) tests\n\n");

  test_free_particle_2d();
  test_free_particle_3d();
  test_ho_coherent_state_2d();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll SOFT checks passed\n");

  return 0;
}
