/*
 * Split-Operator Fourier Transform (SOFT) Method - 2D and 3D
 *
 * schrodinger.c's evolve_tdse_split_step handles 1D. Real wavepacket
 * dynamics, scattering off a 2D barrier, a particle in a 2D/3D trap
 * needs more than one dimension.
 */

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/soft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Split-Operator Fourier Transform: 2D and 3D Wavepackets\n\n");

  // 1. Free-particle Gaussian: spreads, drifts at constant velocity
  printf("   Free-particle 2D Gaussian wavepacket (V=0):\n");
  printf("   %6s  %10s  %10s  %10s  %10s\n", "t", "<x>", "<x>_exact", "Var(x)",
         "Var(x)_exact");
  {
    int Nx = 64, Ny = 64;
    double Lx = 24.0, Ly = 24.0, dx = Lx / Nx, dy = Ly / Ny;
    double hbar = 1.0, mass = 1.0, sigma0 = 1.2;
    double x0 = -2.0, y0 = 0.5, kx0 = 0.8, ky0 = -0.3;

    double *x = malloc(Nx * sizeof(double));
    double *y = malloc(Ny * sizeof(double));
    for (int i = 0; i < Nx; i++) {
      x[i] = -Lx / 2 + i * dx;
    }
    for (int i = 0; i < Ny; i++) {
      y[i] = -Ly / 2 + i * dy;
    }

    cvector_t *psi = cvector_alloc(Nx * Ny);
    double *V = calloc(Nx * Ny, sizeof(double));
    for (int ix = 0; ix < Nx; ix++) {
      for (int iy = 0; iy < Ny; iy++) {
        double gx = exp(-(x[ix] - x0) * (x[ix] - x0) / (4 * sigma0 * sigma0));
        double gy = exp(-(y[iy] - y0) * (y[iy] - y0) / (4 * sigma0 * sigma0));
        double phase = kx0 * x[ix] + ky0 * y[iy];
        psi->data[ix * Ny + iy] =
            c_scale(c_new(cos(phase), sin(phase)), gx * gy);
      }
    }

    double scale = 1.0 / sqrt(grid_norm(psi, dx * dy));
    for (int i = 0; i < Nx * Ny; i++) {
      psi->data[i] = c_scale(psi->data[i], scale);
    }

    double dt = 0.01;
    for (int checkpoint = 0; checkpoint <= 4; checkpoint++) {
      double t = checkpoint * 0.75;
      if (checkpoint > 0) {
        soft_evolve_2d(psi, V, Nx, Ny, dx, dy, dt, (int)(0.75 / dt), hbar,
                       mass);
      }

      double mean_x = 0.0, var_x = 0.0;
      for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
          mean_x += x[ix] * c_abs2(psi->data[ix * Ny + iy]) * dx * dy;
        }
      }

      for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
          var_x += (x[ix] - mean_x) * (x[ix] - mean_x) *
                   c_abs2(psi->data[ix * Ny + iy]) * dx * dy;
        }
      }
      printf("   %6.2f  %10.4f  %10.4f  %10.4f  %10.4f\n", t, mean_x,
             x0 + hbar * kx0 / mass * t, var_x,
             sigma0 * sigma0 *
                 (1.0 + pow(hbar * t / (2 * mass * sigma0 * sigma0), 2)));
    }

    cvector_free(psi);
    free(V);
    free(x);
    free(y);
  }
  printf("\n");

  // 2. HO coherent state: no spreading, exact classical orbit
  printf("   2D harmonic-oscillator coherent state (exact Ehrenfest orbit, "
         "no spreading):\n");
  printf("   %6s  %10s  %10s  %10s  %10s\n", "t", "<x>",
         "x_0 * \\cos(\\omega * t)", "<y>", "y_0 * \\cos(\\omega * t)");
  {
    int Nx = 64, Ny = 64;
    double Lx = 20.0, Ly = 20.0, dx = Lx / Nx, dy = Ly / Ny;
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

    double dt = 0.005;
    for (int checkpoint = 0; checkpoint <= 4; checkpoint++) {
      double t = checkpoint * 1.0;
      if (checkpoint > 0) {
        soft_evolve_2d(psi, V, Nx, Ny, dx, dy, dt, (int)(1.0 / dt), hbar, mass);
      }

      double mean_x = 0.0, mean_y = 0.0;
      for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
          double p = c_abs2(psi->data[ix * Ny + iy]) * dx * dy;
          mean_x += x[ix] * p;
          mean_y += y[iy] * p;
        }
      }

      printf("   %6.2f  %10.4f  %10.4f  %10.4f  %10.4f\n", t, mean_x,
             x0 * cos(omega * t), mean_y, y0 * cos(omega * t));
    }

    printf("     (final norm: %.8f : should stay 1.0 throughout)\n",
           grid_norm(psi, dx * dy));
    cvector_free(psi);
    free(V);
    free(x);
    free(y);
  }

  return 0;
}
