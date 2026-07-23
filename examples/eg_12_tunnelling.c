/*
 * Quantum Tunnelling
 *
 * A Gaussian wave packet hits a rectangular potential barrier.
 * Time evolution via Crank-Nicolson (unconditionally stable, norm-preserving).
 *
 * Initial state:             \phi(x,0) = (2 * \pi * \sigma^2)^(-1/4) *
 * \exp(-(x-x_0)^2/4 * \sigma^2) * \exp(i * k_0x)
 * Barrier:                    V(x) = V_0  for  a <= x <= b,  else 0
 * Transmission probability:   T = \int_{x>b} |\phi(x,T)|^2 dx
 * (In natural units \hbar = m = 1)
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../export/plot.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
  printf(" > Quantum Tunnelling (Crank-Nicolson)\n\n");

  // Grid
  int N = 1001;
  double x_min = -10.0, x_max = 10.0;
  double dx = (x_max - x_min) / (N - 1);
  double *x = linspace(x_min, x_max, N);
  if (!x) {
    return 1;
  }

  // Potential: rectangular barrier
  double V0 = 1.5; // barrier height  (\hbar=m=1 units)
  double a = 1.0;  // barrier left edge
  double b = 2.0;  // barrier right edge
  double *V = calloc(N, sizeof *V);
  for (int i = 0; i < N; i++) {
    V[i] = (x[i] >= a && x[i] <= b) ? V0 : 0.0;
  }

  // Initial Gaussian wave packet
  double x0 = -4.0;   // centre
  double sigma = 0.5; // width
  double k0 = 3.0;    // mean momentum (E_kin = k0^2/2)

  cvector_t *psi = cvector_alloc(N);
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    double xi = x[i];
    double env = exp(-(xi - x0) * (xi - x0) / (4.0 * sigma * sigma));
    double pre = pow(2.0 * M_PI * sigma * sigma, -0.25);
    psi->data[i].re = pre * env * cos(k0 * xi);
    psi->data[i].im = pre * env * sin(k0 * xi);
    norm += (psi->data[i].re * psi->data[i].re +
             psi->data[i].im * psi->data[i].im) *
            dx;
  }

  // Normalize
  double inv = 1.0 / sqrt(norm);
  for (int i = 0; i < N; i++) {
    psi->data[i].re *= inv;
    psi->data[i].im *= inv;
  }

  // Build tridiagonal Hamiltonian
  double hbar_sq_2m = 0.5; // \hbar^2/2m with \hbar=m=1
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc(N * sizeof *offdiag);
  build_tridiagonal_hamiltonian(x, V, N, dx, hbar_sq_2m, diag, offdiag);

  // Time evolution parameters
  double dt = 0.005;
  int steps = 800;      // total time T = steps*dt = 4.0
  int snap = steps / 4; // save snapshot every T/4

  printf("   Grid:    N=%d, x∈[%.1f,%.1f], dx=%.4f\n", N, x_min, x_max, dx);
  printf("   Barrier: V_0=%.2f, x∈[%.1f,%.1f]\n", V0, a, b);
  printf("   Packet:  k_0=%.2f, \\sigma=%.2f, x_0=%.2f\n", k0, sigma, x0);
  printf("   E_kin = k_0^2/2 = %.2f  (barrier height = %.2f)\n", k0 * k0 / 2.0,
         V0);
  printf("   Running %d steps, dt=%.4f, T=%.2f...\n\n", steps, dt, steps * dt);

  // Save initial snapshot
  {
    double *prob = malloc(N * sizeof *prob);
    for (int i = 0; i < N; i++) {
      prob[i] =
          psi->data[i].re * psi->data[i].re + psi->data[i].im * psi->data[i].im;
    }
    save_wavefunction("tunnel_psi_t0.dat", x, psi, N);

    plot_opts_t opts = {0};
    opts.title = "|\\psi(x,t=0)|^2";
    opts.xlabel = "x";
    opts.ylabel = "|\\psi|^2";

    // Overlay barrier as second series
    double *Vscaled = malloc(N * sizeof *Vscaled);
    for (int i = 0; i < N; i++) {
      Vscaled[i] = V[i] * 0.3; // scale for vis
    }

    const double *ys[2] = {prob, Vscaled};
    const char *lbs[2] = {"|\\psi|^2", "V(x)"};
    plot_lines("tunnel_t0", PLOT_FORMAT_PNG, x, ys, 2, N, lbs, &opts);

    free(prob);
    free(Vscaled);
  }

  // Evolve
  for (int step = 0; step < steps; step++) {
    crank_nicolson_step(diag, offdiag, dt, psi);

    if ((step + 1) % snap == 0) {
      int t_idx = (step + 1) / snap;
      double t = (step + 1) * dt;

      // Compute norm (should stay = 1)
      double n2 = 0.0;
      for (int i = 0; i < N; i++) {
        n2 += (psi->data[i].re * psi->data[i].re +
               psi->data[i].im * psi->data[i].im) *
              dx;
      }
      printf("   t=%.2f  norm=%.6f\n", t, n2);

      char fname[64], pname[64];
      snprintf(fname, sizeof fname, "tunnel_psi_t%d.dat", t_idx);
      snprintf(pname, sizeof pname, "tunnel_t%d", t_idx);
      save_wavefunction(fname, x, psi, N);

      double *prob = malloc(N * sizeof *prob);
      for (int i = 0; i < N; i++) {
        prob[i] = psi->data[i].re * psi->data[i].re +
                  psi->data[i].im * psi->data[i].im;
      }

      char title[64];
      snprintf(title, sizeof title, "|\\psi(x,t=%.2f)|^2", t);
      plot_opts_t opts = {0};
      opts.title = title;
      opts.xlabel = "x";
      opts.ylabel = "|\\psi|^2";
      plot_line(pname, PLOT_FORMAT_PNG, x, prob, N, &opts);

      free(prob);
    }
  }

  // Transmission coefficient
  double T_coeff = 0.0, R_coeff = 0.0;
  for (int i = 0; i < N; i++) {
    double p = (psi->data[i].re * psi->data[i].re +
                psi->data[i].im * psi->data[i].im) *
               dx;
    if (x[i] > b) {
      T_coeff += p;
    }

    if (x[i] < a) {
      R_coeff += p;
    }
  }
  printf("\n");
  printf("   Transmission T = %.4f\n", T_coeff);
  printf("   Reflection   R = %.4f\n", R_coeff);
  printf("   T + R          = %.4f (should be ~1)\n", T_coeff + R_coeff);

  /* Analytic WKB transmission (thick barrier, E < V0):
   * T_WKB = \exp(-2 * \kappa * L)  where \kappa = \sqrt(2m(V0-E))/\hbar, L=b-a
   */
  double E_kin = k0 * k0 / 2.0;
  if (E_kin < V0) {
    double kappa = sqrt(2.0 * (V0 - E_kin));
    double T_wkb = exp(-2.0 * kappa * (b - a));
    printf("   WKB approx    = %.4f\n", T_wkb);
  }

  cvector_free(psi);
  free(diag);
  free(offdiag);
  free(x);
  free(V);

  return 0;
}
