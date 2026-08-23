/*
TDSE + TISE solvers
*/
#include "schrodinger.h"
#include "../core/complex.h"
#include "../core/fft/fft.h"       // for split-step
#include "../core/linalg/linalg.h" // for eigen_generic
#include "../core/matrix.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/ode/numerov.h"
#include "../core/vector.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Matrix diagonalization solver
eigen_t *solve_tise_matrix(double *x, int n, double dx, double hbar_sq_2m,
                           potential_fn V, void *params) {
  if (!x || n < 2 || !V) {
    return NULL;
  }
  double coeff = hbar_sq_2m / (dx * dx); // \hbar^2/(2m) / dx^2

  cmatrix_t *H = cmatrix_alloc(n, n);
  if (!H) {
    return NULL;
  }

  // Fill Hamiltonian: H = -\hbar^2/(2m) d^2/dx^2 + V(x)
  for (int i = 0; i < n; i++) {
    double V_i = V(x[i], params);
    CMAT(H, i, i) = c_real(2.0 * coeff + V_i);

    if (i > 0) {
      CMAT(H, i, i - 1) = c_real(-coeff);
    }
    if (i < n - 1) {
      CMAT(H, i, i + 1) = c_real(-coeff);
    }
  }

  // Solve eigenproblem
  eigen_t *eig = cmatrix_eigh_generic(H);
  cmatrix_free(H);

  return eig;
}

// Shooting solver (Numerov)
numerov_solution_t *solve_tise_shoot(numerov_params_t *params, double E_guess,
                                     double E_tol) {
  return numerov_shoot(params, E_guess, E_tol);
}

// Crank-Nicolson time evolution
int evolve_tdse_crank(const double *diag, const double *offdiag, int n,
                      cvector_t *psi, double dt, int steps) {
  if (!diag || !offdiag || !psi || n < 2 || steps < 1) {
    return -1;
  }

  for (int s = 0; s < steps; s++) {
    if (crank_nicolson_step(diag, offdiag, dt, psi) != 0) {
      return -1;
    }
  }

  return 0;
}

// Split-step Fourier
int evolve_tdse_split_step(cvector_t *psi, const double *x, const double *V,
                           int n, double dx, double dt, int steps, double hbar,
                           double mass) {
  if (!psi || !x || !V || n < 2 || steps < 1 || hbar <= 0.0 || mass <= 0.0) {
    return -1;
  }

  double dk = 2.0 * M_PI / (n * dx);
  double hbar_over_2m = hbar / (2.0 * mass);

  for (int s = 0; s < steps; s++) {
    // 1. Half-step in position space: \exp(-i * V * dt / (2 * \hbar))
    for (int i = 0; i < n; i++) {
      double phase = -V[i] * dt / (2.0 * hbar);
      complex_t exp_factor = c_exp(c_imag(phase));

      psi->data[i] = c_mul(psi->data[i], exp_factor);
    }

    // 2. Full step in momentum space: \exp(-i * \hbar k^2 / (2 * m) dt)
    fft(psi); // forward FFT (no normalization)
    for (int i = 0; i < n; i++) {
      double k = (i < n / 2) ? i * dk : (i - n) * dk;
      double phase = -hbar_over_2m * k * k * dt;
      complex_t exp_factor = c_exp(c_imag(phase));

      psi->data[i] = c_mul(psi->data[i], exp_factor);
    }
    ifft(psi); // inverse FFT (with normalization)

    // 3. Half-step in position space: \exp(-i * V * dt / (2 * \hbar))
    for (int i = 0; i < n; i++) {
      double phase = -V[i] * dt / (2.0 * hbar);
      complex_t exp_factor = c_exp(c_imag(phase));

      psi->data[i] = c_mul(psi->data[i], exp_factor);
    }
  }

  return 0;
}
