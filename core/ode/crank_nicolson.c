/*
Crank-Nicolson method for TDSE (complex tridiagonal solver), with support
for time-dependent V(x,t) and complex absorbing potentials (CAP).

TDSE
 - For time evolution `i * \hbar d\phi / dt = H\phi`, Crank-Nicolson is
   unconditionally stable and unitary for Hermitian H. A complex absorbing
   potential (CAP) makes H = H0 - i * \Gamma(x) non-Hermitian by design, so norm
   decays in absorbing region instead of being conserved
 - This reduces to a tridiagonal complex linear system per timestep - solvable
   with Thomas algorithm in O(N)
*/

#include "crank_nicolson.h"
#include "../complex.h"
#include "../vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Complex Thomas algorithm for tridiagonal system.
 *  a: lower diagonal (size N-1),
 *  b: diagonal (size N),
 *  c: upper diagonal (size N-1)
 *  d: RHS (complex) input, solution overwritten in d.
 *
 *  Returns 0 on success.
 */
static int thomas_complex(const complex_t *a, const complex_t *b,
                          const complex_t *c, complex_t *d, int N) {
  if (N < 1) {
    return -1;
  }

  complex_t *cp = malloc((N - 1) * sizeof(complex_t));
  complex_t *dp = malloc(N * sizeof(complex_t));
  if (!cp || !dp) {
    free(cp);
    free(dp);

    return -1;
  }

  cp[0] = c_div(c[0], b[0]);
  dp[0] = c_div(d[0], b[0]);
  for (int i = 1; i < N; i++) {
    complex_t denom = c_sub(b[i], c_mul(a[i - 1], cp[i - 1]));
    if (c_abs(denom) < 1e-15) {
      free(cp);
      free(dp);

      return -1;
    }

    if (i < N - 1)
      cp[i] = c_div(c[i], denom);
    dp[i] = c_div(c_sub(d[i], c_mul(a[i - 1], dp[i - 1])), denom);
  }

  d[N - 1] = dp[N - 1];
  for (int i = N - 2; i >= 0; i--) {
    d[i] = c_sub(dp[i], c_mul(cp[i], d[i + 1]));
  }

  free(cp);
  free(dp);

  return 0;
}

int crank_nicolson_step_general(const complex_t *diag, const double *offdiag,
                                double dt, cvector_t *psi) {
  if (!diag || !offdiag || !psi) {
    return -1;
  }

  int N = psi->n;
  if (N < 2) {
    return -1;
  }

  // Build complex tridiagonal matrix A = I + i * dt / 2 * H, H possibly
  // non-Hermitian (CAP) and/or time-dependent
  complex_t *a = malloc((N - 1) * sizeof(complex_t));
  complex_t *b = malloc(N * sizeof(complex_t));
  complex_t *c = malloc((N - 1) * sizeof(complex_t));
  complex_t *rhs = malloc(N * sizeof(complex_t));
  if (!a || !b || !c || !rhs) {
    free(a);
    free(b);
    free(c);
    free(rhs);

    return -1;
  }

  complex_t i_dt2 = c_imag(dt / 2.0);

  for (int i = 0; i < N; i++) {
    // b[i] = 1 + i * dt / 2 * diag[i]
    b[i] = c_add(c_one(), c_mul(i_dt2, diag[i]));
    if (i < N - 1) {
      // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage)
      complex_t off = c_mul(i_dt2, c_real(offdiag[i]));
      c[i] = off;
      a[i] = off; // symmetric kinetic coupling
    }
  }

  // Compute RHS = (I - i * dt / 2 * H) * \psi
  for (int i = 0; i < N; i++) {
    complex_t sum = c_zero();
    complex_t diag_factor = c_sub(c_one(), c_mul(i_dt2, diag[i]));
    sum = c_add(sum, c_mul(diag_factor, psi->data[i]));
    if (i > 0) {
      complex_t off_factor =
          c_scale(c_mul(i_dt2, c_real(offdiag[i - 1])), -1.0);
      sum = c_add(sum, c_mul(off_factor, psi->data[i - 1]));
    }

    if (i < N - 1) {
      complex_t off_factor = c_scale(c_mul(i_dt2, c_real(offdiag[i])), -1.0);
      sum = c_add(sum, c_mul(off_factor, psi->data[i + 1]));
    }

    rhs[i] = sum;
  }

  if (thomas_complex(a, b, c, rhs, N) != 0) {
    free(a);
    free(b);
    free(c);
    free(rhs);

    return -1;
  }

  for (int i = 0; i < N; i++) {
    psi->data[i] = rhs[i];
  }

  free(a);
  free(b);
  free(c);
  free(rhs);

  return 0;
}

int crank_nicolson_step(const double *diag, const double *offdiag, double dt,
                        cvector_t *psi) {
  if (!diag || !offdiag || !psi) {
    return -1;
  }

  int N = psi->n;
  if (N < 2) {
    return -1;
  }

  complex_t *diag_c = malloc(N * sizeof *diag_c);
  if (!diag_c) {
    return -1;
  }

  for (int i = 0; i < N; i++) {
    diag_c[i] = c_real(diag[i]);
  }

  int rc = crank_nicolson_step_general(diag_c, offdiag, dt, psi);
  free(diag_c);

  return rc;
}

void build_tridiagonal_hamiltonian(const double *x, const double *V, int N,
                                   double dx, double hbar_sq_2m, double *diag,
                                   double *offdiag) {
  double coeff = hbar_sq_2m / (dx * dx);
  for (int i = 0; i < N; i++) {
    diag[i] = 2.0 * coeff + V[i];
    if (i < N - 1) {
      offdiag[i] = -coeff;
    }
  }
}

void build_tridiagonal_hamiltonian_time_dependent(
    const double *x, int N, double dx, double hbar_sq_2m, potential_time_fn V,
    void *params, double t, const double *absorb, complex_t *diag_out,
    double *offdiag_out) {
  if (!x || !V || !diag_out || !offdiag_out || N < 2) {
    return;
  }

  double coeff = hbar_sq_2m / (dx * dx);
  for (int i = 0; i < N; i++) {
    double re = 2.0 * coeff + V(x[i], t, params);
    double im = absorb ? -absorb[i] : 0.0;
    diag_out[i] = c_new(re, im);

    if (i < N - 1) {
      offdiag_out[i] = -coeff;
    }
  }
}

void cap_build_monomial(const double *x, int N, double width, double eta,
                        int power, double *absorb_out) {
  if (!x || !absorb_out || N < 2 || width <= 0.0 || eta < 0.0) {
    return;
  }

  double x_left_start = x[0] + width;
  double x_right_start = x[N - 1] - width;

  for (int i = 0; i < N; i++) {
    if (x[i] < x_left_start) {
      double d = (x_left_start - x[i]) / width;
      absorb_out[i] = eta * pow(d, power);
    } else if (x[i] > x_right_start) {
      double d = (x[i] - x_right_start) / width;
      absorb_out[i] = eta * pow(d, power);
    } else {
      absorb_out[i] = 0.0;
    }
  }
}

int crank_nicolson_evolve_time_dependent(const double *x, int N, double dx,
                                         double hbar_sq_2m, potential_time_fn V,
                                         void *params, const double *absorb,
                                         cvector_t *psi, double t0, double dt,
                                         int steps) {
  if (!x || !V || !psi || N < 2 || steps < 0) {
    return -1;
  }

  complex_t *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);

    return -1;
  }

  double t = t0;
  int rc = 0;

  for (int s = 0; s < steps && rc == 0; s++) {
    // Evaluate H at step's midpoint time for 2nd-order accuracy.
    double t_mid = t + dt / 2.0;
    build_tridiagonal_hamiltonian_time_dependent(
        x, N, dx, hbar_sq_2m, V, params, t_mid, absorb, diag, offdiag);
    rc = crank_nicolson_step_general(diag, offdiag, dt, psi);
    t += dt;
  }

  free(diag);
  free(offdiag);

  return rc;
}
