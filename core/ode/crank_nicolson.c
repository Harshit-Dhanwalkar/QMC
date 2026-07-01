/*
Crank-Nicolson method for TDSE (complex tridiagonal solver)

TDSE
- For time evolution `i hbar d\phi / dt = H\phi`, Crank-Nicolson is
unconditionally stable and unitary
- This reduces to tridiagonal complex linear system per timestep - solvable with
Thomas algorithm in O(N)
*/

#include "crank_nicolson.h"
#include "../complex.h"
#include "../vector.h"
#include "complex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Complex Thomas algorithm for tridiagonal system.
   a: lower diagonal (size N-1),
   b: diagonal (size N),
   c: upper diagonal (size N-1)
   d: RHS (complex) input, solution overwritten in d.
   Returns 0 on success.
*/
static int thomas_complex(const complex_t *a, const complex_t *b,
                          const complex_t *c, complex_t *d, int N) {
  if (N < 1)
    return -1;
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

int crank_nicolson_step(const double *diag, const double *offdiag, double dt,
                        cvector_t *psi) {
  if (!diag || !offdiag || !psi)
    return -1;
  int N = psi->n;
  if (N < 2)
    return -1;

  // Build complex tridiagonal matrix A = I + i*dt/2*H
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

  for (int i = 0; i < N; i++) {
    b[i] = c_add(c_one(), c_mul(c_imag(1.0), c_real(dt / 2.0 * diag[i])));
    if (i < N - 1) {
      c[i] = c_mul(c_imag(1.0), c_real(dt / 2.0 * offdiag[i]));
      a[i] = c[i]; // symmetric
    }
  }

  // Compute RHS = (I - i*dt/2*H) * psi
  for (int i = 0; i < N; i++) {
    complex_t sum = c_zero();
    complex_t diag_factor =
        c_sub(c_one(), c_mul(c_imag(1.0), c_real(dt / 2.0 * diag[i])));
    sum = c_add(sum, c_mul(diag_factor, psi->data[i]));
    if (i > 0) {
      complex_t off_factor =
          c_mul(c_imag(1.0), c_real(-dt / 2.0 * offdiag[i - 1]));
      sum = c_add(sum, c_mul(off_factor, psi->data[i - 1]));
    }
    if (i < N - 1) {
      complex_t off_factor = c_mul(c_imag(1.0), c_real(-dt / 2.0 * offdiag[i]));
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

  for (int i = 0; i < N; i++)
    psi->data[i] = rhs[i];

  free(a);
  free(b);
  free(c);
  free(rhs);
  return 0;
}

void build_tridiagonal_hamiltonian(const double *x, const double *V, int N,
                                   double dx, double hbar_sq_2m, double *diag,
                                   double *offdiag) {
  double coeff = hbar_sq_2m / (dx * dx);
  for (int i = 0; i < N; i++) {
    diag[i] = 2.0 * coeff + V[i];
    if (i < N - 1)
      offdiag[i] = -coeff;
  }
}
