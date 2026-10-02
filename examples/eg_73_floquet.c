/*
 * Floquet theory: quasi-energies of periodically driven two-level systems
 *
 * NOTE: eg_21_driven.c integrates a driven system through time. Floquet theory
 * instead asks for stationary structure of a periodic Hamiltonian
 *   H(t + T) = H(t): solutions \psi(t) = \exp^{-i \epsilon t} u(t)
 * with a T-periodic mode u(t)
 * Quasi-energies \epsilon are only defined modulo \omega = 2 \pi / T, so they
 * are folded into zone [-\omega/2, \omega/2)
 * Two independent solvers are compared here: time-domain one (diagonalize
 * one-period propagator) and Sambe-space one (diagonalize harmonic-extended
 * Hamiltonian)
 *
 * Three textbook effects:
 *  1. Rabi avoided crossing: a circularly polarized drive is exactly solvable
 *     (no RWA); quasi-energies \omega/2 \pm \Omega_R/2 repel, with minimum gap
 *     \Omega at resonance
 *  2. Coherent destruction of tunneling: for a strong high-frequency drive on a
 *     double well, the tunnel splitting becomes \Delta |J_0(A/\omega)| and
 *     vanishes at the zeros of J_0 - particle is frozen in its well
 *  3. Bloch-Siegert shift: for a linearly polarized drive the counter-rotating
 *     term moves the resonance to \omega_0 + \Omega^2 / (4 \omega_0)
 */

#include "../physics/floquet.h"
#include "../physics/variational.h"
#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

static void put(cmatrix_t *m, int i, int j, double re, double im) {
  cmatrix_set(m, i, j, c_new(re, im));
}

static double circ_dist(double a, double b, double omega) {
  double d = fmod(fabs(a - b), omega);

  return fmin(d, omega - d);
}

typedef struct {
  double w0, Om, om;
} circ_t;

static void circular_H(double t, void *p, cmatrix_t *H) {
  const circ_t *q = p;

  put(H, 0, 0, 0.5 * q->w0, 0.0);
  put(H, 1, 1, -0.5 * q->w0, 0.0);
  put(H, 0, 1, 0.5 * q->Om * cos(q->om * t), -0.5 * q->Om * sin(q->om * t));
  put(H, 1, 0, 0.5 * q->Om * cos(q->om * t), 0.5 * q->Om * sin(q->om * t));
}

// Two-level Sambe harmonics (m = -1, 0, +1) for H = diag/offdiag static part +
// diag/offdiag drive amplitudes:
// static = [[sa, sb],[sb, -sa]] and
// \exp^{+-i w t} coefficient is [[da, db],[db, -da]]
static void make_harmonics(cmatrix_t *h[3], double sa, double sb, double da,
                           double db) {
  for (int i = 0; i < 3; i++) {
    h[i] = cmatrix_alloc(2, 2);
  }

  put(h[1], 0, 0, sa, 0.0);
  put(h[1], 1, 1, -sa, 0.0);
  put(h[1], 0, 1, sb, 0.0);
  put(h[1], 1, 0, sb, 0.0);

  for (int i = 0; i < 3; i += 2) {
    put(h[i], 0, 0, da, 0.0);
    put(h[i], 1, 1, -da, 0.0);
    put(h[i], 0, 1, db, 0.0);
    put(h[i], 1, 0, db, 0.0);
  }
}

static void free_harmonics(cmatrix_t *h[3]) {
  for (int i = 0; i < 3; i++) {
    cmatrix_free(h[i]);
  }
}

// quasi-energy gap of linearly driven qubit at drive frequency om
typedef struct {
  double w0, Om;
} bs_t;

static double bs_gap(double om, void *p) {
  const bs_t *q = p;
  cmatrix_t *h[3];

  make_harmonics(h, 0.5 * q->w0, 0.0, 0.0, 0.5 * q->Om);
  floquet_result_t *r = floquet_solve_sambe(2, h, 1, om, 6);
  double gap = INFINITY;
  if (r) {
    gap = circ_dist(r->quasienergies[0], r->quasienergies[1], om);
  }

  floquet_result_free(r);
  free_harmonics(h);

  return gap;
}

static void demo_rabi(void) {
  printf("  === Rabi avoided crossing (circular drive, w0 = 1, \\Omega = 0.3) "
         "===\n");
  printf("     %-8s  %-12s  %-12s  %-12s  %-12s\n", "omega", "time e-",
         "analytic e-", "time e+", "analytic e+");
  const double oms[5] = {0.70, 0.90, 1.00, 1.10, 1.30};
  for (int k = 0; k < 5; k++) {
    circ_t p = {1.0, 0.3, oms[k]};
    floquet_result_t *r = floquet_solve_time(2, circular_H, &p, oms[k], 2000);

    double a = floquet_circular_rabi_quasienergy(p.w0, p.Om, p.om, -1);
    double b = floquet_circular_rabi_quasienergy(p.w0, p.Om, p.om, +1);
    if (a > b) {
      double t = a;
      a = b;
      b = t;
    }

    if (r) {
      printf("     %-8.2f  %-12.8f  %-12.8f  %-12.8f  %-12.8f\n", oms[k],
             r->quasienergies[0], a, r->quasienergies[1], b);
    }

    floquet_result_free(r);
  }

  printf("     (at \\omega = w0 zone-folded gap is exactly \\Omega = 0.3)\n\n");
}

static void demo_cdt(void) {
  double Delta = 0.1;
  double om = 4.0;
  printf("  === Coherent destruction of tunneling (\\Delta = %.2f, \\omega = "
         "%.1f) ===\n",
         Delta, om);
  printf("     %-8s  %-14s  %-14s\n", "A/\\omega", "splitting",
         "\\Delta|J0(A/w)|");
  const double xs[7] = {
      0.0, 0.8, 1.6, 2.404825557695773, 3.2, 4.4, 5.520078110286311};
  for (int k = 0; k < 7; k++) {
    cmatrix_t *h[3];
    make_harmonics(h, 0.0, 0.5 * Delta, 0.25 * xs[k] * om, 0.0);
    floquet_result_t *r = floquet_solve_sambe(2, h, 1, om, 30);
    if (r) {
      printf("     %-8.4f  %-14.8f  %-14.8f\n", xs[k],
             circ_dist(r->quasienergies[0], r->quasienergies[1], om),
             floquet_cdt_splitting(Delta, xs[k] * om, om));
    }

    floquet_result_free(r);
    free_harmonics(h);
  }

  printf("     (tunneling is frozen at A/\\omega = 2.4048 and 5.5201, zeros of "
         "J0)\n\n");
}

static void demo_bloch_siegert(void) {
  printf("  === Bloch-Siegert shift (linear drive, w0 = 1) ===\n");
  printf("     %-8s  %-14s  %-14s  %-10s\n", "\\Omega", "found shift",
         "\\Omega^2/(4w0)", "min gap");
  const double Oms[3] = {0.05, 0.10, 0.15};
  for (int k = 0; k < 3; k++) {
    bs_t q = {1.0, Oms[k]};
    double res = golden_section_minimize(0.99, 1.04, bs_gap, &q, 1e-10);

    printf("     %-8.2f  %-14.6e  %-14.6e  %-10.6f\n", Oms[k], res - q.w0,
           Oms[k] * Oms[k] / (4.0 * q.w0), bs_gap(res, &q));
  }

  printf("     (residual differences are O(\\Omega^4/w0^3) next order)\n");
}

int main(void) {
  printf(" > Floquet quasi-energies of periodically driven qubits\n\n");

  demo_rabi();
  demo_cdt();
  demo_bloch_siegert();

  return 0;
}
