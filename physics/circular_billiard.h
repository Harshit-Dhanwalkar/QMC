/*
 * Quantum particle in a circular well ("quantum corral" / circular drum)
 *
 * Unit disc, hbar = 1, mass 1/2, so Schroedinger equation is
 *  -\nabla \psi = k^2 \psi with \psi = 0 on r = 1 and E = k^2
 *
 * The eigenfunctions are:
 *      \psi(r, \theta) = N * J_m(k r) * \cos(m \theta)
 *   OR \psi(r, \theta) = N * J_m(k r) * \sin(m \theta)
 *
 * with k = j_{m,n} n-th zero of Bessel function J_m
 *
 * Every m >= 1 level is doubly degenerate (\cos and \sin); m = 0 is not
 * degenerate
 * A mode (m, n) has (n - 1) radial nodes and 2 m angular nodal lines
 */
#ifndef CIRCULAR_BILLIARD_H
#define CIRCULAR_BILLIARD_H

typedef struct {
  int m;       /* angular quantum number, 0, 1, 2, ... */
  int n;       /* radial quantum number, 1, 2, ... */
  int sine;    /* 0: \cos(m \theta), 1: \sin(m \theta) (only for m >= 1) */
  double k;    /* j_{m,n}; energy is k^2 */
  double norm; /* normalisation constant: integral of \psi^2 over disc is 1 */
} cb_mode_t;

/* Bessel function J_m(x) of integer order 0 <= m <= 200 for x >= 0; NaN
 * otherwise */
double cb_bessel_j(int m, double x);

/* n-th positive zero of J_m (n >= 1), to about 1e-13; NaN on bad input */
double cb_bessel_zero(int m, int n);

/*
 * All modes with k <= k_max, sorted by k (ties: cos before sin, smaller m
 * first) Writes at most max_modes entries and returns how many modes exist
 * (which may exceed max_modes)
 * Returns -1 on bad input
 */
int cb_modes(double k_max, cb_mode_t *modes, int max_modes);

/* psi at radius r in [0, 1] and angle \theta */
double cb_mode_value(const cb_mode_t *mode, double r, double theta);

/* Weyl estimate of number of modes with k <= K: K^2 / 4 - K / 2 */
double cb_weyl_count(double k);

#endif
