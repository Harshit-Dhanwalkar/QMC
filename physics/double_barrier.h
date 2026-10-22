/*
 * Transmission through a stack of square barriers (transfer matrices), for
 * resonant tunnelling
 *
 * Units \hbar^2 / 2m = 1, so E = k^2. A structure is n layers of potential v[i]
 * and thickness d[i] laid out from x = 0 to L = sum d; leads x < 0 and x > L
 * have V = 0. For a particle incident from left with unit amplitude transfer
 * matrix of (\psi, \psi') through a layer of local wave number q = sqrt(E - V)
 * is
 * [[\cos(qd), \sin(qd) / q], [-q \sin(qd), \cos(qd)]]
 * (\cosh/\sinh for E < V, [[1, d], [0, 1]] at E = V), a real unimodular matrix
 * and
 *
 *     T = 4 / [ (M11 + M22)^2 + (k M12 - M21 / k)^2 ]
 *
 * A single barrier gives T = 1 / (1 + V^2 \sinh^2(\kappa a) / (4 E (V - E)));
 * two identical barriers around a well show Fabry-Perot resonances with T = 1
 * at quasi-bound levels (n \pi / w)^2 of well, whose width shrinks as barriers
 * get opaque and whose dwell time grows as inverse width
 */
#ifndef DOUBLE_BARRIER_H
#define DOUBLE_BARRIER_H

enum { DB_MAX_LAYERS = 16 };

typedef struct {
  double transmission, reflection; /* |t|^2, |r|^2 */
  double r_re, r_im;               /* reflection amplitude at x = 0 */
} db_result_t;

/* Transmission probability by closed transfer-matrix formula; NaN unless E
 * > 0, 1 <= n <= DB_MAX_LAYERS and d >= 0 */
double db_transmission(int n_layers, const double *v, const double *d,
                       double energy);

/* Complex reflection and transmission by matching plane waves; 0 on success, -1
 * on bad input */
int db_solve(int n_layers, const double *v, const double *d, double energy,
             db_result_t *out);

/* |\psi(x)|^2 for unit incident amplitude on `count` points from xmin to xmax
 * (leads included; structure spans 0 .. L). Returns 0, or -1 on bad input
 * (count >= 2, xmin < xmax) */
int db_profile(int n_layers, const double *v, const double *d, double energy,
               double xmin, double xmax, int count, double *psi2);

/* Dwell time: integral of |\psi|^2 over structure divided by incident flux 2k;
 * NaN on bad input */
double db_dwell_time(int n_layers, const double *v, const double *d,
                     double energy);

#endif
