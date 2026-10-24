/*
 * Resistively and Capacitively Shunted Josephson junction (RCSJ model)
 *
 * Superconducting phase difference \phi of a junction driven by a current
 * obeys equation of a damped, driven pendulum, in units where time is 1 /
 * \omega_p (plasma frequency) and voltage is (\hbar \omega_p / 2e):
 *
 *   \phi'' + \phi' / Q + \sin(\phi) = i_dc + i_ac \sin(\Omega t)
 * and
 *   V = <\phi'>
 *
 * Equivalently a particle with friction 1/Q in tilted washboard potential
 *  U = -\cos(\phi) - i_dc \phi
 * Consequences:
 *  - dc drive below critical current (i < 1): particle sits at \phi = \asin(i)
 *    and V = 0 (supercurrent)
 *  - strongly overdamped (Q -> 0): V = Q \sqrt(i^2 - 1) above i = 1 (that is
 *    \sqrt(i^2 - 1) in units of I_c R)
 *  - Underdamped (Q >> 1): hysteresis, with a switch to voltage state at i = 1
 *    and a retrapping current near 4 / (\pi Q)
 *  - ac drive: phase locks to drive and I-V curve has Shapiro steps at
 *    V = n \Omega
 */
#ifndef JOSEPHSON_H
#define JOSEPHSON_H

typedef struct {
  double q;     /* quality factor, damping 1 / q */
  double i_dc;  /* dc bias current / critical current */
  double i_ac;  /* ac amplitude */
  double omega; /* ac frequency in units of plasma frequency */
} jj_params_t;

typedef struct {
  double phi, v; /* phase and d\phi/dt */
  double t;
} jj_state_t;

/* One RK4 step of length dt; 0 on success, -1 on bad input (q > 0, finite
 * values, dt > 0) */
int jj_step(jj_state_t *s, const jj_params_t *p, double dt);

/* Advance by `time` in steps of at most dt_max
 *
 * Returns 0 or -1. */
int jj_advance(jj_state_t *s, const jj_params_t *p, double time, double dt_max);

/* Time-averaged voltage <\phi'> = (\phi(t_end) - \phi(t_start)) / t_avg after a
 transient, starting from and updatingstate s. NaN on bad input. With an ac
 drive average should cover a whole number of drive periods to be exact */
double jj_mean_voltage(jj_state_t *s, const jj_params_t *p, double t_transient,
                       double t_avg, double dt);

/* I-V sweep: n currents from i_start to i_end, carrying state from one point to
 * next (so up and down sweeps show hysteresis). First point starts at rest for
 * an up sweep and in the running state (\phi' = 2 q) for a down sweep v[k]
 * receives the mean voltage at current k
 *
 * Returns 0, or -1 on bad input (n >= 2) */
int jj_sweep(const jj_params_t *p, double i_start, double i_end, int n,
             double *v);

#endif
