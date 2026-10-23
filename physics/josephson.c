/*
 * Numerical integrator for RCSJ junction model
 *
 * State (\phi, v, t) is integrated with classical fixed-step RK4 applied to 2D
 * first-order system
 *
 *     \phi' = v
 *      v'   = i_dc + i_ac \sin(\omega t) - \sin(\phi) - v / q
 *
 * RK4 was chosen over adaptive stepping because RHS is smooth everywhere (no
 * events, no stiffness to track), cost per step is negligible next to settling
 * times used by jj_sweep, and a fixed step makes sweeps reproducible for a
 * given dt
 *
 * NOTE: Integration:
 * - jj_step performs one RK4 step with k1..k4 stages inlined; no Jacobian and
 *   no stage array are needed for a 2D autonomous-ish (time-dependent via ac
 *   drive) system
 * - jj_advance subdivides requested span into n = ceil(time / dt_max)
 *   equal steps, so actual step is dt = time / n <= dt_max and final
 *   step is not a stub.
 * - jj_mean_voltage reports mean winding rate
 *       <\phi'> = (\phi(t_end) - \phi(t_start)) / t_avg
 *   Rather than average of sampled v values. Winding rate is exact for
 *   continuous system (since <v> = <\phi'> and endpoints are primitive of v),
 *   and it avoids bias you get from a Riemann sum of a rapidly oscillating v
 *
 * Sweep heuristics (jj_sweep)
 * I-V curve is traced by carrying state from one bias point to next. 3 choices
 * keep this sane across overdamped / underdamped range:
 *
 *  1. Settling time per point: settle = 20 + 20 * min(q, 20)
 *     - 20*q term matches damping timescale tau = q (from v/q term in ODE:
 *       relaxation rate 1/q)
 *     - Fixed 20 covers residual phase settling near a stable fixed point,
 *       where v decays quickly but \phi still has to stop
 *     - min(q, 20) caps cost, so q = 1e3 costs same 420 time units per point as
 *       q = 20. For q >> 20 up-sweep switch to running state then lags ideal
 *       I_c; that is standard trade when sweeping strongly underdamped
 *       junctions with a fixed budget
 *  2. Initial condition for a down sweep. A down sweep must start on running
 *     branch or first point lands on zero-voltage branch. s.v = 2q is a cheap
 *     kick into running state, adequate for i_start of order 1..3; exact value
 *     is immaterial because settle then relaxes it onto true trajectory. Up
 *     sweep starts at rest, which is physically "just turn on current" case
 *  3. Averaging window. When driven, t_avg is rounded up to a whole number of
 *     drive periods T = 2 \pi / \omega, starting from 60 time units. Skip this
 *     and Shapiro steps appear broadened / shifted because average is taken
 *     over a non-integer number of carrier cycles. t_avg is computed once per
 *     sweep, not per point, so window does not drift along sweep
 *
 * Internal max step handed to jj_mean_voltage from jj_sweep is a fixed 0.02
 * For drive frequencies of order 1e2 that no longer resolves carrier; lower
 * constant (or caller should pre-smooth ac term) before trusting high-\omega
 * sweeps
 */
#include "josephson.h"
#include <math.h>

static void deriv(double phi, double v, double t, const jj_params_t *p,
                  double *dphi, double *dv) {
  *dphi = v;
  *dv = p->i_dc + p->i_ac * sin(p->omega * t) - sin(phi) - v / p->q;
}

static int valid(const jj_params_t *p) {
  return p && p->q > 0.0 && isfinite(p->q) && isfinite(p->i_dc) &&
         isfinite(p->i_ac) && isfinite(p->omega) && p->omega >= 0.0;
}

int jj_step(jj_state_t *s, const jj_params_t *p, double dt) {
  if (!s || !valid(p) || !(dt > 0.0) || !isfinite(dt)) {
    return -1;
  }

  double k1p;
  double k1v;
  double k2p;
  double k2v;
  double k3p;
  double k3v;
  double k4p;
  double k4v;

  deriv(s->phi, s->v, s->t, p, &k1p, &k1v);
  deriv(s->phi + 0.5 * dt * k1p, s->v + 0.5 * dt * k1v, s->t + 0.5 * dt, p,
        &k2p, &k2v);
  deriv(s->phi + 0.5 * dt * k2p, s->v + 0.5 * dt * k2v, s->t + 0.5 * dt, p,
        &k3p, &k3v);
  deriv(s->phi + dt * k3p, s->v + dt * k3v, s->t + dt, p, &k4p, &k4v);

  s->phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p);
  s->v += dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v);
  s->t += dt;

  return 0;
}

int jj_advance(jj_state_t *s, const jj_params_t *p, double time,
               double dt_max) {
  if (!s || !valid(p) || !(time >= 0.0) || !(dt_max > 0.0) || !isfinite(time) ||
      !isfinite(dt_max)) {
    return -1;
  }

  int n = (int)ceil(time / dt_max);
  if (n < 1) {
    return 0;
  }

  double dt = time / n;

  for (int i = 0; i < n; i++) {
    jj_step(s, p, dt);
  }

  return 0;
}

double jj_mean_voltage(jj_state_t *s, const jj_params_t *p, double t_transient,
                       double t_avg, double dt) {
  if (!s || !valid(p) || !(t_transient >= 0.0) || !(t_avg > 0.0) ||
      !(dt > 0.0) || !isfinite(t_transient) || !isfinite(t_avg)) {
    return NAN;
  }

  jj_advance(s, p, t_transient, dt);

  double phi0 = s->phi;

  jj_advance(s, p, t_avg, dt);

  return (s->phi - phi0) / t_avg;
}

int jj_sweep(const jj_params_t *p, double i_start, double i_end, int n,
             double *v) {
  if (!v || n < 2 || !valid(p) || !isfinite(i_start) || !isfinite(i_end)) {
    return -1;
  }

  jj_params_t q = *p;
  jj_state_t s = {0.0, 0.0, 0.0};
  double settle = 20.0 + 20.0 * fmin(p->q, 20.0);

  // a down sweep starts in running state: give it a spin
  if (i_end < i_start) {
    s.v = 2.0 * p->q;
  }

  double t_avg = 60.0;

  if (p->i_ac != 0.0 && p->omega > 0.0) { // a whole number of drive periods
    double period = 2.0 * 3.14159265358979323846 / p->omega;

    t_avg = ceil(60.0 / period) * period;
  }

  for (int k = 0; k < n; k++) {
    q.i_dc = i_start + (i_end - i_start) * k / (n - 1);
    v[k] = jj_mean_voltage(&s, &q, settle, t_avg, 0.02);
  }

  return 0;
}
