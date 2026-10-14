/*
 * Quantum Zeno effect for a resonantly driven two-level system
 *
 * Let |0> be initial state and Omega resonant Rabi frequency. Without
 * measurement drive rotates Bloch vector about x axis, so survival probability
 * in |0> oscillates as
 *
 *     P0(t) = \cos^2(\Omega t / 2)
 *
 * and a full flip |0> -> |1> occurs at t = pi / Omega
 *
 * Pulsed measurement: Split [0, T] into N intervals of length tau = T / N
 * At end of each interval a projective measurement asks "still in |0>?"
 * Conditional on surviving, qubit is reset to |0> and drive phase restarts
 * One interval therefore contributes
 *
 *     p = \cos^2(\Omega \tau / 2)
 *
 * so final survival probability after N measurements is
 *
 *     P_N(T) = p^N = cos^{2N}(Omega T / (2N))
 *
 * For fixed T and large N this behaves as \exp[-(\Omega T)^2 / (4N)] -> 1:
 * more frequent observation freezes evolution, which is quantum Zeno
 * effect
 *
 * zeno_survival_pulsed() returns this final probability
 * zeno_survival_pulsed_at() returns same product at an arbitrary time t
 * inside sequence:
 *
 *     p^k \cos^2(\Omega (t - k \tau) / 2)
 *
 * Where
 *  - k intervals have already been completed
 *  - \tau = T / N
 *
 * zeno_sample_pulsed() generates independent Monte-Carlo trajectories with
 * exactly per-interval survival probability p. alive[k] is number of
 * trajectories that passed first k measurements, so alive[0] = n_traj and
 * alive[N] is sampled final survivor count
 *
 * Continuous measurement / dephasing: Instead of projective kicks, apparatus
 * dephases qubit at rate \Gamma
 * In rotating frame Bloch equations are
 *
 *     dz/dt = -\Omega y
 *     dy/dt =  \Omega z - \Gamma y
 * with z = P0 - P1
 * Decay of z is governed by slowest root of
 *
 *     \lambda^2 + \Gamma \lambda + \Omega^2 = 0
 *
 * For \Gamma <= 2 \Omega decay rate is \Gamma / 2
 * For \Gamma > 2 Omega it is
 *
 *     (\Gamma - \sqrt(\Gamma^2 - 4 \Omega^2)) / 2  ~  \Omega^2 / \Gamma
 *
 * Decay is fastest at \Gamma = 2 \Omega. Both weaker measurement (poor
 * coupling to apparatus) and stronger measurement (Zeno freezing) slow
 * transition. zeno_continuous_rate() returns this slowest positive decay
 * rate. zeno_survival_continuous() integrates Bloch equations by RK4 from
 * z = 1, y = 0 and returns P0(t) = (1 + z(t)) / 2
 */

#include "zeno.h"
#include <math.h>
#include <stddef.h>

double zeno_survival_pulsed(double omega, double t_total, int n_meas) {
  if (n_meas < 1 || !(t_total >= 0.0) || !isfinite(omega) ||
      !isfinite(t_total)) {
    return NAN;
  }

  double c = cos(omega * t_total / (2.0 * n_meas));

  return pow(c * c, n_meas);
}

double zeno_survival_pulsed_at(double omega, double t_total, int n_meas,
                               double t) {
  if (n_meas < 1 || !(t_total >= 0.0) || !isfinite(omega) ||
      !isfinite(t_total) || !(t >= 0.0) || t > t_total * (1.0 + 1e-12)) {
    return NAN;
  }

  double tau = t_total / n_meas;
  double k = floor(t / tau * (1.0 + 1e-12));

  if (k > n_meas) {
    k = n_meas;
  }

  double c = cos(omega * tau / 2.0), p = c * c;
  double s = cos(omega * (t - k * tau) / 2.0);

  return pow(p, k) * s * s;
}

static void bloch_rhs(double omega, double gamma, const double *s, double *d) {
  // s = (z, y): dz/dt = -\omega y, dy/dt = \omega z - \gamma y
  d[0] = -omega * s[1];
  d[1] = omega * s[0] - gamma * s[1];
}

double zeno_survival_continuous(double omega, double gamma, double t) {
  if (!isfinite(omega) || !(gamma >= 0.0) || !isfinite(gamma) || !(t >= 0.0) ||
      !isfinite(t)) {
    return NAN;
  }

  double rate = fmax(fabs(omega), gamma);
  double steps_f = ceil(t * rate * 40.0) + 200.0;
  int steps = steps_f > 2.0e6 ? 2000000 : (int)steps_f;
  double h = t / steps;
  double s[2] = {1.0, 0.0};
  double k1[2];
  double k2[2];
  double k3[2];
  double k4[2];
  double tmp[2];

  for (int i = 0; i < steps; i++) {
    bloch_rhs(omega, gamma, s, k1);
    tmp[0] = s[0] + 0.5 * h * k1[0];
    tmp[1] = s[1] + 0.5 * h * k1[1];

    bloch_rhs(omega, gamma, tmp, k2);
    tmp[0] = s[0] + 0.5 * h * k2[0];
    tmp[1] = s[1] + 0.5 * h * k2[1];

    bloch_rhs(omega, gamma, tmp, k3);
    tmp[0] = s[0] + h * k3[0];
    tmp[1] = s[1] + h * k3[1];

    bloch_rhs(omega, gamma, tmp, k4);
    s[0] += h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    s[1] += h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
  }

  return 0.5 * (1.0 + s[0]);
}

double zeno_continuous_rate(double omega, double gamma) {
  if (!isfinite(omega) || !(gamma >= 0.0) || !isfinite(gamma)) {
    return NAN;
  }

  double w = fabs(omega);
  if (gamma <= 2.0 * w) {
    return gamma / 2.0;
  }

  // HACK: (\gamma - \sqrt(\gamma^2 - 4 w^2)) / 2, written to avoid cancellation
  // for large \gamma
  return 2.0 * w * w / (gamma + sqrt(gamma * gamma - 4.0 * w * w));
}

static unsigned next(unsigned *s) {
  unsigned x = *s;

  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *s = x;

  return x;
}

int zeno_sample_pulsed(double omega, double t_total, int n_meas, int n_traj,
                       unsigned seed, int *alive) {
  if (alive == NULL || n_meas < 1 || n_meas > 100000 || n_traj < 1 ||
      n_traj > 1000000 || !(t_total >= 0.0) || !isfinite(omega) ||
      !isfinite(t_total)) {
    return -1;
  }

  double c = cos(omega * t_total / (2.0 * n_meas));
  double p = c * c;
  unsigned s = seed ? seed : 0x2545f491u;
  int n_alive = n_traj;

  alive[0] = n_alive;

  // each survivor passes each measurement independently with probability p
  for (int k = 1; k <= n_meas; k++) {
    int next_alive = 0;

    for (int i = 0; i < n_alive; i++) {
      if ((next(&s) >> 8) / 16777216.0 < p) {
        next_alive++;
      }
    }

    n_alive = next_alive;
    alive[k] = n_alive;
  }

  return 0;
}
