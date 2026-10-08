#include "landau_zener.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double lz_probability(double omega, double rate) {
  if (!isfinite(omega) || !isfinite(rate) || rate == 0.0) {
    return NAN;
  }

  return exp(-M_PI * omega * omega / (2.0 * fabs(rate)));
}

double lz_upper_population(const double v[3], double omega, double delta) {
  if (!v || !isfinite(omega) || !isfinite(delta)) {
    return NAN;
  }

  double b = sqrt(omega * omega + delta * delta);
  if (!(b > 0.0)) {
    return NAN;
  }

  return 0.5 * (1.0 + (omega * v[0] + delta * v[2]) / b);
}

int lz_lower_state(double v[3], double omega, double delta) {
  if (!v || !isfinite(omega) || !isfinite(delta)) {
    return -1;
  }

  double b = sqrt(omega * omega + delta * delta);
  if (!(b > 0.0)) {
    return -1;
  }

  v[0] = -omega / b;
  v[1] = 0.0;
  v[2] = -delta / b;

  return 0;
}

// dv/dt = b x v with b = (\omega, 0, \delta)
static void deriv(const double v[3], double omega, double delta, double d[3]) {
  d[0] = -delta * v[1];               /* b_y v_z - b_z v_y */
  d[1] = delta * v[0] - omega * v[2]; /* b_z v_x - b_x v_z */
  d[2] = omega * v[1];                /* b_x v_y - b_y v_x */
}

int lz_advance(double v[3], double omega, double delta0, double rate,
               double duration, double dt) {
  if (!v || !isfinite(v[0]) || !isfinite(v[1]) || !isfinite(v[2]) ||
      sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) > 1.0 + 1e-9 ||
      !isfinite(omega) || omega < 0.0 || !isfinite(delta0) || !isfinite(rate) ||
      !isfinite(duration) || duration < 0.0 || !isfinite(dt) || !(dt > 0.0)) {
    return -1;
  }

  double d1 = delta0 + rate * duration;
  double dmax = fmax(fabs(delta0), fabs(d1));
  // RK4 is not norm preserving: keep |b| h small even at the sweep ends
  double hmax = fmin(dt, 0.03 / sqrt(omega * omega + dmax * dmax));
  double nsteps = ceil(duration / hmax);

  if (!(nsteps <= 1e8)) {
    return -1;
  }
  if (nsteps < 1.0) {
    return 0;
  }

  int n = (int)nsteps;
  double h = duration / n;
  double w[3] = {v[0], v[1], v[2]};

  for (int i = 0; i < n; i++) {
    double dl0 = delta0 + rate * (i * h);
    double dlm = dl0 + rate * 0.5 * h;
    double dl1 = dl0 + rate * h;
    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    double y[3];

    deriv(w, omega, dl0, k1);
    for (int c = 0; c < 3; c++) {
      y[c] = w[c] + 0.5 * h * k1[c];
    }

    deriv(y, omega, dlm, k2);
    for (int c = 0; c < 3; c++) {
      y[c] = w[c] + 0.5 * h * k2[c];
    }

    deriv(y, omega, dlm, k3);
    for (int c = 0; c < 3; c++) {
      y[c] = w[c] + h * k3[c];
    }

    deriv(y, omega, dl1, k4);
    for (int c = 0; c < 3; c++) {
      w[c] += h / 6.0 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
    }
  }

  v[0] = w[0];
  v[1] = w[1];
  v[2] = w[2];

  return 0;
}

int lz_sweep(double v[3], double omega, double rate, double amp, double dt,
             int passes) {
  if (!isfinite(rate) || rate == 0.0 || !isfinite(amp) || !(amp > 0.0) ||
      passes < 1) {
    return -1;
  }

  double speed = fabs(rate);
  double duration = 2.0 * amp / speed;
  double w[3];

  if (!v) {
    return -1;
  }

  for (int c = 0; c < 3; c++) {
    w[c] = v[c];
  }

  for (int p = 0; p < passes; p++) {
    double dir = (p % 2 == 0) ? 1.0 : -1.0;

    if (lz_advance(w, omega, -dir * amp, dir * speed, duration, dt) != 0) {
      return -1;
    }
  }

  for (int c = 0; c < 3; c++) {
    v[c] = w[c];
  }

  return 0;
}
