/*
 * Transfer matrix solver for 1D resonant tunneling
 *
 * Analytic formula for a single barrier
 * Implements numerical machinery for arbitrary stacks of layers, with
 * double-barrier (Fabry-Perot) geometry as primary use case
 *
 * Units: \hbar^2 / 2m = 1, so energy is E = k^2 and local wave number in a
 * layer of potential V is q = \sqrt(E - V)
 *
 * Implementation: Core primitive is `layer()`, which returns 2x2 real transfer
 * matrix mapping (\psi, \psi') across a single layer of thickness d
 * It handles 3 regimes explicitly:
 *   - E > V:   oscillatory (cos/sin)
 *   - E < V:   evanescent (cosh/sinh)
 *   - E = V:   linear (1, d; 0, 1)  - handled separately to avoid a 1/0
 *              singularity in q or \kappa coefficients
 *
 * NOTE: Because all transfer matrices are real and 2x2, `total()` composes them
 * with 3 multiply-adds per layer and no dynamic allocation. Transmission
 * formula T = 4 / [(M11 + M22)^2 + (k M12 - M21 / k)^2] is then trivial to
 * evaluate
 *
 * Reflection amplitude: `db_solve()` needs complex reflection amplitude r, not
 * just |r|^2, because `db_profile()` uses it to reconstruct wavefunction
 * everywhere (phase of r determines where standing-wave nodes sit in left
 * lead). Transfer matrix itself is real, so r is recovered by matching boundary
 * conditions at x = 0:
 *
 *     \psi(0)   = 1 + r
 *     \psi'(0)  = i k (1 - r)
 *
 * and solving resulting complex linear equation for r
 * WARN: This is a trick that keeps transfer matrix real while still getting
 * full complex S-matrix
 *
 * Wavefunction profile: `db_profile()` evaluates |\psi(x)|^2 on an arbitrary
 *   grid
 *   - For x <= 0 (left lead), it uses exact plane-wave form
 *     \exp(ikx) + r \exp(-ikx)
 *   - For x > 0, it propagates (\psi, \psi') forward through layers by
 *     applying `layer()` one step at a time, stopping at requested x
 *     Points past right edge of structure (x > L) use free propagation from L
 * NOTE: This avoids building an interpolation table and keeps profile exact at
 * every query point
 *
 * Dwell time: `db_dwell_time()` integrates |\psi|^2 over structure using
 * composite Simpson's rule (400 subintervals per layer), then divides by
 * incident flux 2k. Dwell time diverges at a Fabry-Perot resonance as barrier
 * transparency goes to zero; this numerical integration is accurate to Simpson
 * error (O(h^4)) for any smooth profile
 */

#include "double_barrier.h"
#include <math.h>

typedef struct {
  double m11, m12, m21, m22;
} mat2;

static mat2 layer(double v, double d, double energy) {
  double e = energy - v;

  if (fabs(e) < 1e-14) {
    return (mat2){1.0, d, 0.0, 1.0};
  }

  if (e > 0.0) {
    double q = sqrt(e);
    double c = cos(q * d);
    double s = sin(q * d);

    return (mat2){c, s / q, -q * s, c};
  }

  double kap = sqrt(-e);
  double c = cosh(kap * d);
  double s = sinh(kap * d);

  return (mat2){c, s / kap, kap * s, c};
}

static mat2 mul(mat2 a, mat2 b) {
  return (mat2){a.m11 * b.m11 + a.m12 * b.m21, a.m11 * b.m12 + a.m12 * b.m22,
                a.m21 * b.m11 + a.m22 * b.m21, a.m21 * b.m12 + a.m22 * b.m22};
}

static int valid(int n, const double *v, const double *d, double energy) {
  if (n < 1 || n > DB_MAX_LAYERS || !v || !d || !(energy > 0.0) ||
      !isfinite(energy)) {
    return 0;
  }

  for (int i = 0; i < n; i++) {
    if (!isfinite(v[i]) || !(d[i] >= 0.0) || !isfinite(d[i])) {
      return 0;
    }
  }

  return 1;
}

static mat2 total(int n, const double *v, const double *d, double energy) {
  mat2 m = {1.0, 0.0, 0.0, 1.0};

  for (int i = 0; i < n; i++) {
    m = mul(layer(v[i], d[i], energy), m);
  }

  return m;
}

double db_transmission(int n, const double *v, const double *d, double energy) {
  if (!valid(n, v, d, energy)) {
    return NAN;
  }

  mat2 m = total(n, v, d, energy);
  double k = sqrt(energy);
  double a = m.m11 + m.m22;
  double b = k * m.m12 - m.m21 / k;

  return 4.0 / (a * a + b * b);
}

// Complex helpers
typedef struct {
  double re;
  double im;
} cx;

static cx cadd(cx a, cx b) { return (cx){a.re + b.re, a.im + b.im}; }
static cx cmul(cx a, cx b) {
  return (cx){a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
static cx cscale(cx a, double s) { return (cx){a.re * s, a.im * s}; }
static cx cdiv(cx a, cx b) {
  double dd = b.re * b.re + b.im * b.im;

  return (cx){(a.re * b.re + a.im * b.im) / dd,
              (a.im * b.re - a.re * b.im) / dd};
}

// Reflection amplitude r with \psi(0) = 1 + r, \psi'(0) = i k (1 - r)
static int reflect(int n, const double *v, const double *d, double energy,
                   cx *r) {
  mat2 m = total(n, v, d, energy);
  double k = sqrt(energy);

  // (1 + r)(M21 - i k M11) + (1 - r)(i k M22 + k^2 M12) = 0
  cx a = {m.m21, -k * m.m11};
  cx b = {k * k * m.m12, k * m.m22};
  cx num = {a.re + b.re, a.im + b.im};
  cx den = {a.re - b.re, a.im - b.im};

  if (den.re * den.re + den.im * den.im < 1e-300) {
    return -1;
  }

  cx q = cdiv(num, den);

  *r = (cx){-q.re, -q.im};

  return 0;
}

int db_solve(int n, const double *v, const double *d, double energy,
             db_result_t *out) {
  cx r;

  if (!out || !valid(n, v, d, energy) || reflect(n, v, d, energy, &r) != 0) {
    return -1;
  }

  mat2 m = total(n, v, d, energy);
  double k = sqrt(energy);

  // \psi'(0) = i k (1 - r) = i k (1 - r.re - i r.im) = k r.im + i k (1 - r.re)
  cx psi0 = {1.0 + r.re, r.im};
  cx dpsi0 = {k * r.im, k * (1.0 - r.re)};
  cx psil = cadd(cscale(psi0, m.m11), cscale(dpsi0, m.m12));

  out->transmission = psil.re * psil.re + psil.im * psil.im;
  out->reflection = r.re * r.re + r.im * r.im;
  out->r_re = r.re;
  out->r_im = r.im;

  return 0;
}

int db_profile(int n, const double *v, const double *d, double energy,
               double xmin, double xmax, int count, double *psi2) {
  cx r;

  if (!psi2 || count < 2 || !(xmin < xmax) || !isfinite(xmin) ||
      !isfinite(xmax) || !valid(n, v, d, energy) ||
      reflect(n, v, d, energy, &r) != 0) {
    return -1;
  }

  double k = sqrt(energy);
  double len = 0.0;

  for (int i = 0; i < n; i++) {
    len += d[i];
  }

  for (int p = 0; p < count; p++) {
    double x = xmin + (xmax - xmin) * p / (count - 1);
    cx psi;

    if (x <= 0.0) {
      // \exp(ikx) + r \exp(-ikx)
      cx a = {cos(k * x), sin(k * x)};
      cx b = cmul(r, (cx){cos(k * x), -sin(k * x)});

      psi = cadd(a, b);
    } else {
      cx dpsi;

      // propagate (\psi, \psi') from 0 to x
      psi = (cx){1.0 + r.re, r.im};
      dpsi = (cx){k * r.im, k * (1.0 - r.re)};

      double pos = 0.0;

      for (int i = 0; i < n && pos < x; i++) {
        double step = (x - pos < d[i]) ? x - pos : d[i];
        mat2 m = layer(v[i], step, energy);

        cx np = cadd(cscale(psi, m.m11), cscale(dpsi, m.m12));
        cx nd = cadd(cscale(psi, m.m21), cscale(dpsi, m.m22));

        psi = np;
        dpsi = nd;
        pos += step;
      }

      if (x > len) { // free propagation in right lead
        mat2 m = layer(0.0, x - len, energy);

        psi = cadd(cscale(psi, m.m11), cscale(dpsi, m.m12));
      }
    }

    psi2[p] = psi.re * psi.re + psi.im * psi.im;
  }

  return 0;
}

double db_dwell_time(int n, const double *v, const double *d, double energy) {
  double len = 0.0;
  double sum = 0.0;
  enum { PER = 400 };
  double buf[PER + 1];

  if (!valid(n, v, d, energy)) {
    return NAN;
  }

  for (int i = 0; i < n; i++) {
    len += d[i];
  }

  if (len <= 0.0) {
    return 0.0;
  }

  // Composite Simpson on each layer
  double x0 = 0.0;

  for (int i = 0; i < n; i++) {
    if (d[i] > 0.0 &&
        db_profile(n, v, d, energy, x0, x0 + d[i], PER + 1, buf) == 0) {
      double h = d[i] / PER;
      double s = buf[0] + buf[PER];

      for (int p = 1; p < PER; p++) {
        s += buf[p] * (p % 2 ? 4.0 : 2.0);
      }

      sum += s * h / 3.0;
    }

    x0 += d[i];
  }

  return sum / (2.0 * sqrt(energy));
}
