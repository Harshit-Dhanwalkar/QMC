/*
 * Sequential spin-1/2 measurements and beam evolution
 *
 * Models a beam of spin-1/2 particles (like silver atoms) passing through a
 * sequence of Stern-Gerlach magnets. State of beam is represented by its
 * intensity I and a Bloch vector r (|r| <= 1).
 *
 * Measurement model:
 * - A projective measurement along axis n projects state onto r = +n or r = -n
 *   Beam intensity is split according to P(+) = (1 + r.n)/2
 * - If both paths are allowed through (mode = SG_PASS_BOTH), spatial
 *   wavepackets can recombine. Perfect recombination (coherence = 1) fully
 *   restores original state. Imperfect recombination (c < 1) dephases
 *   components transverse to n, shrinking Bloch vector:
 *   r -> (r.n)n + c(r - (r.n)n)
 *
 * Screen pattern model:
 * - Physical separation between two spots on detector is d
 * - Width of each spatial wavepacket is \sigma
 * - Quantum prediction is two Gaussians separated by d. Classical prediction
 *   (random magnetic moments) is a uniform distribution stretched by Gaussian
 *   packet width
 *
 * All functions are pure and side-effect free. They perform strict input
 * validation (returning -1 or NaN) to protect against invalid numerical inputs
 * (e.g. negative standard deviations, non-finite values, or unnormalised Bloch
 * vectors).
 */

#include "stern_gerlach.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

sg_vec_t sg_axis(double theta, double phi) {
  sg_vec_t n = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
  return n;
}

static double dot(sg_vec_t a, sg_vec_t b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double sg_prob_up(sg_vec_t r, sg_vec_t n) { return 0.5 * (1.0 + dot(r, n)); }

int sg_chain(sg_vec_t r0, double i0, const sg_analyzer_t *a, int n,
             double intensity[], sg_vec_t r[]) {
  if (n < 0 || !intensity || !r || !isfinite(i0) || i0 < 0.0 ||
      !isfinite(r0.x + r0.y + r0.z) || dot(r0, r0) > 1.0 + 1e-9) {
    return -1;
  }
  if (n > 0 && !a) {
    return -1;
  }

  for (int k = 0; k < n; k++) {
    if (!isfinite(a[k].theta) || !isfinite(a[k].phi) ||
        !(a[k].coherence >= 0.0 && a[k].coherence <= 1.0) || a[k].mode < 0 ||
        a[k].mode > 2) {
      return -1;
    }
  }

  intensity[0] = i0;
  r[0] = r0;

  for (int k = 0; k < n; k++) {
    sg_vec_t n_k = sg_axis(a[k].theta, a[k].phi), cur = r[k];
    double s = dot(cur, n_k);
    double i = intensity[k];

    if (a[k].mode == SG_PASS_UP || a[k].mode == SG_PASS_DOWN) {
      double sign = a[k].mode == SG_PASS_UP ? 1.0 : -1.0;

      intensity[k + 1] = i * 0.5 * (1.0 + sign * s);
      if (intensity[k + 1] > 1e-15) {
        r[k + 1].x = sign * n_k.x;
        r[k + 1].y = sign * n_k.y;
        r[k + 1].z = sign * n_k.z;
      } else {
        intensity[k + 1] = 0.0;
        r[k + 1] = cur;
      }
    } else {
      double c = a[k].coherence;

      intensity[k + 1] = i;
      r[k + 1].x = s * n_k.x + c * (cur.x - s * n_k.x);
      r[k + 1].y = s * n_k.y + c * (cur.y - s * n_k.y);
      r[k + 1].z = s * n_k.z + c * (cur.z - s * n_k.z);
    }
  }

  return 0;
}

double sg_coherence(double d, double sigma) {
  if (!(sigma > 0.0) || !isfinite(d) || !isfinite(sigma)) {
    return NAN;
  }
  return exp(-d * d / (8.0 * sigma * sigma));
}

double sg_gauss(double z, double c, double s) {
  if (!(s > 0.0) || !isfinite(z) || !isfinite(c) || !isfinite(s)) {
    return NAN;
  }

  double u = (z - c) / s;

  return exp(-0.5 * u * u) / (s * sqrt(2.0 * M_PI));
}

double sg_screen_density(double z, double pup, double d, double sigma) {
  if (!(sigma > 0.0) || !isfinite(z) || !isfinite(pup) || !isfinite(d) ||
      !isfinite(sigma)) {
    return NAN;
  }

  return pup * sg_gauss(z, 0.5 * d, sigma) +
         (1.0 - pup) * sg_gauss(z, -0.5 * d, sigma);
}

double sg_classical_density(double z, double d, double sigma) {
  if (!(sigma > 0.0) || !isfinite(z) || !isfinite(d) || !isfinite(sigma)) {
    return NAN;
  }

  if (d <= 0.0) {
    return sg_gauss(z, 0.0, sigma);
  }

  double a = (z + 0.5 * d) / (sigma * M_SQRT2);
  double b = (z - 0.5 * d) / (sigma * M_SQRT2);

  return 0.5 * (erf(a) - erf(b)) / d;
}

double sg_width(double sigma0, double t) {
  if (!(sigma0 > 0.0) || !isfinite(t) || !isfinite(sigma0)) {
    return NAN;
  }

  double u = t / (2.0 * sigma0 * sigma0);

  return sigma0 * sqrt(1.0 + u * u);
}
