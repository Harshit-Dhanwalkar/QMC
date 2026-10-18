/*
 * Implementation of circular quantum billiard (quantum corral) modes.
 *
 * 1. cb_bessel_j: Uses Miller's downward recurrence because upward recurrence
 *                 is numerically unstable for x < m
 *                 Normalized using identity: J_0 + 2 * \sum(J_{2k}) = 1
 *                 Re-scales intermediate values to prevent floating-point
 *                 overflow
 * 2. cb_bessel_zero: Scans J_m(x) in steps of 0.05 to bracket sign change,
 *                 then refines root using 80 iterations of bisection method
 * 3. cb_modes: Generates all modes up to k_max, calculates their normalization
 *                 constants, and sorts them using qsort
 *
 * References: Numerical Recipes (Press et al.) for Bessel function algorithms
 * and Weyl's Law for eigenvalue asymptotics
 */

#include "circular_billiard.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Bessel Functions
double cb_bessel_j(int m, double x) {
  if (m < 0 || m > 200 || !(x >= 0.0) || !isfinite(x)) {
    return NAN;
  }
  if (x == 0.0) {
    return m == 0 ? 1.0 : 0.0;
  }

  // NOTE: Miller's backward recurrence, normalised with J_0 + 2 \sum J_{2k} = 1
  int base =
          (int)fmax((double)m, x) + 30 + (int)(1.5 * sqrt(fmax((double)m, x))),
      start = base + base % 2;
  double jp1 = 0.0;
  double j = 1e-30;
  double sum = 0.0;
  double result = 0.0;

  for (int k = start; k >= 1; k--) {
    double jm1 = 2.0 * k / x * j - jp1;

    jp1 = j;
    j = jm1;

    // j = J_{k-1} (unnormalised), jp1 = J_k
    if (k == m + 1) {
      result = j;
    }
    if (k - 1 > 0 && (k - 1) % 2 == 0) {
      sum += 2.0 * j;
    }

    if (fabs(j) > 1e250) {
      j *= 1e-250;
      jp1 *= 1e-250;
      sum *= 1e-250;
      result *= 1e-250;
    }
  }

  if (m == 0) {
    result = j;
  }

  sum += j; // J_0

  return result / sum;
}

// Root Finding
double cb_bessel_zero(int m, int n) {
  if (m < 0 || m > 200 || n < 1 || n > 2000) {
    return NAN;
  }

  const double step = 0.05;
  double a = fmax((double)m, 0.0) + 1e-9;
  double fa = cb_bessel_j(m, a);
  int found = 0;

  for (;;) {
    double b = a + step, fb = cb_bessel_j(m, b);

    if (fa * fb < 0.0) {
      if (++found == n) {
        // bisection to full precision
        double lo = a;
        double hi = b;
        double flo = fa;

        for (int it = 0; it < 80; it++) {
          double mid = 0.5 * (lo + hi);
          double fm = cb_bessel_j(m, mid);

          if (flo * fm <= 0.0) {
            hi = mid;
          } else {
            lo = mid;
            flo = fm;
          }
        }

        return 0.5 * (lo + hi);
      }
    }

    a = b;
    fa = fb;
    if (a > 1e4) {
      return NAN;
    }
  }
}

// Mode Generation & Sorting
static int cmp(const void *pa, const void *pb) {
  const cb_mode_t *a = pa;
  const cb_mode_t *b = pb;

  if (a->k < b->k - 1e-12) {
    return -1;
  }
  if (a->k > b->k + 1e-12) {
    return 1;
  }
  if (a->sine != b->sine) {
    return a->sine - b->sine;
  }

  return a->m - b->m;
}

int cb_modes(double k_max, cb_mode_t *modes, int max_modes) {
  if (modes == NULL || max_modes < 1 || !(k_max > 0.0) || k_max > 200.0) {
    return -1;
  }

  int count = 0;
  int cap = 0;
  cb_mode_t *all = NULL;

  for (int m = 0; m <= (int)k_max; m++) {
    for (int n = 1;; n++) {
      double k = cb_bessel_zero(m, n);

      if (!(k <= k_max)) {
        break;
      }

      // Degeneracy loop: m=0 has only 1 state, m>=1 has \cos and \sin
      for (int s = 0; s <= (m >= 1 ? 1 : 0); s++) {
        if (count == cap) {
          cap = cap ? 2 * cap : 256;
          cb_mode_t *grown = realloc(all, (size_t)cap * sizeof *all);

          if (grown == NULL) {
            free(all);

            return -1;
          }

          all = grown;
        }

        double jm1 = cb_bessel_j(m + 1, k);

        all[count].m = m;
        all[count].n = n;
        all[count].sine = s;
        all[count].k = k;

        // NOTE: Normalization factor for unit disc
        all[count].norm =
            1.0 / sqrt(0.5 * jm1 * jm1 * (m == 0 ? 2.0 * M_PI : M_PI));
        count++;
      }
    }
  }

  if (count > 0) {
    qsort(all, (size_t)count, sizeof *all, cmp);

    for (int i = 0; i < count && i < max_modes; i++) {
      modes[i] = all[i];
    }
  }

  free(all);

  return count;
}

double cb_mode_value(const cb_mode_t *mode, double r, double theta) {
  if (mode == NULL || !(r >= 0.0 && r <= 1.0 + 1e-12)) {
    return NAN;
  }

  // NOTE: Evaluate angular part based on sine flag
  double ang = mode->sine ? sin(mode->m * theta) : cos(mode->m * theta);

  return mode->norm * cb_bessel_j(mode->m, mode->k * r) * ang;
}

double cb_weyl_count(double k) {
  // NOTE: Weyl's law for 2D domain:
  // N(k) \approx (Area / 4\pi) * k^2 - (Perimeter / 4\pi) * k
  // For unit circle: Area = pi, Perimeter = 2pi => k^2 / 4 - k / 2
  return k * k / 4.0 - k / 2.0;
}
