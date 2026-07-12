/*
WKB approximation: quantization and tunneling.
 TODO: Bohr-Sommerfeld
*/

#include "wkb.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>

/* Helper compute integral \int_{x_min}^{x_max} \sqrt(2m(E - V(x))) dx between
   turning points. For a given E, find turning points x1, x2 where V(x) = E.
*/
// NOTE: Simple Riemann sum over given range
static double action_integral(potential_fn V, void *params, double m, double E,
                              double x_min, double x_max, int N) {
  double dx = (x_max - x_min) / (N - 1);
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    double x = x_min + i * dx;
    double Vx = V(x, params);
    if (Vx < E) {
      integral += sqrt(2.0 * m * (E - Vx)) * dx;
    }
  }
  return integral;
}

double wkb_quantization(potential_fn V, void *params, double m, double hbar,
                        int n, double E_min, double E_max, double x_min,
                        double x_max, int N, double tol) {
  if (!V || N < 2 || E_max <= E_min)
    return E_min;

  double target = (n + 0.5) * M_PI * hbar;
  double lo = E_min, hi = E_max;

  const int max_iter = 100;
  for (int iter = 0; iter < max_iter; iter++) {
    double mid = 0.5 * (lo + hi);
    double action = action_integral(V, params, m, mid, x_min, x_max, N);
    if (fabs(action - target) < tol)
      return mid;
    if (action < target)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

double wkb_transmission(potential_fn V, void *params, double m, double hbar,
                        double E, double x1, double x2, int N) {
  if (!V || N < 2)
    return 0.0;
  double dx = (x2 - x1) / (N - 1);
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    double x = x1 + i * dx;
    double Vx = V(x, params);
    if (Vx > E) {
      double kappa = sqrt(2.0 * m * (Vx - E)) / hbar;
      integral += kappa * dx;
    }
  }
  return exp(-2.0 * integral);
}

double wkb_energy_harmonic(int n, double hbar, double omega) {
  double m = 1.0;
  double omega_param = omega;
  double E_max_guess = (n + 3.0) * hbar * omega;
  double x_max = sqrt(2.0 * E_max_guess / (omega * omega)) * 1.2;

  return wkb_quantization(V_harmonic, &omega_param, m, hbar, n, 1e-6,
                          E_max_guess, -x_max, x_max, 2000, 1e-6);
}

double wkb_action_integral_harmonic(double E, double hbar, double omega) {
  (void)hbar;
  double m = 1.0;
  double omega_param = omega;
  double x_turn = sqrt(2.0 * E / (omega * omega));

  double one_way =
      action_integral(V_harmonic, &omega_param, m, E, -x_turn, x_turn, 2000);
  return 2.0 * one_way;
}

double wkb_tunneling_rectangular(double E, double V0, double width,
                                 double hbar_sq_2m) {
  if (hbar_sq_2m <= 0.0 || V0 <= E)
    return 1.0;
  // \kappa = \sqrt(2m(V0-E))/ \hbar = \sqrt((V0-E)/ hbar_sq_2m),
  // where:
  // hbar_sq_2m = \hbar^2/(2m)
  // V(x) = V0 across the barrier,
  double kappa = sqrt((V0 - E) / hbar_sq_2m);
  return exp(-2.0 * kappa * width);
}
