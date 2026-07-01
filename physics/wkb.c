/*
WKB approximation: quantization and tunneling.
 TODO: Bohr-Sommerfeld
*/

#include "wkb.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>

/* Helper compute integral \int \sqrt(2m(E - V(x))) dx between turning points.
   For a given E, find turning points x1, x2 where V(x) = E.
   HACK: function is placeholder
   TODO: additional arguments (x_min, x_max, grid size) or root‑finding routine
*/
static double action_integral(potential_fn V, void *params, double m, double E,
                              double x_min, double x_max, int N) {
  // HACK: Simple Riemann sum over the given range, assuming turning points
  // lie within [x_min, x_max].
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
                        int n, double E_min, double E_max, double tol) {
  /*
     WKB quantization condition:
     \int_{x1}^{x2} \sqrt(2m(E - V(x))) dx = (n + 1/2) \pi \hbar
     Use bisection to find energy E that satisfies.
     HACK: assume turning points lie within a default range.
     TODO: implementation more parameters (x_min, x_max, grid).
  */
  // Stub: returns E_min.  TODO: implement proper bisection with
  // action_integral.
  return E_min;
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
