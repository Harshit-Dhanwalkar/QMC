/*
Generic Cartesian geometry optimizer
*/

#include "geometry_optimization.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define GEOM_OPT_MAX_BACKTRACK 40
#define GEOM_OPT_MIN_STEP 1e-14

static double max_abs_component(const double *v, int n) {
  double m = 0.0;
  for (int i = 0; i < n; i++) {
    double a = fabs(v[i]);

    if (a > m) {
      m = a;
    }
  }

  return m;
}

geometry_optimization_result_t *optimize_geometry_steepest_descent(
    geom_energy_grad_fn func, void *user_data, const double *coords0,
    int n_atoms, double step0, double grad_tol, int max_iter) {
  if (!func || !coords0 || n_atoms <= 0 || step0 <= 0.0 || grad_tol <= 0.0 ||
      max_iter <= 0) {
    return NULL;
  }

  int ndim = n_atoms * 3;

  geometry_optimization_result_t *res =
      malloc(sizeof(geometry_optimization_result_t));
  if (!res) {
    return NULL;
  }

  double *coords = malloc(sizeof(double) * ndim);
  double *grad = malloc(sizeof(double) * ndim);
  double *trial_coords = malloc(sizeof(double) * ndim);
  double *trial_grad = malloc(sizeof(double) * ndim);
  if (!coords || !grad || !trial_coords || !trial_grad) {
    free(coords);
    free(grad);
    free(trial_coords);
    free(trial_grad);
    free(res);

    return NULL;
  }

  memcpy(coords, coords0, sizeof(double) * ndim);

  double energy;
  double grad_norm = 0.0;
  int converged = 0;
  int iter = 0;

  int rc = func(coords, n_atoms, user_data, &energy, grad);
  if (rc != 0) {
    free(coords);
    free(grad);
    free(trial_coords);
    free(trial_grad);
    free(res);

    return NULL;
  }

  grad_norm = max_abs_component(grad, ndim);

  double step = step0;

  for (iter = 0; iter < max_iter; iter++) {
    if (grad_norm < grad_tol) {
      converged = 1;

      break;
    }

    int accepted = 0;
    double trial_energy = 0.0;

    for (int bt = 0; bt < GEOM_OPT_MAX_BACKTRACK; bt++) {
      for (int i = 0; i < ndim; i++) {
        // move downhill along force = -gradient
        trial_coords[i] = coords[i] - step * grad[i];
      }

      int rc2 =
          func(trial_coords, n_atoms, user_data, &trial_energy, trial_grad);

      if (rc2 == 0 && trial_energy < energy) {
        accepted = 1;

        break;
      }

      step *= 0.5;
      if (step < GEOM_OPT_MIN_STEP) {
        break;
      }
    }

    if (!accepted) {
      // NOTE: Can't find a downhill step even after backtracking - stuck
      break;
    }

    memcpy(coords, trial_coords, sizeof(double) * ndim);
    memcpy(grad, trial_grad, sizeof(double) * ndim);
    energy = trial_energy;
    grad_norm = max_abs_component(grad, ndim);

    /* NOTE: Grow step back a little after a successful move, so a step shrunk
     * by backtracking earlier doesn't stay pessimistically small for rest of
     * the optimization */
    step *= 1.2;
  }

  res->coords = coords;
  res->energy = energy;
  res->grad_norm = grad_norm;
  res->iterations = iter;
  res->converged = converged;

  free(grad);
  free(trial_coords);
  free(trial_grad);

  return res;
}

void geometry_optimization_result_free(geometry_optimization_result_t *res) {
  if (!res) {
    return;
  }

  free(res->coords);
  free(res);
}
