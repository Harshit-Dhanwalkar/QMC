/*
Runge–Kutta algorithm
Runge-Kutta 4th order integration.

ODE Integration
- For TDSE time evolution
*/

#include "rk4.h"
#include "../vector.h"
#include <stdlib.h>
#include <string.h>

int rk4_step(double t, double dt, cvector_t *y, ode_func_t f, void *params) {
  if (!y || !f)
    return -1;
  int n = y->n;

  cvector_t *k1 = cvector_alloc(n);
  cvector_t *k2 = cvector_alloc(n);
  cvector_t *k3 = cvector_alloc(n);
  cvector_t *k4 = cvector_alloc(n);
  cvector_t *ytmp = cvector_alloc(n);
  if (!k1 || !k2 || !k3 || !k4 || !ytmp) {
    cvector_free(k1);
    cvector_free(k2);
    cvector_free(k3);
    cvector_free(k4);
    cvector_free(ytmp);
    return -1;
  }

  f(t, y, k1, params);
  for (int i = 0; i < n; i++)
    ytmp->data[i] = c_add(y->data[i], c_scale(k1->data[i], dt / 2.0));
  f(t + dt / 2.0, ytmp, k2, params);

  for (int i = 0; i < n; i++)
    ytmp->data[i] = c_add(y->data[i], c_scale(k2->data[i], dt / 2.0));
  f(t + dt / 2.0, ytmp, k3, params);

  for (int i = 0; i < n; i++)
    ytmp->data[i] = c_add(y->data[i], c_scale(k3->data[i], dt));
  f(t + dt, ytmp, k4, params);

  for (int i = 0; i < n; i++) {
    complex_t sum = c_add(k1->data[i], c_scale(k2->data[i], 2.0));
    sum = c_add(sum, c_scale(k3->data[i], 2.0));
    sum = c_add(sum, k4->data[i]);
    y->data[i] = c_add(y->data[i], c_scale(sum, dt / 6.0));
  }

  cvector_free(k1);
  cvector_free(k2);
  cvector_free(k3);
  cvector_free(k4);
  cvector_free(ytmp);
  return 0;
}
