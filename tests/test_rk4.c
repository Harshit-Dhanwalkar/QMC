#include "../core/complex.h"
#include "../core/ode/rk4.h"
#include "../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * A single complex phase oscillator :
 *   dy/dt = i * \omega * y => y(t) = y0 * \exp(i * \omega * t)
 */
static void phase_rhs(double t, const cvector_t *y, cvector_t *dydt,
                      void *params) {
  (void)t;
  double omega = *(double *)params;
  complex_t iomega = c_imag(omega);
  dydt->data[0] = c_mul(iomega, y->data[0]);
}

// Integrate from t=0 to t=T with fixed step dt
// Return |error| vs analytic
static double run_to(double T, double dt, double omega) {
  cvector_t *y = cvector_alloc(1);
  y->data[0].re = 1.0;
  y->data[0].im = 0.0;

  int steps = (int)(T / dt + 0.5);
  double t = 0.0;
  for (int i = 0; i < steps; i++) {
    rk4_step(t, dt, y, phase_rhs, &omega);
    t += dt;
  }

  double ana_re = cos(omega * T);
  double ana_im = sin(omega * T);
  double dre = y->data[0].re - ana_re;
  double dim = y->data[0].im - ana_im;
  double err = sqrt(dre * dre + dim * dim);

  cvector_free(y);
  return err;
}

int main(void) {
  printf(" > Testing RK4 integrator (complex phase oscillator)...\n");
  int pass = 1;

  double omega = 2.0;
  double T = 1.0;

  // Accuracy at a moderate step size
  double dt = 0.01;
  double err = run_to(T, dt, omega);
  printf("   dt=%.4f  final-state error = %.3e\n", dt, err);
  if (err > 1e-6) {
    printf("   FAIL: error too large for RK4 at this step size\n");
    pass = 0;
  }

  // 4th-order convergence check: halving dt should shrink the error
  // by roughly 2^4 = 16
  double err_full = run_to(T, 0.05, omega);
  double err_half = run_to(T, 0.025, omega);
  double ratio = (err_half > 0.0) ? err_full / err_half : 0.0;
  printf(
      "   Convergence: err(dt)=%.3e, err(dt/2)=%.3e, ratio=%.2f (expect ~16)\n",
      err_full, err_half, ratio);
  if (ratio < 8.0) {
    printf("   FAIL: convergence order looks lower than 4th order\n");
    pass = 0;
  }

  // Norm conservation: |y| should stay 1
  {
    cvector_t *y = cvector_alloc(1);
    y->data[0].re = 1.0;
    y->data[0].im = 0.0;
    double t = 0.0;
    double dtn = 0.02;
    int steps = (int)(T / dtn + 0.5);

    for (int i = 0; i < steps; i++) {
      rk4_step(t, dtn, y, phase_rhs, &omega);
      t += dtn;
    }

    double norm =
        sqrt(y->data[0].re * y->data[0].re + y->data[0].im * y->data[0].im);
    printf("   |y(T)| = %.6f (expected 1.0)\n", norm);

    if (fabs(norm - 1.0) > 1e-6) {
      printf("   FAIL: norm not conserved\n");
      pass = 0;
    }

    cvector_free(y);
  }

  if (pass) {
    printf("   RK4 test passed.\n");
    return 0;
  }
  printf("FAIL\n");

  return 1;
}
