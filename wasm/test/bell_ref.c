/* Native reference for wasm/qmc_bell.c; diffed against bell_check.mjs. */
#include <math.h>
#include <stdio.h>

double *qmc_bl_quantum_buffer(void);
double *qmc_bl_local_buffer(void);
int qmc_bl_set(double a, double a2, double b, double b2, double visibility);
double qmc_bl_chsh(void);
double qmc_bl_local_chsh(void);
double qmc_bl_max(void);
double qmc_bl_e(int i);
double qmc_bl_e_local(int i);
int qmc_bl_curve(int count);
double qmc_bl_sample(int n, int seed);
double qmc_bl_est(int i);
double qmc_bl_sigma(void);

/* Printed to 5 digits so libm last-bit differences cannot show up. */
int main(void) {
  const double pi = 3.14159265358979323846;

  qmc_bl_set(0.0, pi / 2, pi / 4, 3 * pi / 4, 1.0);
  printf("S=%.5f local=%.5f max=%.5f e0=%.5f e1=%.5f el0=%.5f\n", qmc_bl_chsh(), qmc_bl_local_chsh(), qmc_bl_max(),
         qmc_bl_e(0), qmc_bl_e(1), qmc_bl_e_local(0));
  qmc_bl_set(0.0, pi / 2, pi / 4, 3 * pi / 4, 0.6);
  printf("V=0.6 S=%.5f max=%.5f\n", qmc_bl_chsh(), qmc_bl_max());
  qmc_bl_set(0.1, 1.2, 0.7, 2.0, 0.9);
  printf("odd S=%.5f local=%.5f\n", qmc_bl_chsh(), qmc_bl_local_chsh());
  qmc_bl_curve(181);
  const double *q = qmc_bl_quantum_buffer(), *l = qmc_bl_local_buffer();
  printf("curve q0=%.5f q40=%.5f q90=%.5f l40=%.5f l60=%.5f\n", q[0], q[40], q[90], l[40], l[60]);
  qmc_bl_set(0.0, pi / 2, pi / 4, 3 * pi / 4, 1.0);
  double s_est = qmc_bl_sample(5000, 42); /* sample first: printf argument order is unspecified */
  printf("sample S=%.5f e0=%.5f sigma=%.5f\n", s_est, qmc_bl_est(0), qmc_bl_sigma());
  if (qmc_bl_set(0, 0, 0, 0, 1.5) != -1 || qmc_bl_curve(1) != -1 || !isnan(qmc_bl_e(7))) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
