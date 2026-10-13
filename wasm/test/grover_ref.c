#include <stdio.h>

double *qmc_gv_amp_buffer(void);
double *qmc_gv_curve_buffer(void);
double *qmc_gv_exact_buffer(void);
int qmc_gv_set_n(int n);
int qmc_gv_mark(int i, int flag);
int qmc_gv_step(void);
int qmc_gv_goto(int k);
int qmc_gv_curve(int count);
double qmc_gv_success(void);
double qmc_gv_theta(void);
double qmc_gv_exact(void);
int qmc_gv_kopt(void);
int qmc_gv_marks(void);
int qmc_gv_iter(void);

// Printed to 5 digits so libm last-bit differences cannot show up
int main(void) {
  qmc_gv_set_n(10);
  qmc_gv_mark(37, 1);
  const double *a = qmc_gv_amp_buffer();
  printf("M=1 theta=%.5f kopt=%d\n", qmc_gv_theta(), qmc_gv_kopt());

  qmc_gv_goto(10);
  printf("k=%d p=%.5f exact=%.5f a37=%.5f a0=%.5f\n", qmc_gv_iter(),
         qmc_gv_success(), qmc_gv_exact(), a[37], a[0]);

  qmc_gv_goto(25);
  printf("k=%d p=%.5f exact=%.5f a37=%.5f a0=%.5f\n", qmc_gv_iter(),
         qmc_gv_success(), qmc_gv_exact(), a[37], a[0]);

  qmc_gv_mark(100, 1);
  qmc_gv_mark(500, 1);
  qmc_gv_mark(900, 1);
  qmc_gv_goto(qmc_gv_kopt());
  printf("M=%d kopt=%d p=%.5f\n", qmc_gv_marks(), qmc_gv_kopt(),
         qmc_gv_success());

  qmc_gv_curve(60);
  const double *c = qmc_gv_curve_buffer(), *e = qmc_gv_exact_buffer();
  printf("curve c5=%.5f c12=%.5f c40=%.5f e12=%.5f\n", c[5], c[12], c[40],
         e[12]);

  qmc_gv_set_n(4);
  qmc_gv_mark(3, 1);
  printf("step_one_mark=%d\n", qmc_gv_step());

  for (int i = 0; i < 8; i++) {
    qmc_gv_mark(i, 1);
}
  printf("marks=%d (capped below N/2)\n", qmc_gv_marks());

  if (qmc_gv_set_n(1) != -1 || qmc_gv_mark(99, 1) != -1) {
    printf("bad args accepted\n");
    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
