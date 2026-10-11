#include <math.h>
#include <stdio.h>

double *qmc_tq_zz_buffer(void);
double *qmc_tq_rate_buffer(void);
double *qmc_tq_mx_buffer(void);
double *qmc_tq_density_buffer(void);
int qmc_tq_set(double h_i, double h_f);
int qmc_tq_eval(double t);
double qmc_tq_last_mx(void);
double qmc_tq_last_rate(void);
int qmc_tq_curve(double t_max, int count);
double qmc_tq_tcrit(int n);
double qmc_tq_front_speed(void);

// Printed to 5 digits so libm last-bit differences cannot show up
int main(void) {
  qmc_tq_set(0.3, 2.0);
  qmc_tq_eval(0.0);
  const double *z = qmc_tq_zz_buffer();
  printf("t=0 mx=%.5f rate=%.5f zz1=%.5f zz10=%.5f\n", qmc_tq_last_mx(),
         qmc_tq_last_rate(), z[1], z[10]);

  qmc_tq_eval(0.8);
  printf("t=0.8 mx=%.5f rate=%.5f zz1=%.5f zz2=%.5f zz5=%.5f\n",
         qmc_tq_last_mx(), qmc_tq_last_rate(), z[1], z[2], z[5]);
  printf("tcrit %.5f %.5f front %.3f\n", qmc_tq_tcrit(0), qmc_tq_tcrit(1),
         qmc_tq_front_speed());

  qmc_tq_curve(3.0, 61);
  const double *r = qmc_tq_rate_buffer(), *m = qmc_tq_mx_buffer();
  printf("curve r[10]=%.5f r[30]=%.5f r[60]=%.5f m[10]=%.5f m[60]=%.5f\n",
         r[10], r[30], r[60], m[10], m[60]);

  qmc_tq_set(3.0, 0.5);
  qmc_tq_eval(3.0);
  printf("para->ordered zz3=%.5f zz4=%.5f zz30_ok=%d tcrit_nan=%d\n", z[3],
         z[4], fabs(z[30]) < 1e-6, isnan(qmc_tq_tcrit(0)) ? 0 : 1);

  const double *d = qmc_tq_density_buffer();
  printf("density d[0]=%.5f d[60]=%.5f\n", d[0], d[60]);

  return qmc_tq_set(3.5, 1.0) == -1 && qmc_tq_curve(1.0, 1000) == -1 &&
                 qmc_tq_eval(-1.0) == -1
             ? 0
             : 1;
}
