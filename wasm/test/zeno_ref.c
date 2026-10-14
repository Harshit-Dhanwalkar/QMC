/* Native reference for wasm/qmc_zeno.c; diffed against zeno_check.mjs. */
#include <math.h>
#include <stdio.h>

double *qmc_zn_time_buffer(void);
double *qmc_zn_n_buffer(void);
double *qmc_zn_rate_buffer(void);
double *qmc_zn_cont_buffer(void);
double *qmc_zn_alive_buffer(void);
double qmc_zn_pulsed(int n);
int qmc_zn_pulsed_curves(int n);
double qmc_zn_continuous(double gamma);
int qmc_zn_rates(double gmax);
double qmc_zn_rate(double gamma);
double qmc_zn_sample(int n, int n_traj, int seed);

/* Printed to 5 digits so libm last-bit differences cannot show up. */
int main(void) {
  printf("pulsed N=1 %.5f N=4 %.5f N=50 %.5f\n", qmc_zn_pulsed(1), qmc_zn_pulsed(4), qmc_zn_pulsed(50));
  qmc_zn_pulsed_curves(8);
  const double *t = qmc_zn_time_buffer(), *nc = qmc_zn_n_buffer();
  printf("time c0=%.5f c100=%.5f c250=%.5f c400=%.5f\n", t[0], t[100], t[250], t[400]);
  printf("ncurve n1=%.5f n8=%.5f n100=%.5f n400=%.5f\n", nc[0], nc[7], nc[99], nc[399]);
  double end = qmc_zn_continuous(5.0);
  const double *cc = qmc_zn_cont_buffer();
  printf("continuous g5 end=%.5f c200=%.5f\n", end, cc[200]);
  qmc_zn_rates(10.0);
  const double *r = qmc_zn_rate_buffer();
  printf("rates r0=%.5f r80=%.5f r400=%.5f peak %.5f\n", r[0], r[80], r[400], qmc_zn_rate(2.0));
  double s = qmc_zn_sample(10, 5000, 7);
  const double *al = qmc_zn_alive_buffer();
  printf("sample final=%.5f a1=%.5f a5=%.5f\n", s, al[1], al[5]);
  if (qmc_zn_pulsed_curves(0) != -1 || qmc_zn_rates(-1.0) != -1 || !isnan(qmc_zn_continuous(-1.0)) || !isnan(qmc_zn_sample(0, 10, 1))) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
