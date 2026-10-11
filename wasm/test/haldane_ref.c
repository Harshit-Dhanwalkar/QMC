#include <math.h>
#include <stdio.h>

double *qmc_hd_energy_buffer(void);
double *qmc_hd_mean_buffer(void);
double *qmc_hd_prob_buffer(void);
int qmc_hd_set(double t2, double phi, double mass);
int qmc_hd_compute(int n_k);
int qmc_hd_profile(int ik, int n_k, int state);
int qmc_hd_chern(void);
double qmc_hd_limit(void);
double qmc_hd_valley_gap(int valley);
double qmc_hd_vmax(void);
double qmc_hd_cmin(void);

// Printed to 5 digits so libm last-bit differences cannot show up
int main(void) {
  qmc_hd_set(0.1, M_PI / 2, 0.0);
  qmc_hd_compute(25);

  const double *e = qmc_hd_energy_buffer();
  const double *m = qmc_hd_mean_buffer();
  printf("chern=%d limit=%.5f gapK=%.5f gapKp=%.5f vmax=%.5f cmin=%.5f\n",
         qmc_hd_chern(), qmc_hd_limit(), qmc_hd_valley_gap(1),
         qmc_hd_valley_gap(-1), qmc_hd_vmax(), qmc_hd_cmin());

  // k1 = pi is index 12 of 25; states 13 and 14 are zero modes
  printf("k=pi E13=%.5f E14=%.5f\n", fabs(e[12 * 28 + 13]),
         fabs(e[12 * 28 + 14]));
  printf("k=3pi/4 E13=%.5f E14=%.5f m13=%.4f m14=%.4f\n", e[9 * 28 + 13],
         e[9 * 28 + 14], m[9 * 28 + 13], m[9 * 28 + 14]);
  printf("k=0 E0=%.5f E27=%.5f\n", e[0], e[27]);

  qmc_hd_profile(9, 25, 13);
  const double *p = qmc_hd_prob_buffer();
  printf("profile p0=%.5f p1=%.5f p13=%.6f\n", p[0], p[1], p[13]);

  qmc_hd_set(0.1, -M_PI / 2, 1.0);
  printf("trivial chern=%d\n", qmc_hd_chern());

  return qmc_hd_set(0.5, 0.0, 0.0) == -1 && qmc_hd_compute(500) == -1 ? 0 : 1;
}
