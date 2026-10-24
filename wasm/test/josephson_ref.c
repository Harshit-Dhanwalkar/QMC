#include <math.h>
#include <stdio.h>

int qmc_jj_set(double q, double i_dc, double i_ac, double omega);
int qmc_jj_reset(double phi, double v);
int qmc_jj_advance(double time);
double qmc_jj_phi(void);
double qmc_jj_v(void);
double qmc_jj_time(void);
int qmc_jj_iv(double imax);
double *qmc_jj_up_buffer(void);
double *qmc_jj_down_buffer(void);

// Printed to 5 digits (tiny values clamped to 0 so a signed zero cannot show
// up)
static void pr(const char *label, double x) {
  printf("%s=%.5f ", label, fabs(x) < 5e-6 ? 0.0 : x);
}

int main(void) {
  qmc_jj_set(3.0, 1.4, 0.0, 0.5);
  qmc_jj_reset(0.2, 0.0);
  qmc_jj_advance(10.0);
  qmc_jj_advance(10.0);
  qmc_jj_advance(5.0);
  pr("phi", qmc_jj_phi());
  pr("v", qmc_jj_v());
  pr("t", qmc_jj_time());
  printf("\n");

  qmc_jj_set(3.0, 0.0, 0.0, 0.5);
  qmc_jj_iv(2.0);
  const double *up = qmc_jj_up_buffer();
  const double *dn = qmc_jj_down_buffer();
  printf("hyst ");

  pr("up40", up[40]);
  pr("up100", up[100]);
  pr("up160", up[160]);
  pr("dn40", dn[40]);
  pr("dn80", dn[80]);
  pr("dn160", dn[160]);
  printf("\n");

  qmc_jj_set(0.7, 0.0, 1.2, 0.8);
  qmc_jj_iv(2.2);
  printf("shapiro ");

  pr("up60", up[60]);
  pr("up90", up[90]);
  pr("up120", up[120]);
  pr("dn90", dn[90]);
  printf("\n");

  if (qmc_jj_set(0.0, 0, 0, 0) != -1 || qmc_jj_set(1, 9, 0, 0) != -1 ||
      qmc_jj_iv(0.0) != -1 || qmc_jj_advance(-1.0) != -1 ||
      qmc_jj_reset(NAN, 0) != -1) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
