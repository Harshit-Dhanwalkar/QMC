#include <math.h>
#include <stdio.h>

int qmc_zm_curves(double xi, double bmax, int npts);
int qmc_zm_lines(double xi, double b);
double qmc_zm_moment(double xi, double b, int k);
double *qmc_zm_curves_buffer(void);
double *qmc_zm_lines_buffer(void);
double *qmc_zm_mj_buffer(void);

// Printed to 5 digits (tiny values clamped to 0 so a signed zero cannot show
// up)
static void pr(const char *label, double x) {
  printf("%s=%.5f ", label, fabs(x) < 5e-6 ? 0.0 : x);
}

int main(void) {
  qmc_zm_curves(1.0, 5.0, 101);
  const double *c = qmc_zm_curves_buffer(), *mj = qmc_zm_mj_buffer();
  printf("mj ");
  for (int k = 0; k < 8; k++) {
    pr("m", mj[k]);
  }
  printf("\n");

  for (int k = 0; k < 8; k++) {
    printf("curve%d ", k);
    pr("b0", c[k * 101]);
    pr("b1", c[k * 101 + 20]);
    pr("b5", c[k * 101 + 100]);
    printf("\n");
  }

  int n = qmc_zm_lines(1.0, 0.7);
  const double *l = qmc_zm_lines_buffer();
  printf("lines n=%d\n", n);
  for (int j = 0; j < 12; j++) {
    printf("row%d ", j);
    pr("e", l[j * 5]);
    pr("i", l[j * 5 + 1]);
    pr("q", l[j * 5 + 2]);
    printf("\n");
  }

  printf("moment ");
  pr("m0", qmc_zm_moment(1.0, 0.7, 0));
  pr("m4", qmc_zm_moment(1.0, 0.7, 4));
  pr("m5", qmc_zm_moment(1.0, 0.7, 5));
  printf("\n");

  if (qmc_zm_curves(0.0, 5.0, 101) != -1 || qmc_zm_curves(1.0, 5.0, 1) != -1 ||
      qmc_zm_curves(1.0, 5.0, 9999) != -1 || qmc_zm_lines(-1.0, 1.0) != -1 ||
      !isnan(qmc_zm_moment(1.0, 1.0, 7))) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
