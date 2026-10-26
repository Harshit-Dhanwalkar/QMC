#include <math.h>
#include <stdio.h>

int qmc_dw_set(double barrier, double a);
double qmc_dw_splitting(void);
double qmc_dw_wkb(void);
int qmc_dw_split_curve(double bmax);
double qmc_dw_launch(double center, double sigma);
double qmc_dw_evolve(double t);
double qmc_dw_mean_x(void);
double *qmc_dw_e_buffer(void);
double *qmc_dw_density_buffer(void);
double *qmc_dw_psi_buffer(void);
double *qmc_dw_split_buffer(void);
double *qmc_dw_wkb_buffer(void);
double *qmc_dw_v_buffer(void);

// Printed to 5 digits (tiny values clamped to 0 so a signed zero cannot show
// up)
static void pr(const char *label, double x) {
  printf("%s=%.5f ", label, fabs(x) < 5e-6 ? 0.0 : x);
}

int main(void) {
  qmc_dw_set(8.0, 1.5);
  const double *e = qmc_dw_e_buffer(), *v = qmc_dw_v_buffer();
  printf("levels ");
  for (int k = 0; k < 6; k++) {
    pr("e", e[k]);
  }
  printf("\n");

  printf("split ");
  pr("d", qmc_dw_splitting());
  pr("wkb", qmc_dw_wkb());
  pr("v0", v[0]);
  pr("v250", v[250]);
  printf("\n");

  double cap = qmc_dw_launch(-1.5, 0.35);
  printf("launch ");
  pr("cap", cap);
  double d = qmc_dw_splitting();
  for (int i = 0; i < 4; i++) {
    double t = i * 0.6 * 3.141592653589793 / d / 3.0;
    double p = qmc_dw_evolve(t);
    pr("P", p);
    pr("x", qmc_dw_mean_x());
  }
  printf("\n");

  const double *dn = qmc_dw_density_buffer();
  qmc_dw_evolve(1.0);
  printf("dens ");
  pr("d100", dn[100]);
  pr("d150", dn[150]);
  pr("d350", dn[350]);
  printf("\n");

  qmc_dw_split_curve(20.0);
  const double *sp = qmc_dw_split_buffer(), *wk = qmc_dw_wkb_buffer();
  printf("curve ");
  pr("s10", sp[10]);
  pr("s30", sp[30]);
  pr("s60", sp[60]);
  pr("w30", wk[30]);
  pr("w60", wk[60]);
  printf("\n");

  if (qmc_dw_set(0.0, 1.5) != -1 || qmc_dw_set(8.0, 9.0) != -1 ||
      qmc_dw_split_curve(1.0) != -1 || qmc_dw_launch(-1.5, 5.0) != -1.0 ||
      qmc_dw_launch(9.0, 0.3) != -1.0) {
    printf("bad args accepted\n");

    return 1;
  }

  printf("bad args rejected\n");

  return 0;
}
