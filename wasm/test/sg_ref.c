/* Native reference for wasm/qmc_sg.c; diffed against sg_check.mjs. */
#include <math.h>
#include <stdio.h>

int qmc_sg_set(int k, double theta, double phi, int mode, double coherence);
int qmc_sg_run(double src_theta, double src_phi, double pol, int n);
double qmc_sg_screen(double pup, double d, double sigma, int npts, double zmax);
int qmc_sg_shoot(int n, double pup, double d, double sigma, int seed);
double *qmc_sg_intensity_buffer(void);
double *qmc_sg_r_buffer(void);
double *qmc_sg_quantum_buffer(void);
double *qmc_sg_classical_buffer(void);
double *qmc_sg_hits_buffer(void);

/* Printed to 5 digits (tiny values clamped to 0 so a signed zero cannot show up). */
static void pr(const char *label, double x) { printf("%s=%.5f ", label, fabs(x) < 5e-6 ? 0.0 : x); }

int main(void) {
  qmc_sg_set(0, 0.0, 0.0, 0, 1.0);
  qmc_sg_set(1, 1.5707963267948966, 0.0, 0, 1.0);
  qmc_sg_set(2, 0.0, 0.0, 0, 1.0);
  qmc_sg_set(3, 0.8, 0.4, 2, 0.35);
  qmc_sg_run(0.0, 0.0, 1.0, 3);
  const double *in = qmc_sg_intensity_buffer(), *r = qmc_sg_r_buffer();

  printf("zxz ");
  for (int k = 0; k <= 3; k++) {
    pr("i", in[k]);
  }
  printf("\n");
  qmc_sg_run(1.0, 0.5, 0.7, 4);
  printf("mixed ");
  for (int k = 0; k <= 4; k++) {
    pr("i", in[k]);
    pr("rz", r[3 * k + 2]);
  }
  printf("\n");
  double ov = qmc_sg_screen(0.4, 3.0, 0.7, 101, 6.0);
  const double *q = qmc_sg_quantum_buffer(), *c = qmc_sg_classical_buffer();

  printf("screen ");
  pr("ov", ov);
  pr("q20", q[20]);
  pr("q50", q[50]);
  pr("q85", q[85]);
  pr("c50", c[50]);
  pr("c85", c[85]);
  printf("\n");
  int ups = qmc_sg_shoot(500, 0.3, 3.0, 0.7, 7);
  const double *h = qmc_sg_hits_buffer();

  printf("shots ups=%d ", ups);
  pr("z0", h[0]);
  pr("z1", h[2]);
  pr("z499", h[998]);
  printf("\n");
  if (qmc_sg_set(4, 0, 0, 0, 1) != -1 || qmc_sg_set(0, 0, 0, 3, 1) != -1 || qmc_sg_set(0, 0, 0, 0, 2) != -1 || qmc_sg_run(0, 0, 1.5, 1) != -1 ||
      qmc_sg_run(0, 0, 1, 5) != -1 || qmc_sg_screen(0.5, 1, 0.0, 10, 5) != -1.0 || qmc_sg_shoot(99999, 0.5, 1, 1, 1) != -1) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
