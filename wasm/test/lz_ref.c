/* Native reference for wasm/qmc_lz.c; diffed against lz_check.mjs. */
#include <math.h>
#include <stdio.h>

double *qmc_lz_scan_buffer(void);
double qmc_lz_time(void);
double qmc_lz_total_time(void);
double qmc_lz_delta(void);
double qmc_lz_upper(void);
double qmc_lz_exact(void);
int qmc_lz_reset(double omega, double rate, double amp, int passes);
int qmc_lz_step(double duration);
int qmc_lz_scan_rate(double omega, double amp, int n, double rmin, double rmax);
int qmc_lz_scan_amp(double omega, double rate, int n, double a0, double a1);

int main(void) {
  qmc_lz_reset(1.0, 1.0, 12.0, 1);
  int rc = 0;

  while (rc == 0) {
    rc = qmc_lz_step(0.5);
  }
  printf("one pass t=%.3f delta=%.6f upper=%.9f exact=%.9f rc=%d\n",
         qmc_lz_time(), qmc_lz_delta(), qmc_lz_upper(), qmc_lz_exact(), rc);
  qmc_lz_reset(1.0, 1.0, 12.0, 2);
  rc = 0;
  while (rc == 0) {
    rc = qmc_lz_step(0.7);
  }
  printf("two pass t=%.3f delta=%.6f upper=%.9f\n", qmc_lz_time(),
         qmc_lz_delta(), qmc_lz_upper());
  qmc_lz_scan_rate(1.0, 10.0, 5, 0.1, 10.0);
  const double *s = qmc_lz_scan_buffer();
  printf("scan rate %.6f %.6f %.6f %.6f %.6f\n", s[0], s[1], s[2], s[3], s[4]);
  qmc_lz_scan_amp(1.0, 1.0, 4, 10.0, 11.0);
  printf("scan amp %.6f %.6f %.6f %.6f\n", s[0], s[1], s[2], s[3]);

  return qmc_lz_reset(5.0, 1.0, 10.0, 1) == -1 &&
                 qmc_lz_scan_rate(1.0, 10.0, 1, 0.1, 1.0) == -1
             ? 0
             : 1;
}
