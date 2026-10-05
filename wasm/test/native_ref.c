/* Native reference run of the playground scenario; compared against the WASM
 * build by wasm/test/check.mjs. Build: see wasm/test/run.sh */
#include <stdio.h>
double *qmc_potential(void);
float *qmc_density(void);
int qmc_init(int, double, double);
void qmc_build_slits(double, double, int, double, double, double);
int qmc_wavepacket(double, double, double, double, double);
int qmc_step(int);
double qmc_norm(void);
int main(void) {
  if (qmc_init(128, 40.0, 0.02) != 0) return 1;
  qmc_build_slits(2.0, 1.0, 2, 6.0, 1.6, 60.0);
  qmc_wavepacket(-10.0, 0.0, 3.0, 0.0, 2.0);
  printf("norm0 %.12f\n", qmc_norm());
  if (qmc_step(300) != 0) return 2;
  printf("norm300 %.12f\n", qmc_norm());
  float *d = qmc_density();
  double s = 0.0, wsum = 0.0;
  for (int i = 0; i < 128 * 128; i++) s += d[i];
  /* weighted sum to catch layout bugs */
  for (int ix = 0; ix < 128; ix++)
    for (int iy = 0; iy < 128; iy++) wsum += d[ix * 128 + iy] * (ix * 3 + iy);
  printf("sum %.9e\nwsum %.9e\n", s, wsum);
  return 0;
}
