#include <math.h>
#include <stdio.h>

int qmc_dc_klein(double k0, double v0);
int qmc_dc_zitter(void);
int qmc_dc_step(int n);
double qmc_dc_norm(void);
double qmc_dc_position(void);
double qmc_dc_right(void);
double qmc_dc_exact(double k0, double v0);

int main(void) {
  int fails = 0;
  const double heights[] = {0.0, 0.5, 4.0, 10.0};

  for (int h = 0; h < 4; h++) {
    qmc_dc_klein(sqrt(3.0), heights[h]);
    qmc_dc_step(2500);
    double t = qmc_dc_right();
    double ex = qmc_dc_exact(sqrt(3.0), heights[h]);

    printf("klein V0=%.1f norm=%.10f T=%.10f\n", heights[h], qmc_dc_norm(), t);
    if (fabs(t - ex) > 0.01 || fabs(qmc_dc_norm() - 1.0) > 1e-9) {
      fprintf(stderr, "FAIL klein V0=%.1f T=%.4f exact=%.4f\n", heights[h], t,
              ex);
      fails++;
    }
  }

  qmc_dc_zitter();
  qmc_dc_step(39); /* t = 1.56 ~ \pi/2: <x> = -(1 - \cos 2t)/2 ~ -1 */
  double x = qmc_dc_position();
  printf("zitter x=%.10f\n", x);
  if (fabs(x + 0.5 * (1.0 - cos(2.0 * 1.56))) > 0.02) {
    fprintf(stderr, "FAIL zitter x=%.4f\n", x);
    fails++;
  }

  return fails ? 1 : 0;
}
