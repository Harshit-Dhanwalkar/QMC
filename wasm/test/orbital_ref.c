#include <math.h>
#include <stdio.h>

double *qmc_orb_buffer(void);
int qmc_orb_sample(int n, int l, int m, int real_form, int count, int seed);
double qmc_orb_energy_ev(int n);

int main(void) {
  const int orbs[][4] = {
      {1, 0, 0, 0}, {2, 1, 1, 0}, {3, 2, -1, 1}, {4, 3, 2, 1}, {6, 2, 0, 0}};

  for (int k = 0; k < 5; k++) {
    if (qmc_orb_sample(orbs[k][0], orbs[k][1], orbs[k][2], orbs[k][3], 2000,
                       7 + k) != 0) {
      return 1;
    }

    const double *p = qmc_orb_buffer();
    double sr = 0.0, sz = 0.0, sp = 0.0;

    for (int i = 0; i < 2000; i++) {
      sr += sqrt(p[4 * i] * p[4 * i] + p[4 * i + 1] * p[4 * i + 1] +
                 p[4 * i + 2] * p[4 * i + 2]);
      sz += p[4 * i + 2];
      sp += p[4 * i + 3];
    }

    printf("orb %d%d%+d r=%.6f z=%.6f ph=%.6f\n", orbs[k][0], orbs[k][1],
           orbs[k][2], sr, sz, sp);
  }

  printf("E2=%.6f\n", qmc_orb_energy_ev(2));

  return qmc_orb_sample(1, 0, 0, 0, 0, 1) == -1 ? 0 : 1;
}
