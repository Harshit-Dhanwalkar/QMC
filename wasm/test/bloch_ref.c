#include <math.h>
#include <stdio.h>

double *qmc_bl_vec(void);
double qmc_bl_time(void);
int qmc_bl_reset(double x, double y, double z);
int qmc_bl_step(double omega, double delta, double gamma1, double gamma_phi,
                int steps);
double qmc_bl_purity(void);
double qmc_bl_rabi(double t, double omega, double delta);

int main(void) {
  qmc_bl_reset(0.2, 0.9, -0.3);
  qmc_bl_step(1.3, 0.7, 0.3, 0.2, 240);
  const double *v = qmc_bl_vec();
  printf("damped v=(%.9f, %.9f, %.9f) t=%.3f purity=%.9f\n", v[0], v[1], v[2],
         qmc_bl_time(), qmc_bl_purity());

  qmc_bl_reset(0.0, 0.0, 1.0);
  qmc_bl_step(2.0, 0.0, 0.0, 0.0, 157);
  printf("pulse z=%.9f rabi=%.9f\n", v[2],
         1.0 - 2.0 * qmc_bl_rabi(qmc_bl_time(), 2.0, 0.0));

  return qmc_bl_reset(2.0, 0.0, 0.0) == -1 && qmc_bl_step(NAN, 0, 0, 0, 1) == -1
             ? 0
             : 1;
}
