#include <math.h>
#include <stdio.h>

double *qmc_kr_qbuffer(void);
double *qmc_kr_cbuffer(void);
int qmc_kr_time(void);
int qmc_kr_reset(double k, double hbar, int seed);
int qmc_kr_step(int n);
double qmc_kr_qm2(void);
double qmc_kr_cm2(void);
double qmc_kr_norm(void);
double qmc_kr_diffusion(void);
double qmc_kr_edge(void);

int main(void) {
  qmc_kr_reset(5.0, 1.0, 7);
  qmc_kr_step(100);
  qmc_kr_step(100);

  const double *q = qmc_kr_qbuffer();
  const double *c = qmc_kr_cbuffer();

  printf("local t=%d qm2=%.9f norm=%.10f D=%.6f edge_ok=%d\n", qmc_kr_time(),
         qmc_kr_qm2(), qmc_kr_norm(), qmc_kr_diffusion(),
         qmc_kr_edge() < 1e-20);
  printf("q[1024]=%.9f q[1030]=%.9f\n", q[1024], q[1030]);

  /* classical map is chaotic, so compare it only after a few kicks (a last-bit
   * libm difference grows by about \exp(n)) */
  qmc_kr_reset(5.0, 1.0, 7);
  qmc_kr_step(8);
  printf("classical t=%d cm2=%.5f c0=(%.5f,%.5f)\n", qmc_kr_time(),
         qmc_kr_cm2(), c[0], c[1]);

  /* quantum resonance hbar = 4 \pi, K = 5 */
  qmc_kr_reset(5.0, 4.0 * 3.14159265358979323846, 1);
  qmc_kr_step(10);
  double x = 10.0 * 5.0 / (4.0 * 3.14159265358979323846);
  printf("resonance qm2=%.9f exact=%.9f\n", qmc_kr_qm2(), x * x / 2.0);

  return qmc_kr_reset(13.0, 1.0, 1) == -1 && qmc_kr_step(1000) == -1 ? 0 : 1;
}
