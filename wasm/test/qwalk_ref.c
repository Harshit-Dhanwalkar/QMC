#include <math.h>
#include <stdio.h>

double *qmc_qw_buffer(void);
int qmc_qw_time(void);
int qmc_qw_reset(double theta, int coin);
int qmc_qw_step(int n);
double qmc_qw_norm(void);
double qmc_qw_mean(void);
double qmc_qw_variance(void);
double qmc_qw_asymptote(void);

int main(void) {
  qmc_qw_reset(0.7853981633974483, 2);
  qmc_qw_step(300);

  const double *b = qmc_qw_buffer();
  printf("hadamard t=%d var=%.9f ratio=%.9f asym=%.9f norm=%.12f\n",
         qmc_qw_time(), qmc_qw_variance(), qmc_qw_variance() / (300.0 * 300.0),
         qmc_qw_asymptote(), qmc_qw_norm());
  printf("P(0)=%.9f P(212)=%.9f classical P(0)=%.9f P(2)=%.9f\n", b[512],
         b[512 + 212], b[1024 + 512], b[1024 + 514]);

  qmc_qw_reset(0.3, 0);
  qmc_qw_step(120);

  printf("biased mean=%.9f var=%.9f\n", qmc_qw_mean(), qmc_qw_variance());

  return qmc_qw_reset(2.0, 0) == -1 && qmc_qw_step(1000) == -1 ? 0 : 1;
}
