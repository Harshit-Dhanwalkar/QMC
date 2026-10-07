#include <math.h>
#include <stdio.h>

double *qmc_an_buffer(void);
double qmc_an_time(void);
int qmc_an_reset(double w, int seed);
int qmc_an_step(int n);
double qmc_an_norm(void);
double qmc_an_width(void);
double qmc_an_ipr(void);
double qmc_an_energy(void);
double qmc_an_xi(double w, double e);

int main(void) {
  // clean chain: <x^2> = 2 t^2 \tau^2, so width = \sqrt(2) * \tau
  qmc_an_reset(0.0, 1);
  qmc_an_step(1000);
  printf("clean t=%.2f width=%.9f exact=%.9f norm=%.9f\n", qmc_an_time(),
         qmc_an_width(), sqrt(2.0) * qmc_an_time(), qmc_an_norm());

  qmc_an_reset(6.0, 2024);
  qmc_an_step(3000);
  const double *b = qmc_an_buffer();
  printf("disordered t=%.2f width=%.9f ipr=%.9f energy=%.9f norm=%.9f\n",
         qmc_an_time(), qmc_an_width(), qmc_an_ipr(), qmc_an_energy(),
         qmc_an_norm());
  printf("centre=%.9f eps0=%.9f\n", b[400], b[801]);
  printf("xi(W=1,E=1)=%.6f\n", qmc_an_xi(1.0, 1.0));

  return qmc_an_reset(-1.0, 1) == -1 && isnan(qmc_an_xi(1.0, 2.5)) ? 0 : 1;
}
