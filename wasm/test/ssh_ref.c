#include <math.h>
#include <stdio.h>

double *qmc_ssh_energies(void);
double *qmc_ssh_vectors(void);
double *qmc_ssh_sweep_buffer(void);
int qmc_ssh_set(double w, double noise, double m, int seed);
int qmc_ssh_sweep(double w0, double w1, int count);
int qmc_ssh_winding(void);
double qmc_ssh_zak(void);
double qmc_ssh_xi(void);
double qmc_ssh_gap(void);
double qmc_ssh_edge_weight(int s, int edge);

// Print to 6 digits so libm last-bit differences cannot show up
int main(void) {
  const double *e = qmc_ssh_energies();

  qmc_ssh_set(2.0, 0.0, 0.0, 1);
  printf("topo w=2 wind=%d zak=%.6f xi=%.6f gap=%.6f\n", qmc_ssh_winding(),
         qmc_ssh_zak(), qmc_ssh_xi(), qmc_ssh_gap());
  printf("E19=%.9f E20=%.9f E0=%.6f E39=%.6f edge=%.6f\n", e[19], e[20], e[0],
         e[39], qmc_ssh_edge_weight(19, 3));

  qmc_ssh_set(0.5, 0.0, 0.0, 1);
  printf("trivial w=0.5 wind=%d zak=%.6f xi_nan=%d E19=%.6f E20=%.6f\n",
         qmc_ssh_winding(), qmc_ssh_zak(), isnan(qmc_ssh_xi()) ? 1 : 0, e[19],
         e[20]);

  qmc_ssh_set(2.0, 0.4, 0.0, 5);
  printf("noisy chiral E0+E39=%.6f E19=%.6f\n", e[0] + e[39], e[19]);

  qmc_ssh_set(2.0, 0.0, 0.3, 1);
  printf("staggered E19=%.6f E20=%.6f\n", e[19], e[20]);

  qmc_ssh_set(1.5, 0.0, 0.0, 1);
  qmc_ssh_sweep(0.0, 3.0, 61);
  const double *s = qmc_ssh_sweep_buffer();
  printf("sweep w=0 E0=%.6f w=1.5 s[30*40+19]=%.6f w=3 E39=%.6f\n", s[0],
         s[30 * 40 + 19], s[60 * 40 + 39]);

  return qmc_ssh_set(3.5, 0, 0, 1) == -1 && qmc_ssh_sweep(0, 3, 1000) == -1 ? 0
                                                                            : 1;
}
