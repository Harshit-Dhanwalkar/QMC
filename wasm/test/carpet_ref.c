#include <math.h>
#include <stdio.h>

double *qmc_qc_carpet_buffer(void);
double *qmc_qc_slice_buffer(void);
double *qmc_qc_curve_buffer(void);
int qmc_qc_set(double x0, double sigma, double n0);
double qmc_qc_slice(double frac);
double qmc_qc_weight(void);
double qmc_qc_energy(void);

// Printed to 5 digits so libm last-bit differences cannot show up
int main(void) {
  qmc_qc_set(0.3, 0.05, 20.0);
  const double *car = qmc_qc_carpet_buffer(), *sl = qmc_qc_slice_buffer(),
               *cv = qmc_qc_curve_buffer();
  printf("weight=%.5f energy=%.4f\n", qmc_qc_weight(), qmc_qc_energy());
  printf("carpet t0 x57=%.5f x134=%.5f t192 x40=%.5f t383 x57=%.5f\n", car[57],
         car[134], car[192 * 192 + 40], car[383 * 192 + 57]);
  printf("curve c0=%.5f c128=%.5f c256=%.5f c512=%.5f\n", cv[0], cv[128],
         cv[256], cv[512]);

  double r = qmc_qc_slice(0.5);
  printf("slice half r=%.5f s57=%.5f s134=%.5f\n", r, sl[57], sl[134]);

  r = qmc_qc_slice(0.25);
  printf("slice quarter r=%.5f s57=%.5f s134=%.5f\n", r, sl[57], sl[134]);

  qmc_qc_set(0.5, 0.1, 0.0);
  printf("centred weight=%.5f energy=%.4f half=%.5f\n", qmc_qc_weight(),
         qmc_qc_energy(), qmc_qc_slice(0.5));

  if (qmc_qc_set(1.5, 0.05, 1.0) != -1 || qmc_qc_set(0.5, 0.001, 1.0) != -1 ||
      !isnan(qmc_qc_slice(NAN))) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
