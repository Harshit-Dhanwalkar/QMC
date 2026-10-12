#include <math.h>
#include <stdio.h>

double *qmc_wg_w_buffer(void);
double *qmc_wg_x_buffer(void);
double *qmc_wg_p_buffer(void);
double *qmc_wg_pn_buffer(void);
int qmc_wg_state(int kind, double amp, double param, double chi_t);
double qmc_wg_nbar(void);
double qmc_wg_negativity(void);
double qmc_wg_wmin(void);
double qmc_wg_wmax(void);
double qmc_wg_extent(void);

// Printed to 5 digits so libm last-bit differences cannot show up
int main(void) {
  const double *w = qmc_wg_w_buffer(), *mx = qmc_wg_x_buffer(),
               *mp = qmc_wg_p_buffer(), *pn = qmc_wg_pn_buffer();

  printf("extent %.4f\n", qmc_wg_extent());

  qmc_wg_state(0, 2.5, 0.0, 0.0);
  printf("coherent nbar=%.5f neg=%.5f wmax=%.5f mx64=%.5f\n", qmc_wg_nbar(),
         qmc_wg_negativity(), qmc_wg_wmax(), mx[64]);

  qmc_wg_state(0, 2.5, 0.0, 1.5707963267948966);
  printf("kerr cat nbar=%.5f neg=%.5f wmin=%.5f w[64,64]=%.5f p64=%.5f\n",
         qmc_wg_nbar(), qmc_wg_negativity(), qmc_wg_wmin(), w[64 * 128 + 64],
         mp[64]);

  qmc_wg_state(1, 3.0, 0.0, 0.0);
  printf("fock3 nbar=%.5f wmin=%.5f w[63,63]=%.5f pn3=%.5f\n", qmc_wg_nbar(),
         qmc_wg_wmin(), w[63 * 128 + 63], pn[3]);

  qmc_wg_state(2, 2.0, 0.0, 0.0);
  printf("even cat neg=%.5f pn0=%.5f pn1=%.5f pn2=%.5f\n", qmc_wg_negativity(),
         pn[0], pn[1], pn[2]);

  qmc_wg_state(4, 0.8, 0.0, 0.0);
  printf("squeezed nbar=%.5f wmax=%.5f\n", qmc_wg_nbar(), qmc_wg_wmax());

  if (qmc_wg_state(0, 5.0, 0.0, 0.0) != -1 ||
      qmc_wg_state(1, 20.0, 0.0, 0.0) != -1 ||
      qmc_wg_state(9, 1.0, 0.0, 0.0) != -1) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
