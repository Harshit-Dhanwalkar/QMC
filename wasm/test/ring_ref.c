#include <math.h>
#include <stdio.h>

int qmc_ar_spectrum(int n);
int qmc_ar_current(int n, int ne);
double qmc_ar_levels(int n, int ne, double flux);
double qmc_ar_current_at(int n, int ne, double flux);
double qmc_ar_point(int n, int b, double flux, double energy);
double qmc_ar_last_reflection(void);
int qmc_ar_trans_flux(int n, int b, double e);
int qmc_ar_trans_energy(int n, int b, double flux);
int qmc_ar_map(int n, int b);
double *qmc_ar_spec_buffer(void);
double *qmc_ar_cur_buffer(void);
double *qmc_ar_tf_buffer(void);
double *qmc_ar_te_buffer(void);
double *qmc_ar_map_buffer(void);
double *qmc_ar_levels_buffer(void);

// Printed to 5 digits so libm last-bit differences cannot show up
int main(void) {
  qmc_ar_spectrum(12);
  const double *sp = qmc_ar_spec_buffer();
  printf("spec l0 i0=%.5f i50=%.5f i100=%.5f l11 i75=%.5f\n", sp[0], sp[50],
         sp[100], sp[11 * 201 + 75]);

  qmc_ar_current(11, 5);
  const double *cu = qmc_ar_cur_buffer();
  printf("cur i20=%.5f i60=%.5f i120=%.5f i140=%.5f\n", cu[20], cu[60], cu[120],
         cu[140]);

  double e = qmc_ar_levels(12, 6, 0.21);
  const double *lv = qmc_ar_levels_buffer();
  printf("levels E6=%.5f l0=%.5f l5=%.5f l11=%.5f I=%.5f\n", e, lv[0], lv[5],
         lv[11], qmc_ar_current_at(12, 6, 0.21));

  double t = qmc_ar_point(16, 8, 0.15, 0.6);
  printf("point T=%.5f R=%.5f\n", t, qmc_ar_last_reflection());

  qmc_ar_trans_flux(16, 8, 0.4);
  const double *tf = qmc_ar_tf_buffer();
  printf("tflux t0=%.5f t50=%.5f t100=%.5f t200=%.5f\n", tf[0], tf[50], tf[100],
         tf[200]);

  qmc_ar_trans_energy(16, 5, 0.3);
  const double *te = qmc_ar_te_buffer();
  printf("tenergy t0=%.5f t60=%.5f t100=%.5f t180=%.5f\n", te[0], te[60],
         te[100], te[180]);

  qmc_ar_map(10, 5);
  const double *mp = qmc_ar_map_buffer();
  printf("map %.5f %.5f %.5f %.5f %.5f\n", mp[0], mp[20 * 64 + 10],
         mp[32 * 64 + 32], mp[63 * 64 + 40], mp[40 * 64 + 63]);

  if (qmc_ar_spectrum(3) != -1 || qmc_ar_current(8, 9) != -1 ||
      qmc_ar_trans_flux(8, 0, 0.0) != -1 ||
      qmc_ar_trans_flux(8, 4, 2.5) != -1 || qmc_ar_map(40, 3) != -1 ||
      !isnan(qmc_ar_point(8, 8, 0.0, 0.0)) ||
      !isnan(qmc_ar_levels(8, 9, 0.0))) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
