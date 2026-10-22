#include <math.h>
#include <stdio.h>

int qmc_dt_set(double v1, double b1, double w, double b2, double v2, double vw);
double qmc_dt_length(void);
double qmc_dt_well_level(int n);
int qmc_dt_curves(double emax);
double qmc_dt_point(double e);
double qmc_dt_last_reflection(void);
double qmc_dt_dwell(double e);
int qmc_dt_profile(double e, double xmin, double xmax);
double *qmc_dt_te_buffer(void);
double *qmc_dt_dwell_buffer(void);
double *qmc_dt_psi_buffer(void);
double *qmc_dt_pot_buffer(void);

// Printed to 5 significant digits so libm last-bit differences cannot show up
int main(void) {
  qmc_dt_set(30.0, 0.5, 3.0, 0.5, 30.0, 0.0);
  printf("length %.5f well1 %.5f well2 %.5f\n", qmc_dt_length(),
         qmc_dt_well_level(1), qmc_dt_well_level(2));

  qmc_dt_curves(20.0);
  const double *te = qmc_dt_te_buffer();
  const double *dw = qmc_dt_dwell_buffer();
  printf("T i0=%.5f i20=%.5f i100=%.5f i300=%.5f i400=%.5f\n", te[0], te[20],
         te[100], te[300], te[400]);

  printf("dwell i20=%.5f i100=%.5f i300=%.5f\n", dw[20], dw[100], dw[300]);
  double t = qmc_dt_point(2.2);
  printf("point T=%.5f R=%.5f dwell=%.5f\n", t, qmc_dt_last_reflection(),
         qmc_dt_dwell(2.2));

  qmc_dt_profile(2.2, -2.0, 6.0);
  const double *ps = qmc_dt_psi_buffer(), *pt = qmc_dt_pot_buffer();
  printf("\\psi %.5f %.5f %.5f %.5f %.5f | pot %.5f %.5f %.5f\n", ps[0], ps[100],
         ps[200], ps[300], ps[400], pt[0], pt[150], pt[200]);

  qmc_dt_set(20.0, 0.3, 2.0, 0.6, 35.0, -3.0);
  t = qmc_dt_point(4.0);
  printf("asym T=%.5f R=%.5f length %.5f\n", t, qmc_dt_last_reflection(),
         qmc_dt_length());
  if (qmc_dt_set(-1.0, 0.5, 3.0, 0.5, 30.0, 0.0) != -1 ||
      qmc_dt_set(30.0, 0.0, 3.0, 0.5, 30.0, 0.0) != -1 ||
      qmc_dt_curves(0.0) != -1 || !isnan(qmc_dt_point(-1.0)) ||
      qmc_dt_profile(1.0, 2.0, 1.0) != -1) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
