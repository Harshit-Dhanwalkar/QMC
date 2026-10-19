#include <math.h>
#include <stdio.h>

double *qmc_cb_field_buffer(void);
int qmc_cb_mode_count(void);
int qmc_cb_mode_m(int i);
int qmc_cb_mode_n(int i);
int qmc_cb_mode_sine(int i);
double qmc_cb_mode_k(int i);
int qmc_cb_show_mode(int i);
double qmc_cb_launch(double x0, double y0, double s, double px, double py);
double qmc_cb_energy(void);
double qmc_cb_evolve(double t);

// Printed to 5 digits so libm last-bit differences cannot show up
int main(void) {
  const double *f = qmc_cb_field_buffer();

  printf("modes=%d\n", qmc_cb_mode_count() > 150);
  printf("mode0 m=%d n=%d k=%.5f  mode1 m=%d n=%d s=%d k=%.5f  mode9 m=%d n=%d "
         "k=%.5f\n",
         qmc_cb_mode_m(0), qmc_cb_mode_n(0), qmc_cb_mode_k(0), qmc_cb_mode_m(1),
         qmc_cb_mode_n(1), qmc_cb_mode_sine(1), qmc_cb_mode_k(1),
         qmc_cb_mode_m(9), qmc_cb_mode_n(9), qmc_cb_mode_k(9));
  qmc_cb_show_mode(0);

  printf("ground f[40,40]=%.5f f[10,40]=%.5f f[0,0]=%.5f\n", f[40 * 80 + 40],
         f[40 * 80 + 10], f[0]);
  qmc_cb_show_mode(5);

  printf("mode5 f[30,50]=%.5f f[50,30]=%.5f\n", f[30 * 80 + 50],
         f[50 * 80 + 30]);
  double kept = qmc_cb_launch(-0.4, 0.0, 0.1, 8.0, 3.0);
  double energy = qmc_cb_energy();
  printf("launch kept=%.4f energy=%.3f\n", kept, energy);

  double p0 = qmc_cb_evolve(0.0);
  printf("t=0 total=%.4f f[40,24]=%.5f\n", p0, f[40 * 80 + 24]);

  double p1 = qmc_cb_evolve(0.05);
  printf("t=0.05 total=%.4f f[40,24]=%.5f f[44,45]=%.5f\n", p1, f[40 * 80 + 24],
         f[44 * 80 + 45]);

  if (qmc_cb_show_mode(-1) != -1 ||
      qmc_cb_launch(2.0, 0.0, 0.1, 0.0, 0.0) != -1.0 ||
      qmc_cb_launch(0.0, 0.0, 0.1, 50.0, 0.0) != -1.0) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
