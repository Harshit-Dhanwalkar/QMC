/*
 * Josephson junction (RCSJ model)
 *
 * NOTE: Phase difference across a Josephson junction is a damped, driven
 * pendulum in a tilted washboard potential. Below critical current phase is
 * trapped (zero voltage); above it phase runs and a voltage appears. A lightly
 * damped junction is hysteretic (it switches at i = 1 but retraps near 4 / \pi
 * Q), and an ac drive locks running phase to drive, giving flat Shapiro steps
 * at V = n \Omega
 */

#include "../physics/josephson.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  double up[21];
  double down[21];
  double steps[21];
  jj_params_t hyst = {3.0, 0.0, 0.0, 0.5};
  jj_params_t shap = {0.7, 0.0, 1.2, 0.8};

  jj_sweep(&hyst, 0.0, 2.0, 21, up);
  jj_sweep(&hyst, 2.0, 0.0, 21, down);
  jj_sweep(&shap, 0.0, 2.0, 21, steps);

  printf(" > I-V curves (Q = 3 hysteresis; Q = 0.7 with ac drive 1.2 at \\Omega "
         "= 0.8, steps at V = 0.8)\n  %6s %12s %12s %12s\n",
         "i", "V up", "V down", "V shapiro");
  for (int k = 0; k <= 20; k++) {
    printf("  %6.2f %12.4f %12.4f %12.4f\n", 0.1 * k,
           fabs(up[k]) < 5e-5 ? 0.0 : up[k],
           fabs(down[20 - k]) < 5e-5 ? 0.0 : down[20 - k],
           fabs(steps[k]) < 5e-5 ? 0.0 : steps[k]);
  }

  return 0;
}
