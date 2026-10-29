/*
 * Zeeman effect of a p electron: weak field to Paschen-Back
 *
 * NOTE: H = \\xi L.S + b (L_z + 2 S_z) is block diagonal in m_j
 * Below, six 2P levels (in units of \\xi) as field grows from anomalous Zeeman
 * regime (slopes g_J m_j with g = 4/3 and 2/3) to Paschen-Back (slopes m_l + 2
 * m_s) and D-line pattern 2S -> 2P at one field with its polarisations and
 * intensities
 */

#include "../physics/zeeman_full.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  zf_state_t s[6];

  printf(" > 2P levels versus field (units of \\xi)\n\n  %6s", "b/\\xi");
  const char *lab[6] = {"3/2,+3/2", "3/2,+1/2", "3/2,-1/2",
                        "3/2,-3/2", "1/2,+1/2", "1/2,-1/2"};

  for (int k = 0; k < 6; k++) {
    printf(" %9s", lab[k]);
  }
  printf("\n");

  const double bs[] = {0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0};

  for (int i = 0; i < 10; i++) {
    double b = bs[i];

    zf_solve(1.0, b, s);
    printf("  %6.2f", b);
    for (int k = 0; k < 6; k++) {
      printf(" %9.4f", fabs(s[k].energy) < 5e-5 ? 0.0 : s[k].energy);
    }
    printf("\n");
  }

  double b = 0.7;

  zf_solve(1.0, b, s);
  printf("\n > 2S -> 2P lines at b = %.1f \\xi\n\n  %8s %8s %6s %10s\n", b,
         "energy", "weight", "q", "from spin");
  for (int ms = 0; ms < 2; ms++) {
    for (int k = 0; k < 6; k++) {
      int q = (s[k].two_mj - (ms ? -1 : 1)) / 2;
      if (q < -1 || q > 1) {
        continue;
      }

      double w = zf_line_intensity(&s[k], ms, q);
      if (w > 1e-12) {
        printf("  %8.4f %8.4f %6s %10s\n", s[k].energy - zf_s_energy(b, ms), w,
               q == 0 ? "\\pi" : (q > 0 ? "\\sigma+" : "\\sigma-"),
               ms ? "down" : "up");
      }
    }
  }

  return 0;
}
