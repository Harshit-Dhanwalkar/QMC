/*
 * Su-Schrieffer-Heeger chain: a topological insulator in one dimension
 *
 * NOTE: A chain of 20 cells with intra-cell hopping v = 1. For inter-cell
 * hopping w < v chain is trivial and has a gap of 2|v - w|. For w > v winding
 * number is 1, Zak phase is pi and two states appear inside gap, one at each
 * end of chain, with energy +/- (w^2 - v^2)/w (v/w)^n (exponentially close to
 * zero) and weight (1 - r) r^c on cell c, r = (v/w)^2
 * Random hopping leaves them at zero (chiral symmetry); a staggered potential
 * moves them to +/- m
 */

#include "../physics/ssh_chain.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { CELLS = 20, SITES = 2 * CELLS };

int main(void) {
  double bond[SITES - 1];
  double e[SITES];
  double vec[SITES * SITES];
  const double ws[] = {0.5, 0.9, 1.1, 1.5, 2.0, 3.0};

  printf(" > SSH chain, %d cells, v = 1\n\n", CELLS);
  printf("  %5s  %7s  %7s  %9s  %12s  %12s\n", "w", "winding", "(Zak/\\pi)",
         "gap", "E(+edge)", "exact");

  for (int i = 0; i < 6; i++) {
    double w = ws[i];

    ssh_chain_bonds(bond, CELLS, 1.0, w, 0.0, 1);
    ssh_chain_solve(CELLS, bond, NULL, e, vec);

    double exact = w > 1.0 ? (w * w - 1.0) / w * pow(1.0 / w, CELLS) : NAN;

    printf("  %5.2f  %7d  %7.3f  %9.4f  %12.3e  %12.3e\n", w,
           ssh_winding_number(1.0, w), ssh_zak_phase(1.0, w) / M_PI,
           2.0 * fabs(1.0 - w), e[CELLS], exact);
  }

  // hopping disorder does not move edge states
  ssh_chain_bonds(bond, CELLS, 1.0, 2.0, 0.5, 3);
  ssh_chain_solve(CELLS, bond, NULL, e, NULL);
  printf("\n  w = 2 with 50%% random hopping: E = %.2e, %.2e (chiral)\n",
         e[CELLS - 1], e[CELLS]);

  // staggered potential does
  double eps[SITES];

  ssh_chain_bonds(bond, CELLS, 1.0, 2.0, 0.0, 1);
  ssh_chain_onsite(eps, CELLS, 0.3);
  ssh_chain_solve(CELLS, bond, eps, e, NULL);

  printf("  w = 2 with staggered potential m = 0.3: E = %.4f, %.4f\n",
         e[CELLS - 1], e[CELLS]);

  return 0;
}
