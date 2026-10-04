/*
 * Tight-binding band structure and band topology
 *
 * NOTE: Every model here is the same object (a Bloch Hamiltonian H(k) built
 * from a list of hoppings) and the same two tools are applied to all of them:
 * band energies E_n(k) and gauge-invariant Berry phases of band groups
 * (Fukui-Hatsugai-Suzuki link variables)
 *
 *  1. Graphene: massless Dirac cones at K and K', density of states with
 *     van Hove peaks at |E| = t
 *  2. Haldane model: a sublattice mass M opens a trivial gap; complex
 *     second-neighbor hopping t2 e^{i phi} opens a topological one. Scanning
 *     M maps phase diagram C = \sign(t2 \sin(\phi)) for |M| < 3 \sqrt(3) |t2
 *     \sin(\phi)|, 0 beyond it
 *  3. SSH chain: Zak phase jumps 0 -> \pi when the intercell bond overtakes
 *     the intracell one
 *  4. Kagome lattice: the exactly flat band at -2t and a Dirac point
 *  5. Hofstadter butterfly bands: Chern numbers from the same code, including a
 *     band group for flux 1/4 where two bands touch
 */

#include "../physics/tight_binding.h"
#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Avoid printing "-0.0000" for round-off-sized negative values.
static double clean(double x) { return fabs(x) < 5e-5 ? 0.0 : x; }

static void graphene_demo(void) {
  printf("  === Graphene ===\n");

  tb_model_t *g = tb_model_graphene(1.0);
  if (!g) {
    return;
  }

  double e[2];
  tb_bands(g, 0.0, 0.0, e);
  printf("  \\Gamma: E = %+.4f, %+.4f   (band edges +-3t)\n", e[0], e[1]);
  tb_bands(g, 2.0 * M_PI / 3.0, -2.0 * M_PI / 3.0, e);
  printf("  K    : E = %+.4f, %+.4f   (Dirac point, gap %.1e)\n", e[0], e[1],
         e[1] - e[0]);

  enum { NE = 401 };
  static double dos[NE];
  tb_dos(g, 120, -4.0, 4.0, NE, 0.06, dos);
  int imax = NE / 2;
  for (int i = NE / 2; i < NE; i++) {
    if (dos[i] > dos[imax]) {
      imax = i;
    }
  }

  printf("  DOS at E = 0: %.4f   van Hove peak at E = %.3f (height %.3f)\n",
         dos[NE / 2], -4.0 + 8.0 * imax / (NE - 1), dos[imax]);

  tb_model_free(g);
}

static void haldane_demo(void) {
  printf("  === Haldane model: Chern number of the lower band ===\n");

  const double t1 = 1.0;
  const double t2 = 0.1;
  printf("  t1 = %.1f, t2 = %.2f, transition at |M| = %.4f\n", t1, t2,
         3.0 * sqrt(3.0) * t2);
  printf("  %8s | %12s %12s\n", "M", "C(\\phi=+\\pi/2)", "C(\\phi=-\\pi/2)");
  for (double mass = -0.8; mass <= 0.801; mass += 0.2) {
    tb_model_t *a = tb_model_haldane(t1, t2, M_PI / 2.0, mass);
    tb_model_t *b = tb_model_haldane(t1, t2, -M_PI / 2.0, mass);
    double ca = NAN;
    double cb = NAN;
    if (a) {
      tb_chern_number(a, 0, 1, 48, &ca);
    }
    if (b) {
      tb_chern_number(b, 0, 1, 48, &cb);
    }

    printf("  %+8.2f | %12.4f %12.4f\n", mass, clean(ca), clean(cb));

    tb_model_free(a);
    tb_model_free(b);
  }

  tb_model_t *m = tb_model_haldane(t1, t2, M_PI / 2.0, 0.3);
  if (m) {
    double w = NAN;
    tb_wilson_winding(m, 0, 1, 48, 48, &w);
    printf("  Wilson-loop winding at M = 0.3: %.4f (equals the Chern number)\n",
           w);
    for (int v = 1; v >= -1; v -= 2) {
      double k = v * 2.0 * M_PI / 3.0;
      double e[2];
      tb_bands(m, k, -k, e);

      printf("  valley %+d gap: %.4f (analytic %.4f)\n", v, e[1] - e[0],
             tb_haldane_valley_gap(t2, M_PI / 2.0, 0.3, v));
    }

    tb_model_free(m);
  }
}

static void ssh_demo(void) {
  printf("  === SSH chain: Zak phase ===\n");

  const double v = 1.0;
  for (double w = 0.2; w <= 2.01; w += 0.4) {
    tb_model_t *m = tb_model_ssh(v, w);
    double z = NAN;
    if (m) {
      tb_wilson_phase(m, 0, 1, 0.0, 200, &z);
    }
    if (isnan(z)) {
      printf("  v = %.1f, w = %.1f : gapless, Zak phase undefined\n", v, w);
    } else {
      printf("  v = %.1f, w = %.1f : Zak phase = %+.4f pi\n", v, w,
             clean(z / M_PI));
    }

    tb_model_free(m);
  }
}

static void kagome_demo(void) {
  printf("  === Kagome lattice ===\n");

  tb_model_t *m = tb_model_kagome(1.0);
  if (!m) {
    return;
  }

  const double pts[3][2] = {
      {0.0, 0.0}, {1.0, 2.0}, {2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0}};
  for (int i = 0; i < 3; i++) {
    double e[3];
    tb_bands(m, pts[i][0], pts[i][1], e);

    printf("  k = (%.3f, %.3f): E = %+.4f %+.4f %+.4f\n", pts[i][0], pts[i][1],
           e[0], e[1], e[2]);
  }

  printf("  flat band energy -2t = %.4f\n", tb_kagome_flat_band_energy(1.0));

  tb_model_free(m);
}

static void hofstadter_demo(void) {
  printf("  === Hofstadter bands ===\n");
  const int flux[3][2] = {{1, 3}, {2, 5}, {1, 4}};
  for (int f = 0; f < 3; f++) {
    int p = flux[f][0];
    int q = flux[f][1];

    tb_model_t *m = tb_model_hofstadter(p, q, 1.0);
    if (!m) {
      continue;
    }

    printf("  flux %d/%d:", p, q);
    if (q == 4) {
      // bands 1 and 2 touch: treat them as one group
      double c0 = NAN;
      double c12 = NAN;
      double c3 = NAN;

      tb_chern_number(m, 0, 1, 32, &c0);
      tb_chern_number(m, 1, 2, 32, &c12);
      tb_chern_number(m, 3, 1, 32, &c3);

      printf("  C = %+.3f | bands 1+2 together %+.3f | %+.3f\n", c0, c12, c3);
    } else {
      for (int b = 0; b < q; b++) {
        double c = NAN;

        tb_chern_number(m, b, 1, 40, &c);

        printf("  %+.3f", c);
      }

      printf("\n");
    }

    tb_model_free(m);
  }
}

int main(void) {
  printf(" > Tight-binding band structure and topology (hbar = 1)\n");

  graphene_demo();
  haldane_demo();
  ssh_demo();
  kagome_demo();
  hofstadter_demo();

  return 0;
}
