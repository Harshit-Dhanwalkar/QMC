/*
 * Zeeman Effect (Weak-Field / Anomalous)
 *
 * Lande g-factors for hydrogen-like S,P,D terms, cross-checked against a
 * direct Clebsch-Gordan-coupling computation of <Sz>, and a "Zeeman fan"
 * of energy vs B for the 2P_3/2 term's four mj sublevels.
 */

#include "../core/utils.h"
#include "../export/plot.h"
#include "../physics/zeeman.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void print_term(const char *name, int l, int j_2) {
  double g = zeeman_lande_g_factor(l, j_2);
  printf("   %-6s  l=%d j=%d/2   g_J = %.6f\n", name, l, j_2, g);

  for (int mj_2 = -j_2; mj_2 <= j_2; mj_2 += 2) {
    double sz_closed = (g - 1.0) * (mj_2 / 2.0);
    double sz_coupling = zeeman_sz_expect_from_coupling(l, j_2, mj_2);
    printf("      mj=%4.1f   <Sz> closed=%.6f  from-coupling=%.6f\n",
           mj_2 / 2.0, sz_closed, sz_coupling);
  }
  printf("\n");
}

int main(void) {
  printf(" > Zeeman Effect: Lande g-factors and Weak-Field Splitting\n\n");

  print_term("S1/2", 0, 1);
  print_term("P1/2", 1, 1);
  print_term("P3/2", 1, 3);
  print_term("D3/2", 2, 3);
  print_term("D5/2", 2, 5);

  // "Zeeman fan": energy shift vs B for 2P_3/2 (l=1, j=3/2), 4 mj
  // sublevels, \mu_B = 1 (natural units).
  int l = 1, j_2 = 3;
  int n_mj = j_2 + 1;
  int N = 100;
  double *B = linspace(0.0, 2.0, N);
  if (!B) {
    return 1;
  }

  double **E = malloc(n_mj * sizeof *E);
  const double **E_const = malloc(n_mj * sizeof *E_const);
  char **labels = malloc(n_mj * sizeof *labels);
  if (!E || !E_const || !labels) {
    free(B);
    free(E);
    free((void *)E_const);
    free(labels);

    return 1;
  }

  int idx = 0;
  for (int mj_2 = -j_2; mj_2 <= j_2; mj_2 += 2, idx++) {
    E[idx] = malloc(N * sizeof **E);
    for (int i = 0; i < N; i++) {
      E[idx][i] = zeeman_energy_shift(l, j_2, mj_2, B[i], 1.0);
    }
    E_const[idx] = E[idx];

    labels[idx] = malloc(16);
    snprintf(labels[idx], 16, "m_j=%.1f", mj_2 / 2.0);
  }

  plot_opts_t opts = {0};
  opts.title = "2P_{3/2} Zeeman fan";
  opts.xlabel = "B (\\mu_B units)";
  opts.ylabel = "\\Delta E / mu_B";
  opts.tex_text = 1;

  plot_lines("zeeman_fan_2p32", PLOT_FORMAT_PNG, B, E_const, n_mj, N,
             (const char **)labels, &opts);
  printf("   Saved zeeman_fan_2p32.png");

  for (int i = 0; i < n_mj; i++) {
    free(E[i]);
    free(labels[i]);
  }
  free(E);
  free((void *)E_const);
  free(labels);
  free(B);

  return 0;
}
