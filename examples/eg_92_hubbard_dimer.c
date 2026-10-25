/*
 * Hubbard dimer
 *
 * NOTE: Two lattice sites, two electrons, hopping t and on-site repulsion U,
 * solved exactly. Singlet ground state crosses over from a band state (U = 0,
 * half weight on doubly occupied sites) to a Heisenberg pair (U >> t, one
 * electron per site, J = 4 t^2 / U): double occupancy falls, <S1.S2> goes from
 * -3/8 to -3/4 and entanglement entropy of one site moves from charge-and-spin
 * to pure spin (ln 2). Quenching Neel state |up, down> then makes staggered
 * moment oscillate at exchange J
 */

#include "../physics/hubbard_dimer.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
  printf(" > Ground state versus U / t (t = 1)\n  %6s %10s %10s %10s %10s\n",
         "U", "E0", "doublons", "<S1.S2>", "S_site");
  for (int i = 0; i <= 8; i++) {
    double u = 2.0 * i;
    double re[HD_DIM];
    double im[HD_DIM];
    hd_solution_t s;
    int idx[6];

    hd_solve(1.0, u, &s);
    hd_n2_states(&s, idx);
    hd_eigenstate(&s, idx[0], re, im);

    printf("  %6.1f %10.5f %10.5f %10.5f %10.5f\n", u, s.e[idx[0]],
           hd_expect_double_occ(re, im), hd_expect_ss(re, im),
           hd_site_entropy(re, im, 0));
  }

  printf(
      "\n > Neel quench at U = 20: 2<Sz1>(t) against \\cos(J t), J = 4 t^2 / "
      "U\n  %6s %10s %10s\n",
      "time", "2<Sz1>", "\\cos(Jt)");
  {
    hd_solution_t s;
    double re0[HD_DIM];
    double im0[HD_DIM];
    double re[HD_DIM];
    double im[HD_DIM];
    const double u = 20.0;
    const double j = 4.0 / u;

    hd_solve(1.0, u, &s);
    hd_basis_state(9, re0, im0);
    for (int i = 0; i <= 8; i++) {
      double time = 2.0 * M_PI / j * i / 8.0;

      hd_evolve(&s, re0, im0, time, re, im);

      printf("  %6.2f %10.5f %10.5f\n", time, 2 * hd_expect_sz(re, im, 0),
             cos(j * time));
    }
  }

  return 0;
}
