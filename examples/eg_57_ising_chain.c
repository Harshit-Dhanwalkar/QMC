/*
 * Transverse-Field Ising Model: Quantum Phase Transition
 *
 * physics/ising_chain.c implements exact diagonalization for the TFIM (H = -J
 * \sum \sigma^z_i \sigma^z_{i+1} - h \sum \sigma^x_i) plus its exact
 * Jordan-Wigner/Bogoliubov ground-state energy
 *
 * NOTE: physics: at h=0 the ground state is a doubly-degenerate classical
 * ferromagnet (all spins up or all down, spontaneously breaking the Z2
 * spin-flip symmetry); at h -> \infty every spin aligns with the transverse
 * field instead (a unique paramagnetic ground state, Z2 symmetric). These two
 * phases are separated by a continuous quantum phase transition at h_c = J (for
 * N -> \infty chain), visible in finite N as an increasingly sharp crossover in
 * ground-state energy and transverse magnetization order parameter <\sigma^x>
 * as N grows.
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/sparse.h"
#include "../export/plot.h"
#include "../physics/ising_chain.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define J_COUPLING 1.0
#define N_H_POINTS 41
#define H_MIN 0.0
#define H_MAX 2.0

/* <\sigma^x_0> in ground state, from its Lanczos eigenvector.
 *  \sigma^x_0 flips bit 0 of computational basis, so
 *  <\psi|\sigma^x_0|\psi> = \sum_s conj(\psi[s]) * \psi[s^1]. */
static double transverse_magnetization(const cmatrix_t *psi, int dim) {
  complex_t acc = c_zero();

  for (int s = 0; s < dim; s++) {
    complex_t a = CMAT(psi, s, 0);
    complex_t b = CMAT(psi, s ^ 1, 0);

    acc = c_add(acc, c_mul(c_conj(a), b));
  }

  return acc.re;
}

int main(void) {
  printf(" > Transverse-Field Ising Model: quantum phase transition\n\n");

  double h_vals[N_H_POINTS];
  double dh = (H_MAX - H_MIN) / (N_H_POINTS - 1);
  for (int i = 0; i < N_H_POINTS; i++) {
    h_vals[i] = H_MIN + i * dh;
  }

  /* NOTE: Ground energy per site: numerical ED vs the exact JW/Bogoliubov
   * formula, for N=10 should overlap exactly, demonstrating the module against
   * its own independent validation reference. */
  int N_energy = 10;
  double *E0_numeric = malloc((size_t)N_H_POINTS * sizeof *E0_numeric);
  double *E0_exact = malloc((size_t)N_H_POINTS * sizeof *E0_exact);

  printf("   %-8s %-16s %-16s %-10s\n", "h/J", "E0/N (ED)", "E0/N (exact)",
         "diff");
  for (int i = 0; i < N_H_POINTS; i++) {
    sparse_matrix_t *H = ising_hamiltonian(N_energy, J_COUPLING, h_vals[i], 1);
    lanczos_result_t *res = lanczos_eigs(H, 1, 150, 1e-12);

    E0_numeric[i] = res->values[0] / N_energy;
    E0_exact[i] =
        ising_exact_ground_energy_per_site(N_energy, J_COUPLING, h_vals[i]);
    if (i % 8 == 0) {
      printf("   %-8.3f %-16.8f %-16.8f %-10.2e\n", h_vals[i], E0_numeric[i],
             E0_exact[i], fabs(E0_numeric[i] - E0_exact[i]));
    }

    lanczos_free(res);
    sparse_free(H);
  }

  plot_opts_t opts = {0};
  opts.title = "TFIM ground energy per site: ED vs exact solution (N=10)";
  opts.xlabel = "h / J";
  opts.ylabel = "E0 / N";

  const double *ys1[2] = {E0_numeric, E0_exact};
  const char *labels1[2] = {"Lanczos ED", "Exact (Jordan-Wigner+Bogoliubov)"};

  plot_lines("ising_ground_energy", PLOT_FORMAT_PNG, h_vals, ys1, 2, N_H_POINTS,
             labels1, &opts);
  printf("\n   Saved ising_ground_energy.png (two curves should overlap)\n");

  /* Transverse magnetization order parameter across transition, for 2 system
   * sizes to show finite-size crossover sharpening toward true N->infty phase
   * transition at h_c=J. */
  const int Ns[2] = {6, 10};
  double *mag[2];
  for (int n_idx = 0; n_idx < 2; n_idx++) {
    int N = Ns[n_idx];

    mag[n_idx] = malloc((size_t)N_H_POINTS * sizeof *mag[n_idx]);
    for (int i = 0; i < N_H_POINTS; i++) {
      sparse_matrix_t *H = ising_hamiltonian(N, J_COUPLING, h_vals[i], 1);
      int dim = 1 << N;
      lanczos_result_t *res = lanczos_eigs(H, 1, dim < 150 ? dim : 150, 1e-12);

      mag[n_idx][i] = transverse_magnetization(res->vectors, dim);

      lanczos_free(res);
      sparse_free(H);
    }
  }

  plot_opts_t opts2 = {0};

  opts2.title = "Transverse magnetization <\\sigma^x> vs h/J";
  opts2.xlabel = "h / J";
  opts2.ylabel = "<\\sigma^x>";

  const double *ys2[2] = {mag[0], mag[1]};
  char label_n6[32], label_n10[32];
  snprintf(label_n6, sizeof label_n6, "N=%d", Ns[0]);
  snprintf(label_n10, sizeof label_n10, "N=%d", Ns[1]);
  const char *labels2[2] = {label_n6, label_n10};

  plot_lines("ising_magnetization", PLOT_FORMAT_PNG, h_vals, ys2, 2, N_H_POINTS,
             labels2, &opts2);
  printf("   Saved ising_magnetization.png (N=10 crossover near h/J=1 should "
         "be visibly sharper than N=6's)\n");

  /* Mid-chain entanglement entropy across the transition, for N=10
   * NOTE: should peak near h_c=J (entanglement is largest at a critical point,
   * where correlations are long-ranged) and fall off deep in either phase
   * (product-state-like ground states on both the h=0 ferromagnetic side and
   * h->infinity paramagnetic side have low entanglement). */
  int N_ee = 10;
  double *S_mid = malloc((size_t)N_H_POINTS * sizeof *S_mid);
  for (int i = 0; i < N_H_POINTS; i++) {
    sparse_matrix_t *H = ising_hamiltonian(N_ee, J_COUPLING, h_vals[i], 1);
    lanczos_result_t *res = lanczos_eigs(H, 1, 150, 1e-12);

    S_mid[i] = ising_entanglement_entropy(res->vectors, N_ee, N_ee / 2);

    lanczos_free(res);
    sparse_free(H);
  }

  plot_opts_t opts3 = {0};

  opts3.title = "Mid-chain entanglement entropy vs h/J (N=10)";
  opts3.xlabel = "h / J";
  opts3.ylabel = "S (L_A = N/2)";

  plot_line("ising_entanglement", PLOT_FORMAT_PNG, h_vals, S_mid, N_H_POINTS,
            &opts3);
  printf("   Saved ising_entanglement.png (should peak near h/J=1, critical "
         "point, and fall off in both phases)\n");

  free(S_mid);
  free(E0_numeric);
  free(E0_exact);
  free(mag[0]);
  free(mag[1]);

  return 0;
}
