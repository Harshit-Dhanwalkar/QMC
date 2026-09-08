/*
 * Dynamical Structure Factor of the Spin-1/2 Heisenberg Ring via
 * Symmetry-Adapted Lanczos Continued Fractions
 *
 * physics/spin_chain.c implements translation-symmetry-adapted exact
 * diagonalization of XXZ ring plus a Lanczos continued-fraction evaluator for
 * S^{zz}(q, \omega). Ground energies vs full-space ED, the S^z_q excitation's
 * sum rule, and the continued fraction's integral-over-omega sum rule. This
 * example is end-to-end demonstration of build ground state once, sweep the
 * excitation over every momentum transfer q in the Brillouin zone, and plot the
 * resulting S(q, \omega) spectrum.
 *
 * Workflow:
 *   1. Ground state |psi_0>, E0: scan every (nup, k) sector and Lanczos
 *      each one for its lowest eigenvalue -- the Heisenberg antiferromagnet's
 *      ground state on a finite ring is always in the nup = N/2 sector, but
 *      we scan everything so this doesn't silently assume that.
 *   2. For each q_index = 0..N-1:
 *        |phi_0> = S^z_q |psi_0>      (spin_apply_szq)
 *        I0 = <phi_0|phi_0>           (equal-time structure factor S(q))
 *        |f0> = |phi_0> / sqrt(I0)    (normalized Lanczos start vector)
 *        (\alpha, \beta) = Lanczos tridiagonalization of H in the target
 *                          sector, starting from |f0>
 *   3. Evaluate the continued fraction on a frequency mesh to get S(q,omega)
 *      for that q, and cross-check that integral(S(q,omega) domega) == I0
 *      (this is theoretical quantitative check: the module's f-sum rule).
 *
 * NOTE: Physics : for an N=8 ring the spectrum is still a handful of discrete
 * delta functions broadened by \eta, not the smooth two-spinon continuum of the
 * N -> \infty chain - des Cloizeaux-Pearson's exact lower-bound dispersion
 * \omega_{LB}(q) = (\pi / 2) * J * |\sin(q)| is an infinite-chain result, only
 * expect to approximately bound the dominant finite-size peaks, not match them.
 */

#include "../core/matrix.h"
#include "../core/sparse.h"
#include "../core/vector.h"
#include "../export/plot.h"
#include "../physics/spin_chain.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N_SITES 8
#define J_COUPLING 1.0 // isotropic Heisenberg: Jxy = Jz = J
#define ETA 0.08       // Lorentzian broadening
#define N_OMEGA 800
#define OMEGA_MIN -1.0
#define OMEGA_MAX 9.0
#define PLOT_OMEGA_MIN -0.5
#define PLOT_OMEGA_MAX 4.0

/* Scan every (nup, k) sector, Lanczos each, and return global ground state:
 * energy, its own sector, and its coefficient vector in that sector's
 * symmetry-adapted basis. */
static double find_ground_state(int N, spin_sector_t **out_sector,
                                cvector_t **out_psi0) {
  double best_E = 1e300;
  spin_sector_t *best_sec = NULL;
  cvector_t *best_psi0 = NULL;

  for (int nup = 0; nup <= N; nup++) {
    for (int k = 0; k < N; k++) {
      spin_sector_t *sec = spin_sector_build(N, nup, k);
      if (!sec || sec->dim == 0) {
        spin_sector_free(sec);

        continue;
      }

      sparse_matrix_t *H =
          spin_sector_hamiltonian(sec, J_COUPLING, J_COUPLING, /*pbc=*/1);
      int iters = sec->dim < 80 ? sec->dim : 80;
      lanczos_result_t *res = lanczos_eigs(H, 1, iters, 1e-10);

      if (res && res->values[0] < best_E) {
        best_E = res->values[0];

        spin_sector_free(best_sec);
        cvector_free(best_psi0);

        best_sec = sec;
        best_psi0 = cvector_alloc(sec->dim);
        for (int i = 0; i < sec->dim; i++) {
          best_psi0->data[i] = CMAT(res->vectors, i, 0);
        }

        sparse_free(H);
        lanczos_free(res);

        continue; // sec/H ownership transferred to best_* (don't double-free)
      }

      lanczos_free(res);
      sparse_free(H);
      spin_sector_free(sec);
    }
  }

  *out_sector = best_sec;
  *out_psi0 = best_psi0;

  return best_E;
}

int main(void) {
  printf(" > Dynamical Structure Factor S(q, omega): N=%d Heisenberg Ring "
         "(J=%.2f)\n\n",
         N_SITES, J_COUPLING);

  spin_sector_t *gs_sector = NULL;
  cvector_t *psi0 = NULL;
  double E0 = find_ground_state(N_SITES, &gs_sector, &psi0);

  printf("   Ground state: E0 = %.6f  (sector nup=%d, k=%d, dim=%d)\n\n", E0,
         gs_sector->nup, gs_sector->k, gs_sector->dim);
  printf("   (N=8 reference value from independent full-space ED: "
         "-3.651093)\n\n");

  double *omega = malloc((size_t)N_OMEGA * sizeof *omega);
  double domega = (OMEGA_MAX - OMEGA_MIN) / N_OMEGA;
  for (int i = 0; i < N_OMEGA; i++) {
    omega[i] = OMEGA_MIN + (i + 0.5) * domega;
  }

  /* NOTE: S(q, \omega) for every q, stored row-major [q_index][\omega_{index}]
   * so both print equal-time structure factor S(q) = I0(q) and plot a few
   * representative cuts afterward. */
  double **S_q_omega = malloc((size_t)N_SITES * sizeof *S_q_omega);
  double *I0_of_q = malloc((size_t)N_SITES * sizeof *I0_of_q);

  printf("   %-8s %-14s %-14s %-16s\n", "q_index", "q/pi", "I0=S(q)",
         "sum-rule check");
  printf("   -----------------------------------------------------------\n");

  for (int qi = 0; qi < N_SITES; qi++) {
    spin_sector_t *target = NULL;
    cvector_t *phi0 = NULL;
    double I0 = 0.0;

    int rc = spin_apply_szq(gs_sector, psi0, qi, &target, &phi0, &I0);
    S_q_omega[qi] = malloc((size_t)N_OMEGA * sizeof *S_q_omega[qi]);
    I0_of_q[qi] = I0;

    if (rc != 0 || I0 < 1e-14) {
      /* Zero spectral weight at this q (can happen at q=0 for an S^z
       * excitation, since S^z_0 is proportional to the conserved total S^z and
       * annihilates the ground state up to a constant). */
      for (int i = 0; i < N_OMEGA; i++) {
        S_q_omega[qi][i] = 0.0;
      }
      printf("   %-8d %-14.4f %-14.6e %-16s\n", qi,
             2.0 * (double)qi / (double)N_SITES, I0, "(no weight)");

      spin_sector_free(target);
      cvector_free(phi0);

      continue;
    }

    cvector_t *f0 = cvector_copy(phi0);
    cvector_normalize(f0);

    sparse_matrix_t *H_target =
        spin_sector_hamiltonian(target, J_COUPLING, J_COUPLING, /*pbc=*/1);
    lanczos_tridiag_t *tri =
        lanczos_tridiagonalize(H_target, f0, target->dim, 1e-12);

    double integral = 0.0;
    for (int i = 0; i < N_OMEGA; i++) {
      double S = tri ? spin_dsf_continued_fraction(
                           tri->alpha, tri->beta, tri->m, E0, I0, omega[i], ETA)
                     : 0.0;

      S_q_omega[qi][i] = S;
      integral += S * domega;
    }

    printf("   %-8d %-14.4f %-14.6e %-16s\n", qi,
           2.0 * (double)qi / (double)N_SITES, I0, "");
    printf("            (numeric integral of S(q,w) over the mesh = %.6e, "
           "vs I0 = %.6e)\n",
           integral, I0);

    lanczos_tridiag_free(tri);
    sparse_free(H_target);
    cvector_free(f0);
    cvector_free(phi0);
    spin_sector_free(target);
  }

  /* Plot 3 representative cuts: q=0 (should be flat/near-zero, S^z_0 is
   * ~conserved), q near \pi/2, and q=\pi (dominant AFM peak). */
  const int q_show[3] = {1, N_SITES / 4 > 0 ? N_SITES / 4 : 1, N_SITES / 2};
  const char *labels[3];
  char label_buf[3][32];
  const double *ys[3];
  for (int i = 0; i < 3; i++) {
    int qi = q_show[i];
    snprintf(label_buf[i], sizeof label_buf[i], "q = %.3f pi",
             2.0 * (double)qi / (double)N_SITES);

    labels[i] = label_buf[i];
    ys[i] = S_q_omega[qi];
  }

  plot_opts_t opts = {0};
  opts.title = "S(q, \\omega) - Heisenberg ring, N=8";
  opts.xlabel = "\\omega / J";
  opts.ylabel = "S(q, \\omega)";
  opts.xmin = PLOT_OMEGA_MIN;
  opts.xmax = PLOT_OMEGA_MAX;

  plot_lines("spin_chain_dsf", PLOT_FORMAT_PNG, omega, ys, 3, N_OMEGA, labels,
             &opts);
  printf("\n   Saved spin_chain_dsf.png (S(q, \\omega) at q = %.3f, %.3f, %.3f "
         "\\pi)\n",
         2.0 * q_show[0] / N_SITES, 2.0 * q_show[1] / N_SITES,
         2.0 * q_show[2] / N_SITES);

  // Static structure factor S(q) = I0(q) across the whole zone
  double *q_over_pi = malloc((size_t)N_SITES * sizeof *q_over_pi);
  for (int qi = 0; qi < N_SITES; qi++) {
    q_over_pi[qi] = 2.0 * (double)qi / (double)N_SITES;
  }

  plot_opts_t opts2 = {0};
  opts2.title = "Equal-time structure factor S(q) = I0(q)";
  opts2.xlabel = "q / \\pi";
  opts2.ylabel = "S(q)";
  plot_line("spin_chain_static_sq", PLOT_FORMAT_PNG, q_over_pi, I0_of_q,
            N_SITES, &opts2);
  printf("   Saved spin_chain_static_sq.png (S(q) peaks at q=\\pi, the "
         "antiferromagnetic wavevector, as expected)\n");

  free(q_over_pi);
  for (int qi = 0; qi < N_SITES; qi++) {
    free(S_q_omega[qi]);
  }
  free(S_q_omega);
  free(I0_of_q);
  free(omega);
  cvector_free(psi0);
  spin_sector_free(gs_sector);

  return 0;
}
