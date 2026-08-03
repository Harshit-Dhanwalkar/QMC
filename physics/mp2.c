/*
Second-order Moller-Plesset perturbation theory (MP2), s-orbitals-only
restricted.
*/

#include "mp2.h"
#include "hartree_fock.h"
#include <math.h>
#include <stdlib.h>

// (ia|jb) = int u_i(r) u_a(r) Y0_jb(r) dr
static double two_electron_integral(const double *r, int N, double dr,
                                    const double *ui, const double *ua,
                                    const double *uj, const double *ub,
                                    double *Y0_buf) {
  compute_Y0(r, N, dr, uj, ub, Y0_buf);

  double integral = 0.0;
  for (int k = 0; k < N; k++) {
    integral += ui[k] * ua[k] * Y0_buf[k] * dr;
  }

  return integral;
}

mp2_result_t mp2_correlation_energy(const hf_result_t *hf, const double *r,
                                    int N, int n_virtual) {
  mp2_result_t result = {0};

  if (!hf || !r || N < 10 || n_virtual < 1 || n_virtual > hf->n_virtual ||
      hf->n_orbitals < 1) {
    return result;
  }

  double dr = r[1] - r[0];
  if (dr <= 0.0) {
    return result;
  }

  int n_occ = hf->n_orbitals;

  // Extract raw double arrays once
  double **occ = malloc((size_t)n_occ * sizeof *occ);
  double **virt = malloc((size_t)n_virtual * sizeof *virt);
  double *Y0_buf = malloc((size_t)N * sizeof *Y0_buf);
  if (!occ || !virt || !Y0_buf) {
    free(occ);
    free(virt);
    free(Y0_buf);

    return result;
  }

  for (int i = 0; i < n_occ; i++) {
    occ[i] = malloc((size_t)N * sizeof *occ[i]);
    for (int k = 0; k < N; k++) {
      occ[i][k] = hf->orbitals[i]->data[k].re;
    }
  }

  for (int a = 0; a < n_virtual; a++) {
    virt[a] = malloc((size_t)N * sizeof *virt[a]);
    for (int k = 0; k < N; k++) {
      virt[a][k] = hf->virtual_orbitals[a]->data[k].re;
    }
  }

  // Recompute all (ia|jb)-relevant Y0 kernels (j,b) pair inside loop
  double e_mp2 = 0.0;

  for (int i = 0; i < n_occ; i++) {
    for (int j = 0; j < n_occ; j++) {
      for (int a = 0; a < n_virtual; a++) {
        for (int b = 0; b < n_virtual; b++) {
          double iajb = two_electron_integral(r, N, dr, occ[i], virt[a], occ[j],
                                              virt[b], Y0_buf);
          double ibja = two_electron_integral(r, N, dr, occ[i], virt[b], occ[j],
                                              virt[a], Y0_buf);

          double denom = hf->orbital_energies[i] + hf->orbital_energies[j] -
                         hf->virtual_energies[a] - hf->virtual_energies[b];

          if (fabs(denom) < 1e-12) {
            continue;
          }

          e_mp2 += iajb * (2.0 * iajb - ibja) / denom;
        }
      }
    }
  }

  for (int i = 0; i < n_occ; i++) {
    free(occ[i]);
  }
  for (int a = 0; a < n_virtual; a++) {
    free(virt[a]);
  }

  free(occ);
  free(virt);
  free(Y0_buf);

  result.e_hf = hf->total_energy;
  result.e_mp2 = e_mp2;
  result.e_total = hf->total_energy + e_mp2;
  result.n_occ = n_occ;
  result.n_virt = n_virtual;

  return result;
}
