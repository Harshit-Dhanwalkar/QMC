/*
Second-order Moller-Plesset perturbation theory (MP2), s-orbitals-only
restricted.
*/

#include "mp2.h"
#include "hartree_fock.h"
#include "molecular_integrals.h"
#include <math.h>
#include <stdlib.h>

// (ia|jb) = \int u_i(r) u_a(r) Y0_jb(r) dr
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

molecular_mp2_result_t molecular_mp2(int n_basis, const double *eri_mo,
                                     const double *mo_energy, int n_electrons,
                                     int n_frozen_spatial, double e_rhf) {
  molecular_mp2_result_t result = {0};

  if (n_basis <= 0 || !eri_mo || !mo_energy || n_electrons % 2 != 0 ||
      n_frozen_spatial < 0) {
    return result;
  }

  int n_occ_total = n_electrons / 2;
  if (n_occ_total <= n_frozen_spatial || n_occ_total > n_basis) {
    return result;
  }

  int n_virt = n_basis - n_occ_total;
  if (n_virt < 1) {
    return result;
  }

  double e_mp2 = 0.0;

  for (int i = n_frozen_spatial; i < n_occ_total; i++) {
    for (int j = n_frozen_spatial; j < n_occ_total; j++) {
      for (int a = n_occ_total; a < n_basis; a++) {
        for (int b = n_occ_total; b < n_basis; b++) {
          double iajb = MOLINT_ERI(eri_mo, n_basis, i, a, j, b);
          double ibja = MOLINT_ERI(eri_mo, n_basis, i, b, j, a);

          double denom =
              mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b];

          if (fabs(denom) < 1e-12) {
            continue;
          }

          e_mp2 += iajb * (2.0 * iajb - ibja) / denom;
        }
      }
    }
  }

  result.e_rhf = e_rhf;
  result.e_mp2 = e_mp2;
  result.e_total = e_rhf + e_mp2;
  result.n_occ = n_occ_total - n_frozen_spatial;
  result.n_virt = n_virt;

  return result;
}

molecular_ump2_result_t
molecular_ump2(int n_basis, const double *eri_aaaa, const double *eri_bbbb,
               const double *eri_aabb, const double *mo_energy_alpha,
               const double *mo_energy_beta, int n_alpha, int n_beta,
               int n_frozen_alpha, int n_frozen_beta, double e_uhf) {
  molecular_ump2_result_t result = {0};

  if (n_basis <= 0 || !eri_aaaa || !eri_bbbb || !eri_aabb || !mo_energy_alpha ||
      !mo_energy_beta || n_frozen_alpha < 0 || n_frozen_beta < 0 ||
      n_alpha < n_beta) {
    return result;
  }

  if (n_alpha <= n_frozen_alpha || n_alpha > n_basis ||
      (n_beta > 0 && n_beta <= n_frozen_beta) || n_beta > n_basis) {
    return result;
  }

  double e_aa = 0.0;
  for (int i = n_frozen_alpha; i < n_alpha; i++) {
    for (int j = n_frozen_alpha; j < n_alpha; j++) {
      for (int a = n_alpha; a < n_basis; a++) {
        for (int b = n_alpha; b < n_basis; b++) {
          double iajb = MOLINT_ERI(eri_aaaa, n_basis, i, a, j, b);
          double ibja = MOLINT_ERI(eri_aaaa, n_basis, i, b, j, a);
          double diff = iajb - ibja;

          double denom = mo_energy_alpha[i] + mo_energy_alpha[j] -
                         mo_energy_alpha[a] - mo_energy_alpha[b];

          if (fabs(denom) < 1e-12) {
            continue;
          }

          e_aa += 0.25 * diff * diff / denom;
        }
      }
    }
  }

  double e_bb = 0.0;
  for (int i = n_frozen_beta; i < n_beta; i++) {
    for (int j = n_frozen_beta; j < n_beta; j++) {
      for (int a = n_beta; a < n_basis; a++) {
        for (int b = n_beta; b < n_basis; b++) {
          double iajb = MOLINT_ERI(eri_bbbb, n_basis, i, a, j, b);
          double ibja = MOLINT_ERI(eri_bbbb, n_basis, i, b, j, a);
          double diff = iajb - ibja;

          double denom = mo_energy_beta[i] + mo_energy_beta[j] -
                         mo_energy_beta[a] - mo_energy_beta[b];

          if (fabs(denom) < 1e-12) {
            continue;
          }

          e_bb += 0.25 * diff * diff / denom;
        }
      }
    }
  }

  double e_ab = 0.0;
  for (int i = n_frozen_alpha; i < n_alpha; i++) {
    for (int j = n_frozen_beta; j < n_beta; j++) {
      for (int a = n_alpha; a < n_basis; a++) {
        for (int b = n_beta; b < n_basis; b++) {
          double iajb = MOLINT_ERI(eri_aabb, n_basis, i, a, j, b);

          double denom = mo_energy_alpha[i] + mo_energy_beta[j] -
                         mo_energy_alpha[a] - mo_energy_beta[b];

          if (fabs(denom) < 1e-12) {
            continue;
          }

          e_ab += iajb * iajb / denom;
        }
      }
    }
  }

  result.e_uhf = e_uhf;
  result.e_aa = e_aa;
  result.e_bb = e_bb;
  result.e_ab = e_ab;
  result.e_mp2 = e_aa + e_bb + e_ab;
  result.e_total = e_uhf + result.e_mp2;

  return result;
}
