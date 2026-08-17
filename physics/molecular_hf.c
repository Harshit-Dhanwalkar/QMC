#include "molecular_hf.h"
#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "molecular_integrals.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

molecular_hf_result_t *molecular_rhf(basis_function_t **basis, int n_basis,
                                     const molecule_t *mol, int n_electrons,
                                     double tol, int max_iter) {
  if (!basis || n_basis <= 0 || !mol || n_electrons <= 0 ||
      n_electrons % 2 != 0 || n_electrons > 2 * n_basis || max_iter < 1) {
    return NULL;
  }

  int n_occ = n_electrons / 2;
  int n = n_basis;

  cmatrix_t *S = molecular_overlap_matrix(basis, n);
  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, n, mol);
  double *eri = molecular_eri_tensor(basis, n);
  if (!S || !Hcore || !eri) {
    cmatrix_free(S);
    cmatrix_free(Hcore);
    free(eri);

    return NULL;
  }

  // Lowdin symmetric orthogonalization:
  // S^{-1/2} via S's own eigendecomposition
  eigen_t *eig_S = cmatrix_eigh_complex(S);
  cmatrix_t *Shalf = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;

      for (int k = 0; k < n; k++) {
        double inv_sqrt = 1.0 / sqrt(eig_S->eigenvalues[k]);

        sum += CMAT(eig_S->eigenvectors, i, k).re * inv_sqrt *
               CMAT(eig_S->eigenvectors, j, k).re;
      }

      CMAT(Shalf, i, j) = c_real(sum);
    }
  }

  eigen_free(eig_S);

  double *D = calloc((size_t)n * n, sizeof(double));
  double *D_new = malloc((size_t)n * n * sizeof(double));
  double *C_arr = malloc((size_t)n * n * sizeof(double));
  double *orbital_energies = malloc((size_t)n * sizeof(double));

  double E_elec = 0.0, E_elec_old = 0.0;
  int converged = 0;
  int iter = 0;

  for (; iter < max_iter; iter++) {
    cmatrix_t *F = cmatrix_alloc(n, n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double v = CMAT(Hcore, i, j).re;

        for (int k = 0; k < n; k++) {
          for (int m = 0; m < n; m++) {
            v += D[k * n + m] * (MOLINT_ERI(eri, n, i, j, m, k) -
                                 0.5 * MOLINT_ERI(eri, n, i, k, m, j));
          }
        }

        CMAT(F, i, j) = c_real(v);
      }
    }

    // F' = Shalf^T F Shalf
    cmatrix_t *tmp = cmatrix_multiply(Shalf, F);
    cmatrix_t *Fp = cmatrix_multiply(tmp, Shalf);
    cmatrix_free(tmp);
    cmatrix_free(F);

    eigen_t *eig_F = cmatrix_eigh_complex(Fp);
    cmatrix_free(Fp);

    // C = Shalf * C'
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        double v = 0.0;

        for (int p = 0; p < n; p++) {
          v += CMAT(Shalf, i, p).re * CMAT(eig_F->eigenvectors, p, k).re;
        }

        C_arr[i * n + k] = v;
      }
    }

    memcpy(orbital_energies, eig_F->eigenvalues, (size_t)n * sizeof(double));

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double v = 0.0;

        for (int k = 0; k < n_occ; k++) {
          v += 2.0 * C_arr[i * n + k] * C_arr[j * n + k];
        }

        D_new[i * n + j] = v;
      }
    }

    E_elec = 0.0;
    /* NOTE: Recompute F once more at D_new for energy expression, cheaper to
     * just reuse F built above evaluated with D (pre-update), using RHF energy
     * formula E = (1/2) * \sum D_new_ij (Hcore_ij + F_ij) with F built from
     * D_new for a fully self-consistent value
     */
    cmatrix_t *F_new = cmatrix_alloc(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double v = CMAT(Hcore, i, j).re;

        for (int k = 0; k < n; k++) {
          for (int m = 0; m < n; m++) {
            v += D_new[k * n + m] * (MOLINT_ERI(eri, n, i, j, m, k) -
                                     0.5 * MOLINT_ERI(eri, n, i, k, m, j));
          }
        }

        CMAT(F_new, i, j) = c_real(v);
        E_elec += 0.5 * D_new[i * n + j] * (CMAT(Hcore, i, j).re + v);
      }
    }

    cmatrix_free(F_new);

    memcpy(D, D_new, (size_t)n * n * sizeof(double));
    eigen_free(eig_F);

    if (iter > 0 && fabs(E_elec - E_elec_old) < tol) {
      converged = 1;
      iter++;

      break;
    }

    E_elec_old = E_elec;
  }

  double E_nuc = molecule_nuclear_repulsion(mol);

  molecular_hf_result_t *res = malloc(sizeof(molecular_hf_result_t));
  res->n_basis = n;
  res->n_electrons = n_electrons;
  res->electronic_energy = E_elec;
  res->total_energy = E_elec + E_nuc;
  res->iterations = iter;
  res->converged = converged;
  res->orbital_energies = orbital_energies;

  res->C = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(res->C, i, j) = c_real(C_arr[i * n + j]);
    }
  }

  free(D);
  free(D_new);
  free(C_arr);
  cmatrix_free(S);
  cmatrix_free(Hcore);
  cmatrix_free(Shalf);
  free(eri);

  return res;
}

void molecular_hf_result_free(molecular_hf_result_t *res) {
  if (!res) {
    return;
  }

  cmatrix_free(res->C);
  free(res->orbital_energies);
  free(res);
}

molecular_uhf_result_t *molecular_uhf(basis_function_t **basis, int n_basis,
                                      const molecule_t *mol, int n_alpha,
                                      int n_beta, double tol, int max_iter) {
  if (!basis || n_basis <= 0 || !mol || n_alpha < 0 || n_beta < 0 ||
      n_alpha + n_beta <= 0 || n_alpha > n_basis || n_beta > n_basis ||
      max_iter < 1) {
    return NULL;
  }

  int n = n_basis;

  cmatrix_t *S = molecular_overlap_matrix(basis, n);
  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, n, mol);
  double *eri = molecular_eri_tensor(basis, n);
  if (!S || !Hcore || !eri) {
    cmatrix_free(S);
    cmatrix_free(Hcore);
    free(eri);

    return NULL;
  }

  // Lowdin symmetric orthogonalization (shared by both spin channels)
  eigen_t *eig_S = cmatrix_eigh_complex(S);
  cmatrix_t *Shalf = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;

      for (int k = 0; k < n; k++) {
        double inv_sqrt = 1.0 / sqrt(eig_S->eigenvalues[k]);

        sum += CMAT(eig_S->eigenvectors, i, k).re * inv_sqrt *
               CMAT(eig_S->eigenvectors, j, k).re;
      }

      CMAT(Shalf, i, j) = c_real(sum);
    }
  }

  eigen_free(eig_S);

  double *Da = calloc((size_t)n * n, sizeof(double));
  double *Db = calloc((size_t)n * n, sizeof(double));
  double *Da_new = malloc((size_t)n * n * sizeof(double));
  double *Db_new = malloc((size_t)n * n * sizeof(double));
  double *Ca_arr = malloc((size_t)n * n * sizeof(double));
  double *Cb_arr = malloc((size_t)n * n * sizeof(double));
  double *orb_a = malloc((size_t)n * sizeof(double));
  double *orb_b = malloc((size_t)n * sizeof(double));

  double E_elec = 0.0, E_elec_old = 0.0;
  int converged = 0;
  int iter = 0;

  for (; iter < max_iter; iter++) {
    // Coulomb term shared by both spins, built from the total density
    double *Dtot = malloc((size_t)n * n * sizeof(double));

    for (int i = 0; i < n * n; i++) {
      Dtot[i] = Da[i] + Db[i];
    }

    cmatrix_t *Fa = cmatrix_alloc(n, n);
    cmatrix_t *Fb = cmatrix_alloc(n, n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double J = 0.0, Ka = 0.0, Kb = 0.0;

        for (int k = 0; k < n; k++) {
          for (int m = 0; m < n; m++) {
            double eri_coul = MOLINT_ERI(eri, n, i, j, m, k);
            double eri_exch = MOLINT_ERI(eri, n, i, k, m, j);

            J += Dtot[k * n + m] * eri_coul;
            Ka += Da[k * n + m] * eri_exch;
            Kb += Db[k * n + m] * eri_exch;
          }
        }

        CMAT(Fa, i, j) = c_real(CMAT(Hcore, i, j).re + J - Ka);
        CMAT(Fb, i, j) = c_real(CMAT(Hcore, i, j).re + J - Kb);
      }
    }

    free(Dtot);

    // Diagonalize each spin channel: F' = Shalf^T F Shalf, C = Shalf C'
    cmatrix_t *tmpA = cmatrix_multiply(Shalf, Fa);
    cmatrix_t *FpA = cmatrix_multiply(tmpA, Shalf);
    cmatrix_free(tmpA);
    cmatrix_free(Fa);

    cmatrix_t *tmpB = cmatrix_multiply(Shalf, Fb);
    cmatrix_t *FpB = cmatrix_multiply(tmpB, Shalf);
    cmatrix_free(tmpB);
    cmatrix_free(Fb);

    eigen_t *eig_A = cmatrix_eigh_complex(FpA);
    eigen_t *eig_B = cmatrix_eigh_complex(FpB);
    cmatrix_free(FpA);
    cmatrix_free(FpB);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        double va = 0.0, vb = 0.0;

        for (int p = 0; p < n; p++) {
          va += CMAT(Shalf, i, p).re * CMAT(eig_A->eigenvectors, p, k).re;
          vb += CMAT(Shalf, i, p).re * CMAT(eig_B->eigenvectors, p, k).re;
        }

        Ca_arr[i * n + k] = va;
        Cb_arr[i * n + k] = vb;
      }
    }

    memcpy(orb_a, eig_A->eigenvalues, (size_t)n * sizeof(double));
    memcpy(orb_b, eig_B->eigenvalues, (size_t)n * sizeof(double));
    eigen_free(eig_A);
    eigen_free(eig_B);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double va = 0.0, vb = 0.0;

        for (int k = 0; k < n_alpha; k++) {
          va += Ca_arr[i * n + k] * Ca_arr[j * n + k];
        }

        for (int k = 0; k < n_beta; k++) {
          vb += Cb_arr[i * n + k] * Cb_arr[j * n + k];
        }

        Da_new[i * n + j] = va;
        Db_new[i * n + j] = vb;
      }
    }

    // Recompute Fock matrices at the new densities for a fully self-consistent
    // energy expression
    double *Dtot_new = malloc((size_t)n * n * sizeof(double));
    for (int i = 0; i < n * n; i++) {
      Dtot_new[i] = Da_new[i] + Db_new[i];
    }

    E_elec = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double J = 0.0, Ka = 0.0, Kb = 0.0;

        for (int k = 0; k < n; k++) {
          for (int m = 0; m < n; m++) {
            double eri_coul = MOLINT_ERI(eri, n, i, j, m, k);
            double eri_exch = MOLINT_ERI(eri, n, i, k, m, j);

            J += Dtot_new[k * n + m] * eri_coul;
            Ka += Da_new[k * n + m] * eri_exch;
            Kb += Db_new[k * n + m] * eri_exch;
          }
        }

        double Fa_ij = CMAT(Hcore, i, j).re + J - Ka;
        double Fb_ij = CMAT(Hcore, i, j).re + J - Kb;

        E_elec += 0.5 * (Da_new[i * n + j] * (CMAT(Hcore, i, j).re + Fa_ij) +
                         Db_new[i * n + j] * (CMAT(Hcore, i, j).re + Fb_ij));
      }
    }

    free(Dtot_new);

    memcpy(Da, Da_new, (size_t)n * n * sizeof(double));
    memcpy(Db, Db_new, (size_t)n * n * sizeof(double));

    if (iter > 0 && fabs(E_elec - E_elec_old) < tol) {
      converged = 1;
      iter++;

      break;
    }

    E_elec_old = E_elec;
  }

  // <S^2> diagnostic: Sz(Sz+1) + n_beta - sum_{i occ a, j occ b} |<a_i|b_j>|^2
  // (AO overlap S sandwiched between occupied \alpha / \beta MO coefficients)
  double sz = 0.5 * (n_alpha - n_beta);
  double overlap_sum = 0.0;
  for (int i = 0; i < n_alpha; i++) {
    for (int j = 0; j < n_beta; j++) {
      double ov = 0.0;

      for (int p = 0; p < n; p++) {
        for (int q = 0; q < n; q++) {
          ov += Ca_arr[p * n + i] * CMAT(S, p, q).re * Cb_arr[q * n + j];
        }
      }

      overlap_sum += ov * ov;
    }
  }

  double spin_squared = sz * (sz + 1.0) + n_beta - overlap_sum;
  double E_nuc = molecule_nuclear_repulsion(mol);

  molecular_uhf_result_t *res = malloc(sizeof(molecular_uhf_result_t));
  res->n_basis = n;
  res->n_alpha = n_alpha;
  res->n_beta = n_beta;
  res->electronic_energy = E_elec;
  res->total_energy = E_elec + E_nuc;
  res->iterations = iter;
  res->converged = converged;
  res->spin_squared = spin_squared;
  res->orbital_energies_alpha = orb_a;
  res->orbital_energies_beta = orb_b;

  res->C_alpha = cmatrix_alloc(n, n);
  res->C_beta = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(res->C_alpha, i, j) = c_real(Ca_arr[i * n + j]);
      CMAT(res->C_beta, i, j) = c_real(Cb_arr[i * n + j]);
    }
  }

  free(Da);
  free(Db);
  free(Da_new);
  free(Db_new);
  free(Ca_arr);
  free(Cb_arr);
  cmatrix_free(S);
  cmatrix_free(Hcore);
  cmatrix_free(Shalf);
  free(eri);

  return res;
}

void molecular_uhf_result_free(molecular_uhf_result_t *res) {
  if (!res) {
    return;
  }

  cmatrix_free(res->C_alpha);
  cmatrix_free(res->C_beta);
  free(res->orbital_energies_alpha);
  free(res->orbital_energies_beta);
  free(res);
}

void molecular_ao_to_mo(const cmatrix_t *h_ao, const double *eri_ao,
                        const cmatrix_t *C, int n_basis, double *h_mo,
                        double *eri_mo) {
  int n = n_basis;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double v = 0.0;

      for (int p = 0; p < n; p++) {
        for (int q = 0; q < n; q++) {
          v += CMAT(C, p, i).re * CMAT(h_ao, p, q).re * CMAT(C, q, j).re;
        }
      }

      h_mo[i * n + j] = v;
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
          double v = 0.0;

          for (int p = 0; p < n; p++) {
            double cip = CMAT(C, p, i).re;

            if (cip == 0.0) {
              continue;
            }

            for (int q = 0; q < n; q++) {
              double cjq = CMAT(C, q, j).re;

              if (cjq == 0.0) {
                continue;
              }

              for (int r = 0; r < n; r++) {
                double ckr = CMAT(C, r, k).re;

                if (ckr == 0.0) {
                  continue;
                }

                for (int s = 0; s < n; s++) {
                  v += cip * cjq * ckr * CMAT(C, s, l).re *
                       MOLINT_ERI(eri_ao, n, p, q, r, s);
                }
              }
            }
          }

          MOLINT_ERI(eri_mo, n, i, j, k, l) = v;
        }
      }
    }
  }
}

/*
 * NOTE: Analytic RHF nuclear gradient (Reference :Pulay 1969 / Szabo & Ostlund
 * eq. 3.184; Helgaker, Jorgensen & Olsen ch. 10).
 *
 *   dE/dX_A = dE_nuc / dX_A
 *           + \sum_{uv} P_uv dH_uv / dX_A
 *           + (1/2) * \sum_{uvls} \Gamma_{uvls} d(uv|ls) / dX_A
 *           - \sum_{uv} W_{uv} dS_{uv} / dX_A
 *
 * Where
 *  \Gamma_{uvls} = P_{uv} * P_{ls} - 0.5 * P_{ul} * P_{vs} is RHF two-particle
 *      density
 *  P is the one-electron density matrix
 *  W_{uv} = 2 * \sum_{occ} eps_i C_{ui} C_{vi} is the energy-weighted density
 *      matrix (Pulay trick stationarity condition makes the implicit dP/dX
 *      contributions vanish against this term, neither P nor W nor C need to be
 *      re-differentiated wrt nuclear position).
 *
 * NOTE: dH_uv/dX_A (core Hamiltonian = kinetic + nuclear attraction) picks up a
 * contribution from every relevant piece: if AO u is centered on atom A, its
 * own-center kinetic/nuclear-attraction derivative (summed over every atom's
 * nuclear charge, since AOs interact with all nuclei); likewise for v; and
 * separately, for every atom's own nuclear-attraction operator, a
 * "Hellmann-Feynman" contribution from differentiating wrt THAT atom's own
 * charge position (gto_nuclear_attraction_grad_C), present regardless of which
 * atoms u,v happen to be centered on.
 */
double *molecular_rhf_gradient(basis_function_t **basis, int n_basis,
                               const molecule_t *mol, const int *atom_of_basis,
                               const molecular_hf_result_t *scf) {
  if (!basis || n_basis <= 0 || !mol || !atom_of_basis || !scf || !scf->C ||
      !scf->orbital_energies || !scf->converged) {
    return NULL;
  }

  int n = n_basis;
  int n_occ = scf->n_electrons / 2;
  int n_atoms = mol->n_atoms;

  double *grad = calloc((size_t)n_atoms * 3, sizeof(double));
  if (!grad) {
    return NULL;
  }

  /* One-electron density P_{uv} = 2 \sum_occ C_{ui} C_{vi}, and energy-weighted
   * density W_{uv} = 2 * \sum_{occ} eps_i C_ui C_vi.
   * Real orbitals throughout (molecular_rhf always produces real C/eps), so
   * only .re is used. */
  double *P = malloc((size_t)n * n * sizeof(double));
  double *W = malloc((size_t)n * n * sizeof(double));
  if (!P || !W) {
    free(P);
    free(W);
    free(grad);

    return NULL;
  }

  for (int u = 0; u < n; u++) {
    for (int v = 0; v < n; v++) {
      double p_uv = 0.0, w_uv = 0.0;

      for (int i = 0; i < n_occ; i++) {
        double cu = CMAT(scf->C, u, i).re;
        double cv = CMAT(scf->C, v, i).re;

        p_uv += cu * cv;
        w_uv += scf->orbital_energies[i] * cu * cv;
      }

      P[u * n + v] = 2.0 * p_uv;
      W[u * n + v] = 2.0 * w_uv;
    }
  }

  // Nuclear-nuclear repulsion: closed form
  for (int A = 0; A < n_atoms; A++) {
    for (int B = 0; B < n_atoms; B++) {
      if (B == A) {
        continue;
      }

      const double Rvec[3] = {mol->center[A][0] - mol->center[B][0],
                              mol->center[A][1] - mol->center[B][1],
                              mol->center[A][2] - mol->center[B][2]};
      double dist2 = Rvec[0] * Rvec[0] + Rvec[1] * Rvec[1] + Rvec[2] * Rvec[2];
      double dist3 = dist2 * sqrt(dist2);
      double pref = -mol->charge[A] * mol->charge[B] / dist3;

      for (int d = 0; d < 3; d++) {
        grad[3 * A + d] += pref * Rvec[d];
      }
    }
  }

  // One-electron (kinetic + nuclear-attraction) and overlap contributions
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < n; v++) {
      double p_uv = P[u * n + v];
      double w_uv = W[u * n + v];
      if (p_uv == 0.0 && w_uv == 0.0) {
        continue;
      }

      int atom_u = atom_of_basis[u];
      int atom_v = atom_of_basis[v];

      double gT[3], gS[3];
      if (atom_u >= 0) {
        gto_kinetic_grad_a(basis[u], basis[v], gT);
        gto_overlap_grad_a(basis[u], basis[v], gS);

        for (int d = 0; d < 3; d++) {
          grad[3 * atom_u + d] += p_uv * gT[d];
          grad[3 * atom_u + d] -= w_uv * gS[d];
        }

        for (int A = 0; A < n_atoms; A++) {
          double gV[3];
          gto_nuclear_attraction_grad_a(basis[u], basis[v], mol->center[A], gV);

          for (int d = 0; d < 3; d++) {
            grad[3 * atom_u + d] += p_uv * (-mol->charge[A]) * gV[d];
          }
        }
      }

      if (atom_v >= 0) {
        gto_kinetic_grad_a(basis[v], basis[u], gT);
        gto_overlap_grad_a(basis[v], basis[u], gS);

        for (int d = 0; d < 3; d++) {
          grad[3 * atom_v + d] += p_uv * gT[d];
          grad[3 * atom_v + d] -= w_uv * gS[d];
        }

        for (int A = 0; A < n_atoms; A++) {
          double gV[3];
          gto_nuclear_attraction_grad_a(basis[v], basis[u], mol->center[A], gV);

          for (int d = 0; d < 3; d++) {
            grad[3 * atom_v + d] += p_uv * (-mol->charge[A]) * gV[d];
          }
        }
      }

      // Hellmann-Feynman term: every atom's own nuclear-attraction operator
      // contributes here regardless of where u,v are centered
      for (int A = 0; A < n_atoms; A++) {
        double gC[3];
        gto_nuclear_attraction_grad_C(basis[u], basis[v], mol->center[A], gC);

        for (int d = 0; d < 3; d++) {
          grad[3 * A + d] += p_uv * (-mol->charge[A]) * gC[d];
        }
      }
    }
  }

  // Two-electron (ERI) contribution: 1/2 * \Gamma_{uvls} * d(uv|ls) / dX_A,
  // \Gamma_{uvls} = P_{uv} * P_{ls} - 1/2 * P_{ul} * P_{vs} (RHF two-particle
  // density)
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < n; v++) {
      for (int l = 0; l < n; l++) {
        for (int s = 0; s < n; s++) {
          double gamma = 0.5 * (P[u * n + v] * P[l * n + s] -
                                0.5 * P[u * n + l] * P[v * n + s]);
          if (gamma == 0.0) {
            continue;
          }

          const int atoms4[4] = {atom_of_basis[u], atom_of_basis[v],
                                 atom_of_basis[l], atom_of_basis[s]};
          basis_function_t *bfs4[4] = {basis[u], basis[v], basis[l], basis[s]};

          for (int slot = 0; slot < 4; slot++) {
            int A = atoms4[slot];
            if (A < 0) {
              continue;
            }

            double gE[3];
            /* NOTE: Reorder so `slot` becomes the bra-first argument,
             * preserving bra-pair/ket-pair structure. */
            switch (slot) {
            case 0:
              gto_eri_grad_a(bfs4[0], bfs4[1], bfs4[2], bfs4[3], gE);
              break;
            case 1:
              gto_eri_grad_a(bfs4[1], bfs4[0], bfs4[2], bfs4[3], gE);
              break;
            case 2:
              gto_eri_grad_a(bfs4[2], bfs4[3], bfs4[0], bfs4[1], gE);
              break;
            default:
              gto_eri_grad_a(bfs4[3], bfs4[2], bfs4[0], bfs4[1], gE);
              break;
            }

            for (int d = 0; d < 3; d++) {
              grad[3 * A + d] += gamma * gE[d];
            }
          }
        }
      }
    }
  }

  free(P);
  free(W);

  return grad;
}
