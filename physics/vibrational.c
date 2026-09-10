/*
Numerical Hessian (central-difference of the analytic gradient) + harmonic
vibrational frequency analysis via mass-weighting and translation/rotation
projection.
*/

#include "vibrational.h"
#include "../core/complex.h"
#include "../core/linalg/eigen_generic.h"
#include "../core/matrix.h"
#include "physics/molecular_hf.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double *hessian_numerical_generic(molecular_geometry_fn build_system,
                                         molecular_system_free_fn free_system,
                                         void *userdata,
                                         const double geom0[][3], int n_atoms,
                                         double scf_tol, int scf_max_iter,
                                         double h, int use_uhf, int n_beta) {
  if (!build_system || !free_system || n_atoms <= 0 || h <= 0.0) {
    return NULL;
  }

  int dim = 3 * n_atoms;
  double *H = malloc((size_t)dim * dim * sizeof *H);
  if (!H) {
    return NULL;
  }

  double (*geom)[3] = malloc((size_t)n_atoms * sizeof *geom);
  if (!geom) {
    free(H);

    return NULL;
  }

  int ok = 1;

  for (int i = 0; i < dim && ok; i++) {
    int atom_i = i / 3, d_i = i % 3;

    for (int sign = 0; sign < 2 && ok; sign++) {
      memcpy(geom, geom0, (size_t)n_atoms * sizeof *geom);
      geom[atom_i][d_i] += (sign == 0) ? h : -h;

      molecular_system_t *sys = build_system(geom, n_atoms, userdata);
      if (!sys) {
        ok = 0;

        break;
      }

      double *grad = NULL;
      int converged = 0;

      if (!use_uhf) {
        molecular_hf_result_t *scf =
            molecular_rhf(sys->basis, sys->n_basis, sys->mol, sys->n_electrons,
                          scf_tol, scf_max_iter);
        if (scf) {
          converged = scf->converged;
          if (converged) {
            grad = molecular_rhf_gradient(sys->basis, sys->n_basis, sys->mol,
                                          sys->atom_of_basis, scf);
          }
          molecular_hf_result_free(scf);
        }
      } else {
        molecular_uhf_result_t *scf =
            molecular_uhf(sys->basis, sys->n_basis, sys->mol, sys->n_electrons,
                          n_beta, scf_tol, scf_max_iter);
        if (scf) {
          converged = scf->converged;
          if (converged) {
            grad = molecular_uhf_gradient(sys->basis, sys->n_basis, sys->mol,
                                          sys->atom_of_basis, scf);
          }
          molecular_uhf_result_free(scf);
        }
      }

      free_system(sys, userdata);

      if (!converged || !grad) {
        free(grad);
        ok = 0;

        break;
      }

      // Central difference: H[i][j] accumulates +grad[j]/(2h) for +h
      // displacement and -\grad[j]/(2h) for the -h displacement
      double factor = (sign == 0) ? (1.0 / (2.0 * h)) : (-1.0 / (2.0 * h));
      for (int j = 0; j < dim; j++) {
        if (sign == 0) {
          H[i * dim + j] = factor * grad[j];
        } else {
          H[i * dim + j] += factor * grad[j];
        }
      }

      free(grad);
    }
  }

  free(geom);

  if (!ok) {
    free(H);

    return NULL;
  }

  // Symmetrize: central differences leave a small numerical asymmetry
  for (int i = 0; i < dim; i++) {
    for (int j = i + 1; j < dim; j++) {
      double avg = 0.5 * (H[i * dim + j] + H[j * dim + i]);

      H[i * dim + j] = avg;
      H[j * dim + i] = avg;
    }
  }

  return H;
}

double *molecular_rhf_hessian_numerical(molecular_geometry_fn build_system,
                                        molecular_system_free_fn free_system,
                                        void *userdata, const double geom0[][3],
                                        int n_atoms, double scf_tol,
                                        int scf_max_iter, double h) {
  return hessian_numerical_generic(build_system, free_system, userdata, geom0,
                                   n_atoms, scf_tol, scf_max_iter, h, 0, 0);
}

double *molecular_uhf_hessian_numerical(molecular_geometry_fn build_system,
                                        molecular_system_free_fn free_system,
                                        void *userdata, const double geom0[][3],
                                        int n_atoms, int n_beta, double scf_tol,
                                        int scf_max_iter, double h) {
  return hessian_numerical_generic(build_system, free_system, userdata, geom0,
                                   n_atoms, scf_tol, scf_max_iter, h, 1,
                                   n_beta);
}

void molecular_vib_result_free(molecular_vib_result_t *res) {
  if (!res) {
    return;
  }

  free(res->frequencies_cm1);
  free(res);
}

// Standard atomic mass unit -> atomic unit of mass (electron mass) ratio.
#define AMU_TO_ME 1822.888486209
// 1 Hartree in wavenumbers (cm^-1): E(cm^-1) = E(Hartree) * this constant.
#define HARTREE_TO_CM1 219474.6313705

molecular_vib_result_t *
molecular_vibrational_analysis(const double *hessian, int n_atoms,
                               const double masses_amu[],
                               const double geom[][3]) {
  if (!hessian || n_atoms <= 0 || !masses_amu || !geom) {
    return NULL;
  }

  int dim = 3 * n_atoms;

  double *masses_au = malloc((size_t)n_atoms * sizeof *masses_au);
  double *sqm = malloc((size_t)n_atoms * sizeof *sqm);
  double *Mw = malloc((size_t)dim * dim * sizeof *Mw);
  if (!masses_au || !sqm || !Mw) {
    free(masses_au);
    free(sqm);
    free(Mw);

    return NULL;
  }

  for (int a = 0; a < n_atoms; a++) {
    masses_au[a] = masses_amu[a] * AMU_TO_ME;
    sqm[a] = sqrt(masses_au[a]);
  }

  for (int a = 0; a < n_atoms; a++) {
    for (int b = 0; b < n_atoms; b++) {
      double denom = sqm[a] * sqm[b];

      for (int d1 = 0; d1 < 3; d1++) {
        for (int d2 = 0; d2 < 3; d2++) {
          int i = 3 * a + d1, j = 3 * b + d2;

          Mw[i * dim + j] = hessian[i * dim + j] / denom;
        }
      }
    }
  }

  // Mass-weighted center of mass and translation/rotation vectors.
  double total_mass = 0.0, com[3] = {0, 0, 0};
  for (int a = 0; a < n_atoms; a++) {
    total_mass += masses_au[a];
    for (int d = 0; d < 3; d++) {
      com[d] += masses_au[a] * geom[a][d];
    }
  }
  for (int d = 0; d < 3; d++) {
    com[d] /= total_mass;
  }

  double *D = calloc((size_t)dim * 6, sizeof *D); // column-major, 6 columns
  if (!D) {
    free(masses_au);
    free(sqm);
    free(Mw);

    return NULL;
  }

  for (int a = 0; a < n_atoms; a++) {
    double x = geom[a][0] - com[0];
    double y = geom[a][1] - com[1];
    double z = geom[a][2] - com[2];

    // Translation: T_d[3a+d'] = sqm[a] if d'==d else 0.
    D[(3 * a + 0) * 6 + 0] = sqm[a];
    D[(3 * a + 1) * 6 + 1] = sqm[a];
    D[(3 * a + 2) * 6 + 2] = sqm[a];

    // Rotation about x/y/z: mass-weighted (r x e_axis).
    D[(3 * a + 0) * 6 + 3] = sqm[a] * 0.0;
    D[(3 * a + 1) * 6 + 3] = sqm[a] * (-z);
    D[(3 * a + 2) * 6 + 3] = sqm[a] * (y);

    D[(3 * a + 0) * 6 + 4] = sqm[a] * (z);
    D[(3 * a + 1) * 6 + 4] = sqm[a] * 0.0;
    D[(3 * a + 2) * 6 + 4] = sqm[a] * (-x);

    D[(3 * a + 0) * 6 + 5] = sqm[a] * (-y);
    D[(3 * a + 1) * 6 + 5] = sqm[a] * (x);
    D[(3 * a + 2) * 6 + 5] = sqm[a] * 0.0;
  }

  // NOTE: Gram-Schmidt orthonormalize the 6 columns of D, dropping any that are
  // (numerically) linearly dependent on the earlier ones
  double *Q = malloc((size_t)dim * 6 * sizeof *Q);
  int n_indep = 0;
  const double norm_tol = 1e-8;

  for (int c = 0; c < 6; c++) {
    double *col = malloc((size_t)dim * sizeof *col);
    for (int i = 0; i < dim; i++) {
      col[i] = D[i * 6 + c];
    }

    for (int k = 0; k < n_indep; k++) {
      double dot = 0.0;
      for (int i = 0; i < dim; i++) {
        dot += Q[i * 6 + k] * col[i];
      }
      for (int i = 0; i < dim; i++) {
        col[i] -= dot * Q[i * 6 + k];
      }
    }

    double norm2 = 0.0;
    for (int i = 0; i < dim; i++) {
      norm2 += col[i] * col[i];
    }
    double norm = sqrt(norm2);

    if (norm > norm_tol) {
      for (int i = 0; i < dim; i++) {
        Q[i * 6 + n_indep] = col[i] / norm;
      }
      n_indep++;
    }

    free(col);
  }

  free(D);

  // Project: Hproj = (I - Q Q^T) Mw (I - Q Q^T).
  double *Hproj = malloc((size_t)dim * dim * sizeof *Hproj);
  double *tmp = malloc((size_t)dim * dim * sizeof *tmp);
  if (!Q || !Hproj || !tmp) {
    free(masses_au);
    free(sqm);
    free(Mw);
    free(Q);
    free(Hproj);
    free(tmp);

    return NULL;
  }

  // tmp = Mw - Mw * Q * Q^T (apply projector on the right)
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      double v = Mw[i * dim + j];

      for (int k = 0; k < n_indep; k++) {
        double mwq_ik = 0.0;
        for (int m = 0; m < dim; m++) {
          mwq_ik += Mw[i * dim + m] * Q[m * 6 + k];
        }

        v -= mwq_ik * Q[j * 6 + k];
      }

      tmp[i * dim + j] = v;
    }
  }

  // Hproj = tmp - Q*Q^T*tmp (apply projector on the left)
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      double v = tmp[i * dim + j];

      for (int k = 0; k < n_indep; k++) {
        double qtmp_kj = 0.0;
        for (int m = 0; m < dim; m++) {
          qtmp_kj += Q[m * 6 + k] * tmp[m * dim + j];
        }

        v -= Q[i * 6 + k] * qtmp_kj;
      }

      Hproj[i * dim + j] = v;
    }
  }

  free(tmp);
  free(Q);
  free(Mw);
  free(masses_au);
  free(sqm);

  // Diagonalize the (real, symmetric) projected Hessian via generic Hermitian
  // eigensolver (imaginary parts all zero)
  cmatrix_t *Hc = cmatrix_alloc(dim, dim);
  if (!Hc) {
    free(Hproj);

    return NULL;
  }

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      CMAT(Hc, i, j) = c_new(Hproj[i * dim + j], 0.0);
    }
  }

  free(Hproj);

  eigen_t *eig = cmatrix_eigh_generic(Hc);
  cmatrix_free(Hc);
  if (!eig) {
    return NULL;
  }

  int n_trans_rot = n_indep;
  int n_modes = dim - n_trans_rot;
  if (n_modes < 0) {
    n_modes = 0;
  }

  molecular_vib_result_t *result = malloc(sizeof *result);
  double *freqs = malloc((size_t)(n_modes > 0 ? n_modes : 1) * sizeof *freqs);
  if (!result || !freqs) {
    free(result);
    free(freqs);
    eigen_free(eig);

    return NULL;
  }

  // NOTE: eig->eigenvalues includes n_trans_rot values that are (numerically)
  // ~0 from the projected-out subspace; skip the n_trans_rot smallest in
  // magnitude and report the rest, sorted ascending (eigen_t already returns
  // eigenvalues in ascending order for a Hermitian solve, but the near-zero
  // trans/rot ones can land anywhere near the low end depending on residual
  // numerical noise, so explicitly pick by |value| rank)
  int *order = malloc((size_t)dim * sizeof *order);
  for (int i = 0; i < dim; i++) {
    order[i] = i;
  }

  // simple insertion sort by |eigenvalue|, ascending
  for (int i = 1; i < dim; i++) {
    int key = order[i];
    double key_abs = fabs(eig->eigenvalues[key]);
    int j = i - 1;

    while (j >= 0 && fabs(eig->eigenvalues[order[j]]) > key_abs) {
      order[j + 1] = order[j];
      j--;
    }

    order[j + 1] = key;
  }

  // NOTE: The n_trans_rot smallest-|value| eigenvalues are the projected-out
  // modes; the remaining n_modes are genuine vibrations. Collect those, then
  // sort them ascending by signed value (imaginary/negative first).
  double *vib_vals =
      malloc((size_t)(n_modes > 0 ? n_modes : 1) * sizeof *vib_vals);
  for (int i = 0; i < n_modes; i++) {
    vib_vals[i] = eig->eigenvalues[order[n_trans_rot + i]];
  }

  free(order);

  for (int i = 1; i < n_modes; i++) {
    double key = vib_vals[i];
    int j = i - 1;

    while (j >= 0 && vib_vals[j] > key) {
      vib_vals[j + 1] = vib_vals[j];
      j--;
    }

    vib_vals[j + 1] = key;
  }

  for (int i = 0; i < n_modes; i++) {
    double lam = vib_vals[i];

    freqs[i] = (lam >= 0.0) ? sqrt(lam) * HARTREE_TO_CM1
                            : -sqrt(-lam) * HARTREE_TO_CM1;
  }

  free(vib_vals);
  eigen_free(eig);

  result->n_atoms = n_atoms;
  result->n_modes = n_modes;
  result->n_trans_rot = n_trans_rot;
  result->frequencies_cm1 = freqs;

  return result;
}
