/*
Test: Jordan-Wigner fermion-to-qubit mapping.

1. Anticommutation relations {a_i, a_j^\dagger} = \delta_ij*I,
  {a_i, a_j} = 0 must hold exactly for JW operators, for every pair (i,j) in a
  4-mode system : these are the defining algebraic relations of fermionic
  operators, and whole point of Z-string construction is to make bosonic qubit
  tensor products satisfy them.
2. The number operator a_j^dagger a_j must correctly read off mode j's
  occupation on a fixed test state.
3. Cross-validation: a Hamiltonian built by composing
  jw_creation_operator/jw_annihilation_operator (tensor-product construction)
  must exactly match (to machine precision) the same hamiltonian built by
  second_quant_build_hopping_hamiltonian (a completely independent direct
  bit-manipulation + fermionic-sign-counting construction, no tensor products at
  all) and their eigenvalue spectra must agree.
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/second_quant.h"
#include <math.h>
#include <stdio.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);

  if (err > tol) {
    printf("  FAIL: %s\n", label);

    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

static void test_anticommutation(void) {
  printf("test_anticommutation:\n");

  int n_modes = 4;
  int dim = 1 << n_modes;

  cmatrix_t *a[4];
  cmatrix_t *a_dag[4];
  for (int j = 0; j < n_modes; j++) {
    a[j] = jw_annihilation_operator(j, n_modes);
    a_dag[j] = jw_creation_operator(j, n_modes);
  }

  double max_err = 0.0;

  for (int i = 0; i < n_modes; i++) {
    for (int j = 0; j < n_modes; j++) {
      // {a_i, a_j^dagger} = a_i*a_j^dagger + a_j^dagger*a_i
      cmatrix_t *p1 = cmatrix_multiply(a[i], a_dag[j]);
      cmatrix_t *p2 = cmatrix_multiply(a_dag[j], a[i]);
      cmatrix_t *sum1 = cmatrix_add(p1, p2);

      for (int r = 0; r < dim; r++) {
        for (int c = 0; c < dim; c++) {
          complex_t expected = (i == j && r == c) ? c_real(1.0) : c_zero();
          complex_t got = CMAT(sum1, r, c);
          double err = c_abs2(c_sub(got, expected));
          if (err > max_err) {
            max_err = err;
          }
        }
      }

      cmatrix_free(p1);
      cmatrix_free(p2);
      cmatrix_free(sum1);

      // {a_i, a_j} = 0
      cmatrix_t *q1 = cmatrix_multiply(a[i], a[j]);
      cmatrix_t *q2 = cmatrix_multiply(a[j], a[i]);
      cmatrix_t *sum2 = cmatrix_add(q1, q2);

      for (int k = 0; k < dim * dim; k++) {
        double err = c_abs2(sum2->data[k]);
        if (err > max_err) {
          max_err = err;
        }
      }

      cmatrix_free(q1);
      cmatrix_free(q2);
      cmatrix_free(sum2);
    }
  }

  check_close(
      max_err, 0.0, 1e-20,
      "max squared error across all {a_i,a_j^\\dagger} and {a_i,a_j} checks");

  for (int j = 0; j < n_modes; j++) {
    cmatrix_free(a[j]);
    cmatrix_free(a_dag[j]);
  }
}

static void test_number_operator(void) {
  printf("test_number_operator:\n");

  int n_modes = 4;

  // Test state |0011> (modes 2,3 occupied): index 0b0011 = 3
  int idx = 3;

  for (int j = 0; j < n_modes; j++) {
    cmatrix_t *aj = jw_annihilation_operator(j, n_modes);
    cmatrix_t *aj_dag = jw_creation_operator(j, n_modes);
    cmatrix_t *n_op = cmatrix_multiply(aj_dag, aj);

    double expected = (j == 2 || j == 3) ? 1.0 : 0.0;
    char label[32];
    snprintf(label, sizeof label, "<n_%d> for |0011>", j);
    check_close(CMAT(n_op, idx, idx).re, expected, 1e-12, label);

    cmatrix_free(aj);
    cmatrix_free(aj_dag);
    cmatrix_free(n_op);
  }
}

static void test_jw_vs_direct_hamiltonian(void) {
  printf("test_jw_vs_direct_hamiltonian:\n");

  int n_modes = 4;
  int dim = 1 << n_modes;
  double t = 1.0, U = 0.7;
  const double epsilon[4] = {0.3, -0.1, 0.5, 0.2};

  // Build H via JW operator composition
  cmatrix_t *a[4];
  cmatrix_t *a_dag[4];
  for (int j = 0; j < n_modes; j++) {
    a[j] = jw_annihilation_operator(j, n_modes);
    a_dag[j] = jw_creation_operator(j, n_modes);
  }

  cmatrix_t *H_jw = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H_jw->data[i] = c_zero();
  }

  for (int i = 0; i < n_modes; i++) {
    cmatrix_t *n_i = cmatrix_multiply(a_dag[i], a[i]);
    cmatrix_scale(n_i, c_real(epsilon[i]));
    cmatrix_t *tmp = cmatrix_add(H_jw, n_i);
    cmatrix_free(H_jw);
    cmatrix_free(n_i);

    H_jw = tmp;
  }

  for (int i = 0; i + 1 < n_modes; i++) {
    // -t * (a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i)
    cmatrix_t *hop1 = cmatrix_multiply(a_dag[i], a[i + 1]);
    cmatrix_t *hop2 = cmatrix_multiply(a_dag[i + 1], a[i]);
    cmatrix_t *hop_sum = cmatrix_add(hop1, hop2);
    cmatrix_scale(hop_sum, c_real(-t));
    cmatrix_t *tmp = cmatrix_add(H_jw, hop_sum);
    cmatrix_free(hop1);
    cmatrix_free(hop2);
    cmatrix_free(hop_sum);
    cmatrix_free(H_jw);

    H_jw = tmp;

    // U*n_i*n_{i+1}
    cmatrix_t *n_i = cmatrix_multiply(a_dag[i], a[i]);
    cmatrix_t *n_ip1 = cmatrix_multiply(a_dag[i + 1], a[i + 1]);
    cmatrix_t *interaction = cmatrix_multiply(n_i, n_ip1);
    cmatrix_scale(interaction, c_real(U));

    tmp = cmatrix_add(H_jw, interaction);

    cmatrix_free(n_i);
    cmatrix_free(n_ip1);
    cmatrix_free(interaction);
    cmatrix_free(H_jw);

    H_jw = tmp;
  }

  for (int j = 0; j < n_modes; j++) {
    cmatrix_free(a[j]);
    cmatrix_free(a_dag[j]);
  }

  // Build H via the independent direct construction
  cmatrix_t *H_direct =
      second_quant_build_hopping_hamiltonian(n_modes, epsilon, t, U);

  double max_diff = 0.0;
  for (int k = 0; k < dim * dim; k++) {
    double d = c_abs2(c_sub(H_jw->data[k], H_direct->data[k]));

    if (d > max_diff) {
      max_diff = d;
    }
  }

  check_close(max_diff, 0.0, 1e-20,
              "H_JW and H_direct agree exactly (independent constructions)");

  // Eigenvalues should also match exactly
  cmatrix_t *H_jw_copy = cmatrix_copy(H_jw);
  cmatrix_t *H_direct_copy = cmatrix_copy(H_direct);
  eigen_t *eig_jw = cmatrix_eigh_complex(H_jw_copy);
  eigen_t *eig_direct = cmatrix_eigh_complex(H_direct_copy);
  cmatrix_free(H_jw_copy);
  cmatrix_free(H_direct_copy);

  double max_eig_diff = 0.0;
  for (int k = 0; k < dim; k++) {
    double d = fabs(eig_jw->eigenvalues[k] - eig_direct->eigenvalues[k]);
    if (d > max_eig_diff) {
      max_eig_diff = d;
    }
  }

  check_close(max_eig_diff, 0.0, 1e-9,
              "eigenvalue spectra of H_JW and H_direct agree");

  printf("  ground state energy: %.6f\n", eig_direct->eigenvalues[0]);

  eigen_free(eig_jw);
  eigen_free(eig_direct);
  cmatrix_free(H_jw);
  cmatrix_free(H_direct);
}

static void test_invalid_input(void) {
  printf("test_invalid_input:\n");

  check_true(jw_creation_operator(-1, 3) == NULL, "negative mode rejected");
  check_true(jw_creation_operator(3, 3) == NULL, "mode >= n_modes rejected");
  check_true(jw_annihilation_operator(0, 0) == NULL, "n_modes=0 rejected");

  const double eps[2] = {0.0, 0.0};
  check_true(second_quant_build_hopping_hamiltonian(0, eps, 1.0, 0.0) == NULL,
             "n_modes=0 rejected (Hamiltonian builder)");
  check_true(second_quant_build_hopping_hamiltonian(2, NULL, 1.0, 0.0) == NULL,
             "NULL epsilon rejected");
}

int main(void) {
  test_anticommutation();
  test_number_operator();
  test_jw_vs_direct_hamiltonian();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_second_quant checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
