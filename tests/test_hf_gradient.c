/*
 * Test: analytic RHF nuclear gradient (Pulay 1969)
 *
 * 1. Gradient integral primitives (gto_*_grad_a, gto_nuclear_attraction_grad_C)
 *    vs finite difference, isolated from the SCF/Pulay-formula machinery:
 *    validates McMurchie-Davidson "increment/decrement" identity (Reference:
 *    Helgaker, Jorgensen & Olsen eq. 9.3.8) gradient code is built on,
 *    independent of any density-matrix or energy-weighted-density.
 * 2. molecular_rhf_gradient vs finite difference of the total SCF energy:
 *    H2/STO-3G at several bond lengths (translational invariance - sum of
 *    forces on both atoms is exactly zero), and LiH/STO-3G (multi-atom,
 *    multiple basis functions sharing a center, including p-orbitals)
 * 3. H2 equilibrium bond length via bisection on the analytic gradient,
 *    cross-checked against literature (STO-3G/RHF H2 equilibrium ~1.34-1.40
 *    bohr) as an end-to-end sanity check that gradient actually locates a real
 *    energy minimum, not just any zero-crossing.
 */

#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
  printf("  %s: %s\n", label, cond ? "ok" : "FAILED");
  if (!cond) {
    failures++;
  }
}

static void test_boys_large_x_regression(void) {
  printf("test_boys_large_x_regression:\n");
  struct {
    int n;
    double x;
    double expected;
  } cases[] = {
      {0, 201.5, 6.243201867e-02},  {0, 290.15, 5.202757904e-02},
      {0, 515.8, 3.902153152e-02},  {0, 805.98, 3.121640117e-02},
      {1, 290.15, 8.965636452e-05}, {1, 515.8, 3.782621777e-05},
  };

  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    double *F = malloc((size_t)(cases[i].n + 1) * sizeof(double));

    boys_function_array(cases[i].n, cases[i].x, F);

    char label[128];
    snprintf(label, sizeof(label), "F_%d(%.2f) matches independent reference",
             cases[i].n, cases[i].x);

    check_close(F[cases[i].n], cases[i].expected, 1e-8, label);

    free(F);
  }
}

static void test_gradient_integral_primitives(void) {
  printf("test_gradient_integral_primitives:\n");
  double centerA[3] = {0.3, -0.2, 0.55};
  double centerB[3] = {-0.4, 0.1, 0.9};
  double centerC[3] = {0.05, 0.6, -0.3};
  const double alphas_a[1] = {0.9}, coeffs_a[1] = {1.0};
  const double alphas_b[1] = {1.3}, coeffs_b[1] = {1.0};

  basis_function_t *a =
      basis_function_alloc(1, 0, 1, centerA, 1, alphas_a, coeffs_a);
  basis_function_t *b =
      basis_function_alloc(0, 2, 0, centerB, 1, alphas_b, coeffs_b);

  double h = 1e-6;
  double Ax_plus[3] = {centerA[0] + h, centerA[1], centerA[2]};
  double Ax_minus[3] = {centerA[0] - h, centerA[1], centerA[2]};

  basis_function_t *ap =
      basis_function_alloc(1, 0, 1, Ax_plus, 1, alphas_a, coeffs_a);
  basis_function_t *am =
      basis_function_alloc(1, 0, 1, Ax_minus, 1, alphas_a, coeffs_a);

  double grad[3];
  gto_overlap_grad_a(a, b, grad);

  double fd = (gto_overlap(ap, b) - gto_overlap(am, b)) / (2 * h);
  check_close(grad[0], fd, 1e-8, "overlap grad_a matches finite difference");

  gto_kinetic_grad_a(a, b, grad);
  fd = (gto_kinetic(ap, b) - gto_kinetic(am, b)) / (2 * h);
  check_close(grad[0], fd, 1e-8, "kinetic grad_a matches finite difference");

  gto_nuclear_attraction_grad_a(a, b, centerC, grad);
  fd = (gto_nuclear_attraction(ap, b, centerC) -
        gto_nuclear_attraction(am, b, centerC)) /
       (2 * h);
  check_close(grad[0], fd, 1e-8,
              "nuclear-attraction grad_a matches finite difference");

  gto_nuclear_attraction_grad_C(a, b, centerC, grad);
  double Cx_plus[3] = {centerC[0] + h, centerC[1], centerC[2]};
  double Cx_minus[3] = {centerC[0] - h, centerC[1], centerC[2]};
  fd = (gto_nuclear_attraction(a, b, Cx_plus) -
        gto_nuclear_attraction(a, b, Cx_minus)) /
       (2 * h);
  check_close(grad[0], fd, 1e-8,
              "nuclear-attraction grad_C (point charge) matches finite "
              "difference");

  basis_function_t *c =
      basis_function_alloc(0, 0, 0, centerC, 1, alphas_a, coeffs_a);
  double centerD[3] = {-0.1, 0.3, -0.2};

  basis_function_t *d =
      basis_function_alloc(0, 0, 0, centerD, 1, alphas_b, coeffs_b);
  gto_eri_grad_a(a, b, c, d, grad);
  fd = (gto_eri(ap, b, c, d) - gto_eri(am, b, c, d)) / (2 * h);
  check_close(grad[0], fd, 1e-8, "ERI grad_a matches finite difference");

  basis_function_free(a);
  basis_function_free(b);
  basis_function_free(c);
  basis_function_free(d);
  basis_function_free(ap);
  basis_function_free(am);
}

static double h2_energy(double R) {
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};

  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};

  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};

  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);

  double E = res->total_energy;

  molecular_hf_result_free(res);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);

  return E;
}

static double h2_gradient_atom1_z(double R) {
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};

  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};

  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};

  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);

  const int atom_of_basis[2] = {0, 1};
  double *grad = molecular_rhf_gradient(basis, 2, mol, atom_of_basis, res);
  double gz = grad[5];

  free(grad);
  molecular_hf_result_free(res);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);

  return gz;
}

static void test_h2_gradient_vs_finite_difference(void) {
  printf("test_h2_gradient_vs_finite_difference:\n");
  double Rs[] = {1.0, 1.4, 2.0, 3.0};
  double h = 1e-4;
  for (size_t i = 0; i < sizeof(Rs) / sizeof(Rs[0]); i++) {
    double R0 = Rs[i];
    double fd = (h2_energy(R0 + h) - h2_energy(R0 - h)) / (2 * h);
    double analytic = h2_gradient_atom1_z(R0);

    char label[64];
    snprintf(label, sizeof(label), "H2 R=%.2f analytic matches FD", R0);

    check_close(analytic, fd, 1e-6, label);
  }

  // translational invariance: sum of forces on both atoms must vanish
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, 1.4};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};

  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, 1.4}};

  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);

  const int atom_of_basis[2] = {0, 1};
  double *grad = molecular_rhf_gradient(basis, 2, mol, atom_of_basis, res);
  check_true(fabs(grad[2] + grad[5]) < 1e-10,
             "sum of forces on both atoms vanishes (translational invariance)");

  free(grad);
  molecular_hf_result_free(res);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);
}

static void lih_energy_and_grad(double R, double *E_out, double grad_out[6]) {
  double cLi[3] = {0, 0, 0}, cH[3] = {0, 0, R};

  basis_function_t *li[5];
  molint_basis_sto3g_li(cLi, li);
  basis_function_t *h = molint_basis_sto3g_h(cH);
  basis_function_t *basis[6] = {li[0], li[1], li[2], li[3], li[4], h};

  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 6, mol, 4, 1e-14, 500);
  if (E_out) {
    *E_out = res->total_energy;
  }
  if (grad_out) {
    const int atom_of_basis[6] = {0, 0, 0, 0, 0, 1};
    double *grad = molecular_rhf_gradient(basis, 6, mol, atom_of_basis, res);

    for (int i = 0; i < 6; i++) {
      grad_out[i] = grad[i];
    }

    free(grad);
  }

  molecular_hf_result_free(res);
  molecule_free(mol);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
  basis_function_free(h);
}

static void test_lih_gradient_vs_finite_difference(void) {
  printf("test_lih_gradient_vs_finite_difference:\n");
  double Rs[] = {2.0, 3.0, 3.5, 4.0};
  double h = 1e-4;

  for (size_t i = 0; i < sizeof(Rs) / sizeof(Rs[0]); i++) {
    double R0 = Rs[i];
    double Ep, Em;

    lih_energy_and_grad(R0 + h, &Ep, NULL);
    lih_energy_and_grad(R0 - h, &Em, NULL);

    double fd = (Ep - Em) / (2 * h);
    double grad[6];

    lih_energy_and_grad(R0, NULL, grad);

    char label[64];
    snprintf(label, sizeof(label), "LiH R=%.2f analytic matches FD", R0);
    check_close(grad[5], fd, 1e-6, label);

    char label2[96];
    snprintf(label2, sizeof(label2),
             "LiH R=%.2f sum of forces vanishes (translational invariance)",
             R0);
    check_true(fabs(grad[2] + grad[5]) < 1e-9, label2);
  }
}

static void test_h2_equilibrium_bond_length(void) {
  printf("test_h2_equilibrium_bond_length:\n");
  double lo = 1.0, hi = 2.0;

  for (int i = 0; i < 40; i++) {
    double mid = 0.5 * (lo + hi);
    double g = h2_gradient_atom1_z(mid);

    if (g > 0) {
      hi = mid;
    } else {
      lo = mid;
    }
  }

  double Req = 0.5 * (lo + hi);
  printf("  H2/STO-3G equilibrium bond length (via gradient root): %.4f "
         "bohr\n",
         Req);
  // Literature STO-3G/RHF H2 equilibrium is ~1.34-1.40 bohr
  check_true(Req > 1.3 && Req < 1.45,
             "equilibrium bond length matches literature STO-3G/RHF range");
}

int main(void) {
  test_boys_large_x_regression();
  test_gradient_integral_primitives();
  test_h2_gradient_vs_finite_difference();
  test_lih_gradient_vs_finite_difference();
  test_h2_equilibrium_bond_length();

  if (failures > 0) {
    printf("\n%d test_hf_gradient check(s) FAILED.\n", failures);
    return 1;
  }

  printf("\nAll test_hf_gradient checks passed.\n");

  return 0;
}
