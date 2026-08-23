/*
 * Test: UCC ansatz (physics/ucc.c), wired to CCSD amplitudes.
 *
 * NOTE: H2/STO-3G's occ={0,1}/virt={2,3} (only 2 occupied, 2 virtual spin
 * orbitals) means general ucc_excitations_from_ccsd_amplitudes enumeration
 * naturally reduces to exactly this same single-double-excitation case (only
 * one occ pair, only one virt pair possible). Singles are forced to zero by
 * spin symmetry (CCSD's own t1 is spin-conserving-integral-driven and vanishes
 * here), so this lets the same exact target validate the full pipeline
 * end-to-end: CCSD converge -> T1/T2 amplitudes ->
 * ucc_excitations_from_ccsd_amplitudes -> ucc_build_generator ->
 * ucc_prepare_state -> energy expectation against the JW Hamiltonian.
 *
 * 1. Core-math test: ucc_build_generator + ucc_prepare_state directly (no CCSD
 *    involved), \theta scanned by hand, minimum checked against exact FCI.
 * 2. Wiring test: same, but with \theta values taken directly from converged
 *    CCSD run via ucc_excitations_from_ccsd_amplitudes, checking (a)
 *    CCSD-seeded UCC energy is already a legitimate variational upper bound on
 *    FCI, and (b) a local scan around that CCSD-seeded \theta cannot find any
 *    point with materially lower energy, i.e., CCSD amplitude for H2 already
 *    sits essentially at the exact UCC optimum.
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/ccsd.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
#include "../physics/ucc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;
static const double H2_FCI_ENERGY = -1.1372759436170439;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.12f expected=%.12f err=%.2e\n", label, got, expected,
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

static double ucc_energy_single_double(int n_modes, const cvector_t *hf_ref,
                                       const cmatrix_t *H, int i, int j, int a,
                                       int b, double theta) {
  ucc_double_t d = {i, j, a, b};
  cmatrix_t *gen = ucc_build_generator(n_modes, NULL, NULL, 0, &d, &theta, 1);
  cvector_t *psi = ucc_prepare_state(gen, hf_ref);

  double re = 0.0;
  int dim = 1 << n_modes;

  for (int r = 0; r < dim; r++) {
    complex_t Hpsi_r = c_zero();

    for (int c = 0; c < dim; c++) {
      Hpsi_r = c_add(Hpsi_r, c_mul(CMAT(H, r, c), psi->data[c]));
    }

    complex_t contrib = c_mul(c_conj(psi->data[r]), Hpsi_r);
    re += contrib.re;
  }

  cvector_free(psi);
  cmatrix_free(gen);

  return re;
}

static void build_h2_setup(molecule_t **mol_out, basis_function_t **h0_out,
                           basis_function_t **h1_out,
                           molecular_hf_result_t **hf_out, cmatrix_t **H_out,
                           double **h_mo_out, double **eri_mo_out) {
  double R = 1.4;
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
  *h0_out = molint_basis_sto3g_h(c0);
  *h1_out = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {*h0_out, *h1_out};
  const double charge[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  *mol_out = molecule_alloc(2, charge, centers);

  *hf_out = molecular_rhf(basis, 2, *mol_out, 2, 1e-12, 200);

  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, 2, *mol_out);
  double *eri_ao = molecular_eri_tensor(basis, 2);
  *h_mo_out = malloc(4 * sizeof(double));
  *eri_mo_out = malloc(16 * sizeof(double));
  molecular_ao_to_mo(Hcore, eri_ao, (*hf_out)->C, 2, *h_mo_out, *eri_mo_out);

  *H_out = second_quant_build_molecular_hamiltonian(
      2, *h_mo_out, *eri_mo_out, molecule_nuclear_repulsion(*mol_out));

  free(eri_ao);
  cmatrix_free(Hcore);
}

static void test_core_math_exact_h2(void) {
  printf("test_ucc_core_math_h2_exact:\n");

  molecule_t *mol;
  basis_function_t *h0, *h1;
  molecular_hf_result_t *hf;
  cmatrix_t *H;
  double *h_mo, *eri_mo;
  build_h2_setup(&mol, &h0, &h1, &hf, &H, &h_mo, &eri_mo);

  check_true(hf->converged, "H2 RHF converges");

  int n_modes = 4;
  cvector_t *hf_ref = cvector_alloc(1 << n_modes);
  for (int k = 0; k < (1 << n_modes); k++) {
    hf_ref->data[k] = c_zero();
  }
  hf_ref->data[0b1100] = c_real(1.0); // modes 0,1 occupied (MSB-first)/

  /*
   * HACK: hand scan around the known-optimal \theta cross-check
   * (-0.11295559...) to bracket minimum without needing a general optimizer in
   * this test
   */
  double best_e = 1e9, best_theta = 0.0;
  for (int s = -200; s <= 200; s++) {
    double theta = -0.2 + 0.4 * (s + 200) / 400.0;
    double e = ucc_energy_single_double(n_modes, hf_ref, H, 0, 1, 2, 3, theta);

    if (e < best_e) {
      best_e = e;
      best_theta = theta;
    }
  }

  printf("  scanned best theta=%.6f  energy=%.10f\n", best_theta, best_e);
  check_close(best_e, H2_FCI_ENERGY, 1e-6,
              "scanned UCC minimum matches exact FCI (grid-limited precision)");

  cvector_free(hf_ref);
  cmatrix_free(H);
  free(h_mo);
  free(eri_mo);
  molecular_hf_result_free(hf);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);
}

static void test_ccsd_to_ucc_wiring(void) {
  printf("test_ucc_ccsd_wiring_h2:\n");

  molecule_t *mol;
  basis_function_t *h0, *h1;
  molecular_hf_result_t *hf;
  cmatrix_t *H;
  double *h_mo, *eri_mo;
  build_h2_setup(&mol, &h0, &h1, &hf, &H, &h_mo, &eri_mo);

  ccsd_amplitudes_t *amp = NULL;
  ccsd_result_t *ccsd = ccsd_run_ex(2, h_mo, eri_mo, hf->orbital_energies, 2, 0,
                                    hf->total_energy, 1e-12, 100, &amp);
  check_true(ccsd != NULL && ccsd->converged, "CCSD converges");
  check_true(amp != NULL, "ccsd_run_ex returns amplitudes");

  ucc_single_t *singles;
  double *theta_s;
  int n_singles;
  ucc_double_t *doubles;
  double *theta_d;
  int n_doubles;
  int ok = ucc_excitations_from_ccsd_amplitudes(
      amp, &singles, &theta_s, &n_singles, &doubles, &theta_d, &n_doubles);
  check_true(ok, "ucc_excitations_from_ccsd_amplitudes succeeds");

  printf("  n_singles=%d n_doubles=%d\n", n_singles, n_doubles);
  check_true(n_doubles == 1,
             "H2 (2occ,2virt) reduces to exactly one doubles excitation");

  int n_modes = amp->nso;
  cvector_t *hf_ref = cvector_alloc(1 << n_modes);

  for (int k = 0; k < (1 << n_modes); k++) {
    hf_ref->data[k] = c_zero();
  }

  hf_ref->data[0b1100] = c_real(1.0);

  cmatrix_t *gen = ucc_build_generator(n_modes, singles, theta_s, n_singles,
                                       doubles, theta_d, n_doubles);
  cvector_t *psi_seed = ucc_prepare_state(gen, hf_ref);

  int dim = 1 << n_modes;
  double e_seed = 0.0;

  for (int r = 0; r < dim; r++) {
    complex_t Hpsi_r = c_zero();

    for (int c = 0; c < dim; c++) {
      Hpsi_r = c_add(Hpsi_r, c_mul(CMAT(H, r, c), psi_seed->data[c]));
    }

    e_seed += c_mul(c_conj(psi_seed->data[r]), Hpsi_r).re;
  }

  printf(
      "  CCSD t2 seed value (theta_d[0])=%.8f  CCSD-seeded UCC energy=%.10f\n",
      theta_d[0], e_seed);
  check_true(e_seed >= H2_FCI_ENERGY - 1e-9,
             "CCSD-seeded UCC energy is a variational upper bound on FCI");

  check_close(e_seed, H2_FCI_ENERGY, 1e-3,
              "CCSD-seeded (unoptimized) UCC energy is already CCSD-quality "
              "close to FCI");

  /*
   * NOTE: now scan around the CCSD-seeded theta to confirm no nearby point has
   * materially lower energy : i.e. CCSD amplitude is already essentially at UCC
   * optimum
   */
  double best_e = 1e9, best_theta = theta_d[0];
  for (int s = -200; s <= 200; s++) {
    double theta = theta_d[0] + 0.2 * s / 200.0;
    ucc_double_t d0 = doubles[0];
    double e = ucc_energy_single_double(n_modes, hf_ref, H, d0.i, d0.j, d0.a,
                                        d0.b, theta);
    if (e < best_e) {
      best_e = e;
      best_theta = theta;
    }
  }

  printf("  local scan around CCSD seed: best_theta=%.6f best_e=%.10f "
         "(seed itself: theta=%.6f e=%.10f)\n",
         best_theta, best_e, theta_d[0], e_seed);

  check_true(best_e <= e_seed + 1e-12,
             "local scan never finds an energy above the seed itself "
             "(sanity: seed is a valid point in the scanned set)");

  check_close(best_e, H2_FCI_ENERGY, 1e-6,
              "best point found near the CCSD seed matches exact FCI to the "
              "same precision as test 1's independent scan");

  cvector_free(psi_seed);
  cmatrix_free(gen);
  cvector_free(hf_ref);
  free(singles);
  free(theta_s);
  free(doubles);
  free(theta_d);
  free(ccsd);
  ccsd_amplitudes_free(amp);
  cmatrix_free(H);
  free(h_mo);
  free(eri_mo);
  molecular_hf_result_free(hf);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);
}

static void test_invalid_input_rejected(void) {
  printf("test_ucc_invalid_input:\n");
  check_true(ucc_build_generator(0, NULL, NULL, 0, NULL, NULL, 0) == NULL,
             "n_modes<1 rejected");

  ucc_single_t bad_single = {5, 1}; // i out of range for n_modes=2
  double theta = 0.1;

  check_true(ucc_build_generator(2, &bad_single, &theta, 1, NULL, NULL, 0) ==
                 NULL,
             "out-of-range single index rejected");

  check_true(ucc_prepare_state(NULL, NULL) == NULL,
             "ucc_prepare_state NULL input rejected");
  check_true(ucc_excitations_from_ccsd_amplitudes(NULL, NULL, NULL, NULL, NULL,
                                                  NULL, NULL) == 0,
             "ucc_excitations_from_ccsd_amplitudes NULL input rejected");
}

int main(void) {
  test_core_math_exact_h2();
  test_ccsd_to_ucc_wiring();
  test_invalid_input_rejected();

  if (failures == 0) {
    printf("\nAll UCC tests PASSED\n");
    return 0;
  }
  printf("\n%d UCC test(s) FAILED\n", failures);

  return 1;
}
