/*
 * Test: general Gaussian94-format basis-set-file parser (basis_parser.c).
 *
 * Validation strategy: rather than trust the parser in isolation, cross-check
 * it two ways against pre-existing, already-validated hardcoded STO-3G builders
 * (molint_basis_sto3g_h / molint_basis_sto3g_li):
 *   1. Parse a literal STO-3G basis-set-file string (published H/Li parameters
 *      as hardcoded builders, Szabo & Ostlund Table 3.7 / EMSL Basis Set
 *      Exchange) and diff every exponent/coefficient against hardcoded arrays
 *      exactly (they are literature numbers, so this must match to full double
 *      precision).
 *   2. Run the parser's output through the exact downstream pipeline
 *      (molecular_rhf) used by test_molecular_hf.c / test_lih.c and confirm SCF
 *      energies match those tests' already-validated reference values (H2
 *      -1.116714 Hartree at R=1.4 bohr; LiH's RHF energy), which only happens
 *       if every primitive, contraction coefficient, and Cartesian-component
 *       expansion (including the SP-shell split) is correct end to end.
 *
 * Covers: malformed input rejection, D-exponent Fortran notation, comment
 * lines, element lookup, and general (non-minimal) multi-shell-type parsing via
 * a synthetic basis with an explicit D shell.
 */

#include "../physics/basis_parser.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

/*
 * Published STO-3G H + Li text as EMSL/Basis Set Exchange Gaussian94
 * export
 */
static const char *STO3G_TEXT =
    "!  STO-3G  EMSL Basis Set Exchange Library\n"
    "\n"
    "H     0\n"
    "S   3   1.00\n"
    "      3.42525091              0.15432897\n"
    "      0.62391373              0.53532814\n"
    "      0.16885540              0.44463454\n"
    "****\n"
    "LI     0\n"
    "S   3   1.00\n"
    "     16.1195750               0.15432897\n"
    "      2.9362007               0.53532814\n"
    "      0.7946505               0.44463454\n"
    "SP   3   1.00\n"
    "      0.6362897              -0.09996723             0.15591627\n"
    "      0.1478601               0.39951283             0.60768372\n"
    "      0.0480887               0.70011547             0.39195739\n"
    "****\n";

static void test_parse_basic_structure(void) {
  printf("Test: basic parse structure (element/shell counts)\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_TEXT);
  check(bs != NULL, "parse should succeed on well-formed text");
  if (!bs) {
    return;
  }

  check(bs->n_elements == 2, "should find exactly 2 element blocks (H, LI)");

  const basis_element_t *h = basis_set_find_element(bs, "H");
  check(h != NULL, "H element should be found");
  if (h) {
    check(h->n_shells == 1, "H should have exactly 1 shell (S)");
    check(strcmp(h->element, "H") == 0, "H element symbol stored correctly");

    if (h->n_shells == 1) {
      check(h->shells[0].n_prim == 3, "H S shell should have 3 primitives");
      check(strcmp(h->shells[0].type, "S") == 0, "H shell type is S");
    }
  }

  const basis_element_t *li = basis_set_find_element(bs, "li"); /* lowercase */
  check(li != NULL, "Li element should be found case-insensitively");
  if (li) {
    check(li->n_shells == 2, "Li should have exactly 2 shells (S, SP)");

    if (li->n_shells == 2) {
      check(strcmp(li->shells[0].type, "S") == 0, "Li shell 0 is S");
      check(strcmp(li->shells[1].type, "SP") == 0, "Li shell 1 is SP");
      check(li->shells[1].coef2 != NULL, "SP shell must have a coef2 (P part)");
    }
  }

  const basis_element_t *missing = basis_set_find_element(bs, "C");
  check(missing == NULL, "unlisted element should not be found");

  basis_set_free(bs);
}

static void test_parsed_values_match_hardcoded_h(void) {
  printf("Test: parsed H exponents/coefficients exactly match hardcoded "
         "molint_basis_sto3g_h\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_TEXT);
  check(bs != NULL, "parse should succeed");
  if (!bs) {
    return;
  }

  const basis_element_t *h = basis_set_find_element(bs, "H");
  check(h != NULL && h->n_shells == 1, "H parsed correctly");
  if (h && h->n_shells == 1) {
    static const double ref_alphas[3] = {3.42525091, 0.62391373, 0.16885540};
    static const double ref_coeffs[3] = {0.15432897, 0.53532814, 0.44463454};

    for (int i = 0; i < 3; i++) {
      check_close(h->shells[0].exps[i], ref_alphas[i], 1e-12,
                  "H exponent matches hardcoded literature value exactly");
      check_close(h->shells[0].coef1[i], ref_coeffs[i], 1e-12,
                  "H raw coefficient matches hardcoded literature value "
                  "exactly");
    }
  }

  basis_set_free(bs);
}

static void test_build_atom_h_matches_hardcoded(void) {
  printf("Test: basis_set_build_atom(H) produces a normalized basis "
         "function matching molint_basis_sto3g_h's overlap-normalized "
         "output\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_TEXT);
  check(bs != NULL, "parse should succeed");
  if (!bs) {
    return;
  }

  const double center[3] = {0.0, 0.0, 0.7};
  const basis_element_t *h = basis_set_find_element(bs, "H");
  check(h != NULL, "H found");

  basis_function_t **parsed_fns = NULL;
  int n = basis_set_build_atom(h, center, &parsed_fns);
  check(n == 1, "H should expand to exactly 1 basis function (1s, S shell)");

  basis_function_t *ref_fn = molint_basis_sto3g_h(center);
  check(ref_fn != NULL, "reference builder should succeed");

  if (n == 1 && parsed_fns && ref_fn) {
    /* Self-overlap must be 1.0 (normalization), and parsed function's overlap
     * with reference hardcoded function must also be exactly 1.0 : identical
     * functions have unit overlap with themselves. */
    double self_overlap = gto_overlap(parsed_fns[0], parsed_fns[0]);
    check_close(self_overlap, 1.0, 1e-10,
                "parsed H function is properly normalized (self-overlap=1)");

    double cross_overlap = gto_overlap(parsed_fns[0], ref_fn);
    check_close(cross_overlap, 1.0, 1e-10,
                "parsed H function has unit overlap with hardcoded "
                "reference (i.e. they are the same function)");
  }

  if (ref_fn) {
    basis_function_free(ref_fn);
  }

  basis_set_free_functions(parsed_fns, n);
  basis_set_free(bs);
}

static void test_build_atom_li_matches_hardcoded(void) {
  printf("Test: basis_set_build_atom(Li) SP-shell split produces 5 "
         "functions (1s,2s,2px,2py,2pz) matching molint_basis_sto3g_li\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_TEXT);
  check(bs != NULL, "parse should succeed");
  if (!bs) {
    return;
  }

  const double center[3] = {0.0, 0.0, 0.0};
  const basis_element_t *li = basis_set_find_element(bs, "LI");
  check(li != NULL, "Li found");

  basis_function_t **parsed_fns = NULL;
  int n = basis_set_build_atom(li, center, &parsed_fns);
  check(n == 5,
        "Li should expand to exactly 5 functions (1s + SP->2s,2px,2py,2pz)");

  basis_function_t *ref_fns[5];
  int ok = molint_basis_sto3g_li(center, ref_fns);
  check(ok == 1, "reference builder should succeed");

  if (n == 5 && parsed_fns && ok) {
    const char *names[5] = {"1s", "2s", "2px", "2py", "2pz"};

    for (int i = 0; i < 5; i++) {
      double cross = gto_overlap(parsed_fns[i], ref_fns[i]);
      char msg[128];

      snprintf(msg, sizeof(msg),
               "parsed Li %s matches hardcoded reference (unit overlap)",
               names[i]);

      check_close(cross, 1.0, 1e-9, msg);
    }
  }

  if (ok) {
    for (int i = 0; i < 5; i++) {
      basis_function_free(ref_fns[i]);
    }
  }

  basis_set_free_functions(parsed_fns, n);
  basis_set_free(bs);
}

static void test_build_molecule_h2_rhf_energy(void) {
  printf("Test: end-to-end H2/STO-3G RHF energy from parsed basis-set-file "
         "text matches test_molecular_hf.c's reference (-1.116714 Hartree @ "
         "R=1.4 bohr)\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_TEXT);
  check(bs != NULL, "parse should succeed");
  if (!bs) {
    return;
  }

  double R = 1.4;
  double centerA[3] = {0.0, 0.0, 0.0};
  double centerB[3] = {0.0, 0.0, R};
  const char *symbols[2] = {"H", "H"};
  double centers[2][3];

  memcpy(centers[0], centerA, sizeof(centerA));
  memcpy(centers[1], centerB, sizeof(centerB));

  basis_function_t **basis = NULL;
  int n_basis = basis_set_build_molecule(bs, symbols, centers, 2, &basis);
  check(n_basis == 2, "H2/STO-3G should have exactly 2 basis functions");

  const double charges[2] = {1.0, 1.0};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  check(mol != NULL, "molecule_alloc should succeed");

  if (n_basis == 2 && basis && mol) {
    molecular_hf_result_t *scf =
        molecular_rhf(basis, n_basis, mol, 2, 1e-10, 200);

    check(scf != NULL, "RHF should converge");
    if (scf) {
      check(scf->converged, "RHF should report converged");

      check_close(scf->total_energy, -1.116714, 1e-5,
                  "H2/STO-3G total RHF energy from parsed basis matches known "
                  "reference");
      molecular_hf_result_free(scf);
    }
  }

  if (mol) {
    molecule_free(mol);
  }

  basis_set_free_functions(basis, n_basis);
  basis_set_free(bs);
}

static void test_build_molecule_lih_rhf_energy(void) {
  printf("Test: end-to-end LiH/STO-3G RHF energy from parsed basis-set-file "
         "text matches test_lih.c's reference (-7.862009272 Hartree @ R=3.015 "
         "bohr)\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_TEXT);
  check(bs != NULL, "parse should succeed");
  if (!bs) {
    return;
  }

  double R = 3.015;
  const char *symbols[2] = {"LI", "H"};
  double centers[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, R}};

  basis_function_t **basis = NULL;
  int n_basis = basis_set_build_molecule(bs, symbols, centers, 2, &basis);

  check(n_basis == 6, "LiH/STO-3G should have exactly 6 basis functions "
                      "(Li: 1s,2s,2px,2py,2pz + H: 1s)");

  const double charges[2] = {3.0, 1.0};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  check(mol != NULL, "molecule_alloc should succeed");

  if (n_basis == 6 && basis && mol) {
    molecular_hf_result_t *scf =
        molecular_rhf(basis, n_basis, mol, 4, 1e-10, 200);

    check(scf != NULL, "RHF should converge");
    if (scf) {
      check(scf->converged, "RHF should report converged");

      check_close(scf->total_energy, -7.862009272, 1e-5,
                  "LiH/STO-3G total RHF energy from parsed basis (incl. the"
                  "SP-shell split path) matches known reference");

      molecular_hf_result_free(scf);
    }
  }

  if (mol) {
    molecule_free(mol);
  }

  basis_set_free_functions(basis, n_basis);
  basis_set_free(bs);
}

static void test_malformed_input_rejected(void) {
  printf("Test: malformed basis-set text is rejected (returns NULL), not "
         "silently mis-parsed\n");

  // Shell header claims 3 primitives but only provides 1 line before EOF
  const char *truncated = "H     0\nS   3   1.00\n   3.42525091   0.1543\n";
  basis_set_t *bs1 = basis_set_parse_string(truncated);

  check(bs1 == NULL, "truncated shell body should be rejected");

  basis_set_free(bs1);

  // Unrecognized shell-type letter
  const char *bad_type = "H     0\nQ   1   1.00\n   1.0   1.0\n****\n";
  basis_set_t *bs2 = basis_set_parse_string(bad_type);

  check(bs2 == NULL, "unrecognized shell type letter should be rejected");

  basis_set_free(bs2);

  check(basis_set_parse_string(NULL) == NULL, "NULL input should return NULL");
}

static void test_fortran_d_exponent_notation(void) {
  printf("Test: Fortran-style 'D' exponent notation parses identically to "
         "'E'/plain notation\n");

  const char *plain_text = "H     0\nS   1   1.00\n   1.5   1.0\n****\n";
  const char *fortran_text =
      "H     0\nS   1   1.00\n   1.5D+00   1.0D+00\n****\n";

  basis_set_t *bs_plain = basis_set_parse_string(plain_text);
  basis_set_t *bs_fortran = basis_set_parse_string(fortran_text);
  check(bs_plain != NULL && bs_fortran != NULL,
        "both variants should parse successfully");

  if (bs_plain && bs_fortran) {
    const basis_element_t *h_plain = basis_set_find_element(bs_plain, "H");
    const basis_element_t *h_fortran = basis_set_find_element(bs_fortran, "H");

    check(h_plain && h_fortran && h_plain->n_shells == 1 &&
              h_fortran->n_shells == 1,
          "both should have 1 shell");

    if (h_plain && h_fortran) {
      check_close(h_plain->shells[0].exps[0], h_fortran->shells[0].exps[0],
                  1e-14, "D-notation exponent equals plain-notation exponent");
      check_close(h_plain->shells[0].coef1[0], h_fortran->shells[0].coef1[0],
                  1e-14,
                  "D-notation coefficient equals plain-notation coefficient");
    }
  }

  basis_set_free(bs_plain);
  basis_set_free(bs_fortran);
}

static void test_general_d_shell_expansion(void) {
  printf("Test: a synthetic D shell expands into 6 Cartesian components, "
         "each individually normalized\n");

  /* Not a real published basis set : just a single D-shell carbon-like
   * exponent, used purely to exercise l=2 Cartesian expansion beyond the STO-3G
   * s/p-only cases above. */
  const char *text = "C     0\nD   1   1.00\n   0.8   1.0\n****\n";

  basis_set_t *bs = basis_set_parse_string(text);
  check(bs != NULL, "D-shell text should parse");
  if (!bs) {
    return;
  }

  const basis_element_t *c = basis_set_find_element(bs, "C");
  check(c != NULL && c->n_shells == 1, "C element with 1 D shell parsed");

  if (c && c->n_shells == 1) {
    check(c->shells[0].l == 2, "D shell angular momentum is 2");

    const double center[3] = {0.0, 0.0, 0.0};
    basis_function_t **fns = NULL;

    int n = basis_set_build_atom(c, center, &fns);
    check(n == 6, "D shell should expand to 6 Cartesian components "
                  "(xx,xy,xz,yy,yz,zz)");

    if (n == 6 && fns) {
      int total_l_sum = 0;

      for (int i = 0; i < 6; i++) {
        double self_overlap = gto_overlap(fns[i], fns[i]);

        check_close(self_overlap, 1.0, 1e-9,
                    "each D component is individually normalized");

        total_l_sum += fns[i]->l + fns[i]->m + fns[i]->n;

      }

      check(total_l_sum == 12, "all 6 components carry total angular "
                               "momentum l=2 each (sum=12)");
    }

    basis_set_free_functions(fns, n);
  }

  basis_set_free(bs);
}

static void test_build_molecule_missing_element_fails(void) {
  printf("Test: basis_set_build_molecule fails cleanly when an atom's "
         "element is absent from the basis set\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_TEXT); // only has H, Li
  check(bs != NULL, "parse should succeed");
  if (!bs) {
    return;
  }

  const char *symbols[1] = {"O"}; /* not present */
  double centers[1][3] = {{0.0, 0.0, 0.0}};
  basis_function_t **basis = NULL;

  int n = basis_set_build_molecule(bs, symbols, centers, 1, &basis);
  check(n == 0, "should return 0 for a molecule with a missing element");
  check(basis == NULL, "output array should be NULL on failure");

  basis_set_free(bs);
}

int main(void) {
  printf("=== Basis-set-file parser tests ===\n\n");

  test_parse_basic_structure();
  test_parsed_values_match_hardcoded_h();
  test_build_atom_h_matches_hardcoded();
  test_build_atom_li_matches_hardcoded();
  test_build_molecule_h2_rhf_energy();
  test_build_molecule_lih_rhf_energy();
  test_malformed_input_rejected();
  test_fortran_d_exponent_notation();
  test_general_d_shell_expansion();
  test_build_molecule_missing_element_fails();

  printf("\n=== %s ===\n", failures == 0 ? "ALL TESTS PASSED" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
