/*
 * General Basis-Set-File Parser: Water (H2O) via a Published Gaussian94
 * Basis-Set Text
 *
 * NOTE: Every prior molecular example (eg_40, eg_42..eg_46) is built from
 * hardcoded convenience builders molint_basis_sto3g_h/li in
 * molecular_integrals.c, which only cover hydrogen and lithium. Water needs
 * oxygen, which those builders don't provide.
 *
 * This example instead parses a real Gaussian94-format basis-set-file
 * (plain-text format Basis Set Exchange, https://www.basissetexchange.org,
 * exports) at runtime via basis_parser.c, covering hydrogen and oxygen's
 * published STO-3G parameters (Reference: Hehre, Stewart & Pople 1969). It
 * demonstrates general path any new element/basis combination would take from
 * here: write/download a Gaussian94-format .gbs text file, parse it once, then
 * build any molecule made of elements it contains.
 */

#include "../physics/basis_parser.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * NOTE: Published STO-3G basis-set text (Gaussian94 format) for H and O. The
 * contraction coefficients (0.15432897/0.53532814/0.44463454 for every S shell;
 * -0.09996723/0.39951283/0.70011547 and 0.15591627/0.60768372/ 0.39195739 for
 * every SP shell's s-/p-parts respectively) are same "universal" STO-3G values
 * used across the periodic table; only per-element exponents change. Same
 * numeric family already validated in tests/test_basis_parser.c against the
 * hardcoded H/Li builders.
 */
static const char *STO3G_H_O_TEXT =
    "!  STO-3G  (Hehre, Stewart & Pople 1969)\n"
    "\n"
    "H     0\n"
    "S   3   1.00\n"
    "      3.42525091              0.15432897\n"
    "      0.62391373              0.53532814\n"
    "      0.16885540              0.44463454\n"
    "****\n"
    "O     0\n"
    "S   3   1.00\n"
    "    130.70932000              0.15432897\n"
    "     23.80886100              0.53532814\n"
    "      6.44360830              0.44463454\n"
    "SP   3   1.00\n"
    "      5.03315130             -0.09996723             0.15591627\n"
    "      1.16959610              0.39951283             0.60768372\n"
    "      0.38038900              0.70011547             0.39195739\n"
    "****\n";

int main(void) {
  printf(" > General Basis-Set-File Parser: Water (H2O), STO-3G\n\n");

  printf("Step 1: parse the STO-3G Gaussian94-format text (H + O)\n\n");

  basis_set_t *bs = basis_set_parse_string(STO3G_H_O_TEXT);
  if (!bs) {
    fprintf(stderr, "Failed to parse basis-set text.\n");
    return 1;
  }

  printf("  Parsed %d element(s):\n", bs->n_elements);
  for (int i = 0; i < bs->n_elements; i++) {
    const basis_element_t *e = &bs->elements[i];

    printf("    %-3s : %d shell(s) [", e->element, e->n_shells);
    for (int s = 0; s < e->n_shells; s++) {
      printf("%s%s", e->shells[s].type, s + 1 < e->n_shells ? "," : "");
    }

    printf("]\n");
  }

  printf("\n");

  /* Near-experimental water geometry: r(O-H) = 0.9584 Angstrom = 1.8111 bohr,
   * HOH angle = 104.45 degrees, O at the origin with the bisector along -z. */
  printf("Step 2: build the water molecule (near-experimental geometry: "
         "r(O-H)=1.8111 bohr, angle=104.45 deg)\n\n");

  double r_oh = 1.8111;
  double half_angle_rad = (104.45 / 2.0) * M_PI / 180.0;
  double centers[3][3] = {
      {0.0, 0.0, 0.0},                                                // O
      {r_oh * sin(half_angle_rad), 0.0, -r_oh * cos(half_angle_rad)}, // H1
      {-r_oh * sin(half_angle_rad), 0.0, -r_oh * cos(half_angle_rad)} // H2
  };

  const char *symbols[3] = {"O", "H", "H"};
  const double charges[3] = {8.0, 1.0, 1.0};

  basis_function_t **basis = NULL;
  int n_basis = basis_set_build_molecule(bs, symbols, centers, 3, &basis);
  if (n_basis == 0) {
    fprintf(stderr, "Failed to build the molecular basis.\n");

    basis_set_free(bs);

    return 1;
  }

  printf("  Total basis functions: %d (O: 1s,2s,2px,2py,2pz = 5, H x2: 1s each "
         "= 2)\n\n",
         n_basis);

  molecule_t *mol = molecule_alloc(3, charges, centers);

  printf("Step 3: general N-basis RHF (10 valence+core electrons: O has 8, "
         "each H has 1)\n\n");

  molecular_hf_result_t *scf =
      molecular_rhf(basis, n_basis, mol, 10, 1e-9, 300);
  if (!scf) {
    fprintf(stderr, "RHF failed to allocate.\n");
  } else {
    printf("  Water/STO-3G RHF total energy: %.6f Hartree (converged=%d, %d "
           "iterations)\n"
           "  (Reference: STO-3G water near this geometry is documented in "
           "theoretically at approximately -74.94 Hartree. The result, exact "
           "digits depend on the precise geometry used)\n\n",
           scf->total_energy, scf->converged, scf->iterations);

    printf("  Occupied orbital energies (Hartree):\n");
    for (int i = 0; i < 5; i++) {
      printf("    MO %d: %.6f\n", i, scf->orbital_energies[i]);
    }

    molecular_hf_result_free(scf);
  }

  molecule_free(mol);
  basis_set_free_functions(basis, n_basis);
  basis_set_free(bs);

  return 0;
}
