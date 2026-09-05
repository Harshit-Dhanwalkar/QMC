/*
 * Full Configuration Interaction (FCI): Comparing the Quantum-Chemistry Method
 * Hierarchy on H2 and LiH
 *
 * NOTE: physics/fci.c formalizes exact diagonalization within a fixed-electron-
 * number Fock-space sector
 */

#include "../core/matrix.h"
#include "../physics/ccsd.h"
#include "../physics/ccsd_t.h"
#include "../physics/fci.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <stdio.h>
#include <stdlib.h>

static void run_method_ladder(const char *label, basis_function_t **basis,
                              int n_basis, const molecule_t *mol,
                              int n_electrons) {
  printf(" > %s\n\n", label);

  molecular_hf_result_t *hf =
      molecular_rhf(basis, n_basis, mol, n_electrons, 1e-12, 200);
  printf("  RHF:      %.8f Hartree\n", hf->total_energy);

  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, n_basis, mol);
  double *eri_ao = molecular_eri_tensor(basis, n_basis);
  double *h_mo = malloc((size_t)n_basis * n_basis * sizeof(double));
  double *eri_mo =
      malloc((size_t)n_basis * n_basis * n_basis * n_basis * sizeof(double));
  molecular_ao_to_mo(h_ao, eri_ao, hf->C, n_basis, h_mo, eri_mo);
  double Enuc = molecule_nuclear_repulsion(mol);

  ccsd_result_t *ccsd = ccsd_run(n_basis, h_mo, eri_mo, hf->orbital_energies,
                                 n_electrons, 0, hf->total_energy, 1e-10, 100);
  printf("  CCSD:     %.8f Hartree\n", ccsd->total_energy);

  ccsdt_result_t *ccsdt =
      ccsdt_run(n_basis, h_mo, eri_mo, hf->orbital_energies, n_electrons, 0,
                hf->total_energy, 1e-10, 100);
  printf("  CCSD(T):  %.8f Hartree\n", ccsdt->total_energy);

  fci_result_t *fci = fci_solve(n_basis, h_mo, eri_mo, Enuc, n_electrons);
  printf("  FCI:      %.8f Hartree  (exact, for this basis set)\n\n",
         fci->ground_energy);

  double rhf_corr = hf->total_energy - fci->ground_energy;
  double ccsd_corr = ccsd->total_energy - fci->ground_energy;
  double ccsdt_corr = ccsdt->total_energy - fci->ground_energy;
  printf("  Error vs. exact FCI:\n");
  printf("    RHF:     %.6f Hartree (%.2f%% of correlation energy missed)\n",
         rhf_corr, 100.0);
  printf("    CCSD:    %.6f Hartree (%.4f%% of RHF's error remains)\n",
         ccsd_corr, 100.0 * ccsd_corr / rhf_corr);
  printf("    CCSD(T): %.6f Hartree (%.4f%% of RHF's error remains)\n\n",
         ccsdt_corr, 100.0 * ccsdt_corr / rhf_corr);

  free(ccsd);
  free(ccsdt);
  fci_result_free(fci);
  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(h_ao);
  molecular_hf_result_free(hf);
}

int main(void) {
  printf(" > FCI and the Quantum-Chemistry Method Hierarchy\n\n");

  {
    double R = 1.4;
    double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
    basis_function_t *h0 = molint_basis_sto3g_h(c0);
    basis_function_t *h1 = molint_basis_sto3g_h(c1);
    basis_function_t *basis[2] = {h0, h1};
    const double charges[2] = {1.0, 1.0};
    double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
    molecule_t *mol = molecule_alloc(2, charges, centers);

    run_method_ladder("H2/STO-3G, R=1.4 bohr", basis, 2, mol, 2);

    basis_function_free(h0);
    basis_function_free(h1);
    molecule_free(mol);
  }

  {
    double R = 3.015;
    double c_li[3] = {0, 0, 0}, c_h[3] = {0, 0, R};
    basis_function_t *li_orbs[5];
    molint_basis_sto3g_li(c_li, li_orbs);
    basis_function_t *h_orb = molint_basis_sto3g_h(c_h);
    basis_function_t *basis[6] = {li_orbs[0], li_orbs[1], li_orbs[2],
                                  li_orbs[3], li_orbs[4], h_orb};
    const double charges[2] = {3.0, 1.0};
    double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
    molecule_t *mol = molecule_alloc(2, charges, centers);

    run_method_ladder("LiH/STO-3G, R=3.015 bohr", basis, 6, mol, 4);

    for (int i = 0; i < 5; i++) {
      basis_function_free(li_orbs[i]);
    }
    basis_function_free(h_orb);
    molecule_free(mol);
  }

  printf("For H2 (a 2-electron system), CCSD is exactly equivalent to FCI "
         "error vs. FCI above is essentially zero. LiH shows CCSD(T)'s usual "
         "pattern: recovering vast majority of RHF's missing correlation "
         "energy that plain CCSD leaves behind.\n");

  return 0;
}
