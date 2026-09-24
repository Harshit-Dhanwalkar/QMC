/*
 * CASSCF: Simultaneous Orbital + CI Optimization
 *
 * NOTE: eg_68_cisd.c's CISD diagonalizes a fixed-orbital Hamiltonian
 * (RHF orbitals never change). CASSCF goes further: it simultaneously
 * re-optimizes the molecular orbitals themselves alongside a full CI
 * within a chosen "active space" of orbitals, recovering the multi-
 * reference (static) correlation that single-reference methods like RHF,
 * MP2, and CCSD can struggle with. This example shows the orbitals
 * actually rotating: each iteration's active-space CI energy improves as
 * the orbital-optimization loop proceeds.
 */

#include "../physics/casscf.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <stdio.h>

int main(void) {
  printf(" > CASSCF: Simultaneous Orbital + CI Optimization\n\n");

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

  molecular_hf_result_t *hf = molecular_rhf(basis, 6, mol, 4, 1e-12, 200);

  printf("  LiH/STO-3G, R=%.3f bohr\n", R);
  printf("  RHF energy:            %.8f Hartree\n\n", hf->total_energy);

  printf("  CASSCF(2e,2o), freezing the Li 1s core orbital:\n");
  casscf_result_t *r1 = casscf_run(basis, 6, mol, hf->C, 1, 2, 2, 1e-6, 60);
  printf("    Converged in %d iterations (final |grad|=%.2e)\n",
         r1->iterations, r1->grad_norm);
  printf("    CASSCF(2,2) energy:  %.8f Hartree\n",
         r1->total_energy);
  printf("    Correlation recovered vs RHF: %.6f Hartree\n\n",
         hf->total_energy - r1->total_energy);

  printf("  CASSCF(4e,3o), all 4 electrons active, no frozen core:\n");
  casscf_result_t *r2 = casscf_run(basis, 6, mol, hf->C, 0, 3, 4, 1e-6, 60);
  printf("    Converged in %d iterations (final |grad|=%.2e)\n",
         r2->iterations, r2->grad_norm);
  printf("    CASSCF(4,3) energy:  %.8f Hartree\n",
         r2->total_energy);
  printf("    Correlation recovered vs RHF: %.6f Hartree\n\n",
         hf->total_energy - r2->total_energy);

  casscf_result_free(r1);
  casscf_result_free(r2);
  molecular_hf_result_free(hf);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  basis_function_free(h_orb);
  molecule_free(mol);

  return 0;
}
