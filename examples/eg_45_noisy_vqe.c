/*
 * H2 Noisy VQE: T1/T2 Decoherence and Zero-Noise Extrapolation
 *
 * Takes the same H2/STO-3G 4-qubit Hamiltonian and VQE-optimized ansatz
 * parameters as eg_42_h2_vqe.c, but instead of evaluating ansatz as a pure
 * state, evaluates it as a density matrix that picks up T1 (amplitude-damping)
 * and T2 (dephasing) noise after every gate, which drives lindblad.c's
 * already-validated GKSL machinery gate-by-gate rather than reimplementing
 * open-system dynamics.
 *
 * This demonstrates two things a real NISQ device faces that the noiseless
 * VQE examples don't show at all:
 *   1. Decoherence during the circuit measurably degrades the energy estimate,
 *     growing with noise strength; even though the ansatz *parameters* found by
 *     noiseless VQE are unchanged.
 *   2. Zero-noise extrapolation (run the same circuit at amplified noise
 *     levels, then extrapolate back to zero) recovers a meaningfully better
 *     estimate than the raw noisy result, using nothing more than a linear fit
 *     across noise scales.
 *
 * NOTE: VQE parameters used here are still the *noiseless*-optimal ones. This
 * example is about evaluating a fixed circuit under noise and mitigating that
 * noise, not about re-optimizing parameters in the presence of noise.
 * TODO: noise-aware VQE  optimization
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
#include "../physics/vqe.h"
#include "../physics/vqe_noisy.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > H2 Noisy VQE: T1/T2 Decoherence and Zero-Noise Extrapolation\n\n");

  double R_bond = 1.4;
  double A[3] = {0, 0, 0}, B[3] = {0, 0, R_bond};
  basis_function_t *funcs[2] = {molint_basis_sto3g_h(A),
                                molint_basis_sto3g_h(B)};
  const double charge[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R_bond}};
  molecule_t *mol = molecule_alloc(2, charge, centers);

  printf("Step 1: STO-3G H2 RHF, AO -> MO transform, Jordan-Wigner "
         "Hamiltonian (4 qubits, 16x16) -- same pipeline as eg_42\n\n");
  molecular_hf_result_t *hf = molecular_rhf(funcs, 2, mol, 2, 1e-10, 200);
  cmatrix_t *Hcore = molecular_core_hamiltonian(funcs, 2, mol);
  double *eri_ao = molecular_eri_tensor(funcs, 2);
  double nuclear_repulsion = molecule_nuclear_repulsion(mol);
  double h_mo[4];
  double *eri_mo = malloc(16 * sizeof(double));
  molecular_ao_to_mo(Hcore, eri_ao, hf->C, 2, h_mo, eri_mo);
  cmatrix_t *H = second_quant_build_molecular_hamiltonian(2, h_mo, eri_mo,
                                                          nuclear_repulsion);

  cmatrix_t *H_copy = cmatrix_alloc(16, 16);
  for (int i = 0; i < 16 * 16; i++) {
    H_copy->data[i] = H->data[i];
  }

  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  double E_fci = eig->eigenvalues[0];
  printf("  RHF: %.6f Hartree   FCI (exact): %.6f Hartree\n\n",
         hf->total_energy, E_fci);

  printf("Step 2: Noiseless VQE (4 qubits, hardware-efficient ansatz) -- "
         "same as eg_42\n\n");
  vqe_result_t vqe_res = vqe_run(4, 3, H, 8, 0.6, 20260810ULL);
  printf("  Noiseless VQE: %.6f Hartree  (error vs FCI: %.4f mHartree)\n\n",
         vqe_res.energy, 1000.0 * fabs(vqe_res.energy - E_fci));

  printf("Step 3: Evaluate the SAME ansatz parameters under increasing T1/T2 "
         "noise (density-matrix simulation, physics/vqe_noisy.c)\n\n");
  double gate_time = 0.05; // natural units; a nominal per-gate duration
  const double noise_levels[4] = {0.0, 0.002, 0.005,
                                  0.008}; // \gamma1 (T1 rate)
  printf("  %-12s %-16s %-16s\n", "gamma1", "noisy energy",
         "error vs FCI (mH)");
  for (int i = 0; i < 4; i++) {
    double gamma1 = noise_levels[i];
    double gamma2 = gamma1 * 0.6; // T2 typically faster than T1 in practice
    double e_noisy =
        vqe_noisy_energy(4, 3, vqe_res.theta_opt, H, gamma1, gamma2, gate_time);
    printf("  %-12.3f %-16.6f %-16.4f\n", gamma1, e_noisy,
           1000.0 * fabs(e_noisy - E_fci));
  }
  printf("\n  Energy error grows with noise strength, even though the circuit "
         "parameters never changed.\n\n");

  printf("Step 4: Zero-noise extrapolation at a representative noise level "
         "(\\gamma1=0.005, \\gamma2=0.003)\n\n");
  double gamma1 = 0.005, gamma2 = 0.003;
  double raw_c1 = 0.0;
  double e_zne = vqe_noisy_zne_energy(4, 3, vqe_res.theta_opt, H, gamma1,
                                      gamma2, gate_time, 3, &raw_c1);
  printf("  Raw noisy energy (c=1):         %.6f Hartree  (error: %.4f "
         "mHartree)\n",
         raw_c1, 1000.0 * fabs(raw_c1 - E_fci));
  printf("  ZNE-extrapolated energy (c->0): %.6f Hartree  (error: %.4f "
         "mHartree)\n\n",
         e_zne, 1000.0 * fabs(e_zne - E_fci));

  printf("=== Summary ===\n\n");
  printf("  %-30s %-14.6f\n", "RHF (mean-field)", hf->total_energy);
  printf("  %-30s %-14.6f\n", "FCI (exact)", E_fci);
  printf("  %-30s %-14.6f\n", "VQE, noiseless", vqe_res.energy);
  printf("  %-30s %-14.6f\n", "VQE, noisy (raw, c=1)", raw_c1);
  printf("  %-30s %-14.6f\n", "VQE, noisy + ZNE", e_zne);
  printf("\n  ZNE recovered %.4f mHartree of the %.4f mHartree that noise cost "
         "the raw readout at this noise level, using only a polynomial "
         "extrapolation across 3 amplified noise scales.\n",
         1000.0 * (fabs(raw_c1 - E_fci) - fabs(e_zne - E_fci)),
         1000.0 * fabs(raw_c1 - E_fci));

  free(vqe_res.theta_opt);
  eigen_free(eig);
  cmatrix_free(H_copy);
  cmatrix_free(H);
  cmatrix_free(Hcore);
  free(eri_ao);
  free(eri_mo);
  molecular_hf_result_free(hf);
  molecule_free(mol);
  basis_function_free(funcs[0]);
  basis_function_free(funcs[1]);

  return 0;
}
