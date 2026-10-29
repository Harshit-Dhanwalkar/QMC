/*
 * QMC demo driver: runs through examples in sequence (or one at a time).
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

typedef struct {
  const char *path;
  const char *description;
} demo_entry_t;

static const demo_entry_t demos[] = {
    {"build/eg_01_particle_box", "Infinite square well (particle in a box)"},
    {"build/eg_02_harmonic", "Harmonic oscillator"},
    {"build/eg_03_finite_well", "Finite square well"},
    {"build/eg_04_infinite_well", "Infinite well, numerical vs analytic"},
    {"build/eg_05_uncertainty", "Heisenberg uncertainty (Gaussian state)"},
    {"build/eg_06_hydrogen", "Hydrogen atom radial solve"},
    {"build/eg_07_central_potential", "General 3D central potential"},
    {"build/eg_08_helium", "Helium ground state (variational)"},
    {"build/eg_09_identical_particles", "Slater determinants / identical particles"},
    {"build/eg_10_perturbation", "Perturbation theory"},
    {"build/eg_11_wkb", "WKB approximation"},
    {"build/eg_12_tunnelling", "Barrier tunnelling"},
    {"build/eg_13_scattering", "Phase shifts and Born approximation"},
    {"build/eg_14_rabi", "Rabi oscillations"},
    {"build/eg_15_angular_coupling", "Angular momentum coupling (CG coefficients)"},
    {"build/eg_16_finestructure", "Hydrogen fine structure"},
    {"build/eg_17_dirac", "1D Dirac equation"},
    {"build/eg_18_qubits", "Multi-qubit states and entanglement"},
    {"build/eg_19_lindblad", "Lindblad open-system evolution"},
    {"build/eg_20_hartree_fock", "Hartree-Fock (closed-shell)"},
    {"build/eg_21_driven", "Driven two-level systems / Landau-Zener"},
    {"build/eg_22_soft", "Split-operator (SOFT) wavepacket evolution"},
    {"build/eg_23_fermi_golden_rule", "Fermi's golden rule"},
    {"build/eg_24_central_potential_3d", "3D central potential (HO, finite well)"},
    {"build/eg_25_boson_sampling", "Boson sampling / Hong-Ou-Mandel"},
    {"build/eg_26_zeeman", "Zeeman effect"},
    {"build/eg_27_cap_tdse", "Complex absorbing potential + time-dependent V(x,t)"},
    {"build/eg_28_vmc_helium", "Variational Monte Carlo (helium)"},
    {"build/eg_29_dmc_helium", "Diffusion Monte Carlo (helium)"},
    {"build/eg_30_pimc_helium", "Path Integral Monte Carlo (helium)"},
    {"build/eg_31_vqe", "Variational Quantum Eigensolver"},
    {"build/eg_32_mp2", "MP2 correlation energy"},
    {"build/eg_33_lattice", "Tight-binding lattice models (bands, Anderson, SSH)"},
    {"build/eg_34_quantum_info", "Teleportation, superdense coding, Bell's inequality"},
    {"build/eg_35_qec", "3-qubit quantum error correction (bit-flip, phase-flip)"},
    {"build/eg_36_shor_code", "9-qubit Shor code: correcting arbitrary single-qubit errors (X, Y, Z)"},
    {"build/eg_37_qec5", "5-qubit perfect error correcting code"},
    {"build/eg_38_qec_surface17", "Distance-3 rotated surface code"},
    {"build/eg_39_quantum_algorithms", "Deutsch-Jozsa and Grover's search algorithm"},
    {"build/eg_40_second_quant", "Second quantization and Jordan-Wigner transformation"},
    {"build/eg_41_landau_levels", "Landau levels on a lattice via Peierls substitution"},
    {"build/eg_42_openmp_qmc", "OpenMP-parallel VMC/DMC/PIMC replica runners (wall-clock speedup)"},
    {"build/eg_44_molecular_integrals", "General molecular integrals (McMurchie-Davidson): H2/STO-3G restricted Hartree-Fock"},
    {"build/eg_44_dft_atoms", "Kohn-Sham LDA DFT for atoms (He, Be), compared against Hartree-Fock"},
    {"build/eg_45_h2_vqe", "Real H2 VQE: AO integrals -> Jordan-Wigner -> variational quantum eigensolver (VQE)"},
    {"build/eg_46_h4_vqe", "H4 chain: general N-basis RHF, FCI, and VQE at 8 qubits"},
    {"build/eg_47_lih_vqe", "LiH: p-orbitals, frozen-core reduction, and VQE"},
    {"build/eg_48_noisy_vqe", "H2 Noisy VQE: T1/T2 decoherence and ZNE"},
    {"build/eg_49_geometry_optimization", "Molecular geometry optimization via HF/DFT gradients"},
    {"build/eg_50_ccsd", "Coupled-Cluster Singles and Doubles (CCSD)"},
    {"build/eg_51_ccsd_t", "Coupled-Cluster Singles, Doubles, and Triples CCSD(T)"},
    {"build/eg_52_general_basis_parser", "General Gaussian basis set parser and integral evaluator"},
    {"build/eg_53_ucc", "Unitary Coupled-Cluster (UCCSD) ansatz generation for VQE"},
    {"build/eg_54_qpe", "Quantum Phase Estimation (QPE) for molecular Hamiltonians"},
    {"build/eg_55_molecular_dft", "Kohn-Sham molecular DFT solver (LDA/VWN)"},
    {"build/eg_56_molecular_pbe_dft", "Molecular DFT with GGA functionals (PBE)"},
    {"build/eg_57_molecular_mp2", "General molecular MP2 (beyond s-orbitals): LiH with/without frozen core"},
    {"build/eg_58_tddft", "Linear-response TDDFT excitation energies"},
    {"build/eg_59_schrodinger_solvers", "Comprehensive 1D/3D Schrödinger eigensolvers comparison"},
    {"build/eg_60_fci", "Full Configuration Interaction (FCI) exact diagonalization"},
    {"build/eg_61_cisd", "Configuration interaction singles and doubles (CISD)"},
    {"build/eg_62_spin_chain_dsf", "Dynamic structure factor for 1D Heisenberg spin chains"},
    {"build/eg_63_ising_chain", "Transverse-field Ising chain exact diagonalization and Jordan-Wigner"},
    {"build/eg_64_uhf_openshell", "Unrestricted Hartree-Fock: open-shell radicals (H, Li doublet, Li+)"},
    {"build/eg_65_ump2", "UMP2: post-HF correlation energy for open-shell radicals (Li doublet)"},
    {"build/eg_66_uhf_geometry_optimization", "Analytic UHF nuclear gradient: open-shell geometry optimization (LiH+)"},
    {"build/eg_67_vqe_pauli_hamiltonian", "VQE with a general Pauli-string Hamiltonian: the XXZ chain"},
    {"build/eg_68_vibrational_analysis", "Harmonic vibrational frequency analysis via numerical Hessian (H2)"},
    {"build/eg_69_dmrg", "Infinite-system DMRG: Heisenberg antiferromagnetic chain convergence"},
    {"build/eg_70_finite_dmrg", "Finite-system DMRG sweeps"},
    {"build/eg_71_casscf", "Complete Active Space SCF (CASSCF)"},
    {"build/eg_72_qaoa", "Quantum Approximate Optimization Algorithm (QAOA)"},
    {"build/eg_73_floquet", "Floquet theory for periodically driven systems"},
    {"build/eg_74_td_perturbation", "Time-dependent perturbation theory"},
    {"build/eg_75_tight_binding", "Generic multi-band tight-binding models"},
    {"build/eg_76_dirac_evolve", "1D Dirac equation time evolution"},
    {"build/eg_77_orbital_sample", "Monte Carlo sampling of hydrogen orbitals"},
    {"build/eg_78_bloch", "Single-qubit Bloch-vector picture"},
    {"build/eg_79_anderson", "Anderson localisation on a 1D chain"},
    {"build/eg_80_quantum_walk", "Discrete-time quantum walk vs classical random walk"},
    {"build/eg_81_landau_zener", "Landau-Zener transitions at an avoided crossing"},
    {"build/eg_82_kicked_rotor", "Quantum kicked rotor: dynamical localisation"},
    {"build/eg_83_ssh_chain", "Su-Schrieffer-Heeger chain topological insulator"},
    {"build/eg_84_tfim_quench", "Exact quench dynamics of the transverse-field Ising chain"},
    {"build/eg_85_haldane_ribbon", "Haldane Chern insulator ribbon and phase diagram"},
    {"build/eg_86_phase_space", "Phase-space picture of a bosonic mode"},
    {"build/eg_87_grover_multi", "Grover search with several marked items"},
    {"build/eg_88_bell_chsh", "CHSH Bell test on a Werner state"},
    {"build/eg_89_box_revival", "Particle in a box: wavefunction revivals"},
    {"build/eg_90_zeno", "Quantum Zeno effect"},
    {"build/eg_91_circular_billiard", "Quantum billiard in a disk"},
    {"build/eg_92_hubbard_dimer", "Two-site Hubbard model exact diagonalisation"},
    {"build/eg_93_ab_ring", "Aharonov-Bohm ring: persistent currents and transmission"},
    {"build/eg_94_double_barrier", "Resonant tunneling through a double barrier"},
    {"build/eg_95_josephson", "RCSJ Josephson junction phase dynamics"},
    {"build/eg_96_double_well", "Symmetric double-well tunnelling and parity doublet"},
    {"build/eg_97_zeeman_full", "Zeeman effect from weak field to Paschen-Back"},
    {"build/eg_latex_gen", "Automated LaTeX matrix and equation code generator"},
};

static const int n_demos = (int)(sizeof(demos) / sizeof(demos[0]));

static int binary_exists(const char *path) {
  struct stat st;

  return stat(path, &st) == 0;
}

static void run_one(int idx) {
  if (idx < 0 || idx >= n_demos) {
    printf("Invalid selection.\n");

    return;
  }

  const demo_entry_t *d = &demos[idx];
  if (!binary_exists(d->path)) {
    printf("  [skipped] %s not built yet - run `make` first.\n", d->path);

    return;
  }

  printf("\n========================================\n");
  printf(" %2d. %s\n", idx + 1, d->description);
  printf("========================================\n");
  int rc = system(d->path);
  if (rc != 0) {
    printf("  (exited with status %d)\n", rc);
  }
}

static void run_all(void) {
  for (int i = 0; i < n_demos; i++) {
    run_one(i);
  }
}

static void print_menu(void) {
  printf("\nQMC demo driver (plot backend: %s)\n", QMC_PLOT_BACKEND_NAME);
  printf("--------------------------------------------------\n");
  for (int i = 0; i < n_demos; i++) {
    printf("  %2d. %s\n", i + 1, demos[i].description);
  }
  printf("   a. Run all\n");
  printf("   q. Quit\n");
}

int main(void) {
  char line[64];

  for (;;) {
    print_menu();
    printf("\nSelect: ");
    if (!fgets(line, sizeof line, stdin)) {
      break;
    }

    if (line[0] == 'q' || line[0] == 'Q') {
      break;
    }

    if (line[0] == 'a' || line[0] == 'A') {
      run_all();
      continue;
    }

    int choice = atoi(line);
    if (choice >= 1 && choice <= n_demos) {
      run_one(choice - 1);
    } else {
      printf("Not a valid choice.\n");
    }
  }

  return 0;
}
