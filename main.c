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
    {"build/eg_15_angular_coupling",
     "Angular momentum coupling (CG coefficients)"},
    {"build/eg_16_finestructure", "Hydrogen fine structure"},
    {"build/eg_17_dirac", "1D Dirac equation"},
    {"build/eg_18_qubits", "Multi-qubit states and entanglement"},
    {"build/eg_19_lindblad", "Lindblad open-system evolution"},
    {"build/eg_20_hartree_fock", "Hartree-Fock (closed-shell)"},
    {"build/eg_21_driven", "Driven two-level systems / Landau-Zener"},
    {"build/eg_22_soft", "Split-operator (SOFT) wavepacket evolution"},
    {"build/eg_23_fermi_golden_rule", "Fermi's golden rule"},
    {"build/eg_24_central_potential_3d",
     "3D central potential (HO, finite well)"},
    {"build/eg_25_boson_sampling", "Boson sampling / Hong-Ou-Mandel"},
    {"build/eg_26_zeeman", "Zeeman effect"},
    {"build/eg_27_cap_tdse",
     "Complex absorbing potential + time-dependent V(x,t)"},
    {"build/eg_28_vmc_helium", "Variational Monte Carlo (helium)"},
    {"build/eg_29_dmc_helium", "Diffusion Monte Carlo (helium)"},
    {"build/eg_30_pimc_helium", "Path Integral Monte Carlo (helium)"},
    {"build/eg_31_vqe", "Variational Quantum Eigensolver"},
    {"build/eg_32_mp2", "MP2 correlation energy"},
    {"build/eg_33_lattice", "Tight-binding lattice models (bands, Anderson, SSH)"},
    {"build/eg_34_quantum_info", "Teleportation, superdense coding, Bell's inequality"},
    {"build/eg_35_qec", "3-qubit quantum error correction (bit-flip, phase-flip)"},
    {"build/eg_36_quantum_algorithms", "Deutsch-Jozsa and Grover's search algorithm"},
    {"build/eg_37_second_quant", "Second quantization and Jordan-Wigner transformation"},
    {"build/eg_38_landau_levels", "Landau levels on a lattice via Peierls substitution (Hofstadter model)"},
    {"build/eg_39_openmp_qmc", "OpenMP-parallel VMC/DMC/PIMC replica runners (wall-clock speedup)"},
    {"build/eg_40_molecular_integrals", "General molecular integrals (McMurchie-Davidson): H2/STO-3G restricted Hartree-Fock"},
    {"build/eg_41_dft_atoms", "Kohn-Sham LDA DFT for atoms (He, Be), compared against Hartree-Fock"},
    {"build/eg_42_h2_vqe", "Real H2 VQE: AO integrals -> Jordan-Wigner -> variational quantum eigensolver"},
    {"build/eg_43_h4_vqe", "H4 chain: general N-basis RHF, FCI, and VQE at 8 qubits"},
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
