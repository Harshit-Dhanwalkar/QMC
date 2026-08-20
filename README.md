# Quantum Mechanics in C

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-PolyForm__Noncommercial__1.0.0-blue.svg" alt="License">
  </a>
  <a href="https://github.com/Harshit-Dhanwalkar/QMC/actions/workflows/docs.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/Harshit-Dhanwalkar/QMC/docs.yml?branch=main&label=docs" alt="Docs Deploy Status">
  </a>
  <a href="https://github.com/Harshit-Dhanwalkar/QMC/actions/workflows/ci.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/Harshit-Dhanwalkar/QMC/ci.yml?branch=main&label=build_and_tests" alt="CI Status">
  </a>
  <img src="https://img.shields.io/badge/c-C99-blue.svg" alt="C Standard">
  <img src="https://img.shields.io/badge/parallel-OpenMP-informational.svg" alt="OpenMP">
  <img src="https://img.shields.io/badge/sanitizers-ASan%20%7C%20UBSan-green.svg" alt="Sanitizers">
  <img src="https://img.shields.io/badge/platform-Linux-lightgrey.svg" alt="OS Linux">
</p>

A pure-C library and simulation engine for numerical and semi-analytic quantum mechanics, covering undergraduate through early-graduate physics: wavefunctions, eigensolvers, time evolution, perturbation theory, scattering, angular momentum coupling, identical particles, open quantum systems, relativistic wave equations (Klein-Gordon, Dirac), Hartree-Fock self-consistent field theory, and quantum Monte Carlo (variational and diffusion).

> **Units Convention:** Natural/atomic units ($\hbar = m = 1$) are used throughout, except `hydrogen.c` and `fine_structure.c`, which work in SI units and take explicit physical parameters ($\hbar, m_e, e$, etc.).

## Prerequisites

- C compiler with C99 and OpenMP support (`gcc` or `clang`)
- `make` and `cmake`
- `grplot` : GR Framework (Default plot backend, built via Git submodule or system libs)
- `gnuplot` : GNU plot (Optional, fallback plot backend)
- MIP :
- LAPACK :

## Building

```bash
# Clone repository with submodules
git clone --recursive https://github.com/Harshit-Dhanwalkar/QMC.git
cd QMC

# Build GR submodule (if not using system GR)
cd third_party/gr
cmake -B build -DGR_INSTALL=ON \
  -DCMAKE_INSTALL_PREFIX=$(pwd)/install \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt6=ON \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt5=ON \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt4=ON
cmake --build build -j
cmake --install build
cd ../..

# Build QMC demo driver and binaries
make all

# Run interactive demo selector
make demo

# Run test suite
make run-tests
# All example/test output is written under `output/` (`QMC_OUTPUT_DIR`)

# Cleanup
make clean
```

### Optional: LAPACK/BLAS acceleration

Dense Hermitian diagonalization (used by `hartree_fock.c`, `dft.c`, `molecular_hf.c`, exact diagonalization of `second_quant_build_molecular_hamiltonian`'s output, and `vqe.c`) uses a hand-rolled QR/real-embedding solver by default. This becomes the real bottleneck as system size grows (e.g. a radial-grid DFT run at larger N, or an 8$\pm$qubit molecular Hamiltonian's exact diagonalization).

If `liblapacke-dev` and a BLAS backend (e.g. `libopenblas-dev`) are installed, opt in to LAPACK-backed solvers:

```bash
# Debian/Ubuntu:
sudo apt-get install liblapacke-dev libopenblas-dev libtmglib-dev

# make USE_LAPACK=1 run-tests
PATH=/usr/bin:$PATH make USE_LAPACK=1 run-tests

# make USE_LAPACK=1 run-examples
PATH=/usr/bin:$PATH make USE_LAPACK=1 run-examples
```

This swaps in LAPACK's `dsyev`/`zheev` for both the real-symmetric and complex-Hermitian solvers for the complex case.

### Optional: disabling AddressSanitizer

AddressSanitizer (`-fsanitize=address`) is on by default. It must be turned OFF (for `EXTRA_CFLAGS` to be cleared and that hook is additive and doesn't remove flags already in `CFLAGS`) when running binaries under Valgrind, since Valgrind cannot run an ASan-instrumented binary correctly:

```bash
# install headers + libraries to a prefix (default /usr/local)
make SANITIZE=0 install-lib PREFIX=/usr/local
```

## Benchmarks

Standalone benchmark programs under `benchmarks/`, built the same way as examples/tests:

```bash
# Build and run everything except the MPI hybrid benchmark
make run-benchmarks

# Just build them (no run)
make benchmarks
```

- **`bench_accuracy`** : VMC/DMC ground-state energy vs. published values for the He isoelectronic two-electron sequence (He, Li+, Be2+).
- **`bench_eigensolver`** : hand-rolled (Jacobi) vs. LAPACK complex-Hermitian eigensolver - these are compared by building twice, once with `USE_LAPACK=0` and once `USE_LAPACK=1` (see `benchmarks/README.md`).
- **`bench_vmc_convergence`** : VMC statistical convergence vs. sample count.
- **`bench_reference_vs_blocked`** / **`instrument_iteration_counts`** : cache-blocking benchmark for the tridiagonal QL eigensolver -built and run standalone with `gcc`/`valgrind`+`cachegrind`, not through `make`.

`run-benchmarks` writes each benchmark's stdout to `output/bench_<name>.dat` in addition to printing it, for later plotting/comparison.

See [[benchmarks/README.md]] for build/run details, validation notes, and the reasoning behind each benchmark's design (e.g. why `bench_accuracy` uses Be2+ instead of neutral Be, why `bench_eigensolver` needs two separate builds instead of a runtime flag).

---

## Features & Modules

- Core Math : Complex vectors/matrices, special functions (Hermite, Laguerre, Legendre, Spherical Harmonics, Bessel)
  - ODE solver : Numerov integrator, RK4, Crank-Nicolson propagator, FFT (1D, 2D, 3D), and custom eigensolvers (`tridiag_eigh`, `complex_eigh`).
- 1D & 3D Potentials: Particle in a box, finite square well, harmonic oscillator, tunneling barriers, radial hydrogen, general 3D central potentials.
- Quantum Dynamics: Time-dependent Schr$\ddot{o}$dinger equation, Split-Operator Fourier Transform (SOFT 2D/3D), driven two-level systems, Rabi oscillations, and Complex Absorbing Potentials (CAP).
- Many-Body & Chemistry: Identical particles (Slater determinants), closed-shell Hartree-Fock SCF, and Quantum Monte Carlo (VMC and DMC for two-electron atoms/ions with a Slater-Jastrow ansatz, and finite-temperature PIMC for helium with Kelbg-regularized Coulomb pair actions and bisection sampling).
- Open Quantum Systems & Quantum Info: Lindblad master equation for dissipative density matrices, multi-qubit state evolution, entanglement, Boson sampling, and a Variational Quantum Eigensolver (hardware-efficient ansatz + classical coordinate-descent optimizer).
- Relativistic QM: 1D/3D Dirac equation and Klein-Gordon solver.
- Visualization: Multi-backend plot abstraction supporting GR Framework (default), GNUplot pipe, or Matplotlib (Python subprocess) pipe.
- Optimised linear algebra (LAPACK/BLAS backend).

<details>
<summary>Examples for each topic.</summary>
<ul>
  <li><code>eg_01_particle_box.c</code></li>
  <li><code>eg_02_harmonic.c</code></li>
  <li><code>eg_03_finite_well.c</code></li>
  <li><code>eg_04_infinite_well.c</code></li>
  <li><code>eg_05_uncertainty.c</code></li>
  <li><code>eg_06_hydrogen.c</code></li>
  <li><code>eg_07_central_potential.c</code></li>
  <li><code>eg_08_helium.c</code></li>
  <li><code>eg_09_identical_particles.c</code></li>
  <li><code>eg_10_perturbation.c</code></li>
  <li><code>eg_11_wkb.c</code></li>
  <li><code>eg_12_tunnelling.c</code></li>
  <li><code>eg_13_scattering.c</code></li>
  <li><code>eg_14_rabi.c</code></li>
  <li><code>eg_15_angular_coupling.c</code></li>
  <li><code>eg_16_finestructure.c</code></li>
  <li><code>eg_17_dirac.c</code></li>
  <li><code>eg_18_qubits.c</code></li>
  <li><code>eg_19_lindblad.c</code></li>
  <li><code>eg_20_hartree_fock.c</code></li>
  <li><code>eg_21_driven.c</code></li>
  <li><code>eg_22_soft.c</code></li>
  <li><code>eg_23_fermi_golden_rule.c</code></li>
  <li><code>eg_24_central_potential_3d.c</code></li>
  <li><code>eg_25_boson_sampling.c</code></li>
  <li><code>eg_26_zeeman.c</code></li>
  <li><code>eg_27_cap_tdse.c</code></li>
  <li><code>eg_28_vmc_helium.c</code></li>
  <li><code>eg_29_dmc_helium.c</code></li>
  <li><code>eg_30_pimc_helium.c</code></li>
  <li><code>eg_31_vqe.c</code></li>
  <li><code>eg_32_mp2.c</code></li>
  <li><code>eg_33_lattice.c</code></li>
  <li><code>eg_34_quantum_info.c</code></li>
  <li><code>eg_35_qec.c</code></li>
  <li><code>eg_36_quantum_algorithms.c</code></li>
  <li><code>eg_37_second_quant.c</code></li>
  <li><code>eg_38_landau_levels.c</code></li>
  <li><code>eg_39_openmp_qmc</code></li>
  <li><code>eg_40_molecular_integrals.c</code></li>
  <li><code>eg_41_dft_atoms.c</code></li>
  <li><code>eg_42_h2_vqe.c</code></li>
  <li><code>eg_43_h4_vqe.c</code></li>
  <li><code>eg_44_lih_vqe.c</code></li>
  <li><code>eg_45_vqe_nosiy.c</code></li>
  <li><code>eg_46_geometry_optimization.c</code></li>
</ul>
</details>

<details>
<summary>Unit tests for each topic.</summary>
<ul>
  <li><code>test_angular_coupling.c</code></li>
  <li><code>test_boson_sampling.c</code></li>
  <li><code>test_ccsd.c</code></li>
  <li><code>test_central_potential.c</code></li>
  <li><code>test_central_potential_3d.c</code></li>
  <li><code>test_complex.c</code></li>
  <li><code>test_complex_eigh.c</code></li>
  <li><code>test_crank_nicolson.c</code></li>
  <li><code>test_dirac.c</code></li>
  <li><code>test_dmc.c</code></li>
  <li><code>test_driven.c</code></li>
  <li><code>test_fermi_golden_rule.c</code></li>
  <li><code>test_fft.c</code></li>
  <li><code>test_fine_structure.c</code></li>
  <li><code>test_grplot.c</code></li>
  <li><code>test_hartree_fock.c</code></li>
  <li><code>test_helium.c</code></li>
  <li><code>test_hydrogen.c</code></li>
  <li><code>test_identical.c</code></li>
  <li><code>test_klein_gordon.c</code></li>
  <li><code>test_lanczos.c</code></li>
  <li><code>test_lattice.c</code></li>
  <li><code>test_lindblad.c</code></li>
  <li><code>test_matrix.c</code></li>
  <li><code>test_mp2.c</code></li>
  <li><code>test_numerov.c</code></li>
  <li><code>test_perturbation.c</code></li>
  <li><code>test_pimc.c</code></li>
  <li><code>test_potentials.c</code></li>
  <li><code>test_quantum_algorithms.c</code></li>
  <li><code>test_qec.c</code></li>
  <li><code>test_qubits.c</code></li>
  <li><code>test_quantum_info.c</code></li>
  <li><code>test_rabi.c</code></li>
  <li><code>test_random.c</code></li>
  <li><code>test_rk4.c</code></li>
  <li><code>test_scattering.c</code></li>
  <li><code>test_second_quant.c</code></li>
  <li><code>test_molecular_integrals.c</code></li>
  <li><code>test_dft.c</code></li>
  <li><code>test_h2_vqe.c</code></li>
  <li><code>test_molecular_hf.c</code></li>
  <li><code>test_hf_gradient.c</code></li>
  <li><code>test_uhf.c</code></li>
  <li><code>test_lih.c</code></li>
  <li><code>test_soft.c</code></li>
  <li><code>test_tridiag.c</code></li>
  <li><code>test_tridiag_eigvals.c</code></li>
  <li><code>test_vmc.c</code></li>
  <li><code>test_vqe.c</code></li>
  <li><code>test_vqe_noisy.c</code></li>
  <li><code>test_wkb.c</code></li>
  <li><code>test_zeeman.c</code></li>
</ul>
</details>

## Reference texts

- Sakurai & Napolitano, Bransden & Joachain, Szabo & Ostlund, Thijssen's
- Computational Physics, Breuer & Petruccione, Feynman & Hibbs, Taylor's
- Scattering Theory, Greiner's Relativistic Quantum Mechanics.
- Quantum Theory of Many-Particle Systems, Fetter & Walecka.
- Quantum Computation and Quantum Information, Michael A. Nielsen and Isaac L. Chuang

# License and Commercial Use

This project is licensed under the [PolyForm Noncommercial License 1.0.0](LICENSE).
It is free for personal use, education, research, and non-profit organizations.

## Commercial Use Requires Permission

If you wish to use this software, or any part of its source code, for commercial or monetary purposes, you **must** obtain a commercial license.

To discuss commercial licensing, please contact me directly at:

- **Email:** [Contact Me via Email](mailto:harshitpd1729@gmail.com)
- **GitHub:** [Open a Licensing Issue](https://github.com/Harshit-Dhanwalkar/QMC/issues)
