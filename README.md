# Quantum Mechanics in C

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-PolyForm_Noncommercial_1.0.0-blue.svg" alt="License">
  </a>
  <a href="https://github.com/Harshit-Dhanwalkar/QMC/actions/workflows/docs.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/Harshit-Dhanwalkar/QMC/docs.yml?branch=main&label=docs" alt="Docs Deploy Status">
  </a>
  <a href="https://github.com/Harshit-Dhanwalkar/QMC/actions/workflows/ci.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/Harshit-Dhanwalkar/QMC/ci.yml?branch=main&label=build%20%26%20tests" alt="CI Status">
  </a>
</p>

[![SPDX-License-Identifier](https://img.shields.io/badge/License-PolyForm_Noncommercial_1.0.0-blue.svg)](LICENSE)

A pure-C library and simulation engine for numerical and semi-analytic quantum mechanics, covering undergraduate through early-graduate physics: wavefunctions, eigensolvers, time evolution, perturbation theory, scattering, angular momentum coupling, identical particles, open quantum systems, relativistic wave equations (Klein-Gordon, Dirac), Hartree-Fock self-consistent field theory, and quantum Monte Carlo (variational and diffusion).

> **Units Convention:** Natural/atomic units ($\hbar = m = 1$) are used throughout, except `hydrogen.c` and `fine_structure.c`, which work in SI units and take explicit physical parameters ($\hbar, m_e, e$, etc.).

## Prerequisites

- C compiler with C99 and OpenMP support (`gcc` or `clang`)
- `make` and `cmake`
- `grplot` : GR Framework (Default plot backend, built via Git submodule or system libs)
- `gnuplot` : GNU plot (Optional, fallback plot backend)

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
# All example/test output is written under `output/` (`QMC_OUTPUT_DIR`);

# Cleanup
make clean
```

## Features & Modules

- Core Math : Complex vectors/matrices, special functions (Hermite, Laguerre, Legendre, Spherical Harmonics, Bessel)
  - ODE solver : Numerov integrator, RK4, Crank-Nicolson propagator, FFT (1D, 2D, 3D), and custom eigensolvers (`tridiag_eigh`, `complex_eigh`).
- 1D & 3D Potentials: Particle in a box, finite square well, harmonic oscillator, tunneling barriers, radial hydrogen, general 3D central potentials.
- Quantum Dynamics: Time-dependent Schr$\ddot{o}$dinger equation, Split-Operator Fourier Transform (SOFT 2D/3D), driven two-level systems, Rabi oscillations, and Complex Absorbing Potentials (CAP).
- Many-Body & Chemistry: Identical particles (Slater determinants), closed-shell Hartree-Fock SCF, and Quantum Monte Carlo (VMC and DMC for two-electron atoms/ions with a Slater-Jastrow ansatz, and finite-temperature PIMC for helium with Kelbg-regularized Coulomb pair actions and bisection sampling).
- Open Quantum Systems & Quantum Info: Lindblad master equation for dissipative density matrices, multi-qubit state evolution, entanglement, Boson sampling, and a Variational Quantum Eigensolver (hardware-efficient ansatz + classical coordinate-descent optimizer).
- Relativistic QM: 1D/3D Dirac equation and Klein-Gordon solver.
- Visualization: Multi-backend plot abstraction supporting GR Framework (default), GNUplot pipe, or Matplotlib (Python subprocess) pipe.
- [ ] Optimised linear algebra (LAPACK/BLAS backend).

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
</ul>
</details>

<details>
<summary>Unit tests for each topic.</summary>
<ul>
  <li><code>test_angular_coupling.c</code></li>
  <li><code>test_boson_sampling.c</code></li>
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
  <li><code>test_lindblad.c</code></li>
  <li><code>test_matrix.c</code></li>
  <li><code>test_numerov.c</code></li>
  <li><code>test_perturbation.c</code></li>
  <li><code>test_pimc.c</code></li>
  <li><code>test_potentials.c</code></li>
  <li><code>test_qubits.c</code></li>
  <li><code>test_rabi.c</code></li>
  <li><code>test_random.c</code></li>
  <li><code>test_rk4.c</code></li>
  <li><code>test_scattering.c</code></li>
  <li><code>test_soft.c</code></li>
  <li><code>test_tridiag.c</code></li>
  <li><code>test_tridiag_eigvals.c</code></li>
  <li><code>test_vmc.c</code></li>
  <li><code>test_vqe.c</code></li>
  <li><code>test_wkb.c</code></li>
  <li><code>test_zeeman.c</code></li>
</ul>
</details>

## Reference texts

- Sakurai & Napolitano, Bransden & Joachain, Szabo & Ostlund, Thijssen's
- Computational Physics, Breuer & Petruccione, Feynman & Hibbs, Taylor's
- Scattering Theory, Greiner's Relativistic Quantum Mechanics.

# License and Commercial Use

This project is licensed under the [PolyForm Noncommercial License 1.0.0](LICENSE).
It is free for personal use, education, research, and non-profit organizations.

## Commercial Use Requires Permission

If you wish to use this software, or any part of its source code, for commercial or monetary purposes, you **must** obtain a commercial license.

To discuss commercial licensing, please contact me directly at:

- **Email:** [Contact Me via Email](mailto:harshitpd1729@gmail.com)
- **GitHub:** [Open a Licensing Issue](https://github.com/Harshit-Dhanwalkar/QMC/issues)
