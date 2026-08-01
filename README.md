# Quantum Mechanics in C

A pure-C library for numerical and semi-analytic quantum mechanics, covering undergraduate through early-graduate topics: wavefunctions, eigensolvers, time evolution, perturbation theory, scattering, angular momentum coupling, identical particles, open quantum systems, relativistic wave equations (Klein-Gordon, Dirac), Hartree-Fock self-consistent field theory, and quantum Monte Carlo (variational and diffusion).

> Natural/atomic units ($\hbar = m = 1$) are used throughout, except `hydrogen.c` and `fine_structure.c`, which work in SI units and take `hbar`, `mass`, etc. as explicit parameters.

## Building

- Build with `make`.
- Run demo with `make demo`. Which will provide you with options for demo individual implementations.
- Run all tests with `make run-tests`.
- Run all examples with `make run-examples`.
  All example/test output is written under `output/` (`QMC_OUTPUT_DIR`);
- `make clean` removes it.

## Features

- Core: complex numbers, vectors, matrices, special functions (Hermite, Laguerre, Legendre, Bessel).
- ODE solvers: Numerov, RK4, Crank-Nicolson.
- 1D potentials: infinite well, finite well, harmonic, step, barrier, Coulomb, Yukawa, Morse.
- Quantum Monte Carlo: Variational (VMC) and Diffusion (DMC) Monte Carlo for the helium ground state, Slater-Jastrow trial wavefunction.
- Plotting via GR framework, GNUplot pipe, and matplotlib pipe for Python.
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
  <li><code>test_lindblad.c</code></li>
  <li><code>test_matrix.c</code></li>
  <li><code>test_numerov.c</code></li>
  <li><code>test_perturbation.c</code></li>
  <li><code>test_potentials.c</code></li>
  <li><code>test_qubits.c</code></li>
  <li><code>test_rabi.c</code></li>
  <li><code>test_random.c</code></li>
  <li><code>test_rk4.c</code></li>
  <li><code>test_scattering.c</code></li>
  <li><code>test_soft.c</code></li>
  <li><code>test_tridiag.c</code></li>
  <li><code>test_vmc.c</code></li>
  <li><code>test_wkb.c</code></li>
  <li><code>test_zeeman.c</code></li>
  <li><code>test_tridiag_eigvals.c</code></li>
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
