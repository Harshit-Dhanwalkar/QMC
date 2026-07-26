# Quantum Mechanics in C

A pure-C library for numerical and semi-analytic quantum mechanics, covering undergraduate through early-graduate topics: wavefunctions, eigensolvers, time evolution, perturbation theory, scattering, angular momentum coupling, identical particles, open quantum systems, and relativistic wave equations.

> Natural/atomic units ($\hbar = m = 1$) are used throughout, except `hydrogen.c` and `fine_structure.c`, which work in SI units and take `hbar`, `mass`, etc. as explicit parameters.

## Building

- Build with `make`.
- Run all tests with `make run-tests`.
- Run all examples with `make run-examples`.
  All example/test output is written under `output/` (`QMC_OUTPUT_DIR`);
- `make clean` removes it.

## Features

- Core: complex numbers, vectors, matrices, special functions (Hermite, Laguerre, Legendre, Bessel).
- ODE solvers: Numerov, RK4, Crank‑Nicolson.
- 1D potentials: infinite well, harmonic, barrier, Coulomb.
- Plotting via GR framework, GNUplot pipe and matplotlib pipe for python.
- [ ] Optimised linear algebra (LAPACK/BLAS backend).

- <details>
  <summary>Unit examples for each topic.</summary>
  <ul>
    <li>`eg_01_particle_box.c`</lo>
    <li> `eg_02_harmonic.c`</li>
    <li> `eg_03_finite_well.c` </li>
    <li> `eg_04_infinite_well.c` </li>
    <li> `eg_05_uncertainty.c` </li>
    <li> `eg_06_hydrogen.c` </li>
    <li> `eg_07_central_potential.c` </li>
    <li> `eg_08_helium.c` </li>
    <li> `eg_09_identical_particles.c` </li>
    <li> `eg_10_perturbation.c` </li>
    <li> `eg_11_wkb.c` </li>
    <li> `eg_12_tunnelling.c` </li>
    <li> `eg_13_scattering.c` </li>
    <li> `eg_14_rabi.c` </li>
    <li> `eg_15_angular_coupling.c` </li>
    <li> `eg_16_finestructure.c` </li>
    <li> `eg_17_dirac.c` </li>
    <li> `eg_18_qubits.c` </li>
    <li> `eg_19_lindblad.c` </li>
    <li> `eg_20_hartree_fock.c` </li>
    <li> `eg_21_driven.c` </li>
    <li> `eg_22_soft.c` </li>
    <li> `eg_23_fermi_golden_rule.c` </li>
    <li> `eg_24_central_potential_3d.c` </li>
    <li> `eg_25_boson_sampling.c` </li>
  </ul>
- <details>
    <summary>Unit tests for each topic.</summary>
    <ul>
    <li> `test_angular_coupling.c`</li>
    <li> `test_boson_sampling.c`</li>
    <li> `test_central_potential.c`</li>
    <li> `test_central_potential_3d.c`</li>
    <li> `test_complex.c`</li>
    <li> `test_complex_eigh.c`</li>
    <li> `test_crank_nicolson.c`</li>
    <li> `test_dirac.c`</li>
    <li> `test_driven.c`</li>
    <li> `test_fermi_golden_rule.c`</li>
    <li> `test_fft.c`</li>
    <li> `test_fine_structure.c`</li>
    <li> `test_grplot.c`</li>
    <li> `test_hartree_fock.c`</li>
    <li> `test_helium.c`</li>
    <li> `test_hydrogen.c`</li>
    <li> `test_identical.c`</li>
    <li> `test_lindblad.c`</li>
    <li> `test_matrix.c`</li>
    <li> `test_numerov.c`</li>
    <li> `test_perturbation.c`</li>
    <li> `test_potentials.c`</li>
    <li> `test_qubits.c`</li>
    <li> `test_rabi.c`</li>
    <li> `test_rk4.c`</li>
    <li> `test_scattering.c`</li>
    <li> `test_soft.c`</li>
    <li> `test_tridiag.c`</li>
    <li> `test_wkb.c`</li>
  </ul>
  </details>

## Reference texts

- Sakurai & Napolitano, Bransden & Joachain, Szabo & Ostlund, Thijssen's
- _Computational Physics_, Breuer & Petruccione, Feynman & Hibbs, Taylor's
- _Scattering Theory_, Greiner's _Relativistic Quantum Mechanics_.

# License and Commercial Use

This project is licensed under the **PolyForm Noncommercial License 1.0.0**. 
It is free for personal use, education, research, and non-profit organizations.

## Commercial Use Requires Permission
If you wish to use this software, or any part of its source code, for commercial or  monetary purposes, you **must** obtain a commercial license. 

To discuss commercial licensing, please contact me directly at:
* **Email:** [Contact Me via Email](mailto:harshitpd1729@://gmail.com)
* **GitHub:** [Open a Licensing Issue](https://github.com)
