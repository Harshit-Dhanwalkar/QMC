# Quantum Mechanics in C

Quantum Mechanics implementations : C library and executables for solving quantum mechanical problems.

## Features

- Core: complex numbers, vectors, matrices, special functions (Hermite, Laguerre, Legendre, Bessel).
- ODE solvers: Numerov, RK4, Crank‑Nicolson.
- 1D potentials: infinite well, harmonic, barrier, Coulomb.
- Plotting via GR (2D), GNUplot pipe and matplotlib pipe for python.

## Planned / TODO

- [ ] Full 1D TISE solver
  - [x] Diagonalisation
  - [ ] Shooting
- [x] 2D/3D systems and density of states.
- [x] Angular momentum (spherical harmonics, CG coefficients).
- [x] Hydrogen atom (radial wavefunctions).
- [x] Perturbation theory (static & time‑dependent).
- [x] Scattering theory (phase shifts, Born).
- [x] WKB approximation.
- [x] Identical particles
- [x] Slater determinants.
- [ ] Relativistic quantum mechanics
  - [ ] Dirac 1D : complex-Hermitian solver
  - [ ] Klein‑Gordon
- [ ] 3D plotting (via [GR](https://github.com/sciapp/gr.git)).
- [ ] LaTeX annotation on plots (via Pango/GR or external generation).
- [ ] Optimised linear algebra (LAPACK/BLAS backend).
- [ ] Unit examples for each topic.
  - [x] `eg_01_particle_box.c`
  - [x] `eg_02_harmonic.c`
  - [x] `eg_03_hydrogen.c`
  - [x] `eg_04_perturbation.c`
  - [x] `eg_05_tunnelling.c`
  - [x] `eg_06_finite_well.c`
  - [x] `eg_07_infinite_well.c`
- [ ] Unit tests for each topic.
  - [x] `test_complex.c`
  - [x] `test_crank_nicolson.c`
  - [x] `test_fft.c`
  - [x] `test_grplot.c`
  - [x] `test_hydrogen.c`
  - [x] `test_matrix.c`
  - [x] `test_numerov.c`
  - [x] `test_perturbation.c`
  - [x] `test_potentials.c`
  - [x] `test_rk4.c`
  - [x] `test_wkb.c`

# LICENSE

[QMC](https://github.com/Harshit-Dhanwalkar/QMC) currently is under [GNU General Public License v3.0](https://github.com/Harshit-Dhanwalkar/QMC/blob/main/LICENSE)
// TODO: Update to 'non-commerial' License.
