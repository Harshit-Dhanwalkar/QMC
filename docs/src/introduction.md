# qmc

Quantum mechanics simulation and visualization engine written in pure C.

This book documents both the physics and the implementation - derivations
sit alongside the C code that computes them. The goal is a self-contained
engine that runs from 1D potentials up to relativistic wave equations, with
publication-quality plots where math annotations render as actual symbols,
not ASCII approximations.

## What this covers

The physics spans a standard graduate QM curriculum:

- Wave functions, normalization, expectation values
- 1D, 2D, and 3D potentials - infinite well, finite well, harmonic oscillator,
  tunneling barriers, Kronig-Penney, Morse, double well
- Time evolution via Crank-Nicolson (TDSE)
- Angular momentum, spin, spherical harmonics
- Hydrogen atom - radial wave functions, energy levels, orbital plots
- Perturbation theory - Stark effect, Zeeman effect, anharmonic oscillator
- WKB approximation and Bohr-Sommerfeld quantization
- Scattering theory - Born approximation, phase shifts, cross-sections
- Relativistic wave equations - Klein-Gordon and Dirac
- Identical particles - Slater determinants, Pauli exclusion

## How this book is organized

**Physics** sections derive the equations first, then show the C code that
implements them and the plots it produces. Each section corresponds directly
to a module in `src/physics/`.

**Internals** sections cover the numerical methods in isolation - the
Numerov integrator, tridiagonal eigensolver, Crank-Nicolson propagator,
and GR plotting wrapper. These are the building blocks reused across every
physics module.

## Status

Early stage. Foundation is being laid. Pages will be filled in as each
module is implemented.
