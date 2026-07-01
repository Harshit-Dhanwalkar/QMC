# 1D Potentials

Common 1D potentials implemented in `src/physics/potentials.c` and used throughout the library.

## Implemented Potentials

All potentials are functions of position `x` and fill an array `V` of length `N`.

- **Infinite Square Well**
  `potential_infinite_well(V, x, N, L)` – zero inside `[0,L]`, infinite (1e30) outside.

- **Finite Square Well**
  `potential_finite_well(V, x, N, V0, L)` – `-V0` inside `[-L/2, L/2]`, zero outside.

- **Harmonic Oscillator**
  `potential_harmonic(V, x, N, m, omega)` – `0.5*m*omega^2*x^2`.

- **Morse Potential**
  `potential_morse(V, x, N, D_e, a, x0)` – `D_e*(1 - exp(-a*(x-x0)))^2`.

- **Double Well**
  `potential_double_well(V, x, N, V0, a)` – `-V0*((x/a)^2 - 1)^2`.

- **Kronig-Penney (periodic delta)**  
  Implemented separately for band structure calculations.

## Usage Example

```c
int N = 1000;
double *x = linspace(-10, 10, N);
double *V = malloc(N * sizeof(double));
potential_harmonic(V, x, N, 1.0, 1.0);
```

All potentials are used with the Numerov integrator or matrix diagonalization to solve the Schrödinger equation.
