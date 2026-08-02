# SOFT: Split-Operator Fourier Transform (2D/3D)

The 2D/3D generalization of split-step time evolution, for an arbitrary (precomputed, grid-sampled) potential $V(x,y[,z])$. Builds on [`fft2d`/`fft3d`](../internals/linalg.md) for the kinetic term. Implemented in `physics/soft.h`.

## Method

Each step applies the Strang (symmetric) splitting of the propagator:

$$
\psi(t+dt) \approx e^{-iV\,dt/(2\hbar)}\cdot e^{-i\hbar k^2 dt/(2m)}\cdot e^{-iV\,dt/(2\hbar)}\cdot\psi(t)
$$

The potential half-steps are applied pointwise in position space, and the kinetic full-step is applied in momentum space via `fft2d`/`fft3d`, round-tripping through the FFT once per step. This is second-order accurate in `dt` due to the symmetric (half-step, full-step, half-step) splitting, compared to a naive first-order split.

## Grid convention

Grids are **flat, row-major** arrays, not multi-dimensional C arrays:

- 2D: `psi[ix * Ny + iy]`, length $N_x \times N_y$
- 3D: `psi[(ix * Ny + iy) * Nz + iz]`, length $N_x \times N_y \times N_z$

`V` is a precomputed flat array of the same length as `psi`, `V[grid index] = V(x,y[,z])` at that grid point - there's no potential function-pointer here (unlike [1D Potentials](potentials_1d.md)'s `potential_fn` convention); build the array before calling.

## Implementation

```c
int soft_evolve_2d(cvector_t *psi, const double *V, int Nx, int Ny, double dx,
                   double dy, double dt, int steps, double hbar, double mass);

int soft_evolve_3d(cvector_t *psi, const double *V, int Nx, int Ny, int Nz,
                   double dx, double dy, double dz, double dt, int steps,
                   double hbar, double mass);
```

```c
int Nx = 128, Ny = 128;
cvector_t *psi = // initial 2D wavepacket, flattened
double *V = malloc(Nx * Ny * sizeof(double));
// fill V[ix * Ny + iy] = V(x[ix], y[iy]) ...

soft_evolve_2d(psi, V, Nx, Ny, dx, dy, dt, steps, hbar, mass);
```

## Norm checking

```c
/* Total probability int|\psi|^2 over the grid, with cell_volume = dx*dy (2D) or dx*dy*dz (3D). */
double grid_norm(const cvector_t *psi, double cell_volume);
```

Since SOFT is unitary in principle, `grid_norm` staying close to 1 over the course of a long propagation is a cheap sanity check for numerical drift - particularly useful near grid boundaries where an absorbing layer (see [Crank-Nicolson's CAP support](../internals/crank_nicolson.md), if ported to this 2D/3D case) would otherwise be needed to prevent wraparound artifacts from the periodic boundary conditions implicit in using an FFT.

## Running the Example

```sh
./build/eg_22_soft
```
