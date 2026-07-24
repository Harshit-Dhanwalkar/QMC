#ifndef QMC_SOFT_H
#define QMC_SOFT_H

#include "../core/vector.h"

/*
 * Split-Operator Fourier Transform (SOFT) method: 2D/3D generalization to
 * evolve TDSE split step, for an arbitrary (precomputed, grid-sampled)
 * potential V(x,y[,z]).
 *
 * Each step applies the Strang (symmetric) splitting of the propagator:
 *   \psi(t+dt) ~ \exp(-i V dt/(2 \hbar)) . \exp(-i hbar k^2 dt/(2m))
 *               . \exp(-i V dt/(2 \hbar)) . \psi(t)
 * (kinetic term applied in momentum space via fft2d/fft3d, potential terms
 * applied pointwise in position space).
 *
 * Grids are flat, row-major arrays: \psi[ix * Ny + iy] for 2D (length Nx * Ny),
 * \psi[(ix * Ny + iy) * Nz + iz] for 3D (length Nx * Ny * Nz)
 * V is a precomputed flat array of same length as \psi, V[grid index] =
 * V(x,y[,z]) at that grid point.
 */

int soft_evolve_2d(cvector_t *psi, const double *V, int Nx, int Ny, double dx,
                   double dy, double dt, int steps, double hbar, double mass);

int soft_evolve_3d(cvector_t *psi, const double *V, int Nx, int Ny, int Nz,
                   double dx, double dy, double dz, double dt, int steps,
                   double hbar, double mass);

/*
 * Norm (\int |psi|^2 over the grid, i.e. total probability) for flat 2D or 3D
 * grid wavefunction with grid spacing product cell_volume = dx * dy (2D) or dx
 * * dy * dz (3D).
 */
double grid_norm(const cvector_t *psi, double cell_volume);

#endif // QMC_SOFT_H
