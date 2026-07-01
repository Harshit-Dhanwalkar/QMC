#ifndef QMC_CRANK_NICOLSON_H
#define QMC_CRANK_NICOLSON_H

#include "../matrix.h"
#include "../vector.h"

/* Crank-Nicolson time evolution for 1D Schrödinger equation.
   Hamiltonian must be tridiagonal (kinetic + diagonal potential).
   Solves: (I + i*dt/2*H) ψ_{n+1} = (I - i*dt/2*H) ψ_n
   using a complex tridiagonal Thomas algorithm.
*/

/* Evolve wavefunction psi for one time step dt.
   diag: diagonal of H (real, size N)
   offdiag: off-diagonal of H (real, size N-1)
   psi: input/output complex vector
   Returns 0 on success, -1 on failure.
*/
int crank_nicolson_step(const double *diag, const double *offdiag, double dt,
                        cvector_t *psi);

/* Build tridiagonal Hamiltonian for 1D system:
   H = -hbar^2/(2m) * d^2/dx^2 + V(x)
   discretized: diag[i] = 2*coeff + V[i], offdiag[i] = -coeff,
   where coeff = hbar^2/(2m) / dx^2.
*/
void build_tridiagonal_hamiltonian(const double *x, const double *V, int N,
                                   double dx, double hbar_sq_2m, double *diag,
                                   double *offdiag);

#endif
