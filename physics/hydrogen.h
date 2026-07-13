#ifndef QMC_HYDROGEN_H
#define QMC_HYDROGEN_H

#include "../core/matrix.h"
#include "../core/vector.h"
#include "wavefn.h"

/* Radial Schrodinger equation for hydrogen:
   [ -\hbar^2 / (2m) * (d^2 / dr^2 - l(l+1) / r^2) - e^2 / (4 * \pi * \epsilon_0
   r) ] R(r) = E R(r)
   Returns eigen_t with eigenvalues (energies) and eigenvectors (radial
   functions).
   NOTE: The radial grid r[0..N-1] must be given (typically logarithmic). l is
   the angular momentum quantum number.
*/
eigen_t *hydrogen_radial_solve(double *r, int N, int l, double hbar,
                               double mass, double e_charge, double eps0);

/* Analytic hydrogen energy levels (SI units) */
double hydrogen_energy_level(int n);

/* Hydrogen radial wavefunction (analytic, unnormalized) for given n, l.
   Returns new cvector_t* with values on the grid r.
   r: grid, N: size, n: principal quantum number, l: angular momentum.
   NOTE: The function uses associated Laguerre polynomials.
*/
cvector_t *hydrogen_radial_wavefunction(double *r, int N, int n, int l);

#endif
