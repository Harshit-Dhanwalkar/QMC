#ifndef QMC_POTENTIALS_H
#define QMC_POTENTIALS_H

/* Potential function pointer */
typedef double (*potential_fn)(double x, void *params);

/* Specific potentials */

/* Infinite square well: V=0 for |x|<a, V=inf elsewhere */
double V_infinite_well(double x, void *params);
/* params: double *a (half-width) */

/* Finite square well: V=-V0 for |x|<a, V=0 elsewhere */
double V_finite_well(double x, void *params);
/* params: struct {double a, V0} */

/* Harmonic oscillator: V = (1/2)*m*omega^2*x^2 */
double V_harmonic(double x, void *params);
/* params: double *omega */

/* Step potential: V=0 for x<0, V=V0 for x>=0 */
double V_step(double x, void *params);
/* params: double *V0 */

/* Rectangular barrier: V=V0 for 0<x<a, V=0 elsewhere */
double V_barrier(double x, void *params);
/* params: struct {double a, V0} */

/* Coulomb potential: V = -k*e^2/r (for hydrogen) */
double V_coulomb(double r, void *params);

/* Yukawa potential: V = -g*exp(-mu*r)/r (screening) */
double V_yukawa(double r, void *params);

/* Morse potential: V = D*[1 - exp(-a(x-x0))]^2 (molecules) */
double V_morse(double x, void *params);

/* Generic: evaluate array of V values */
void potential_array(double *x, int n, potential_fn V, void *params,
                     double *V_out);

#endif
