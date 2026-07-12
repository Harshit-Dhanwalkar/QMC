#ifndef QMC_WKB_H
#define QMC_WKB_H

#include "potentials.h"

/* WKB quantization for bound states (Bohr-Sommerfeld):
   \int_{x1}^{x2} \sqrt(2m(E - V(x))) dx = (n + 1/2) \pi \hbar
   Returns energy for given n.
*/
double wkb_quantization(potential_fn V, void *params, double m, double hbar,
                        int n, double E_min, double E_max, double x_min,
                        double x_max, int N, double tol);

/* WKB transmission coefficient for a barrier:
   T = \exp(-2 \int_{x1}^{x2} \sqrt(2m(V(x) - E)) dx / \hbar)
*/
double wkb_transmission(potential_fn V, void *params, double m, double hbar,
                        double E, double x1, double x2, int N);

/* Wrapper: WKB energy for level n of the harmonic oscillator
   V(x) = 0.5*omega^2*x^2
*/
double wkb_energy_harmonic(int n, double hbar, double omega);

/* Wrapper: full round-trip Bohr-Sommerfeld action integral
   \oint p dx = 2 * \int_{-x_turn}^{x_turn} * \sqrt(2m(E - V(x))) dx
   for the harmonic oscillator at given energy E (m=1).
   Quantized levels satisfy action (2 * \pi * \hbar) = n + 1/2.
*/
double wkb_action_integral_harmonic(double E, double hbar, double omega);

/* Wrapper: WKB tunneling probability through rectangular barrier of
   height V0 and width, closed-form :
   T = \exp(-2 * \kappa * width), \kappa = \sqrt((V0-E) / hbar_sq_2m)
   where hbar_sq_2m = \hbar^2/(2m),
   Returns 1.0 (fully transmitting) if V0 <= E
*/
double wkb_tunneling_rectangular(double E, double V0, double width,
                                 double hbar_sq_2m);

#endif
