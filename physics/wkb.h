#ifndef QMC_WKB_H
#define QMC_WKB_H

#include "potentials.h"

/* WKB quantization for bound states (Bohr-Sommerfeld):
   \int_{x1}^{x2} \sqrt(2m(E - V(x))) dx = (n + 1/2) \pi \hbar
   Returns energy for given n.
*/
double wkb_quantization(potential_fn V, void *params, double m, double hbar,
                        int n, double E_min, double E_max, double tol);

/* WKB transmission coefficient for a barrier:
   T = \exp(-2 \int_{x1}^{x2} \sqrt(2m(V(x) - E)) dx / \hbar)
*/
double wkb_transmission(potential_fn V, void *params, double m, double hbar,
                        double E, double x1, double x2, int N);

#endif
