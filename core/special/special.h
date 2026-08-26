#ifndef QMC_SPECIAL_H
#define QMC_SPECIAL_H

#include "../complex.h"
#include <math.h>

/* Hermite polynomials */
double hermite(int n, double x);
void hermite_array(int n, const double *x, int m, double *H);
double hermite_deriv(int n, double x);
double hermite_zeros(int n, int k);
int hermite_zeros_all(int n, double *zeros);

/* Numerically stable normalized Hermite function, via direct recurrence on
 * normalized function */
double hermite_function_stable(int n, double xi);

/* Laguerre polynomials */
double laguerre(int n, double alpha, double x);
void laguerre_array(int n, double alpha, const double *x, int N, double *L);

/* Legendre polynomials */
double legendre(int l, double x);
double assoc_legendre(int l, int m, double x);
void legendre_array(int l, const double *x, int N, double *P);

/* Spherical harmonics */
double spherical_harmonic_real(int l, int m, double theta, double phi);
double spherical_harmonic_imag(int l, int m, double theta, double phi);
complex_t spherical_harmonic(int l, int m, double theta, double phi);

/* Spherical Bessel functions */
double sph_bessel_j(int l, double x);
double sph_bessel_y(int l, double x);
double sph_bessel_j_deriv(int l, double x);
void sph_bessel_array(int lmax, double x, double *j, double *y);

/* Riccati-Bessel functions */
double riccati_bessel_j(int l, double x);
double riccati_bessel_y(int l, double x);
double riccati_bessel_j_deriv(int l, double x);

/* Utilities */
double factorial(int n);
double log_factorial(int n);

#endif
