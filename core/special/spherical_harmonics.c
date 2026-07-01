/*
Spherical harmonics
*/
#include "../complex.h"
#include "special.h"
#include <math.h>

/* Real part of spherical harmonic Y_l^m(\theta, \phi) */
double spherical_harmonic_real(int l, int m, double theta, double phi) {
  if (m < 0) {
    // Use Y_l^{-m} = (-1)^m conj(Y_l^m)
    // Real part = (-1)^m * real(Y_l^m) for m>0? Actually formula: Y_l^{-m} =
    // (-1)^m (Y_l^m)^* So real = (-1)^m * real(Y_l^m)
    double val = spherical_harmonic_real(l, -m, theta, phi);
    if (m % 2)
      return -val;
    return val;
  }
  double norm = sqrt((2.0 * l + 1.0) / (4.0 * M_PI) *
                     exp(log_factorial(l - m) - log_factorial(l + m)));
  double P = assoc_legendre(l, m, cos(theta));
  return norm * P * cos(m * phi);
}

double spherical_harmonic_imag(int l, int m, double theta, double phi) {
  if (m < 0) {
    // Imag part = -(-1)^m * imag(Y_l^m)
    double val = spherical_harmonic_imag(l, -m, theta, phi);
    if (m % 2)
      return val; // since -(-1)^m = (-1)^{m+1}? Let's derive: Y_l^{-m} = (-1)^m
                  // conj(Y_l^m), so imag = (-1)^m * (-imag(Y_l^m)) = -(-1)^m
                  // imag(Y_l^m)
    return -val;
  }
  double norm = sqrt((2.0 * l + 1.0) / (4.0 * M_PI) *
                     exp(log_factorial(l - m) - log_factorial(l + m)));
  double P = assoc_legendre(l, m, cos(theta));
  return norm * P * sin(m * phi);
}

complex_t spherical_harmonic(int l, int m, double theta, double phi) {
  return c_new(spherical_harmonic_real(l, m, theta, phi),
               spherical_harmonic_imag(l, m, theta, phi));
}
