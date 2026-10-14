/*
 * Quantum revivals of a wave packet in a box
 *
 * NOTE: A Gaussian packet bounces around an infinite square well
 * Because energies are E_n = n^2 whole state is periodic with T = 2 \pi: at T
 * it is exactly back, at T/2 it is its own mirror image, at T/4 it is a
 * superposition of packet and its mirror image (a cat), and in between it forms
 * intricate "quantum carpet" Return probability |<\psi_0|\psi(t)>|^2 shows
 * fractional revivals
 */

#include "../physics/box_revival.h"
#include "complex.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
  complex_t c[200];
  double kept = box_gaussian_coefficients(c, 200, 0.3, 0.04, 25.0);

  printf(" > Gaussian packet x0 = 0.3, \\sigma = 0.04, n0 = 25 (weight in 200 "
         "modes: %.6f)\n\n",
         kept);
  printf("  %8s  %14s\n", "t / T", "return prob.");

  for (int i = 0; i <= 16; i++) {
    double f = i / 16.0;

    printf("  %8.4f  %14.6f\n", f,
           c_abs2(box_overlap(c, 200, f * BOX_REVIVAL_TIME)));
  }

  return 0;
}
