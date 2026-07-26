/*
 * Boson Sampling: Hong-Ou-Mandel Interference and a 3-Mode Example
 *
 * Two-photon Hong-Ou-Mandel (HOM) dip through a 50:50 beam splitter:
 * indistinguishable bosons bunch together at the output, with two "one photon
 * in each output" amplitude paths interfering destructively.
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../physics/boson_sampling.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Boson Sampling: Hong-Ou-Mandel Interference\n\n");

  cmatrix_t *U = beam_splitter_50_50();
  int input_modes[2] = {0, 1};

  printf("   Two photons into a 50:50 beam splitter (modes 0,1):\n");
  printf("   output      P\n");
  printf("   --------   --------\n");

  int out_00[2] = {0, 0};
  int out_11[2] = {1, 1};
  int out_01[2] = {0, 1};

  double p00 = boson_sampling_probability(U, 2, input_modes, out_00, 2);
  double p11 = boson_sampling_probability(U, 2, input_modes, out_11, 2);
  double p01 = boson_sampling_probability(U, 2, input_modes, out_01, 2);

  printf("   (2,0)      %.6f  (expect 0.5, both photons bunch at mode 0)\n",
         p00);
  printf("   (0,2)      %.6f  (expect 0.5, both photons bunch at mode 1)\n",
         p11);
  printf(
      "   (1,1)      %.6f  (expect 0.0, HOM dip: destructive interference)\n",
      p01);
  printf("   sum = %.6f (expect 1.0)\n\n", p00 + p11 + p01);

  cmatrix_free(U);

  printf("   3-mode DFT network, 2 photons in modes (0,1):\n");
  printf("   output      P\n");
  printf("   --------   --------\n");

  cmatrix_t *D = dft_unitary(3);
  int dft_input[2] = {0, 1};
  double total = 0.0;

  for (int a = 0; a < 3; a++) {
    for (int b = a; b < 3; b++) {
      int out[2] = {a, b};
      double p = boson_sampling_probability(D, 3, dft_input, out, 2);
      printf("   (%d,%d)      %.6f\n", a, b, p);
      total += p;
    }
  }
  printf("   sum = %.6f (expect 1.0 - unitarity check)\n", total);

  cmatrix_free(D);

  return 0;
}
