/*
 * Test: Boson sampling
 *
 * 1. Hong-Ou-Mandel dip: two photons through a 50:50 beam splitter must give
 *    P(1,1)=0 exactly (destructive interference/bunching) and
 *    P(2,0)=P(0,2)=0.5, summing to 1.
 * 2. Unitarity: for an arbitrary unitary network (DFT_3), summing probabilities
 *    over all output configurations for a fixed input must equal 1.
 * 3. Permutation invariance: probability must not depend on order photons are
 *    listed in input_modes (only on multiset).
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../physics/boson_sampling.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static int test_hom_dip(void) {
  cmatrix_t *U = beam_splitter_50_50();
  const int in[2] = {0, 1};
  const int out_20[2] = {0, 0};
  const int out_02[2] = {1, 1};
  const int out_11[2] = {0, 1};

  double p20 = boson_sampling_probability(U, 2, in, out_20, 2);
  double p02 = boson_sampling_probability(U, 2, in, out_02, 2);
  double p11 = boson_sampling_probability(U, 2, in, out_11, 2);

  cmatrix_free(U);

  int fail = check_close(p20, 0.5, 1e-10, "P(2,0)");
  fail |= check_close(p02, 0.5, 1e-10, "P(0,2)");
  fail |= check_close(p11, 0.0, 1e-10, "P(1,1) HOM dip");
  fail |= check_close(p20 + p02 + p11, 1.0, 1e-10, "sum");

  return fail;
}

static int test_unitarity_dft(void) {
  cmatrix_t *D = dft_unitary(3);
  const int in[2] = {0, 1};
  double total = 0.0;

  for (int a = 0; a < 3; a++) {
    for (int b = a; b < 3; b++) {
      const int out[2] = {a, b};

      total += boson_sampling_probability(D, 3, in, out, 2);
    }
  }

  cmatrix_free(D);

  return check_close(total, 1.0, 1e-10, "sum over all outputs (DFT_3)");
}

static int test_permutation_invariance(void) {
  cmatrix_t *D = dft_unitary(3);
  const int in_a[2] = {0, 1};
  const int in_b[2] = {1, 0}; // same input multiset, different order
  const int out[2] = {2, 2};

  double pa = boson_sampling_probability(D, 3, in_a, out, 2);
  double pb = boson_sampling_probability(D, 3, in_b, out, 2);

  cmatrix_free(D);

  return check_close(pa, pb, 1e-12, "P invariant under input photon order");
}

int main(void) {
  int failed = 0;

  printf("Hong-Ou-Mandel dip (2 photons, 50:50 beam splitter):\n");
  failed += test_hom_dip();

  printf("Unitarity check (DFT_3 network, sum over all outputs):\n");
  failed += test_unitarity_dft();

  printf("Permutation invariance:\n");
  failed += test_permutation_invariance();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
