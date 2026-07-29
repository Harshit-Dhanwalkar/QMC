/*
Tests for core/random.c: xoshiro256** PRNG and derived distributions.

Validates: seed reproducibility, stream divergence under distinct seeds,
uniform range/moments, a chi-square goodness-of-fit check on rng_uniform,
and mean/variance of rng_gaussian.
*/

#include "../core/random.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

static int failures = 0;

#define CHECK(cond, msg)                                                       \
  do {                                                                         \
    if (!(cond)) {                                                             \
      printf("  FAIL: %s\n", msg);                                             \
      failures++;                                                              \
    }                                                                          \
  } while (0)

// Reproducibility: same seed -> identical raw stream
static void test_seed_reproducibility(void) {
  printf("test_seed_reproducibility...\n");

  rng_state_t a, b;
  rng_seed(&a, 12345ULL);
  rng_seed(&b, 12345ULL);

  int identical = 1;
  for (int i = 0; i < 1000; i++) {
    uint64_t xa = rng_next_u64(&a);
    uint64_t xb = rng_next_u64(&b);
    if (xa != xb) {
      identical = 0;
      break;
    }
  }
  CHECK(identical, "identical seeds must produce identical u64 streams");
}

// Distinct seeds should (with overwhelming probability) diverge immediately
static void test_distinct_seeds_diverge(void) {
  printf("test_distinct_seeds_diverge...\n");

  rng_state_t a, b;
  rng_seed(&a, 1ULL);
  rng_seed(&b, 2ULL);

  int any_diff = 0;
  for (int i = 0; i < 16; i++) {
    if (rng_next_u64(&a) != rng_next_u64(&b)) {
      any_diff = 1;
      break;
    }
  }

  CHECK(any_diff, "distinct seeds must diverge within first few draws");
}

// State must actually evolve
static void test_state_evolves(void) {
  printf("test_state_evolves...\n");

  rng_state_t r;
  rng_seed(&r, 42ULL);

  uint64_t first = rng_next_u64(&r);
  int saw_different = 0;
  for (int i = 0; i < 100; i++) {
    if (rng_next_u64(&r) != first) {
      saw_different = 1;

      break;
    }
  }

  CHECK(saw_different, "successive draws must not be constant");
}

// rng_uniform: bounds, sample mean, sample variance
static void test_uniform_moments(void) {
  printf("test_uniform_moments...\n");

  rng_state_t r;
  rng_seed(&r, 777ULL);

  const int N = 2000000;
  double sum = 0.0, sum2 = 0.0;
  int in_bounds = 1;

  for (int i = 0; i < N; i++) {
    double u = rng_uniform(&r);
    if (u < 0.0 || u >= 1.0) {
      in_bounds = 0;
    }

    sum += u;
    sum2 += u * u;
  }

  double mean = sum / N;
  double var = sum2 / N - mean * mean;

  CHECK(in_bounds, "rng_uniform must stay within [0,1)");
  // Analytic: E[U] = 0.5, Var[U] = 1/12 ~= 0.08333
  CHECK(fabs(mean - 0.5) < 0.001, "rng_uniform sample mean near 0.5");
  CHECK(fabs(var - (1.0 / 12.0)) < 0.001,
        "rng_uniform sample variance near 1/12");
}

// rng_uniform_range: bounds and mean for non-trivial [a,b)
static void test_uniform_range(void) {
  printf("test_uniform_range...\n");

  rng_state_t r;
  rng_seed(&r, 555ULL);

  const double a = -3.0, b = 7.0;
  const int N = 1000000;
  double sum = 0.0;
  int in_bounds = 1;

  for (int i = 0; i < N; i++) {
    double x = rng_uniform_range(&r, a, b);
    if (x < a || x >= b) {
      in_bounds = 0;
    }

    sum += x;
  }

  double mean = sum / N;
  double expected_mean = 0.5 * (a + b);

  CHECK(in_bounds, "rng_uniform_range must stay within [a,b)");
  CHECK(fabs(mean - expected_mean) < 0.01,
        "rng_uniform_range sample mean near (a+b)/2");
}

// Chi-square goodness-of-fit on rng_uniform against uniform histogram
static void test_uniform_chi_square(void) {
  printf("test_uniform_chi_square...\n");

  rng_state_t r;
  rng_seed(&r, 2024ULL);

  const int K = 10;     // bins
  const int N = 200000; // samples
  int counts[10] = {0};

  for (int i = 0; i < N; i++) {
    double u = rng_uniform(&r);
    int bin = (int)(u * K);
    if (bin < 0) {
      bin = 0;
    }
    if (bin >= K) {
      bin = K - 1;
    }

    counts[bin]++;
  }

  double expected = (double)N / K;
  double chi2 = 0.0;
  for (int i = 0; i < K; i++) {
    double d = counts[i] - expected;
    chi2 += d * d / expected;
  }

  // dof = K-1 = 9. Critical value at \alpha=0.001 is ~27.9 (chi2 table).
  CHECK(chi2 < 30.0, "\\chi^2 statistic for uniform bins within bound");
  printf("  (\\chi^2 = %.3f, dof = %d)\n", chi2, K - 1);
}

// rng_gaussian: mean ~0, variance ~1
static void test_gaussian_moments(void) {
  printf("test_gaussian_moments...\n");

  rng_state_t r;
  rng_seed(&r, 9999ULL);

  const int N = 2000000;
  double sum = 0.0, sum2 = 0.0;

  for (int i = 0; i < N; i++) {
    double g = rng_gaussian(&r);
    sum += g;
    sum2 += g * g;
  }

  double mean = sum / N;
  double var = sum2 / N - mean * mean;

  CHECK(fabs(mean) < 0.003, "rng_gaussian sample mean near 0");
  CHECK(fabs(var - 1.0) < 0.005, "rng_gaussian sample variance near 1");
}

// rng_gaussian_scaled: mean/variance track (mean, \sigma^2)
static void test_gaussian_scaled(void) {
  printf("test_gaussian_scaled...\n");

  rng_state_t r;
  rng_seed(&r, 31415ULL);

  const double mu = 5.0, sigma = 2.0;
  const int N = 1000000;
  double sum = 0.0, sum2 = 0.0;

  for (int i = 0; i < N; i++) {
    double g = rng_gaussian_scaled(&r, mu, sigma);
    sum += g;
    sum2 += g * g;
  }

  double mean = sum / N;
  double var = sum2 / N - mean * mean;

  CHECK(fabs(mean - mu) < 0.01, "rng_gaussian_scaled sample mean near \\mu");
  CHECK(fabs(var - sigma * sigma) < 0.02,
        "rng_gaussian_scaled sample variance near \\sigma^2");
}

int main(void) {
  test_seed_reproducibility();
  test_distinct_seeds_diverge();
  test_state_evolves();
  test_uniform_moments();
  test_uniform_range();
  test_uniform_chi_square();
  test_gaussian_moments();
  test_gaussian_scaled();

  if (failures == 0) {
    printf("\nAll test_random checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
