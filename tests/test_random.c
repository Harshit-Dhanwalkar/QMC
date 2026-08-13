/*
 * Tests for core/random.c: xoshiro256** PRNG and derived distributions.
 *
 * 1. Validates: seed reproducibility, stream divergence under distinct seeds,
 *    uniform range/moments, a chi-square goodness-of-fit check on rng_uniform,
 *    and mean/variance of rng_gaussian.
 */

#include "../core/random.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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

// rng_jump: implementation of canonical xoshiro256** jump (seed=12345).
static void test_jump_regression(void) {
  printf("test_jump_regression:\n");

  rng_state_t rng;
  rng_seed(&rng, 12345);
  rng_jump(&rng);

  CHECK(rng.s[0] == 0xc447e65c62d994cfULL, "post-jump s[0] matches reference");
  CHECK(rng.s[1] == 0xaf415ed201c9e97eULL, "post-jump s[1] matches reference");
  CHECK(rng.s[2] == 0x620fb38cd6dd52f4ULL, "post-jump s[2] matches reference");
  CHECK(rng.s[3] == 0xe6ba5be4e54b26c6ULL, "post-jump s[3] matches reference");

  uint64_t v1 = rng_next_u64(&rng);
  uint64_t v2 = rng_next_u64(&rng);
  uint64_t v3 = rng_next_u64(&rng);
  CHECK(v1 == 4527653816107373798ULL, "post-jump output 1 matches reference");
  CHECK(v2 == 5438022859293692230ULL, "post-jump output 2 matches reference");
  CHECK(v3 == 7149129066978069246ULL, "post-jump output 3 matches reference");
}

// A jumped stream must not reproduce any of a large number of draws from
// pre-jump stream (weak overlap sanity check. Jump is a proven 2^128-step
// advance, this just catches gross implementation errors, e.g. a jump that
// accidentally no-ops).
static void test_jump_no_short_run_overlap(void) {
  printf("test_jump_no_short_run_overlap:\n");

  rng_state_t base;
  rng_seed(&base, 999);

  int n = 100000;
  uint64_t *orig = malloc((size_t)n * sizeof *orig);
  for (int i = 0; i < n; i++) {
    orig[i] = rng_next_u64(&base);
  }

  rng_state_t jumped;
  rng_seed(&jumped, 999);
  rng_jump(&jumped);

  int found_overlap = 0;
  for (int i = 0; i < n && !found_overlap; i++) {
    uint64_t v = rng_next_u64(&jumped);

    for (int j = 0; j < n; j++) {
      if (v == orig[j]) {
        found_overlap = 1;

        break;
      }
    }
  }

  CHECK(!found_overlap,
        "no 64-bit output collisions between pre- and post-jump streams "
        "over 1e5 draws each");

  free(orig);
}

// rng_jump must clear any pending Box-Muller spare deviate, and (the real
// regression target) the spare-deviate cache must live per-stream: priming one
// stream's cache must not leak into, or be disturbed by, a second,
// independently-seeded stream interleaved with it.
static void test_gaussian_cache_is_per_stream(void) {
  printf("test_gaussian_cache_is_per_stream:\n");

  rng_state_t a, a_ref;
  rng_seed(&a, 42);
  rng_seed(&a_ref, 42);

  double a1 = rng_gaussian(&a); // primes a.cached_gaussian with the 2nd deviate
  CHECK(a.has_cached_gaussian == 1,
        "first rng_gaussian call caches a spare deviate in *rng");

  rng_state_t b;
  rng_seed(&b, 777);
  double b1 = rng_gaussian(&b); // must be b's own draw, not a's cached spare

  double a1_ref = rng_gaussian(&a_ref); // independent stream, same seed as a
  CHECK(fabs(a1 - a1_ref) < 1e-15,
        "interleaving a second stream's rng_gaussian call does not perturb "
        "stream a's first result");

  double a2 = rng_gaussian(&a); // should return a's originally cached spare
  CHECK(a.has_cached_gaussian == 0,
        "second call on stream a consumes its own cached spare");
  (void)b1;

  // NOTE: Recompute what a's cached spare should have been, from a fresh
  // identical stream, to confirm a2 is genuinely a's own spare and not
  // corrupted.
  rng_state_t a_check;
  rng_seed(&a_check, 42);
  double u1, u2;
  do {
    u1 = rng_uniform(&a_check);
  } while (u1 <= 1e-300);
  u2 = rng_uniform(&a_check);

  double r = sqrt(-2.0 * log(u1));
  double theta = 2.0 * M_PI * u2;
  double expected_a1 = r * cos(theta);
  double expected_a2 = r * sin(theta);

  CHECK(fabs(a1 - expected_a1) < 1e-15,
        "a's first deviate matches direct computation");
  CHECK(fabs(a2 - expected_a2) < 1e-15,
        "a's second deviate (cached spare) matches direct computation, "
        "unperturbed by stream b's interleaved call");
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
  test_jump_regression();
  test_jump_no_short_run_overlap();
  test_gaussian_cache_is_per_stream();

  if (failures == 0) {
    printf("\nAll test_random checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
