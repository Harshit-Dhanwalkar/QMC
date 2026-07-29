/*
xoshiro256** PRNG and derived distributions.
*/

#include "random.h"
#include <math.h>
#include <stdint.h>

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

/* SplitMix64, used only to expand single seed into 4 well-mixed
 * xoshiro256** state words. */
static uint64_t splitmix64_next(uint64_t *state) {
  uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;

  return z ^ (z >> 31);
}

void rng_seed(rng_state_t *rng, uint64_t seed) {
  if (!rng) {
    return;
  }

  uint64_t sm_state = seed;
  for (int i = 0; i < 4; i++) {
    rng->s[i] = splitmix64_next(&sm_state);
  }
}

uint64_t rng_next_u64(rng_state_t *rng) {
  const uint64_t result = rotl(rng->s[1] * 5, 7) * 9;
  const uint64_t t = rng->s[1] << 17;

  rng->s[2] ^= rng->s[0];
  rng->s[3] ^= rng->s[1];
  rng->s[1] ^= rng->s[2];
  rng->s[0] ^= rng->s[3];
  rng->s[2] ^= t;
  rng->s[3] = rotl(rng->s[3], 45);

  return result;
}

double rng_uniform(rng_state_t *rng) {
  // Top 53 bits -> exact double in [0,1), no bias
  return (double)(rng_next_u64(rng) >> 11) * (1.0 / 9007199254740992.0);
}

double rng_uniform_range(rng_state_t *rng, double a, double b) {
  return a + (b - a) * rng_uniform(rng);
}

double rng_gaussian(rng_state_t *rng) {
  static _Thread_local int have_cached = 0;
  static _Thread_local double cached = 0.0;

  if (have_cached) {
    have_cached = 0;

    return cached;
  }

  double u1, u2;
  do {
    u1 = rng_uniform(rng);
  } while (u1 <= 1e-300); // avoid log(0)
  u2 = rng_uniform(rng);

  double r = sqrt(-2.0 * log(u1));
  double theta = 2.0 * M_PI * u2;

  cached = r * sin(theta);
  have_cached = 1;

  return r * cos(theta);
}

double rng_gaussian_scaled(rng_state_t *rng, double mean, double sigma) {
  return mean + sigma * rng_gaussian(rng);
}
