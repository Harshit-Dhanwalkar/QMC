/*
xoshiro256** PRNG and derived distributions.
*/

#include "random.h"
#include <math.h>
#include <stdint.h>

static inline uint64_t rotl(const uint64_t val, int shift) {
  return (val << shift) | (val >> (64 - shift));
}

/*
 * SplitMix64, used only to expand single seed into 4 well-mixed xoshiro256**
 * state words.
 */
static uint64_t splitmix64_next(uint64_t *state) {
  uint64_t state_val = (*state += 0x9E3779B97F4A7C15ULL);
  state_val = (state_val ^ (state_val >> 30)) * 0xBF58476D1CE4E5B9ULL;
  state_val = (state_val ^ (state_val >> 27)) * 0x94D049BB133111EBULL;

  return state_val ^ (state_val >> 31);
}

void rng_seed(rng_state_t *rng, uint64_t seed) {
  if (!rng) {
    return;
  }

  uint64_t sm_state = seed;
  for (int i = 0; i < 4; i++) {
    rng->s[i] = splitmix64_next(&sm_state);
  }

  rng->has_cached_gaussian = 0;
  rng->cached_gaussian = 0.0;
}

uint64_t rng_next_u64(rng_state_t *rng) {
  const uint64_t result = rotl(rng->s[1] * 5, 7) * 9;
  const uint64_t temp_val = rng->s[1] << 17;

  rng->s[2] ^= rng->s[0];
  rng->s[3] ^= rng->s[1];
  rng->s[1] ^= rng->s[2];
  rng->s[0] ^= rng->s[3];
  rng->s[2] ^= temp_val;
  rng->s[3] = rotl(rng->s[3], 45);

  return result;
}

double rng_uniform(rng_state_t *rng) {
  // Top 53 bits -> exact double in [0,1), no bias
  return (double)(rng_next_u64(rng) >> 11) * (1.0 / 9007199254740992.0);
}

double rng_uniform_range(rng_state_t *rng, double min_val, double max_val) {
  return min_val + (max_val - min_val) * rng_uniform(rng);
}

double rng_gaussian(rng_state_t *rng) {
  if (rng->has_cached_gaussian) {
    rng->has_cached_gaussian = 0;

    return rng->cached_gaussian;
  }

  double rand1;
  double rand2;
  do {
    rand1 = rng_uniform(rng);
  } while (rand1 <= 1e-300); // avoid log(0)
  rand2 = rng_uniform(rng);

  double radius = sqrt(-2.0 * log(rand1));
  double theta = 2.0 * M_PI * rand2;

  rng->cached_gaussian = radius * sin(theta);
  rng->has_cached_gaussian = 1;

  return radius * cos(theta);
}

double rng_gaussian_scaled(rng_state_t *rng, double mean, double sigma) {
  return mean + sigma * rng_gaussian(rng);
}

void rng_jump(rng_state_t *rng) {
  static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
                                  0xa9582618e03fc9aa, 0x39abdc4529b1661c};

  uint64_t state0 = 0;
  uint64_t state1 = 0;
  uint64_t state2 = 0;
  uint64_t state3 = 0;
  for (int i = 0; i < 4; i++) {
    for (int bit_idx = 0; bit_idx < 64; bit_idx++) {
      if (JUMP[i] & (UINT64_C(1) << bit_idx)) {
        state0 ^= rng->s[0];
        state1 ^= rng->s[1];
        state2 ^= rng->s[2];
        state3 ^= rng->s[3];
      }

      rng_next_u64(rng);
    }
  }

  rng->s[0] = state0;
  rng->s[1] = state1;
  rng->s[2] = state2;
  rng->s[3] = state3;

  /* NOTE: A stream that gets jumped mid-sequence should not carry a Box-Muller
   * spare deviate computed from pre-jump stream into the post-jump stream. */
  rng->has_cached_gaussian = 0;
}

void rng_long_jump(rng_state_t *rng) {
  static const uint64_t LONG_JUMP[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3,
                                       0x77710069854ee241, 0x39109bb02acbe635};

  uint64_t state0 = 0;
  uint64_t state1 = 0;
  uint64_t state2 = 0;
  uint64_t state3 = 0;
  for (int i = 0; i < 4; i++) {
    for (int bit_idx = 0; bit_idx < 64; bit_idx++) {
      if (LONG_JUMP[i] & (UINT64_C(1) << bit_idx)) {
        state0 ^= rng->s[0];
        state1 ^= rng->s[1];
        state2 ^= rng->s[2];
        state3 ^= rng->s[3];
      }

      rng_next_u64(rng);
    }
  }

  rng->s[0] = state0;
  rng->s[1] = state1;
  rng->s[2] = state2;
  rng->s[3] = state3;

  rng->has_cached_gaussian = 0;
}
