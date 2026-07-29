#ifndef QMC_RANDOM_H
#define QMC_RANDOM_H

#include <stdint.h>

/*
 * xoshiro256** PRNG (Blackman & Vigna, 2018 https://arxiv.org/abs/1805.01407)
 * 256-bit state, period 2^256-1, passes BigCrush/PractRand
 * Not cryptographically secure for Monte Carlo sampling only
 */
typedef struct {
  uint64_t s[4];
} rng_state_t;

/* Seed from a single 64-bit value via SplitMix64 state expansion */
void rng_seed(rng_state_t *rng, uint64_t seed);

/* Raw 64-bit output */
uint64_t rng_next_u64(rng_state_t *rng);

/* Uniform double in [0, 1) */
double rng_uniform(rng_state_t *rng);

/* Uniform double in [a, b) */
double rng_uniform_range(rng_state_t *rng, double a, double b);

/* Standard normal N(0,1) via Box-Muller (one deviate per call; 2nd one
 * generated internally is cached and returned on next call) */
double rng_gaussian(rng_state_t *rng);

/* N(mean, \sigma^2) */
double rng_gaussian_scaled(rng_state_t *rng, double mean, double sigma);

#endif
