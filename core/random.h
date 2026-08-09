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
  int has_cached_gaussian; /* rng_gaussian's Box-Muller spare-deviate cache,
                            * stored per-stream (not per-thread/global) so that
                            * interleaving calls across multiple independent
                            * rng_state_t streams on same thread */
  double cached_gaussian;
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
 * generated internally is cached in *rng and returned on next call) */
double rng_gaussian(rng_state_t *rng);

/* N(mean, \sigma^2) */
double rng_gaussian_scaled(rng_state_t *rng, double mean, double sigma);

/*
 * Advance *rng by 2^128 steps in-place: canonical xoshiro256** jump (Reference
 * : Blackman & Vigna). Produces a stream guaranteed non-overlapping with
 * pre-jump stream for up to 2^128 draws from each. Generate independent
 * parallel RNG streams from one seed, used here so that e.g. thread/walker k's
 * stream is rng_jump()'d k times from a single master seed rather than seeded
 * independently (which relies only on SplitMix64's mixing quality, not a proven
 * non-overlap guarantee).
 */
void rng_jump(rng_state_t *rng);

/*
 * Advance *rng by 2^192 steps in-place: coarser xoshiro256** long-jump, for
 * generating independent *jump* starting points (e.g. one per MPI rank, each of
 * which then uses rng_jump() internally to fan out to its own threads/walkers).
 * 2^64 calls to rng_long_jump() partition the period into 2^64 subsequences of
 * length 2^192, each of which contains 2^64 non-overlapping rng_jump()
 * subsequences.
 */
void rng_long_jump(rng_state_t *rng);

#endif
