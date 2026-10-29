/*
 * qmc_sg.c - WebAssembly front end for Stern-Gerlach demo
 *
 * Physics: physics/stern_gerlach.c (sequential spin-1/2 analyzers, partial
 * recombination, screen pattern). Shot generator below is a seeded xorshift
 * stream, so a run is reproducible and native and WASM builds give identical
 * hits
 */

#include "../physics/stern_gerlach.h"
#include <math.h>
#include <stdint.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { SG_MAX_ANALYZERS = 4, SG_MAX_POINTS = 1024, SG_MAX_SHOTS = 4000 };

static sg_analyzer_t g_an[SG_MAX_ANALYZERS];
static double g_inten[SG_MAX_ANALYZERS + 1];
static double g_r[(SG_MAX_ANALYZERS + 1) * 3];
static double g_z[SG_MAX_POINTS];
static double g_quantum[SG_MAX_POINTS];
static double g_classical[SG_MAX_POINTS];
static double g_hits[SG_MAX_SHOTS * 2];
static uint64_t g_state;

QMC_EXPORT(qmc_sg_max_analyzers)
int qmc_sg_max_analyzers(void) { return SG_MAX_ANALYZERS; }
QMC_EXPORT(qmc_sg_max_points)
int qmc_sg_max_points(void) { return SG_MAX_POINTS; }
QMC_EXPORT(qmc_sg_max_shots)
int qmc_sg_max_shots(void) { return SG_MAX_SHOTS; }

// Beam intensity after source and each analyzer (n + 1 values) and its
// polarisation vector (3 per entry)
QMC_EXPORT(qmc_sg_intensity_buffer)
double *qmc_sg_intensity_buffer(void) { return g_inten; }
QMC_EXPORT(qmc_sg_r_buffer)
double *qmc_sg_r_buffer(void) { return g_r; }
QMC_EXPORT(qmc_sg_z_buffer)
double *qmc_sg_z_buffer(void) { return g_z; }
QMC_EXPORT(qmc_sg_quantum_buffer)
double *qmc_sg_quantum_buffer(void) { return g_quantum; }
QMC_EXPORT(qmc_sg_classical_buffer)
double *qmc_sg_classical_buffer(void) { return g_classical; }
// Pairs (z, outcome) with outcome 1 for up beam and 0 for down
QMC_EXPORT(qmc_sg_hits_buffer)
double *qmc_sg_hits_buffer(void) { return g_hits; }

// Analyzer k in [0, 4): polar angle \theta, azimuth \phi, mode 0 pass up / 1
// pass down / 2 pass both, recombination coherence
QMC_EXPORT(qmc_sg_set)
int qmc_sg_set(int k, double theta, double phi, int mode, double coherence) {
  if (k < 0 || k >= SG_MAX_ANALYZERS || !isfinite(theta) || !isfinite(phi) ||
      mode < 0 || mode > 2 || !(coherence >= 0.0 && coherence <= 1.0)) {
    return -1;
  }

  g_an[k].theta = theta;
  g_an[k].phi = phi;
  g_an[k].mode = mode;
  g_an[k].coherence = coherence;

  return 0;
}

// Runs first n analyzers for a source of polarisation magnitude pol along
// (\theta, \phi)
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_sg_run)
int qmc_sg_run(double src_theta, double src_phi, double pol, int n) {
  if (n < 0 || n > SG_MAX_ANALYZERS || !(pol >= 0.0 && pol <= 1.0) ||
      !isfinite(src_theta) || !isfinite(src_phi)) {
    return -1;
  }

  sg_vec_t r0 = sg_axis(src_theta, src_phi), r[SG_MAX_ANALYZERS + 1];
  r0.x *= pol;
  r0.y *= pol;
  r0.z *= pol;
  if (sg_chain(r0, 1.0, g_an, n, g_inten, r) != 0) {
    return -1;
  }

  for (int k = 0; k <= n; k++) {
    g_r[3 * k] = r[k].x;
    g_r[3 * k + 1] = r[k].y;
    g_r[3 * k + 2] = r[k].z;
  }

  return 0;
}

// Screen pattern on npts points over [-zmax, zmax]
// Returns packet overlap \exp(-d^2 / 8 \sigma^2) or -1 on bad input
QMC_EXPORT(qmc_sg_screen)
double qmc_sg_screen(double pup, double d, double sigma, int npts,
                     double zmax) {
  if (!(pup >= 0.0 && pup <= 1.0 && d >= 0.0 && d < 1e3 && sigma > 1e-3 &&
        sigma < 1e3 && zmax > 0.0 && zmax < 1e4) ||
      npts < 2 || npts > SG_MAX_POINTS) {
    return -1.0;
  }

  for (int i = 0; i < npts; i++) {
    double z = -zmax + 2.0 * zmax * i / (npts - 1);

    g_z[i] = z;
    g_quantum[i] = sg_screen_density(z, pup, d, sigma);
    g_classical[i] = sg_classical_density(z, d, sigma);
  }

  return sg_coherence(d, sigma);
}

static double next_uniform(void) {
  g_state ^= g_state << 13;
  g_state ^= g_state >> 7;
  g_state ^= g_state << 17;

  return ((g_state >> 11) + 0.5) / 9007199254740992.0;
}

// Draws n atoms: outcome up with probability pup, then z from Gaussian of that
// beam
// Returns number of ups, or -1 on bad input
QMC_EXPORT(qmc_sg_shoot)
int qmc_sg_shoot(int n, double pup, double d, double sigma, int seed) {
  if (n < 0 || n > SG_MAX_SHOTS ||
      !(pup >= 0.0 && pup <= 1.0 && d >= 0.0 && sigma > 1e-3 && sigma < 1e3)) {
    return -1;
  }

  g_state = 0x9E3779B97F4A7C15ull ^
            ((uint64_t)(uint32_t)seed * 0xD1342543DE82EF95ull);
  if (!g_state) {
    g_state = 1;
  }

  int ups = 0;

  for (int i = 0; i < n; i++) {
    int up = next_uniform() < pup;
    double u1 = next_uniform(), u2 = next_uniform(),
           g = sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);

    ups += up;
    g_hits[2 * i] = (up ? 0.5 : -0.5) * d + sigma * g;
    g_hits[2 * i + 1] = up;
  }

  return ups;
}
