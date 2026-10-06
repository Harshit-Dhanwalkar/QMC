/*
 * qmc_orbital.c - WebAssembly front end for hydrogen orbital viewer
 *
 * Physics is physics/orbital_sample.c (inverse-CDF radial sampling of
 * library's analytic R_nl, rejection-sampled spherical harmonics). This file
 * only owns point buffer JavaScript reads
 */

#include "../core/constants.h"
#include "../physics/hydrogen.h"
#include "../physics/orbital_sample.h"

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { ORB_MAX_POINTS = 80000 };

static double g_pts[4 * ORB_MAX_POINTS]; /* x, y, z, phase per point (a0, rad) */

QMC_EXPORT(qmc_orb_max_points)
int qmc_orb_max_points(void) { return ORB_MAX_POINTS; }

// Pointer to point buffer: 4 doubles per point (x, y, z, phase)
QMC_EXPORT(qmc_orb_buffer)
double *qmc_orb_buffer(void) { return g_pts; }

// Draw `count` points of orbital (n, l, m). Returns 0 on success
QMC_EXPORT(qmc_orb_sample)
int qmc_orb_sample(int n, int l, int m, int real_form, int count, int seed) {
  if (count < 1 || count > ORB_MAX_POINTS) {
    return -1;
  }

  return hydrogen_orbital_sample(n, l, m, real_form, count,
                                 (unsigned long long)(unsigned int)seed, g_pts);
}

// Bohr energy of shell n in eV (0 for invalid n)
QMC_EXPORT(qmc_orb_energy_ev)
double qmc_orb_energy_ev(int n) { return hydrogen_energy_level(n) / E_CHARGE; }

QMC_EXPORT(qmc_orb_max_n)
int qmc_orb_max_n(void) { return ORBITAL_MAX_N; }
