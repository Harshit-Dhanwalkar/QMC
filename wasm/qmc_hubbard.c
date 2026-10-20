/*
 * qmc_hubbard.c - WebAssembly front end for the Hubbard-dimer demo.
 *
 * Physics: physics/hubbard_dimer.c (exact diagonalisation of two Hubbard sites). The hopping is the unit of energy
 * (t = 1, hbar = 1); the page varies U / t. Half filling (N = 2) is shown: the spectrum and ground-state observables
 * against U, the Neel quench, one-site probabilities and canonical thermal averages.
 */

#include "../physics/hubbard_dimer.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { HB_CURVE = 201, HB_TIME = 401, HB_NEEL = 9 }; /* Neel: site 1 up (bit 0) + site 2 down (bit 3) */

static double g_spec[6 * HB_CURVE], g_gs_d[HB_CURVE], g_gs_ss[HB_CURVE], g_gs_s[HB_CURVE];
static double g_q_d[HB_TIME], g_q_m[HB_TIME], g_q_s[HB_TIME], g_q_ss[HB_TIME];
static double g_th_d[HB_CURVE], g_th_ss[HB_CURVE], g_th_s[HB_CURVE], g_th_e[HB_CURVE];
static double g_probs[8], g_levels[6];

QMC_EXPORT(qmc_hb_curve_len)
int qmc_hb_curve_len(void) { return HB_CURVE; }

QMC_EXPORT(qmc_hb_time_len)
int qmc_hb_time_len(void) { return HB_TIME; }

/* Spectrum: level l (ascending at each U) at index l * curve_len + i, U = umax i / (curve_len - 1). */
QMC_EXPORT(qmc_hb_spec_buffer)
double *qmc_hb_spec_buffer(void) { return g_spec; }

/* Ground-state double occupancy (sum over both sites), <S1.S2> and one-site entropy against U. */
QMC_EXPORT(qmc_hb_gs_d_buffer)
double *qmc_hb_gs_d_buffer(void) { return g_gs_d; }
QMC_EXPORT(qmc_hb_gs_ss_buffer)
double *qmc_hb_gs_ss_buffer(void) { return g_gs_ss; }
QMC_EXPORT(qmc_hb_gs_s_buffer)
double *qmc_hb_gs_s_buffer(void) { return g_gs_s; }

/* Neel quench: doublons, 2<Sz1>, one-site entropy and <S1.S2> at t = tmax i / (time_len - 1). */
QMC_EXPORT(qmc_hb_q_d_buffer)
double *qmc_hb_q_d_buffer(void) { return g_q_d; }
QMC_EXPORT(qmc_hb_q_m_buffer)
double *qmc_hb_q_m_buffer(void) { return g_q_m; }
QMC_EXPORT(qmc_hb_q_s_buffer)
double *qmc_hb_q_s_buffer(void) { return g_q_s; }
QMC_EXPORT(qmc_hb_q_ss_buffer)
double *qmc_hb_q_ss_buffer(void) { return g_q_ss; }

/* Thermal curves at T = tmax (i + 1) / curve_len: doublons, <S1.S2>, entropy, energy. */
QMC_EXPORT(qmc_hb_th_d_buffer)
double *qmc_hb_th_d_buffer(void) { return g_th_d; }
QMC_EXPORT(qmc_hb_th_ss_buffer)
double *qmc_hb_th_ss_buffer(void) { return g_th_ss; }
QMC_EXPORT(qmc_hb_th_s_buffer)
double *qmc_hb_th_s_buffer(void) { return g_th_s; }
QMC_EXPORT(qmc_hb_th_e_buffer)
double *qmc_hb_th_e_buffer(void) { return g_th_e; }

/* Occupation probabilities (empty, up, down, double) of site 1 then site 2. */
QMC_EXPORT(qmc_hb_probs_buffer)
double *qmc_hb_probs_buffer(void) { return g_probs; }

/* The six half-filling energies at one U, ascending. */
QMC_EXPORT(qmc_hb_levels_buffer)
double *qmc_hb_levels_buffer(void) { return g_levels; }

static int u_ok(double u) { return u >= 0.0 && u <= 40.0; }

static void gs_state(const hd_solution_t *s, double *re, double *im) {
  int idx[6];

  hd_n2_states(s, idx);
  hd_eigenstate(s, idx[0], re, im);
}

/* Curves over U in [0, umax]. Returns 0, or -1 on bad input. */
QMC_EXPORT(qmc_hb_curves)
int qmc_hb_curves(double umax) {
  if (!(umax > 0.0 && umax <= 40.0)) {
    return -1;
  }
  for (int i = 0; i < HB_CURVE; i++) {
    hd_solution_t s;
    int idx[6];
    double re[HD_DIM], im[HD_DIM];

    hd_solve(1.0, umax * i / (HB_CURVE - 1), &s);
    hd_n2_states(&s, idx);
    for (int l = 0; l < 6; l++) {
      g_spec[l * HB_CURVE + i] = s.e[idx[l]];
    }
    hd_eigenstate(&s, idx[0], re, im);
    g_gs_d[i] = hd_expect_double_occ(re, im);
    g_gs_ss[i] = hd_expect_ss(re, im);
    g_gs_s[i] = hd_site_entropy(re, im, 0);
  }

  return 0;
}

/* Levels and ground energy at one U; returns the ground energy, NaN on bad input. */
QMC_EXPORT(qmc_hb_levels)
double qmc_hb_levels(double u) {
  hd_solution_t s;
  int idx[6];

  if (!u_ok(u)) {
    return NAN;
  }
  hd_solve(1.0, u, &s);
  hd_n2_states(&s, idx);
  for (int l = 0; l < 6; l++) {
    g_levels[l] = s.e[idx[l]];
  }

  return g_levels[0];
}

/* Neel quench at interaction u over [0, tmax]. Returns 0, or -1 on bad input. */
QMC_EXPORT(qmc_hb_quench)
int qmc_hb_quench(double u, double tmax) {
  hd_solution_t s;
  double re0[HD_DIM], im0[HD_DIM], re[HD_DIM], im[HD_DIM];

  if (!u_ok(u) || !(tmax > 0.0 && tmax <= 500.0)) {
    return -1;
  }
  hd_solve(1.0, u, &s);
  hd_basis_state(HB_NEEL, re0, im0);
  for (int i = 0; i < HB_TIME; i++) {
    hd_evolve(&s, re0, im0, tmax * i / (HB_TIME - 1), re, im);
    g_q_d[i] = hd_expect_double_occ(re, im);
    g_q_m[i] = 2.0 * hd_expect_sz(re, im, 0);
    g_q_s[i] = hd_site_entropy(re, im, 0);
    g_q_ss[i] = hd_expect_ss(re, im);
  }

  return 0;
}

/* Site occupation probabilities. mode 0: ground state, mode 1: Neel state evolved for `time`. Returns 0, or -1 on bad input. */
QMC_EXPORT(qmc_hb_probs)
int qmc_hb_probs(double u, int mode, double time) {
  hd_solution_t s;
  double re[HD_DIM], im[HD_DIM], re0[HD_DIM], im0[HD_DIM];

  if (!u_ok(u) || (mode != 0 && mode != 1) || !isfinite(time)) {
    return -1;
  }
  hd_solve(1.0, u, &s);
  if (mode == 0) {
    gs_state(&s, re, im);
  } else {
    hd_basis_state(HB_NEEL, re0, im0);
    hd_evolve(&s, re0, im0, time, re, im);
  }
  hd_site_probs(re, im, 0, g_probs);
  hd_site_probs(re, im, 1, g_probs + 4);

  return 0;
}

/* Thermal curves at interaction u over T in (0, tmax]. Returns 0, or -1 on bad input. */
QMC_EXPORT(qmc_hb_thermal)
int qmc_hb_thermal(double u, double tmax) {
  hd_solution_t s;
  hd_thermal_t th;

  if (!u_ok(u) || !(tmax > 0.0 && tmax <= 100.0)) {
    return -1;
  }
  hd_solve(1.0, u, &s);
  for (int i = 0; i < HB_CURVE; i++) {
    hd_thermal_n2(&s, tmax * (i + 1) / HB_CURVE, &th);
    g_th_d[i] = th.double_occ;
    g_th_ss[i] = th.ss;
    g_th_s[i] = th.entropy;
    g_th_e[i] = th.energy;
  }

  return 0;
}
