/*
Test + demonstration: driven two-level systems (general time-dependent H(t),
Landau-Zener sweeps, and lab-frame driving beyond RWA).

1. Constant \Delta, \Omega cross-check: driven_two_level_evolve with
  time_fn_constant for both must reproduce rabi_evolve_exact()'s closed-form
  result, validates the RK4 integrator and rotating-frame RHS independent of any
  time-dependence-specific code.
2. Landau-Zener: sweep \Delta(t) = \alpha * t from deep negative to deep
  positive time, starting in the diabatic ground state; compare final diabatic
  population to the closed-form landau_zener_probability(), in both
  fast/diabatic (large alpha) and slow/adiabatic (small alpha) regimes.
3. Lab-frame vs RWA (Bloch-Siegert): in weak-driving regime (\Omega_0 <<
  \omega_0), full non-RWA lab-frame simulation must agree with
  rabi_excited_probability() evaluated at RWA-equivalent parameters (\Delta =
  \omega_L - \omega_0, \Omega = \Omega_0/2) to within a small tolerance.
  Bloch-Siegert-type correction from counter-rotating term RWA drops.
*/

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/driven.h"
#include "../physics/rabi.h"
#include <math.h>
#include <stdio.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);
  return err > tol;
}

// Test 1: constant Delta/Omega must match rabi_evolve_exact
static int test_constant_matches_rabi_exact(void) {
  double Omega = 2.0, Delta = 1.0;
  double T = 3.0, dt = 1e-4;
  int steps = (int)(T / dt);

  cvector_t *psi_driven = cvector_alloc(2);
  psi_driven->data[0] = c_real(1.0);
  cvector_t *psi_exact = cvector_alloc(2);
  psi_exact->data[0] = c_real(1.0);

  int fail =
      driven_two_level_evolve(psi_driven, time_fn_constant, &Delta,
                              time_fn_constant, &Omega, 0.0, dt, steps) != 0;
  fail |= rabi_evolve_exact(psi_exact, T, Omega, Delta) != 0;

  fail |= check_close(c_abs2(psi_driven->data[1]), c_abs2(psi_exact->data[1]),
                      1e-6, "P_excited: driven(constant) vs rabi_evolve_exact");
  fail |= check_close(psi_driven->data[0].re, psi_exact->data[0].re, 1e-5,
                      "Re(amp ground): driven(constant) vs exact");
  fail |= check_close(psi_driven->data[1].im, psi_exact->data[1].im, 1e-5,
                      "Im(amp excited): driven(constant) vs exact");

  cvector_free(psi_driven);
  cvector_free(psi_exact);

  return fail;
}

// Test 2: Landau-Zener, fast (diabatic) and slow (adiabatic) regimes.
// NOTE: T and dt differ between the two regimes on purpose.
// The crossing integrated out to |\Delta(T)| = \alpha * T >> \Omega for
// t->\pm\inf asymptotic formula to apply, but a large alpha reaches that regime
// (and starts oscillating at frequency ~\alpha * T) far sooner than a small
// alpha does. Using one shared (T, dt) for both regimes either under-resolves
// fast case or wastes enormous step counts on slow one. The residual error at
// finite T is a physical effect (Stuckelberg oscillations : finite-time
// interference b/w diabatic and adiabatic paths)
static int test_landau_zener(void) {
  int fail = 0;
  double Omega = 1.0;

  // Fast sweep: large alpha -> near-diabatic passage -> P(stay in state0) -> 1
  {
    double alpha = 20.0;
    double T = 10.0,
           dt = 1e-4; // \Delta(T)=200: fine dt (fast local oscillation)
    int steps = (int)(2.0 * T / dt);
    cvector_t *psi = cvector_alloc(2);
    psi->data[0] = c_real(1.0);
    driven_two_level_evolve(psi, time_fn_linear_ramp, &alpha, time_fn_constant,
                            &Omega, -T, dt, steps);
    double p_diabatic = c_abs2(psi->data[0]);
    double p_lz = landau_zener_probability(Omega, alpha);
    fail |= check_close(p_diabatic, p_lz, 0.02,
                        "Landau-Zener fast sweep (alpha=20, near-diabatic)");
    cvector_free(psi);
  }

  // Slow sweep: small alpha -> near-adiabatic passage -> P(stay in state0) -> 0
  {
    double alpha = 0.3;
    double T = 40.0,
           dt = 2e-3; // \Delta(T)=12: large T (slow to reach asymptote)
    int steps = (int)(2.0 * T / dt);
    cvector_t *psi = cvector_alloc(2);
    psi->data[0] = c_real(1.0);
    driven_two_level_evolve(psi, time_fn_linear_ramp, &alpha, time_fn_constant,
                            &Omega, -T, dt, steps);
    double p_diabatic = c_abs2(psi->data[0]);
    double p_lz = landau_zener_probability(Omega, alpha);
    fail |= check_close(p_diabatic, p_lz, 0.02,
                        "Landau-Zener slow sweep (alpha=0.3, near-adiabatic)");
    cvector_free(psi);
  }

  return fail;
}

// Test 3: weak-driving lab-frame simulation vs RWA prediction
static int test_lab_frame_weak_driving_matches_rwa(void) {
  double omega0 = 50.0, Omega0 = 1.0, omega_L = 50.0; // resonant, weak drive
  double T = 3.0, dt = 1e-3;
  int steps = (int)(T / dt);

  cvector_t *psi = cvector_alloc(2);
  psi->data[0] = c_real(1.0);
  driven_two_level_evolve_lab_frame(psi, omega0, Omega0, omega_L, 0.0, 0.0, dt,
                                    steps);
  double p_lab_frame = c_abs2(psi->data[1]);

  // NOTE: with convention : H_lab = (\omega_0/2) * \sigma_z + \Omega_0 *
  // cos(...) * \sigma_x, H_RWA = (\Delta/2) * sigma_z + (\Omega/2) * sigma_x -
  // RWA reduction maps \Omega_rwa = \Omega_0 (no extra factor of 2): the "1/2"
  // from dropping counter-rotating term is already same "1/2" built into
  // rabi.h's own (Omega/2) convention.
  double Delta_rwa = omega_L - omega0;
  double Omega_rwa = Omega0;
  double p_rwa = rabi_excited_probability(T, Omega_rwa, Delta_rwa);

  int fail = check_close(p_lab_frame, p_rwa, 0.01,
                         "P_excited: lab-frame (no RWA) vs RWA prediction, "
                         "weak driving");

  cvector_free(psi);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Constant Delta/Omega matches rabi_evolve_exact:\n");
  failed += test_constant_matches_rabi_exact();

  printf("Landau-Zener sweep (fast=diabatic, slow=adiabatic):\n");
  failed += test_landau_zener();

  printf("Lab-frame (beyond RWA) vs RWA prediction, weak driving:\n");
  failed += test_lab_frame_weak_driving_matches_rwa();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
