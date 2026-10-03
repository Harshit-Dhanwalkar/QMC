#include "td_perturbation.h"
#include "complex.h"
#include "matrix.h"
#include "vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Result management
void tdpt_result_free(tdpt_result_t *result) {
  if (!result) {
    return;
  }

  if (result->orders) {
    for (int k = 0; k <= result->max_order; k++) {
      free(result->orders[k].data);
    }
  }
  free(result->orders);
  free(result->energies);
  free(result);
}

static tdpt_result_t *tdpt_result_alloc(int n, int max_order, double t0,
                                        double t1, const double *energies) {
  tdpt_result_t *res = calloc(1, sizeof *res);
  if (!res) {
    return NULL;
  }

  res->n = n;
  res->max_order = max_order;
  res->t0 = t0;
  res->t1 = t1;
  res->energies = malloc((size_t)n * sizeof *res->energies);
  res->orders = calloc((size_t)max_order + 1, sizeof *res->orders);
  if (!res->energies || !res->orders) {
    tdpt_result_free(res);

    return NULL;
  }

  memcpy(res->energies, energies, (size_t)n * sizeof *energies);
  for (int k = 0; k <= max_order; k++) {
    res->orders[k].n = n;
    res->orders[k].data = calloc((size_t)n, sizeof *res->orders[k].data);
    if (!res->orders[k].data) {
      tdpt_result_free(res);

      return NULL;
    }
  }

  return res;
}

// Dyson hierarchy integration
// Right-hand side of the hierarchy at time t:
//   dy_0 = 0,   dy_k = -i P(t) V(t) P(t)^dagger y_{k-1}   (k >= 1),
//
// Where
//   P(t) = diag(e^{i E_n t})
// y is laid out as y[k * n + i]. `ph` and `z` are length-n scratch arrays.
static void tdpt_rhs(int n, int K, const double *E, double t,
                     const cmatrix_t *V, const complex_t *y, complex_t *dy,
                     complex_t *ph, complex_t *z) {
  for (int i = 0; i < n; i++) {
    ph[i] = (complex_t){cos(E[i] * t), sin(E[i] * t)};
    dy[i] = c_zero();
  }

  for (int k = 1; k <= K; k++) {
    const complex_t *prev = y + (size_t)(k - 1) * n;
    complex_t *out = dy + (size_t)k * n;
    for (int m = 0; m < n; m++) {
      z[m] = c_mul(c_conj(ph[m]), prev[m]);
    }

    for (int f = 0; f < n; f++) {
      double re = 0.0;
      double im = 0.0;

      for (int m = 0; m < n; m++) {
        complex_t v = V->data[(size_t)f * n + m];
        re += v.re * z[m].re - v.im * z[m].im;
        im += v.re * z[m].im + v.im * z[m].re;
      }

      complex_t pw = c_mul(ph[f], (complex_t){re, im});
      out[f] = (complex_t){pw.im, -pw.re}; // multiply by -i
    }
  }
}

static void tdpt_axpy(size_t len, const complex_t *y, double scale,
                      const complex_t *k, complex_t *dst) {
  for (size_t i = 0; i < len; i++) {
    dst[i].re = y[i].re + scale * k[i].re;
    dst[i].im = y[i].im + scale * k[i].im;
  }
}

static int tdpt_all_finite(const complex_t *a, size_t len) {
  for (size_t i = 0; i < len; i++) {
    if (!isfinite(a[i].re) || !isfinite(a[i].im)) {
      return 0;
    }
  }

  return 1;
}

// Integrate hierarchy y over [t0, t1]; y holds initial condition on entry and
// final state on exit. Returns 0 on success, -1 on allocation failure
static int tdpt_integrate(int n, int K, const double *E,
                          tdpt_perturbation_fn Vfn, void *params, double t0,
                          double t1, int steps, complex_t *y) {
  size_t len = (size_t)(K + 1) * (size_t)n;
  size_t nn = (size_t)n * (size_t)n;

  complex_t *work = calloc(5 * len + 2 * (size_t)n, sizeof *work);
  cmatrix_t *Va = cmatrix_alloc(n, n);
  cmatrix_t *Vm = cmatrix_alloc(n, n);
  cmatrix_t *Vb = cmatrix_alloc(n, n);

  if (!work || !Va || !Vm || !Vb) {
    free(work);
    cmatrix_free(Va);
    cmatrix_free(Vm);
    cmatrix_free(Vb);

    return -1;
  }

  complex_t *k1 = work;
  complex_t *k2 = work + len;
  complex_t *k3 = work + 2 * len;
  complex_t *k4 = work + 3 * len;
  complex_t *tmp = work + 4 * len;
  complex_t *ph = work + 5 * len;
  complex_t *z = ph + n;

  double h = (t1 - t0) / (double)steps;
  memset(Va->data, 0, nn * sizeof *Va->data);

  Vfn(t0, params, Va);
  for (int s = 0; s < steps; s++) {
    double t = t0 + (double)s * h;
    double t_end = (s == steps - 1) ? t1 : t0 + (double)(s + 1) * h;

    memset(Vm->data, 0, nn * sizeof *Vm->data);
    memset(Vb->data, 0, nn * sizeof *Vb->data);

    Vfn(t + 0.5 * h, params, Vm);
    Vfn(t_end, params, Vb);

    tdpt_rhs(n, K, E, t, Va, y, k1, ph, z);
    tdpt_axpy(len, y, 0.5 * h, k1, tmp);
    tdpt_rhs(n, K, E, t + 0.5 * h, Vm, tmp, k2, ph, z);
    tdpt_axpy(len, y, 0.5 * h, k2, tmp);
    tdpt_rhs(n, K, E, t + 0.5 * h, Vm, tmp, k3, ph, z);
    tdpt_axpy(len, y, h, k3, tmp);
    tdpt_rhs(n, K, E, t + h, Vb, tmp, k4, ph, z);

    for (size_t i = 0; i < len; i++) {
      y[i].re +=
          h / 6.0 * (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re);
      y[i].im +=
          h / 6.0 * (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im);
    }

    // V at end of this step is V at start of next one
    cmatrix_t *swap = Va;
    Va = Vb;
    Vb = swap;
  }

  free(work);
  cmatrix_free(Va);
  cmatrix_free(Vm);
  cmatrix_free(Vb);

  return 0;
}

tdpt_result_t *tdpt_dyson(int n, const double *energies, tdpt_perturbation_fn V,
                          void *params, const cvector_t *psi0, double t0,
                          double t1, int steps, int max_order) {
  if (n < 1 || n > TDPT_MAX_DIM || !energies || !V || !psi0 || !psi0->data ||
      psi0->n != n || !isfinite(t0) || !isfinite(t1) || !(t1 > t0) ||
      steps < 4 || max_order < 0 || max_order > TDPT_MAX_ORDER) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    if (!isfinite(energies[i])) {
      return NULL;
    }
  }

  if (!tdpt_all_finite(psi0->data, (size_t)n)) {
    return NULL;
  }

  size_t len = (size_t)(max_order + 1) * (size_t)n;
  complex_t *y = calloc(len, sizeof *y);
  if (!y) {
    return NULL;
  }

  memcpy(y, psi0->data, (size_t)n * sizeof *y); // order 0, constant

  if (max_order > 0 && (tdpt_integrate(n, max_order, energies, V, params, t0,
                                       t1, steps, y) != 0 ||
                        !tdpt_all_finite(y, len))) {
    free(y);

    return NULL;
  }

  tdpt_result_t *res = tdpt_result_alloc(n, max_order, t0, t1, energies);
  if (!res) {
    free(y);

    return NULL;
  }

  for (int k = 0; k <= max_order; k++) {
    memcpy(res->orders[k].data, y + (size_t)k * n, (size_t)n * sizeof *y);
  }

  free(y);

  return res;
}

// Series assembly
static int tdpt_valid(const tdpt_result_t *r, int order) {
  return r && r->orders && r->energies && order >= 0 && order <= r->max_order;
}

cvector_t *tdpt_total_amplitude(const tdpt_result_t *result, int order) {
  if (!tdpt_valid(result, order)) {
    return NULL;
  }

  cvector_t *out = cvector_alloc(result->n);
  if (!out) {
    return NULL;
  }

  for (int i = 0; i < result->n; i++) {
    complex_t s = c_zero();
    for (int k = 0; k <= order; k++) {
      s = c_add(s, result->orders[k].data[i]);
    }

    out->data[i] = s;
  }

  return out;
}

cvector_t *tdpt_schrodinger_amplitude(const tdpt_result_t *result, int order) {
  cvector_t *out = tdpt_total_amplitude(result, order);
  if (!out) {
    return NULL;
  }

  for (int i = 0; i < result->n; i++) {
    double phase = -result->energies[i] * result->t1;

    out->data[i] = c_mul(out->data[i], (complex_t){cos(phase), sin(phase)});
  }

  return out;
}

double tdpt_probability(const tdpt_result_t *result, int order, int state) {
  if (!tdpt_valid(result, order) || state < 0 || state >= result->n) {
    return NAN;
  }

  complex_t s = c_zero();
  for (int k = 0; k <= order; k++) {
    s = c_add(s, result->orders[k].data[state]);
  }

  return c_abs2(s);
}

double tdpt_norm(const tdpt_result_t *result, int order) {
  cvector_t *amp = tdpt_total_amplitude(result, order);
  if (!amp) {
    return NAN;
  }

  double norm = 0.0;
  for (int i = 0; i < result->n; i++) {
    norm += c_abs2(amp->data[i]);
  }

  cvector_free(amp);

  return norm;
}

// Closed-form first-order building blocks
complex_t tdpt_phase_integral(double delta, double t) {
  if (!isfinite(delta) || !isfinite(t)) {
    return (complex_t){NAN, NAN};
  }

  double x = delta * t;
  if (fabs(x) < 1e-2) {
    // t * sum_k (i x)^k / (k + 1)!
    complex_t term = c_one();
    complex_t sum = c_one();
    for (int k = 1; k <= 8; k++) {
      // term_k = term_{k-1} * (i x) / (k + 1)
      complex_t next = {-term.im * x / (double)(k + 1),
                        term.re * x / (double)(k + 1)};
      term = next;
      sum = c_add(sum, term);
    }

    return (complex_t){t * sum.re, t * sum.im};
  }

  // (\exp^{i x} - 1) / (i \delta) = (\sin(x) + i (1 - \cos(x))) / \delta
  return (complex_t){sin(x) / delta, (1.0 - cos(x)) / delta};
}

complex_t tdpt_constant_first_order(complex_t Vfi, double omega_fi, double t) {
  complex_t integral = tdpt_phase_integral(omega_fi, t);
  complex_t v_int = c_mul(Vfi, integral);

  return (complex_t){v_int.im, -v_int.re}; // times -i
}

complex_t tdpt_harmonic_first_order(complex_t Vfi, double omega_fi,
                                    double Omega, double t) {
  complex_t sum = c_add(tdpt_phase_integral(omega_fi + Omega, t),
                        tdpt_phase_integral(omega_fi - Omega, t));
  complex_t v_int = c_mul(Vfi, sum);

  return (complex_t){0.5 * v_int.im, -0.5 * v_int.re}; // times -i/2
}

complex_t tdpt_gaussian_first_order(complex_t Vfi, double omega_fi,
                                    double tau) {
  if (!isfinite(omega_fi) || !isfinite(tau) || !(tau > 0.0)) {
    return (complex_t){NAN, NAN};
  }

  double mag =
      tau * sqrt(2.0 * M_PI) * exp(-0.5 * omega_fi * omega_fi * tau * tau);

  return (complex_t){mag * Vfi.im, -mag * Vfi.re}; // -i * mag * Vfi
}
