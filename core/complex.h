#ifndef QMC_COMPLEX_H
#define QMC_COMPLEX_H

#include <math.h>

typedef struct {
  double re, im;
} complex_t;

/* Basic operations */
static inline complex_t c_new(double re, double im) {
  return (complex_t){re, im};
}
static inline complex_t c_one(void) { return (complex_t){1.0, 0.0}; }
static inline complex_t c_zero(void) { return (complex_t){0.0, 0.0}; }
static inline complex_t c_neg(complex_t z) { return (complex_t){-z.re, -z.im}; }
static inline complex_t c_real(double x) { return (complex_t){x, 0.0}; }
static inline complex_t c_imag(double y) { return (complex_t){0.0, y}; }
static inline complex_t c_conj(complex_t z) { return (complex_t){z.re, -z.im}; }
// NOLINTNEXTLINE(clang-analyzer-core.UndefinedBinaryOperatorResult)
static inline double c_abs(complex_t z) { return sqrt(z.re * z.re + z.im * z.im);}
static inline double c_abs2(complex_t z) { return z.re * z.re + z.im * z.im; }
static inline double c_arg(complex_t z) { return atan2(z.im, z.re); }

/* Arithmetic */
static inline complex_t c_add(complex_t a, complex_t b) {
  return (complex_t){a.re + b.re, a.im + b.im};
}
static inline complex_t c_sub(complex_t a, complex_t b) {
  return (complex_t){a.re - b.re, a.im - b.im};
}
static inline complex_t c_mul(complex_t a, complex_t b) {
  return (complex_t){a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
static inline complex_t c_div(complex_t a, complex_t b) {
  double denom = c_abs2(b);
  return (complex_t){(a.re * b.re + a.im * b.im) / denom,
                     (a.im * b.re - a.re * b.im) / denom};
}
static inline complex_t c_scale(complex_t z, double s) {
  return (complex_t){z.re * s, z.im * s};
}
static inline complex_t c_exp(complex_t z) {
  double exp_re = exp(z.re);
  return (complex_t){exp_re * cos(z.im), exp_re * sin(z.im)};
}

/* Polar coorindates */
static inline double c_abs_sqr(complex_t z) {
  return z.re * z.re + z.im * z.im;
}
static inline complex_t c_from_polar(double r, double theta) {
  return (complex_t){r * cos(theta), r * sin(theta)};
}

#endif
