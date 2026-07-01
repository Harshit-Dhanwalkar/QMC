#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

cvector_t *cvector_alloc(int n) {
  cvector_t *v = malloc(sizeof(cvector_t));
  if (!v)
    return NULL;

  v->data = calloc(n, sizeof(complex_t));
  if (!v->data) {
    free(v);
    return NULL;
  }
  v->n = n;
  return v;
}

void cvector_free(cvector_t *v) {
  if (!v)
    return;
  free(v->data);
  free(v);
}

cvector_t *cvector_copy(const cvector_t *v) {
  cvector_t *copy = cvector_alloc(v->n);
  if (!copy)
    return NULL;

  memcpy(copy->data, v->data, v->n * sizeof(complex_t));
  return copy;
}

complex_t cvector_dot(const cvector_t *a, const cvector_t *b) {
  if (a->n != b->n) {
    fprintf(stderr, "Error: vector size mismatch in dot product\n");
    return c_zero();
  }

  complex_t result = c_zero();
  for (int i = 0; i < a->n; i++) {
    result = c_add(result, c_mul(c_conj(a->data[i]), b->data[i]));
  }
  return result;
}

double cvector_norm(const cvector_t *v) {
  double sum = 0.0;
  for (int i = 0; i < v->n; i++) {
    sum += c_abs2(v->data[i]);
  }
  return sqrt(sum);
}

void cvector_normalize(cvector_t *v) {
  double norm = cvector_norm(v);
  if (norm > 1e-15) {
    for (int i = 0; i < v->n; i++) {
      v->data[i] = c_scale(v->data[i], 1.0 / norm);
    }
  }
}

double cvector_expect(const cvector_t *psi, const cvector_t *op_psi) {
  complex_t result = cvector_dot(psi, op_psi);
  return result.re;
}

void cvector_print(const cvector_t *v, const char *label) {
  if (label)
    printf("%s:\n", label);
  for (int i = 0; i < v->n; i++) {
    printf("  [%d] = %.6e + %.6e i\n", i, v->data[i].re, v->data[i].im);
  }
}

void cvector_fill(cvector_t *v, complex_t val) {
  for (int i = 0; i < v->n; i++) {
    v->data[i] = val;
  }
}
