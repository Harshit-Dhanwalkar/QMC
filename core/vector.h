#ifndef QMC_VECTOR_H
#define QMC_VECTOR_H

#include "complex.h"

typedef struct {
  complex_t *data;
  int n;
} cvector_t;

/* Allocate/Free */
cvector_t *cvector_alloc(int n);
void cvector_free(cvector_t *v);
cvector_t *cvector_copy(const cvector_t *v);

/* Operations */
complex_t cvector_dot(const cvector_t *a, const cvector_t *b);
double cvector_norm(const cvector_t *v);
void cvector_normalize(cvector_t *v);
double cvector_expect(const cvector_t *psi, const cvector_t *op_psi);

/* Utilities */
void cvector_print(const cvector_t *v, const char *label);
void cvector_fill(cvector_t *v, complex_t val);

#endif
