/*
 * Benchmark: hand-rolled (Jacobi, via cmatrix_eigh_generic) vs LAPACK zheev
 * complex-Hermitian eigensolver, both reached through the single entry point
 * core/linalg/complex_eigh.c:cmatrix_eigh_complex().
 *
 * NOTE: cmatrix_eigh_complex() picks its backend at COMPILE time (#ifdef
 * USE_LAPACK), not runtime, so a single binary can only exercise one path. This
 * benchmark is therefore built twice : once with USE_LAPACK=0, once with
 * USE_LAPACK=1. And each run writes a tagged results line; combine the two
 * output files for the side-by-side comparison.
 * See benchmarks/README.md for the exact build/run recipe.
 *
 * Accuracy is checked with a residual test independent of any reference solver:
 *  for each returned eigenpair : (lambda_i, v_i),
 *  compute : max_i ||H v_i - lambda_i v_i|| / ||v_i||
 * directly against the input matrix, rather than diffing against a
 * simultaneously-computed "other" solver (which this compile-time split makes
 * impossible within one run anyway).
 *
 * USAGE:
 *   gcc/USE_LAPACK=0 build:
 *     make PLOT_BACKEND=NONE SANITIZE=0 build/bench_eigensolver
 *   USE_LAPACK=1 build (separate build dir to avoid stale objects):
 *     make PLOT_BACKEND=NONE SANITIZE=0 USE_LAPACK=1 build/bench_eigensolver
 *
 *   ./build/bench_eigensolver [size1 size2 ...]   (defaults below if none
 * given)
 */

#include "../core/complex.h"
#include "../core/linalg/linalg.h"
#include "../core/matrix.h"
#include "../core/random.h"
#include "../core/vector.h"
#include <bits/time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef USE_LAPACK
#define BACKEND_TAG "LAPACK"
#else
#define BACKEND_TAG "hand-rolled"
#endif

static double wall_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);

  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Random n x n complex Hermitian matrix: A_ij ~ N(0,1) + i N(0,1) off-diagonal
 * (with A_ji = conj(A_ij)), A_ii ~ N(0,1) real. */
static cmatrix_t *random_hermitian(int n, rng_state_t *rng) {
  cmatrix_t *H = cmatrix_alloc(n, n);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    CMAT(H, i, i) = c_real(rng_gaussian(rng));

    for (int j = i + 1; j < n; j++) {
      complex_t z = c_new(rng_gaussian(rng), rng_gaussian(rng));

      CMAT(H, i, j) = z;
      CMAT(H, j, i) = c_conj(z);
    }
  }

  return H;
}

/* max_i ||H v_i - lambda_i v_i|| / ||v_i|| over all n returned eigenpairs. */
static double max_residual(const cmatrix_t *H, const eigen_t *eig) {
  int n = eig->n;
  double worst = 0.0;

  cvector_t *v = cvector_alloc(n);
  if (!v) {
    return -1.0;
  }

  for (int k = 0; k < n; k++) {
    for (int i = 0; i < n; i++) {
      v->data[i] = CMAT(eig->eigenvectors, i, k);
    }

    double vnorm = cvector_norm(v);
    if (vnorm < 1e-300) {
      continue;
    }

    cvector_t *Hv = cmatrix_mv(H, v);
    if (!Hv) {
      cvector_free(v);

      return -1.0;
    }

    double lambda = eig->eigenvalues[k];
    double resid2 = 0.0;
    for (int i = 0; i < n; i++) {
      complex_t r = c_sub(Hv->data[i], c_scale(v->data[i], lambda));
      resid2 += c_abs2(r);
    }

    double resid = sqrt(resid2) / vnorm;
    if (resid > worst) {
      worst = resid;
    }

    cvector_free(Hv);
  }

  cvector_free(v);

  return worst;
}

int main(int argc, char **argv) {
  int default_sizes[] = {10, 20, 50, 100};
  int *sizes = default_sizes;
  int n_sizes = 4;

  if (argc > 1) {
    n_sizes = argc - 1;
    sizes = malloc(n_sizes * sizeof(int));

    for (int i = 0; i < n_sizes; i++) {
      sizes[i] = atoi(argv[i + 1]);
    }
  }

  printf("# Eigensolver benchmark -- backend: %s\n", BACKEND_TAG);
  printf("# Size | Time (ms) | Max residual ||Hv - lambda*v|| / ||v||\n");

  rng_state_t rng;
  rng_seed(&rng, 42ULL);

  for (int idx = 0; idx < n_sizes; idx++) {
    int n = sizes[idx];
    cmatrix_t *H = random_hermitian(n, &rng);

    if (!H) {
      fprintf(stderr, "allocation failed for n=%d\n", n);

      continue;
    }

    double t0 = wall_seconds();
    eigen_t *eig = cmatrix_eigh_complex(H);
    double elapsed_ms = (wall_seconds() - t0) * 1000.0;

    if (!eig) {
      fprintf(stderr, "cmatrix_eigh_complex failed for n=%d\n", n);

      cmatrix_free(H);

      continue;
    }

    double resid = max_residual(H, eig);

    printf("%d %.3f %.3e\n", n, elapsed_ms, resid);
    fflush(stdout);

    eigen_free(eig);
    cmatrix_free(H);
  }

  if (sizes != default_sizes) {
    free(sizes);
  }

  return 0;
}
