# `shift_invert`

> Generated from `core/shift_invert.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Shift-invert eigensolver for k interior eigenvalues and eigenvectors of a Hermitian sparse matrix.

**Category:** Core: linear algebra

```c
#include <qmc/core/shift_invert.h>
```

**Start with:** [`shift_invert_eigs`](#shift_invert_eigs)

- Guide: [internals/sparse_solvers.md](../internals/sparse_solvers.md)
- Test (every call with expected results): `tests/test_shift_invert.c`
- Runtime help: `qmc_help(stdout, "shift_invert");` or `qmc-help shift_invert`

## Types

### `shift_invert_result_t`

```text
Shift-invert eigensolver: finds k eigenvalues (and eigenvectors) of a
Hermitian sparse matrix nearest an arbitrary shift `\sigma`, via inverse
iteration on (A - \sigma * I)

NOTE: Unlike lanczos_eigs (sparse.h), which only reaches extremal
eigenvalues, this targets interior eigenvalues near any \sigma - at cost of
one dense LU factorization of (A - \sigma * I), so it's meant for matrices
small/dense enough that an n x n dense factorization is affordable, not for
very large sparse systems
```

```c
typedef struct {
  int n;              /* no. of eigenpairs actually returned (== k passed  to
                       * shift_invert_eigs); values/vectors are sized for
                       * this many entries */
  double *values;     /* eigenvalues nearest \sigma, ascending */
  cmatrix_t *vectors; /* n_A x k, eigenvectors as columns */
} shift_invert_result_t;
```

## Functions (2)

[`shift_invert_eigs`](#shift_invert_eigs), [`shift_invert_result_free`](#shift_invert_result_free)

### `shift_invert_eigs`

```c
shift_invert_result_t *shift_invert_eigs(const sparse_matrix_t *A, double sigma,
                                         int k, int max_iter, double tol);
```

```text
A        : Hermitian sparse matrix
sigma    : shift point; results are k eigenvalues closest to \sigma
k        : no. of eigenpairs wanted (1 <= k <= A->nrows)
max_iter : inverse-iteration steps per eigenpair (>= 1)
tol      : convergence threshold on successive-iterate overlap,
           |1 - |<v_prev|v>|| < tol

Returns NULL on invalid input, allocation failure, or if (A - \sigma * I)
is exactly singular (\sigma is exactly an eigenvalue of A - shift-invert
is undefined there; perturb sigma slightly and retry)

Later eigenpairs are found by Gram-Schmidt deflation against all earlier
ones, so they converge to progressively farther eigenvalues from sigma
rather than repeating same one
```

### `shift_invert_result_free`

```c
void shift_invert_result_free(shift_invert_result_t *res);
```

_No description in header yet._
