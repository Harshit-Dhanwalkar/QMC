# Linear Algebra Core

The library provides a set of linear algebra routines for complex matrices and vectors, optimized for quantum mechanics.

## Matrix and Vector Types

- `cvector_t` - complex vector (size `n`, data `complex_t*`)
- `cmatrix_t` - complex matrix (row-major, `nrows`x`ncols`)

## Basic Operations

- Allocation/free: `cmatrix_alloc`, `cmatrix_free`
- Copy: `cmatrix_copy`
- Access: `CMAT(m,i,j)` macro
- Multiplication, transpose, adjoint, scaling

## Linear Solvers

- LU decomposition: `cmatrix_lu_decomp`, `cmatrix_solve` (uses `lu.c`)
- Complex tridiagonal solver (Thomas algorithm) used in Crank-Nicolson.

## Eigenvalue Solvers

The main function is `cmatrix_eigh()` which computes all eigenvalues and eigenvectors of a Hermitian matrix.

### Implementation

Currently uses the **Jacobi eigenvalue algorithm** for real symmetric matrices (after assuming real Hermitian). It iteratively rotates off-diagonal elements until the matrix is diagonal.

```c
eigen_t *cmatrix_eigh(cmatrix_t *A) {
    // Copy A, initialize eigenvector matrix to identity
    // Sweep over all off-diagonal pairs, apply Jacobi rotations
    // Sort eigenvalues and eigenvectors in ascending order
}
```

## Sparse Matrices

For large systems, sparse matrix support (`sparse.h`) is provided, but not yet fully integrated. // TODO

## FFT

FFT routines are in `core/fft/` for Fourier transforms (momentum space).

## Usage Example

```c
cmatrix_t *H = cmatrix_alloc(N, N);
// Fill H ...
eigen_t *eig = cmatrix_eigh(H);
for (int i=0; i<N; i++)
    printf("E_%d = %f\n", i, eig->eigenvalues[i]);
eigen_free(eig);
cmatrix_free(H);
```

TODO: Later implement LAPACK
