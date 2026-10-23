# `molecular_integrals`

> Generated from `physics/molecular_integrals.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Gaussian-type-orbital molecular integrals by McMurchie-Davidson, with built-in STO-3G data.

**Category:** Quantum chemistry (molecules)

```c
#include <qmc/physics/molecular_integrals.h>
```

**Start with:** [`boys_function`](#boys_function), [`boys_function_array`](#boys_function_array), [`basis_function_alloc`](#basis_function_alloc), [`molint_normalize_contraction`](#molint_normalize_contraction), [`basis_function_value`](#basis_function_value)

- Example: `examples/eg_43_molecular_integrals.c`
- Test (every call with expected results): `tests/test_molecular_integrals.c`
- Runtime help: `qmc_help(stdout, "molecular_integrals");` or `qmc-help molecular_integrals`

## Overview

```text
General Gaussian-type-orbital (GTO) molecular integrals via
McMurchie-Davidson (MD) Hermite-expansion scheme (Refreneces: McMurchie &
Davidson, J. Comput. Phys. 26, 218 (1974); Helgaker, Jorgensen & Olsen,
"Molecular Electronic-Structure Theory", ch. 9)
NOTE: This is lets second_quant.c / vqe.c build a real molecular
electronic-structure Hamiltonian (e.g. for H2 VQE) instead of only the toy
fermionic-lattice Hamiltonians second_quant_build_hopping_hamiltonian()
provides

Supports arbitrary Cartesian angular momentum (s, p, d, f, ...) and arbitrary
contraction (number of primitives), so this is a general basis-set /
general-molecule engine, not restricted to any one element or minimal basis

See basis_parser.c/.h for a general Gaussian94-format basis-set-file parser
built on the same basis_function_t representation

Units: atomic units throughout (bohr for length, Hartree for energy)
```

## Types

### `gto_primitive_t`

```text
A single primitive Cartesian Gaussian, unnormalized:
  g(r) = x^l y^m z^n \exp(-\alpha |r - center|^2)

l+m+n is this primitive's total angular momentum (0=s, 1=p, 2=d, ...)
Individual (l,m,n) triples distinguish e.g. px/py/pz (1,0,0)/(0,1,0)/
(0,0,1), or the six Cartesian d components
```

```c
typedef struct {
  int l, m, n;
  double alpha;
} gto_primitive_t;
```

### `basis_function_t`

```text
A contracted basis function: a fixed linear combination of primitives sharing
one center and one (l,m,n). coefficients[] are expected to already include
each primitive's own normalization and the overall contraction normalization
```

```c
typedef struct {
  int l, m, n;
  double center[3];
  int n_primitives;
  double *exponents;    /* length n_primitives */
  double *coefficients; /* length n_primitives */
} basis_function_t;
```

### `molecule_t`

```text
A molecule: point nuclear charges (atomic units, so charge is in units of
|e|, e.g. 1.0 for hydrogen) at fixed 3D centers. Same molecule_t can be
paired with any basis set
```

```c
typedef struct {
  int n_atoms;
  double *charge;      /* length n_atoms */
  double (*center)[3]; /* length n_atoms */
} molecule_t;
```

## Macros and constants

```c
#define MOLINT_ERI(eri, n, i, j, k, l)                                         \
  (eri)[(((size_t)(i) * (n) + (j)) * (n) + (k)) * (n) + (l)]
```

## Functions (25)

[`boys_function`](#boys_function), [`boys_function_array`](#boys_function_array), [`basis_function_alloc`](#basis_function_alloc), [`basis_function_free`](#basis_function_free), [`molint_normalize_contraction`](#molint_normalize_contraction), [`basis_function_value`](#basis_function_value), [`basis_function_gradient`](#basis_function_gradient), [`molecule_alloc`](#molecule_alloc), [`molecule_free`](#molecule_free), [`molecule_nuclear_repulsion`](#molecule_nuclear_repulsion), [`gto_overlap`](#gto_overlap), [`gto_kinetic`](#gto_kinetic), [`gto_nuclear_attraction`](#gto_nuclear_attraction), [`gto_eri`](#gto_eri), [`gto_overlap_grad_a`](#gto_overlap_grad_a), [`gto_kinetic_grad_a`](#gto_kinetic_grad_a), [`gto_nuclear_attraction_grad_a`](#gto_nuclear_attraction_grad_a), [`gto_nuclear_attraction_grad_C`](#gto_nuclear_attraction_grad_c), [`gto_eri_grad_a`](#gto_eri_grad_a), [`molecular_overlap_matrix`](#molecular_overlap_matrix), [`molecular_kinetic_matrix`](#molecular_kinetic_matrix), [`molecular_core_hamiltonian`](#molecular_core_hamiltonian), [`molecular_eri_tensor`](#molecular_eri_tensor), [`molint_basis_sto3g_h`](#molint_basis_sto3g_h), [`molint_basis_sto3g_li`](#molint_basis_sto3g_li)

### `boys_function`

```c
double boys_function(int n, double x);
```

```text
Boys function F_n(x) = \int^1 t^{2n} * \exp(-x * t^2) dt, universal building
block of every GTO-based Coulomb-type integral (nuclear attraction, ERIs)
Computed via the \exp(-x)*convergent-power-series form (all terms same sign,
no cancellation) for F_nmax, then the stable downward recursion :
  F_{n-1}(x) = (2x * F_n(x) + \exp(-x)) / (2n - 1) for n < nmax

(WARN: upward recursion is numerically unstable and must never be used)

Where
 boys_function      : single value F_n(x).
 boys_function_array: fills F[0..nmax] in one downward-recursion pass
```

### `boys_function_array`

```c
void boys_function_array(int nmax, double x, double *F  );
```

_No description in header yet._

### `basis_function_alloc`

```c
basis_function_t *basis_function_alloc(int l, int m, int n,
                                       const double center[3], int n_primitives,
                                       const double *exponents,
                                       const double *raw_coefficients);
```

```text
Allocate a basis function and copy in exponents/raw (unnormalized)
literature contraction coefficients

Returns NULL on allocation failure or n_primitives <= 0
```

### `basis_function_free`

```c
void basis_function_free(basis_function_t *bf);
```

_No description in header yet._

### `molint_normalize_contraction`

```c
void molint_normalize_contraction(basis_function_t *bf);
```

```text
In place:
  - folds each primitive's own Cartesian-Gaussian normalization constant into
    coefficients[i]
  - rescales the whole contraction so self-overlap integral(bf * bf) = 1
```

### `basis_function_value`

```c
double basis_function_value(const basis_function_t *bf, const double r[3]);
```

```text
Evaluate a contracted Cartesian GTO basis function's value at an arbitrary
point r (not just the integrals between two basis functions). Needed for any
real-space numerical-grid method (e.g. Becke-grid DFT)

Where
  - density n(r) = \sum_{pq} P_{pq} * \phi_{p}(r) * \phi_{q}(r) must be
    evaluated pointwise on a quadrature grid
```

### `basis_function_gradient`

```c
void basis_function_gradient(const basis_function_t *bf, const double r[3],
                             double grad[3]);
```

```text
Analytic gradient (d/dx, d/dy, d/dz) of a contracted Cartesian GTO basis
function at point r. Needed for any GGA-level (gradient-dependent) exchange-
correlation functional on a real-space grid, since the density gradient
 \nabla n(r) = 2 * \sum_{pq} D_pq * \nabla(\phi_p)(r) * \phi_q(r) (D
symmetric)
requires each basis function's own gradient, not just its value
Writes result into \nabla[3]
```

### `molecule_alloc`

```c
molecule_t *molecule_alloc(int n_atoms, const double *charge,
                           const double center[][3]);
```

_No description in header yet._

### `molecule_free`

```c
void molecule_free(molecule_t *mol);
```

_No description in header yet._

### `molecule_nuclear_repulsion`

```c
double molecule_nuclear_repulsion(const molecule_t *mol);
```

```text
Classical nuclear-nuclear repulsion energy :
 \Sum_{A<B} Z_A Z_B / |R_A - R_B|

piece of total electronic-structure energy that has nothing to do with basis
set or electrons
```

### `gto_overlap`

```c
double gto_overlap(const basis_function_t *a, const basis_function_t *b);
```

```text
<a|b>
```

### `gto_kinetic`

```c
double gto_kinetic(const basis_function_t *a, const basis_function_t *b);
```

```text
<a|-1/2 \nabla^2|b>
```

### `gto_nuclear_attraction`

```c
double gto_nuclear_attraction(const basis_function_t *a,
                              const basis_function_t *b,
                              const double center[3]);
```

```text
NOTE: Positive convention : integral itself, not yet multiplied by a nuclear
charge or a minus sign; molecular_nuclear_matrix below applies -Z_A per atom
and sums over atoms, matching one-electron Hamiltonian h = T - \sum_A Z_A *
this_integral)
```

### `gto_eri`

```c
double gto_eri(const basis_function_t *a, const basis_function_t *b,
               const basis_function_t *c, const basis_function_t *d);
```

```text
Full four-center, four-arbitrary-angular-momentum two-electron repulsion
integral
```

### `gto_overlap_grad_a`

```c
void gto_overlap_grad_a(const basis_function_t *a, const basis_function_t *b,
                        double grad[3]);
```

```text
Analytic derivative (gradient) integrals wrt a basis-function center or a
nuclear point-charge center, via the McMurchie-Davidson "increment/decrement"
identity (Helgaker, Jorgensen & Olsen eq. 9.3.8): reuses the integral
routines above with one center's angular momentum shifted by +-1, rather than
a separate derivative-Hermite-coefficient implementation

Every integral above is symmetric under swapping which argument comes first
(real GTOs, Hermitian/self-adjoint operators: <a|O|b> = <b|O|a>, (ab|cd) =
(ba|cd) = (cd|ab)), so the derivative wrt any other center is obtained by
simply reordering the arguments to put that center first:
  d<a|b>/dB   = gto_overlap_grad_a(b, a)
  d(ab|cd)/dB = gto_eri_grad_a(b, a, c, d)
  d(ab|cd)/dC = gto_eri_grad_a(c, d, a, b)
  d(ab|cd)/dD = gto_eri_grad_a(d, c, a, b)

NOTE: bra-pair order and ket-pair order may each be freely swapped, and bra
pair may swap with ket pair as a whole; a bra index may never be paired with
a ket index directly, since that changes which two centers' exponents sum
into the "p" of the bra vs "q" of the ket
```

### `gto_kinetic_grad_a`

```c
void gto_kinetic_grad_a(const basis_function_t *a, const basis_function_t *b,
                        double grad[3]);
```

_No description in header yet._

### `gto_nuclear_attraction_grad_a`

```c
void gto_nuclear_attraction_grad_a(const basis_function_t *a,
                                   const basis_function_t *b,
                                   const double center[3], double grad[3]);
```

_No description in header yet._

### `gto_nuclear_attraction_grad_C`

```c
void gto_nuclear_attraction_grad_C(const basis_function_t *a,
                                   const basis_function_t *b,
                                   const double center[3], double grad[3]);
```

```text
Derivative wrt the point-charge center itself (Hellmann-Feynman term):
present regardless of which atoms a/b happen to be centered on, since
moving that nucleus changes the 1/|r-C| operator felt by every AO pair
```

### `gto_eri_grad_a`

```c
void gto_eri_grad_a(const basis_function_t *a, const basis_function_t *b,
                    const basis_function_t *c, const basis_function_t *d,
                    double grad[3]);
```

_No description in header yet._

### `molecular_overlap_matrix`

```c
cmatrix_t *molecular_overlap_matrix(basis_function_t **basis, int n_basis);
```

```text
Whole-basis matrix/tensor builders: loop the above over every pair (or
quadruple) of functions in a basis_function_t** array of length n_basis

Returned cmatrix_t are real-valued (stored as complex_t with zero imaginary
part)
```

### `molecular_kinetic_matrix`

```c
cmatrix_t *molecular_kinetic_matrix(basis_function_t **basis, int n_basis);
```

_No description in header yet._

### `molecular_core_hamiltonian`

```c
cmatrix_t *molecular_core_hamiltonian(basis_function_t **basis, int n_basis,
                                      const molecule_t *mol);
```

```text
h_core = T - \sum_A Z_A * gto_nuclear_attraction(., ., R_A),
NOTE: full one-electron Hamiltonian matrix (kinetic + nuclear attraction to
every atom in mol, not just one center)
```

### `molecular_eri_tensor`

```c
double *molecular_eri_tensor(basis_function_t **basis, int n_basis);
```

```text
Flat row-major n_basis^4 array of (ij|kl) chemists'-notation ERIs, index (((i
* n_basis + j) * n_basis + k) * n_basis + l)
NOTE: Exploits 8-fold permutational symmetry (ij|kl) = (ji|kl) = (ij|lk) =
(kl|ij) = ... by computing each symmetry-distinct integral once and copying
to every equivalent slot, rather than calling gto_eri n_basis^4 times
independently

Returns NULL on n_basis<=0 or allocation failure
```

### `molint_basis_sto3g_h`

```c
basis_function_t *molint_basis_sto3g_h(const double center[3]);
```

```text
Builder: STO-3G hydrogen 1s (Refrenece: Szabo & Ostlund Table 3.7 / published
EMSL/Gaussian94 STO-3G parameters), already normalized

Returns NULL on allocation failure
```

### `molint_basis_sto3g_li`

```c
int molint_basis_sto3g_li(const double center[3], basis_function_t *out[5]);
```

```text
Published STO-3G lithium parameters (Refreneces: Hehre, Stewart & Pople 1969;
same numeric values as widely-used psi4/EMSL Basis Set Exchange sto-3g.gbs).
minimal-basis Li needs 5 basis functions : 1s (core), 2s (valence), and a
full 2p shell (2px, 2py, 2pz), even though 2p is unoccupied in the atomic
ground state, because STO-3G is defined as a full valence-shell minimal
basis. The 2s and 2p functions share the same 3 primitive exponents (a single
"SP shell", the minimal-basis-set convention) but have different contraction
coefficients

Fills out[0..4] with newly allocated, already-normalized basis functions in
the order 1s, 2s, 2px, 2py, 2pz

Returns 1 on success, 0 on allocation failure
```
