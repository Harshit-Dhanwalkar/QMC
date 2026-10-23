# `qaoa`

> Generated from `physics/qaoa.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

QAOA for MaxCut and Ising cost functions: simulator, closed-form p=1 expectation and parameter optimisation.

**Category:** Quantum computing & information

```c
#include <qmc/physics/qaoa.h>
```

**Start with:** [`qaoa_problem_maxcut`](#qaoa_problem_maxcut), [`qaoa_problem_ising`](#qaoa_problem_ising), [`qaoa_problem_qubits`](#qaoa_problem_qubits), [`qaoa_problem_cost_table`](#qaoa_problem_cost_table), [`qaoa_problem_best_cost`](#qaoa_problem_best_cost)

- Example: `examples/eg_72_qaoa.c`
- Test (every call with expected results): `tests/test_qaoa.c`
- Runtime help: `qmc_help(stdout, "qaoa");` or `qmc-help qaoa`

## Overview

```text
Quantum Approximate Optimization Algorithm (QAOA): hybrid quantum-classical
solver for diagonal (classical) cost functions such as MaxCut and Ising
models (Reference: Farhi, Goldstone, Gutmann 2014)

NOTE: Depth-p ansatz, starting from uniform superposition |+>^n:
  |\psi(\gamma, \beta)> = \prod_{l=p-1}^{0} [ e^{-i \beta_l B} e^{-i \gamma_l
                          C} ] |+>^n,   B = \sum_i X_i
Where
 - C is a cost function that is diagonal in computational basis

Both layers are applied without ever forming a 2^n x 2^n matrix:
 - phase separator \exp(-i \gamma C) multiplies amplitude k by
   \exp(-i \gamma C[k])
 - mixer \exp(-i \beta B) = \prod_i RX_i(2 \beta) is n single-qubit gates
   (RX(\theta) = \exp(-i \theta X / 2))
Memory is 2^n complex amplitudes plus 2^n-double cost table, so problems
up to QAOA_MAX_QUBITS qubits are practical (unlike vqe.c's dense-matrix
Hamiltonians, which stop around 10-12 qubits)

NOTE: Bit convention matches every other qstate_* function: qubit 0 is most
significant bit of basis-state index, and a bit b maps to Ising spin
z = 1 - 2b (bit 0 -> z = +1, bit 1 -> z = -1)

NOTE: Every problem *maximizes* or *minimizes* its cost as recorded in
`maximize`, optimizer and approximation ratio respect that. MaxCut is a
maximization problem; Ising energies are minimization problems

NOTE: Optimizer, per restart: coordinate descent over 2p angles (each swept
via golden_section_minimize over a window around its current value) to get
into right basin, then a BFGS polish with finite-difference gradients
Coordinate descent alone stalls in narrow, strongly correlated gamma-beta
valleys of QAOA landscape once p >= 2. Restart 0 is a deterministic
linear-ramp (annealing-inspired) schedule; further restarts are random
Landscape is non-convex, so this finds a good local optimum, not a guaranteed
global one; more restarts raise odds

Reference: Farhi et al., arXiv:1411.4028; Wang, Hadfield, Jiang, Rieffel,
PRA 97, 022304 (2018) for closed-form p = 1 MaxCut expectation
```

## Types

### `qaoa_problem_t`

```c
typedef struct qaoa_problem qaoa_problem_t;
```

### `qaoa_result_t`

```c
typedef struct {
  double expectation;  // <\psi|C|\psi> at optimized angles
  double approx_ratio; // (<C> - C_worst) / (C_best - C_worst), in [0, 1]
  double optimum_prob; // total probability of measuring an optimal bitstring
  int best_bitstring;  // most probable basis state, MSB-first (qubit 0 = MSB)
  int depth;           // p
  double *gamma;       // p phase-separator angles (heap-allocated)
  double *beta;        // p mixer angles (heap-allocated)
} qaoa_result_t;
```

## Macros and constants

```c
#define QAOA_MAX_QUBITS 24
```

```text
Largest supported qubit count (2^24 amplitudes = 256 MiB of complex128)
```

## Functions (13)

[`qaoa_problem_maxcut`](#qaoa_problem_maxcut), [`qaoa_problem_ising`](#qaoa_problem_ising), [`qaoa_problem_free`](#qaoa_problem_free), [`qaoa_problem_qubits`](#qaoa_problem_qubits), [`qaoa_problem_cost_table`](#qaoa_problem_cost_table), [`qaoa_problem_best_cost`](#qaoa_problem_best_cost), [`qaoa_problem_worst_cost`](#qaoa_problem_worst_cost), [`qaoa_prepare_state`](#qaoa_prepare_state), [`qaoa_expectation`](#qaoa_expectation), [`qaoa_maxcut_p1_expectation`](#qaoa_maxcut_p1_expectation), [`qaoa_run`](#qaoa_run), [`qaoa_result_free`](#qaoa_result_free), [`qaoa_approximation_ratio`](#qaoa_approximation_ratio)

### `qaoa_problem_maxcut`

```c
qaoa_problem_t *qaoa_problem_maxcut(int n_qubits, const int *u, const int *v,
                                    const double *w, int n_edges);
```

```text
Build a weighted MaxCut problem on `n_qubits` vertices:
  C(z) = \sum_{e = (u,v)} w_e (1 - z_u z_v) / 2
i.e. total weight of edges cut by partition z
Maximization

u, v    : arrays of n_edges endpoint indices in [0, n_qubits), u[e] != v[e]
w       : array of n_edges edge weights, or NULL for unit weights
n_edges : number of edges (>= 1). Repeated edges are allowed and add

Returns NULL for invalid input (n_qubits outside [1, QAOA_MAX_QUBITS], NULL
endpoint arrays, self-loops, out-of-range or non-finite data) or allocation
failure. Free with qaoa_problem_free
```

### `qaoa_problem_ising`

```c
qaoa_problem_t *qaoa_problem_ising(int n_qubits, const double *h, const int *ci,
                                   const int *cj, const double *J,
                                   int n_couplings, double offset);
```

```text
Build a classical Ising problem:
  C(z) = offset + \sum_i h_i z_i + \sum_{c} J_c z_{i_c} z_{j_c}
Minimization (a ground-state energy problem)

h       : array of n_qubits local fields, or NULL for none
ci, cj  : arrays of n_couplings endpoint indices (ci[c] != cj[c])
J       : array of n_couplings coupling strengths
n_couplings may be 0 (then ci, cj, J may be NULL)

Returns NULL on invalid input or allocation failure
```

### `qaoa_problem_free`

```c
void qaoa_problem_free(qaoa_problem_t *problem);
```

_No description in header yet._

### `qaoa_problem_qubits`

```c
int qaoa_problem_qubits(const qaoa_problem_t *problem);
```

```text
Number of qubits, or 0 if `problem` is NULL
```

### `qaoa_problem_cost_table`

```c
const double *qaoa_problem_cost_table(const qaoa_problem_t *problem);
```

```text
Read-only view of the length-2^n cost table C[k] (k = basis-state index)

Returns NULL if `problem` is NULL. Valid until qaoa_problem_free
```

### `qaoa_problem_best_cost`

```c
double qaoa_problem_best_cost(const qaoa_problem_t *problem, int *best_index);
```

```text
Exact best/worst cost over all 2^n bitstrings, by a full scan of table (so
exact for every supported n). "Best" is maximum for maximization problems and
minimum for minimization problems; "worst" is other extreme. `*best_index`
(optional, may be NULL) receives lowest-index optimal bitstring

Both return 0.0 if `problem` is NULL
```

### `qaoa_problem_worst_cost`

```c
double qaoa_problem_worst_cost(const qaoa_problem_t *problem);
```

_No description in header yet._

### `qaoa_prepare_state`

```c
cvector_t *qaoa_prepare_state(const qaoa_problem_t *problem, int p,
                              const double *gamma, const double *beta);
```

```text
Prepare |\psi(\gamma, \beta)> for depth p >= 1 with p angles in each array

Returns an allocated normalized cvector_t, or NULL on invalid input
```

### `qaoa_expectation`

```c
double qaoa_expectation(const qaoa_problem_t *problem, int p,
                        const double *gamma, const double *beta);
```

```text
<\psi(\gamma, \beta)| C |\psi(\gamma, \beta)>

Returns 0.0 on invalid input
```

### `qaoa_maxcut_p1_expectation`

```c
double qaoa_maxcut_p1_expectation(const qaoa_problem_t *problem, double gamma,
                                  double beta);
```

```text
Closed-form depth-1 MaxCut expectation (Wang et al. 2018), valid for
UNWEIGHTED MaxCut problems (all weights == 1) of any graph, triangles
included:
  <C_uv> = 1/2 + 1/4 sin(4 \beta) \sin(\gamma)
                                   (\cos^d \gamma + \cos^e(\gamma))
           - 1/4 sin^2(2 \beta) \cos^{d + e - 2f}(\gamma)
                                   (1 - \cos^f(2 \gamma))

with d = deg(u) - 1, e = deg(v) - 1, f = number of triangles on edge (u,v);
<C> is sum over edges. Costs O(E * deg) rather than O(2^n), so it is usable
far beyond QAOA_MAX_QUBITS and doubles as an exact cross-check of
state-vector simulator

Returns NaN if `problem` is NULL, is not a MaxCut problem, or has any edge
weight other than 1 (or parallel edges), where formula does not apply
```

### `qaoa_run`

```c
qaoa_result_t qaoa_run(const qaoa_problem_t *problem, int p, int n_restarts,
                       int n_sweeps, double window, uint64_t seed);
```

```text
Optimize the 2p angles of a depth-p QAOA

p          : depth (>= 1)
n_restarts : number of starting points (>= 1); restart 0 is a linear ramp,
             rest are uniform random in gamma in [-\pi, \pi] and \beta in
             [-\pi/2, \pi/2] (fixed by `seed`)
n_sweeps   : coordinate-descent sweeps per restart (>= 1)
window     : half-width of each 1D golden-section search (> 0), radians

Returns best restart. On invalid input, returns a zeroed result with `\gamma
== beta == NULL`
```

### `qaoa_result_free`

```c
void qaoa_result_free(qaoa_result_t *result);
```

_No description in header yet._

### `qaoa_approximation_ratio`

```c
double qaoa_approximation_ratio(const qaoa_problem_t *problem,
                                double expectation);
```

```text
Approximation ratio (<C> - C_worst) / (C_best - C_worst) of an arbitrary
expectation value, which is offset-invariant and lies in [0, 1] for any
achievable <C>

Returns 1.0 if best == worst (constant cost) and 0.0 if problem` is NULL
```
