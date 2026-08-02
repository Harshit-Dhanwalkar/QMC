# Lindblad Master Equation (Open Quantum Systems)

The Lindblad (GKSL) equation describes the evolution of an open quantum system's density matrix under both unitary (Hamiltonian) dynamics and dissipative coupling to an environment. Natural units ($\hbar=m=1$). Implemented in `physics/lindblad.h`.

$$
$$\frac{d\rho}{dt} = -\frac{i}{\hbar} [H, \rho] + \sum_k \gamma_k \left( L_k \rho L_k^\dagger - \frac{1}{2} \{ L_k^\dagger L_k, \rho \} \right)$$
$$

<!-- $$ -->
<!-- \frac{d\rho}{dt} = -i[H,\rho] + \sum_k \mathcal{D}[L_k](\rho), \qquad -->
<!-- \mathcal{D}[L](\rho) = L\rho L^\dagger - \frac12\{L^\dagger L, \rho\} -->
<!-- $$ -->

Where:
* $H$: $N \times N$ System Hamiltonian
- $\rho$ : $N\times N$ Hermitian density matrix (`cmatrix_t`, $\text{Tr}(\rho)=1$);
* Each $L_k$: $N \times N$ Jump operators modeling interaction with environment i.e. jump operator encoding one dissipation channel.
* $\gamma_k$: Relaxation / dissipation rates
* $\{A, B\} = AB + BA$: Anti-commutator

## Density matrix basics

```c
// \rho = |\psi><\psi|, an NxN density matrix from a pure state vector
cmatrix_t *density_from_pure_state(const cvector_t *psi);

// Purity Tr(\rho^2). 1 for pure, <1 for mixed, 1/N for maximally mixed
double density_purity(const cmatrix_t *rho);

/* Von Neumann entropy S = -Tr(\rho \log2 \rho), in bits.
   0 for pure,\log2(N) for maximally mixed */
double density_von_neumann_entropy(cmatrix_t *rho);
```

## The Lindblad RHS and time-stepping

```c
/* -i[H,\rho] + \sum_k(L_k \rho L_k^\dagger - (1/2) * {L_k^\dagger L_k, \rho}).
   n_ops=0 (L may be NULL) gives closed-system unitary evolution. */
cmatrix_t *lindblad_rhs(const cmatrix_t *H, const cmatrix_t *rho, cmatrix_t **L,
                        int n_ops);

// One RK4 step of size dt applied to lindblad_rhs.
int lindblad_step_rk4(cmatrix_t *rho, const cmatrix_t *H, cmatrix_t **L,
                      int n_ops, double dt);

// `steps` repeated calls to lindblad_step_rk4.
int lindblad_evolve(cmatrix_t *rho, const cmatrix_t *H, cmatrix_t **L,
                    int n_ops, double dt, int steps);
```

Setting `n_ops = 0` recovers ordinary unitary (closed-system) evolution as a special case of the same machinery - useful for validating the RK4 integration itself against known unitary dynamics before trusting the dissipative terms.

## Building jump operators for qubit systems

```c
/* Embeds a 2x2 single-qubit operator (row-major: op[0]=<0|op|0>, op[1]=<0|op|1>,
   op[2]=<1|op|0>, op[3]=<1|op|1>) acting on qubit `target` into the full
   2^n_qubits x 2^n_qubits Hilbert space (identity on every other qubit). */
cmatrix_t *embed_single_qubit_op(const complex_t op[4], int n_qubits,
                                 int target);
```

Three standard single-qubit noise channels, already embedded and ready to use as an `L_k`:

```c
/* T1 decay (amplitude damping): L = \sqrt(\gamma) * \sigma_minus,
   \sigma_minus = |0><1|, population decays |1> -> |0> */
cmatrix_t *lindblad_amplitude_damping_op(int n_qubits, int target, double gamma);

// T2 (pure dephasing): L = \sqrt(\gamma/2) * \sigma_z
cmatrix_t *lindblad_dephasing_op(int n_qubits, int target, double gamma);

// Bit flip: L = \sqrt(\gamma) * \sigma_x
cmatrix_t *lindblad_bitflip_op(int n_qubits, int target, double gamma);
```

`gamma` is the channel rate (1/time, natural units).

```c
cmatrix_t *rho = density_from_pure_state(psi);
cmatrix_t *L_ops[1] = { lindblad_amplitude_damping_op(1, 0, gamma) };
lindblad_evolve(rho, H, L_ops, 1, dt, steps);
```

## Measurement

```c
/* Samples an outcome in [0, \rho->nrows) with probability \rho[i][i].re,
   using caller-supplied uniform u in [0,1). Collapses \rho in place to the
   pure-state projector |outcome><outcome|. Returns outcome, or -1 on
   invalid input (\rho not square). */
int density_measure_computational_basis(cmatrix_t *rho, double u);
```

Same caller-owns-the-RNG convention as [`qstate_measure`](qubits.md#measurement), for reproducible tests.

## Running the Example

```sh
./build/eg_19_lindblad
```
