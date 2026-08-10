# Rabi Oscillations

A two-level system driven near resonance, treated in the rotating-wave approximation (RWA). Implemented in `physics/rabi.h`.

## The Hamiltonian

In the rotating frame, in the $\{|1\rangle\text{ (ground)}, |2\rangle\text{ (excited)}\}$ basis:

$$
H = \frac{\hbar}{2}\begin{pmatrix} \Delta & \Omega \\ \Omega & -\Delta \end{pmatrix}
$$

Where:

- $\Delta = \omega_L - \omega_0$: Detuning between drive frequency $\omega_L$ and transition frequency $\omega_0$
- $\Omega$ is the (real, angular-frequency-units) Rabi coupling strength i.e. Rabi frequency
- $H$ is time-independent in the rotating frame, which is what makes an exact closed-form solution possible

## Transition Probability

Starting in ground state $\vert{}0\rangle$, the probability of occupying excited state $\vert{}1\rangle$ at time $t$ is:

$$
P_{0 \to 1}(t) = \frac{\vert{}\Omega_R\vert{}^2}{\Omega_R^2 + \Delta^2} \sin^2\left( \frac{\sqrt{\vert{}\Omega_R\vert{}^2 + \Delta^2}}{2} t \right)
$$

## Exact solutions

```c
/* P_e(t) = (\Omega^2 / \Omega_R^2) * \sin^2(\Omega_R * t/2),
   \Omega_R = \sqrt(\Omega^2 + \Delta^2) is the generalized Rabi frequency
*/
double rabi_excited_probability(double t, double Omega, double Delta);
```

Starting from the ground state at $t=0$, this is the exact excited-state population - the textbook Rabi flopping formula.

```c
/* U(t) = \cos(\Omega_R t/2) I - i*\sin(\Omega_R t/2) / \Omega_R * [[\Delta, \Omega], [\Omega, -\Delta]]
   Overwrites \psi in place.
   Returns 0 on success, -1 on invalid input.
*/
int rabi_evolve_exact(cvector_t *psi, double t, double Omega, double Delta);
```

Full unitary propagation of a 2-component spinor ($\psi[0]$ = ground amplitude, $\psi[1]$ = excited amplitude) using the closed-form evolution operator - not just the excited-state population, the complete complex state at time $t$.

```c
cvector_t *psi = cvector_alloc(2);
psi->data[0] = c_one();            // start in ground state
psi->data[1] = c_zero();
rabi_evolve_exact(psi, t, Omega, Delta);
```

## Numerical cross-check

```c
/* Treats H as a 2x2 real symmetric (trivially tridiagonal) energy
   Hamiltonian: diag = [\hbar * \Delta/2, -\hbar * \Delta/2], offdiag = [\hbar * \Omega/2] */
int rabi_evolve_numerical(cvector_t *psi, double hbar, double Omega,
                          double Delta, double dt, int steps);
```

Propagates the same 2-level system via the general Crank-Nicolson machinery (see [Crank-Nicolson Solver](../internals/crank_nicolson.md)) instead of the closed-form operator - useful as a validation path for the time-evolution infrastructure itself, since the exact answer is known and cheap to compute via `rabi_evolve_exact` for comparison.

## Running the Example

```sh
./build/eg_14_rabi
```
