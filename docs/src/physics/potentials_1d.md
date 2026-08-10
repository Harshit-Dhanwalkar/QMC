# 1D Potentials

Potentials implemented in `src/physics/potentials.c` / `physics/potentials.h`.

## API shape

Unlike an array-filling function per potential, QMC uses single-point
evaluators plus a generic array filler:

```c
typedef double (*potential_fn)(double x, void *params);

void potential_array(double *x, int n, potential_fn V, void *params,
                     double *V_out);
```

Every potential below has this signature - `V(x, params)` - and gets
evaluated across a grid via `potential_array`, rather than each potential
having its own array-filling function.

## Implemented Potentials

- **Infinite Square Well** - `V_infinite_well(x, params)`
  `params`: `double *a` (**half-width**). Well is `|x| < a`, centered at the
  origin - not `[0, L]`.

- **Finite Square Well** - `V_finite_well(x, params)`
  `params`: `struct { double a, V0; }`. `-V0` inside `|x| < a`, zero outside.

- **Harmonic Oscillator** - `V_harmonic(x, params)`
  `params`: `double *omega`. Computes `0.5*omega^2*x^2` - **no mass term**,
  consistent with the project's natural-units convention ($\hbar = m = 1$). If
  you need a physical mass, build the array manually instead of going through
  this function.

- **Step Potential** - `V_step(x, params)`
  `params`: `double *V0`. `0` for `x < 0`, `V0` for `x >= 0`.

- **Rectangular Barrier** - `V_barrier(x, params)`
  `params`: `struct { double a, V0; }`. `V0` inside `0 < x < a`, zero outside.

- **Coulomb Potential** - `V_coulomb(r, params)`
  For hydrogen-like systems. _`params` isn't documented in the header -
  needs confirming against `potentials.c` (unlike every other entry here, it
  has no accompanying params comment)._

- **Yukawa Potential** - `V_yukawa(r, params)`
  Screened Coulomb, $V = -g\,e^{-\mu r}/r$. _`params` also undocumented in
  the header - likely `struct { double g, mu; }` by analogy with the others,
  but unconfirmed._

- **Morse Potential** - `V_morse(x, params)`
  $D_e[1 - e^{-a(x-x_0)}]^2$, for molecular vibration. Params struct not
  documented in the header.

> **Not currently in `potentials.h`:** a double-well potential and a > Kronig-Penney periodic potential were both referenced in the previous version of this page but have no declaration here. Either they live in a different header I haven't seen yet, or they were aspirational and haven't landed - worth confirming either way before this page claims they exist.

## Usage Example

```c
int N = 1000;
double *x = linspace(-10, 10, N);
double *V = malloc(N * sizeof(double));

double omega = 1.0;
potential_array(x, N, V_harmonic, &omega, V);
```

Finite well, using the anonymous params struct inline:

```c
struct { double a, V0; } fw = { .a = 2.0, .V0 = 5.0 };
potential_array(x, N, V_finite_well, &fw, V);
```

All potentials are used with the Numerov integrator (see
[Numerov Integrator](../internals/numerov.md)) or tridiagonal matrix
diagonalization (see [Linear Algebra Core](../internals/linalg.md)) to solve
the Schrödinger equation.
