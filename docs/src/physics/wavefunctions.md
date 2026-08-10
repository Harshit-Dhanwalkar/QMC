# Wave Functions

Properties of quantum mechanical wavefunctions, including normalization, expectation values, and basis transformations.

## The `wavefunction_t` type

The physics layer wraps a grid-aware wavefunction in `physics/wavefn.h`,
rather than passing a raw `cvector_t` and a separate `dx` around everywhere:

```c
typedef struct {
  cvector_t *psi; // wavefunction values on grid
  double *x;      // position grid
  double dx;      // grid spacing
  int n;          // number of grid points
} wavefunction_t;

wavefunction_t *wavefunction_alloc(int n);
void wavefunction_free(wavefunction_t *wf);
wavefunction_t *wavefunction_copy(const wavefunction_t *wf);
```

> The header comments `x` as "position grid (owned)" with your own inline `// NOTE: what is 'owned'?` still sitting next to it. Presumably this means `wavefunction_t` allocates and owns its own copy of the grid (so `wavefunction_free` frees it too) rather than borrowing a caller-supplied pointer - but that's a guess pending confirmation, since it's flagged as an open question in header.

## Normalization

The wavefunction $\psi(x)$ must be normalized:

$$
\int_{-\infty}^{\infty} |\psi(x)|^2 dx = 1
$$

```c
void wavefunction_normalize(wavefunction_t *wf);
```

This is the physics-layer entry point and correctly has access to `wf->dx` internally, so there's no ambiguity here (unlike the lower-level `core/utils.h::normalize_wavefunction(psi)`, which takes no `dx` - see the [Linear Algebra Core](../internals/linalg.md) note on that discrepancy; whether one calls the other, or they're independent, isn't clear without the `.c` file).

## Probability Density

```c
double *wavefunction_prob_density(const wavefunction_t *wf); // |psi|^2, per grid point
double wavefunction_prob_in_interval(const wavefunction_t *wf, double a, double b);
```

`wavefunction_prob_in_interval` integrates $|\psi|^2$ between `a` and `b` - useful for e.g. tunneling probability past a barrier.

## Expectation Values

For an operator $\hat{O}$, the expectation value in state $\psi$ is:

$$
\langle \hat{O} \rangle = \int \psi^*(x)\, \hat{O}\psi(x)\, dx
$$

```c
double wavefunction_expect_x(const wavefunction_t *wf);   // <x>, position space
double wavefunction_expect_x2(const wavefunction_t *wf);  // <x^2>
double wavefunction_expect_p(const wavefunction_t *wf);   // <p>, via FFT
double wavefunction_expect_p2(const wavefunction_t *wf);  // <p^2>
```

`wavefunction_expect_p` and `wavefunction_expect_p2` go through the FFT interface (`core/fft/fft.h`) to reach momentum space - see `position_to_momentum` in `core/utils.h`.

## Uncertainty

```c
double wavefunction_delta_x(const wavefunction_t *wf);
double wavefunction_delta_p(const wavefunction_t *wf);
double wavefunction_uncertainty_product(const wavefunction_t *wf);
```

These are thin wrappers over the expectation values above ($\Delta x = \sqrt{\langle x^2\rangle - \langle x\rangle^2}$, etc.). For the combined result as a single struct, see [Uncertainty Principle](uncertainty.md) and `compute_uncertainties()`.

## I/O

```c
void wavefunction_save(const wavefunction_t *wf, const char *filename);
void wavefunction_save_prob(const wavefunction_t *wf, const char *filename);
```

## Example (Gaussian wavepacket)

For the minimum-uncertainty Gaussian wavepacket:

$$
\phi(x)=\frac{1}{(2\pi\sigma^2)^{1/4}} \exp\left(-\frac{(x-x_0)^2}{4\sigma^2} + \frac{ip_0 x}{\hbar}\right)
$$

Properties:

- $\langle x \rangle = x_0$
- $\langle p \rangle = p_0$
- $\Delta x = \sigma$
- $\Delta p = \hbar/(2\sigma)$
