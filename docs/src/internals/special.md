# Special Functions

A collection of commonly used special functions in quantum mechanics: orthogonal polynomials, spherical harmonics, spherical Bessel functions, and factorials. Implemented in `core/special/` (headers in `core/special/special.h`).

These are used throughout the physics modules - e.g., Hermite polynomials for the harmonic oscillator, Laguerre polynomials for hydrogen, Legendre polynomials and spherical harmonics for angular momentum, and spherical Bessel functions for scattering.

---

## Hermite Polynomials

```c
double hermite(int n, double x);
void hermite_array(int n, double *x, int m, double *H);
double hermite_deriv(int n, double x);
double hermite_zeros(int n, int k);
```

- `hermite(n, x)` - returns $H_n(x)$ using the recurrence:
  $$H_0=1,  H_1=2x,  H_{n+1}=2xH_n − 2nH_{n−1}$$
- `hermite_array` - fills an array $H[0..m-1]$ with $H_n(x_i)$ for a list of points.
- `hermite_deriv` - derivative $Hn′(x)=2nH_{n−1}(x)$.
- `hermite_zeros` - returns the _k_-th zero of $H_n$ (0‑based), computed via Newton refinement from an asymptotic approximation. Useful for Gauss‑Hermite quadrature.

## Laguerre Polynomials

```c
double laguerre(int n, double alpha, double x);
void laguerre_array(int n, double alpha, double *x, int N, double *L);
```

