# `scattering`

> Generated from `physics/scattering.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Scattering theory: phase shifts, Born approximation, partial waves and the S-matrix.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/scattering.h>
```

**Start with:** [`phase_shift`](#phase_shift), [`born_amplitude`](#born_amplitude), [`born_cross_section`](#born_cross_section), [`s_matrix_element`](#s_matrix_element), [`partial_wave_amplitude`](#partial_wave_amplitude)

- Guide: [physics/scattering.md](../physics/scattering.md)
- Example: `examples/eg_13_scattering.c`
- Test (every call with expected results): `tests/test_scattering.c`
- Runtime help: `qmc_help(stdout, "scattering");` or `qmc-help scattering`

## Functions (8)

[`phase_shift`](#phase_shift), [`born_amplitude`](#born_amplitude), [`born_cross_section`](#born_cross_section), [`s_matrix_element`](#s_matrix_element), [`partial_wave_amplitude`](#partial_wave_amplitude), [`total_cross_section_partial_waves`](#total_cross_section_partial_waves), [`phase_shift_lippmann_schwinger`](#phase_shift_lippmann_schwinger), [`born_series_phase_shift`](#born_series_phase_shift)

### `phase_shift`

```c
double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N, double hbar_sq_2m);
```

```text
phase_shift: scattering phase shift for partial wave l, via Numerov
integration of radial equation + asymptotic matching.

l         : partial wave (angular momentum)
k         : wavenumber, E = hbar_sq_2m * k^2
V, params : central potential V(r) and its parameters
r_min     : small nonzero starting radius
r_max     : extend well beyond range of V
N         : number of radial grid points
hbar_sq_2m: \hbar^2/(2m)

Returns phase shift \delta_l via \atan2, in (-\pi, \pi].
```

### `born_amplitude`

```c
complex_t born_amplitude(potential_fn V, void *params, double k, double theta,
                         double r_max, int N, double hbar_sq_2m);
```

```text
Born approximation scattering amplitude f(\theta):
f(\theta) = - (1/hbar_sq_2m) * (1/q) * \int_0^\infty r * V(r) * \sin(q * r)
dr, q = 2k * \sin(\theta / 2)
```

### `born_cross_section`

```c
double born_cross_section(complex_t f_theta);
```

```text
Differential cross section: d\sigma / d\Omega = |f(\theta)|^2
```

### `s_matrix_element`

```c
complex_t s_matrix_element(double delta_l);
```

```text
S-matrix element for partial wave l: S_l = \exp(2 * i * \delta_l)
```

### `partial_wave_amplitude`

```c
complex_t partial_wave_amplitude(const double *delta_l, int l_max, double k,
                                 double theta);
```

```text
Full scattering amplitude from a partial-wave phase-shift set:
  f(\theta) = (1/k) * \sum_{l=0}^{l_max} (2l+1) * \exp^{i * \delta_l} *
\sin(\delta_l) P_l(\cos(\theta)) * \delta_l
must be an array of length l_max+1
```

### `total_cross_section_partial_waves`

```c
double total_cross_section_partial_waves(const double *delta_l, int l_max,
                                         double k);
```

```text
Total cross section from a partial-wave phase-shift set:
  \sigma_{tot} = (4 * \pi / k^2) * \sum_l (2l+1) \sin^2(\delta_l)
Cross-check two against optical theorem: \sigma_{tot} = (4 * \pi / k) *
Im[f(0)]
```

### `phase_shift_lippmann_schwinger`

```c
double phase_shift_lippmann_schwinger(int l, double k, potential_fn V,
                                      void *params, double r_min, double r_max,
                                      int N, double hbar_sq_2m);
```

```text
Exact (non-perturbative) phase shift via direct linear solve of discretized
Lippmann-Schwinger equation
```

### `born_series_phase_shift`

```c
double born_series_phase_shift(int l, double k, potential_fn V, void *params,
                               double r_min, double r_max, int N,
                               double hbar_sq_2m, int order);
```

```text
Phase shift via truncated Neumann/Born series to a given order:
  order = 1 reduces exactly to standard first-Born partial-wave phase shift
          (no iteration, u = hat_j_l(kr)).
  order >= 2 self-consistently iterates u through the LS kernel; with more
          decimal places.
  As order -> \infty this converges to phase_shift_lippmann_schwinger()'s
          exact result (for potentials where Born  series converges at all)
```
