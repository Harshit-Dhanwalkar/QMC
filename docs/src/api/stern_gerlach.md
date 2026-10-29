# `stern_gerlach`

> Generated from `physics/stern_gerlach.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Stern-Gerlach analyzers for spin-1/2: sequential measurements, partial coherence and screen patterns.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/stern_gerlach.h>
```

**Start with:** [`sg_chain`](#sg_chain), [`sg_screen_density`](#sg_screen_density), [`sg_prob_up`](#sg_prob_up)

- Runtime help: `qmc_help(stdout, "stern_gerlach");` or `qmc-help stern_gerlach`

## Overview

```text
Stern-Gerlach analyzers for spin-1/2: sequential measurements, partial
coherence and screen pattern

Spin state of beam is its Bloch (polarisation) vector r (|r| <= 1, a mixed
state when < 1), carried together with beam intensity I. An analyzer along
unit vector n has two outputs with intensities I (1 \pm r.n) / 2. Blocking
one output is a projective measurement: transmitted beam has r = \pm n
Passing both outputs through a reversed gradient recombines beam; with
perfect recombination (coherence 1) spin state is untouched and with
imperfect recombination (coherence c in [0, 1], overlap of two spatial
packets) only component of r along n survives intact: r -> (r.n) n + c (r -
(r.n) n)

Screen pattern: a magnet of gradient G sends up and down components to z =
\pm d/2 with Gaussian profiles of standard deviation \sigma (packet width at
screen). Packets overlap by <f|g> = \exp(-d^2 / (8 \sigma^2)); a classical
magnetic moment with random orientation would give a uniform spread over
[-d/2, d/2] instead of two spots
Atomic Units \hbar = m = 1
```

## Types

### `sg_vec_t`

```c
typedef struct {
  double x, y, z;
} sg_vec_t;
```

### `?`

```c
enum { SG_PASS_UP = 0, SG_PASS_DOWN = 1, SG_PASS_BOTH = 2 };
```

### `sg_analyzer_t`

```c
typedef struct {
  double theta,
      phi;          /* analyzer axis, polar angle from +z and azimuth from +x */
  int mode;         /* SG_PASS_UP, SG_PASS_DOWN or SG_PASS_BOTH */
  double coherence; /* recombination quality for SG_PASS_BOTH, in [0, 1] */
} sg_analyzer_t;
```

## Functions (8)

[`sg_axis`](#sg_axis), [`sg_prob_up`](#sg_prob_up), [`sg_chain`](#sg_chain), [`sg_coherence`](#sg_coherence), [`sg_gauss`](#sg_gauss), [`sg_screen_density`](#sg_screen_density), [`sg_classical_density`](#sg_classical_density), [`sg_width`](#sg_width)

### `sg_axis`

```c
sg_vec_t sg_axis(double theta, double phi);
```

```text
Unit vector (\sin(\theta) \cos(\phi) \sin(\theta) \sin(\phi) \cos(\theta))
```

### `sg_prob_up`

```c
double sg_prob_up(sg_vec_t r, sg_vec_t n);
```

```text
Probability of + outcome along n for polarisation r: (1 + r.n) / 2
```

### `sg_chain`

```c
int sg_chain(sg_vec_t r0, double i0, const sg_analyzer_t *a, int n,
             double intensity[], sg_vec_t r[]);
```

```text
Sends a beam of intensity i0 and polarisation r0 through analyzers in turn
Fills intensity[0..n] and r[0..n] (entry 0 is source, entry k beam after
analyzer k; after an extinguished beam r keeps its last value)
Returns 0, or -1 for non-finite input, |r0| > 1, a negative i0, n < 0 or a
coherence outside [0, 1]
```

### `sg_coherence`

```c
double sg_coherence(double d, double sigma);
```

```text
Overlap \exp(-d^2 / (8 \sigma^2)) of 2 outgoing packets (separation d,
density width \sigma)
```

### `sg_gauss`

```c
double sg_gauss(double z, double c, double s);
```

```text
Gaussian density of centre c and standard deviation s
```

### `sg_screen_density`

```c
double sg_screen_density(double z, double pup, double d, double sigma);
```

```text
Pattern on screen:
 pup N(z; +d/2, \sigma) + (1 - pup) N(z; -d/2, \sigma)
```

### `sg_classical_density`

```c
double sg_classical_density(double z, double d, double sigma);
```

```text
Classical smear: uniform orientation, z = (d/2) \cos(angle), convolved with
packet width
```

### `sg_width`

```c
double sg_width(double sigma0, double t);
```

```text
Width after free spreading for time t (\hbar = m = 1):
 \sigma_0 sqrt(1 + (t / (2 \sigma_0^2))^2)
```
