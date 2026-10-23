# `wkb`

> Generated from `physics/wkb.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

WKB approximation: Bohr-Sommerfeld quantisation and tunnelling transmission.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/wkb.h>
```

**Start with:** [`wkb_quantization`](#wkb_quantization), [`wkb_transmission`](#wkb_transmission), [`wkb_energy_harmonic`](#wkb_energy_harmonic), [`wkb_action_integral_harmonic`](#wkb_action_integral_harmonic), [`wkb_tunneling_rectangular`](#wkb_tunneling_rectangular)

- Guide: [physics/wkb.md](../physics/wkb.md)
- Example: `examples/eg_11_wkb.c`
- Test (every call with expected results): `tests/test_wkb.c`
- Runtime help: `qmc_help(stdout, "wkb");` or `qmc-help wkb`

## Functions (5)

[`wkb_quantization`](#wkb_quantization), [`wkb_transmission`](#wkb_transmission), [`wkb_energy_harmonic`](#wkb_energy_harmonic), [`wkb_action_integral_harmonic`](#wkb_action_integral_harmonic), [`wkb_tunneling_rectangular`](#wkb_tunneling_rectangular)

### `wkb_quantization`

```c
double wkb_quantization(potential_fn V, void *params, double m, double hbar,
                        int n, double E_min, double E_max, double x_min,
                        double x_max, int N, double tol);
```

```text
WKB quantization for bound states (Bohr-Sommerfeld):
\int_{x1}^{x2} \sqrt(2m(E - V(x))) dx = (n + 1/2) \pi \hbar
Returns energy for given n.
```

### `wkb_transmission`

```c
double wkb_transmission(potential_fn V, void *params, double m, double hbar,
                        double E, double x1, double x2, int N);
```

```text
WKB transmission coefficient for a barrier:
T = \exp(-2 \int_{x1}^{x2} \sqrt(2m(V(x) - E)) dx / \hbar)
```

### `wkb_energy_harmonic`

```c
double wkb_energy_harmonic(int n, double hbar, double omega);
```

```text
Wrapper: WKB energy for level n of the harmonic oscillator
V(x) = 0.5*omega^2*x^2
```

### `wkb_action_integral_harmonic`

```c
double wkb_action_integral_harmonic(double E, double hbar, double omega);
```

```text
Wrapper: full round-trip Bohr-Sommerfeld action integral
\oint p dx = 2 * \int_{-x_turn}^{x_turn} * \sqrt(2m(E - V(x))) dx
for the harmonic oscillator at given energy E (m=1).
Quantized levels satisfy action (2 * \pi * \hbar) = n + 1/2.
```

### `wkb_tunneling_rectangular`

```c
double wkb_tunneling_rectangular(double E, double V0, double width,
                                 double hbar_sq_2m);
```

```text
Wrapper: WKB tunneling probability through rectangular barrier of
height V0 and width, closed-form :
T = \exp(-2 * \kappa * width), \kappa = \sqrt((V0-E) / hbar_sq_2m)
where hbar_sq_2m = \hbar^2/(2m),
Returns 1.0 (fully transmitting) if V0 <= E
```
