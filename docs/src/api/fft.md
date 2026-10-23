# `fft`

> Generated from `core/fft/fft.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

1D radix-2 Cooley-Tukey FFT and inverse for power-of-two lengths.

**Category:** Core: FFT

```c
#include <qmc/core/fft/fft.h>
```

**Start with:** [`fft`](#fft), [`ifft`](#ifft), [`fft_normalized`](#fft_normalized), [`ifft_normalized`](#ifft_normalized), [`fft_shift`](#fft_shift)

- Guide: [internals/fft.md](../internals/fft.md)
- Test (every call with expected results): `tests/test_fft.c`
- Runtime help: `qmc_help(stdout, "fft");` or `qmc-help fft`

## Overview

```text
Fast Fourier Transform (Cooley‑Tukey, radix‑2, iterative)
Works for power‑of‑two lengths.
```

## Functions (6)

[`fft`](#fft), [`ifft`](#ifft), [`fft_normalized`](#fft_normalized), [`ifft_normalized`](#ifft_normalized), [`fft_shift`](#fft_shift), [`position_to_momentum`](#position_to_momentum)

### `fft`

```c
void fft(cvector_t *x);
```

```text
In‑place forward FFT (no normalization)
```

### `ifft`

```c
void ifft(cvector_t *x);
```

```text
In‑place inverse FFT (no normalization)
```

### `fft_normalized`

```c
void fft_normalized(cvector_t *x);
```

```text
Normalized forward FFT (1 / \sqrt(N) factor)
```

### `ifft_normalized`

```c
void ifft_normalized(cvector_t *x);
```

```text
Normalized inverse FFT (1 / \sqrt(N) factor)
```

### `fft_shift`

```c
void fft_shift(cvector_t *x);
```

```text
Utility: shift zero‑frequency component to centre of array
```

### `position_to_momentum`

```c
cvector_t *position_to_momentum(const cvector_t *psi_x, double dx);
```

```text
Convert position‑space wavefunction to momentum space.
dx : grid spacing in position space.
The output array is ordered as k = 2 * \pi ** n / (N * dx) with n=0..N-1.
Normalisation is such that \int |\phi(x)|^2 dx = \int|\phi(k)|^2 dk.

Returns a newly allocated cvector_t in momentum space.
```
