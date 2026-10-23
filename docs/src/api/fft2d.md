# `fft2d`

> Generated from `core/fft/fft2d.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

2D FFT by the separable row-column algorithm.

**Category:** Core: FFT

```c
#include <qmc/core/fft/fft2d.h>
```

**Start with:** [`fft2d`](#fft2d), [`ifft2d`](#ifft2d)

- Guide: [internals/fft.md](../internals/fft.md)
- Runtime help: `qmc_help(stdout, "fft2d");` or `qmc-help fft2d`

## Overview

```text
2D FFT built from 1D radix-2 fft()/ifft() via standard row-column (separable)
algorithm: FFT every row, then FFT every column
```

## Functions (2)

[`fft2d`](#fft2d), [`ifft2d`](#ifft2d)

### `fft2d`

```c
void fft2d(cvector_t *psi, int Nx, int Ny);
```

```text
Forward 2D FFT, unnormalized
```

### `ifft2d`

```c
void ifft2d(cvector_t *psi, int Nx, int Ny);
```

```text
Inverse 2D FFT, normalized by 1/(Nx*Ny)
```
