# `fft3d`

> Generated from `core/fft/fft3d.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

3D FFT by separable row-column-depth decomposition.

**Category:** Core: FFT

```c
#include <qmc/core/fft/fft3d.h>
```

**Start with:** [`fft3d`](#fft3d), [`ifft3d`](#ifft3d)

- Guide: [internals/fft.md](../internals/fft.md)
- Runtime help: `qmc_help(stdout, "fft3d");` or `qmc-help fft3d`

## Overview

```text
3D FFT via row-column-depth (separable) decomposition; transform every line
along z, then every line along y, then every line along x.
```

## Functions (2)

[`fft3d`](#fft3d), [`ifft3d`](#ifft3d)

### `fft3d`

```c
void fft3d(cvector_t *psi, int Nx, int Ny, int Nz);
```

_No description in the header yet._

### `ifft3d`

```c
void ifft3d(cvector_t *psi, int Nx, int Ny, int Nz);
```

_No description in the header yet._
