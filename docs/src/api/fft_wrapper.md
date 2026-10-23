# `fft_wrapper`

> Generated from `core/fft/fft_wrapper.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

FFT wrapper that dispatches to the built-in implementation or to FFTW.

**Category:** Core: FFT

```c
#include <qmc/core/fft/fft_wrapper.h>
```

**Start with:** [`fft_wrapper`](#fft_wrapper), [`ifft_wrapper`](#ifft_wrapper), [`fft_wrapper_normalized`](#fft_wrapper_normalized), [`ifft_wrapper_normalized`](#ifft_wrapper_normalized)

- Guide: [internals/fft.md](../internals/fft.md)
- Runtime help: `qmc_help(stdout, "fft_wrapper");` or `qmc-help fft_wrapper`

## Functions (4)

[`fft_wrapper`](#fft_wrapper), [`ifft_wrapper`](#ifft_wrapper), [`fft_wrapper_normalized`](#fft_wrapper_normalized), [`ifft_wrapper_normalized`](#ifft_wrapper_normalized)

### `fft_wrapper`

```c
void fft_wrapper(cvector_t *x);
```

```text
FFT wrapper that uses either internal or external (FFTW) library
```

### `ifft_wrapper`

```c
void ifft_wrapper(cvector_t *x);
```

_No description in the header yet._

### `fft_wrapper_normalized`

```c
void fft_wrapper_normalized(cvector_t *x);
```

_No description in the header yet._

### `ifft_wrapper_normalized`

```c
void ifft_wrapper_normalized(cvector_t *x);
```

_No description in the header yet._
