# `boson_sampling`

> Generated from `physics/boson_sampling.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Boson sampling: amplitudes and probabilities for N photons through an interferometer.

**Category:** Quantum computing & information

```c
#include <qmc/physics/boson_sampling.h>
```

**Start with:** [`boson_sampling_amplitude`](#boson_sampling_amplitude), [`boson_sampling_probability`](#boson_sampling_probability), [`beam_splitter_50_50`](#beam_splitter_50_50), [`dft_unitary`](#dft_unitary)

- Guide: [physics/boson_sampling.md](../physics/boson_sampling.md)
- Example: `examples/eg_25_boson_sampling.c`
- Test (every call with expected results): `tests/test_boson_sampling.c`
- Runtime help: `qmc_help(stdout, "boson_sampling");` or `qmc-help boson_sampling`

## Functions (4)

[`boson_sampling_amplitude`](#boson_sampling_amplitude), [`boson_sampling_probability`](#boson_sampling_probability), [`beam_splitter_50_50`](#beam_splitter_50_50), [`dft_unitary`](#dft_unitary)

### `boson_sampling_amplitude`

```c
complex_t boson_sampling_amplitude(const cmatrix_t *U, int M,
                                   const int *input_modes,
                                   const int *output_modes, int N);
```

```text
Boson sampling: transition amplitude/probability for N indistinguishable
photons through an MxM linear-optical network described by unitary U,
from an input Fock configuration to an output Fock configuration.

bosonic_permanent_value alone gives Perm(U_ST) / \sqrt(N! * prod_i s_i!).
Photon-number Fock-state transition amplitudes need normalization by both
input and output mode occupation factorials:
A = Perm(U_ST) / \sqrt(prod_i * s_i! * prod_j t_j!)

input_modes:  array of N input mode indices (0..M-1), one entry per
              photon. Repeated values = multiply-occupied input mode.
output_modes: array of N output mode indices (0..M-1), one entry per
              photon. Repeated values = multiply-occupied output mode.
N:            number of photons.
U:            MxM unitary.
```

### `boson_sampling_probability`

```c
double boson_sampling_probability(const cmatrix_t *U, int M,
                                  const int *input_modes,
                                  const int *output_modes, int N);
```

```text
|amplitude|^2 - transition probability.
```

### `beam_splitter_50_50`

```c
cmatrix_t *beam_splitter_50_50(void);
```

```text
50:50 beam-splitter unitary (2x2, Hadamard convention):
  U = (1 / \sqrt2) * [[1, 1], [1, -1]]
The theoretical Hong-Ou-Mandel test case.
Returns matrix
```

### `dft_unitary`

```c
cmatrix_t *dft_unitary(int M);
```

```text
MxM discrete Fourier transform matrix:
U_jk = (1 / \sqrt(M)) * \exp(2 * \pi * i * j * k / M)
Manifestly unitary; a
network for boson-sampling demonstrations/tests without needing a general
Haar-random-unitary generator. Returns matrix
```
