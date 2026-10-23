# `utils`

> Generated from `core/utils.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Grid and array utilities: linspace, logspace, range, wavefunction normalisation and norms.

**Category:** Core: data types & utilities

```c
#include <qmc/core/utils.h>
```

**Start with:** [`linspace`](#linspace), [`logspace`](#logspace), [`range`](#range), [`normalize_wavefunction`](#normalize_wavefunction), [`compute_norm_squared`](#compute_norm_squared)

- Runtime help: `qmc_help(stdout, "utils");` or `qmc-help utils`

## Functions (17)

[`linspace`](#linspace), [`logspace`](#logspace), [`range`](#range), [`normalize_wavefunction`](#normalize_wavefunction), [`compute_norm_squared`](#compute_norm_squared), [`expectation_value`](#expectation_value), [`expectation_position`](#expectation_position), [`expectation_position_squared`](#expectation_position_squared), [`position_to_momentum`](#position_to_momentum), [`expectation_momentum`](#expectation_momentum), [`mean`](#mean), [`variance`](#variance), [`std_dev`](#std_dev), [`cvector_from_matrix_column`](#cvector_from_matrix_column), [`save_wavefunction`](#save_wavefunction), [`save_potential`](#save_potential), [`save_eigenvalues`](#save_eigenvalues)

### `linspace`

```c
double *linspace(double start, double end, int num_points);
```

```text
Array/Grid utilities
```

### `logspace`

```c
double *logspace(double start, double end, int num_points);
```

_No description in the header yet._

### `range`

```c
int *range(int start, int end);
```

_No description in the header yet._

### `normalize_wavefunction`

```c
void normalize_wavefunction(cvector_t *psi);
```

```text
Wave function utilities
```

### `compute_norm_squared`

```c
double compute_norm_squared(const cvector_t *psi, double grid_dx);
```

_No description in the header yet._

### `expectation_value`

```c
double expectation_value(const cvector_t *psi, const cvector_t *op_psi,
                         double grid_dx);
```

```text
Expectation values
```

### `expectation_position`

```c
double expectation_position(const cvector_t *psi, const double *pos_grid,
                            double grid_dx);
```

_No description in the header yet._

### `expectation_position_squared`

```c
double expectation_position_squared(const cvector_t *psi,
                                    const double *pos_grid, double grid_dx);
```

_No description in the header yet._

### `position_to_momentum`

```c
cvector_t *position_to_momentum(const cvector_t *psi_x, double grid_dx);
```

```text
Momentum space (via FFT interface)
```

### `expectation_momentum`

```c
double expectation_momentum(const cvector_t *psi_k, const double *k_grid,
                            double grid_dk);
```

_No description in the header yet._

### `mean`

```c
double mean(const double *data, int n);
```

```text
Statistics
```

### `variance`

```c
double variance(const double *data, int n);
```

_No description in the header yet._

### `std_dev`

```c
double std_dev(const double *data, int n);
```

_No description in the header yet._

### `cvector_from_matrix_column`

```c
cvector_t *cvector_from_matrix_column(const cmatrix_t *m, int col);
```

```text
Matrix column extraction
```

### `save_wavefunction`

```c
void save_wavefunction(const char *filename, const double *pos_grid,
                       const cvector_t *psi, int num_points);
```

```text
File I/O utilities
```

### `save_potential`

```c
void save_potential(const char *filename, const double *x, const double *V,
                    int n);
```

_No description in the header yet._

### `save_eigenvalues`

```c
void save_eigenvalues(const char *filename, const double *eigenvals, int n);
```

_No description in the header yet._
