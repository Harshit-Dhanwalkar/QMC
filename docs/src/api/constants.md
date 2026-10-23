# `constants`

> Generated from `core/constants.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Fundamental physical constants in SI units (HBAR, M_ELECTRON, E_CHARGE, EPSILON_0, ...).

**Category:** Core: data types & utilities

```c
#include <qmc/core/constants.h>
```

- Runtime help: `qmc_help(stdout, "constants");` or `qmc-help constants`

## Macros and constants

```c
#define HBAR 1.054571817e-34        /* Planck constant / 2π (J·s) */
```

```text
Fundamental constants (SI units)
```

```c
#define M_ELECTRON 9.1093837015e-31 /* Electron mass (kg) */
```

```c
#define M_PROTON 1.67262192369e-27  /* Proton mass (kg) */
```

```c
#define E_CHARGE 1.602176634e-19    /* Elementary charge (C) */
```

```c
#define K_COULOMB 8.9875517923e9    /* Coulomb constant (N·m²/C²) */
```

```c
#define EPSILON_0 8.8541878128e-12  /* Vacuum permittivity (F/m) */
```

```c
#define C_LIGHT 299792458.0         /* Speed of light (m/s) */
```

```c
#define AU_LENGTH 5.29177210903e-11   /* Bohr radius (m) */
```

```text
Atomic units (useful for QM calculations)
```

```c
#define AU_ENERGY 4.3597447222071e-18 /* Hartree (J) */
```

```c
#define AU_TIME 2.4188843265857e-17   /* Atomic time unit (s) */
```

```c
#define HBAR_SQ (HBAR * HBAR)
```

```text
Convenient combinations
```

```c
#define HBAR_2M (HBAR_SQ / (2.0 * M_ELECTRON))
```
