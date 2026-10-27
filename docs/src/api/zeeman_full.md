# `zeeman_full`

> Generated from `physics/zeeman_full.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Exact Zeeman effect for a p-electron across all field strengths (weak to Paschen-Back), with line intensities.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/zeeman_full.h>
```

**Start with:** [`zf_solve`](#zf_solve), [`zf_line_intensity`](#zf_line_intensity), [`zf_magnetic_moment`](#zf_magnetic_moment)

- Runtime help: `qmc_help(stdout, "zeeman_full");` or `qmc-help zeeman_full`

## Overview

```text
Zeeman effect of a p electron at any field: weak-field (anomalous Zeeman) to
Paschen-Back, with line intensities

 H = \xi L.S + b (L_z + 2 S_z) for l = 1, s = 1/2

Where:
 - \xi is spin-orbit constant
 - b = \mu_B B (both in same energy unit)

H conserves m_j = m_l + m_s, so it is block diagonal: two 1 x 1 blocks (m_j =
+-3/2) and two 2 x 2 blocks (m_j = +-1/2), solved in closed form. In  basis
(|1,dn>, |0,up>) of m_j = +1/2 block, H = [[-xi/2, xi / \sqrt(2)], [xi /
\sqrt(2), b]] and for m_j = -1/2, in basis (|0,dn>, |-1,up>), H = [[-b, xi /
\sqrt(2)], [xi / \sqrt(2), -xi/2]]

Limits:
 - b = 0 gives fine-structure levels \xi/2 (2P3/2, four-fold) and -\xi
   (2P1/2, two-fold)
 - b << \xi gives E = E_j + g_J m_j b (g = 4/3, 2/3)
 - b >> \xi gives Paschen-Back levels b (m_l + 2 m_s) + \xi m_l m_s

2S states are simply E = 2 m_s b = \pm b. An electric-dipole line from 2S
(spin m_s) into a P state with polarisation q (q = m_l of final state: -1, 0,
+1 for \sigma_-, \pi, \sigma_+) has relative intensity |<m_l = q, m_s | P
state>|^2: dipole operator does not act on spin. At b = 0, D2 (3/2) lines
carry twice weight of D1 (1/2) lines

Basis ordering of amplitudes: index (m_l + 1) * 2 + (1 if spin down else 0)
```

## Types

### `zf_state_t`

```c
typedef struct {
  double energy;
  int two_mj;    /* 2 m_j */
  int j_branch;  /* 3 for branch that is 2P3/2 at b = 0, 1 for 2P1/2
                    (labels 2j) */
  double amp[6]; /* amplitudes in |m_l, m_s> basis */
} zf_state_t;
```

## Functions (4)

[`zf_solve`](#zf_solve), [`zf_s_energy`](#zf_s_energy), [`zf_line_intensity`](#zf_line_intensity), [`zf_magnetic_moment`](#zf_magnetic_moment)

### `zf_solve`

```c
int zf_solve(double xi, double b, zf_state_t s[6]);
```

```text
6 2P states ordered (3/2, +3/2), (3/2, +1/2), (3/2, -1/2), (3/2, -3/2),
(1/2, +1/2), (1/2, -1/2); 0, or -1 for non-finite input
```

### `zf_s_energy`

```c
double zf_s_energy(double b, int ms_down);
```

```text
2S energy of spin projection ms_down = 0 (up) or 1 (down): -+... E = b for
up, -b for down
```

### `zf_line_intensity`

```c
double zf_line_intensity(const zf_state_t *s, int ms_down, int q);
```

```text
Relative intensity of line from 2S (spin ms_down) into state s with
polarisation q in {-1, 0, 1}: |amp(q, ms)|^2
```

### `zf_magnetic_moment`

```c
double zf_magnetic_moment(const zf_state_t *s);
```

```text
<L_z + 2 S_z> in a state (field derivative of its energy, per unit b)
```
