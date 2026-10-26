# `ccsd_t`

> Generated from `physics/ccsd_t.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Perturbative triples correction CCSD(T) on top of CCSD.

**Category:** Quantum chemistry (molecules)

```c
#include <qmc/physics/ccsd_t.h>
```

**Start with:** [`ccsdt_run`](#ccsdt_run)

- Example: `examples/eg_51_ccsd_t.c`
- Test (every call with expected results): `tests/test_ccsd_t.c`
- Runtime help: `qmc_help(stdout, "ccsd_t");` or `qmc-help ccsd_t`

## Overview

```text
Perturbative triples correction, CCSD(T) (Reference: Raghavachari, Trucks,
Pople & Head-Gordon, Chem. Phys. Lett. 157, 479 (1989)), converged T1/T2
amplitudes (via ccsd_run_ex)

Spin-orbital formulas (Crawford & Schaefer's "An Introduction to Coupled
Cluster Theory", project #6):

  t3d_ijkabc = P(i / jk) P(a / bc) [ t1_ia * <jk||bc> ]
  t3c_ijkabc = P(i / jk) P(a / bc) [ \sum_e t2_jkae * <ei||bc>
                                   - \sum_m t2_imbc * <majk> ]
  E_(T) = (1/36) * \sum_ijkabc t3c_ijkabc * (t3c_ijkabc + t3d_ijkabc)
                                             / D_ijkabc

Where :
  P(p/qr) f(p,q,r) = f(p,q,r) - f(q,p,r) - f(r,q,p)
3-term cyclic antisymmetrizer over index being singled out) and
D_ijkabc is usual triples energy denominator
f_ii + f_jj + f_kk - f_aa - f_bb - f_cc.
```

## Types

### `ccsdt_result_t`

```c
typedef struct {
  double ccsd_correlation_energy; /* E_CCSD - E_RHF, (same as ccsd_result_t) */
  double perturbative_correction; /* (T) term */
  double total_energy;            /* E_RHF + CCSD corr + (T) */
  int ccsd_converged;
  int ccsd_iterations;
} ccsdt_result_t;
```

## Functions (1)

[`ccsdt_run`](#ccsdt_run)

### `ccsdt_run`

```c
ccsdt_result_t *ccsdt_run(int n_spatial, const double *h_mo,
                          const double *eri_mo, const double *mo_energy,
                          int n_electrons, int n_frozen_spatial, double e_rhf,
                          double conv_tol, int max_iter);
```

```text
Same inputs as ccsd_run : runs CCSD to convergence, then adds perturbative
(T) correction using converged amplitudes

Returns NULL on same invalid-input conditions as ccsd_run, or if underlying
CCSD does not converge within max_iter

 HACK: ccsdt_result_t's ccsd_converged field lets caller distinguish "didn't
converge" from "invalid input" if a non-NULL-on-non-convergence variant is
ever needed : for now this returns NULL in both cases, matching ccsd_run's
contract of not returning results from an unconverged calculation
```
