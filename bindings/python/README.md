# QMC Python Bindings (ctypes)

A small, scoped `ctypes` wrapper around QMC's compiled shared library:
**not** a full Python API for the library. It binds two representative
functions to demonstrate the pattern (load the `.so`, declare
`argtypes`/`restype`, mirror any returned struct's exact layout, wrap
with a Pythonic function) so it's straightforward to extend to more
of `physics/*.h` following the same shape:

- `hydrogen_energy_level(n)` - simple scalar in/out (`physics/hydrogen.h`)
- `dmrg_run(n_target, jz, jxy, m_max)` - returns a struct pointer, freed
  automatically (`physics/dmrg.h`)

## Prerequisites

Build the shared library first, from the repository root:

```sh
make shared SANITIZE=0 USE_LAPACK=0
```

`SANITIZE=0` matters here specifically: a plain `python3` process does not
have AddressSanitizer preloaded, so loading an ASan-instrumented `.so` via
`ctypes.CDLL` will fail or abort at import time.
( `USE_LAPACK` is choice; nothing currently bound depends on it.)

## Usage

```sh
python3 bindings/python/example.py
```

Or from script:

```python
from qmc import hydrogen_energy_level, dmrg_run

print(hydrogen_energy_level(1))  # -2.178960e-18 (Joules; -13.6 eV)

result = dmrg_run(n_target=40, jz=1.0, jxy=1.0, m_max=20)
print(result.energy_per_site)
```

By default, the wrapper looks for `build/libqmc.so` relative to the
repository root. If you built or installed elsewhere (e.g. after
`make install`), point it at the right file with an environment variable:

```sh
QMC_LIBRARY=/usr/local/lib/libqmc.so python3 bindings/python/example.py
```

---

# TODO
## Extending to more functions

Same three steps for any additional `physics/*.h` function:

1. Declare `_lib.<function_name>.argtypes` and `.restype` to match the C signature exactly.
2. If it returns a struct pointer, mirror the struct's fields in same order and with the same C types as a `ctypes.Structure` -- `ctypes` does not check this against the header.
3. If the C API says to `free()` the result, do so (see `dmrg_run`'s wrapper for the pattern).
