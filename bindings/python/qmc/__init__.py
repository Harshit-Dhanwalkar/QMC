# WARN: This is not a full Python wrapper for the library 
"""
Bound:
  - hydrogen_energy_level (physics/hydrogen.h)   - simple scalar in/out
  - dmrg_run              (physics/dmrg.h)       - returns a struct pointer

Prerequisites: build the shared library first:
    make shared SANITIZE=0 USE_LAPACK=0
"""

import ctypes
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
_DEFAULT_LIB_PATH = os.path.join(_REPO_ROOT, "build", "libqmc.so")
_LIB_PATH = os.environ.get("QMC_LIBRARY", _DEFAULT_LIB_PATH)

try:
    _lib = ctypes.CDLL(_LIB_PATH)
except OSError as exc:
    raise OSError(
        f"Could not load {_LIB_PATH!r}. Build it first with "
        "'make shared SANITIZE=0 USE_LAPACK=0', or set the QMC_LIBRARY "
        "environment variable to an existing libqmc.so (e.g. after "
        "'make install', to $PREFIX/lib/libqmc.so)."
    ) from exc


# hydrogen_energy_level(int n) -> double  (physics/hydrogen.h)
_lib.hydrogen_energy_level.argtypes = [ctypes.c_int]
_lib.hydrogen_energy_level.restype = ctypes.c_double


def hydrogen_energy_level(n: int) -> float:
    """Analytic hydrogen energy level E_n, in Joules (SI units).

    Mirrors physics/hydrogen.h's hydrogen_energy_level exactly -- see that
    header for the underlying physics.
    """
    return _lib.hydrogen_energy_level(n)


# dmrg_run(...) -> dmrg_result_t*  (physics/dmrg.h)
class _DmrgResult(ctypes.Structure):
    # Field order/types must match physics/dmrg.h's dmrg_result_t exactly --
    # ctypes does not check this against the C header for you.
    _fields_ = [
        ("energy", ctypes.c_double),
        ("energy_per_site", ctypes.c_double),
        ("N_reached", ctypes.c_int),
        ("truncation_dim", ctypes.c_int),
        ("truncation_error", ctypes.c_double),
    ]


_lib.dmrg_run.argtypes = [
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_int,
]
_lib.dmrg_run.restype = ctypes.POINTER(_DmrgResult)

_libc = ctypes.CDLL(None)
_libc.free.argtypes = [ctypes.c_void_p]
_libc.free.restype = None


class DmrgResult:
    """Plain-Python copy of dmrg_result_t's fields (safe to keep around
    after the underlying C struct is freed)."""

    def __init__(self, c_result):
        self.energy = c_result.energy
        self.energy_per_site = c_result.energy_per_site
        self.n_reached = c_result.N_reached
        self.truncation_dim = c_result.truncation_dim
        self.truncation_error = c_result.truncation_error

    def __repr__(self):
        return (
            f"DmrgResult(energy={self.energy!r}, "
            f"energy_per_site={self.energy_per_site!r}, "
            f"n_reached={self.n_reached!r}, "
            f"truncation_dim={self.truncation_dim!r}, "
            f"truncation_error={self.truncation_error!r})"
        )


def dmrg_run(n_target: int, jz: float, jxy: float, m_max: int) -> DmrgResult:
    """Infinite-system DMRG ground state of the open spin-1/2 XXZ chain.

    Mirrors physics/dmrg.h's dmrg_run. Raises RuntimeError for the same
    invalid-input cases the C function returns NULL for (n_target < 2,
    m_max < 1, non-finite jz/jxy) or on allocation failure.
    """
    ptr = _lib.dmrg_run(n_target, jz, jxy, m_max)
    if not ptr:
        raise RuntimeError(
            "dmrg_run returned NULL (invalid input or allocation failure) "
            "-- see physics/dmrg.h for the exact validity conditions"
        )
    try:
        return DmrgResult(ptr.contents)
    finally:
        _libc.free(ptr)


__all__ = ["hydrogen_energy_level", "dmrg_run", "DmrgResult"]
