#!/usr/bin/env python3
"""Demonstrates the bound subset of QMC's C API from Python

USAGE: Run from the repository root after building the shared library:
    make shared SANITIZE=0 USE_LAPACK=0
    python3 bindings/python/example.py
"""

import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])

from qmc import dmrg_run, hydrogen_energy_level  # noqa: E402


def main():
    print("hydrogen_energy_level:")
    for n in (1, 2, 3, 4):
        e_joules = hydrogen_energy_level(n)
        e_ev = e_joules / 1.602176634e-19
        print(f"  n={n}: E = {e_joules:.6e} J  ({e_ev:.4f} eV)")

    print("\ndmrg_run (Heisenberg AFM chain, N=40, m=20):")
    result = dmrg_run(n_target=40, jz=1.0, jxy=1.0, m_max=20)
    print(f"  {result}")
    print(f"  energy per site: {result.energy_per_site:.8f}")


if __name__ == "__main__":
    main()
