#!/bin/sh
# Check that WebAssembly build of playground solver reproduces
# native C build bit-for-bit on key observables, check physics
#
#   make wasm-test        (or ./wasm/test/run.sh after ./wasm/build.sh)
set -eu
cd "$(dirname "$0")/../.."

TMP=${TMPDIR:-/tmp}
gcc -O2 -Wall -Wextra -I. wasm/test/native_ref.c wasm/qmc_wasm.c core/vector.c \
  core/fft/fft.c core/fft/fft2d.c core/fft/fft3d.c physics/soft.c -lm \
  -o "$TMP/qmc_native_ref"

"$TMP/qmc_native_ref" | grep -E '^norm' > "$TMP/qmc_native.txt"
node wasm/test/check.mjs| grep -E '^norm' > "$TMP/qmc_wasm.txt"

python3 - "$TMP/qmc_native.txt" "$TMP/qmc_wasm.txt" <<'PY'
import sys
native, wasm = [dict(l.split() for l in open(p)) for p in sys.argv[1:3]]
for key in native:
    n, w = float(native[key]), float(wasm[key])
    err = abs(n - w)
    if err > 1e-10:
        print(f"FAIL: {key}: native={n:.15g} wasm={w:.15g} diff={err:.3e}")
        sys.exit(1)
    print(f"ok   {key}: {n:.15g} (wasm {w:.15g}, diff {err:.2e})")
PY

echo "--- double-slit fringe check ---"
node wasm/test/fringes.mjs | grep -E 'peak|gaps|predicted'

echo "--- Hofstadter butterfly: band edges, TKNN vs library Chern numbers ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics -Icore/linalg \
  wasm/test/butterfly_check.c wasm/qmc_butterfly.c physics/tight_binding.c \
  core/matrix.c core/vector.c core/linalg/*.c -lm -o "$TMP/qmc_bf_check"
"$TMP/qmc_bf_check" > "$TMP/qmc_bf_native.txt"   # exits non-zero on any failure
cat "$TMP/qmc_bf_native.txt"
node wasm/test/butterfly_check.mjs| grep digest > "$TMP/qmc_bf_wasm.txt"
if grep digest "$TMP/qmc_bf_native.txt" | diff -u - "$TMP/qmc_bf_wasm.txt"; then
  echo "OK: butterfly WASM matches native"
else
  echo "FAIL: butterfly WASM and native results differ" >&2
  exit 1
fi

echo "--- Dirac: Klein step vs closed form, Zitterbewegung, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics wasm/test/dirac_ref.c wasm/qmc_dirac.c physics/dirac_evolve.c \
  core/vector.c core/fft/fft.c -lm -o "$TMP/qmc_dirac_ref"
"$TMP/qmc_dirac_ref" > "$TMP/qmc_dirac_native.txt"   # exits non-zero on physics failure
cat "$TMP/qmc_dirac_native.txt"
node wasm/test/dirac_check.mjs > "$TMP/qmc_dirac_wasm.txt"
if diff -u "$TMP/qmc_dirac_native.txt" "$TMP/qmc_dirac_wasm.txt"; then
  echo "OK: Dirac WASM matches native"
else
  echo "FAIL: Dirac WASM and native results differ" >&2
  exit 1
fi

echo "--- Hydrogen orbitals: native vs WASM sampling ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics -Icore/linalg wasm/test/orbital_ref.c wasm/qmc_orbital.c \
  physics/orbital_sample.c physics/hydrogen.c physics/central_potential.c physics/potentials.c \
  physics/wavefn.c core/special/*.c core/vector.c core/matrix.c core/utils.c core/fft/fft.c core/linalg/*.c core/ode/*.c \
  -lm -o "$TMP/qmc_orbital_ref"
"$TMP/qmc_orbital_ref" > "$TMP/qmc_orbital_native.txt"
cat "$TMP/qmc_orbital_native.txt"
node wasm/test/orbital_check.mjs > "$TMP/qmc_orbital_wasm.txt"
if diff -u "$TMP/qmc_orbital_native.txt" "$TMP/qmc_orbital_wasm.txt"; then
  echo "OK: orbital WASM matches native"
else
  echo "FAIL: orbital WASM and native results differ" >&2
  exit 1
fi

echo "--- Bloch sphere: Lindblad qubit, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics -Icore/linalg wasm/test/bloch_ref.c wasm/qmc_bloch.c \
  physics/bloch.c physics/lindblad.c physics/rabi.c core/matrix.c core/vector.c core/utils.c \
  core/fft/fft.c core/linalg/*.c core/ode/*.c -lm -o "$TMP/qmc_bloch_ref"
"$TMP/qmc_bloch_ref" > "$TMP/qmc_bloch_native.txt"
cat "$TMP/qmc_bloch_native.txt"
node wasm/test/bloch_check.mjs > "$TMP/qmc_bloch_wasm.txt"
if diff -u "$TMP/qmc_bloch_native.txt" "$TMP/qmc_bloch_wasm.txt"; then
  echo "OK: Bloch WASM matches native"
else
  echo "FAIL: Bloch WASM and native results differ" >&2
  exit 1
fi

echo "--- Anderson localisation: tight-binding chain, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics wasm/test/anderson_ref.c wasm/qmc_anderson.c physics/anderson.c \
  core/vector.c -lm -o "$TMP/qmc_anderson_ref"
"$TMP/qmc_anderson_ref" > "$TMP/qmc_anderson_native.txt"
cat "$TMP/qmc_anderson_native.txt"
node wasm/test/anderson_check.mjs > "$TMP/qmc_anderson_wasm.txt"
if diff -u "$TMP/qmc_anderson_native.txt" "$TMP/qmc_anderson_wasm.txt"; then
  echo "OK: Anderson WASM matches native"
else
  echo "FAIL: Anderson WASM and native results differ" >&2
  exit 1
fi

echo "--- Quantum walk: coin + shift, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics wasm/test/qwalk_ref.c wasm/qmc_qwalk.c physics/quantum_walk.c \
  core/vector.c -lm -o "$TMP/qmc_qwalk_ref"
"$TMP/qmc_qwalk_ref" > "$TMP/qmc_qwalk_native.txt"
cat "$TMP/qmc_qwalk_native.txt"
node wasm/test/qwalk_check.mjs > "$TMP/qmc_qwalk_wasm.txt"
if diff -u "$TMP/qmc_qwalk_native.txt" "$TMP/qmc_qwalk_wasm.txt"; then
  echo "OK: quantum-walk WASM matches native"
else
  echo "FAIL: quantum-walk WASM and native results differ" >&2
  exit 1
fi

echo "--- Landau-Zener: swept two-level system, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics wasm/test/lz_ref.c wasm/qmc_lz.c physics/landau_zener.c \
  -lm -o "$TMP/qmc_lz_ref"
"$TMP/qmc_lz_ref" > "$TMP/qmc_lz_native.txt"
cat "$TMP/qmc_lz_native.txt"
node wasm/test/lz_check.mjs > "$TMP/qmc_lz_wasm.txt"
if diff -u "$TMP/qmc_lz_native.txt" "$TMP/qmc_lz_wasm.txt"; then
  echo "OK: Landau-Zener WASM matches native"
else
  echo "FAIL: Landau-Zener WASM and native results differ" >&2
  exit 1
fi

echo "--- Kicked rotor: quantum map + standard map, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics wasm/test/rotor_ref.c wasm/qmc_rotor.c physics/kicked_rotor.c \
  core/vector.c core/fft/fft.c -lm -o "$TMP/qmc_rotor_ref"
"$TMP/qmc_rotor_ref" > "$TMP/qmc_rotor_native.txt"
cat "$TMP/qmc_rotor_native.txt"
node wasm/test/rotor_check.mjs > "$TMP/qmc_rotor_wasm.txt"
if diff -u "$TMP/qmc_rotor_native.txt" "$TMP/qmc_rotor_wasm.txt"; then
  echo "OK: kicked-rotor WASM matches native"
else
  echo "FAIL: kicked-rotor WASM and native results differ" >&2
  exit 1
fi

echo "--- SSH chain: edge states, winding and Zak phase, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics -Icore/linalg wasm/test/ssh_ref.c wasm/qmc_ssh.c physics/ssh_chain.c \
  core/matrix.c core/vector.c core/utils.c core/fft/fft.c core/linalg/*.c core/ode/*.c -lm -o "$TMP/qmc_ssh_ref"
"$TMP/qmc_ssh_ref" > "$TMP/qmc_ssh_native.txt"
cat "$TMP/qmc_ssh_native.txt"
node wasm/test/ssh_check.mjs > "$TMP/qmc_ssh_wasm.txt"
if diff -u "$TMP/qmc_ssh_native.txt" "$TMP/qmc_ssh_wasm.txt"; then
  echo "OK: SSH-chain WASM matches native"
else
  echo "FAIL: SSH-chain WASM and native results differ" >&2
  exit 1
fi

echo "--- Ising quench: free-fermion correlations and Loschmidt rate, native vs WASM ---"
gcc -O2 -Wall -Wextra -I. -Icore -Iphysics wasm/test/tfim_ref.c wasm/qmc_tfim.c physics/tfim_quench.c \
  -lm -o "$TMP/qmc_tfim_ref"
"$TMP/qmc_tfim_ref" > "$TMP/qmc_tfim_native.txt"
cat "$TMP/qmc_tfim_native.txt"
node wasm/test/tfim_check.mjs > "$TMP/qmc_tfim_wasm.txt"
if diff -u "$TMP/qmc_tfim_native.txt" "$TMP/qmc_tfim_wasm.txt"; then
  echo "OK: Ising-quench WASM matches native"
else
  echo "FAIL: Ising-quench WASM and native results differ" >&2
  exit 1
fi
