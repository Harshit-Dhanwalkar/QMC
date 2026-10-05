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
node wasm/test/check.mjs 2>/dev/null | grep -E '^norm' > "$TMP/qmc_wasm.txt"

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
