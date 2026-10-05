#!/bin/sh
# Check that the WebAssembly build of the playground solver reproduces the
# native C build bit-for-bit on the key observables, then check the physics.
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

if diff -u "$TMP/qmc_native.txt" "$TMP/qmc_wasm.txt"; then
  echo "OK: WASM and native norms agree to 12 digits"
else
  echo "FAIL: WASM and native results differ" >&2
  exit 1
fi

echo "--- double-slit fringe check ---"
node wasm/test/fringes.mjs | grep -E 'peak|gaps|predicted'
