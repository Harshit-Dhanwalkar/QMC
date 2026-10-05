#!/bin/sh
# Build the QMC WebAssembly playground module.
#
#   ./wasm/build.sh            # uses emcc if installed, otherwise zig cc
#
# Output: docs/src/playground/qmc.wasm (picked up by the mdBook docs deploy).
#
# Toolchains (either one works):
#   * Emscripten:  https://emscripten.org/docs/getting_started/downloads.html
#   * Zig (no install needed, just pip):  pip install ziglang
set -eu

cd "$(dirname "$0")/.."
OUT=docs/src/playground
mkdir -p "$OUT"

SRCS="wasm/qmc_wasm.c core/vector.c core/fft/fft.c core/fft/fft2d.c \
core/fft/fft3d.c physics/soft.c"

EXPORTS="qmc_init qmc_n qmc_box qmc_potential qmc_clear_potential \
qmc_fill_rect qmc_build_slits qmc_wavepacket qmc_step qmc_norm qmc_density"

if command -v emcc >/dev/null 2>&1; then
  echo "[wasm] using emcc"
  EXPORT_LIST=$(printf '"_%s",' $EXPORTS | sed 's/,$//')
  # shellcheck disable=SC2086
  emcc -O3 $SRCS -o "$OUT/qmc.wasm" \
    -sSTANDALONE_WASM=1 -sALLOW_MEMORY_GROWTH=1 --no-entry \
    -sEXPORTED_FUNCTIONS="[$EXPORT_LIST]"
elif command -v zig >/dev/null 2>&1; then
  echo "[wasm] using zig"
  # shellcheck disable=SC2086
  zig cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
    -mexec-model=reactor -Wl,--no-entry -Wl,--strip-all \
    $SRCS -lm -o "$OUT/qmc.wasm"
elif python3 -c 'import ziglang' 2>/dev/null; then
  echo "[wasm] using python -m ziglang"
  # shellcheck disable=SC2086
  python3 -m ziglang cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
    -mexec-model=reactor -Wl,--no-entry -Wl,--strip-all \
    $SRCS -lm -o "$OUT/qmc.wasm"
else
  echo "[wasm] ERROR: no suitable toolchain found." >&2
  echo "Install one of:" >&2
  echo "  - emcc (Emscripten):  https://emscripten.org/docs/getting_started/downloads.html" >&2
  echo "  - zig:                https://ziglang.org/download/" >&2
  echo "  - python3 -m pip install ziglang" >&2
  exit 1
fi

ls -l "$OUT/qmc.wasm"
