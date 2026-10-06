#!/bin/sh
# Build the QMC WebAssembly modules used by the browser demos.
#
#   ./wasm/build.sh            # uses emcc if installed, otherwise zig cc
#
# Output:
#   docs/src/playground/qmc.wasm            2D double-slit playground
#   docs/src/playground/qmc_butterfly.wasm  Hofstadter butterfly explorer
#   docs/src/playground/qmc_dirac.wasm      Dirac Klein-paradox demo
#   docs/src/playground/qmc_orbital.wasm    hydrogen orbital viewer
#   docs/src/playground/qmc_bloch.wasm      Bloch-sphere qubit demo
#
# Toolchains:
#   * Emscripten:  https://emscripten.org/docs/getting_started/downloads.html
#   * Zig:         https://ziglang.org/download/
#   * Zig via pip: pip install ziglang
set -eu

cd "$(dirname "$0")/.."
OUT=docs/src/playground
mkdir -p "$OUT"

# Pick a toolchain once.
if command -v emcc >/dev/null 2>&1; then
  echo "[wasm] using emcc"
  TOOL=emcc
elif command -v zig >/dev/null 2>&1; then
  echo "[wasm] using zig"
  TOOL=zig
elif python3 -c 'import ziglang' 2>/dev/null; then
  echo "[wasm] using python -m ziglang"
  TOOL=ziglang
else
  echo "[wasm] ERROR: no suitable toolchain found." >&2
  echo "Install one of:" >&2
  echo "  - emcc (Emscripten):  https://emscripten.org/docs/getting_started/downloads.html" >&2
  echo "  - zig:                https://ziglang.org/download/" >&2
  echo "  - python3 -m pip install ziglang" >&2
  exit 1
fi

# build_module <output.wasm> "<exported functions>" <source files...>
build_module() {
  out=$1; exports=$2; shift 2
  case "$TOOL" in
    emcc)
      export_list=$(printf '"_%s",' $exports | sed 's/,$//')
      # shellcheck disable=SC2086
      emcc -O3 "$@" -o "$out" -Icore -Icore/linalg -Iphysics \
        -sSTANDALONE_WASM=1 -sALLOW_MEMORY_GROWTH=1 --no-entry \
        -sEXPORTED_FUNCTIONS="[$export_list]"
      ;;
    zig)
      # shellcheck disable=SC2086
      zig cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
        -mexec-model=reactor -Wl,--no-entry -Wl,--max-memory=268435456 \
        -Wl,--strip-all -Icore -Icore/linalg -Iphysics "$@" -lm -o "$out"
      ;;
    ziglang)
      # shellcheck disable=SC2086
      python3 -m ziglang cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
        -mexec-model=reactor -Wl,--no-entry -Wl,--max-memory=268435456 \
        -Wl,--strip-all -Icore -Icore/linalg -Iphysics "$@" -lm -o "$out"
      ;;
  esac
  ls -l "$out"
}

build_module "$OUT/qmc.wasm" \
  "qmc_init qmc_n qmc_box qmc_potential qmc_clear_potential qmc_fill_rect \
qmc_build_slits qmc_wavepacket qmc_step qmc_norm qmc_density" \
  wasm/qmc_wasm.c core/vector.c core/fft/fft.c core/fft/fft2d.c \
  core/fft/fft3d.c physics/soft.c

build_module "$OUT/qmc_butterfly.wasm" \
  "qmc_bf_buffer qmc_bf_max_q qmc_bf_edges qmc_bf_tknn qmc_bf_chern" \
  wasm/qmc_butterfly.c physics/tight_binding.c core/matrix.c core/vector.c \
  core/linalg/*.c

build_module "$OUT/qmc_dirac.wasm" \
  "qmc_dc_n qmc_dc_dx qmc_dc_dt qmc_dc_time qmc_dc_buffer qmc_dc_density \
qmc_dc_klein qmc_dc_zitter qmc_dc_step qmc_dc_norm qmc_dc_position \
qmc_dc_right qmc_dc_exact" \
  wasm/qmc_dirac.c physics/dirac_evolve.c core/vector.c core/fft/fft.c

build_module "$OUT/qmc_orbital.wasm" \
  "qmc_orb_max_points qmc_orb_buffer qmc_orb_sample qmc_orb_energy_ev qmc_orb_max_n" \
  wasm/qmc_orbital.c physics/orbital_sample.c physics/hydrogen.c \
  physics/central_potential.c physics/potentials.c physics/wavefn.c \
  core/special/*.c core/vector.c core/matrix.c core/utils.c core/fft/fft.c core/linalg/*.c core/ode/*.c

build_module "$OUT/qmc_bloch.wasm" \
  "qmc_bl_dt qmc_bl_vec qmc_bl_time qmc_bl_reset qmc_bl_step qmc_bl_purity qmc_bl_rabi" \
  wasm/qmc_bloch.c physics/bloch.c physics/lindblad.c physics/rabi.c \
  core/matrix.c core/vector.c core/utils.c core/fft/fft.c core/linalg/*.c core/ode/*.c
