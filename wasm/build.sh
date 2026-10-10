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
#   docs/src/playground/qmc_anderson.wasm   Anderson-localisation demo
#   docs/src/playground/qmc_qwalk.wasm      quantum-walk demo
#   docs/src/playground/qmc_lz.wasm         Landau-Zener / Stueckelberg demo
#   docs/src/playground/qmc_rotor.wasm      kicked-rotor demo
#   docs/src/playground/qmc_ssh.wasm        SSH topological-chain demo
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
    # zig)
    #   shellcheck disable=SC2086
    #   zig cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
    #     -mexec-model=reactor -Wl,--no-entry -Wl,--max-memory=0x10000000 \
    #     -Wl,--strip-all -Icore -Icore/linalg -Iphysics "$@" -lm -o "$out"
    #   ;;
    # ziglang)
    #   # shellcheck disable=SC2086
    #   python3 -m ziglang cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
    #     -mexec-model=reactor -Wl,--no-entry -Wl,--max-memory=0x10000000 \
    #     -Wl,--strip-all -Icore -Icore/linalg -Iphysics "$@" -lm -o "$out"
    #   ;;
    zig)
      zig cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
        -mexec-model=reactor -Wl,--no-entry \
        -Wl,--strip-all -Icore -Icore/linalg -Iphysics "$@" -lm -o "$out"
      ;;
    ziglang)
      # shellcheck disable=SC2086
      python3 -m ziglang cc -target wasm32-wasi -O2 -fno-sanitize=undefined \
        -mexec-model=reactor -Wl,--no-entry \
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

build_module "$OUT/qmc_anderson.wasm" \
  "qmc_an_n qmc_an_dt qmc_an_time qmc_an_buffer qmc_an_density qmc_an_reset \
qmc_an_step qmc_an_norm qmc_an_width qmc_an_ipr qmc_an_energy qmc_an_xi" \
  wasm/qmc_anderson.c physics/anderson.c core/vector.c

build_module "$OUT/qmc_qwalk.wasm" \
  "qmc_qw_n qmc_qw_tmax qmc_qw_time qmc_qw_buffer qmc_qw_reset qmc_qw_step \
qmc_qw_norm qmc_qw_mean qmc_qw_variance qmc_qw_asymptote" \
  wasm/qmc_qwalk.c physics/quantum_walk.c core/vector.c

build_module "$OUT/qmc_lz.wasm" \
  "qmc_lz_vec qmc_lz_scan_buffer qmc_lz_scan_max qmc_lz_time qmc_lz_total_time qmc_lz_delta \
qmc_lz_upper qmc_lz_exact qmc_lz_reset qmc_lz_step qmc_lz_scan_rate qmc_lz_scan_amp" \
  wasm/qmc_lz.c physics/landau_zener.c

build_module "$OUT/qmc_rotor.wasm" \
  "qmc_kr_n qmc_kr_np qmc_kr_tmax qmc_kr_time qmc_kr_qbuffer qmc_kr_cbuffer qmc_kr_reset \
qmc_kr_step qmc_kr_qm2 qmc_kr_cm2 qmc_kr_norm qmc_kr_diffusion qmc_kr_edge" \
  wasm/qmc_rotor.c physics/kicked_rotor.c core/vector.c core/fft/fft.c

build_module "$OUT/qmc_ssh.wasm" \
  "qmc_ssh_cells qmc_ssh_sweep_max qmc_ssh_energies qmc_ssh_vectors qmc_ssh_bonds qmc_ssh_sweep_buffer \
qmc_ssh_set qmc_ssh_sweep qmc_ssh_winding qmc_ssh_zak qmc_ssh_xi qmc_ssh_gap qmc_ssh_bulk qmc_ssh_edge_weight" \
  wasm/qmc_ssh.c physics/ssh_chain.c core/matrix.c core/vector.c core/utils.c \
  core/fft/fft.c core/linalg/*.c core/ode/*.c
