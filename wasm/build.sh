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
#   docs/src/playground/qmc_tfim.wasm       Ising-quench demo
#   docs/src/playground/qmc_haldane.wasm    Haldane Chern-insulator demo
#   docs/src/playground/qmc_wigner.wasm     Wigner-function / Kerr cat-state demo
#   docs/src/playground/qmc_grover.wasm     Grover-search demo
#   docs/src/playground/qmc_bell.wasm       CHSH Bell-test on a Werner state demo
#   docs/src/playground/qmc_carpet.wasm     quantum-carpet / revival demo
#   docs/src/playground/qmc_zeno.wasm       quantum-Zeno demo
#   docs/src/playground/qmc_corral.wasm     quantum-corral (circular well) demo
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
qmc_ssh_set qmc_ssh_sweep qmc_ssh_winding qmc_ssh_zak qmc_ssh_xi qmc_ssh_gap qmc_ssh_bulk \
qmc_ssh_edge_weight" \
  wasm/qmc_ssh.c physics/ssh_chain.c core/matrix.c core/vector.c core/utils.c \
  core/fft/fft.c core/linalg/*.c core/ode/*.c

build_module "$OUT/qmc_tfim.wasm" \
  "qmc_tq_rmax qmc_tq_curve_max qmc_tq_modes qmc_tq_zz_buffer qmc_tq_rate_buffer qmc_tq_mx_buffer \
qmc_tq_density_buffer qmc_tq_set qmc_tq_eval qmc_tq_last_mx qmc_tq_last_rate qmc_tq_curve qmc_tq_tcrit \
qmc_tq_front_speed" \
  wasm/qmc_tfim.c physics/tfim_quench.c

build_module "$OUT/qmc_haldane.wasm" \
  "qmc_hd_width qmc_hd_states qmc_hd_nk_max qmc_hd_energy_buffer qmc_hd_mean_buffer qmc_hd_prob_buffer \
qmc_hd_set qmc_hd_compute qmc_hd_profile qmc_hd_chern qmc_hd_limit qmc_hd_valley_gap qmc_hd_vmax \
qmc_hd_cmin" \
  wasm/qmc_haldane.c physics/tight_binding.c core/matrix.c core/vector.c core/utils.c \
  core/fft/fft.c core/linalg/*.c core/ode/*.c

build_module "$OUT/qmc_wigner.wasm" \
  "qmc_wg_n qmc_wg_nmax qmc_wg_extent qmc_wg_w_buffer qmc_wg_x_buffer qmc_wg_p_buffer qmc_wg_pn_buffer \
qmc_wg_nbar qmc_wg_negativity qmc_wg_wmin qmc_wg_wmax qmc_wg_state" \
  wasm/qmc_wigner.c physics/phase_space.c

build_module "$OUT/qmc_grover.wasm" \
  "qmc_gv_max_n qmc_gv_curve_max qmc_gv_amp_buffer qmc_gv_marked_buffer qmc_gv_curve_buffer \
  qmc_gv_exact_buffer qmc_gv_n qmc_gv_iter qmc_gv_marks qmc_gv_set_n qmc_gv_mark qmc_gv_reset \
  qmc_gv_step qmc_gv_goto qmc_gv_success qmc_gv_theta qmc_gv_kopt qmc_gv_exact qmc_gv_curve" \
  wasm/qmc_grover.c physics/grover_multi.c

build_module "$OUT/qmc_bell.wasm" \
  "qmc_bl_curve_max qmc_bl_quantum_buffer qmc_bl_local_buffer qmc_bl_set qmc_bl_chsh qmc_bl_local_chsh \
  qmc_bl_max qmc_bl_e qmc_bl_e_local qmc_bl_curve qmc_bl_sample qmc_bl_est qmc_bl_sigma" \
  wasm/qmc_bell.c physics/bell_chsh.c

build_module "$OUT/qmc_carpet.wasm" \
  "qmc_qc_nx qmc_qc_nt qmc_qc_curve_len qmc_qc_carpet_buffer qmc_qc_slice_buffer qmc_qc_curve_buffer \
qmc_qc_weight qmc_qc_energy qmc_qc_set qmc_qc_slice" \
  wasm/qmc_carpet.c physics/box_revival.c

build_module "$OUT/qmc_zeno.wasm" \
  "qmc_zn_curve_len qmc_zn_max_n qmc_zn_time_buffer qmc_zn_n_buffer qmc_zn_rate_buffer \
qmc_zn_cont_buffer qmc_zn_alive_buffer qmc_zn_total_time qmc_zn_pulsed qmc_zn_pulsed_curves \ 
qmc_zn_continuous qmc_zn_rates qmc_zn_rate qmc_zn_sample" \
  wasm/qmc_zeno.c physics/zeno.c

build_module "$OUT/qmc_corral.wasm" \
  "qmc_cb_n qmc_cb_mode_count qmc_cb_field_buffer qmc_cb_mode_m qmc_cb_mode_n qmc_cb_mode_sine qmc_cb_mode_k \
qmc_cb_show_mode qmc_cb_launch qmc_cb_energy qmc_cb_evolve" \
  wasm/qmc_corral.c physics/circular_billiard.c
