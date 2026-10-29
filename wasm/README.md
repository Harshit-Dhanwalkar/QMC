# QMC in browser (WebAssembly)

`qmc_wasm.c` is a thin wrapper around the library's own 2D split-operator solver (`physics/soft.c` -> `core/fft/fft2d.c`). The browser page in `docs/src/playground/index.html` calls it every animation frame.

`qmc_butterfly.c` does the same for the tight-binding module: the Hofstadter butterfly page (`docs/src/playground/butterfly.html`) diagonalises `tb_model_hofstadter` for every flux p/q up to q = 50, labels each gap with its Chern number (TKNN equation) and, on click, recomputes that number with the library's Fukui-Hatsugai-Suzuki `tb_chern_number`.

`qmc_dirac.c` wraps the Dirac time-evolution module (`physics/dirac_evolve.c`): `docs/src/playground/dirac.html` fires a positive-energy packet at a sharp step (Klein paradox, compared live with `dirac_step_transmission`) or shows Zitterbewegung of a packet built from both energy signs.

`qmc_orbital.c` wraps the orbital sampler (`physics/orbital_sample.c`): `docs/src/playground/orbitals.html` draws any hydrogen orbital (n up to 8, real or complex) as a rotatable 3D point cloud sampled from |psi_nlm|^2.

`qmc_bloch.c` wraps the Bloch-vector module (`physics/bloch.c`, RK4 Lindblad evolution): `docs/src/playground/bloch.html` drives and damps a qubit and plots the state on the Bloch sphere next to the closed-form Rabi curve.

`qmc_anderson.c` wraps the Anderson-localisation module (`physics/anderson.c`): `docs/src/playground/anderson.html` evolves a particle on a disordered tight-binding chain and shows it spreading, then freezing.

`qmc_qwalk.c` wraps the quantum walk (`physics/quantum_walk.c`): `docs/src/playground/quantum_walk.html` runs a coin-and-shift walk next to the classical random walk.

`qmc_lz.c` wraps the Landau-Zener integrator (`physics/landau_zener.c`): `docs/src/playground/landau_zener.html` sweeps a qubit through an avoided crossing and compares the jump probability with exp(-pi Omega^2 / 2 rate).

`qmc_rotor.c` wraps the kicked rotor (`physics/kicked_rotor.c`): `docs/src/playground/kicked_rotor.html` runs the quantum map next to the classical standard map and shows dynamical localisation.

`qmc_ssh.c` wraps the Su-Schrieffer-Heeger chain (`physics/ssh_chain.c`): `docs/src/playground/ssh_chain.html` diagonalises a 20-cell chain while you tune the hopping and shows the topological edge states, the winding number and the Zak phase.

`qmc_tfim.c` wraps the exact free-fermion quench of the transverse-field Ising chain (`physics/tfim_quench.c`): `docs/src/playground/ising_quench.html` shows the correlation light cone, the transverse magnetisation and the Loschmidt rate function with its dynamical-phase-transition kinks.

`qmc_haldane.c` wraps the Haldane model and its ribbon solver (`physics/tight_binding.c`): `docs/src/playground/haldane.html` shows the chiral edge states, the Chern number and the phase diagram. This module imports 3 WASI calls (`fd_write`, `fd_seek`, `fd_close`) because the complex eigensolver carries stdio error paths; none of them runs on valid input and the page passes stubs.

`qmc_wigner.c` wraps the Fock-basis states, the Kerr evolution and the Wigner function (`physics/phase_space.c`): `docs/src/playground/wigner.html` shows a coherent, Fock, cat or squeezed state in phase space and lets you apply the Kerr phase.

`qmc_grover.c` wraps the multi-target Grover search (`physics/grover_multi.c`): `docs/src/playground/grover.html` lets you mark items in a database of up to 1024 and step the amplitude amplification.

`qmc_bell.c` wraps the CHSH Bell test (`physics/bell_chsh.c`): `docs/src/playground/bell.html` compares a noisy singlet with a local-hidden-variable model.

`qmc_zeno.c` wraps the quantum Zeno effect (`physics/zeno.c`): `docs/src/playground/zeno.html` compares pulsed and continuous measurement of a driven qubit.

`qmc_corral.c` wraps the circular-well eigenmodes (`physics/circular_billiard.c`): `docs/src/playground/corral.html` shows the Bessel modes and a packet bouncing in the disc.

`qmc_ring.c` wraps the tight-binding Aharonov-Bohm ring (`physics/ab_ring.c`): `docs/src/playground/ring.html` shows the flux-dependent spectrum, persistent current and two-terminal transmission.

`qmc_tunnel.c` wraps the transfer-matrix double barrier (`physics/double_barrier.c`): `docs/src/playground/tunnel.html` shows resonant transmission, the dwell time and the wavefunction in the well.

`qmc_josephson.c` wraps the RCSJ Josephson junction (`physics/josephson.c`): `docs/src/playground/josephson.html` animates the phase in the washboard potential and sweeps the I-V curve with hysteresis and Shapiro steps.

`qmc_dwell.c` wraps the double-well eigen-solver (`physics/double_well.c`): `docs/src/playground/doublewell.html` shows the parity doublet, the exponentially small splitting and a packet tunnelling between the wells.

`qmc_zeeman.c` wraps the full-field Zeeman solver (`physics/zeeman_full.c`): `docs/src/playground/zeeman.html` shows the 2P levels, the 2S to 2P lines with polarisation and intensity, from the anomalous Zeeman regime to Paschen-Back.

## Build

```sh
# Option A: Zig as the C compiler, no system install needed
pip install ziglang
./wasm/build.sh            # or: make wasm

# Option B: Emscripten (https://emscripten.org), picked automatically if present
./wasm/build.sh
```

Output: `docs/src/playground/qmc.wasm`. The docs workflow (mdBook) copies the whole `playground/` folder to
`https://harshit-dhanwalkar.github.io/QMC/playground/`.

## Run locally

Browsers refuse to `fetch()` a `.wasm` file from `file://`, so serve the folder:

```sh
cd docs/src/playground && python3 -m http.server 8000
# open http://localhost:8000/
```

## Test

```sh
make wasm-test    # native vs WASM agreement + double-slit fringe spacing
```

## Regenerate README GIF

```sh
pip install pillow numpy
node wasm/tools/dump_frames.mjs /tmp/qmc_frames.bin
python3 wasm/tools/make_gif.py /tmp/qmc_frames.bin docs/src/playground/preview.gif
```

## JavaScript API (exports of `qmc.wasm`)

| function                                                | purpose                                                        |
| ------------------------------------------------------- | -------------------------------------------------------------- |
| `qmc_init(n, box, dt)`                                  | allocate an `n x n` grid (n a power of two, <= 256)            |
| `qmc_build_slits(wall_x, thick, slits, sep, width, v0)` | vertical wall with 1 or 2 slits                                |
| `qmc_fill_rect(x0, x1, y0, y1, v0)`                     | add a rectangular potential                                    |
| `qmc_potential()`                                       | pointer to the `double` potential array (write to paint walls) |
| `qmc_wavepacket(x0, y0, kx, ky, sigma)`                 | Gaussian packet with momentum                                  |
| `qmc_step(steps)`                                       | advance with `soft_evolve_2d`                                  |
| `qmc_density()`                                         | pointer to `float` \|psi\|^2                                   |
| `qmc_norm()`                                            | probability left in the box                                    |

Butterfly module (`qmc_butterfly.wasm`, needs three no-op WASI imports for
stdio: `fd_write`, `fd_seek`, `fd_close`; see the page for a 10-line shim):

| function                     | purpose                                                                                |
| ---------------------------- | -------------------------------------------------------------------------------------- |
| `qmc_bf_edges(p, q)`         | band intervals at flux p/q, written to `qmc_bf_buffer()` as (lo, hi) pairs             |
| `qmc_bf_tknn(p, q, r)`       | Hall number of the gap above the lowest r bands (TKNN); `INT_MIN` if the gap is closed |
| `qmc_bf_chern(p, q, r, n_k)` | the same number computed numerically by the library (NaN if closed or grid too coarse) |

The numerical Chern number is an integer only when the k-grid resolves the Berry curvature;
for large Hall numbers use `n_k` of at least about `2|t| + 8`.

Dirac module (`qmc_dirac.wasm`, no imports;
one global simulation, N = 2048, dx = 0.1, dt = 0.04, units hbar = m = c = 1):

| function                                               | purpose                                                                 |
| ------------------------------------------------------ | ----------------------------------------------------------------------- |
| `qmc_dc_klein(k0, v0)`                                 | positive-energy packet at x = -50 aimed at a step of height v0 at x = 0 |
| `qmc_dc_zitter()`                                      | wide packet of spinor $(1, i)/\sqrt{2}$ at rest, no potential           |
| `qmc_dc_step(n)`                                       | advance n steps with `dirac_evolve_1d`                                  |
| `qmc_dc_density()` / `qmc_dc_buffer()`                 | fill / locate scratch buffer: \|upper\|^2, \|lower\|^2, V               |
| `qmc_dc_right()`, `qmc_dc_position()`, `qmc_dc_norm()` | transmitted probability, mean position, total probability               |
| `qmc_dc_exact(k0, v0)`                                 | closed-form transmission (`dirac_step_transmission`)                    |

Orbital module (`qmc_orbital.wasm`, no imports): `qmc_orb_sample(n, l, m, real, count, seed)` draws points from $\|\psi_{nlm}\|^2$ via `hydrogen_orbital_sample` into `qmc_orb_buffer()` (x, y, z, phase per point in Bohr radii); `qmc_orb_energy_ev(n)` gives the Bohr energy.

Bloch module (`qmc_bloch.wasm`, no imports): `qmc_bl_reset(x, y, z)` sets the state, `qmc_bl_step(omega, delta, gamma1, gamma_phi, steps)` advances `steps` RK4 steps of 0.01 with `bloch_evolve`, `qmc_bl_vec()` points at the Bloch vector, `qmc_bl_rabi(t, omega, delta)` is the closed-form Rabi curve.

Anderson module (`qmc_anderson.wasm`, no imports): `qmc_an_reset(w, seed)` draws new disorder of strength `w` and puts the particle on the middle of the 801-site chain, `qmc_an_step(n)` advances `n` RK4 steps of 0.01 (hop t = 1) with `anderson_evolve`, `qmc_an_buffer()` holds $|\psi|^2$ then the on-site energies (refreshed by every step), `qmc_an_width()`, `qmc_an_ipr()`, `qmc_an_energy()` and `qmc_an_norm()` are the observables, and `qmc_an_xi(w, e)` is the weak-disorder localisation length.

Quantum-walk module (`qmc_qwalk.wasm`, no imports): `qmc_qw_reset(theta, coin)` restarts the walk (coin 0 = up, 1 = down, 2 = symmetric), `qmc_qw_step(n)` advances `n` steps (at most 480 in total), `qmc_qw_buffer()` holds the quantum then the classical probability of each of the 1024 sites, and `qmc_qw_variance()` / `qmc_qw_mean()` / `qmc_qw_asymptote()` give the observables and the ($1 - \sin(\theta)$ limit of variance / $t^2$.

Landau-Zener module (`qmc_lz.wasm`, no imports): `qmc_lz_reset(omega, rate, amp, passes)` starts a sweep of the detuning between -amp and +amp in the lower adiabatic state, `qmc_lz_step(duration)` advances it (returns 1 when finished), `qmc_lz_vec()` points at the Bloch vector, `qmc_lz_delta()` / `qmc_lz_upper()` / `qmc_lz_exact()` give the detuning, upper-level population and closed-form jump probability, and `qmc_lz_scan_rate` / `qmc_lz_scan_amp` fill `qmc_lz_scan_buffer()` with whole-sweep scans.

Kicked-rotor module (`qmc_rotor.wasm`, no imports): `qmc_kr_reset(k, hbar, seed)` restarts the quantum state (momentum eigenstate $m = 0$ on 2048 angle points) and a classical ensemble of 3000 points, `qmc_kr_step(n)` applies `n` kicks (at most 400 in total), `qmc_kr_qbuffer()` holds the quantum momentum probabilities (index $j$ is $m = j - 1024$), `qmc_kr_cbuffer()` the classical (theta, p) pairs, and `qmc_kr_qm2()` / `qmc_kr_cm2()` the two $<m^2>$ values.

SSH-chain module (`qmc_ssh.wasm`, no imports): `qmc_ssh_set(w, noise, m, seed)` builds the 20-cell chain (v = 1, inter-cell hopping w, random hopping `noise`, staggered potential m) and diagonalises it; `qmc_ssh_energies()` and `qmc_ssh_vectors()` hold the 40 eigenvalues and eigenvectors, `qmc_ssh_sweep(w0, w1, n)` fills `qmc_ssh_sweep_buffer()` with the spectrum for n values of w, and `qmc_ssh_winding()`, `qmc_ssh_zak()`, `qmc_ssh_xi()` and `qmc_ssh_gap()` give the topological invariants and length scales.

Ising-quench module (`qmc_tfim.wasm`, no imports): `qmc_tq_set(h_i, h_f)` fixes the two fields (J = 1), `qmc_tq_eval(t)` fills `qmc_tq_zz_buffer()` with <sz_0 sz_r> for r = 0..40 on a 256-site ring and `qmc_tq_last_mx()` / `qmc_tq_last_rate()` with the magnetisation and Loschmidt rate, `qmc_tq_curve(t_max, n)` fills the rate and magnetisation curve buffers, and `qmc_tq_tcrit(n)` gives the critical times of the dynamical transition.

Haldane module (`qmc_haldane.wasm`, imports three WASI stubs): `qmc_hd_set(t2, phi, mass)` builds the model (t1 = 1), `qmc_hd_compute(n_k)` diagonalises a 14-cell ribbon at `n_k` wave numbers into `qmc_hd_energy_buffer()` and `qmc_hd_mean_buffer()` (position of each state across the ribbon, -1 to +1), `qmc_hd_profile(ik, n_k, state)` fills `qmc_hd_prob_buffer()` with one state's weight per cell, and `qmc_hd_chern()` / `qmc_hd_limit()` give the Chern number and phase boundary.

Wigner module (`qmc_wigner.wasm`, no imports): `qmc_wg_state(kind, amp, param, chi_t)` builds a state (0 coherent, 1 Fock, 2 even cat, 3 odd cat, 4 squeezed), applies the Kerr phase `exp(-i chi_t n(n-1))` and fills `qmc_wg_w_buffer()` (128 x 128, row = $p$, column = $x$, window `qmc_wg_extent()`), the marginals `qmc_wg_x_buffer()` / `qmc_wg_p_buffer()` and the photon-number distribution `qmc_wg_pn_buffer()`. `qmc_wg_nbar()`, `qmc_wg_negativity()`, `qmc_wg_wmin()` and `qmc_wg_wmax()` are the readouts.

Grover module (`qmc_grover.wasm`, no imports): `qmc_gv_set_n(n)` picks $2^n$ items (n up to 10), `qmc_gv_mark(i, flag)` marks or unmarks one (fewer than $N/2$ in total), `qmc_gv_step()` / `qmc_gv_goto(k)` apply Grover iterations to `qmc_gv_amp_buffer()`, and `qmc_gv_curve(count)` fills the simulated and closed-form success curves. `qmc_gv_success()`, `qmc_gv_theta()` and `qmc_gv_kopt()` are the readouts.

Bell module (`qmc_bell.wasm`, no imports): `qmc_bl_set(a, a2, b, b2, V)` fixes the four detector angles (radians) and the Werner-state visibility, `qmc_bl_chsh()` / `qmc_bl_local_chsh()` give S for the quantum state and for the hidden-angle local model, `qmc_bl_max()` the best S over a $\pi/4$ grid, `qmc_bl_curve(n)` fills the correlation-versus-angle buffers and `qmc_bl_sample(n, seed)` simulates n pairs per setting (`qmc_bl_est(i)`, `qmc_bl_sigma()`).

Zeno module (`qmc_zeno.wasm`, no imports): the drive is $\Omega = 1$ and the run lasts $T = \pi$. `qmc_zn_pulsed_curves(n)` fills the survival-versus-time curve for $n$ projective measurements and the survival-versus-N curve, `qmc_zn_continuous(gamma)` the survival curve under continuous dephasing, `qmc_zn_rates(gmax)` the decay rate versus dephasing, and `qmc_zn_sample(n, n_traj, seed)` simulates trajectories (`qmc_zn_alive_buffer()`).

Corral module (`qmc_corral.wasm`, no imports): the first call builds a table of every mode of the unit disc with k <= 26 on an 80 x 80 grid (about 10 MB, ~0.3 s). `qmc_cb_show_mode(i)` fills `qmc_cb_field_buffer()` with mode `i` (quantum numbers from `qmc_cb_mode_m/n/sine/k`), `qmc_cb_launch(x0, y0, s, px, py)` expands a Gaussian packet in the modes and `qmc_cb_evolve(t)` writes its density at time t into the same buffer.

Ring module (`qmc_ring.wasm`, no imports): $t = 1$, rings of 4 .. 32 sites with leads at sites 0 and b. `qmc_ar_spectrum(n)` fills the sorted levels against flux in $[-1, 1]$, `qmc_ar_current(n, ne)` the persistent current of ne electrons, `qmc_ar_trans_flux(n, b, e)` and `qmc_ar_trans_energy(n, b, flux)` the transmission curves, `qmc_ar_map(n, b)` a 64 x 64 transmission map over flux and energy and `qmc_ar_point(n, b, flux, e)` one value (`qmc_ar_last_reflection` gives R).

Tunnel module (`qmc_tunnel.wasm`, no imports): units $\hbar^2/2m = 1$. `qmc_dt_set(v1, b1, w, b2, v2, vw)` defines barrier, well and barrier, `qmc_dt_curves(emax)` fills T(E) and the dwell time, `qmc_dt_point(e)` returns T at one energy (`qmc_dt_last_reflection` gives R) and `qmc_dt_profile(e, xmin, xmax)` fills $|\psi|^2$ and the potential along x.

Josephson module (`qmc_josephson.wasm`, no imports): time in units of 1/omega_p. `qmc_jj_set(q, i_dc, i_ac, omega)` sets the junction, `qmc_jj_reset(phi, v)` and `qmc_jj_advance(t)` step the live pendulum (`qmc_jj_phi/v/time` read it), and `qmc_jj_iv(imax)` fills the up and down sweep I-V curves.

Double-well module (`qmc_dwell.wasm`, no imports): $\hbar = m = 1$, a 500-point grid on \[-5, 5\]. `qmc_dw_set(B, a)` solves the lowest eight states of $V = B(1 - x^2/a^2)^2$ (`qmc_dw_e_buffer`, `qmc_dw_psi_buffer`, `qmc_dw_splitting`, `qmc_dw_wkb`), `qmc_dw_launch(center, sigma)` projects a Gaussian and `qmc_dw_evolve(t)` fills the density buffer and returns P(x < 0); `qmc_dw_split_curve(bmax)` gives the splitting and WKB curves against barrier height.

Zeeman module (`qmc_zeeman.wasm`, no imports): $H = \xi L.S + b (L_z + 2 S_z)$ for a p electron, energies in the unit of $\xi$ and b. `qmc_zm_curves(xi, bmax, npts)` fills `qmc_zm_curves_buffer` with the six 2P levels and the two 2S levels (level $k$ at $k \* npts$, `qmc_zm_mj_buffer` holds $2 m_j$), `qmc_zm_lines(xi, b)` fills `qmc_zm_lines_buffer` with twelve rows \[energy, intensity, q, P state, spin\] for the 2S to 2P lines and returns the number with non-zero weight and `qmc_zm_moment(xi, b, k)` returns $dE/db$ of state $k$.

Regenerate the butterfly images:

```sh
node wasm/tools/dump_butterfly.mjs 50 /tmp/butterfly.json
python3 wasm/tools/make_butterfly_png.py /tmp/butterfly.json \
    docs/src/playground/butterfly.png docs/src/playground/butterfly_grow.gif
```
