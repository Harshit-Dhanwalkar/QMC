# QMC in browser (WebAssembly)

`qmc_wasm.c` is a thin wrapper around the library's own 2D split-operator solver (`physics/soft.c` -> `core/fft/fft2d.c`). The browser page in `docs/src/playground/index.html` calls it every animation frame.

`qmc_butterfly.c` does the same for the tight-binding module: the Hofstadter butterfly page (`docs/src/playground/butterfly.html`) diagonalises `tb_model_hofstadter` for every flux p/q up to q = 50, labels each gap with its Chern number (TKNN equation) and, on click, recomputes that number with the library's Fukui-Hatsugai-Suzuki `tb_chern_number`.

`qmc_dirac.c` wraps the Dirac time-evolution module (`physics/dirac_evolve.c`): `docs/src/playground/dirac.html` fires a positive-energy packet at a sharp step (Klein paradox, compared live with `dirac_step_transmission`) or shows Zitterbewegung of a packet built from both energy signs.

`qmc_orbital.c` wraps the orbital sampler (`physics/orbital_sample.c`): `docs/src/playground/orbitals.html` draws any hydrogen orbital (n up to 8, real or complex) as a rotatable 3D point cloud sampled from |psi_nlm|^2.

`qmc_bloch.c` wraps the Bloch-vector module (`physics/bloch.c`, RK4 Lindblad evolution): `docs/src/playground/bloch.html` drives and damps a qubit and plots the state on the Bloch sphere next to the closed-form Rabi curve.

`qmc_anderson.c` wraps the Anderson-localisation module (`physics/anderson.c`): `docs/src/playground/anderson.html` evolves a particle on a disordered tight-binding chain and shows it spreading, then freezing.

`qmc_qwalk.c` wraps the quantum walk (`physics/quantum_walk.c`): `docs/src/playground/quantum_walk.html` runs a coin-and-shift walk next to the classical random walk.

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

The numerical Chern number is an integer only when the k-grid resolves the Berry
curvature; for large Hall numbers use `n_k` of at least about `2|t| + 8`.

Dirac module (`qmc_dirac.wasm`, no imports; one global simulation, N = 2048,
dx = 0.1, dt = 0.04, units hbar = m = c = 1):

| function                                               | purpose                                                                 |
| ------------------------------------------------------ | ----------------------------------------------------------------------- |
| `qmc_dc_klein(k0, v0)`                                 | positive-energy packet at x = -50 aimed at a step of height v0 at x = 0 |
| `qmc_dc_zitter()`                                      | wide packet of spinor (1, i)/sqrt(2) at rest, no potential              |
| `qmc_dc_step(n)`                                       | advance n steps with `dirac_evolve_1d`                                  |
| `qmc_dc_density()` / `qmc_dc_buffer()`                 | fill / locate scratch buffer: \|upper\|^2, \|lower\|^2, V               |
| `qmc_dc_right()`, `qmc_dc_position()`, `qmc_dc_norm()` | transmitted probability, mean position, total probability               |
| `qmc_dc_exact(k0, v0)`                                 | closed-form transmission (`dirac_step_transmission`)                    |

Orbital module (`qmc_orbital.wasm`, no imports): `qmc_orb_sample(n, l, m, real, count, seed)` draws points from |psi_nlm|^2 via `hydrogen_orbital_sample` into `qmc_orb_buffer()` (x, y, z, phase per point in Bohr radii); `qmc_orb_energy_ev(n)` gives the Bohr energy.

Bloch module (`qmc_bloch.wasm`, no imports): `qmc_bl_reset(x, y, z)` sets the state, `qmc_bl_step(omega, delta, gamma1, gamma_phi, steps)` advances `steps` RK4 steps of 0.01 with `bloch_evolve`, `qmc_bl_vec()` points at the Bloch vector, `qmc_bl_rabi(t, omega, delta)` is the closed-form Rabi curve.

Anderson module (`qmc_anderson.wasm`, no imports): `qmc_an_reset(w, seed)` draws new disorder of strength `w` and puts the particle on the middle of the 801-site chain, `qmc_an_step(n)` advances `n` RK4 steps of 0.01 (hop t = 1) with `anderson_evolve`, `qmc_an_buffer()` holds |psi|^2 then the on-site energies (refreshed by every step), `qmc_an_width()`, `qmc_an_ipr()`, `qmc_an_energy()` and `qmc_an_norm()` are the observables, and `qmc_an_xi(w, e)` is the weak-disorder localisation length.

Quantum-walk module (`qmc_qwalk.wasm`, no imports): `qmc_qw_reset(theta, coin)` restarts the walk (coin 0 = up, 1 = down, 2 = symmetric), `qmc_qw_step(n)` advances `n` steps (at most 480 in total), `qmc_qw_buffer()` holds the quantum then the classical probability of each of the 1024 sites, and `qmc_qw_variance()` / `qmc_qw_mean()` / `qmc_qw_asymptote()` give the observables and the (1 - sin theta) limit of variance / t^2.

Regenerate the butterfly images:

```sh
node wasm/tools/dump_butterfly.mjs 50 /tmp/butterfly.json
python3 wasm/tools/make_butterfly_png.py /tmp/butterfly.json \
    docs/src/playground/butterfly.png docs/src/playground/butterfly_grow.gif
```
