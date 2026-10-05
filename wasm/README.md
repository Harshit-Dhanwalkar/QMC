# QMC in browser (WebAssembly)

`qmc_wasm.c` is a thin wrapper around the library's own 2D split-operator
solver (`physics/soft.c` -> `core/fft/fft2d.c`). The browser page in
`docs/src/playground/index.html` calls it every animation frame.

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

| function                                                | purpose                                                    |
| ------------------------------------------------------- | ---------------------------------------------------------- |
| `qmc_init(n, box, dt)`                                  | allocate an `n x n` grid (n a power of two, <= 256)        |
| `qmc_build_slits(wall_x, thick, slits, sep, width, v0)` | vertical wall with 1 or 2 slits                            |
| `qmc_fill_rect(x0, x1, y0, y1, v0)`                     | add a rectangular potential                                |
| `qmc_potential()`                                       | pointer to `double` potential array (write to paint walls) |
| `qmc_wavepacket(x0, y0, kx, ky, sigma)`                 | Gaussian packet with momentum                              |
| `qmc_step(steps)`                                       | advance with `soft_evolve_2d`                              |
| `qmc_density()`                                         | pointer to `float` \|psi\|^2                               |
| `qmc_norm()`                                            | probability left in box                                    |
