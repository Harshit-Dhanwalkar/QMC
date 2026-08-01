# Dependencies

QMC is engineered to run with minimal external dependencies. The core mathematical routines, physics solvers, ODE integrators, and linear algebra backends are written in pure standard C (C99).

Below are the toolchains and libraries required for compilation, visualization, and documentation generation.

---

## 1. System Toolchain

- **C Compiler**: `gcc` or `clang` supporting **C99** standard and **OpenMP** (`-fopenmp`) for multi-threading.
- **Build Systems**: `make` (GNU Make) and `cmake` (v3.16+ required for building submodules).

On Debian/Ubuntu :

```bash
sudo apt update
sudo apt install build-essential cmake libgomp1
```

## 2. Visualization & Plotting Backends

QMC supports multiple plot backends configured via the `PLOT_BACKEND` Makefile flag.

| Backend                  | Flag                      | Requirements             | Notes                                                                  |
| ------------------------ | ------------------------- | ------------------------ | ---------------------------------------------------------------------- |
| GR Framework _(default)_ | `PLOT_BACKEND=GR`         | `libGR.so`, `libqhull_r` | Built‑in via Git submodule (`third_party/gr`). Renders plots directly. |
| GNUplot                  | `PLOT_BACKEND=GNUPLOT`    | `gnuplot` in `PATH`      | Automatic fallback if GR is unavailable. Pipe‑based driver.            |
| Matplotlib               | `PLOT_BACKEND=MATPLOTLIB` | `python3` + `matplotlib` | Pipe‑based driver invoking Python.                                     |
| None                     | `PLOT_BACKEND=NONE`       | None                     | Headless / No plotting.                                                |

### GR Dependencies (Linux Headless / Build Requirements)

To compile GR from source with standard features and X11/OpenGL output support:

```bash
sudo apt install -y \
  libjpeg-dev libpng-dev libfreetype6-dev libfontconfig1-dev \
  libqhull-dev zlib1g-dev libx11-dev libxt-dev libxrender-dev libxext-dev \
  libgl1-mesa-dev libglu1-mesa-dev libglfw3-dev
```

### Compiling the GR Plotting Submodule

By default, the Makefile expects GR binaries installed under `third_party/gr/install`. Build GR headlessly using CMake:

```bash
cd third_party/gr
cmake -B build -DGR_INSTALL=ON \
  -DCMAKE_INSTALL_PREFIX=$(pwd)/install \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt6=ON \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt5=ON \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt4=ON
cmake --build build --parallel
cmake --install build
cd ../..
```

Verify that `third_party/gr/install/lib/libGR.so` exists before proceeding.

## 3. $\LaTeX$ Rendering (Optional)

If $\LaTeX$ rendering is enabled for plot annotations, the build system checks for `pdflatex` (or `lualatex`) and `pdftoppm`:

```bash
sudo apt install texlive-latex-base texlive-extra-utils poppler-utils
```
