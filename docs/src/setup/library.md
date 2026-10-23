# Using QMC as a Library

QMC installs like any C library: static and shared libraries, headers, and a
`pkg-config` file. Nothing here requires the repository to be present at build
time of your project.

## Install

```sh
make install PREFIX=$HOME/.local          # or /usr/local (needs sudo)
make install DESTDIR=/tmp/pkg PREFIX=/usr # staged, for packagers
```

`make install` always builds a **release** library (no AddressSanitizer) in a
separate `build-install/` directory, so your normal sanitized development build
is never installed. Options such as `USE_LAPACK=1`, `USE_HDF5=1` and
`PLOT_BACKEND=...` are honoured, and the resulting link flags are recorded in
`qmc.pc`, so consumers do not need to know how QMC was configured.

Installed layout:

| Path                                               | Contents                                    |
| -------------------------------------------------- | ------------------------------------------- |
| `$PREFIX/include/qmc/qmc.h`                        | umbrella header                             |
| `$PREFIX/bin/qmc-help`                             | command-line module catalog                 |
| `$PREFIX/include/qmc/{core,physics,export,latex}/` | per-module headers                          |
| `$PREFIX/lib/libqmc.a`                             | static library                              |
| `$PREFIX/lib/libqmc.so.<MAJOR>.<MINOR>.<PATCH>`    | shared library (SONAME `libqmc.so.<MAJOR>`) |
| `$PREFIX/lib/pkgconfig/qmc.pc`                     | `pkg-config` metadata                       |

## Use it

```c
#include <qmc/qmc.h>      // everything, or include single modules:
// #include <qmc/physics/hydrogen.h>

#include <stdio.h>

int main(void) {
  printf("QMC %s: E(1s) = %g J\n", qmc_version(), hydrogen_energy_level(1));
  return 0;
}
```

```sh
export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig
cc demo.c $(pkg-config --cflags --libs qmc) -o demo
LD_LIBRARY_PATH=$HOME/.local/lib ./demo   # not needed if PREFIX is on the loader path
```

For a fully static link, pass the archive and the remaining flags yourself:

```sh
cc demo.c $(pkg-config --cflags qmc) $PREFIX/lib/libqmc.a \
   $(pkg-config --libs qmc | sed 's/-lqmc//') -o demo
```

`qmc.h` deliberately leaves out `export/hdf5_writer.h`, `export/netcdf_writer.h`
and `export/plot.h`, which need optional third-party libraries; include those
directly when you built QMC with the matching option.

## Finding your way around

QMC has more than a hundred modules, so it can describe itself. From code:

```c
#include <qmc/qmc.h>

qmc_help(stdout, NULL);             // every category and module
qmc_help(stdout, "tight_binding");  // what it does, #include line, where to start,
                                    // guide, example, test, API page
qmc_help(stdout, "chemistry");      // all modules in a category
qmc_help(stdout, "chern");          // keyword search
```

From a shell (installed with the library):

```sh
qmc-help                    # overview
qmc-help tight_binding      # one module
qmc-help dmrg vqe           # several topics
qmc-help --list             # module names, one per line
```

The same catalog is available as data (`qmc_module_count()`, `qmc_module_at()`,
`qmc_module_find()`, `qmc_module_search()`, `qmc_category_at()`) if you want to
build a menu or a lookup on top of it.

The [API Reference](../api/index.md) lists every type, constant and function of
every module, with the comments from the headers and the example and test that
show each module in use. It, the catalog and the include list of `qmc.h` are all
generated from the public headers by `tools/gen_api.py`
(`make api-docs`); CI fails if they are stale, and a new header must be given a
category and a one-line summary in `tools/api_meta.py`.

## Versioning

`core/qmc_version.h` is the single source of the version. It provides
`QMC_VERSION_MAJOR/MINOR/PATCH`, `QMC_VERSION_STRING` and `QMC_VERSION_NUMBER`
(`MAJOR*10000 + MINOR*100 + PATCH`) at compile time, and `qmc_version()` /
`qmc_version_number()` for the library you actually linked:

```c
#if QMC_VERSION_NUMBER >= 100          /* QMC >= 0.1.0 */
#endif
if (qmc_version_number() != QMC_VERSION_NUMBER) { /* headers and .so differ */ }
```

While QMC is at `0.x`, minor releases may change the API; the SONAME only
promises compatibility within one `MAJOR`.

## Naming caveat

Most public functions carry a module prefix (`tb_`, `qaoa_`, `json_`, ...), but
a handful of older ones are generic (`fft`, `linspace`, `factorial`, `sigma_x`,
`V_harmonic`, ...). If your project defines the same name you will get a
duplicate-symbol error at link time. A `qmc_` prefix for these is planned.

## Checking an install

`make check-install` installs into a scratch prefix, compiles every installed
header on its own, builds and runs a small consumer against it (shared and
static, through `pkg-config`) and audits `libqmc.so` (SONAME, no sanitizer
runtime, no exported symbol missing from the public headers). CI runs it for
gcc, clang and a `USE_LAPACK=1 USE_HDF5=1` configuration.
