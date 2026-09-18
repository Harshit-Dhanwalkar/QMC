# Platform Support

QMC is developed and primarily tested on Linux. This page summarizes what
works elsewhere and what doesn't yet.

## Linux

Fully supported. This is the reference platform for the whole test suite,
all four `SANITIZE` x `USE_LAPACK` build combinations, `cppcheck`, `clang-tidy`,
and `valgrind`.

## macOS

Partially supported, not regularly tested on real hardware.

The Makefile already auto-detects Homebrew (`brew --prefix`) and, when
present, prefers a Homebrew-installed GCC over Apple's `clang`-based `gcc`
shim, and adds a Homebrew LAPACK's include/lib paths automatically so
`make USE_LAPACK=1` can find it (see the top of the `Makefile`).

## Windows

Not natively supported; **use WSL or MSYS2/MinGW**, which provide a POSIX
shell and toolchain and make QMC behave exactly like the Linux case above.
