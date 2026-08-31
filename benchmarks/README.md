# Benchmarks

Standalone benchmark programs, built the same way as `eg_%` examples and `test_%` tests
From the repo root:

```bash
make run-benchmarks   # builds + runs bench_accuracy, bench_eigensolver,
                      # bench_vmc_convergence; writes output/bench_<name>.dat
make benchmarks       # just builds them, doesn't run
```

Add `PLOT_BACKEND=NONE` to any of these if you don't have the GR submodule built and don't want the fallback chain to try GNUPLOT/MATPLOTLIB first.

Add `SANITIZE=0` to avoid paying the AddressSanitizer overhead for what are meant to be timing runs.

## bench_accuracy: VMC/DMC vs. published energies

Runs VMC (`vmc_run`, optimizing the Jastrow parameter `b` via `vmc_optimize_b`) and DMC (`dmc_run`) for the He isoelectronic two-electron sequence : $He$, $Li+$, $Be2+$ and then compares against published Pekeris-type clamped-nucleus nonrelativistic ground-state energies. `physics/vmc.c` and `physics/dmc.c` only support two-electron systems, so this uses Be2+ rather than neutral Be (4 electrons, unsupported), and omits $H-$ (its variational optimum `Zeff` is far enough from the bare nuclear charge that this single-screened-orbital ansatz makes it a misleading benchmark case rather than an informative one).

```bash
./build/bench_accuracy
```

Typical run: DMC lands within 0.002-0.01 Hartree of the published values across all three systems, VMC sits above them by 0.03-0.05 Hartree as the variational principle requires (`examples/eg_28_vmc_helium.c` / `eg_29_dmc_helium.c`'s He numbers).

## bench_eigensolver: hand-rolled vs LAPACK complex-Hermitian eigensolver

`cmatrix_eigh_complex()` (`core/linalg/complex_eigh.c`) picks its backend at **compile time** via `#ifdef USE_LAPACK`, not runtime, so unlike `bench_reference_vs_blocked.c`'s runtime mode switch, this benchmark has to be built twice to compare the two paths:

```bash
make PLOT_BACKEND=NONE SANITIZE=0 build/bench_eigensolver              # hand-rolled (Jacobi)
./build/bench_eigensolver 10 20 50 100 > eigensolver_handrolled.dat

rm -rf build && mkdir -p build/core/linalg build/core/fft build/core/ode \
  build/core/special build/physics build/export build/latex
make PLOT_BACKEND=NONE SANITIZE=0 USE_LAPACK=1 build/bench_eigensolver  # LAPACK zheev
./build/bench_eigensolver 10 20 50 100 > eigensolver_lapack.dat
```

Accuracy is checked with a residual test independent of any reference solver (`max_i ||H v_i - lambda_i v_i|| / ||v_i||` against the input matrix directly), not a diff against a second solver computed in the same run: the compile-time split makes that impossible within one binary anyway.

Hand-rolled `cmatrix_eigh` is cyclic Jacobi on the `2n x 2n` real embedding, so it gets slow well before $n=500$; keep default sizes modest unless you specifically want to watch that curve.

## bench_vmc_convergence: VMC statistical convergence, He atom

Sweeps `n_samples` (100 through 50000) at fixed Jastrow `b=0.15`, `Zeff=Z=2.0`, and reports `vmc_run`'s block-averaged mean/error/wall-time against the exact He ground state. Isolates pure `1/sqrt(N)` statistical convergence from optimizer noise by not re-optimizing `b` per sample count.

Block size is chosen to guarantee at least ~5 blocks even at the smallest sample count.

```bash
./build/bench_vmc_convergence
```

## bench_reference_vs_blocked: tridiagonal QL eigenvector cache-optimization

Not included in `make benchmarks`/`run-benchmarks` (it needs explicit `<N> <mode>` arguments, so running it with none as the aggregate loop would just prints a usage message and exits). Build and run it directly:

```bash
make PLOT_BACKEND=NONE SANITIZE=0 build/bench_reference_vs_blocked
./build/bench_reference_vs_blocked <N> <mode>
#   mode = 0       -> reference (unblocked) tridiagonal QL implementation
#   mode = <block> -> column-blocked variant with that block size (e.g. 64)
```

Compares wall-clock time between a reference (unblocked) and a column-blocked tridiagonal QL eigenvector computation - both mathematically identical, differing only in memory-access pattern (the blocked variant replays recorded Givens rotations against the eigenvector matrix one column-block at a time, so the working set per pass is `block_size*n` instead of `n*n`). See the file's own header comment for suggested `cachegrind`-based diagnostics.

## instrument_iteration_counts: QL sweep/rotation-count instrumentation

A standalone diagnostic tool (not part of the `bench_%` Makefile pattern - build it directly with `gcc`, as its own header comment documents), counting outer QL sweeps and inner rotation steps on the same clustered-diagonal test matrix `bench_reference_vs_blocked` uses, across increasing `N`. Useful for checking whether the QL algorithm's iteration count itself scales as expected (roughly `O(N)` sweeps) independent of any wall-clock/cache effects.

```bash
gcc -O2 -o instrument_iteration_counts benchmarks/instrument_iteration_counts.c -lm
for N in 400 800 1600 3200; do ./instrument_iteration_counts $N; done
```

---

OUTPUT:

1. Batch size 3200 0

QMC/benchmarks$ valgrind --tool=cachegrind --cache-sim=yes ./bench 3200 0

```
==15613== Cachegrind, a cache and branch-prediction profiler
==15613== Copyright (C) 2002-2017, and GNU GPL'd, by Nicholas Nethercote et al.
==15613== Using Valgrind-3.18.1 and LibVEX; rerun with -h for copyright info
==15613== Command: ./bench 3200 0
==15613==
--15613-- warning: L3 cache found, using its data for the LL simulation.

N= 3200 mode=0      time=1400.6608s
==15613==
==15613== I   refs:      451,632,030,099
==15613== I1  misses:              1,617
==15613== LLi misses:              1,610
==15613== I1  miss rate:            0.00%
==15613== LLi miss rate:            0.00%
==15613==
==15613== D   refs:      150,474,848,053  (90,310,710,442 rd   + 60,164,137,611 wr)
==15613== D1  misses:      5,155,874,685  ( 5,134,638,202 rd   +     21,236,483 wr)
==15613== LLd misses:      3,738,817,336  ( 3,737,529,944 rd   +      1,287,392 wr)
==15613== D1  miss rate:             3.4% (           5.7%     +            0.0%  )
==15613== LLd miss rate:             2.5% (           4.1%     +            0.0%  )
==15613==
==15613== LL refs:         5,155,876,302  ( 5,134,639,819 rd   +     21,236,483 wr)
==15613== LL misses:       3,738,818,946  ( 3,737,531,554 rd   +      1,287,392 wr)
==15613== LL miss rate:              0.6% (           0.7%     +            0.0%  )
```

2. Batch size 3200 128

```
QMC/benchmarks$ valgrind --tool=cachegrind --cache-sim=yes ./bench 3200 128
==28115== Cachegrind, a cache and branch-prediction profiler
==28115== Copyright (C) 2002-2017, and GNU GPL'd, by Nicholas Nethercote et al.
==28115== Using Valgrind-3.18.1 and LibVEX; rerun with -h for copyright info
==28115== Command: ./bench 3200 128
==28115==
--28115-- warning: L3 cache found, using its data for the LL simulation.
N= 3200 mode=128    time=1252.1613s
==28115==
==28115== I   refs:      456,279,498,205
==28115== I1  misses:              1,661
==28115== LLi misses:              1,654
==28115== I1  miss rate:            0.00%
==28115== LLi miss rate:            0.00%
==28115==
==28115== D   refs:      151,197,958,504  (90,996,177,934 rd   + 60,201,780,570 wr)
==28115== D1  misses:      4,086,809,145  ( 4,081,999,908 rd   +      4,809,237 wr)
==28115== LLd misses:         94,281,678  (    89,477,090 rd   +      4,804,588 wr)
==28115== D1  miss rate:             2.7% (           4.5%     +            0.0%  )
==28115== LLd miss rate:             0.1% (           0.1%     +            0.0%  )
==28115==
==28115== LL refs:         4,086,810,806  ( 4,082,001,569 rd   +      4,809,237 wr)
==28115== LL misses:          94,283,332  (    89,478,744 rd   +      4,804,588 wr)
==28115== LL miss rate:              0.0% (           0.0%     +            0.0%  )
```

QMC/benchmarks$ cg_annotate cachegrind_ref.out | head -30

```
--------------------------------------------------------------------------------
I1 cache:         32768 B, 64 B, 8-way associative
D1 cache:         49152 B, 64 B, 12-way associative
LL cache:         8388608 B, 64 B, 8-way associative
Command:          ./bench 3200 0
Data file:        cachegrind_ref.out
Events recorded:  Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
Events shown:     Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
Event sort order: Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
Thresholds:       0.1 100 100 100 100 100 100 100 100
Include dirs:
User annotated:
Auto-annotation:  on

--------------------------------------------------------------------------------
Ir                       I1mr           ILmr           Dr                      D1mr                   DLmr                   Dw                      D1mw                DLmw
--------------------------------------------------------------------------------
451,632,030,099 (100.0%) 1,617 (100.0%) 1,610 (100.0%) 90,310,710,442 (100.0%) 5,134,638,202 (99.63%) 3,737,528,340 (100.0%) 60,164,070,026 (100.0%) 21,234,984 (99.99%) 1,285,938 (99.89%)  ???:main
```

QMC/benchmarks$ cg_annotate cachegrind_blocked.out | head -30

```
--------------------------------------------------------------------------------
I1 cache:         32768 B, 64 B, 8-way associative
D1 cache:         49152 B, 64 B, 12-way associative
LL cache:         8388608 B, 64 B, 8-way associative
Command:          ./bench 3200 128
Data file:        cachegrind_blocked.out
Events recorded:  Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
Events shown:     Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
Event sort order: Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
Thresholds:       0.1 100 100 100 100 100 100 100 100
Include dirs:
User annotated:
Auto-annotation:  on

--------------------------------------------------------------------------------
Ir                       I1mr        ILmr        Dr                      D1mr                   DLmr                Dw                      D1mw               DLmw
--------------------------------------------------------------------------------
456,279,498,205 (100.0%) 1,661 (100.0%) 1,654 (100.0%) 90,996,177,934 (100.0%) 4,081,999,908 (100.0%) 89,477,090 (100.0%) 60,201,614,285 (100.0%) 4,806,149 (99.94%) 4,801,580 (99.94%)  PROGRAM TOTALS

--------------------------------------------------------------------------------
Ir                       I1mr        ILmr        Dr                      D1mr                   DLmr                Dw                      D1mw               DLmw
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
```
