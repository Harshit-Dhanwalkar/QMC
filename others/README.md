# Benchmark

Benchmark: reference (unblocked) vs column-blocked tridiagonal QL eigenvector computation. (for performance comparison)

```bash
sudo apt install valgrind

gcc -O2 -o bench bench_reference_vs_blocked.c -lm

valgrind --tool=cachegrind --cache-sim=yes ./bench 3200 0
mv cachegrind.out.* cachegrind_ref.out

valgrind --tool=cachegrind --cache-sim=yes ./bench 3200 128
mv cachegrind.out.* cachegrind_blocked.out

cg_annotate cachegrind_ref.out | head -30
cg_annotate cachegrind_blocked.out | head -30
```

---

OUTPUT:

1. Batch size 3200 0

```
QMC/others$ valgrind --tool=cachegrind --cache-sim=yes ./bench 3200 0

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
QMC/others$ valgrind --tool=cachegrind --cache-sim=yes ./bench 3200 128
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

harshitpd@pop-os:~/Documents/GITHUB/QMC/others$ cg_annotate cachegrind_ref.out | head -30

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
451,632,030,099 (100.0%) 1,617 (100.0%) 1,610 (100.0%) 90,310,710,442 (100.0%) 5,134,638,202 (100.0%) 3,737,529,944 (100.0%) 60,164,137,611 (100.0%) 21,236,483 (100.0%) 1,287,392 (100.0%)  PROGRAM TOTALS

--------------------------------------------------------------------------------
Ir                       I1mr        ILmr        Dr                      D1mr                   DLmr                   Dw                      D1mw                DLmw                file:function
--------------------------------------------------------------------------------
451,410,326,805 (99.95%) 26 ( 1.61%) 26 ( 1.61%) 90,282,448,047 (99.97%) 5,115,863,423 (99.63%) 3,737,528,340 (100.0%) 60,164,070,026 (100.0%) 21,234,984 (99.99%) 1,285,938 (99.89%)  ???:main
```

harshitpd@pop-os:~/Documents/GITHUB/QMC/others$ cg_annotate cachegrind_blocked.out | head -30

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
Ir                       I1mr           ILmr           Dr                      D1mr                   DLmr                Dw                      D1mw               DLmw
--------------------------------------------------------------------------------
456,279,498,205 (100.0%) 1,661 (100.0%) 1,654 (100.0%) 90,996,177,934 (100.0%) 4,081,999,908 (100.0%) 89,477,090 (100.0%) 60,201,780,570 (100.0%) 4,809,237 (100.0%) 4,804,588 (100.0%)  PROGRAM TOTALS

--------------------------------------------------------------------------------
Ir                       I1mr        ILmr        Dr                      D1mr                   DLmr                Dw                      D1mw               DLmw                file:function
--------------------------------------------------------------------------------
456,057,693,051 (99.95%) 32 ( 1.93%) 32 ( 1.93%) 90,967,816,516 (99.97%) 4,081,994,402 (100.0%) 89,475,436 (100.0%) 60,201,614,285 (100.0%) 4,806,149 (99.94%) 4,801,580 (99.94%)  ???:main
```

---
