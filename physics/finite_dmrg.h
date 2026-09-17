#ifndef QMC_FINITE_DMRG_H
#define QMC_FINITE_DMRG_H

#include "dmrg.h"

/*
 * Finite-system Density Matrix Renormalization Group sweeps (Reference: White,
 * Phys. Rev. B 48, 10345 (1993)) for the same open-boundary spin-1/2 XXZ chain
 * as dmrg.h/dmrg.c:
 *
 *   H = \sum_{j=0}^{N-2} [ Jz S^z_j S^z_{j+1}
 *                        + (Jxy/2)(S^+_j S^-_{j+1} + S^-_j S^+_{j+1}) ]
 *
 * NOTE: dmrg.c's infinite algorithm only ever compares an enlarged block
 * against its own mirror image, so every intermediate block is built using an
 * environment that is itself still small and approximate. Finite-size sweeping
 * fixes the total chain length N up front, keeps a cache of every block
 * size 1..N-1 seen so far (reflection symmetry means one cache serves as both
 * left- and right-block storage) and repeatedly re-diagononalizes superblocks
 * built from one growing side plus a previously-cached environment side of
 * complementary size, re-truncating only growing side each step. Sweeping the
 * growing side back and forth from one end of the chain to the other lets every
 * block eventually be rebuilt using a much better (larger effective bond
 * dimension / already-optimized) environment than the infinite algorithm ever
 * had access to, which is  reason finite sweeps improve on the infinite result
 * at fixed m
 *
 * Algorithm:
 *   1. Warmup: run the infinite algorithm (mirrored superblocks) to build
 *      block cache entries for every size 1..N/2 but keep every intermediate
 *      truncated block instead of discarding all but the last, and truncate the
 *      final size-N/2 block too
 *   2. Sweep: pick a growing-side size l starting at 2, form the superblock
 *      from enlarge(cache[l-1]) (system) and cache[N-l] (environment, read from
 *      cache only, not modified), find its ground state, truncate the enlarged
 *      side to at most m_max states via its reduced density matrix, and store
 *      the result back into cache[l]. Advance l by one and repeat until l
 *      reaches N-1
 *   3. Reverse direction (shrunk side now becomes the side that grows, using
 *      the fresh blocks the other side just produced as its improving
 *      environment) and sweep back down to l=1. One outward-and-back pass is
 *      one full sweep; repeat for n_sweeps full sweeps
 *   4. Report the ground energy at the most balanced bipartition (l = N/2)
 *      using the final cache contents, plus the central-bond energy recorded
 *      after each individual sweep for convergence diagnostics
 *
 * N is rounded up internally to nearest even number, exactly as in dmrg_run, so
 * that N/2 is always an integer valid cache index
 */

typedef struct {
  double energy;           /* ground energy at final l=N/2 bipartition */
  double energy_per_site;  /* energy / N_reached */
  int N_reached;           /* total chain length actually used (even) */
  int m_max;               /* requested cap on retained states per block */
  int n_sweeps;            /* number of full sweeps actually performed */
  double truncation_error; /* discarded RDM weight on the final l=N/2 step */
  double *sweep_energy;    /* malloc'd, length n_sweeps: l=N/2 energy recorded
                             after each full sweep, for convergence diagnostics */
} finite_dmrg_result_t;

/*
 * Run n_sweeps full finite-size sweeps for a chain of N_target sites (rounded
 * up to even), truncating every block to at most m_max states
 *
 * Returns a malloc'd finite_dmrg_result_t, or NULL for invalid input (N_target
 * < 4, m_max < 1, n_sweeps < 1, Jz/Jxy non-finite) or allocation failure
 */
finite_dmrg_result_t *finite_dmrg_run(int N_target, double Jz, double Jxy,
                                      int m_max, int n_sweeps);

void finite_dmrg_result_free(finite_dmrg_result_t *r);

#endif
