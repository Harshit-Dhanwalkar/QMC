#ifndef QMC_DMRG_H
#define QMC_DMRG_H

#include "../core/matrix.h"

/*
 * Infinite-system Density Matrix Renormalization Group (Reference: White, Phys.
 * Rev. Lett. 69, 2863 (1992)) for open-boundary spin-1/2 XXZ chain:
 *
 *   H = \sum_{j=0}^{N-2} [ Jz S^z_j S^z_{j+1}
 *                        + (Jxy/2)(S^+_j S^-_{j+1} + S^-_j S^+_{j+1}) ]
 *
 * Site-operator convention basis order [up, down], S^z_up = +1/2, S^z_down =
 * -1/2, S^+ raises down -> up), but for a finite open chain (no periodic
 * wraparound bond), which is natural setting for block-growth
 * construction below rather than periodic ring spin_chain.h targets
 *
 * Algorithm (infinite-system sweep only):
 *   1. Start with a 1-site "block" (dim 2, H=0)
 *   2. Enlarge block by one site: new block = old block (x) new site,
 *      with coupling bond between them added into H
 *   3. Build a superblock = enlarged block (x) mirror-image enlarged block
 *      (reflection symmetry stands in for an explicit environment block),
 *      and find its ground state by dense diagonalization
 *   4. Diagonalize ground state's reduced density matrix on enlarged-block
 *      subspace; keep m eigenvectors of largest eigenvalue ("density-matrix
 *      truncation", method's namesake) These become new, truncated block
 *      basis; go to 2
 *   5. Stop once superblock's total site count reaches target N
 *
 * Each step grows total chain by 2 sites (one on each mirrored side), so N
 * is reached exactly when even, or overshot by 1 when odd
 *
 * WARN: this is infinite-system algorithm only there is no finite-size sweeping
 * back over previously-built (now-truncated) blocks yet. For a gapped chain
 * this still reaches machine precision once m is large enough (as above), but
 * for a fixed, insufficient m on a critical/gapless chain (e.g. isotropic
 * Heisenberg, Jz=Jxy), infinite-algorithm truncation error is somewhat larger
 * than a fully swept finite-system calculation would give for same m.
 * Finite-size sweeps are a natural follow-on extension, not implemented here
 */

typedef struct {
  int dim;           /* current (possibly truncated) block basis dimension */
  cmatrix_t *H;      /* block Hamiltonian, dim x dim */
  cmatrix_t *Sz_end; /* S^z acting on block's open (most recently added)
                       site, in block's basis, dim x dim */
  cmatrix_t *Sp_end; /* S^+ acting on block's open site, dim x dim */
} dmrg_block_t;

/* Single-site block: dim=2, H=0, Sz_end/Sp_end = bare site operators.
 * Returns NULL on allocation failure */
dmrg_block_t *dmrg_block_init(void);
void dmrg_block_free(dmrg_block_t *b);

/* Result of an infinite-algorithm DMRG run */
typedef struct {
  double energy;          /* superblock ground-state energy at N_reached */
  double energy_per_site; /* energy / N_reached */
  int N_reached;          /* total chain length actually reached (even) */
  int truncation_dim; /* m actually used on final truncation step  performed (<=
                         m_max requested; may be less if pre-truncation
                         dimension was smaller) */
  double truncation_error; /* discarded reduced-density-matrix weight (1 - sum
                              of kept eigenvalues) on that step; 0 if N_target
                              was reached before any  truncation was needed */
} dmrg_result_t;

/*
 * Run infinite-system algorithm until superblock's total site count is >=
 * N_target, truncating to at most m_max states at every intermediate step.
 * N_target is rounded up internally to nearest even number (each step grows
 * chain by 2 sites, one per mirrored side)
 *
 * Returns a malloc'd dmrg_result_t (free with free()), or NULL for invalid
 * input (N_target < 2, m_max < 1, Jz/Jxy non-finite) or allocation failure
 */
dmrg_result_t *dmrg_run(int N_target, double Jz, double Jxy, int m_max);

#endif
