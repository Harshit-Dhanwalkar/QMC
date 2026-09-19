#ifndef QMC_GEOMETRY_OPTIMIZATION_H
#define QMC_GEOMETRY_OPTIMIZATION_H

/*
 * Generic Cartesian geometry optimizer, driven by an energy+analytic-gradient
 * callback (Not tied to RHF/UHF or to any one molecule's basis-set
 * construction)
 * NOTE: eg_49_geometry_optimization.c and eg_66_uhf_geometry_optimization.c
 * each hand-roll a scalar (bond-length-only) steepest-descent loop that
 * rebuilds the whole basis/molecule per force evaluation; this generalizes that
 * idea to a full n_atoms*3 Cartesian coordinate vector so it also covers
 * non-collinear (bent, polyatomic) geometries, while keeping the same "rebuild
 * from these coordinates and hand back energy+gradient" shape those examples
 * already use
 */

/*
 * Evaluate energy and its Cartesian gradient at a given geometry.
 *  coords    : input, length n_atoms*3 (coords[3 * A + d], d=0,1,2 for x,y,z)
 *  n_atoms   : number of atoms (fixed for the lifetime of the optimization)
 *  user_data : opaque pointer forwarded from the optimizer call, typically
 *              holding the basis-set kind, charges, electron count, and any SCF
 *              settings needed to rebuild the basis at the new geometry
 * energy_out : callee writes the total (electronic + nuclear repulsion) energy
 *              at coords
 * grad_out   : preallocated length n_atoms*3 by the caller; callee writes
 *              dE/dR_A[d] into grad_out[3*A+d]
 *
 * Returns 0 on success (e.g. SCF converged), nonzero on failure (e.g. SCF did
 * not converge at this geometry) - optimizer aborts on a nonzero return and
 * reports it as not converged
 */
typedef int (*geom_energy_grad_fn)(const double *coords, int n_atoms,
                                   void *user_data, double *energy_out,
                                   double *grad_out);

typedef struct {
  double *coords;   /* final geometry, length n_atoms*3, caller frees */
  double energy;    /* energy at the final geometry */
  double grad_norm; /* max |component| of the gradient at the final geometry */
  int iterations;
  int converged;
} geometry_optimization_result_t;

/*
 * Steepest descent on the analytic gradient, with simple backtracking (halving
 * the step) whenever a step would raise the energy - the same "move downhill
 * along the force", made robust to a poorly-scaled initial step and generalized
 * to n_atoms*3 coordinates instead of one scalar bond length
 *
 * func/user_data: see geom_energy_grad_fn above
 * coords0       : starting geometry, length n_atoms*3 (not modified)
 * n_atoms       : number of atoms
 * step0         : initial step length along -gradient
 * grad_tol      : converged when every |grad component| < grad_tol
 * max_iter      : outer iteration cap (each outer iteration may backtrack the
 *                 step internally without consuming an outer iteration)
 *
 * Returns NULL on invalid input (n_atoms <= 0, step0 <= 0, grad_tol <= 0,
 * max_iter <= 0) or allocation failure. Otherwise always returns a result
 * (check -> converged)
 */
geometry_optimization_result_t *
optimize_geometry_steepest_descent(geom_energy_grad_fn func, void *user_data,
                                   const double *coords0, int n_atoms,
                                   double step0, double grad_tol, int max_iter);

void geometry_optimization_result_free(geometry_optimization_result_t *res);

#endif
