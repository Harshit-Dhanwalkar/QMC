#ifndef QMC_VIBRATIONAL_H
#define QMC_VIBRATIONAL_H

#include "molecular_hf.h"
#include "molecular_integrals.h"

/*
 * Numerical Hessian + harmonic vibrational frequency analysis, built on top of
 * the analytic RHF/UHF nuclear gradients (molecular_rhf_gradient /
 * molecular_uhf_gradient) rather than analytic second derivatives - this has no
 * analytic Hessian, so the Hessian here is a central-difference of the gradient
 * at displaced geometries (semi-numerical Hessian, analytic second derivatives
 * aren't implemented)
 *
 * Because computing the Hessian means re-running a full SCF+gradient at
 * 6*n_atoms displaced geometries, and building the basis set for a given
 * geometry is inherently molecule-specific (which basis-builder function to
 * call, how many atoms, etc.), caller supplies two small callbacks:
 *  - molecular_geometry_fn: given a new geometry (n_atoms x 3 Cartesian
 *    coordinates, bohr), builds and returns a molecular_system_t (basis
 *    functions, atom_of_basis map, and molecule_t) for that geometry
 *  - molecular_system_free_fn: frees everything molecular_geometry_fn allocated
 *    (basis functions, mol, atom_of_basis array, the molecular_system_t struct
 *    itself)
 */

typedef struct {
  basis_function_t **basis;
  int n_basis;
  molecule_t *mol;
  int *atom_of_basis; /* length n_basis, maps basis[i] to its atom index */
  int n_electrons;
} molecular_system_t;

typedef molecular_system_t *(*molecular_geometry_fn)(const double geom[][3],
                                                     int n_atoms,
                                                     void *userdata);
typedef void (*molecular_system_free_fn)(molecular_system_t *sys,
                                         void *userdata);

/*
 * Computes the RHF nuclear Hessian by central-differencing
 * molecular_rhf_gradient() at displaced geometries:
 *  H[i][j] = (dE/dR_j at R_i+h - dE/dR_j at R_i-h) / (2h)
 *
 * then symmetrized ((H+H^T)/2) to remove the small numerical asymmetry central
 * differences leave behind
 * - geom0: n_atoms x 3 equilibrium (or any reference) Cartesian geometry, bohr
 *  - h   : displacement step (bohr); 1e-3 to 1e-2 is a reasonable range, too
 *          small loses precision to floating-point cancellation, too large
 *          loses accuracy to higher-order curvature terms
 * - scf_tol, scf_max_iter: passed through to molecular_rhf() at each displaced
 *          geometry
 *
 * Returns a newly allocated flat (3 * n_atoms)^2 array (row-major, H[3 * i +
 * d1][3 * j + d2] at index (3 * i + d1) * (3 * n_atoms) + (3 * j + d2))
 * Returns NULL on invalid input (n_atoms<=0, h<=0, NULL callbacks) or if any
 * displaced-geometry SCF fails to converge.
 */
double *molecular_rhf_hessian_numerical(molecular_geometry_fn build_system,
                                        molecular_system_free_fn free_system,
                                        void *userdata, const double geom0[][3],
                                        int n_atoms, double scf_tol,
                                        int scf_max_iter, double h);

/* Same as molecular_rhf_hessian_numerical, but for open-shell systems via
 * molecular_uhf_gradient(). The geometry-builder callback's returned
 * molecular_system_t additionally needs n_alpha/n_beta - reuse n_electrons for
 * n_alpha and pass n_beta separately here (molecular_uhf's convention: n_alpha
 * >= n_beta ) */
double *molecular_uhf_hessian_numerical(molecular_geometry_fn build_system,
                                        molecular_system_free_fn free_system,
                                        void *userdata, const double geom0[][3],
                                        int n_atoms, int n_beta, double scf_tol,
                                        int scf_max_iter, double h);

typedef struct {
  int n_atoms;
  int n_modes;             /* 3*n_atoms - (5 or 6): genuine vibrations only,
                            * translation/rotation projected out */
  int n_trans_rot;         /* 5 (linear) or 6 (nonlinear); 3*n_atoms -
                            * n_modes */
  double *frequencies_cm1; /* length n_modes, ascending; negative values
                            * denote imaginary frequencies (a saddle point,
                            * not a minimum - sqrt of a negative
                            * mass-weighted-Hessian eigenvalue) */
} molecular_vib_result_t;

void molecular_vib_result_free(molecular_vib_result_t *res);

/*
 * Harmonic vibrational frequency analysis: mass-weights the given Hessian,
 * projects out the 6 (or 5, for a linear molecule) translation/rotation
 * directions via Gram-Schmidt on the mass-weighted translation/rotation vectors
 * (built from the center of mass and, for rotation, r x e_axis for each
 * Cartesian axis), diagonalizes the projected mass-weighted Hessian, and
 * converts eigenvalues to wavenumbers (cm^-1) via the standard 1 Hartree =
 * 219474.6313705 cm^-1 conversion (frequency = \sqrt(|eigenvalue in
 * Hartree / (bohr^2 * m_e)|) * that constant, sign preserved to flag imaginary
 * frequencies).
 *
 * hessian   : flat (3*n_atoms)^2 array (as returned by
 *             molecular_*_hessian_numerical), Hartree/bohr^2
 * geom      : n_atoms x 3 Cartesian geometry, bohr (same geometry the Hessian
 *             was computed at)
 * masses_amu: length n_atoms, atomic mass units (use each element's atomic
 *             weight - e.g. 1.00782503207 for H-1, not 1.0; this matters at the
 *             ~0.5% level for light atoms)
 *
 * Returns NULL on invalid input (n_atoms<=0, NULL arrays) or allocation failure
 */
molecular_vib_result_t *
molecular_vibrational_analysis(const double *hessian, int n_atoms,
                               const double masses_amu[],
                               const double geom[][3]);

#endif
