#ifndef QMC_BASIS_PARSER_H
#define QMC_BASIS_PARSER_H

/*
 * General basis-set-file parser: reads Gaussian94-format basis-set text
 * (plain-text format exported by Basis Set Exchange, (Refernce:
 * https://www.basissetexchange.org), and used by Gaussian/psi4/NWChem's
 * "Gaussian94" export option) and builds basis_function_t arrays for arbitrary
 * elements and arbitrary basis sets.
 *
 * Format (one element block per element, blocks separated by a "****" line):
 *
 *   H     0
 *   S   3   1.00
 *         3.42525091              0.15432897
 *         0.62391373              0.53532814
 *         0.16885540              0.44463454
 *   ****
 *   LI     0
 *   S   3   1.00
 *        16.1195750              0.15432897
 *         2.9362007              0.53532814
 *         0.7946505              0.44463454
 *   SP   3   1.00
 *         0.6362897             -0.09996723             0.15591627
 *         0.1478601              0.39951283             0.60768372
 *         0.0480887              0.70011547             0.39195739
 *   ****
 *
 * NOTE: Shell types S,P,D,F,G (angular momentum 0..4) are supported, each
 * expanded into its full set of (l+1)(l+2)/2 Cartesian components (e.g. P ->
 * px,py,pz; D -> the 6 Cartesian d components xx,xy,xz,yy,yz,zz). An "SP" shell
 * (shared exponents, separate S and P contraction coefficients :
 * minimal-basis-set convention used by e.g. STO-3G/3-21G/6-31G for valence s+p
 * shell) is split into one S basis_function_t and three P basis_function_t.
 *
 * Every basis_function_t returned is already run through
 * molint_normalize_contraction (own-primitive normalization + whole-contraction
 * normalization to unit self-overlap).
 */

#include "molecular_integrals.h"

/* One shell as parsed from file: n_prim primitives shared by every Cartesian
 * component of this shell. coef1 is the contraction for the shell's own
 * angular momentum (S,P,D,... or the S-half of an SP shell); coef2 is the
 * P-contraction of an SP shell, or NULL for any non-SP shell. */
typedef struct {
  char type[3]; /* "S","P","D","F","G","SP", NUL-terminated */
  int l;        /* angular momentum of coef1's function (0 for SP's S part) */
  int n_prim;
  double *exps;
  double *coef1;
  double *coef2; /* NULL unless type=="SP" (then P-part coefficients) */
} basis_shell_t;

typedef struct {
  char element[4]; /* upper-case element symbol, e.g. "H", "LI" */
  int n_shells;
  basis_shell_t *shells;
} basis_element_t;

typedef struct {
  int n_elements;
  basis_element_t *elements;
} basis_set_t;

/* Parse basis-set text already in memory (NUL-terminated).
 *
 * Returns NULL on malformed input or allocation failure.
 */
basis_set_t *basis_set_parse_string(const char *text);

/* Parse a basis-set file from disk.
 *
 * Returns NULL if the file cannot be opened/read or the contents are malformed.
 */
basis_set_t *basis_set_parse_file(const char *path);

void basis_set_free(basis_set_t *bs);

/* Case-insensitive element lookup.
 *
 *Returns NULL if not present in bs.
 */
const basis_element_t *basis_set_find_element(const basis_set_t *bs,
                                              const char *symbol);

/* Expand one element's shells, centered at `center`, into an array of
 * newly-allocated, normalized basis_function_t. Writes the array to *out
 *
 * Returns the count, or returns 0 and sets *out=NULL on allocation failure.
 */
int basis_set_build_atom(const basis_element_t *elem, const double center[3],
                         basis_function_t ***out);

/* Build the full concatenated basis for a molecule in one call.
 * symbols[i]/centers[i] give element symbol and 3D center of atom i.
 *
 * Returns total basis function count (0 on failure, including any symbol
 * missing from bs), writing the flat array to *out (ready to hand straight to
 * molecular_overlap_matrix / molecular_core_hamiltonian / etc.).
 */
int basis_set_build_molecule(const basis_set_t *bs, const char *const symbols[],
                             const double centers[][3], int n_atoms,
                             basis_function_t ***out);

/* Frees an array of basis_function_t* (and the array itself) as returned by
 * basis_set_build_atom / basis_set_build_molecule.
 */
void basis_set_free_functions(basis_function_t **funcs, int n);

#endif
