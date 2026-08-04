#ifndef QMC_QUANTUM_INFO_H
#define QMC_QUANTUM_INFO_H

#include "../core/complex.h"

/*
 * Quantum information protocols : quantum teleportation, superdense coding, and
 * CHSH form of Bell's inequality.
 */

typedef struct {
  int m1, m2; // two classical bits Alice's Bell measurement produces

  // Bob's qubit after correction \alpha * |0> + \beta * |1>, should exactly
  // equal teleported (\alpha, \beta) up to global phase
  complex_t bob_alpha;
  complex_t bob_beta;

} qi_teleport_result_t;

/*
 * Quantum teleportation: teleports single-qubit state \alpha * |0> + \beta *
 * |1> (not required to be pre-normalized; the state is normalized internally)
 * from Alice to Bob via a shared Bell pair and two classical bits.
 *
 * NOTE: u1, u2 in [0,1) select which of the 4 equally-likely (25% each)
 * Bell-measurement outcomes occurs, using same caller-supplied-randomness
 * convention as qstate_measure_qubit (so callers can drive this with an RNG for
 * a physical simulation, or with fixed values to force/inspect a specific
 * outcome deterministically; teleportation succeeds for every outcome after
 * Bob's correction, which is exactly what makes the protocol work.
 */
qi_teleport_result_t qi_teleport(complex_t alpha, complex_t beta, double u1,
                                 double u2);

typedef struct {
  int decoded_bit1, decoded_bit2;
} qi_superdense_result_t;

/*
 * Superdense coding: encodes two classical bits (bit1, bit2) onto Alice's half
 * of shared Bell pair (I/X/Z/XZ depending on bit pattern), then decodes them
 * exactly from Bob's side after he receives Alice's qubit.
 *
 * NOTE:  Fully deterministic (100% success by construction. After encoding,
 * two-qubit state is exactly one of four mutually orthogonal Bell states, and
 * Bob's decoding circuit rotates each to a distinct, unambiguous computational
 * basis state), so no measurement randomness is needed here.*/
qi_superdense_result_t qi_superdense(int bit1, int bit2);

/*
 * Exact quantum-mechanical correlator :
 *   <\Phi+| A(\theta) (x) B(phi) |\Phi+>
 * for Bell state |\Phi+> = (|00> + |11>) / \sqrt2
 * Where
 *  A(\theta) = \cos(\theta) * \sigma_z + \sin(\theta) * \sigma_x
 * (a spin measurement in the x-z plane at angle \theta from the z-axis), and
 * B(\phi) likewise on second qubit. Computed exactly (dense 4x4 matrix
 * expectation value), not via measurement sampling. Known closed form is
 * \cos(\theta - \phi).
 */
double qi_chsh_correlator(double theta, double phi);

/*
 * CHSH combination S = E(a,b) - E(a,b') + E(a',b) + E(a',b'), using
 * qi_chsh_correlator for each term. Local hidden-variable (classical) theories
 * are bounded by |S| <= 2; quantum mechanics allows up to Tsirelson bound 2 *
 * \sqrt(2) =~ 2.8284, saturated by the Bell state with angles a=0, a'=\pi/2,
 * b=\pi/4, b'=3*\pi/4.
 */
double qi_chsh_S(double a, double a_prime, double b, double b_prime);

#endif
