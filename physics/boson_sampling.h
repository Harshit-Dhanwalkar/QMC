#ifndef QMC_BOSON_SAMPLING_H
#define QMC_BOSON_SAMPLING_H

#include "../core/complex.h"
#include "../core/matrix.h"

/*
 * Boson sampling: transition amplitude/probability for N indistinguishable
 * photons through an MxM linear-optical network described by unitary U,
 * from an input Fock configuration to an output Fock configuration.
 *
 * bosonic_permanent_value alone gives Perm(U_ST) / \sqrt(N! * prod_i s_i!).
 * Photon-number Fock-state transition amplitudes need normalization by both
 * input and output mode occupation factorials:
 * A = Perm(U_ST) / \sqrt(prod_i * s_i! * prod_j t_j!)
 *
 * input_modes:  array of N input mode indices (0..M-1), one entry per
 *               photon. Repeated values = multiply-occupied input mode.
 * output_modes: array of N output mode indices (0..M-1), one entry per
 *               photon. Repeated values = multiply-occupied output mode.
 * N:            number of photons.
 * U:            MxM unitary.
 */
complex_t boson_sampling_amplitude(const cmatrix_t *U, int M,
                                   const int *input_modes,
                                   const int *output_modes, int N);

/* |amplitude|^2 - transition probability. */
double boson_sampling_probability(const cmatrix_t *U, int M,
                                  const int *input_modes,
                                  const int *output_modes, int N);

/*
 * 50:50 beam-splitter unitary (2x2, Hadamard convention):
 *   U = (1 / \sqrt2) * [[1, 1], [1, -1]]
 * The theoretical Hong-Ou-Mandel test case.
 * Returns matrix
 */
cmatrix_t *beam_splitter_50_50(void);

/*
 * MxM discrete Fourier transform matrix:
 * U_jk = (1 / \sqrt(M)) * \exp(2 * \pi * i * j * k / M)
 * Manifestly unitary; a
 * network for boson-sampling demonstrations/tests without needing a general
 * Haar-random-unitary generator. Returns matrix
 */
cmatrix_t *dft_unitary(int M);

#endif
