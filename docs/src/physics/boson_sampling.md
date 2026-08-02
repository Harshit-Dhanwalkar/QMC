# Boson Sampling

Transition amplitudes/probabilities for $N$ indistinguishable photons passing through an $M\times M$ linear-optical network described by a unitary $U$, from an input Fock configuration to an output Fock configuration. Builds directly on the bosonic permanent from [Identical Particles](identical_particles.md). Implemented in `physics/boson_sampling.h`.

## From permanent to transition amplitude

`bosonic_permanent_value` alone gives $\text{Perm}(U_{ST})/\sqrt{N!\prod_i s_i!}$. A genuine Fock-state transition amplitude needs normalization by **both** input and output mode occupation factorials:

$$
A = \frac{\text{Perm}(U_{ST})}{\sqrt{\prod_i s_i!\ \prod_j t_j!}}
$$

```c
complex_t boson_sampling_amplitude(const cmatrix_t *U, int M,
                                   const int *input_modes,
                                   const int *output_modes, int N);

// |amplitude|^2 - transition probability
double boson_sampling_probability(const cmatrix_t *U, int M,
                                  const int *input_modes,
                                  const int *output_modes, int N);
```

`input_modes`/`output_modes` are length-`N` arrays of mode indices ($0..M-1$), one entry per photon - a repeated value means that mode is multiply occupied (bosonic bunching). `U` is the $M\times M$ unitary describing the linear-optical network.

```c
int input_modes[2]  = {0, 1};  // one photon in each of modes 0 and 1
int output_modes[2] = {0, 1};
double P = boson_sampling_probability(U, M, input_modes, output_modes, 2);
```

## Test networks

Two manifestly-unitary networks are provided for demonstrations/tests without needing a general Haar-random-unitary generator:

```c
/* 50:50 beam-splitter, Hadamard convention:
   U = (1 / \sqrt2) * [[1,1],[1,-1]] - the Hong-Ou-Mandel test case. */
cmatrix_t *beam_splitter_50_50(void);

// MxM DFT unitary: U_jk = (1 / \sqrt(M)) * \exp(2 * \pi * i * j * k / M)
cmatrix_t *dft_unitary(int M);
```

`beam_splitter_50_50` is specifically useful for reproducing the **Hong-Ou-Mandel effect**: two indistinguishable photons entering the two input ports of a 50:50 beam splitter always exit together from the same output port (bunching), never one from each - a purely bosonic-interference result that vanishes for distinguishable particles. Computing `boson_sampling_probability` for the "one photon each output" configuration against this beam splitter should give exactly 0.

## Running the Example

```sh
./build/eg_25_boson_sampling
```
