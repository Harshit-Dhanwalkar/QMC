# FFT Module

The FFT module provides discrete Fourier transforms for wavefunctions, used
throughout the library to move between position and momentum representations

- most visibly in the split-operator propagator ([SOFT](../physics/soft.md))
  and in [uncertainty](../physics/uncertainty.md)-type position/momentum
  comparisons.

## 1D FFT

Defined in `core/fft/fft.h`:

```c
void fft(cvector_t *x);              // in-place forward FFT, unnormalized
void ifft(cvector_t *x);             // in-place inverse FFT, unnormalized
void fft_normalized(cvector_t *x);   // forward, unitary 1/sqrt(N) convention
void ifft_normalized(cvector_t *x);  // inverse, unitary 1/sqrt(N) convention
void fft_shift(cvector_t *x);        // move zero-frequency component to center
```

`fft`/`ifft` implement the standard iterative Cooley-Tukey radix-2 algorithm
(bit-reversal permutation, then $\log_2 N$ butterfly passes), which requires
a power-of-two length. **If `x->n` is not a power of two, `fft`/`ifft` fall
back to a direct $O(N^2)$ DFT sum instead of failing** - correct for any
length, but silently much slower for large non-power-of-two grids. Pad to a
power of two if performance matters.

`fft_normalized`/`ifft_normalized` are `fft`/`ifft` followed by a $1/\sqrt{N}$
rescale, giving the unitary DFT convention where forward and inverse are
each other's adjoint and Parseval's theorem holds without an extra factor.

`fft_shift` matches NumPy's `fftshift` convention exactly (`numpy.roll(x,
n // 2)`) for both even and odd `n` - e.g. for $n=6$: $[0,1,2,3,4,5]
\rightarrow [3,4,5,0,1,2]$; for $n=5$: $[0,1,2,3,4] \rightarrow
[3,4,0,1,2]$. Useful for plotting a momentum-space wavefunction with $k=0$
in the middle of the array instead of at the edges.

### FFTW wrapper

Defined in `core/fft/fft_wrapper.h`:

```c
void fft_wrapper(cvector_t *x);
void ifft_wrapper(cvector_t *x);
void fft_wrapper_normalized(cvector_t *x);
void ifft_wrapper_normalized(cvector_t *x);
```

Same four signatures as above, but built with `USE_FFTW` (see
[Building from Source](../setup/building.md)) they call out to FFTW's
`fftw_plan_dft_1d` (`FFTW_ESTIMATE`) instead of the built-in radix-2/naive-DFT
implementation - useful for large grids where FFTW's real speed advantage
matters. Without `USE_FFTW`, all four are thin pass-throughs to
`fft`/`ifft`/`fft_normalized`/`ifft_normalized`, so code written against the
wrapper works identically either way.

## 2D and 3D FFT

Defined in `core/fft/fft2d.h` and `core/fft/fft3d.h`:

```c
void fft2d(cvector_t *psi, int Nx, int Ny);
void ifft2d(cvector_t *psi, int Nx, int Ny);   // normalized by 1/(Nx*Ny)

void fft3d(cvector_t *psi, int Nx, int Ny, int Nz);
void ifft3d(cvector_t *psi, int Nx, int Ny, int Nz);
```

Both are separable transforms built directly from the 1D `fft`/`ifft`:
`fft2d` runs the 1D FFT along every row, then along every column; `fft3d`
adds a third pass along the depth axis. `psi` is a flattened `cvector_t` of
length `Nx*Ny` (or `Nx*Ny*Nz`) in row-major order - the same layout used by
[SOFT](../physics/soft.md)'s 2D/3D split-operator propagator, which calls
these directly rather than going through `fft_wrapper`.

## Position -> Momentum Convenience Function

Defined in `core/utils.h`:

```c
cvector_t *position_to_momentum(const cvector_t *psi_x, double deltax);
```

Copies `psi_x`, applies the raw (unnormalized) forward `fft`, then rescales
by $\Delta x / \sqrt{2\pi}$ so that

$$
\int |\psi(x)|^2\,dx = \int |\phi(k)|^2\,dk
$$

is preserved - i.e. the output is a properly normalized momentum-space
wavefunction on a grid $k_n = 2\pi n / (N \Delta x)$, $n = 0 \ldots N-1$,
_not_ run through `fft_shift` (frequencies wrap around past the Nyquist
point rather than being centered). Returns a newly allocated `cvector_t`;
the input is left untouched.

## Example

```c
int N = 256;
cvector_t *psi_x = /* position-space wavefunction on a grid of spacing dx */;

cvector_t *psi_k = position_to_momentum(psi_x, dx);
fft_shift(psi_k); // so k=0 sits in the middle of the array for plotting

/* ... use psi_k, e.g. to compute <p> or <p^2> ... */

cvector_free(psi_k);
```
