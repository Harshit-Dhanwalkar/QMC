# Anderson Localisation

Why a wave stops spreading in a disordered medium, on a 1D tight-binding chain. Natural units ($\hbar = 1$). Implemented in `physics/anderson.h`.

## The model

$$
H = -t \sum_j \left( |j\rangle\langle j+1| + \text{h.c.} \right) + \sum_j \epsilon_j\, |j\rangle\langle j|
$$

The chain has $n$ sites with open ends and hopping $t > 0$. The on-site energies $\epsilon_j$ are drawn independently and uniformly from $[-W/2, W/2]$, so $W$ is the disorder strength. `anderson_disorder` fills them from a deterministic generator: the same seed gives the same sequence on every platform.

## Clean chain: ballistic spreading

With $W = 0$ a packet started on site $j_0$ evolves exactly as

$$
\psi_j(\tau) = i^{\,j-j_0}\, J_{j-j_0}(2t\tau), \qquad \langle x^2\rangle = 2t^2\tau^2
$$

so the width grows linearly in time, at a front speed $2t$.

## Disordered chain: localisation

In one dimension any disorder, however weak, localises every eigenstate. The packet stops spreading and keeps an exponential profile $\exp(-|x-x_0|/\xi)$. For weak disorder the inverse localisation length (Lyapunov exponent) of the amplitude at energy $E$ inside the band $|E| < 2t$ is, to leading order,

$$
\frac{1}{\xi} = \frac{W^2}{24\,(4t^2 - E^2)}
$$

This estimate fails near the band centre $E = 0$ (the band-centre anomaly) and at the band edges. `anderson_lyapunov` computes the exponent from a transfer-matrix recursion on any set of on-site energies; use a long chain (at least $10^5$ sites) because it is self-averaging.

## Functions

| Function                                          | Purpose                                                                     |
| ------------------------------------------------- | --------------------------------------------------------------------------- |
| `anderson_disorder(eps, n, w, seed)`              | uniform on-site energies in $[-W/2, W/2]$                                   |
| `anderson_evolve(psi, eps, hop, dt, steps)`       | RK4 propagation under $H$                                                   |
| `anderson_norm`, `anderson_width`, `anderson_ipr` | $\sum\lvert\psi_j\rvert^2$, rms width in sites, inverse participation ratio |
| `anderson_energy(psi, eps, hop)`                  | $\langle\psi\vert H\vert\psi\rangle / \langle\psi\vert\psi\rangle$          |
| `anderson_lyapunov(eps, n, hop, energy)`          | inverse localisation length from transfer matrices                          |
| `anderson_lyapunov_weak(w, hop, energy)`          | the weak-disorder formula above                                             |

The inverse participation ratio is 1 for a state on one site and $1/n$ for a state spread evenly over $n$ sites, so its inverse counts the sites the state occupies. RK4 is not exactly unitary, but with $t = 1$ and $\Delta\tau = 0.01$ the norm drifts by about $5\times10^{-10}$ over 2000 steps on a clean chain.

## Example

```c
enum { N = 801 };
double eps[N];
cvector_t *psi = cvector_alloc(N);

anderson_disorder(eps, N, 3.0, 2024);       /* W = 3 */
cvector_fill(psi, c_new(0.0, 0.0));
psi->data[N / 2] = c_new(1.0, 0.0);         /* one site, in the middle */

anderson_evolve(psi, eps, 1.0, 0.01, 2000); /* \tau = 20 */
double w20 = anderson_width(psi);
anderson_evolve(psi, eps, 1.0, 0.01, 4000); /* \tau = 60 */
printf("%.2f %.2f %.4f\n", w20, anderson_width(psi), anderson_ipr(psi));
```

Running this for a clean chain ($W = 0$) and for $W = 3$ gives

| $W$ | width at $\tau = 20$ | width at $\tau = 60$ | IPR at $\tau = 60$ |
| --- | -------------------- | -------------------- | ------------------ |
| 0   | 28.28                | 84.85                | 0.0075             |
| 3   | 6.93                 | 6.77                 | 0.1408             |

The clean widths are exactly $\sqrt2\,\tau$ ($28.28$ and $84.85$); the disordered packet has stopped spreading. The full program is `examples/eg_79_anderson.c`.

## Verification

`tests/test_anderson.c` checks the Bessel-function spreading and $\langle x^2\rangle = 2t^2\tau^2$ on a clean chain with norm and energy conserved, the disorder generator (determinism, range, mean and variance $W^2/12$), that a strongly disordered packet stops spreading while a clean one keeps going, and the transfer-matrix exponent against the weak-disorder formula.
