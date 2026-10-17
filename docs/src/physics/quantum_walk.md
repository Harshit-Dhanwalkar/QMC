# Quantum Walk on a Line

A discrete-time quantum walk with a two-state coin. Implemented in `physics/quantum_walk.h`.

## The walk

A walker with a coin (up, down) sits at integer sites $x = j - n/2$. One step applies the coin to every site, then shifts "up" one site to the right and "down" one site to the left:

$$
\begin{pmatrix} a' \\ b' \end{pmatrix} =
\begin{pmatrix} \cos\theta & \sin\theta \\ \sin\theta & -\cos\theta \end{pmatrix}
\begin{pmatrix} a \\ b \end{pmatrix},
\qquad a_\text{new}(x) = a'(x-1), \quad b_\text{new}(x) = b'(x+1)
$$

$\theta = \pi/4$ is the Hadamard walk. The map is unitary, so the total probability is conserved to round-off.

## Ballistic against diffusive

The classical walk with the same coin flip, probability $1/2$ per step, has $\langle x^2\rangle = t$ and a width that grows as $\sqrt t$. The quantum walk spreads ballistically. For an unbiased starting coin

$$
\langle x^2\rangle \;\to\; \left(1 - \sin\theta\right) t^2 \qquad (t \gg 1)
$$

Probability piles up near the edges $x = \pm t\cos\theta$ of the light cone (the "fronts") and is low in the middle; for the Hadamard walk the fronts sit near $\pm t/\sqrt2$. Starting in $|\text{up}\rangle$ instead biases the drift: $\langle x\rangle \to (1 - 1/\sqrt2)\,t$.

## Functions

| Function                                                | Purpose                                                                      |
| ------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `qwalk_init(up, down, coin)`                            | walker on the centre site; `QWALK_COIN_UP`, `_DOWN` or `_SYM`                |
| `qwalk_step(up, down, theta, steps)`                    | advance `steps` coin-and-shift steps                                         |
| `qwalk_probability`, `qwalk_norm`                       | $P(x_j) = \lvert a_j\rvert^2 + \lvert b_j\rvert^2$ and its sum              |
| `qwalk_mean`, `qwalk_variance`                          | mean position and variance in sites                                          |
| `qwalk_classical(prob, n, t)`                           | binomial distribution of the classical walk                                  |
| `qwalk_asymptotic_variance(theta)`                      | the limit $1 - \sin\theta$                                                   |

The arrays wrap around, so keep the number of steps below $n/2$; the walker moves at most one site per step.

## Example

```c
enum { N = 1024, T = 400 };
cvector_t *up = cvector_alloc(N), *down = cvector_alloc(N);

qwalk_init(up, down, QWALK_COIN_SYM);       /* symmetric start */
qwalk_step(up, down, M_PI / 4.0, T);        /* Hadamard coin   */

double var = qwalk_variance(up, down);
printf("%.5f vs %.5f\n", var / (T * (double)T),
       qwalk_asymptotic_variance(M_PI / 4.0));
```

After $t = 400$ Hadamard steps the norm is $1$ to 12 digits and $\langle x^2\rangle/t^2 = 0.29290$, against the limit $1 - \sin(\pi/4) = 0.29289$. The classical walk has $\langle x^2\rangle = 400$, a width of $20$ sites, while the quantum walk's width is about $0.54\,t \approx 216$ sites. The full program is `examples/eg_80_quantum_walk.c`.

## Verification

`tests/test_quantum_walk.c` checks one Hadamard step from $|\text{up}\rangle$ ($P(\pm1) = 1/2$), unitarity over hundreds of steps for several coin angles, $\langle x^2\rangle/t^2 \to 1 - \sin\theta$ for three angles, the symmetry and light-cone structure of the distribution (peaks near $\pm t/\sqrt2$, $\lvert x\rvert \le t$, low centre), the drift from a biased start, and the classical binomial probabilities.
