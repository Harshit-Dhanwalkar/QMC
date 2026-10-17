# Landau–Zener Sweeps

A two-level system whose detuning is swept through resonance, on the Bloch sphere. Natural units ($\hbar = 1$). Implemented in `physics/landau_zener.h`; the lab-frame and beyond-RWA cases are in [Driven Two-Level Systems](driven.md).

## The Hamiltonian

$$
H(t) = \tfrac12\left[\Omega\,\sigma_x + \delta(t)\,\sigma_z\right], \qquad
\frac{d\mathbf v}{dt} = \mathbf b \times \mathbf v, \quad \mathbf b = (\Omega, 0, \delta)
$$

The detuning is swept linearly, $\delta(t) = \delta_0 + \text{rate}\cdot t$. The two diabatic levels cross at $\delta = 0$; the adiabatic levels $\pm\tfrac12\sqrt{\Omega^2 + \delta^2}$ avoid each other with a minimum gap $\Omega$.

## The Landau–Zener formula

A sweep from $\delta = -\infty$ to $+\infty$ that starts in the lower adiabatic state ends in the upper one, a diabatic jump in which the state keeps its diabatic character, with probability

$$
P_\text{LZ} = \exp\!\left(-\frac{\pi\,\Omega^2}{2\,|\text{rate}|}\right)
$$

A slow sweep follows the adiabatic level ($P_\text{LZ} \to 0$); a fast one jumps ($P_\text{LZ} \to 1$). With no coupling the jump is certain.

## Stückelberg interference

Sweeping out and back makes two paths through the crossing interfere. The upper-level population after two passes oscillates, as the sweep amplitude (and so the accumulated phase) changes, between about $0$ and $4P_\text{LZ}(1 - P_\text{LZ})$.

## Functions

| Function                                           | Purpose                                                                                      |
| -------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `lz_probability(omega, rate)`                      | the closed form $\exp(-\pi\Omega^2/2\lvert\text{rate}\rvert)$                                |
| `lz_lower_state(v, omega, delta)`                  | put $\mathbf v$ in the lower adiabatic state, $\mathbf v = -\mathbf b/\lvert\mathbf b\rvert$ |
| `lz_upper_population(v, omega, delta)`             | $(1 + \mathbf v\cdot\mathbf b/\lvert\mathbf b\rvert)/2$                                      |
| `lz_advance(v, omega, delta0, rate, duration, dt)` | RK4 with linearly moving detuning                                                            |
| `lz_sweep(v, omega, rate, amp, dt, passes)`        | zigzag the detuning between $-\text{amp}$ and $+\text{amp}$                                  |

The step is capped at $0.03/\sqrt{\Omega^2 + \delta_\text{max}^2}$, which keeps $\lvert\mathbf v\rvert$ conserved to about $10^{-9}$ even though RK4 is not exactly norm preserving.

## Example

```c
const double omega = 1.0, amp = 12.0;
double v[3];

lz_lower_state(v, omega, -amp);            /* start on the lower level   */
lz_sweep(v, omega, 2.0, amp, 0.01, 1);     /* one pass, rate = 2         */
printf("%.5f vs %.5f\n", lz_upper_population(v, omega, amp),
       lz_probability(omega, 2.0));
```

For $\Omega = 1$ and a finite sweep from $-12$ to $+12$:

| rate | simulated | $\exp(-\pi\Omega^2/2\,\text{rate})$ |
| ---- | --------- | ----------------------------------- |
| 0.5  | 0.04319   | 0.04321                             |
| 2.0  | 0.45555   | 0.45594                             |
| 8.0  | 0.82066   | 0.82172                             |

The small differences come from the sweep covering a finite range rather than $\pm\infty$. The full program is `examples/eg_81_landau_zener.c`.

## Verification

`tests/test_landau_zener.c` checks the closed form against RK4 integration for nine $(\Omega, \text{rate})$ pairs, the adiabatic, sudden and no-coupling limits with $\lvert\mathbf v\rvert$ conserved, that reversing the sweep gives the same probability, the Stückelberg oscillation between about $0$ and $4P(1-P)$, and that stepping with `lz_advance` agrees with one `lz_sweep`.
