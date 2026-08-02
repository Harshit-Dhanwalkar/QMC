# Driven Two-Level Systems

General time-dependent two-level Hamiltonians, beyond the exact closed-form [Rabi](rabi.md) solutions: Landau-Zener sweeps and lab-frame driving beyond the rotating-wave approximation (RWA). Natural units ($\hbar=m=1$). Implemented in `physics/driven.h`.

## Time-dependent parameters

```c
typedef double (*time_fn)(double t, void *params);

double time_fn_constant(double t, void *params);      /* fixed value */
double time_fn_linear_ramp(double t, void *params);   /* \alpha * t */
```

`time_fn` is the general caller-supplied hook for a time-dependent scalar (detuning, coupling, pulse envelope, chirp, ...); the two provided implementations cover the constant and linearly-swept cases directly.

## Rotating-frame evolution (general $\Delta(t)$, $\Omega(t)$)

```c
int driven_two_level_evolve(cvector_t *psi, time_fn Delta, void *delta_params,
                            time_fn Omega, void *omega_params, double t0,
                            double dt, int steps);
```

Integrates $i\,d\psi/dt = H(t)\psi$ with $H(t) = \frac12\begin{pmatrix}\Delta(t) & \Omega(t)\\ \Omega(t) & -\Delta(t)\end{pmatrix}$ via RK4 (the same rotating-frame form [Rabi](rabi.md) uses, but here $\Delta$ and $\Omega$ can vary arbitrarily in time rather than being fixed). `psi` evolves in place from `t0` to `t0 + steps*dt`.

```c
double alpha = 0.5;
driven_two_level_evolve(psi, time_fn_linear_ramp, &alpha,
                        time_fn_constant, &Omega, t0, dt, steps);
```

## Landau-Zener

For a linearly-swept detuning $\Delta(t) = \alpha t$ through a fixed coupling $\Omega$ (an avoided crossing), the closed-form probability of ending up in the diabatic (bare) state rather than adiabatically following the instantaneous ground state is:

```c
/* P_LZ = \exp(-\pi * \Omega^2 / (2 * \alpha)) */
double landau_zener_probability(double Omega, double alpha);
```

$$
P_{\text{LZ}} = \exp\left(-\frac{\pi\Omega^2}{2\alpha}\right)
$$

`alpha > 0` is the sweep rate. Fast sweeps ($\alpha$ large) push $P_{\text{LZ}} \to 1$ (diabatic limit - no time to react to the crossing); slow sweeps ($\alpha$ small) push $P_{\text{LZ}} \to 0$ (adiabatic limit - follows the instantaneous eigenstate across the crossing). This closed form is a natural validation target for `driven_two_level_evolve` with a `time_fn_linear_ramp` detuning: numerically integrate and compare the final diabatic-state population against `landau_zener_probability`.

## Lab-frame driving (beyond RWA)

```c
/* H(t) = (\omega0 / 2) * \sigma_z + \Omega0 * \cos(\omega_L * t + phase) * \sigma_x */
int driven_two_level_evolve_lab_frame(cvector_t *psi, double omega0,
                                      double Omega0, double omega_L,
                                      double phase, double t0, double dt,
                                      int steps);
```

Integrates the actual lab-frame Hamiltonian via RK4, with `omega0` the bare transition frequency, `Omega0` the physical drive amplitude, `omega_L` the drive frequency, and `phase` a drive phase offset - no rotating-wave approximation is made, so this captures counter-rotating terms that the RWA (and hence `rabi.c` / `driven_two_level_evolve` above) discards. Useful for checking how good the RWA actually is in a given coupling regime by comparing against the RWA solution at the same nominal parameters.

## Running the Example

```sh
./build/eg_21_driven
```
