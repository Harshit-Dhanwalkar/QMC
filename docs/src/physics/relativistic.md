# Relativistic Quantum Mechanics

Relativistic quantum mechanics extends quantum theory to particles moving at speeds comparable to the speed of light.

## Klein-Gordon Equation

The relativistic wave equation for spin-0 particles:

$$
\left(\Box + \frac{m^2c^2}{\hbar^2}\right)\psi = 0
$$

Where $\Box = \frac{1}{c^2}\frac{\partial^2}{\partial t^2} - \nabla^2$.

In energy form: $E^2 = p^2 c^2 + m^2 c^4$.

This equation has both positive and negative energy solutions, leading to the concept of antiparticles.

## Dirac Equation

The relativistic wave equation for spin-1/2 particles:

$$
i\hbar\frac{\partial\psi}{\partial t} = \left(c\boldsymbol{\alpha}\cdot\mathbf{p} + \beta mc^2\right)\psi
$$

Where $\boldsymbol{\alpha}$ and $\beta$ are 4$\times$4 matrices:

$$
\alpha_i = \begin{pmatrix} 0 & \sigma_i \\ \sigma_i & 0 \end{pmatrix}, \quad
\beta = \begin{pmatrix} I & 0 \\ 0 & -I \end{pmatrix}
$$

The spinor has four components (two particle spin states, two antiparticle).

## Implementation

```c
// Gamma matrices
typedef struct {
    complex_t data[4][4];
} gamma_matrix_t;

gamma_matrix_t gamma_0 = {{
    {{1,0},{0,0},{0,0},{0,0}},
    {{0,0},{1,0},{0,0},{0,0}},
    {{0,0},{0,0},{-1,0},{0,0}},
    {{0,0},{0,0},{0,0},{-1,0}}
}};

// Dirac spinor
typedef struct {
    complex_t u1; // particle spin up
    complex_t u2; // particle spin down
    complex_t v1; // antiparticle spin up
    complex_t v2; // antiparticle spin down
} dirac_spinor_t;

dirac_spinor_t dirac_plane_wave(double E, double p, double m, double c) {
    dirac_spinor_t u = {0};
    double N = 1.0 / sqrt(2.0 * E * (E + m*c*c));
    u.u1 = c_real(N * (E + m*c*c));
    u.u2 = c_real(N * 0.0);
    u.v1 = c_real(N * c * p);
    u.v2 = c_real(N * 0.0);
    return u;
}
```

## Non-Relativistic Limit

In the limit $v \ll c$, the Dirac equation reduces to the Schr$\"{o}$dinger equation with spin-orbit coupling.

### Applications

1. Hydrogen Fine Structure

The Dirac equation gives the fine structure of hydrogen:

TODO:

$$
E_{nj} = mc^2 \left[ 1 + \left( \frac{\alpha}{n - (j + 1/2) + \sqrt{(j + 1/2)^2 - \alpha^2}} \right)^2 \right]^{-1/2}
$$

Where $\alpha$ is the fine structure constant.
