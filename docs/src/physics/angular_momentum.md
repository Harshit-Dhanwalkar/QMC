# Angular Momentum & Spin

Quantum mechanical treatment of angular momentum including orbital angular momentum, spin, and their coupling.

## Orbital Angular Momentum

The orbital angular momentum operator:

$$
\mathbf{L} = \mathbf{r} \times \mathbf{p}
$$

In spherical coordinates:

$$
L_z = -i\hbar\frac{\partial}{\partial\phi}
$$

$$
L^2 = -\hbar^2\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}\right]
$$

### Eigenvalues

$$
L^2|l,m\rangle = \hbar^2 l(l+1)|l,m\rangle
$$

$$
L_z|l,m\rangle = \hbar m|l,m\rangle
$$

where $l = 0, 1, 2, \ldots$ and $m = -l, -l+1, \ldots, l$.

## Spherical Harmonics

The spherical harmonics are the eigenfunctions of $L^2$ and $L_z$:

$$
Y\_{lm}(\theta, \phi) = \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}} P_l^m(\cos\theta) e^{i m \phi}
$$

### Implementation

```c
complex_t spherical_harmonic(int l, int m, double theta, double phi) {
    double prefactor = sqrt((2*l+1)/(4*M_PI) *
                           factorial(l-m)/factorial(l+m));
    double Plm = associated_legendre(l, m, cos(theta));
    return c_scale(c_exp(c_imag(m * phi)), prefactor * Plm);
}
```

## Spin Angular Momentum

Spin is an intrinsic angular momentum with $s=1/2$ for electrons:

$$
S_z|\uparrow\rangle = +\frac{\hbar}{2}|\uparrow\rangle, \quad
S_z|\downarrow\rangle = -\frac{\hbar}{2}|\downarrow\rangle
$$

The Pauli matrices represent the spin operators:

$$
\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

## Spin-Orbit Coupling

The spin-orbit interaction:

$$
H_{SO} = \frac{1}{2m^2 c^2}\frac{1}{r}\frac{dV}{dr}\mathbf{L}\cdot\mathbf{S}
$$

This splits energy levels with different total $j = l \pm 1/2$.

### Implementation

```c
// Pauli matrices
typedef struct {
    complex_t data[2][2];
} spin_matrix_t;

spin_matrix_t sigma_x = {{
    {{0,0}, {1,0}},
    {{1,0}, {0,0}}
}};

spin_matrix_t sigma_y = {{
    {{0,0}, {0,-1}},
    {{0,1}, {0,0}}
}};

spin_matrix_t sigma_z = {{
    {{1,0}, {0,0}},
    {{0,0}, {-1,0}}
}};

// Spin states
typedef struct {
    complex_t up;
    complex_t down;
} spinor_t;

spinor_t spin_up() {
    return (spinor_t){{1,0}, {0,0}};
}

spinor_t spin_down() {
    return (spinor_t){{0,0}, {1,0}};
}

// Expectation value of spin
double spin_expectation_z(spinor_t s) {
    return HBAR/2 * (c_abs2(s.up) - c_abs2(s.down));
}
```

## Coupled Angular Momentum

For two angular momenta $J_1$ and $J_2$, the total angular momentum is:

$$
J=J_1 + J_2
$$

The possible values of _j_ are:

$$
j=|j_1−j_2|,|j_1−j_2| + 1, \cdots, j_1 + j_2
$$

The Clebsch-Gordan coefficients connect the uncoupled and coupled bases.
