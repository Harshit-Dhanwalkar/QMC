# Wave Functions

Properties of quantum mechanical wavefunctions, including normalization, expectation values, and basis transformations.

## Normalization

The wavefunction $\psi(x)$ must be normalized:

$$
\int_{-\infty}^{\infty} |\psi(x)|^2 dx = 1
$$

### Implementation

```c
double cvector_norm_sq(const cvector_t *v, double dx) {
    double norm = 0.0;
    for (int i = 0; i < v->n; i++) {
        norm += c_abs2(v->data[i]);
    }
    return norm * dx;
}

void cvector_normalize(cvector_t *v, double dx) {
    double norm = sqrt(cvector_norm_sq(v, dx));
    if (norm > 1e-15) {
        for (int i = 0; i < v->n; i++) {
            v->data[i] = c_scale(v->data[i], 1.0/norm);
        }
    }
}
```

## Expectation Values
For an operator $\hat{O}$ :
TODO

Position
TODO

Momentum
TODO


Properties:
- $\langle x \rangle = x_0$
- $\langle p \rangle = p_0$
- $\Delta x = \sigma$
- $\Delta p = \hbar/(2\sigma)$
