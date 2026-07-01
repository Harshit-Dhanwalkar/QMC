#ifndef QMC_CONSTANTS_H
#define QMC_CONSTANTS_H

#include <math.h>

/* Fundamental constants (SI units) */
#define HBAR 1.054571817e-34        /* Planck constant / 2π (J·s) */
#define M_ELECTRON 9.1093837015e-31 /* Electron mass (kg) */
#define M_PROTON 1.67262192369e-27  /* Proton mass (kg) */
#define E_CHARGE 1.602176634e-19    /* Elementary charge (C) */
#define K_COULOMB 8.9875517923e9    /* Coulomb constant (N·m²/C²) */
#define EPSILON_0 8.8541878128e-12  /* Vacuum permittivity (F/m) */
#define C_LIGHT 299792458.0         /* Speed of light (m/s) */

/* Atomic units (useful for QM calculations) */
#define AU_LENGTH 5.29177210903e-11   /* Bohr radius (m) */
#define AU_ENERGY 4.3597447222071e-18 /* Hartree (J) */
#define AU_TIME 2.4188843265857e-17   /* Atomic time unit (s) */

/* Convenient combinations */
#define HBAR_SQ (HBAR * HBAR)
#define HBAR_2M (HBAR_SQ / (2.0 * M_ELECTRON))

#endif
