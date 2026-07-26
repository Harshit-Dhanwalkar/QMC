// Master header for Physics
#ifndef QMC_PHYSICS_H
#define QMC_PHYSICS_H

/* Standard Libs */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Core headers */
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/fft/fft.h"
#include "../core/fft/fft2d.h"
#include "../core/fft/fft3d.h"
#include "../core/linalg/linalg.h"
#include "../core/matrix.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/ode/numerov.h"
#include "../core/special/special.h"
#include "../core/utils.h"
#include "../core/vector.h"

/* Physics modules */
#include "angular.h"
#include "central_potential.h"
#include "driven.h"
#include "fine_structure.h"
#include "hartree_fock.h"
#include "helium.h"
#include "hydrogen.h"
#include "identical.h"
#include "lindblad.h"
#include "perturbation.h"
#include "potentials.h"
#include "qubits.h"
#include "rabi.h"
#include "relativistic.h"
#include "scattering.h"
#include "schrodinger.h"
#include "boson_sampling.h"
#include "zeeman.h"
#include "soft.h"
#include "uncertainty.h"
#include "variational.h"
#include "wavefn.h"
#include "wkb.h"

#endif
