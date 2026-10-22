/* Downstream-consumer smoke test.
 *
 * Built by tests/consumer/run.sh against an *installed* QMC (not the source
 * tree), using only `pkg-config --cflags --libs qmc` and <qmc/...> includes,
 * exactly as an external project would.
 */
#include <math.h>
#include <stdio.h>

#include <qmc/core/matrix.h>
#include <qmc/core/random.h>
#include <qmc/physics/hydrogen.h>

int main(void) {
  // physics: hydrogen ground state is -13.6 eV = -2.17896e-18 J
  double e1 = hydrogen_energy_level(1);
  if (fabs(e1 - (-2.17896e-18)) > 1e-22) {
    fprintf(stderr, "hydrogen_energy_level(1) = %g\n", e1);
    return 1;
  }

  // core: RNG is deterministic for a given seed
  rng_state_t a, b;
  rng_seed(&a, 12345);
  rng_seed(&b, 12345);

  for (int i = 0; i < 8; i++) {
    double x = rng_uniform(&a);
    if (x != rng_uniform(&b) || x < 0.0 || x >= 1.0) {
      fprintf(stderr, "rng mismatch at draw %d\n", i);
      return 1;
    }
  }

  // core: allocation round trip (zero-initialised)
  cmatrix_t *m = cmatrix_alloc(2, 2);
  if (!m) {
    fprintf(stderr, "cmatrix_alloc failed\n");
    return 1;
  }

  cmatrix_free(m);

  puts("qmc consumer OK");

  return 0;
}
