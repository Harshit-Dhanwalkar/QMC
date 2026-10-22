/* Downstream-consumer smoke test
 *
 * Built by tests/consumer/run.sh against an *installed* QMC (not source tree),
 * using only `pkg-config --cflags --libs qmc` and <qmc/...> includes, exactly
 * as an external project would
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <qmc/qmc.h> // umbrella header

int main(void) {
  if (qmc_version_number() != QMC_VERSION_NUMBER ||
      strcmp(qmc_version(), QMC_VERSION_STRING) != 0) {
    fprintf(stderr, "header/library version mismatch: %s vs %s\n",
            QMC_VERSION_STRING, qmc_version());
    return 1;
  }
  printf("qmc version %s\n", qmc_version());

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

  // exporters: an in-memory-free smoke call through umbrella
  char ts[EXPORT_TIMESTAMP_SIZE];
  export_timestamp_utc(ts, sizeof ts);
  if (ts[4] != '-' || ts[10] != 'T') {
    fprintf(stderr, "bad timestamp: %s\n", ts);
    return 1;
  }

  puts("qmc consumer OK");

  return 0;
}
