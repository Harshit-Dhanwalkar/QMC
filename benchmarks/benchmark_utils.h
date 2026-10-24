/* malloc that aborts whole MPI job on failure, so callers never see NULL */

#include <stdio.h>
#include <stdlib.h>

static void *xmalloc(size_t bytes) {
  void *p = malloc(bytes);

  // if (!p) {
  //   fprintf(stderr, "bench_openmp_parallel: out of memory (%zu bytes)\n",
  //           bytes);
  //   MPI_Abort(MPI_COMM_WORLD, 1);
  // }

  if (!p) {
    fprintf(stderr, "out of memory (%zu bytes)\n", bytes);
    exit(EXIT_FAILURE);
  }

  return p;
}
