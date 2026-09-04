/*
Pyhton Matplotlib subprocess pipe
 */

#include "matplotlib_pipe.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

matplotlib_t *matplotlib_open(void) {
  matplotlib_t *mp = malloc(sizeof *mp);
  if (!mp) {
    return NULL;
  }

  // -u : unbuffered stdout/stderr from child.
  // -  : read script from stdin
  mp->pipe = popen("python3 -u -", "w");
  if (!mp->pipe) {
    free(mp);

    return NULL;
  }

  matplotlib_cmd(mp, "import matplotlib");
  matplotlib_cmd(mp, "matplotlib.use('Agg')");
  matplotlib_cmd(mp, "import matplotlib.pyplot as plt");

  return mp;
}

void matplotlib_close(matplotlib_t *mp) {
  if (!mp) {
    return;
  }

  if (mp->pipe) {
    pclose(mp->pipe);
  }

  free(mp);
}

void matplotlib_cmd(matplotlib_t *mp, const char *cmd, ...) {
  if (!mp || !mp->pipe) {
    return;
  }

  va_list args;
  va_start(args, cmd);
  // NOLINTNEXTLINE(clang-analyzer-valist.Uninitialized)
  vfprintf(mp->pipe, cmd, args);
  va_end(args);
  fprintf(mp->pipe, "\n");
  fflush(mp->pipe);
}
