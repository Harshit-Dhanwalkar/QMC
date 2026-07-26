/*
 GNU plot pipe
*/

#include "gnuplot_pipe.h"
#include "../../core/vector.h"
#include "../plot.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

gnuplot_t *gnuplot_open(void) {
  gnuplot_t *gp = malloc(sizeof *gp);
  if (!gp) {
    return NULL;
  }

  gp->pipe = popen("gnuplot", "w");
  if (!gp->pipe) {
    free(gp);

    return NULL;
  }

  return gp;
}

void gnuplot_close(gnuplot_t *gp) {
  if (!gp) {
    return;
  }

  if (gp->pipe) {
    pclose(gp->pipe);
  }

  free(gp);
}

void gnuplot_cmd(gnuplot_t *gp, const char *cmd, ...) {
  if (!gp || !gp->pipe) {
    return;
  }

  va_list args;
  va_start(args, cmd);
  vfprintf(gp->pipe, cmd, args);
  va_end(args);
  fprintf(gp->pipe, "\n");
  fflush(gp->pipe);
}
