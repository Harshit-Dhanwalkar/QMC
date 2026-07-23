#ifndef QMC_MATPLOTLIB_PIPE_H
#define QMC_MATPLOTLIB_PIPE_H

#include <stdio.h>

typedef struct {
  FILE *pipe;
} matplotlib_t;

/*
 * Opens `python3 -u -` subprocess and pipes preamble:
 *   import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
 */
matplotlib_t *matplotlib_open(void);

void matplotlib_close(matplotlib_t *mp);

void matplotlib_cmd(matplotlib_t *mp, const char *cmd, ...);

#endif
