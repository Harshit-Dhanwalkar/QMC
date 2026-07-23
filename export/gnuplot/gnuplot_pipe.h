#ifndef QMC_GNUPLOT_PIPE_H
#define QMC_GNUPLOT_PIPE_H

#include <stdio.h>

typedef struct {
  FILE *pipe;
} gnuplot_t;

/* Lifetime */
gnuplot_t *gnuplot_open(void);
void gnuplot_close(gnuplot_t *gp);

/* Commands */
void gnuplot_cmd(gnuplot_t *gp, const char *cmd, ...);

/* Plotting */
// void gnuplot_plot_file(gnuplot_t *gp, const char *filename,
//                          const char *title, const char *xlabel,
//                          const char *ylabel, int using_cols);
//
// void gnuplot_plot_wavefunction(gnuplot_t *gp, const char *wavefunction_file,
//                                  const char *potential_file);
//
// void gnuplot_multiplot(gnuplot_t *gp, int rows, int cols);
//
// void gnuplot_end_multiplot(gnuplot_t *gp);

#endif
