#ifndef QMC_CONFIG_H
#define QMC_CONFIG_H

/* Plot backend (informational only)
 * Backend (GR/GNUPLOT/MATPLOTLIB/NONE) is selected at build time by Makefile,
 * which compiles exactly one of export/plot_gr.c plot_gnuplot.c,
 * plot_matplotlib.c or plot_none.c.
 */
#ifndef QMC_PLOT_BACKEND_NAME
#define QMC_PLOT_BACKEND_NAME "UNKNOWN (built outside Makefile)"
#endif

/* Complex Absorbing Potential (CAP) defaults
 * For reasonable CAP without picking width/eta/power (see
 * core/ode/crank_nicolson.c: cap_build_monomial).
 */
#define QMC_CAP_DEFAULT_WIDTH 4.0 /* absorbing layer width at each edge */
#define QMC_CAP_DEFAULT_ETA 5.0   /* absorption strength at outermost point */
#define QMC_CAP_DEFAULT_POWER 2   /* monomial ramp power (2 or 3 typical) */

/* Set to 1 to enable extra diagnostic printf's in code, checks flag */
#ifndef QMC_VERBOSE
#define QMC_VERBOSE 0
#endif

#endif /* QMC_CONFIG_H */
