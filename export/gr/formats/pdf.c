/*
 * PDF output
 * standard print  page-size presets (A4, Letter)

 * NOTE: NO implementation for multi-page PDF composition.
 */

#include "../gr_plot.h"
#include "../gr_plot_internal.h"
#include <stdlib.h>
#include <stdio.h>

int gr_plot_pdf(const char *filename, const double *x, const double **ys,
                int n_series, int n_pts, const char **labels,
                const gr_plot_opt_t *opts) {
  char path[512];
  gr_build_sized_path(path, sizeof path, filename, ".pdf", opts, 1000, 750);
  setenv("GKS_WSTYPE", "nul", 0);

  return gr_emit_single_page(path, x, ys, n_series, n_pts, labels, opts);
}

typedef enum {
  GR_PDF_PAGE_DEFAULT, // opts->width/height or 1000x750
  GR_PDF_PAGE_A4,      // 595 x 842 pt
  GR_PDF_PAGE_LETTER,  // 612 x 792 pt
} gr_pdf_page_size_t;

int gr_plot_pdf_sized(const char *filename, const double *x, const double **ys,
                      int n_series, int n_pts, const char **labels,
                      const gr_plot_opt_t *opts, gr_pdf_page_size_t page_size) {
  int w, h;
  switch (page_size) {
  case GR_PDF_PAGE_A4:
    w = 595;
    h = 842;
    break;
  case GR_PDF_PAGE_LETTER:
    w = 612;
    h = 792;
    break;
  case GR_PDF_PAGE_DEFAULT:
  default:
    w = (opts && opts->width > 0) ? opts->width : 1000;
    h = (opts && opts->height > 0) ? opts->height : 750;
    break;
  }

  char path[512];
  snprintf(path, sizeof path, "%s/%s.pdf|%dx%d", QMC_OUTPUT_DIR, filename, w,
           h);
  setenv("GKS_WSTYPE", "nul", 0);

  return gr_emit_single_page(path, x, ys, n_series, n_pts, labels, opts);
}
