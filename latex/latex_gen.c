#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Render a LaTeX math expression to a PNG file.
   NOTE: Requires pdflatex and pdftoppm (poppler-utils) installed.
   TODO: Add lualatex support
 */
int latex_render_to_png(const char *expr, const char *outpath) {
  FILE *f = fopen("/tmp/qmc_eq.tex", "w");
  if (!f)
    return -1;
  fprintf(f,
          "\\documentclass[12pt,preview]{standalone}\n"
          "\\usepackage{amsmath,amssymb,physics,bm}\n"
          "\\begin{document}\n"
          "\\[ %s \\]\n"
          "\\end{document}\n",
          expr);
  fclose(f);

  int r = system(
      "cd /tmp && pdflatex -interaction=batchmode qmc_eq.tex > /dev/null 2>&1");
  if (r != 0)
    return -2;

  char cmd[512];
  snprintf(cmd, sizeof cmd,
           "pdftoppm -r 200 -png /tmp/qmc_eq.pdf /tmp/qmc_eq_out "
           "&& mv /tmp/qmc_eq_out-1.png %s",
           outpath);
  return system(cmd);
}

// Write a full LaTeX document containing a figure + equation
int latex_write_figure(const char *texpath, const char *figure_pdf,
                       const char *caption, const char *equation) {
  FILE *f = fopen(texpath, "w");
  if (!f)
    return -1;
  fprintf(f,
          "\\documentclass{article}\n"
          "\\usepackage{amsmath,amssymb,physics,graphicx,bm}\n"
          "\\begin{document}\n"
          "\\begin{figure}[h]\\centering\n"
          "  \\includegraphics[width=0.8\\textwidth]{%s}\n"
          "  \\caption{%s}\n"
          "\\end{figure}\n"
          "\\[ %s \\]\n"
          "\\end{document}\n",
          figure_pdf, caption, equation);
  fclose(f);
  return 0;
}
