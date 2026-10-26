/*
 * LaTeX Report Generation
 */

#include "../latex/latex_gen.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > LaTeX Report Generation\n");

  printf("  Step 1: inline and display math snippets\n\n");
  char *inline_eq = latex_inline_math("E = -1.116714\\,\\text{Hartree}");
  char *display_eq =
      latex_display_math("\\hat{H}\\psi = \\left[-\\frac{1}{2}\\nabla^2 + "
                         "V(r)\\right]\\psi = E\\psi");
  printf("  Inline:  %s\n", inline_eq);
  printf("  Display: %s\n\n", display_eq);

  printf("  Step 2: a results table\n\n");
  const char *row0[3] = {"Method", "Basis", "Energy (Hartree)"};
  const char *row1[3] = {"RHF", "STO-3G", "-1.116714"};
  const char *row2[3] = {"CCSD", "STO-3G", "-1.137276"};
  const char *const *table_data[3] = {row0, row1, row2};

  const char *table_path = "/tmp/qmc_report_table.tex";
  latex_generate_table(table_path, table_data, 3, 3, "|l|l|r|",
                       "H$_2$ energies at different levels of theory",
                       "tab:h2_energies", "htbp");
  printf("  Table written to %s\n\n", table_path);

  printf("  Step 3: a matrix (H2/STO-3G overlap matrix)\n\n");
  const char *mrow0[2] = {"1.000000", "0.659318"};
  const char *mrow1[2] = {"0.659318", "1.000000"};
  const char *const *matrix_data[2] = {mrow0, mrow1};
  char *S_matrix = latex_matrix("pmatrix", matrix_data, 2, 2);
  printf("  S = %s\n\n", S_matrix);

  printf("  Step 4: assemble a full article from pieces above\n\n");
  char section1[512];
  snprintf(section1, sizeof(section1),
           "\\section{Results}\n"
           "RHF ground-state energy is %s. Molecular Hamiltonian is\n%s\n"
           "Table~\\ref{tab:h2_energies} compares methods and overlap matrix "
           "is\n\\[ S = %s \\]",
           inline_eq, display_eq, S_matrix);

  // Build a list of key takeaways
  char *energy_fmt = latex_format_double(-1.116714, 6);
  char bullet1[128];
  snprintf(bullet1, sizeof(bullet1), "RHF total energy: %s", energy_fmt);
  const char *bullets[] = {bullet1, "CCSD shows significant correlation energy",
                           "STO-3G is a minimal basis set", NULL};

  char list_tex[512];
  latex_string_t items_s;
  latex_string_init(&items_s);
  latex_string_append(&items_s, "\\begin{itemize}\n");
  for (int i = 0; bullets[i]; i++) {
    latex_string_appendf(&items_s, "  \\item %s\n", bullets[i]);
  }
  latex_string_append(&items_s, "\\end{itemize}\n");
  snprintf(list_tex, sizeof(list_tex), "%s", items_s.data);
  latex_string_free(&items_s);

  char *verbatim_snippet = "\\begin{verbatim}\n"
                           "Energy convergence log:\n"
                           "Iteration 1: -1.1000\n"
                           "Iteration 2: -1.1167\n"
                           "\\end{verbatim}\n";

  const char *sections[] = {section1, list_tex, verbatim_snippet, NULL};

  const char *article_path = "/tmp/qmc_report_article.tex";
  int rc = latex_generate_article(
      article_path, "H$_2$ Electronic Structure Report", "QMC", NULL,
      "A short automatically-generated summary of H$_2$/STO-3G "
      "calculations.",
      sections);
  printf("  Article written to %s (rc=%d)\n\n", article_path, rc);

  printf("  Step 5: compile to PDF and render one equation to PNG\n\n");
  if (!latex_tools_available()) {
    printf("  pdflatex/pdftoppm not found on PATH in this environment: "
           "skipping actual compilation. Generated .tex files above are still "
           "valid LaTeX source; you can compile them yourself with a LaTeX "
           "distribution installed.\n");
  } else {
    // NOTE: latex_generate_article's output isn't a self-contained standalone
    // equation, so compile it with latex_compile directly (working in /tmp,
    // where file already lives)
    rc = latex_compile(article_path, NULL, "/tmp");
    printf("  latex_compile(%s) -> %s\n", article_path,
           rc == 0 ? "success" : "failed");

    const char *png_path = "/tmp/qmc_report_equation.png";
    rc = latex_render_to_png("\\hat{H}\\psi = \\left[-\\frac{1}{2}\\nabla^2 + "
                             "V(r)\\right]\\psi = E\\psi",
                             png_path);
    printf("  latex_render_to_png -> %s (%s)\n", rc == 0 ? "success" : "failed",
           png_path);

    if (rc == 0) {
      remove(png_path);
    }
  }

  printf("\n  Cleaning up temp files created during this run...\n");
  remove(table_path);
  remove(article_path);
  remove("/tmp/qmc_report_article.pdf");
  remove("/tmp/qmc_report_article.aux");
  remove("/tmp/qmc_report_article.log");

  free(inline_eq);
  free(display_eq);
  free(S_matrix);
  free(energy_fmt);

  return 0;
}
