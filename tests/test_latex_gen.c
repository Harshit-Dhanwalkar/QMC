/*
 * Test: latex/latex_gen.c.
 *
 *  1. Pure string-generation functions (latex_inline_math, latex_display_math,
 *     latex_matrix, and the generate-article/table and write-figure file
 *     writers) need no external LaTeX toolchain and always run.
 *  2. Functions that actually invoke pdflatex/pdftoppm
 *     (latex_render_to_png/pdf, latex_compile) only run if
 *     latex_tools_available() reports the tools are present : this project
 *     builds and runs its full test suite fine on systems with no LaTeX
 *     distribution installed (latex_gen.c has no LaTeX-time dependency, only a
 *     run-time one), so unconditionally requiring pdflatex here would make the
 *     test suite non-portable for no reason.
 */

#include "../latex/latex_gen.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static int file_exists(const char *path) {
  struct stat st;

  return stat(path, &st) == 0;
}

static long file_size(const char *path) {
  struct stat st;
  if (stat(path, &st) != 0) {
    return -1;
  }

  return (long)st.st_size;
}

static void test_inline_display_math(void) {
  printf("Test: latex_inline_math / latex_display_math\n");

  char *inline_res = latex_inline_math("x^2+1");
  check(inline_res != NULL, "latex_inline_math allocates a result");
  check(inline_res && strcmp(inline_res, "$x^2+1$") == 0,
        "latex_inline_math wraps in $...$");

  free(inline_res);

  char *display_res = latex_display_math("\\int_0^1 f(x)\\,dx");
  check(display_res != NULL, "latex_display_math allocates a result");
  check(display_res && strcmp(display_res, "\\[ \\int_0^1 f(x)\\,dx \\]") == 0,
        "latex_display_math wraps in \\[ ... \\]");

  free(display_res);

  check(latex_inline_math(NULL) == NULL, "NULL input returns NULL cleanly");
  check(latex_display_math(NULL) == NULL, "NULL input returns NULL cleanly");
}

static void test_matrix_generation(void) {
  printf("Test: latex_matrix, including a cell longer than the old "
         "implementation's rough 10-char/cell size estimate");

  const char *row0[2] = {"1", "0"};
  const char *row1[2] = {"0", "1"};
  const char *const *data[2] = {row0, row1};

  char *m = latex_matrix("pmatrix", data, 2, 2);
  check(m != NULL, "latex_matrix allocates a result");
  check(m && strcmp(m, "\\begin{pmatrix}1 & 0 \\\\ 0 & 1\\end{pmatrix}") == 0,
        "2x2 identity matrix renders correctly");

  free(m);

  const char *long_cell = "\\frac{\\partial^2 \\psi}{\\partial x^2}";
  const char *lrow0[1] = {long_cell};
  const char *const *ldata[1] = {lrow0};
  char *lm = latex_matrix("bmatrix", ldata, 1, 1);
  check(lm != NULL, "latex_matrix with a long cell allocates a result");
  check(lm && strstr(lm, long_cell) != NULL,
        "long cell content is NOT truncated");

  free(lm);

  check(latex_matrix("bmatrix", data, 0, 2) == NULL,
        "rows<=0 returns NULL cleanly");
  check(latex_matrix("bmatrix", NULL, 2, 2) == NULL,
        "NULL data returns NULL cleanly");

  char *default_type = latex_matrix(NULL, data, 2, 2);
  check(default_type && strstr(default_type, "bmatrix") != NULL,
        "NULL type defaults to bmatrix");

  free(default_type);
}

static void test_generate_table(void) {
  printf("Test: latex_generate_table writes a well-formed .tex file\n");

  const char *row0[2] = {"a", "b"};
  const char *row1[2] = {"c", NULL}; // NULL cell -> empty
  const char *const *data[2] = {row0, row1};

  const char *path = "/tmp/qmc_test_table.tex";
  int rc = latex_generate_table(path, data, 2, 2, NULL, "A caption", "tab:test",
                                "htbp");
  check(rc == 0, "latex_generate_table returns 0 on success");
  check(file_exists(path), "output file was created");

  FILE *f = fopen(path, "r");
  check(f != NULL, "output file can be reopened");
  if (f) {
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';

    check(strstr(buf, "\\begin{table}[htbp]") != NULL,
          "placement argument appears in output");
    check(strstr(buf, "A caption") != NULL, "caption appears in output");
    check(strstr(buf, "tab:test") != NULL, "label appears in output");
    check(strstr(buf, "a & b") != NULL, "first row renders correctly");

    fclose(f);
  }

  remove(path);

  check(latex_generate_table(path, data, 0, 2, NULL, NULL, NULL, NULL) == -1,
        "rows<=0 returns an error instead of writing a malformed file");
  check(latex_generate_table(NULL, data, 2, 2, NULL, NULL, NULL, NULL) == -1,
        "NULL texpath returns an error");
}

static void test_generate_article(void) {
  printf("Test: latex_generate_article writes a well-formed .tex file\n");

  const char *sections[] = {"\\section{Intro}\nHello.", NULL};
  const char *path = "/tmp/qmc_test_article.tex";
  int rc = latex_generate_article(path, "Title", "Author", "2026", "Abstract.",
                                  sections);
  check(rc == 0, "latex_generate_article returns 0 on success");

  FILE *f = fopen(path, "r");
  check(f != NULL, "output file can be reopened");
  if (f) {
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';

    check(strstr(buf, "\\title{Title}") != NULL, "title appears");
    check(strstr(buf, "\\begin{abstract}") != NULL, "abstract appears");
    check(strstr(buf, "\\section{Intro}") != NULL, "section body appears");

    fclose(f);
  }

  remove(path);
}

static void test_write_figure(void) {
  printf("Test: latex_write_figure writes a well-formed .tex file, even "
         "with NULL optional fields\n");

  const char *path = "/tmp/qmc_test_figure.tex";
  int rc = latex_write_figure(path, "plot.pdf", "A caption", "E = mc^2");

  check(rc == 0, "returns 0 on success");
  check(file_exists(path), "output file created");

  remove(path);

  /* NULL optional args should not crash and should still produce valid
   * (if sparse) output rather than dereferencing NULL in fprintf's %s. */
  rc = latex_write_figure(path, NULL, NULL, NULL);
  check(rc == 0, "NULL optional args handled without crashing");

  remove(path);

  check(latex_write_figure(NULL, "x", "y", "z") == -1,
        "NULL texpath returns an error");
}

static void test_unsafe_paths_rejected(void) {
  printf("Test: shell-metacharacter-containing arguments are rejected");

  if (!latex_tools_available()) {
    printf("  (skipped: pdflatex/pdftoppm not available in this environment");
    return;
  }

  check(latex_render_to_png("x", "/tmp/evil; rm -rf /tmp/nonexistent") == -5,
        "outpath with a shell metacharacter is rejected");
  check(latex_compile("foo.tex; echo pwned", NULL, NULL) == -5,
        "texfile with a shell metacharacter is rejected");
  check(latex_compile("foo.tex", "pdflatex && echo pwned", NULL) == -5,
        "compiler with a shell metacharacter is rejected");
}

static void test_actual_rendering(void) {
  printf("Test: latex_render_to_pdf / latex_compile actually invoke "
         "pdflatex and produce real output\n");

  if (!latex_tools_available()) {
    printf(
        "  (skipped: pdflatex/pdftoppm not available in this environment)\n");
    return;
  }

  const char *pdf_out = "/tmp/qmc_test_render.pdf";
  int rc = latex_render_to_pdf("E=mc^2", pdf_out);

  check(rc == 0, "latex_render_to_pdf succeeds when tools are available");
  check(file_exists(pdf_out), "output PDF file was created");
  check(file_size(pdf_out) > 0, "output PDF is non-empty");

  remove(pdf_out);

  const char *png_out = "/tmp/qmc_test_render.png";
  rc = latex_render_to_png("E=mc^2", png_out);

  check(rc == 0, "latex_render_to_png succeeds when tools are available");
  check(file_exists(png_out), "output PNG file was created");
  check(file_size(png_out) > 0, "output PNG is non-empty");

  remove(png_out);

  /* Two renders in the same process must not clobber each other's
   * intermediate files (regression test for the old hardcoded
   * "/tmp/qmc_eq.*" basename collision). */
  const char *pdf_out2 = "/tmp/qmc_test_render2.pdf";
  int rc1 = latex_render_to_pdf("\\alpha+\\beta", pdf_out);
  int rc2 = latex_render_to_pdf("\\gamma+\\delta", pdf_out2);

  check(rc1 == 0 && rc2 == 0, "two back-to-back renders both succeed");
  check(file_exists(pdf_out) && file_exists(pdf_out2),
        "both output files exist independently");

  remove(pdf_out);
  remove(pdf_out2);

  latex_clean_temp();
}

int main(void) {
  test_inline_display_math();
  test_matrix_generation();
  test_generate_table();
  test_generate_article();
  test_write_figure();
  test_unsafe_paths_rejected();
  test_actual_rendering();

  if (failures == 0) {
    printf("\nAll test_latex_gen checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
