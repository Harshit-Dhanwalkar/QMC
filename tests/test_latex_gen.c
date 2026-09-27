/*
 * Test: latex/latex_gen.c.
 *
 *  1. Pure string-generation functions (latex_inline_math, latex_display_math,
 *     latex_matrix, and generate-article/table and write-figure file writers)
 *     need no external LaTeX toolchain and always run
 *  2. Functions that actually invoke pdflatex/pdftoppm
 *     (latex_render_to_png/pdf, latex_compile) only run if
 *     latex_tools_available() reports tools are present : this project
 *     builds and runs its full test suite fine on systems with no LaTeX
 *     distribution installed (latex_gen.c has no LaTeX-time dependency, only a
 *     run-time one), so unconditionally requiring pdflatex here would make test
 *     suite non-portable for no reason
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
  printf("   === Test latex_inline_math / latex_display_math ===\n");

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
  printf("  === Test: latex_matrix, including a cell longer than rough "
         "10-char/cell size estimate ===\n");

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
        "long cell content is not truncated");

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
  printf(
      "  === Test: latex_generate_table writes a well-formed .tex file ===\n");

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
  printf(
      "  === Test latex_generate_article writes a well-formed .tex file ===\n");

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

static void test_escaping(void) {
  printf("  === Test latex_escape_text / latex_write_escaped_text / "
         "latex_validate_label ===\n");

  char *e = latex_escape_text("100% & $5 for #1_thing {x} ~y ^z \\slash");
  check(e != NULL, "latex_escape_text allocates a result");
  check(e && strcmp(e, "100\\% \\& \\$5 for \\#1\\_thing \\{x\\} "
                       "\\textasciitilde{}y \\textasciicircum{}z "
                       "\\textbackslash{}slash") == 0,
        "every special character is escaped correctly");

  free(e);

  char *plain = latex_escape_text("nothing special here");
  check(plain && strcmp(plain, "nothing special here") == 0,
        "text with nothing to escape is returned unchanged");

  free(plain);

  char *empty = latex_escape_text("");
  check(empty != NULL && empty[0] == '\0',
        "empty string escapes to an empty (non-NULL) string");

  free(empty);

  check(latex_escape_text(NULL) == NULL, "NULL input returns NULL cleanly");

  const char *path = "/tmp/qmc_test_escaped_write.txt";
  FILE *f = fopen(path, "w");
  check(f != NULL, "test file opens for latex_write_escaped_text");
  if (f) {
    int rc = latex_write_escaped_text(f, "50% & rising");
    check(rc == 0, "latex_write_escaped_text returns 0 on success");

    fclose(f);

    f = fopen(path, "r");
    char buf[128] = {0};
    if (f) {
      size_t n = fread(buf, 1, sizeof buf - 1, f);
      buf[n] = '\0';

      fclose(f);
    }

    check(strcmp(buf, "50\\% \\& rising") == 0,
          "latex_write_escaped_text writes the escaped form to the file");
  }

  remove(path);

  check(latex_write_escaped_text(NULL, "x") == -1,
        "latex_write_escaped_text(NULL file) returns -1");

  check(latex_validate_label("fig:energy_convergence-01.v2") == 1,
        "a conventional label with ':', '_', '-', '.' validates");
  check(latex_validate_label("") == 0, "empty label is invalid");
  check(latex_validate_label(NULL) == 0, "NULL label is invalid");
  check(latex_validate_label("bad label") == 0,
        "a label containing a space is invalid");
  check(latex_validate_label("bad{label}") == 0,
        "a label containing LaTeX-special characters is invalid");
}

static void test_string_builder(void) {
  printf("  === Test latex_string_t dynamic string builder ===\n");

  latex_string_t s;
  latex_string_init(&s);
  check(s.data == NULL && s.length == 0,
        "freshly initialized builder is empty");

  check(latex_string_append(&s, "Hello, ") == 0, "append returns 0");
  check(latex_string_append(&s, "world") == 0, "second append returns 0");
  check(latex_string_append_char(&s, '!') == 0, "append_char returns 0");
  check(latex_string_appendf(&s, " (%d/%d)", 1, 2) == 0, "appendf returns 0");

  check(strcmp(s.data, "Hello, world! (1/2)") == 0,
        "accumulated content is exactly as appended, in order");
  check(s.length == strlen("Hello, world! (1/2)"),
        "length tracks the accumulated content exactly");

  // Force at least one internal reallocation (initial capacity is small)
  for (int i = 0; i < 200; i++) {
    check(latex_string_append(&s, "x") == 0, "repeated append succeeds");
  }
  check(s.length == strlen("Hello, world! (1/2)") + 200,
        "length is still correct after many small appends (reallocation "
        "didn't corrupt anything)");

  char *detached = latex_string_detach(&s);
  check(detached != NULL, "detach returns the accumulated string");
  check(s.data == NULL && s.length == 0 && s.capacity == 0,
        "detach resets the builder to the freshly-initialized state");

  free(detached);

  // Detaching a builder that never had anything appended still returns a valid
  // (empty) string, not NULL
  latex_string_t empty_s;
  latex_string_init(&empty_s);
  char *empty_detached = latex_string_detach(&empty_s);
  check(empty_detached != NULL && empty_detached[0] == '\0',
        "detaching a never-appended-to builder returns \"\", not NULL");

  free(empty_detached);

  check(latex_string_append(NULL, "x") == -1,
        "append on a NULL builder returns -1, not a crash");
  check(latex_string_append(&s, NULL) == -1,
        "append of a NULL string returns -1, not a crash");

  latex_string_t s2;
  latex_string_init(&s2);
  latex_string_append(&s2, "kept");
  latex_string_free(&s2);
  check(s2.data == NULL && s2.length == 0 && s2.capacity == 0,
        "latex_string_free resets the builder to the empty state");
}

static void test_matrix_type_validation(void) {
  printf(
      "  === Test latex_matrix rejects an unrecognized environment name ===\n");

  const char *row0[1] = {"1"};
  const char *const *data[1] = {row0};

  char *valid = latex_matrix("Vmatrix", data, 1, 1);
  check(valid != NULL, "a valid non-default type ('Vmatrix') still works");

  free(valid);

  check(latex_matrix("array; rm -rf /", data, 1, 1) == NULL,
        "an unrecognized/unsafe environment name is rejected rather than "
        "emitted verbatim into \\begin{...}");
}

static void test_render_options(void) {
  printf("  === Test latex_render_options_t / latex_render_to_png_ex / "
         "latex_render_to_pdf_ex ===\n");

  latex_render_options_t opts = latex_render_options_default();
  check(opts.compiler == NULL && opts.dpi == 200 && opts.use_amsmath == 1 &&
            opts.use_amssymb == 1 && opts.use_physics == 1 &&
            opts.use_bm == 1 && opts.keep_temp_files == 0,
        "latex_render_options_default() matches the documented defaults");

  if (!latex_tools_available()) {
    printf("  (skipped remainder: pdflatex/pdftoppm not available in this "
           "environment)\n");

    return;
  }

  // Passing NULL options must behave exactly like legacy function
  const char *pdf_out = "/tmp/qmc_test_render_ex_default.pdf";
  check(latex_render_to_pdf_ex("x^2", pdf_out, NULL) == 0,
        "latex_render_to_pdf_ex(NULL options) succeeds");
  check(file_exists(pdf_out) && file_size(pdf_out) > 0,
        "latex_render_to_pdf_ex(NULL options) produces a real PDF");

  remove(pdf_out);

  // A custom low DPI actually changes produced PNG's size
  const char *png_lo = "/tmp/qmc_test_render_lodpi.png";
  const char *png_hi = "/tmp/qmc_test_render_hidpi.png";
  latex_render_options_t lo = latex_render_options_default();
  lo.dpi = 50;
  latex_render_options_t hi = latex_render_options_default();
  hi.dpi = 300;

  int rc_lo = latex_render_to_png_ex("x^2+y^2=z^2", png_lo, &lo);
  int rc_hi = latex_render_to_png_ex("x^2+y^2=z^2", png_hi, &hi);
  check(rc_lo == 0 && rc_hi == 0, "both DPI variants render successfully");
  if (rc_lo == 0 && rc_hi == 0) {
    check(file_size(png_hi) > file_size(png_lo),
          "a 300 DPI render is larger than a 50 DPI render of the same "
          "expression (dpi option actually took effect)");
  }

  remove(png_lo);
  remove(png_hi);

  // keep_temp_files=1 leaves job directory behind; latex_last_error()
  // reports its path, and no crash occurs on cleanup-less repeated calls
  const char *pdf_kept = "/tmp/qmc_test_render_keep.pdf";
  latex_render_options_t keep = latex_render_options_default();
  keep.keep_temp_files = 1;
  check(latex_render_to_pdf_ex("a+b", pdf_kept, &keep) == 0,
        "render with keep_temp_files=1 still succeeds");

  const char *msg = latex_last_error();
  check(strstr(msg, "qmc-latex-") != NULL,
        "latex_last_error names the kept temp directory");

  remove(pdf_kept);

  // An explicit but nonexistent compiler is reported as unavailable, not a
  // crash or a misleading error code
  latex_render_options_t bad_compiler = latex_render_options_default();
  bad_compiler.compiler = "definitely-not-a-real-latex-compiler";
  check(latex_render_to_pdf_ex("a", "/tmp/qmc_test_bad_compiler.pdf",
                               &bad_compiler) == -4,
        "an unavailable explicit compiler returns -4");
}

static void test_table_builder(void) {
  printf("  === Test latex_table_t builder (classic and booktabs) ===\n");

  latex_table_t *t = latex_table_create(2, 3);
  check(t != NULL, "latex_table_create succeeds");
  if (!t) {
    return;
  }

  check(latex_table_set_header(t, 0, "Method") == 0, "set_header col 0");
  check(latex_table_set_header(t, 1, "Basis") == 0, "set_header col 1");
  check(latex_table_set_header(t, 2, "E (Ha)") == 0, "set_header col 2");

  check(latex_table_set(t, 0, 0, "RHF", 1) == 0, "set(0,0) escaped text");
  check(latex_table_set(t, 0, 1, "STO-3G", 1) == 0, "set(0,1) escaped text");
  check(latex_table_set(t, 0, 2, "$-1.1167$", 0) == 0,
        "set(0,2) raw LaTeX (escape=0) preserves math delimiters");
  check(latex_table_set(t, 1, 0, "CCSD", 1) == 0, "set(1,0)");
  check(latex_table_set(t, 1, 1, "STO-3G", 1) == 0, "set(1,1)");
  check(latex_table_set(t, 1, 2, "$-1.1373$", 0) == 0, "set(1,2)");

  check(latex_table_set_caption(t, "H$_2$ energies") == 0, "set_caption");
  check(latex_table_set_label(t, "tab:h2") == 0,
        "set_label with a valid label");
  check(latex_table_set_label(t, "bad label") == -1,
        "set_label rejects an invalid label");
  check(latex_table_set_placement(t, "htbp") == 0, "set_placement");

  const char *path = "/tmp/qmc_test_table_classic.tex";
  check(latex_table_write(t, path) == 0, "table_write (classic) succeeds");

  FILE *f = fopen(path, "r");
  char buf[2048] = {0};
  if (f) {
    size_t n = fread(buf, 1, sizeof buf - 1, f);
    buf[n] = '\0';

    fclose(f);
  }

  check(strstr(buf, "\\begin{table}[htbp]") != NULL,
        "classic table: placement appears");
  check(strstr(buf, "\\hline") != NULL, "classic table: uses \\hline rules");
  check(strstr(buf, "Method & Basis") != NULL,
        "classic table: header row present");
  check(strstr(buf, "\\$-1.1167\\$") == NULL &&
            strstr(buf, "$-1.1167$") != NULL,
        "classic table: escape=0 cell is written raw, not escaped");
  check(strstr(buf, "RHF & STO-3G") != NULL, "classic table: data row present");
  check(strstr(buf, "\\caption{H\\$\\_2\\$ energies}") != NULL,
        "classic table: caption is escaped as literal text (caption was given "
        "as \"H$_2$ energies\" - since a caption is documented as literal "
        "text, not raw LaTeX, '$' and '_' get escaped too, same as any other "
        "special character)");
  check(strstr(buf, "\\label{tab:h2}") != NULL,
        "classic table: label present, written verbatim");

  remove(path);

  check(latex_table_set_style(t, LATEX_TABLE_STYLE_BOOKTABS) == 0,
        "set_style(BOOKTABS)");

  const char *path2 = "/tmp/qmc_test_table_booktabs.tex";
  check(latex_table_write(t, path2) == 0, "table_write (booktabs) succeeds");

  FILE *f2 = fopen(path2, "r");
  char buf2[2048] = {0};
  if (f2) {
    size_t n = fread(buf2, 1, sizeof buf2 - 1, f2);
    buf2[n] = '\0';

    fclose(f2);
  }

  check(strstr(buf2, "\\toprule") != NULL &&
            strstr(buf2, "\\midrule") != NULL &&
            strstr(buf2, "\\bottomrule") != NULL,
        "booktabs table: uses \\toprule/\\midrule/\\bottomrule");
  check(strstr(buf2, "\\hline") == NULL,
        "booktabs table: does not also emit \\hline");
  check(strstr(buf2, "\\textbf{Method}") != NULL,
        "booktabs table: header cells are bolded");

  remove(path2);

  latex_table_free(t);

  check(latex_table_create(0, 3) == NULL, "rows=0 returns NULL");
  check(latex_table_create(3, 0) == NULL, "cols=0 returns NULL");

  latex_table_t *bad = latex_table_create(2, 2);
  check(latex_table_set(bad, 5, 0, "x", 1) == -1,
        "set() with an out-of-range row is rejected");
  check(latex_table_set(bad, 0, 5, "x", 1) == -1,
        "set() with an out-of-range col is rejected");
  latex_table_free(bad);

  latex_table_free(NULL);
}

static void test_document_builder(void) {
  printf("  === Test latex_document_t builder ===\n");

  latex_document_t *doc = latex_document_create("Results & Discussion");
  check(doc != NULL, "latex_document_create succeeds");
  if (!doc) {
    return;
  }

  check(latex_document_set_author(doc, "J. Doe") == 0, "set_author");
  check(latex_document_set_date(doc, "2026") == 0, "set_date");
  check(latex_document_set_abstract(doc, "50% faster convergence.") == 0,
        "set_abstract");
  check(latex_document_set_class(doc, "article", "11pt") == 0, "set_class");
  check(latex_document_add_package(doc, "amsmath", NULL) == 0,
        "add_package without options");
  check(latex_document_add_package(doc, "geometry", "margin=1in") == 0,
        "add_package with options");

  check(latex_document_add_section(doc, "Introduction") == 0, "add_section");
  check(latex_document_add_paragraph(doc, "Energy at 50% of baseline.") == 0,
        "add_paragraph (escaped)");
  check(latex_document_add_subsection(doc, "Method") == 0, "add_subsection");
  check(latex_document_add_equation(doc, "\\hat{H}\\psi = E\\psi") == 0,
        "add_equation (raw LaTeX, unescaped)");
  check(latex_document_add_raw(doc, "\\noindent Custom raw LaTeX here.") == 0,
        "add_raw");

  latex_figure_t fig = {
      .path = "energy.pdf",
      .caption = "Energy vs. iteration (n=100%)",
      .label = "fig:energy",
      .width_fraction = 0.6,
      .placement = "htbp",
  };

  check(latex_document_add_figure_ex(doc, &fig) == 0, "add_figure_ex");
  check(latex_document_add_figure(doc, "overlap.pdf", "Overlap matrix",
                                  "fig:overlap") == 0,
        "add_figure (simple 3-arg form)");
  check(latex_document_add_figure_ex(doc, NULL) == -1,
        "add_figure_ex(NULL figure) is rejected");

  latex_table_t *tbl = latex_table_create(1, 2);
  latex_table_set_header(tbl, 0, "Step");
  latex_table_set_header(tbl, 1, "Energy");
  latex_table_set(tbl, 0, 0, "1", 1);
  latex_table_set(tbl, 0, 1, "-2.90", 1);
  check(latex_document_add_table(doc, tbl) == 0,
        "add_table succeeds and takes ownership of table");

  const char *path = "/tmp/qmc_test_document.tex";
  check(latex_document_write(doc, path) == 0, "document_write succeeds");

  FILE *f = fopen(path, "r");
  check(f != NULL, "document output file can be reopened");
  char buf[8192] = {0};
  if (f) {
    size_t n = fread(buf, 1, sizeof buf - 1, f);
    buf[n] = '\0';

    fclose(f);
  }

  check(strstr(buf, "\\documentclass[11pt]{article}") != NULL,
        "document class + options appear");
  check(strstr(buf, "\\usepackage{amsmath}") != NULL,
        "package without options appears");
  check(strstr(buf, "\\usepackage[margin=1in]{geometry}") != NULL,
        "package with options appears");
  check(strstr(buf, "\\title{Results \\& Discussion}") != NULL,
        "title is escaped ('&' becomes '\\&')");
  check(strstr(buf, "\\author{J. Doe}") != NULL, "author appears");
  check(strstr(buf, "\\maketitle") != NULL,
        "maketitle emitted since title/author/date were set");
  check(strstr(buf, "\\begin{abstract}") != NULL &&
            strstr(buf, "50\\% faster") != NULL,
        "abstract appears and is escaped");
  check(strstr(buf, "\\section{Introduction}") != NULL, "section appears");
  check(strstr(buf, "Energy at 50\\% of baseline") != NULL,
        "paragraph text is escaped");
  check(strstr(buf, "\\subsection{Method}") != NULL, "subsection appears");
  check(strstr(buf, "\\[ \\hat{H}\\psi = E\\psi \\]") != NULL,
        "equation is written raw/unescaped (braces and backslashes intact)");
  check(strstr(buf, "\\noindent Custom raw LaTeX here.") != NULL,
        "add_raw content appears verbatim");
  check(strstr(buf, "\\includegraphics[width=0.6\\textwidth]{energy.pdf}") !=
            NULL,
        "figure_ex width_fraction/path appear");
  check(strstr(buf, "n=100\\%") != NULL, "figure caption is escaped");
  check(strstr(buf, "\\label{fig:energy}") != NULL, "figure label appears");
  check(strstr(buf, "\\includegraphics[width=0.8\\textwidth]{overlap.pdf}") !=
            NULL,
        "simple add_figure defaults width_fraction to 0.8");
  check(strstr(buf, "Step & Energy") != NULL,
        "embedded table's header row appears in document body");
  check(strstr(buf, "\\end{document}") != NULL, "document is closed properly");

  remove(path);

  latex_document_free(doc);

  latex_document_t *throwaway = latex_document_create(NULL);
  check(throwaway != NULL,
        "latex_document_create(NULL title) still succeeds (untitled doc)");
  latex_document_free(throwaway);

  latex_document_t *untitled = latex_document_create(NULL);
  check(latex_document_write(untitled, "/tmp/qmc_test_untitled.tex") == 0,
        "an untitled document still writes successfully");

  FILE *uf = fopen("/tmp/qmc_test_untitled.tex", "r");
  char ubuf[512] = {0};
  if (uf) {
    size_t n = fread(ubuf, 1, sizeof ubuf - 1, uf);
    ubuf[n] = '\0';

    fclose(uf);
  }

  check(strstr(ubuf, "\\maketitle") == NULL,
        "no \\maketitle is emitted when no title/author/date was ever set");

  remove("/tmp/qmc_test_untitled.tex");

  latex_document_free(untitled);

  latex_document_free(NULL);
}

static void test_document_compiles_for_real(void) {
  printf("  === Test a document assembled with the new builders actually "
         "compiles with pdflatex ===\n");

  if (!latex_tools_available()) {
    printf(
        "  (skipped: pdflatex/pdftoppm not available in this environment)\n");

    return;
  }

  latex_document_t *doc = latex_document_create("A Real Compile Test");
  latex_document_set_author(doc, "Test Suite");
  latex_document_add_package(doc, "amsmath", NULL);
  latex_document_add_section(doc, "Section One");
  latex_document_add_paragraph(doc, "100% real text with & and # and $.");
  latex_document_add_equation(doc, "E = mc^2");

  latex_table_t *tbl = latex_table_create(2, 2);
  latex_table_set_style(tbl, LATEX_TABLE_STYLE_BOOKTABS);
  latex_table_set_header(tbl, 0, "x");
  latex_table_set_header(tbl, 1, "x^2");
  latex_table_set(tbl, 0, 0, "2", 1);
  latex_table_set(tbl, 0, 1, "$4$", 0);
  latex_table_set(tbl, 1, 0, "3", 1);
  latex_table_set(tbl, 1, 1, "$9$", 0);
  latex_document_add_package(doc, "booktabs", NULL);
  latex_document_add_table(doc, tbl);

  const char *texpath = "/tmp/qmc_test_real_compile.tex";
  check(latex_document_write(doc, texpath) == 0,
        "document with escaped text + booktabs table writes successfully");
  latex_document_free(doc);

  check(latex_compile(texpath, NULL, "/tmp") == 0,
        "pdflatex actually accepts and compiles the generated document "
        "with no errors");

  remove(texpath);
  remove("/tmp/qmc_test_real_compile.pdf");
  remove("/tmp/qmc_test_real_compile.aux");
  remove("/tmp/qmc_test_real_compile.log");
}

static void test_write_figure(void) {
  printf("  === Test Function latex_write_figure writes a well-formed .tex "
         "file, even with NULL optional fields\n");

  const char *path = "/tmp/qmc_test_figure.tex";
  int rc = latex_write_figure(path, "plot.pdf", "A caption", "E = mc^2");

  check(rc == 0, "returns 0 on success");
  check(file_exists(path), "output file created");

  remove(path);

  // NULL optional args should not crash and should still produce valid (if
  // sparse) output rather than dereferencing NULL in fprintf's %s
  rc = latex_write_figure(path, NULL, NULL, NULL);
  check(rc == 0, "NULL optional args handled without crashing");

  remove(path);

  check(latex_write_figure(NULL, "x", "y", "z") == -1,
        "NULL texpath returns an error");
}

static void test_unsafe_paths_rejected(void) {
  printf(
      "  === Test Shell-metacharacter-containing arguments are rejected ===\n");

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

  // Two renders in same process must not clobber each other's intermediate
  // files
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
  printf(" > Testing LaTeX Generation:\n");

  test_inline_display_math();
  test_matrix_generation();
  test_escaping();
  test_string_builder();
  test_matrix_type_validation();
  test_render_options();
  test_table_builder();
  test_document_builder();
  test_document_compiles_for_real();
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
