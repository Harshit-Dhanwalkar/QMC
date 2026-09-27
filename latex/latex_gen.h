#ifndef QMC_LATEX_GEN_H
#define QMC_LATEX_GEN_H

#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Returns 1 if the configured LaTeX compiler (pdflatex/lualatex) is on PATH */
int latex_compiler_available(void);

/* Returns 1 if pdftoppm (poppler-utils) is on PATH */
int latex_png_converter_available(void);

/* Both of above */
int latex_tools_available(void);

/* Returns a message describing recent failure, Overwritten on each failing call
 * NOTE: Empty string if no failure has occurred yet in this process */
const char *latex_last_error(void);

/* Render a LaTeX math expression to PNG
 *  - outpath: destination PNG file path
 *
 * Returns 0 on success, <0 on error
 *   -1 file/argument error, -2 compilation failed, -3 conversion failed,
 *   -4 LaTeX tools not available, -5 unsafe path/argument rejected
 */
int latex_render_to_png(const char *expr, const char *outpath);

/* Render a LaTeX math expression to PDF
 *  - outpath: destination PDF file path
 *
 * Returns 0 on success, <0 on error
 */
int latex_render_to_pdf(const char *expr, const char *outpath);

/* Write a full LaTeX document with a figure and an equation */
int latex_write_figure(const char *texpath, const char *figure_pdf,
                       const char *caption, const char *equation);

/* Compile a .tex file to PDF
 * texfile: path to the .tex file
 *  - compiler   : "pdflatex" or "lualatex" (defaults to build's configured
 *                 QMC_LATEX_COMPILER, itself defaulting to "pdflatex", if NULL)
 *  - working_dir: directory to run compiler in (can be NULL for current dir)
 *
 * Returns 0 on success, <0 on error: -1 argument error, -2 compilation failure,
 *   -4 LaTeX tools not available, -5 unsafe path/argument rejected
 * (texfile/compiler/working_dir containing shell metacharacters)
 */
int latex_compile(const char *texfile, const char *compiler,
                  const char *working_dir);

/* Generate a full LaTeX article and write it to a file
 * title, author, date (can be NULL), abstract (can be NULL), sections:
 * NULL‑terminated array of char* (each is a LaTeX body)
 *
 * Returns 0 on success, <0 on file error
 */
int latex_generate_article(const char *texpath, const char *title,
                           const char *author, const char *date,
                           const char *abstract, const char **sections);

/* Generate a LaTeX table from a 2D array of strings
 *  - data      : rows×cols matrix of strings (NULL cell == empty cell)
 *  - col_format: e.g., "|c|c|c|" or "lcr". If NULL, "|c|...c|" is used
 *  - caption   : table caption (can be NULL)
 *  - label     : \label value (can be NULL)
 *  - placement : e.g., "htbp" (can be NULL, defaults to "h")
 *
 *  Returns 0 on success, <0 on file error
 */
int latex_generate_table(const char *texpath, const char *const *const *data,
                         int rows, int cols, const char *col_format,
                         const char *caption, const char *label,
                         const char *placement);

/* Wrap an expression in inline math ($...$) or display math (\[...\])
 *
 * Returns NULL on allocation failure or if expr is NULL
 */
char *latex_inline_math(const char *expr);
char *latex_display_math(const char *expr);

/* Create a bmatrix/pmatrix etc. from a 2D array of strings
 * - type: "bmatrix", "pmatrix", "Bmatrix", "vmatrix", "Vmatrix"
 *
 * Returns a dynamically allocated string, or NULL on error (including rows<=0,
 * cols<=0, or a NULL data pointer)
 */
char *latex_matrix(const char *type, const char *const *const *data, int rows,
                   int cols);

/* Clean up temporary files created during PNG/PDF rendering */
void latex_clean_temp(void);

/* Text escaping */
/* Escape LaTeX special characters (& % $ # _ { } ~ ^ \) in literal, user-facing
 * text - e.g. a title, caption, or paragraph caller does not intend to be
 * interpreted as latex markup
 * WARN: Do not use this on an equation or any other string that is already
 * intentional LaTeX source (e.g. "\\frac{1}{2}mv^2"). API never guesses which
 * kind of string it was given; caller must pick right function Returns a newly
 * malloc'd escaped copy (caller frees), or NULL if text is NULL or on
 * allocation failure
 */
char *latex_escape_text(const char *text);

/* Same as latex_escape_text(), but writes directly to `f` instead of allocating
 *
 * Returns 0 on success, -1 on a NULL argument or write failure
 */
int latex_write_escaped_text(FILE *f, const char *text);

/* Validate a \label{...}/\ref{...} identifier: non-empty, and every
 * character is alphanumeric, ':', '_', '-', or '.' - covers conventional forms
 * ("fig:energy", "tab:results", "eq:hamiltonian") while rejecting characters
 * that have special meaning to LaTeX or would otherwise be unsafe to place
 * inside \label{}
 *
 * Returns 1 if valid, 0 otherwise (including NULL/empty)
 */
int latex_validate_label(const char *label);

/* Dynamic string builder */
/* Growable string buffer used internally by table/document builders below, and
 * available directly for callers assembling their own LaTeX fragments
 * (matrices, custom snippets, etc.) without fixed-estimate
 * allocate-then-snprintf pattern latex_matrix() uses internally
 */
typedef struct {
  char *data;
  size_t length;   /* bytes used, not counting NUL terminator */
  size_t capacity; /* bytes allocated */
} latex_string_t;

/* Zero-initializes `s`. Always succeeds (no allocation happens until first
 * append)
 */
void latex_string_init(latex_string_t *s);

/* Append raw bytes/a NUL-terminated string/a single character/a printf-style
 * formatted fragment. Every latex_string_append* function returns 0 on success,
 * -1 on a NULL argument or allocation failure (in which case `s` is left
 * unchanged - failed append is not partially applied)
 */
int latex_string_append(latex_string_t *s, const char *text);
int latex_string_append_char(latex_string_t *s, char c);
int latex_string_appendf(latex_string_t *s, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 2, 3)))
#endif
    ;

/* Transfers ownership of accumulated string to caller (who must free() it) and
 * resets `s` to freshly-initialized empty state
 *
 * Returns a NUL-terminated string (possibly "" if nothing was ever appended,
 * never NULL unless `s` itself is NULL)
 */
char *latex_string_detach(latex_string_t *s);

/* Frees buffer and resets `s` to freshly-initialized empty state. Safe to call
 * on an already-empty/never-appended-to builder
 */
void latex_string_free(latex_string_t *s);

/* Rendering options */
typedef struct {
  const char *compiler; /* "pdflatex"/"lualatex"/etc; NULL -> build's configured
                           QMC_LATEX_COMPILER default */
  int dpi; /* PNG rendering resolution; ignored by PDF path. Default 200 */
  int use_amsmath;
  int use_amssymb;
  int use_physics;
  int use_bm;
  int keep_temp_files; /* 1 = don't delete render job's temp directory - useful
                          for diagnosing a compile failure by hand */
} latex_render_options_t;

/* Returns default options: compiler=NULL (build default), dpi=200,
 * amsmath/amssymb/physics/bm all enabled, keep_temp_files=0
 */
latex_render_options_t latex_render_options_default(void);

/* latex_render_to_png()/latex_render_to_pdf() with configurable options
 * instead of compiler/DPI/package set
 * latex_render_to_png()/latex_render_to_pdf() remain exactly as before,
 * implemented as thin wrappers passing latex_render_options_default()
 */
int latex_render_to_png_ex(const char *expr, const char *outpath,
                           const latex_render_options_t *options);
int latex_render_to_pdf_ex(const char *expr, const char *outpath,
                           const latex_render_options_t *options);

/* Table builder */
typedef enum {
  LATEX_TABLE_STYLE_CLASSIC, /* current default: |c|c|...|, \hline rules */
  LATEX_TABLE_STYLE_BOOKTABS /* \usepackage{booktabs}: ccc..., \toprule/
                              * \midrule/\bottomrule */
} latex_table_style_t;

typedef struct latex_table latex_table_t;

/* Allocates a rows*cols table, every cell initially empty
 *
 * Returns NULL on invalid dimensions (rows==0 or cols==0) or allocation failure
 */
latex_table_t *latex_table_create(size_t rows, size_t cols);
void latex_table_free(latex_table_t *table);

/* Set the (row, col) cell's content (0-indexed). If `escape` is nonzero,
 * `value` is treated as literal text and LaTeX-escaped on write; if zero,
 * `value` is written verbatim as raw LaTeX (e.g. an already-formatted number,
 * or inline math like "$3.14$")
 *
 * Returns 0 on success, -1 on an out-of-range row/col or a NULL table/value
 */
int latex_table_set(latex_table_t *table, size_t row, size_t col,
                    const char *value, int escape);

/* Set a header cell (row 0 of the emitted table, rendered bold under
 * LATEX_TABLE_STYLE_BOOKTABS). Optional - if never called, the table has no
 * header row. `text` is always escaped (headers are conventionally literal
 * column names, not LaTeX fragments)
 *
 * Returns 0 on success, -1 on an out-of-range col or a NULL table/text
 */
int latex_table_set_header(latex_table_t *table, size_t col, const char *text);

/* All optional; `caption`/`label`/`col_format` may be passed NULL to omit/reset
 * to default. `label`, if non-NULL, is validated with latex_validate_label()
 * and rejected (-1) if invalid. `col_format` overrides the style-driven default
 * (e.g. "lcr" or "|l|c|r|"). Every setter returns 0 on success, -1 on a NULL
 * table or (for the label setter) an invalid label
 */
int latex_table_set_caption(latex_table_t *table, const char *caption);
int latex_table_set_label(latex_table_t *table, const char *label);
int latex_table_set_style(latex_table_t *table, latex_table_style_t style);
int latex_table_set_col_format(latex_table_t *table, const char *col_format);
int latex_table_set_placement(latex_table_t *table, const char *placement);

/* Writes the table (as a complete \begin{table}...\end{table} block) to
 * `texpath`. Returns 0 on success, -1 on a NULL argument or file error */
int latex_table_write(latex_table_t *table, const char *texpath);

/* Figure description (for latex_document_add_figure_ex) */
typedef struct {
  const char *path;    /* image/PDF path, required */
  const char *caption; /* escaped text; NULL to omit */
  const char *label;   /* validated with latex_validate_label(); NULL to omit */
  double width_fraction; /* fraction of \textwidth, e.g. 0.8; <= 0 means "use
                            default of 0.8" */
  const char *placement; /* e.g. "htbp"; NULL means "h" */
} latex_figure_t;

/* Document builder */
typedef struct latex_document latex_document_t;

/* Creates a document with document class "article" (no size option) and no
 * packages
 * `title` may be NULL for an untitled document
 *
 * Returns NULL on allocation failure
 */
latex_document_t *latex_document_create(const char *title);
void latex_document_free(latex_document_t *doc);

/* All escaped as literal text except set_class's `document_class`/
 * `options` and add_package's `package`/`options`, which are identifiers
 * written verbatim (e.g. "article", "11pt", "amsmath", "margin=1in").
 * Every setter here returns 0 on success, -1 on a NULL doc
 */
int latex_document_set_author(latex_document_t *doc, const char *author);
int latex_document_set_date(latex_document_t *doc, const char *date);
int latex_document_set_abstract(latex_document_t *doc, const char *abstract);

/* `options` is the documentclass's bracketed option list, e.g. "11pt" or
 * "12pt,twocolumn"; NULL omits the brackets entirely
 * Default (if never called) is document_class="article" with no options
 */
int latex_document_set_class(latex_document_t *doc, const char *document_class,
                             const char *options);

/* Adds \usepackage{package}, or \usepackage[options]{package} if options
 * is non-NULL. A package added more than once is written more than once
 * (LaTeX tolerates this; the builder does not deduplicate)
 */
int latex_document_add_package(latex_document_t *doc, const char *package,
                               const char *options);

/* Body content, appended in the order called. Returns 0 on success, -1
 * on a NULL doc/argument, allocation failure, or (add_figure/add_table)
 * an invalid label
 */
int latex_document_add_section(latex_document_t *doc, const char *title);
int latex_document_add_subsection(latex_document_t *doc, const char *title);
int latex_document_add_paragraph(latex_document_t *doc, const char *text);
int latex_document_add_equation(latex_document_t *doc, const char *equation);
int latex_document_add_figure(latex_document_t *doc, const char *path,
                              const char *caption, const char *label);
int latex_document_add_figure_ex(latex_document_t *doc,
                                 const latex_figure_t *figure);
/* Appends the table's current content and frees `table` (ownership transfers in
 * on success; on failure -1 is returned and `table` is left for the caller to
 * free) */
int latex_document_add_table(latex_document_t *doc, latex_table_t *table);
/* Escape hatch: appends `latex` verbatim, exactly as given, with no  escaping
 * and no interpretation */
int latex_document_add_raw(latex_document_t *doc, const char *latex);

/* Assembles and writes the complete .tex file (documentclass, packages,
 * title/author/date, \begin{document}, \maketitle if a title/author/date was
 * set, abstract if set, body content in the order added, \end{document})
 *
 * Returns 0 on success, -1 on a NULL argument or file error 
 */
int latex_document_write(latex_document_t *doc, const char *texpath);

#ifdef __cplusplus
}
#endif

#endif /* QMC_LATEX_GEN_H */
