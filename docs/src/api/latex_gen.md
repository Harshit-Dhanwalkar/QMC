# `latex_gen`

> Generated from `latex/latex_gen.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

LaTeX generation: escaping, tables, matrices and documents, rendered to PDF or PNG when pdflatex is available.

**Category:** Export, plotting & LaTeX

```c
#include <qmc/latex/latex_gen.h>
```

**Start with:** [`latex_tools_available`](#latex_tools_available), [`latex_render_to_pdf`](#latex_render_to_pdf), [`latex_document_create`](#latex_document_create)

- Example: `examples/eg_latex_gen.c`
- Test (every call with expected results): `tests/test_latex_gen.c`
- Runtime help: `qmc_help(stdout, "latex_gen");` or `qmc-help latex_gen`

## Types

### `latex_string_t`

```text
Growable string buffer used internally by table/document builders below, and
available directly for callers assembling their own LaTeX fragments
(matrices, custom snippets, etc.) without fixed-estimate
allocate-then-snprintf pattern latex_matrix() uses internally
```

```c
typedef struct {
  char *data;
  size_t length;   /* bytes used, not counting NUL terminator */
  size_t capacity; /* bytes allocated */
} latex_string_t;
```

### `latex_render_options_t`

```text
Rendering options
```

```c
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
```

### `latex_table_style_t`

```text
Table builder
```

```c
typedef enum {
  LATEX_TABLE_STYLE_CLASSIC, /* current default: |c|c|...|, \hline rules */
  LATEX_TABLE_STYLE_BOOKTABS /* \usepackage{booktabs}: ccc..., \toprule/
                              * \midrule/\bottomrule */
} latex_table_style_t;
```

### `latex_table_t`

```c
typedef struct latex_table latex_table_t;
```

### `latex_figure_t`

```text
Figure description (for latex_document_add_figure_ex)
```

```c
typedef struct {
  const char *path;    /* image/PDF path, required */
  const char *caption; /* escaped text; NULL to omit */
  const char *label;   /* validated with latex_validate_label(); NULL to omit */
  double width_fraction; /* fraction of \textwidth, e.g. 0.8; <= 0 means "use
                            default of 0.8" */
  const char *placement; /* e.g. "htbp"; NULL means "h" */
} latex_figure_t;
```

### `latex_document_t`

```text
Document builder
```

```c
typedef struct latex_document latex_document_t;
```

## Functions (55)

[`latex_compiler_available`](#latex_compiler_available), [`latex_png_converter_available`](#latex_png_converter_available), [`latex_tools_available`](#latex_tools_available), [`latex_last_error`](#latex_last_error), [`latex_render_to_png`](#latex_render_to_png), [`latex_render_to_pdf`](#latex_render_to_pdf), [`latex_write_figure`](#latex_write_figure), [`latex_compile`](#latex_compile), [`latex_generate_article`](#latex_generate_article), [`latex_generate_table`](#latex_generate_table), [`latex_inline_math`](#latex_inline_math), [`latex_display_math`](#latex_display_math), [`latex_matrix`](#latex_matrix), [`latex_format_double`](#latex_format_double), [`latex_escape_text`](#latex_escape_text), [`latex_write_escaped_text`](#latex_write_escaped_text), [`latex_validate_label`](#latex_validate_label), [`latex_string_init`](#latex_string_init), [`latex_string_append`](#latex_string_append), [`latex_string_append_char`](#latex_string_append_char), [`latex_string_appendf`](#latex_string_appendf), [`latex_string_detach`](#latex_string_detach), [`latex_string_free`](#latex_string_free), [`latex_render_options_default`](#latex_render_options_default), [`latex_render_to_png_ex`](#latex_render_to_png_ex), [`latex_render_to_pdf_ex`](#latex_render_to_pdf_ex), [`latex_table_create`](#latex_table_create), [`latex_table_free`](#latex_table_free), [`latex_table_set`](#latex_table_set), [`latex_table_set_header`](#latex_table_set_header), [`latex_table_set_caption`](#latex_table_set_caption), [`latex_table_set_label`](#latex_table_set_label), [`latex_table_set_style`](#latex_table_set_style), [`latex_table_set_col_format`](#latex_table_set_col_format), [`latex_table_set_placement`](#latex_table_set_placement), [`latex_table_write`](#latex_table_write), [`latex_document_create`](#latex_document_create), [`latex_document_free`](#latex_document_free), [`latex_document_set_author`](#latex_document_set_author), [`latex_document_set_date`](#latex_document_set_date), [`latex_document_set_abstract`](#latex_document_set_abstract), [`latex_document_set_class`](#latex_document_set_class), [`latex_document_add_package`](#latex_document_add_package), [`latex_document_add_section`](#latex_document_add_section), [`latex_document_add_subsection`](#latex_document_add_subsection), [`latex_document_add_paragraph`](#latex_document_add_paragraph), [`latex_document_add_equation`](#latex_document_add_equation), [`latex_document_add_itemize`](#latex_document_add_itemize), [`latex_document_add_enumerate`](#latex_document_add_enumerate), [`latex_document_add_verbatim`](#latex_document_add_verbatim), [`latex_document_add_figure`](#latex_document_add_figure), [`latex_document_add_figure_ex`](#latex_document_add_figure_ex), [`latex_document_add_table`](#latex_document_add_table), [`latex_document_add_raw`](#latex_document_add_raw), [`latex_document_write`](#latex_document_write)

### `latex_compiler_available`

```c
int latex_compiler_available(void);
```

```text
Returns 1 if configured LaTeX compiler (pdflatex/lualatex) is on PATH
```

### `latex_png_converter_available`

```c
int latex_png_converter_available(void);
```

```text
Returns 1 if pdftoppm (poppler-utils) is on PATH
```

### `latex_tools_available`

```c
int latex_tools_available(void);
```

```text
Both of above
```

### `latex_last_error`

```c
const char *latex_last_error(void);
```

```text
Returns a message describing recent failure, Overwritten on each failing call
NOTE: Empty string if no failure has occurred yet in this process
```

### `latex_render_to_png`

```c
int latex_render_to_png(const char *expr, const char *outpath);
```

```text
Render a LaTeX math expression to PNG
 - outpath: destination PNG file path

Returns 0 on success, <0 on error
  -1 file/argument error, -2 compilation failed, -3 conversion failed,
  -4 LaTeX tools not available, -5 unsafe path/argument rejected
```

### `latex_render_to_pdf`

```c
int latex_render_to_pdf(const char *expr, const char *outpath);
```

```text
Render a LaTeX math expression to PDF
 - outpath: destination PDF file path

Returns 0 on success, <0 on error
```

### `latex_write_figure`

```c
int latex_write_figure(const char *texpath, const char *figure_pdf,
                       const char *caption, const char *equation);
```

```text
Write a full LaTeX document with a figure and an equation
```

### `latex_compile`

```c
int latex_compile(const char *texfile, const char *compiler,
                  const char *working_dir);
```

```text
Compile a .tex file to PDF
texfile: path to .tex file
 - compiler   : "pdflatex" or "lualatex" (defaults to build's configured
                QMC_LATEX_COMPILER, itself defaulting to "pdflatex", if NULL)
 - working_dir: directory to run compiler in (can be NULL for current dir)

Returns 0 on success, <0 on error: -1 argument error, -2 compilation failure,
  -4 LaTeX tools not available, -5 unsafe path/argument rejected
(texfile/compiler/working_dir containing shell metacharacters)
```

### `latex_generate_article`

```c
int latex_generate_article(const char *texpath, const char *title,
                           const char *author, const char *date,
                           const char *abstract, const char **sections);
```

```text
Generate a full LaTeX article and write it to a file
title, author, date (can be NULL), abstract (can be NULL), sections:
NULL-terminated array of char* (each is a LaTeX body)

Returns 0 on success, <0 on file error
```

### `latex_generate_table`

```c
int latex_generate_table(const char *texpath, const char *const *const *data,
                         int rows, int cols, const char *col_format,
                         const char *caption, const char *label,
                         const char *placement);
```

```text
Generate a LaTeX table from a 2D array of strings
 - data      : rows×cols matrix of strings (NULL cell == empty cell)
 - col_format: e.g., "|c|c|c|" or "lcr". If NULL, "|c|...c|" is used
 - caption   : table caption (can be NULL)
 - label     : \label value (can be NULL)
 - placement : e.g., "htbp" (can be NULL, defaults to "h")

Returns 0 on success, <0 on file error
```

### `latex_inline_math`

```c
char *latex_inline_math(const char *expr);
```

```text
Wrap an expression in inline math ($...$) or display math (\[...\])

Returns NULL on allocation failure or if expr is NULL
```

### `latex_display_math`

```c
char *latex_display_math(const char *expr);
```

_No description in header yet._

### `latex_matrix`

```c
char *latex_matrix(const char *type, const char *const *const *data, int rows,
                   int cols);
```

```text
Create a bmatrix/pmatrix etc. from a 2D array of strings
- type: "bmatrix", "pmatrix", "Bmatrix", "vmatrix", "Vmatrix"

Returns a dynamically allocated string, or NULL on error (including rows<=0,
cols<=0, or a NULL data pointer)
```

### `latex_format_double`

```c
char *latex_format_double(double value, int precision);
```

```text
Formats a double to a string suitable for math mode (e.g. "$1.23$")
- value    : double to format
- precision: number of significant digits (e.g. 4)

Returns a dynamically allocated string, or NULL on allocation failure
```

### `latex_escape_text`

```c
char *latex_escape_text(const char *text);
```

```text
Escape LaTeX special characters (& % $ # _ { } ~ ^ \) in literal, user-facing
text - e.g. a title, caption, or paragraph caller does not intend to be
interpreted as latex markup
WARN: Do not use this on an equation or any other string that is already
intentional LaTeX source (e.g. "\\frac{1}{2}mv^2"). API never guesses which
kind of string it was given; caller must pick right function Returns a newly
malloc'd escaped copy (caller frees), or NULL if text is NULL or on
allocation failure
```

### `latex_write_escaped_text`

```c
int latex_write_escaped_text(FILE *f, const char *text);
```

```text
Same as latex_escape_text(), but writes directly to `f` instead of allocating

Returns 0 on success, -1 on a NULL argument or write failure
```

### `latex_validate_label`

```c
int latex_validate_label(const char *label);
```

```text
Validate a \label{...}/\ref{...} identifier: non-empty, and every
character is alphanumeric, ':', '_', '-', or '.' - covers conventional forms
("fig:energy", "tab:results", "eq:hamiltonian") while rejecting characters
that have special meaning to LaTeX or would otherwise be unsafe to place
inside \label{}

Returns 1 if valid, 0 otherwise (including NULL/empty)
```

### `latex_string_init`

```c
void latex_string_init(latex_string_t *s);
```

```text
Zero-initializes `s`. Always succeeds (no allocation happens until first
append)
```

### `latex_string_append`

```c
int latex_string_append(latex_string_t *s, const char *text);
```

```text
Append raw bytes/a NUL-terminated string/a single character/a printf-style
formatted fragment. Every latex_string_append* function returns 0 on success,
-1 on a NULL argument or allocation failure (in which case `s` is left
unchanged - failed append is not partially applied)
```

### `latex_string_append_char`

```c
int latex_string_append_char(latex_string_t *s, char c);
```

_No description in header yet._

### `latex_string_appendf`

```c
int latex_string_appendf(latex_string_t *s, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 2, 3)))
#endif;
```

_No description in header yet._

### `latex_string_detach`

```c
char *latex_string_detach(latex_string_t *s);
```

```text
Transfers ownership of accumulated string to caller (who must free() it) and
resets `s` to freshly-initialized empty state

Returns a NUL-terminated string (possibly "" if nothing was ever appended,
never NULL unless `s` itself is NULL)
```

### `latex_string_free`

```c
void latex_string_free(latex_string_t *s);
```

```text
Frees buffer and resets `s` to freshly-initialized empty state. Safe to call
on an already-empty/never-appended-to builder
```

### `latex_render_options_default`

```c
latex_render_options_t latex_render_options_default(void);
```

```text
Returns default options: compiler=NULL (build default), dpi=200,
amsmath/amssymb/physics/bm all enabled, keep_temp_files=0
```

### `latex_render_to_png_ex`

```c
int latex_render_to_png_ex(const char *expr, const char *outpath,
                           const latex_render_options_t *options);
```

```text
latex_render_to_png()/latex_render_to_pdf() with configurable options
instead of compiler/DPI/package set
latex_render_to_png()/latex_render_to_pdf() remain exactly as before,
implemented as thin wrappers passing latex_render_options_default()
```

### `latex_render_to_pdf_ex`

```c
int latex_render_to_pdf_ex(const char *expr, const char *outpath,
                           const latex_render_options_t *options);
```

_No description in header yet._

### `latex_table_create`

```c
latex_table_t *latex_table_create(size_t rows, size_t cols);
```

```text
Allocates a rows*cols table, every cell initially empty

Returns NULL on invalid dimensions (rows==0 or cols==0) or allocation failure
```

### `latex_table_free`

```c
void latex_table_free(latex_table_t *table);
```

_No description in header yet._

### `latex_table_set`

```c
int latex_table_set(latex_table_t *table, size_t row, size_t col,
                    const char *value, int escape);
```

```text
Set (row, col) cell's content (0-indexed). If `escape` is nonzero,
`value` is treated as literal text and LaTeX-escaped on write; if zero,
`value` is written verbatim as raw LaTeX (e.g. an already-formatted number,
or inline math like "$3.14$")

Returns 0 on success, -1 on an out-of-range row/col or a NULL table/value
```

### `latex_table_set_header`

```c
int latex_table_set_header(latex_table_t *table, size_t col, const char *text);
```

```text
Set a header cell (row 0 of emitted table, rendered bold under
LATEX_TABLE_STYLE_BOOKTABS). Optional - if never called, table has no header
row. `text` is always escaped (headers are conventionally literal column
names, not LaTeX fragments)

Returns 0 on success, -1 on an out-of-range col or a NULL table/text
```

### `latex_table_set_caption`

```c
int latex_table_set_caption(latex_table_t *table, const char *caption);
```

```text
All optional; `caption`/`label`/`col_format` may be passed NULL to omit/reset
to default. `label`, if non-NULL, is validated with latex_validate_label()
and rejected (-1) if invalid. `col_format` overrides style-driven default
(e.g. "lcr" or "|l|c|r|"). Every setter returns 0 on success, -1 on a NULL
table or (for label setter) an invalid label
```

### `latex_table_set_label`

```c
int latex_table_set_label(latex_table_t *table, const char *label);
```

_No description in header yet._

### `latex_table_set_style`

```c
int latex_table_set_style(latex_table_t *table, latex_table_style_t style);
```

_No description in header yet._

### `latex_table_set_col_format`

```c
int latex_table_set_col_format(latex_table_t *table, const char *col_format);
```

_No description in header yet._

### `latex_table_set_placement`

```c
int latex_table_set_placement(latex_table_t *table, const char *placement);
```

_No description in header yet._

### `latex_table_write`

```c
int latex_table_write(latex_table_t *table, const char *texpath);
```

```text
Writes table (as a complete \begin{table}...\end{table} block) to
`texpath`

Returns 0 on success, -1 on a NULL argument or file error
```

### `latex_document_create`

```c
latex_document_t *latex_document_create(const char *title);
```

```text
Creates a document with document class "article" (no size option) and no
packages
`title` may be NULL for an untitled document

Returns NULL on allocation failure
```

### `latex_document_free`

```c
void latex_document_free(latex_document_t *doc);
```

_No description in header yet._

### `latex_document_set_author`

```c
int latex_document_set_author(latex_document_t *doc, const char *author);
```

```text
All escaped as literal text except set_class's `document_class`/
`options` and add_package's `package`/`options`, which are identifiers
written verbatim (e.g. "article", "11pt", "amsmath", "margin=1in").
Every setter here returns 0 on success, -1 on a NULL doc
```

### `latex_document_set_date`

```c
int latex_document_set_date(latex_document_t *doc, const char *date);
```

_No description in header yet._

### `latex_document_set_abstract`

```c
int latex_document_set_abstract(latex_document_t *doc, const char *abstract);
```

_No description in header yet._

### `latex_document_set_class`

```c
int latex_document_set_class(latex_document_t *doc, const char *document_class,
                             const char *options);
```

```text
`options` is documentclass's bracketed option list, e.g. "11pt" or
"12pt,twocolumn"; NULL omits brackets entirely
Default (if never called) is document_class="article" with no options
```

### `latex_document_add_package`

```c
int latex_document_add_package(latex_document_t *doc, const char *package,
                               const char *options);
```

```text
Adds \usepackage{package}, or \usepackage[options]{package} if options
is non-NULL. A package added more than once is written more than once
(LaTeX tolerates this; builder does not deduplicate)
```

### `latex_document_add_section`

```c
int latex_document_add_section(latex_document_t *doc, const char *title);
```

```text
Body content, appended in order called

Returns 0 on success, -1 on a NULL doc/argument, allocation failure, or
(add_figure/add_table) an invalid label
```

### `latex_document_add_subsection`

```c
int latex_document_add_subsection(latex_document_t *doc, const char *title);
```

_No description in header yet._

### `latex_document_add_paragraph`

```c
int latex_document_add_paragraph(latex_document_t *doc, const char *text);
```

_No description in header yet._

### `latex_document_add_equation`

```c
int latex_document_add_equation(latex_document_t *doc, const char *equation);
```

_No description in header yet._

### `latex_document_add_itemize`

```c
int latex_document_add_itemize(latex_document_t *doc, const char *const *items);
```

_No description in header yet._

### `latex_document_add_enumerate`

```c
int latex_document_add_enumerate(latex_document_t *doc,
                                 const char *const *items);
```

_No description in header yet._

### `latex_document_add_verbatim`

```c
int latex_document_add_verbatim(latex_document_t *doc, const char *text);
```

_No description in header yet._

### `latex_document_add_figure`

```c
int latex_document_add_figure(latex_document_t *doc, const char *path,
                              const char *caption, const char *label);
```

_No description in header yet._

### `latex_document_add_figure_ex`

```c
int latex_document_add_figure_ex(latex_document_t *doc,
                                 const latex_figure_t *figure);
```

_No description in header yet._

### `latex_document_add_table`

```c
int latex_document_add_table(latex_document_t *doc, latex_table_t *table);
```

```text
Appends table's current content and frees `table` (ownership transfers in on
success; on failure -1 is returned and `table` is left for caller to free)
```

### `latex_document_add_raw`

```c
int latex_document_add_raw(latex_document_t *doc, const char *latex);
```

```text
Escape hatch: appends `latex` verbatim, exactly as given, with no  escaping
and no interpretation
```

### `latex_document_write`

```c
int latex_document_write(latex_document_t *doc, const char *texpath);
```

```text
Assembles and writes complete .tex file (documentclass, packages,
title/author/date, \begin{document}, \maketitle if a title/author/date was
set, abstract if set, body content in order added, \end{document})

Returns 0 on success, -1 on a NULL argument or file error
```
