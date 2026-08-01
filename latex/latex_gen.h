#ifndef QMC_LATEX_GEN_H
#define QMC_LATEX_GEN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Render a LaTeX math expression to PNG.
   outpath: destination PNG file path.
   Returns 0 on success, <0 on error. */
int latex_render_to_png(const char *expr, const char *outpath);

/* Render a LaTeX math expression to PDF.
   outpath: destination PDF file path.
   Returns 0 on success, <0 on error. */
int latex_render_to_pdf(const char *expr, const char *outpath);

/* Write a full LaTeX document with a figure and an equation. */
int latex_write_figure(const char *texpath, const char *figure_pdf,
                       const char *caption, const char *equation);

/* Compile a .tex file to PDF.
   texfile: path to the .tex file.
   compiler: "pdflatex" or "lualatex" (defaults to "pdflatex" if NULL).
   working_dir: directory to run the compiler in (can be NULL for current dir).
   Returns 0 on success, <0 on compilation failure. */
int latex_compile(const char *texfile, const char *compiler,
                  const char *working_dir);

/* Generate a full LaTeX article and write it to a file.
   title, author, date (can be NULL), abstract (can be NULL),
   sections: NULL‑terminated array of char* (each is a LaTeX body).
   Returns 0 on success, <0 on file error. */
int latex_generate_article(const char *texpath, const char *title,
                           const char *author, const char *date,
                           const char *abstract, const char **sections);

/* Generate a LaTeX table from a 2D array of strings.
   data: rows×cols matrix of strings (NULL cell == empty cell).
   col_format: e.g., "|c|c|c|" or "lcr". If NULL, "|c|...c|" is used.
   caption: table caption (can be NULL).
   label: \label value (can be NULL).
   placement: e.g., "htbp" (can be NULL, defaults to "h").
   Returns 0 on success, <0 on file error. */
int latex_generate_table(const char *texpath, const char *const *const *data,
                         int rows, int cols, const char *col_format,
                         const char *caption, const char *label,
                         const char *placement);

/* Wrap an expression in inline math ($...$) or display math (\[...\]).
   The caller must free the returned string. */
char *latex_inline_math(const char *expr);
char *latex_display_math(const char *expr);

/* Create a bmatrix/pmatrix etc. from a 2D array of strings.
   type: "bmatrix", "pmatrix", "Bmatrix", "vmatrix", "Vmatrix".
   Returns a dynamically allocated string, or NULL on error.
   Caller must free. */
char *latex_matrix(const char *type, const char *const *const *data, int rows,
                   int cols);

/* Clean up temporary files created during PNG/PDF rendering.
   This is optional; the system usually removes /tmp files on reboot. */
void latex_clean_temp(void);

#ifdef __cplusplus
}
#endif

#endif /* QMC_LATEX_GEN_H */
