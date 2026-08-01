#include "latex_gen.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// NOTE: Requires pdflatex and pdftoppm (poppler-utils) installed.

/* Helpers */
// Path join
static void make_tmp_path(const char *basename, char *buf, size_t bufsz) {
  snprintf(buf, bufsz, "/tmp/%s", basename);
}

/* Rendering functions */
int latex_render_to_png(const char *expr, const char *outpath) {
  // Write .tex file
  FILE *f = fopen("/tmp/qmc_eq.tex", "w");
  if (!f) {
    return -1;
  }
  fprintf(f,
          "\\documentclass[12pt,preview]{standalone}\n"
          "\\usepackage{amsmath,amssymb,physics,bm}\n"
          "\\begin{document}\n"
          "\\[ %s \\]\n"
          "\\end{document}\n",
          expr);
  fclose(f);

  // Compile
  int r = system("cd /tmp && pdflatex -interaction=batchmode qmc_eq.tex "
                 "> /dev/null 2>&1");
  if (r != 0) {
    return -2;
  }

  // Convert to PNG
  char cmd[1024];
  snprintf(cmd, sizeof(cmd),
           "pdftoppm -r 200 -png /tmp/qmc_eq.pdf /tmp/qmc_eq_out "
           "&& mv /tmp/qmc_eq_out-1.png %s",
           outpath);
  r = system(cmd);

  return (r == 0) ? 0 : -3;
}

int latex_render_to_pdf(const char *expr, const char *outpath) {
  FILE *f = fopen("/tmp/qmc_eq.tex", "w");
  if (!f) {
    return -1;
  }
  fprintf(f,
          "\\documentclass[12pt,preview]{standalone}\n"
          "\\usepackage{amsmath,amssymb,physics,bm}\n"
          "\\begin{document}\n"
          "\\[ %s \\]\n"
          "\\end{document}\n",
          expr);
  fclose(f);

  int r = system("cd /tmp && pdflatex -interaction=batchmode qmc_eq.tex "
                 "> /dev/null 2>&1");
  if (r != 0) {
    return -2;
  }

  char cmd[512];
  snprintf(cmd, sizeof(cmd), "mv /tmp/qmc_eq.pdf %s", outpath);
  r = system(cmd);

  return (r == 0) ? 0 : -3;
}

int latex_write_figure(const char *texpath, const char *figure_pdf,
                       const char *caption, const char *equation) {
  FILE *f = fopen(texpath, "w");
  if (!f) {
    return -1;
  }
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

/* Compilation helper  */
int latex_compile(const char *texfile, const char *compiler,
                  const char *working_dir) {
  if (!texfile) {
    return -1;
  }
  if (!compiler) {
    compiler = "pdflatex"; // default
                           // compiler = "lualatex";
  }

  char cmd[1024];
  if (working_dir && working_dir[0]) {
    snprintf(cmd, sizeof(cmd),
             "cd %s && %s -interaction=batchmode %s "
             "> /dev/null 2>&1",
             working_dir, compiler, texfile);
  } else {
    snprintf(cmd, sizeof(cmd),
             "%s -interaction=batchmode %s "
             "> /dev/null 2>&1",
             compiler, texfile);
  }
  return system(cmd) == 0 ? 0 : -2;
}

/* Document generation  */
int latex_generate_article(const char *texpath, const char *title,
                           const char *author, const char *date,
                           const char *abstract, const char **sections) {
  FILE *f = fopen(texpath, "w");
  if (!f) {
    return -1;
  }

  fprintf(f, "\\documentclass{article}\n");
  fprintf(f, "\\usepackage{amsmath,amssymb,physics,graphicx,bm}\n");
  if (title) {
    fprintf(f, "\\title{%s}\n", title);
  }
  if (author) {
    fprintf(f, "\\author{%s}\n", author);
  }
  if (date) {
    fprintf(f, "\\date{%s}\n", date);
  }
  fprintf(f, "\\begin{document}\n");
  fprintf(f, "\\maketitle\n");

  if (abstract && abstract[0]) {
    fprintf(f, "\\begin{abstract}\n%s\n\\end{abstract}\n", abstract);
  }

  if (sections) {
    for (int i = 0; sections[i] != NULL; i++) {
      fprintf(f, "%s\n\n", sections[i]);
    }
  }

  fprintf(f, "\\end{document}\n");
  fclose(f);

  return 0;
}

/* Table generation */
int latex_generate_table(const char *texpath, const char *const *const *data,
                         int rows, int cols, const char *col_format,
                         const char *caption, const char *label,
                         const char *placement) {
  FILE *f = fopen(texpath, "w");
  if (!f) {
    return -1;
  }

  // If no format given, create "|c|c|...|"
  char format[256] = {0};
  if (!col_format) {
    char *p = format;
    *p++ = '|';
    for (int i = 0; i < cols; i++) {
      *p++ = 'c';
      *p++ = '|';
    }

    *p = '\0';
  }

  fprintf(f, "\\begin{table}[%s]\n\\centering\n", placement ? placement : "h");
  fprintf(f, "\\begin{tabular}{%s}\n\\hline\n",
          col_format ? col_format : format);

  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      fprintf(f, "%s", data[r][c] ? data[r][c] : "");
      if (c < cols - 1) {
        fprintf(f, " & ");
      }
    }

    fprintf(f, " \\\\ \\hline\n");
  }

  fprintf(f, "\\end{tabular}\n");
  if (caption) {
    fprintf(f, "\\caption{%s}\n", caption);
  }

  if (label) {
    fprintf(f, "\\label{%s}\n", label);
  }
  fprintf(f, "\\end{table}\n");

  fclose(f);
  return 0;
}

/* Math helpers */

char *latex_inline_math(const char *expr) {
  size_t len = strlen(expr) + 4; /* $...$ */
  char *res = malloc(len);
  if (!res) {
    return NULL;
  }

  snprintf(res, len, "$%s$", expr);

  return res;
}

char *latex_display_math(const char *expr) {
  size_t len = strlen(expr) + 8; /* \[ ... \] */
  char *res = malloc(len);
  if (!res) {
    return NULL;
  }

  snprintf(res, len, "\\[ %s \\]", expr);

  return res;
}

char *latex_matrix(const char *type, const char *const *const *data, int rows,
                   int cols) {
  if (!type) {
    type = "bmatrix";
  }
  // HACK: Rough estimate: 10 chars per cell + braces
  size_t est = 20 + strlen(type) + 2 + rows * cols * 20;
  char *res = malloc(est);
  if (!res) {
    return NULL;
  }
  size_t pos = 0;

  pos += snprintf(res + pos, est - pos, "\\begin{%s}", type);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      pos +=
          snprintf(res + pos, est - pos, "%s", data[i][j] ? data[i][j] : "0");
      if (j < cols - 1) {
        pos += snprintf(res + pos, est - pos, " & ");
      }
    }

    if (i < rows - 1) {
      pos += snprintf(res + pos, est - pos, " \\\\ ");
    }
  }

  pos += snprintf(res + pos, est - pos, "\\end{%s}", type);

  return res;
}

void latex_clean_temp(void) {
  // Remove temporary files created
  remove("/tmp/qmc_eq.tex");
  remove("/tmp/qmc_eq.aux");
  remove("/tmp/qmc_eq.log");
  remove("/tmp/qmc_eq.pdf");
  remove("/tmp/qmc_eq_out-1.png");
  // WARN: there may be more temp files from multiple runs
}
