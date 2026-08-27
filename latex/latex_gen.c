#include "latex_gen.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

// NOTE: Requires pdflatex (or lualatex) and pdftoppm (poppler-utils)
//       installed and on PATH.

#ifndef QMC_LATEX_COMPILER
#define QMC_LATEX_COMPILER "pdflatex"
#endif

#define QMC_LATEX_MAX_TRACKED_TEMP 256

static char g_tracked_temp_basenames[QMC_LATEX_MAX_TRACKED_TEMP][64];
static int g_tracked_temp_count = 0;

/* Helpers */
// Path join
static void make_tmp_path(const char *basename, const char *ext, char *buf,
                          size_t bufsz) {
  snprintf(buf, bufsz, "/tmp/%s.%s", basename, ext);
}

static void make_unique_basename(char *buf, size_t bufsz) {
  static long counter = 0;
  long id = (long)getpid() * 1000003L + counter++;
  snprintf(buf, bufsz, "qmc_eq_%ld", id);

  if (g_tracked_temp_count < QMC_LATEX_MAX_TRACKED_TEMP) {
    snprintf(g_tracked_temp_basenames[g_tracked_temp_count++],
             sizeof(g_tracked_temp_basenames[0]), "%s", buf);
  }
}

// Runs a shell command and reports success/failure correctly
static int run_system(const char *cmd) {
  int status = system(cmd);
  if (status == -1) {
    return -1;
  }

  if (WIFEXITED(status)) {
    return WEXITSTATUS(status) == 0 ? 0 : -1;
  }

  return -1; // terminated by signal, etc.
}

static int has_unsafe_shell_chars(const char *s) {
  if (!s) {
    return 0;
  }

  return strpbrk(s, ";&|`$()<>\"'\\\n") != NULL;
}

static int check_tool_on_path(const char *tool) {
  char cmd[256];
  snprintf(cmd, sizeof(cmd), "command -v %s >/dev/null 2>&1", tool);

  return run_system(cmd) == 0;
}

int latex_tools_available(void) {
  return check_tool_on_path(QMC_LATEX_COMPILER) &&
         check_tool_on_path("pdftoppm");
}

// Shared .tex-writing logic for latex_render_to_png/pdf
static int write_standalone_equation_tex(const char *texpath,
                                         const char *expr) {
  FILE *f = fopen(texpath, "w");
  if (!f) {
    return -1;
  }

  int rc = fprintf(f,
                   "\\documentclass[12pt,preview]{standalone}\n"
                   "\\usepackage{amsmath,amssymb,physics,bm}\n"
                   "\\begin{document}\n"
                   "\\[ %s \\]\n"
                   "\\end{document}\n",
                   expr);

  int close_rc = fclose(f);

  return (rc >= 0 && close_rc == 0) ? 0 : -1;
}

/* Rendering functions */
int latex_render_to_png(const char *expr, const char *outpath) {
  if (!expr || !outpath) {
    return -1;
  }

  if (has_unsafe_shell_chars(outpath)) {
    return -5;
  }

  if (!latex_tools_available()) {
    return -4;
  }

  char basename[64];
  make_unique_basename(basename, sizeof(basename));

  char texpath[128], pdfpath[128], pngprefix[128], pngout[136];
  make_tmp_path(basename, "tex", texpath, sizeof(texpath));
  make_tmp_path(basename, "pdf", pdfpath, sizeof(pdfpath));

  snprintf(pngprefix, sizeof(pngprefix), "/tmp/%s_out", basename);
  snprintf(pngout, sizeof(pngout), "%s-1.png", pngprefix);

  if (write_standalone_equation_tex(texpath, expr) != 0) {
    return -1;
  }

  char cmd[512];
  snprintf(cmd, sizeof(cmd),
           "cd /tmp && %s -interaction=batchmode %s > /dev/null 2>&1",
           QMC_LATEX_COMPILER, texpath);
  if (run_system(cmd) != 0) {
    return -2;
  }

  snprintf(cmd, sizeof(cmd), "pdftoppm -r 200 -png %s %s", pdfpath, pngprefix);
  if (run_system(cmd) != 0) {
    return -3;
  }

  if (rename(pngout, outpath) != 0) {
    return -3;
  }

  return 0;
}

int latex_render_to_pdf(const char *expr, const char *outpath) {
  FILE *f = fopen("/tmp/qmc_eq.tex", "w");
  if (!expr || !outpath) {
    return -1;
  }

  if (has_unsafe_shell_chars(outpath)) {
    return -5;
  }

  if (!latex_tools_available()) {
    return -4;
  }

  char basename[64];
  make_unique_basename(basename, sizeof(basename));

  char texpath[128], pdfpath[128];
  make_tmp_path(basename, "tex", texpath, sizeof(texpath));
  make_tmp_path(basename, "pdf", pdfpath, sizeof(pdfpath));

  if (write_standalone_equation_tex(texpath, expr) != 0) {
    return -1;
  }

  char cmd[512];
  snprintf(cmd, sizeof(cmd),
           "cd /tmp && %s -interaction=batchmode %s > /dev/null 2>&1",
           QMC_LATEX_COMPILER, texpath);
  if (run_system(cmd) != 0) {
    return -2;
  }

  if (rename(pdfpath, outpath) != 0) {
    return -3;
  }

  return 0;
}

int latex_write_figure(const char *texpath, const char *figure_pdf,
                       const char *caption, const char *equation) {
  if (!texpath) {
    return -1;
  }

  FILE *f = fopen(texpath, "w");
  if (!f) {
    return -1;
  }

  int rc = fprintf(f,
                   "\\documentclass{article}\n"
                   "\\usepackage{amsmath,amssymb,physics,graphicx,bm}\n"
                   "\\begin{document}\n"
                   "\\begin{figure}[h]\\centering\n"
                   "  \\includegraphics[width=0.8\\textwidth]{%s}\n"
                   "  \\caption{%s}\n"
                   "\\end{figure}\n"
                   "\\[ %s \\]\n"
                   "\\end{document}\n",
                   figure_pdf ? figure_pdf : "", caption ? caption : "",
                   equation ? equation : "");
  int close_rc = fclose(f);

  return (rc >= 0 && close_rc == 0) ? 0 : -1;
}

/* Compilation helper  */
int latex_compile(const char *texfile, const char *compiler,
                  const char *working_dir) {
  if (!texfile) {
    return -1;
  }

  if (!compiler) {
    compiler = QMC_LATEX_COMPILER;
  }

  if (has_unsafe_shell_chars(texfile) || has_unsafe_shell_chars(compiler) ||
      has_unsafe_shell_chars(working_dir)) {
    return -5;
  }

  if (!check_tool_on_path(compiler)) {
    return -4;
  }

  char cmd[1024];
  int written;
  if (working_dir && working_dir[0]) {
    written = snprintf(cmd, sizeof(cmd),
                       "cd %s && %s -interaction=batchmode %s "
                       "> /dev/null 2>&1",
                       working_dir, compiler, texfile);
  } else {
    written = snprintf(cmd, sizeof(cmd),
                       "%s -interaction=batchmode %s "
                       "> /dev/null 2>&1",
                       compiler, texfile);
  }

  if (written < 0 || (size_t)written >= sizeof(cmd)) {
    return -1; // command too long for the buffer
  }

  return run_system(cmd) == 0 ? 0 : -2;
}

/* Document generation  */
int latex_generate_article(const char *texpath, const char *title,
                           const char *author, const char *date,
                           const char *abstract, const char **sections) {
  if (!texpath) {
    return -1;
  }

  FILE *f = fopen(texpath, "w");
  if (!f) {
    return -1;
  }

  int ok = 1;
  ok &= fprintf(f, "\\documentclass{article}\n") >= 0;
  ok &= fprintf(f, "\\usepackage{amsmath,amssymb,physics,graphicx,bm}\n") >= 0;

  if (title) {
    ok &= fprintf(f, "\\title{%s}\n", title) >= 0;
  }
  if (author) {
    ok &= fprintf(f, "\\author{%s}\n", author) >= 0;
  }
  if (date) {
    ok &= fprintf(f, "\\date{%s}\n", date) >= 0;
  }

  ok &= fprintf(f, "\\begin{document}\n") >= 0;
  ok &= fprintf(f, "\\maketitle\n") >= 0;

  if (abstract && abstract[0]) {
    ok &= fprintf(f, "\\begin{abstract}\n%s\n\\end{abstract}\n", abstract) >= 0;
  }

  if (sections) {
    for (int i = 0; sections[i] != NULL; i++) {
      ok &= fprintf(f, "%s\n\n", sections[i]) >= 0;
    }
  }

  ok &= fprintf(f, "\\end{document}\n") >= 0;
  int close_rc = fclose(f);

  return (ok && close_rc == 0) ? 0 : -1;
}

/* Table generation */
int latex_generate_table(const char *texpath, const char *const *const *data,
                         int rows, int cols, const char *col_format,
                         const char *caption, const char *label,
                         const char *placement) {
  if (!texpath || rows <= 0 || cols <= 0 || !data) {
    return -1;
  }

  FILE *f = fopen(texpath, "w");
  if (!f) {
    return -1;
  }

  // If no format given, create "|c|c|...|"
  char format[256] = {0};
  if (!col_format) {
    char *p = format;
    *p++ = '|';
    for (int i = 0; i < cols && (size_t)(p - format) < sizeof(format) - 2;
         i++) {
      *p++ = 'c';
      *p++ = '|';
    }

    *p = '\0';
  }

  int ok = 1;
  ok &= fprintf(f, "\\begin{table}[%s]\n\\centering\n",
                placement ? placement : "h") >= 0;
  ok &= fprintf(f, "\\begin{tabular}{%s}\n\\hline\n",
                col_format ? col_format : format) >= 0;

  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      ok &= fprintf(f, "%s", data[r][c] ? data[r][c] : "") >= 0;
      if (c < cols - 1) {
        ok &= fprintf(f, " & ") >= 0;
      }
    }

    ok &= fprintf(f, " \\\\ \\hline\n") >= 0;
  }

  ok &= fprintf(f, "\\end{tabular}\n") >= 0;
  if (caption) {
    ok &= fprintf(f, "\\caption{%s}\n", caption) >= 0;
  }

  if (label) {
    ok &= fprintf(f, "\\label{%s}\n", label) >= 0;
  }
  ok &= fprintf(f, "\\end{table}\n") >= 0;

  int close_rc = fclose(f);

  return (ok && close_rc == 0) ? 0 : -1;
}

/* Math helpers */
char *latex_inline_math(const char *expr) {
  if (!expr) {
    return NULL;
  }

  size_t len = strlen(expr) + 4; // $...$\0
  char *res = malloc(len);
  if (!res) {
    return NULL;
  }

  snprintf(res, len, "$%s$", expr);

  return res;
}

char *latex_display_math(const char *expr) {
  if (!expr) {
    return NULL;
  }

  size_t len = strlen(expr) + 8; // \[ ... \]\0
  char *res = malloc(len);
  if (!res) {
    return NULL;
  }

  snprintf(res, len, "\\[ %s \\]", expr);

  return res;
}

char *latex_matrix(const char *type, const char *const *const *data, int rows,
                   int cols) {
  if (!data || rows <= 0 || cols <= 0) {
    return NULL;
  }

  if (!type) {
    type = "bmatrix";
  }

  size_t est = 2 * strlen(type) + 20; // \begin{type} ... \end{type}
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      est += data[i][j] ? strlen(data[i][j]) : 1;
      est += 4; // " & " or " \\ " separator
    }
  }

  char *res = malloc(est);
  if (!res) {
    return NULL;
  }
  size_t pos = 0;

  pos += (size_t)snprintf(res + pos, est - pos, "\\begin{%s}", type);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      pos += (size_t)snprintf(res + pos, est - pos, "%s",
                              data[i][j] ? data[i][j] : "0");
      if (j < cols - 1) {
        pos += (size_t)snprintf(res + pos, est - pos, " & ");
      }
    }

    if (i < rows - 1) {
      pos += (size_t)snprintf(res + pos, est - pos, " \\\\ ");
    }
  }

  snprintf(res + pos, est - pos, "\\end{%s}", type);

  return res;
}

void latex_clean_temp(void) {
  const char *exts[] = {"tex", "aux", "log", "pdf"};
  for (int i = 0; i < g_tracked_temp_count; i++) {
    char path[192];
    for (size_t e = 0; e < sizeof(exts) / sizeof(exts[0]); e++) {
      snprintf(path, sizeof(path), "/tmp/%.63s.%s", g_tracked_temp_basenames[i],
               exts[e]);
      remove(path);
    }

    snprintf(path, sizeof(path), "/tmp/%.63s_out-1.png",
             g_tracked_temp_basenames[i]);
    remove(path);
  }

  g_tracked_temp_count = 0;

  remove("/tmp/qmc_eq.tex");
  remove("/tmp/qmc_eq.aux");
  remove("/tmp/qmc_eq.log");
  remove("/tmp/qmc_eq.pdf");
  remove("/tmp/qmc_eq_out-1.png");
}
