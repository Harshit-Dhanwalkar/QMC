#include "latex_gen.h"
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Platform setup
#ifdef _WIN32

#include <direct.h>
#include <fcntl.h>
#include <io.h>
#include <process.h> // _getpid()
#include <sys/stat.h>

// #define getpid _getpid
// #define QMC_PATH_SEP "\\"
#define QMC_SEP '\\'
#define QMC_SEP_STR "\\"
// #define QMC_DEVNULL_REDIRECT "> NUL 2>&1"
// #define QMC_CD_CMD "cd /d" // /d also switches drive letter
// #define QMC_CHECK_TOOL_FMT "where %s >NUL 2>&1"
#define qmc_mkdir(p) _mkdir(p)
#define qmc_rmdir(p) _rmdir(p)

#else

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h> // WIFEXITED/WEXITSTATUS
#include <unistd.h>   // getpid()

// #define QMC_PATH_SEP "/"
#define QMC_SEP '/'
#define QMC_SEP_STR "/"
// #define QMC_DEVNULL_REDIRECT "> /dev/null 2>&1"
// #define QMC_CD_CMD "cd"
// #define QMC_CHECK_TOOL_FMT "command -v %s >/dev/null 2>&1"
#define qmc_mkdir(p) mkdir((p), 0755)
#define qmc_rmdir(p) rmdir(p)

#endif

// NOTE: Requires pdflatex (or lualatex) and pdftoppm (poppler-utils)
//       installed and on PATH.

#ifndef QMC_LATEX_COMPILER
#define QMC_LATEX_COMPILER "pdflatex"
#endif

#define QMC_PDFTOPPM "pdftoppm"
// #define QMC_LATEX_MAX_TRACKED_TEMP 256
//
// static char g_tracked_temp_basenames[QMC_LATEX_MAX_TRACKED_TEMP][64];
// static int g_tracked_temp_count = 0;

// Error state
static char g_last_error[2048];

static void set_last_error(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(g_last_error, sizeof g_last_error, fmt, ap);
  va_end(ap);
}

const char *latex_last_error(void) { return g_last_error; }

// Path Helpers
static const char *tmp_base(void) {
#ifdef _WIN32

  // NOTE: Windows has no fixed /tmp; conventional per-user scratch dir is
  // whatever %TEMP% (or %TMP%) points to, falling back to current directory if
  // neither is set
  static const char *tmp_dir_prefix(void) {
    const char *tmp_dir = getenv("TEMP");
    if (!tmp_dir) {
      tmp_dir = getenv("TMP");
    }

    if (!tmp_dir) {
      tmp_dir = ".";
    }

    return tmp_dir;
  }

#else

  return "/tmp";

#endif
}

// static void make_tmp_path(const char *basename, const char *ext, char *buf,
//                           size_t bufsz) {
//   snprintf(buf, bufsz, "%s" QMC_PATH_SEP "%s.%s", tmp_dir_prefix(),
//   basename,
//            ext);
// }
//
// static void make_unique_basename(char *buf, size_t bufsz) {
//   static long counter = 0;
//   long id = (long)getpid() * 1000003L + counter++;
//   snprintf(buf, bufsz, "qmc_eq_%ld", id);
//
//   if (g_tracked_temp_count < QMC_LATEX_MAX_TRACKED_TEMP) {
//     snprintf(g_tracked_temp_basenames[g_tracked_temp_count++],
//              sizeof(g_tracked_temp_basenames[0]), "%s", buf);
//   }
// }

// // Runs a shell command and reports success/failure correctly
// static int run_system(const char *cmd) {
//   int status = system(cmd);
//   if (status == -1) {
//     return -1;
//   }
//
// #ifdef _WIN32
//
//   // cmd.exe's system() return value IS the child's exit code directly - no
//   // POSIX-style wait-status encoding to unpack.
//   return status == 0 ? 0 : -1;
//
// #else
//
//   if (WIFEXITED(status)) {
//     return WEXITSTATUS(status) == 0 ? 0 : -1;
//   }
//
//   return -1; // terminated by signal, etc.
//
// #endif
// }

static int has_unsafe_shell_metachars(const char *s) {
  if (!s) {
    return 0;
  }

  return strpbrk(s, ";&|`$()<>\"'\\\n") != NULL;
}

/* Runs `command -v <tool>` (POSIX) or `where <tool>` (Windows)
 *
 * Returns 1 if the tool is on PATH */
static int check_tool_on_path(const char *tool) {
  if (!tool) {
    return 0;
  }

  char cmd[256];
  // snprintf(cmd, sizeof(cmd), QMC_CHECK_TOOL_FMT, tool);
  //
  // return run_system(cmd) == 0;

#ifdef _WIN32

  snprintf(cmd, sizeof cmd, "where %s >NUL 2>&1", tool);

#else

  snprintf(cmd, sizeof cmd, "command -v %s >/dev/null 2>&1", tool);

#endif

  /* NOTE: Tool name is a compile-time constant, never user-supplied, so shell
   * here is not a security boundary */
  return system(cmd) == 0;
}

// int latex_compiler_available(void) {
//   return check_tool_on_path(QMC_LATEX_COMPILER) &&
//          check_tool_on_path("pdftoppm");
// }
int latex_compiler_available(void) {
  return check_tool_on_path(QMC_LATEX_COMPILER);
}

int latex_png_converter_available(void) {
  return check_tool_on_path(QMC_PDFTOPPM);
}

int latex_tools_available(void) {
  return latex_compiler_available() && latex_png_converter_available();
}

/*
 * Process execution
 *
 * argv must be NULL-terminated; argv[0] is tool name
 * cwd, stdout_path, stderr_path may be NULL

 * Returns 0 if the process exited with status 0, -1 otherwise
 */
static int run_tool(const char *tool, char *const argv[], const char *cwd,
                    const char *stdout_path, const char *stderr_path) {
#ifdef _WIN32

  char saved_cwd[1024];
  int have_cwd = 0;
  if (cwd && *cwd) {
    if (!_getcwd(saved_cwd, sizeof saved_cwd)) {
      return -1;
    }
    if (_chdir(cwd) != 0) {
      return -1;
    }

    have_cwd = 1;
  }

  int saved_out = -1;
  int saved_err = -1;
  if (stdout_path) {
    saved_out = _dup(_fileno(stdout));
    if (!freopen(stdout_path, "w", stdout)) {
      if (have_cwd) {
        _chdir(saved_cwd);
      }

      return -1;
    }
  }

  if (stderr_path) {
    saved_err = _dup(_fileno(stderr));
    if (!freopen(stderr_path, "w", stderr)) {
      if (saved_out >= 0) {
        _dup2(saved_out, _fileno(stdout));
        _close(saved_out);
      }

      if (have_cwd) {
        _chdir(saved_cwd);
      }

      return -1;
    }
  }

  intptr_t rc = _spawnvp(_P_WAIT, tool, (const char *const *)argv);

  if (saved_out >= 0) {
    _dup2(saved_out, _fileno(stdout));
    _close(saved_out);
  }
  if (saved_err >= 0) {
    _dup2(saved_err, _fileno(stderr));
    _close(saved_err);
  }
  if (have_cwd) {
    _chdir(saved_cwd);
  }

  return (rc == 0) ? 0 : -1;

#else

  pid_t pid = fork();
  if (pid < 0) {
    return -1;
  }

  if (pid == 0) {
    if (cwd && *cwd && chdir(cwd) != 0) {
      _exit(127);
    }

    if (stdout_path) {
      int fd = open(stdout_path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
      if (fd < 0 || dup2(fd, STDOUT_FILENO) < 0) {
        _exit(127);
      }

      close(fd);
    }

    if (stderr_path) {
      int fd = open(stderr_path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
      if (fd < 0 || dup2(fd, STDERR_FILENO) < 0) {
        _exit(127);
      }

      close(fd);
    }

    execvp(tool, argv);
    _exit(127); // exec failed
  }

  int status;
  if (waitpid(pid, &status, 0) < 0) {
    return -1;
  }
  if (!WIFEXITED(status)) {
    return -1;
  }

  return (WEXITSTATUS(status) == 0) ? 0 : -1;

#endif
}

// Per-render temporary job directory
typedef struct {
  char root[512]; /* e.g. /tmp/qmc-latex-XXXXXX */
  char tex[600];
  char pdf[600];
  char log[600];
  char png_prefix[512];
} latex_job_t;

static int job_create(latex_job_t *job) {
  memset(job, 0, sizeof *job);

#ifdef _WIN32

  static long counter = 0;
  long id = (long)_getpid() * 1000003L + counter++;

  snprintf(job->root, sizeof job->root, "%s\\qmc-latex-%ld", tmp_base(), id);
  if (qmc_mkdir(job->root) != 0 && errno != EEXIST) {
    return -1;
  }

#else

  snprintf(job->root, sizeof job->root, "%s/qmc-latex-XXXXXX", tmp_base());
  if (!mkdtemp(job->root)) {
    return -1;
  }

#endif

  snprintf(job->tex, sizeof job->tex, "%s%sequation.tex", job->root,
           QMC_SEP_STR);
  snprintf(job->pdf, sizeof job->pdf, "%s%sequation.pdf", job->root,
           QMC_SEP_STR);
  snprintf(job->log, sizeof job->log, "%s%sequation.log", job->root,
           QMC_SEP_STR);
  snprintf(job->png_prefix, sizeof job->png_prefix, "%s%sequation", job->root,
           QMC_SEP_STR);

  return 0;
}

static void job_destroy(latex_job_t *job) {
  if (!job || !job->root[0]) {
    return;
  }

  remove(job->tex);
  remove(job->pdf);
  remove(job->log);

  char aux[600];
  snprintf(aux, sizeof aux, "%s%sequation.aux", job->root, QMC_SEP_STR);
  remove(aux);

  char out1[600];
  snprintf(out1, sizeof out1, "%s-1.png", job->png_prefix);
  remove(out1);

  qmc_rmdir(job->root);
}

// .tex writing
static int write_file(const char *path, const char *content) {
  FILE *f = fopen(path, "w");
  if (!f) {
    return -1;
  }

  int rc = fputs(content, f);
  int close_rc = fclose(f);

  return (rc >= 0 && close_rc == 0) ? 0 : -1;
}

// Shared .tex-writing logic for latex_render_to_png/pdf
static int write_standalone_equation_tex(const char *texpath,
                                         const char *expr) {
  char buf[8192];

  snprintf(buf, sizeof buf,
           "\\documentclass[12pt,preview]{standalone}\n"
           "\\usepackage{amsmath,amssymb,physics,bm}\n"
           "\\begin{document}\n"
           "\\[ %s \\]\n"
           "\\end{document}\n",
           expr);

  return write_file(texpath, buf);
}

// Diagnostics - slurp a compiler log into g_last_error
static void capture_log(const char *log_path, const char *tag) {
  FILE *f = fopen(log_path, "r");
  if (!f) {
    set_last_error("%s failed (no log available)", tag);

    return;
  }

  char buf[1600];
  size_t n = fread(buf, 1, sizeof buf - 1, f);
  buf[n] = '\0';

  fclose(f);

  set_last_error("%s failed:\n%s", tag, buf);
}

// Rendering functions (compiler + pdftoppm)
int latex_render_to_png(const char *expr, const char *outpath) {
  if (!expr || !outpath) {
    set_last_error("latex_render_to_png: NULL argument");

    return -1;
  }

  if (has_unsafe_shell_metachars(outpath)) {
    set_last_error("latex_render_to_png: unsafe output path");

    return -5;
  }

  if (!latex_compiler_available()) {
    set_last_error("latex_render_to_png: compiler not available");

    return -4;
  }

  if (!latex_png_converter_available()) {
    set_last_error("latex_render_to_png: pdftoppm not available");

    return -4;
  }

  latex_job_t job;
  if (job_create(&job) != 0) {
    set_last_error("latex_render_to_png: temp dir creation failed");

    return -1;
  }

  // char basename[64];
  // make_unique_basename(basename, sizeof(basename));
  //
  // char texpath[128], pdfpath[128], pngprefix[128], pngout[136];
  // make_tmp_path(basename, "tex", texpath, sizeof(texpath));
  // make_tmp_path(basename, "pdf", pdfpath, sizeof(pdfpath));
  //
  // snprintf(pngprefix, sizeof(pngprefix), "%s" QMC_PATH_SEP "%s_out",
  //          tmp_dir_prefix(), basename);
  // snprintf(pngout, sizeof(pngout), "%s-1.png", pngprefix);
  //
  // if (rename(pngout, outpath) != 0) {
  //   return -3;
  // }

  if (write_standalone_equation_tex(job.tex, expr) != 0) {
    set_last_error("latex_render_to_png: could not write .tex");
    job_destroy(&job);

    return -1;
  }

  // char cmd[512];
  // snprintf(cmd, sizeof(cmd),
  //          QMC_CD_CMD
  //          " %s && %s -interaction=batchmode %s " QMC_DEVNULL_REDIRECT,
  //          tmp_dir_prefix(), QMC_LATEX_COMPILER, texpath);
  // if (run_system(cmd) != 0) {
  //   return -2;
  // }
  {
    char *argv[] = {(char *)QMC_LATEX_COMPILER, "-interaction=batchmode",
                    "-halt-on-error", job.tex, NULL};
    if (run_tool(QMC_LATEX_COMPILER, argv, job.root, NULL, NULL) != 0) {
      capture_log(job.log, "latex_render_to_png: compile");
      job_destroy(&job);

      return -2;
    }
  }

  // snprintf(cmd, sizeof(cmd), "pdftoppm -r 200 -png %s %s", pdfpath,
  // pngprefix); if (run_system(cmd) != 0) {
  //   return -3;
  // }
  //
  {
    char *argv[] = {(char *)QMC_PDFTOPPM, "-r", "200", "-png", job.pdf,
                    job.png_prefix,       NULL};
    if (run_tool(QMC_PDFTOPPM, argv, NULL, NULL, NULL) != 0) {
      set_last_error("latex_render_to_png: pdftoppm failed");

      job_destroy(&job);

      return -3;
    }
  }

  {
    char produced[600];
    snprintf(produced, sizeof produced, "%s-1.png", job.png_prefix);
    if (rename(produced, outpath) != 0) {
      set_last_error("latex_render_to_png: rename failed: %s", strerror(errno));
      job_destroy(&job);

      return -3;
    }
  }

  job_destroy(&job);

  return 0;
}

int latex_render_to_pdf(const char *expr, const char *outpath) {
  if (!expr || !outpath) {
    set_last_error("latex_render_to_pdf: NULL argument");

    return -1;
  }

  if (has_unsafe_shell_metachars(outpath)) {
    set_last_error("latex_render_to_pdf: unsafe output path");

    return -5;
  }

  if (!latex_compiler_available()) {
    set_last_error("latex_render_to_pdf: compiler not available");

    return -4;
  }

  // char basename[64];
  // make_unique_basename(basename, sizeof(basename));
  latex_job_t job;
  if (job_create(&job) != 0) {
    set_last_error("latex_render_to_pdf: temp dir creation failed");

    return -1;
  }

  // char texpath[128];
  // char pdfpath[128];
  // make_tmp_path(basename, "tex", texpath, sizeof(texpath));
  // make_tmp_path(basename, "pdf", pdfpath, sizeof(pdfpath));
  //
  // if (write_standalone_equation_tex(texpath, expr) != 0) {
  //   return -1;
  // }
  //
  if (write_standalone_equation_tex(job.tex, expr) != 0) {
    set_last_error("latex_render_to_pdf: could not write .tex");

    job_destroy(&job);

    return -1;
  }

  // char cmd[512];
  // snprintf(cmd, sizeof(cmd),
  //          QMC_CD_CMD
  //          " %s && %s -interaction=batchmode %s " QMC_DEVNULL_REDIRECT,
  //          tmp_dir_prefix(), QMC_LATEX_COMPILER, texpath);
  // if (run_system(cmd) != 0) {
  //   return -2;
  // }
  {
    char *argv[] = {(char *)QMC_LATEX_COMPILER, "-interaction=batchmode",
                    "-halt-on-error", job.tex, NULL};
    if (run_tool(QMC_LATEX_COMPILER, argv, job.root, NULL, NULL) != 0) {
      capture_log(job.log, "latex_render_to_pdf: compile");

      job_destroy(&job);

      return -2;
    }
  }

  // if (rename(pdfpath, outpath) != 0) {
  //   return -3;
  // }
  if (rename(job.pdf, outpath) != 0) {
    set_last_error("latex_render_to_pdf: rename failed: %s", strerror(errno));

    job_destroy(&job);

    return -3;
  }

  job_destroy(&job);

  return 0;
}

int latex_compile(const char *texfile, const char *compiler,
                  const char *working_dir) {
  if (!texfile) {
    set_last_error("latex_compile: NULL texfile");

    return -1;
  }

  if (!compiler) {
    compiler = QMC_LATEX_COMPILER;
  }

  if (has_unsafe_shell_metachars(texfile) ||
      has_unsafe_shell_metachars(compiler) ||
      has_unsafe_shell_metachars(working_dir)) {
    set_last_error("latex_compile: unsafe argument");

    return -5;
  }

  if (!check_tool_on_path(compiler)) {
    set_last_error("latex_compile: compiler not found: %s", compiler);

    return -4;
  }

  char *argv[] = {(char *)compiler, "-interaction=batchmode", "-halt-on-error",
                  (char *)texfile, NULL};

  if (run_tool(compiler, argv, working_dir, NULL, NULL) != 0) {
    set_last_error("latex_compile: compilation failed");

    return -2;
  }

  return 0;
}

// Document / table / equation writers
int latex_write_figure(const char *texpath, const char *figure_pdf,
                       const char *caption, const char *equation) {
  if (!texpath) {
    return -1;
  }

  char buf[8192];
  snprintf(buf, sizeof buf,
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

  return write_file(texpath, buf);
}

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

  char format[256] = {0};
  if (!col_format) {
    char *p = format;
    *p++ = '|';
    for (int i = 0; i < cols && (size_t)(p - format) < sizeof format - 2; i++) {
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

// Math helpers
char *latex_inline_math(const char *expr) {
  if (!expr) {
    return NULL;
  }

  size_t len = strlen(expr) + 4;
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

  size_t len = strlen(expr) + 8;
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

  size_t est = 2 * strlen(type) + 20;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      est += data[i][j] ? strlen(data[i][j]) : 1;
      est += 4;
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

// Cleanup
void latex_clean_temp(void) {
  // const char *exts[] = {"tex", "aux", "log", "pdf"};
  // for (int i = 0; i < g_tracked_temp_count; i++) {
  //   char path[192];
  //
  //   for (size_t e = 0; e < sizeof(exts) / sizeof(exts[0]); e++) {
  //     snprintf(path, sizeof(path), "%s" QMC_PATH_SEP "%.63s.%s",
  //              tmp_dir_prefix(), g_tracked_temp_basenames[i], exts[e]);
  //
  //     remove(path);
  //   }
  //
  //   snprintf(path, sizeof(path), "%s" QMC_PATH_SEP "%.63s_out-1.png",
  //            tmp_dir_prefix(), g_tracked_temp_basenames[i]);
  //
  //   remove(path);
  // }
  //
  // g_tracked_temp_count = 0;
  //
  // char path[192];
  // const char *leftover_exts[] = {"tex", "aux", "log", "pdf"};
  // for (size_t e = 0; e < sizeof(leftover_exts) / sizeof(leftover_exts[0]);
  //      e++) {
  //   snprintf(path, sizeof(path), "%s" QMC_PATH_SEP "qmc_eq.%s",
  //            tmp_dir_prefix(), leftover_exts[e]);
  //
  //   remove(path);
  // }
  //
  // snprintf(path, sizeof(path), "%s" QMC_PATH_SEP "qmc_eq_out-1.png",
  //          tmp_dir_prefix());
  //
  // remove(path);

  /* Every render job now removes its own temporary directory on
   * completion, success or failure */
  // WARN: This function is retained for backward source compatibility only
}
