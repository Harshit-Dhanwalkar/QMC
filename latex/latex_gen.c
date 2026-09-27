#include "latex_gen.h"
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
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

/*
 * Base path buffers hold "<root>/<name>": root is [512], plus a separator,
 * plus longest base name ("equation.tex" = 12 chars). 640 is comfortably
 * above that maximum (512 + 1 + 12 + 1 = 526)
 */
#define LATEX_TMP_PATH_MAX 640

/*
 * Derived path buffers hold "<base><suffix>": png prefix plus a  suffix such as
 * "-1.png" (6 chars). 32 bytes of slack covers every current  case ("-1.png",
 * ".aux", ".aux"'s directory-suffixed form, etc)
 */
#define LATEX_TMP_DERIVED_MAX (LATEX_TMP_PATH_MAX + 32)

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

static char *dupstr(const char *s) {
  if (!s) {
    return NULL;
  }

  size_t len = strlen(s) + 1;
  char *copy = malloc(len);
  if (copy) {
    memcpy(copy, s, len);
  }

  return copy;
}

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
  char tex[LATEX_TMP_PATH_MAX];
  char pdf[LATEX_TMP_PATH_MAX];
  char log[LATEX_TMP_PATH_MAX];
  char png_prefix[LATEX_TMP_PATH_MAX];
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
  int n = snprintf(job->png_prefix, sizeof job->png_prefix, "%s%sequation",
                   job->root, QMC_SEP_STR);
  if (n < 0 || (size_t)n >= sizeof job->png_prefix) {
    return -1; /* job_create failure */
  }

  return 0;
}

static void job_destroy(latex_job_t *job) {
  if (!job || !job->root[0]) {
    return;
  }

  remove(job->tex);
  remove(job->pdf);
  remove(job->log);

  char aux[LATEX_TMP_DERIVED_MAX];
  snprintf(aux, sizeof aux, "%s%sequation.aux", job->root, QMC_SEP_STR);
  remove(aux);

  char out1[LATEX_TMP_DERIVED_MAX];
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

latex_render_options_t latex_render_options_default(void) {
  latex_render_options_t o = {
      .compiler = NULL,
      .dpi = 200,
      .use_amsmath = 1,
      .use_amssymb = 1,
      .use_physics = 1,
      .use_bm = 1,
      .keep_temp_files = 0,
  };

  return o;
}

// Shared .tex-writing logic for latex_render_to_png/pdf
static int
write_standalone_equation_tex_ex(const char *texpath, const char *expr,
                                 const latex_render_options_t *options) {
  latex_string_t s;
  latex_string_init(&s);

  int ok = 1;
  ok &=
      latex_string_append(
          &s, "\\documentclass[12pt,preview]{standalone}\n\\usepackage{") == 0;

  const char *pkgs[4];
  int npkgs = 0;
  if (options->use_amsmath) {
    pkgs[npkgs++] = "amsmath";
  }
  if (options->use_amssymb) {
    pkgs[npkgs++] = "amssymb";
  }
  if (options->use_physics) {
    pkgs[npkgs++] = "physics";
  }
  if (options->use_bm) {
    pkgs[npkgs++] = "bm";
  }

  for (int i = 0; i < npkgs; i++) {
    if (i > 0) {
      ok &= latex_string_append_char(&s, ',') == 0;
    }

    ok &= latex_string_append(&s, pkgs[i]) == 0;
  }

  ok &= latex_string_append(&s, "}\n\\begin{document}\n\\[ ") == 0;
  ok &= latex_string_append(&s, expr) == 0;
  ok &= latex_string_append(&s, " \\]\n\\end{document}\n") == 0;

  if (!ok) {
    latex_string_free(&s);

    return -1;
  }

  char *content = latex_string_detach(&s);
  if (!content) {
    return -1;
  }

  int rc = write_file(texpath, content);
  free(content);

  return rc;
}

// Diagnostics - slurp a compiler log into g_last_error
static void capture_log(const char *log_path, const char *tag) {
  FILE *f = fopen(log_path, "r");
  if (!f) {
    set_last_error("%s failed (no log available)", tag);

    return;
  }

  char buf[LATEX_TMP_PATH_MAX];
  size_t n = fread(buf, 1, sizeof buf - 1, f);
  buf[n] = '\0';

  fclose(f);

  set_last_error("%s failed:\n%s", tag, buf);
}

// Rendering functions (compiler + pdftoppm)
int latex_render_to_png_ex(const char *expr, const char *outpath,
                           const latex_render_options_t *options) {
  if (!expr || !outpath) {
    set_last_error("latex_render_to_png: NULL argument");

    return -1;
  }

  if (has_unsafe_shell_metachars(outpath)) {
    set_last_error("latex_render_to_png: unsafe output path");

    return -5;
  }

  latex_render_options_t default_opts = latex_render_options_default();
  if (!options) {
    options = &default_opts;
  }

  const char *compiler =
      options->compiler ? options->compiler : QMC_LATEX_COMPILER;
  int dpi = (options->dpi > 0) ? options->dpi : 200;

  if (!check_tool_on_path(compiler)) {
    set_last_error("latex_render_to_png: compiler not available: %s", compiler);

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

  if (write_standalone_equation_tex_ex(job.tex, expr, options) != 0) {
    set_last_error("latex_render_to_png: could not write .tex");
    if (!options->keep_temp_files) {
      job_destroy(&job);
    }

    return -1;
  }

  {
    char *argv[] = {(char *)compiler, "-interaction=batchmode",
                    "-halt-on-error", job.tex, NULL};
    if (run_tool(compiler, argv, job.root, NULL, NULL) != 0) {
      capture_log(job.log, "latex_render_to_png: compile");
      if (!options->keep_temp_files) {
        job_destroy(&job);
      }

      return -2;
    }
  }

  {
    char dpi_str[16];
    snprintf(dpi_str, sizeof dpi_str, "%d", dpi);

    char *argv[] = {(char *)QMC_PDFTOPPM, "-r", dpi_str, "-png", job.pdf,
                    job.png_prefix,       NULL};
    if (run_tool(QMC_PDFTOPPM, argv, NULL, NULL, NULL) != 0) {
      set_last_error("latex_render_to_png: pdftoppm failed");

      if (!options->keep_temp_files) {
        job_destroy(&job);
      }

      return -3;
    }
  }

  {
    char produced[LATEX_TMP_DERIVED_MAX];
    snprintf(produced, sizeof produced, "%s-1.png", job.png_prefix);
    if (rename(produced, outpath) != 0) {
      set_last_error("latex_render_to_png: rename failed: %s", strerror(errno));
      if (!options->keep_temp_files) {
        job_destroy(&job);
      }

      return -3;
    }
  }

  if (!options->keep_temp_files) {
    job_destroy(&job);
  } else {
    set_last_error("latex_render_to_png: succeeded; temp dir kept at %s",
                   job.root);
  }

  return 0;
}

int latex_render_to_png(const char *expr, const char *outpath) {
  return latex_render_to_png_ex(expr, outpath, NULL);
}

int latex_render_to_pdf_ex(const char *expr, const char *outpath,
                           const latex_render_options_t *options) {
  if (!expr || !outpath) {
    set_last_error("latex_render_to_pdf: NULL argument");

    return -1;
  }

  if (has_unsafe_shell_metachars(outpath)) {
    set_last_error("latex_render_to_pdf: unsafe output path");

    return -5;
  }

  latex_render_options_t default_opts = latex_render_options_default();
  if (!options) {
    options = &default_opts;
  }

  const char *compiler =
      options->compiler ? options->compiler : QMC_LATEX_COMPILER;

  if (!check_tool_on_path(compiler)) {
    set_last_error("latex_render_to_pdf: compiler not available: %s", compiler);

    return -4;
  }

  latex_job_t job;
  if (job_create(&job) != 0) {
    set_last_error("latex_render_to_pdf: temp dir creation failed");

    return -1;
  }

  if (write_standalone_equation_tex_ex(job.tex, expr, options) != 0) {
    set_last_error("latex_render_to_pdf: could not write .tex");

    if (!options->keep_temp_files) {
      job_destroy(&job);
    }

    return -1;
  }

  {
    char *argv[] = {(char *)compiler, "-interaction=batchmode",
                    "-halt-on-error", job.tex, NULL};
    if (run_tool(compiler, argv, job.root, NULL, NULL) != 0) {
      capture_log(job.log, "latex_render_to_pdf: compile");

      if (!options->keep_temp_files) {
        job_destroy(&job);
      }

      return -2;
    }
  }

  if (rename(job.pdf, outpath) != 0) {
    set_last_error("latex_render_to_pdf: rename failed: %s", strerror(errno));

    if (!options->keep_temp_files) {
      job_destroy(&job);
    }

    return -3;
  }

  if (!options->keep_temp_files) {
    job_destroy(&job);
  } else {
    set_last_error("latex_render_to_pdf: succeeded; temp dir kept at %s",
                   job.root);
  }

  return 0;
}

int latex_render_to_pdf(const char *expr, const char *outpath) {
  return latex_render_to_pdf_ex(expr, outpath, NULL);
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

// Dynamic string builder
void latex_string_init(latex_string_t *s) {
  if (!s) {
    return;
  }

  s->data = NULL;
  s->length = 0;
  s->capacity = 0;
}

static int latex_string_reserve(latex_string_t *s, size_t extra) {
  size_t needed = s->length + extra + 1; // +1 for NUL terminator
  if (needed <= s->capacity) {
    return 0;
  }

  size_t new_cap = s->capacity ? s->capacity * 2 : 64;
  while (new_cap < needed) {
    new_cap *= 2;
  }

  char *nd = realloc(s->data, new_cap);
  if (!nd) {
    return -1;
  }

  s->data = nd;
  s->capacity = new_cap;

  return 0;
}

int latex_string_append(latex_string_t *s, const char *text) {
  if (!s || !text) {
    return -1;
  }

  size_t tlen = strlen(text);
  if (latex_string_reserve(s, tlen) != 0) {
    return -1;
  }

  memcpy(s->data + s->length, text, tlen);
  s->length += tlen;
  s->data[s->length] = '\0';

  return 0;
}

int latex_string_append_char(latex_string_t *s, char c) {
  if (!s) {
    return -1;
  }

  if (latex_string_reserve(s, 1) != 0) {
    return -1;
  }

  s->data[s->length++] = c;
  s->data[s->length] = '\0';

  return 0;
}

int latex_string_appendf(latex_string_t *s, const char *fmt, ...) {
  if (!s || !fmt) {
    return -1;
  }

  va_list ap;
  va_start(ap, fmt);
  va_list ap2;
  va_copy(ap2, ap);
  int needed = vsnprintf(NULL, 0, fmt, ap);
  va_end(ap);

  if (needed < 0) {
    va_end(ap2);

    return -1;
  }

  if (latex_string_reserve(s, (size_t)needed) != 0) {
    va_end(ap2);

    return -1;
  }

  vsnprintf(s->data + s->length, (size_t)needed + 1, fmt, ap2);
  va_end(ap2);
  s->length += (size_t)needed;

  return 0;
}

char *latex_string_detach(latex_string_t *s) {
  if (!s) {
    return NULL;
  }

  char *result = s->data;
  if (!result) {
    result = malloc(1);
    if (result) {
      result[0] = '\0';
    }
  }

  s->data = NULL;
  s->length = 0;
  s->capacity = 0;

  return result;
}

void latex_string_free(latex_string_t *s) {
  if (!s) {
    return;
  }

  free(s->data);
  s->data = NULL;
  s->length = 0;
  s->capacity = 0;
}

// Text escaping
char *latex_escape_text(const char *text) {
  if (!text) {
    return NULL;
  }

  latex_string_t s;
  latex_string_init(&s);

  for (const char *p = text; *p; p++) {
    int rc;

    switch (*p) {
    case '&':
      rc = latex_string_append(&s, "\\&");
      break;
    case '%':
      rc = latex_string_append(&s, "\\%");
      break;
    case '$':
      rc = latex_string_append(&s, "\\$");
      break;
    case '#':
      rc = latex_string_append(&s, "\\#");
      break;
    case '_':
      rc = latex_string_append(&s, "\\_");
      break;
    case '{':
      rc = latex_string_append(&s, "\\{");
      break;
    case '}':
      rc = latex_string_append(&s, "\\}");
      break;
    case '~':
      rc = latex_string_append(&s, "\\textasciitilde{}");
      break;
    case '^':
      rc = latex_string_append(&s, "\\textasciicircum{}");
      break;
    case '\\':
      rc = latex_string_append(&s, "\\textbackslash{}");
      break;
    default:
      rc = latex_string_append_char(&s, *p);
      break;
    }

    if (rc != 0) {
      latex_string_free(&s);

      return NULL;
    }
  }

  return latex_string_detach(&s);
}

int latex_write_escaped_text(FILE *f, const char *text) {
  if (!f || !text) {
    return -1;
  }

  char *escaped = latex_escape_text(text);
  if (!escaped) {
    return -1;
  }

  int rc = fputs(escaped, f);
  free(escaped);

  return (rc >= 0) ? 0 : -1;
}

int latex_validate_label(const char *label) {
  if (!label || !label[0]) {
    return 0;
  }

  for (const char *p = label; *p; p++) {
    if (!(isalnum((unsigned char)*p) || *p == ':' || *p == '_' || *p == '-' ||
          *p == '.')) {
      return 0;
    }
  }

  return 1;
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
  } else if (strcmp(type, "bmatrix") != 0 && strcmp(type, "pmatrix") != 0 &&
             strcmp(type, "Bmatrix") != 0 && strcmp(type, "vmatrix") != 0 &&
             strcmp(type, "Vmatrix") != 0) {
    return NULL;
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

// Table builder
typedef struct {
  char *text;
  int escape;
} latex_table_cell_t;

struct latex_table {
  size_t rows;
  size_t cols;
  latex_table_cell_t *cells; // rows*cols, row-major
  char **headers;            // cols entries; only rendered if has_header
  int has_header;
  char *caption;
  char *label;
  char *col_format;
  char *placement;
  latex_table_style_t style;
};

latex_table_t *latex_table_create(size_t rows, size_t cols) {
  if (rows == 0 || cols == 0 || rows > SIZE_MAX / cols) {
    return NULL;
  }

  latex_table_t *t = calloc(1, sizeof *t);
  if (!t) {
    return NULL;
  }

  t->rows = rows;
  t->cols = cols;
  t->style = LATEX_TABLE_STYLE_CLASSIC;

  t->cells = calloc(rows * cols, sizeof *t->cells);
  t->headers = calloc(cols, sizeof *t->headers);
  if (!t->cells || !t->headers) {
    latex_table_free(t);

    return NULL;
  }

  return t;
}

void latex_table_free(latex_table_t *t) {
  if (!t) {
    return;
  }

  if (t->cells) {
    for (size_t i = 0; i < t->rows * t->cols; i++) {
      free(t->cells[i].text);
    }

    free(t->cells);
  }

  if (t->headers) {
    for (size_t c = 0; c < t->cols; c++) {
      free(t->headers[c]);
    }

    free(t->headers);
  }

  free(t->caption);
  free(t->label);
  free(t->col_format);
  free(t->placement);
  free(t);
}

int latex_table_set(latex_table_t *t, size_t row, size_t col, const char *value,
                    int escape) {
  if (!t || !value || row >= t->rows || col >= t->cols) {
    return -1;
  }

  char *copy = dupstr(value);
  if (!copy) {
    return -1;
  }

  latex_table_cell_t *cell = &t->cells[row * t->cols + col];
  free(cell->text);
  cell->text = copy;
  cell->escape = escape;

  return 0;
}

int latex_table_set_header(latex_table_t *t, size_t col, const char *text) {
  if (!t || !text || col >= t->cols) {
    return -1;
  }

  char *copy = dupstr(text);
  if (!copy) {
    return -1;
  }

  free(t->headers[col]);
  t->headers[col] = copy;
  t->has_header = 1;

  return 0;
}

int latex_table_set_caption(latex_table_t *t, const char *caption) {
  if (!t) {
    return -1;
  }

  char *copy = caption ? dupstr(caption) : NULL;
  if (caption && !copy) {
    return -1;
  }

  free(t->caption);

  t->caption = copy;

  return 0;
}

int latex_table_set_label(latex_table_t *t, const char *label) {
  if (!t) {
    return -1;
  }

  if (label && !latex_validate_label(label)) {
    return -1;
  }

  char *copy = label ? dupstr(label) : NULL;
  if (label && !copy) {
    return -1;
  }

  free(t->label);

  t->label = copy;

  return 0;
}

int latex_table_set_style(latex_table_t *t, latex_table_style_t style) {
  if (!t) {
    return -1;
  }

  t->style = style;

  return 0;
}

int latex_table_set_col_format(latex_table_t *t, const char *col_format) {
  if (!t) {
    return -1;
  }

  char *copy = col_format ? dupstr(col_format) : NULL;
  if (col_format && !copy) {
    return -1;
  }

  free(t->col_format);

  t->col_format = copy;

  return 0;
}

int latex_table_set_placement(latex_table_t *t, const char *placement) {
  if (!t) {
    return -1;
  }

  char *copy = placement ? dupstr(placement) : NULL;
  if (placement && !copy) {
    return -1;
  }

  free(t->placement);

  t->placement = copy;

  return 0;
}

// NOTE: Shared by latex_table_write() and latex_document_add_table() - appends
// table's LaTeX to `out` rather than writing a file directly, so a table can be
// embedded in a document's body
static int latex_table_render(latex_table_t *t, latex_string_t *out) {
  if (!t || !out) {
    return -1;
  }

  int ok = 1;
  int booktabs = (t->style == LATEX_TABLE_STYLE_BOOKTABS);

  ok &= latex_string_appendf(out, "\\begin{table}[%s]\n\\centering\n",
                             t->placement ? t->placement : "h") == 0;

  char generated_format[256];
  const char *col_format = t->col_format;
  if (!col_format) {
    char *p = generated_format;
    const char *end = generated_format + sizeof generated_format - 2;

    if (!booktabs) {
      *p++ = '|';
    }

    for (size_t i = 0; i < t->cols && p < end; i++) {
      *p++ = 'c';
      if (!booktabs) {
        *p++ = '|';
      }
    }

    *p = '\0';

    col_format = generated_format;
  }

  ok &= latex_string_appendf(out, "\\begin{tabular}{%s}\n", col_format) == 0;
  ok &= latex_string_append(out, booktabs ? "\\toprule\n" : "\\hline\n") == 0;

  if (t->has_header) {
    for (size_t c = 0; c < t->cols; c++) {
      if (c > 0) {
        ok &= latex_string_append(out, " & ") == 0;
      }

      char *escaped = latex_escape_text(t->headers[c] ? t->headers[c] : "");
      if (!escaped) {
        return -1;
      }

      if (booktabs) {
        ok &= latex_string_appendf(out, "\\textbf{%s}", escaped) == 0;
      } else {
        ok &= latex_string_append(out, escaped) == 0;
      }

      free(escaped);
    }

    ok &= latex_string_append(out, " \\\\\n") == 0;
    ok &= latex_string_append(out, booktabs ? "\\midrule\n" : "\\hline\n") == 0;
  }

  for (size_t r = 0; r < t->rows; r++) {
    for (size_t c = 0; c < t->cols; c++) {
      if (c > 0) {
        ok &= latex_string_append(out, " & ") == 0;
      }

      latex_table_cell_t *cell = &t->cells[r * t->cols + c];
      const char *text = cell->text ? cell->text : "";

      if (cell->escape) {
        char *escaped = latex_escape_text(text);
        if (!escaped) {
          return -1;
        }

        ok &= latex_string_append(out, escaped) == 0;

        free(escaped);
      } else {
        ok &= latex_string_append(out, text) == 0;
      }
    }

    ok &= latex_string_append(out, " \\\\\n") == 0;
    if (!booktabs) {
      ok &= latex_string_append(out, "\\hline\n") == 0;
    }
  }

  if (booktabs) {
    ok &= latex_string_append(out, "\\bottomrule\n") == 0;
  }

  ok &= latex_string_append(out, "\\end{tabular}\n") == 0;

  if (t->caption) {
    char *escaped = latex_escape_text(t->caption);
    if (!escaped) {
      return -1;
    }

    ok &= latex_string_appendf(out, "\\caption{%s}\n", escaped) == 0;

    free(escaped);
  }

  if (t->label) {
    ok &= latex_string_appendf(out, "\\label{%s}\n", t->label) == 0;
  }

  ok &= latex_string_append(out, "\\end{table}\n") == 0;

  return ok ? 0 : -1;
}

int latex_table_write(latex_table_t *t, const char *texpath) {
  if (!t || !texpath) {
    return -1;
  }

  latex_string_t s;
  latex_string_init(&s);

  if (latex_table_render(t, &s) != 0) {
    latex_string_free(&s);

    return -1;
  }

  char *content = latex_string_detach(&s);
  if (!content) {
    return -1;
  }

  int rc = write_file(texpath, content);

  free(content);

  return rc;
}

// Document builder
struct latex_document {
  char *document_class;
  char *class_options;
  char *title;
  char *author;
  char *date;
  char *abstract;
  latex_string_t packages;
  latex_string_t body;
};

latex_document_t *latex_document_create(const char *title) {
  latex_document_t *doc = calloc(1, sizeof *doc);
  if (!doc) {
    return NULL;
  }

  doc->document_class = dupstr("article");
  if (!doc->document_class) {
    free(doc);

    return NULL;
  }

  if (title) {
    doc->title = dupstr(title);
    if (!doc->title) {
      free(doc->document_class);
      free(doc);

      return NULL;
    }
  }

  latex_string_init(&doc->packages);
  latex_string_init(&doc->body);

  return doc;
}

void latex_document_free(latex_document_t *doc) {
  if (!doc) {
    return;
  }

  free(doc->document_class);
  free(doc->class_options);
  free(doc->title);
  free(doc->author);
  free(doc->date);
  free(doc->abstract);
  latex_string_free(&doc->packages);
  latex_string_free(&doc->body);
  free(doc);
}

int latex_document_set_author(latex_document_t *doc, const char *author) {
  if (!doc) {
    return -1;
  }

  char *copy = author ? dupstr(author) : NULL;
  if (author && !copy) {
    return -1;
  }

  free(doc->author);

  doc->author = copy;

  return 0;
}

int latex_document_set_date(latex_document_t *doc, const char *date) {
  if (!doc) {
    return -1;
  }

  char *copy = date ? dupstr(date) : NULL;
  if (date && !copy) {
    return -1;
  }

  free(doc->date);

  doc->date = copy;

  return 0;
}

int latex_document_set_abstract(latex_document_t *doc, const char *abstract) {
  if (!doc) {
    return -1;
  }

  char *copy = abstract ? dupstr(abstract) : NULL;
  if (abstract && !copy) {
    return -1;
  }

  free(doc->abstract);

  doc->abstract = copy;

  return 0;
}

int latex_document_set_class(latex_document_t *doc, const char *document_class,
                             const char *options) {
  if (!doc || !document_class) {
    return -1;
  }

  char *dc = dupstr(document_class);
  if (!dc) {
    return -1;
  }

  char *opts = options ? dupstr(options) : NULL;
  if (options && !opts) {
    free(dc);

    return -1;
  }

  free(doc->document_class);
  free(doc->class_options);

  doc->document_class = dc;
  doc->class_options = opts;

  return 0;
}

int latex_document_add_package(latex_document_t *doc, const char *package,
                               const char *options) {
  if (!doc || !package) {
    return -1;
  }

  if (options) {
    return latex_string_appendf(&doc->packages, "\\usepackage[%s]{%s}\n",
                                options, package);
  }

  return latex_string_appendf(&doc->packages, "\\usepackage{%s}\n", package);
}

int latex_document_add_section(latex_document_t *doc, const char *title) {
  if (!doc || !title) {
    return -1;
  }

  char *escaped = latex_escape_text(title);
  if (!escaped) {
    return -1;
  }

  int rc = latex_string_appendf(&doc->body, "\\section{%s}\n", escaped);

  free(escaped);

  return rc;
}

int latex_document_add_subsection(latex_document_t *doc, const char *title) {
  if (!doc || !title) {
    return -1;
  }

  char *escaped = latex_escape_text(title);
  if (!escaped) {
    return -1;
  }

  int rc = latex_string_appendf(&doc->body, "\\subsection{%s}\n", escaped);

  free(escaped);

  return rc;
}

int latex_document_add_paragraph(latex_document_t *doc, const char *text) {
  if (!doc || !text) {
    return -1;
  }

  char *escaped = latex_escape_text(text);
  if (!escaped) {
    return -1;
  }

  int rc = latex_string_appendf(&doc->body, "%s\n\n", escaped);

  free(escaped);

  return rc;
}

int latex_document_add_equation(latex_document_t *doc, const char *equation) {
  if (!doc || !equation) {
    return -1;
  }

  return latex_string_appendf(&doc->body, "\\[ %s \\]\n\n", equation);
}

int latex_document_add_raw(latex_document_t *doc, const char *latex) {
  if (!doc || !latex) {
    return -1;
  }

  return latex_string_appendf(&doc->body, "%s\n", latex);
}

int latex_document_add_figure_ex(latex_document_t *doc,
                                 const latex_figure_t *figure) {
  if (!doc || !figure || !figure->path) {
    return -1;
  }

  if (figure->label && !latex_validate_label(figure->label)) {
    return -1;
  }

  double width = (figure->width_fraction > 0) ? figure->width_fraction : 0.8;
  const char *placement = figure->placement ? figure->placement : "h";

  int ok = 1;
  ok &= latex_string_appendf(&doc->body,
                             "\\begin{figure}[%s]\\centering\n"
                             "  \\includegraphics[width=%.3g\\textwidth]{%s}\n",
                             placement, width, figure->path) == 0;

  if (figure->caption) {
    char *escaped = latex_escape_text(figure->caption);
    if (!escaped) {
      return -1;
    }

    ok &= latex_string_appendf(&doc->body, "  \\caption{%s}\n", escaped) == 0;

    free(escaped);
  }

  if (figure->label) {
    ok &=
        latex_string_appendf(&doc->body, "  \\label{%s}\n", figure->label) == 0;
  }

  ok &= latex_string_append(&doc->body, "\\end{figure}\n\n") == 0;

  return ok ? 0 : -1;
}

int latex_document_add_figure(latex_document_t *doc, const char *path,
                              const char *caption, const char *label) {
  latex_figure_t fig = {
      .path = path,
      .caption = caption,
      .label = label,
      .width_fraction = 0.0,
      .placement = NULL,
  };

  return latex_document_add_figure_ex(doc, &fig);
}

int latex_document_add_table(latex_document_t *doc, latex_table_t *table) {
  if (!doc || !table) {
    return -1;
  }

  int rc = latex_table_render(table, &doc->body);
  if (rc == 0) {
    latex_string_append(&doc->body, "\n");
    latex_table_free(table);
  }

  return rc;
}

int latex_document_write(latex_document_t *doc, const char *texpath) {
  if (!doc || !texpath) {
    return -1;
  }

  latex_string_t s;
  latex_string_init(&s);

  int ok = 1;

  if (doc->class_options) {
    ok &= latex_string_appendf(&s, "\\documentclass[%s]{%s}\n",
                               doc->class_options, doc->document_class) == 0;
  } else {
    ok &= latex_string_appendf(&s, "\\documentclass{%s}\n",
                               doc->document_class) == 0;
  }

  if (doc->packages.length > 0) {
    ok &= latex_string_append(&s, doc->packages.data) == 0;
  }

  if (doc->title) {
    char *escaped = latex_escape_text(doc->title);
    if (!escaped) {
      latex_string_free(&s);

      return -1;
    }

    ok &= latex_string_appendf(&s, "\\title{%s}\n", escaped) == 0;

    free(escaped);
  }

  if (doc->author) {
    char *escaped = latex_escape_text(doc->author);
    if (!escaped) {
      latex_string_free(&s);

      return -1;
    }

    ok &= latex_string_appendf(&s, "\\author{%s}\n", escaped) == 0;

    free(escaped);
  }

  if (doc->date) {
    char *escaped = latex_escape_text(doc->date);
    if (!escaped) {
      latex_string_free(&s);

      return -1;
    }

    ok &= latex_string_appendf(&s, "\\date{%s}\n", escaped) == 0;

    free(escaped);
  }

  ok &= latex_string_append(&s, "\\begin{document}\n") == 0;

  if (doc->title || doc->author || doc->date) {
    ok &= latex_string_append(&s, "\\maketitle\n") == 0;
  }

  if (doc->abstract) {
    char *escaped = latex_escape_text(doc->abstract);
    if (!escaped) {
      latex_string_free(&s);

      return -1;
    }

    ok &= latex_string_appendf(&s, "\\begin{abstract}\n%s\n\\end{abstract}\n",
                               escaped) == 0;

    free(escaped);
  }

  if (doc->body.length > 0) {
    ok &= latex_string_append(&s, doc->body.data) == 0;
  }

  ok &= latex_string_append(&s, "\\end{document}\n") == 0;

  if (!ok) {
    latex_string_free(&s);

    return -1;
  }

  char *content = latex_string_detach(&s);
  if (!content) {
    return -1;
  }

  int rc = write_file(texpath, content);

  free(content);

  return rc;
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
