#include "basis_parser.h"
#include "molecular_integrals.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INITIAL_CAP 16

void basis_set_free_functions(basis_function_t **funcs, int num_funcs) {
  if (!funcs) {
    return;
  }

  for (int i = 0; i < num_funcs; i++) {
    basis_function_free(funcs[i]);
  }

  free((void *)funcs); // NOLINT: explicit cast for multi-level pointer
}

/* ---------------------------------------------------------------------
 * Small line/token helpers
 * ------------------------------------------------------------------- */

/* strtod that also accepts Fortran-style 'D'/'d' exponent markers ("1.0D+00")
 * by rewriting a private copy before handing off to strtod. */
static double parse_double_token(const char *tok) {
  char buf[64];
  size_t i;
  for (i = 0; i < sizeof(buf) - 1 && tok[i] != '\0'; i++) {
    char c = tok[i];

    buf[i] = (c == 'D' || c == 'd') ? 'e' : c;
  }

  buf[i] = '\0';

  return strtod(buf, NULL);
}

static void to_upper_str(char *s) {
  for (; *s; s++) {
    *s = (char)toupper((unsigned char)*s);
  }
}

static int is_blank_or_comment(const char *line) {
  while (*line == ' ' || *line == '\t' || *line == '\r') {
    line++;
  }

  return (*line == '\0' || *line == '\n' || *line == '!' || *line == '#');
}

static int is_separator_line(const char *line) {
  while (*line == ' ' || *line == '\t') {
    line++;
  }

  return (line[0] == '*');
}

// Angular momentum letter -> l. Returns -1 if unrecognized.
static int shell_letter_to_l(char c) {
  switch (toupper((unsigned char)c)) {
  case 'S':
    return 0;
  case 'P':
    return 1;
  case 'D':
    return 2;
  case 'F':
    return 3;
  case 'G':
    return 4;
  default:
    return -1;
  }
}

/* ---------------------------------------------------------------------
 * Growable arrays used only during parsing
 * ------------------------------------------------------------------- */

typedef struct {
  basis_element_t *items;
  int count, cap;
} elem_vec_t;

typedef struct {
  basis_shell_t *items;
  int count, cap;
} shell_vec_t;

static int elem_vec_push(elem_vec_t *vec, basis_element_t elem) {
  if (vec->count == vec->cap) {
    int new_cap = vec->cap == 0 ? 8 : vec->cap * 2;
    basis_element_t *tmp = realloc(vec->items, (size_t)new_cap * sizeof(*tmp));

    if (!tmp) {
      return 0;
    }

    vec->items = tmp;
    vec->cap = new_cap;
  }

  vec->items[vec->count++] = elem;

  return 1;
}

static int shell_vec_push(shell_vec_t *vec, basis_shell_t shell) {
  if (vec->count == vec->cap) {
    int new_cap = vec->cap == 0 ? 4 : vec->cap * 2;
    basis_shell_t *tmp = realloc(vec->items, (size_t)new_cap * sizeof(*tmp));

    if (!tmp) {
      return 0;
    }

    vec->items = tmp;
    vec->cap = new_cap;
  }

  vec->items[vec->count++] = shell;

  return 1;
}

static void shell_free(basis_shell_t *shell) {
  if (!shell) {
    return;
  }
  free(shell->exps);
  free(shell->coef1);
  free(shell->coef2);
}

/* ---------------------------------------------------------------------
 * Parser
 * ------------------------------------------------------------------- */

/*
 * Splits `line` into up to max_tok whitespace-separated tokens.
 *
 * Returns token count.
 */
static int split_tokens(char *line, char **tok, int max_tok) {
  int n = 0;
  char *p = strtok(line, " \t\r\n");

  while (p && n < max_tok) {
    tok[n++] = p;
    p = strtok(NULL, " \t\r\n");
  }

  return n;
}

basis_set_t *basis_set_parse_string(const char *text) {
  if (!text) {
    return NULL;
  }

  // Work on a mutable copy, line by line
  char *buf = strdup(text);
  if (!buf) {
    return NULL;
  }

  elem_vec_t elems = {0};
  int ok = 1;

  char *saveptr_lines = NULL;
  const char *line = strtok_r(buf, "\n", &saveptr_lines);

  char cur_symbol[4] = {0};
  shell_vec_t cur_shells = {0};
  int have_open_element = 0;

  while (line != NULL && ok) {
    if (is_blank_or_comment(line)) {
      line = strtok_r(NULL, "\n", &saveptr_lines);

      continue;
    }

    if (is_separator_line(line)) {
      // "****": close out any open element block
      if (have_open_element) {
        basis_element_t e;
        memset(&e, 0, sizeof(e));
        snprintf(e.element, sizeof(e.element), "%s", cur_symbol);

        e.n_shells = cur_shells.count;
        e.shells = cur_shells.items;
        if (!elem_vec_push(&elems, e)) {
          ok = 0;
        }

        cur_shells.items = NULL;
        cur_shells.count = cur_shells.cap = 0;
        have_open_element = 0;
      }

      line = strtok_r(NULL, "\n", &saveptr_lines);

      continue;
    }

    // Tokenize this line (on a scratch copy; strtok mutates)
    char scratch[256];
    strncpy(scratch, line, sizeof(scratch) - 1);
    scratch[sizeof(scratch) - 1] = '\0';

    char *tok[8];
    int ntok = split_tokens(scratch, tok, 8);
    if (ntok == 0) {
      line = strtok_r(NULL, "\n", &saveptr_lines);

      continue;
    }

    if (!have_open_element) {
      // Expect an element header: "SYMBOL 0"
      if (ntok < 1) {
        ok = 0;
        break;
      }

      memset(cur_symbol, 0, sizeof(cur_symbol));
      strncpy(cur_symbol, tok[0], sizeof(cur_symbol) - 1);

      to_upper_str(cur_symbol);

      have_open_element = 1;
      cur_shells.items = NULL;
      cur_shells.count = cur_shells.cap = 0;
      line = strtok_r(NULL, "\n", &saveptr_lines);

      continue;
    }

    // Otherwise this is a shell header: "TYPE nprim scale"
    if (ntok < 2) {
      ok = 0;
      break;
    }

    char shell_type[3] = {0};
    strncpy(shell_type, tok[0], 2);

    to_upper_str(shell_type);

    int n_prim = (int)strtod(tok[1], NULL);
    if (n_prim <= 0) {
      ok = 0;
      break;
    }

    int is_sp = (strcmp(shell_type, "SP") == 0);
    int l = is_sp ? 0 : shell_letter_to_l(shell_type[0]);
    if (!is_sp && (l < 0 || shell_type[1] != '\0')) {
      // Single-letter shell types only (beyond the special-cased "SP")
      ok = 0;
      break;
    }

    double *exps = malloc((size_t)n_prim * sizeof(double));
    double *coef1 = malloc((size_t)n_prim * sizeof(double));
    double *coef2 = is_sp ? malloc((size_t)n_prim * sizeof(double)) : NULL;
    if (!exps || !coef1 || (is_sp && !coef2)) {
      free(exps);
      free(coef1);
      free(coef2);

      ok = 0;
      break;
    }

    for (int i = 0; i < n_prim && ok; i++) {
      line = strtok_r(NULL, "\n", &saveptr_lines);
      while (line != NULL && is_blank_or_comment(line)) {
        line = strtok_r(NULL, "\n", &saveptr_lines);
      }

      if (!line) {
        ok = 0;
        break;
      }

      char pscratch[256];
      strncpy(pscratch, line, sizeof(pscratch) - 1);
      pscratch[sizeof(pscratch) - 1] = '\0';

      char *ptok[4];
      int pntok = split_tokens(pscratch, ptok, 4);
      if (pntok < (is_sp ? 3 : 2)) {
        ok = 0;
        break;
      }

      exps[i] = parse_double_token(ptok[0]);
      coef1[i] = parse_double_token(ptok[1]);
      if (is_sp) {
        coef2[i] = parse_double_token(ptok[2]);
      }
    }

    if (!ok) {
      free(exps);
      free(coef1);
      free(coef2);

      break;
    }

    basis_shell_t sh;
    memset(&sh, 0, sizeof(sh));
    snprintf(sh.type, sizeof(sh.type), "%s", shell_type);

    sh.l = l;
    sh.n_prim = n_prim;
    sh.exps = exps;
    sh.coef1 = coef1;
    sh.coef2 = coef2;
    if (!shell_vec_push(&cur_shells, sh)) {
      shell_free(&sh);
      ok = 0;
      break;
    }

    line = strtok_r(NULL, "\n", &saveptr_lines);
  }

  // Trailing element block with no terminating "****" line
  if (ok && have_open_element) {
    basis_element_t e;
    memset(&e, 0, sizeof(e));
    snprintf(e.element, sizeof(e.element), "%s", cur_symbol);

    e.n_shells = cur_shells.count;
    e.shells = cur_shells.items;
    if (!elem_vec_push(&elems, e)) {
      ok = 0;
    }

    cur_shells.items = NULL;
  }

  free(buf);

  if (!ok) {
    for (int i = 0; i < cur_shells.count; i++) {
      shell_free(&cur_shells.items[i]);
    }

    if (cur_shells.items) {
      free(cur_shells.items);
    }

    for (int i = 0; i < elems.count; i++) {
      for (int j = 0; j < elems.items[i].n_shells; j++) {
        shell_free(&elems.items[i].shells[j]);
      }

      free(elems.items[i].shells);
    }

    free(elems.items);

    return NULL;
  }

  basis_set_t *bs = malloc(sizeof(*bs));
  if (!bs) {
    for (int i = 0; i < elems.count; i++) {
      for (int j = 0; j < elems.items[i].n_shells; j++) {
        shell_free(&elems.items[i].shells[j]);
      }

      free(elems.items[i].shells);
    }

    free(elems.items);

    return NULL;
  }

  bs->n_elements = elems.count;
  bs->elements = elems.items;

  elems.items = NULL;
  elems.count = elems.cap = 0;

  return bs;
}

basis_set_t *basis_set_parse_file(const char *path) {
  FILE *fp = fopen(path, "rb");
  if (!fp) {
    return NULL;
  }

  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);

    return NULL;
  }

  long size = ftell(fp);
  if (size < 0) {
    fclose(fp);

    return NULL;
  }

  rewind(fp);

  char *buf = malloc((size_t)size + 1);
  if (!buf) {
    fclose(fp);

    return NULL;
  }

  size_t nread = fread(buf, 1, (size_t)size, fp);
  fclose(fp);

  buf[nread] = '\0';
  basis_set_t *bs = basis_set_parse_string(buf);

  free(buf);

  return bs;
}

void basis_set_free(basis_set_t *bs) {
  if (!bs) {
    return;
  }

  for (int i = 0; i < bs->n_elements; i++) {
    for (int j = 0; j < bs->elements[i].n_shells; j++) {
      shell_free(&bs->elements[i].shells[j]);
    }

    free(bs->elements[i].shells);
  }

  free(bs->elements);
  free(bs);
}

const basis_element_t *basis_set_find_element(const basis_set_t *bs,
                                              const char *symbol) {
  if (!bs || !symbol) {
    return NULL;
  }

  char up[4] = {0};
  strncpy(up, symbol, sizeof(up) - 1);
  to_upper_str(up);
  for (int i = 0; i < bs->n_elements; i++) {
    if (strcmp(bs->elements[i].element, up) == 0) {
      return &bs->elements[i];
    }
  }

  return NULL;
}

/* ---------------------------------------------------------------------
 * Shell -> basis_function_t expansion
 * ------------------------------------------------------------------- */

/* Enumerate the (l+1)(l+2)/2 Cartesian (lx,ly,lz) components of angular
 * momentum l, in decreasing-lx canonical order (matches the convention already
 * used by molint_basis_sto3g_li's explicit px,py,pz ordering for l=1). Writes
 * into out[][3], returns the count.
 */
static int cartesian_components(int l_quantum, int out[][3]) {
  int n = 0;
  for (int lx = l_quantum; lx >= 0; lx--) {
    for (int ly = l_quantum - lx; ly >= 0; ly--) {
      int lz = l_quantum - lx - ly;

      out[n][0] = lx;
      out[n][1] = ly;
      out[n][2] = lz;

      n++;
    }
  }

  return n;
}

/* Appends one shell's expanded, normalized basis functions into a growable
 * basis_function_t* array. Returns 1 on success, 0 on failure (out params left
 * unmodified on failure other than partial growth).
 */
static int append_shell_functions(const basis_shell_t *sh,
                                  const double center[3],
                                  basis_function_t ***arr, int *count,
                                  int *cap) {
  int comps[15][3]; // enough for up to l=4 (15 components)

  if (strcmp(sh->type, "SP") == 0) {
    // S part (l=0) from coef1, P part (l=1) from coef2
    basis_function_t *s_fn =
        basis_function_alloc(0, 0, 0, center, sh->n_prim, sh->exps, sh->coef1);
    if (!s_fn) {
      return 0;
    }

    molint_normalize_contraction(s_fn);

    if (*count == *cap) {
      int new_cap = *cap == 0 ? 8 : *cap * 2;

      basis_function_t **tmp = realloc(*arr, (size_t)new_cap * sizeof(*tmp));
      if (!tmp) {
        basis_function_free(s_fn);
        return 0;
      }

      *arr = tmp;
      *cap = new_cap;
    }

    (*arr)[(*count)++] = s_fn;

    int nc = cartesian_components(1, comps);
    for (int i = 0; i < nc; i++) {
      basis_function_t *p_fn =
          basis_function_alloc(comps[i][0], comps[i][1], comps[i][2], center,
                               sh->n_prim, sh->exps, sh->coef2);
      if (!p_fn) {
        return 0;
      }

      molint_normalize_contraction(p_fn);
      if (*count == *cap) {
        int new_cap = *cap == 0 ? 8 : *cap * 2;

        basis_function_t **tmp = realloc(*arr, (size_t)new_cap * sizeof(*tmp));
        if (!tmp) {
          basis_function_free(p_fn);
          return 0;
        }

        *arr = tmp;
        *cap = new_cap;
      }

      (*arr)[(*count)++] = p_fn;
    }

    return 1;
  }

  int nc = cartesian_components(sh->l, comps);
  for (int i = 0; i < nc; i++) {
    basis_function_t *fn =
        basis_function_alloc(comps[i][0], comps[i][1], comps[i][2], center,
                             sh->n_prim, sh->exps, sh->coef1);
    if (!fn) {
      return 0;
    }

    molint_normalize_contraction(fn);
    if (*count == *cap) {
      int new_cap = *cap == 0 ? 8 : *cap * 2;

      basis_function_t **tmp = realloc(*arr, (size_t)new_cap * sizeof(*tmp));
      if (!tmp) {
        basis_function_free(fn);
        return 0;
      }

      *arr = tmp;
      *cap = new_cap;
    }

    (*arr)[(*count)++] = fn;
  }

  return 1;
}
int basis_set_build_atom(const basis_element_t *elem, const double center[3],
                         basis_function_t ***out) {
  if (!elem || !out) {
    if (out) {
      *out = NULL;
    }

    return 0;
  }

  basis_function_t **arr = NULL;
  int count = 0;
  int cap = 0;

  for (int shell_idx = 0; shell_idx < elem->n_shells; shell_idx++) {
    if (!append_shell_functions(&elem->shells[shell_idx], center, &arr, &count,
                                &cap)) {
      basis_set_free_functions(arr, count);
      *out = NULL;

      return 0;
    }
  }

  *out = arr;

  return count;
}

static int append_atom_functions(basis_function_t ***all_fns, int *total,
                                 int *cap, basis_function_t **atom_fns,
                                 int n_atom_fns) {
  for (int i = 0; i < n_atom_fns; i++) {
    if (*total == *cap) {
      int new_cap = *cap == 0 ? INITIAL_CAP : *cap * 2;
      void *realloc_ptr = realloc((void *)*all_fns,
                                  (size_t)new_cap * sizeof(basis_function_t *));
      if (!realloc_ptr) {
        return 0;
      }

      *all_fns = (basis_function_t **)realloc_ptr;
      *cap = new_cap;
    }

    (*all_fns)[(*total)++] = atom_fns[i];
  }
  return 1;
}

int basis_set_build_molecule(const basis_set_t *basis_set,
                             const char *const symbols[],
                             const double centers[][3], int n_atoms,
                             basis_function_t ***out) {
  if (!basis_set || !symbols || !centers || n_atoms <= 0 || !out) {
    if (out) {
      *out = NULL;
    }

    return 0;
  }

  basis_function_t **all = NULL;
  int total = 0;
  int cap = 0;

  for (int atom_idx = 0; atom_idx < n_atoms; atom_idx++) {
    const basis_element_t *elem =
        basis_set_find_element(basis_set, symbols[atom_idx]);
    if (!elem) {
      basis_set_free_functions(all, total);
      *out = NULL;

      return 0;
    }

    basis_function_t **atom_fns = NULL;
    int n_atom_fns = basis_set_build_atom(elem, centers[atom_idx], &atom_fns);

    if (n_atom_fns == 0 && elem->n_shells > 0 && !atom_fns) {
      /* NOTE: Distinguish "genuinely zero shells" (shouldn't happen for a valid
       * element block) from allocation failure inside basis_set_build_atom:
       * atom_fns is NULL on failure. */
      basis_set_free_functions(all, total);
      *out = NULL;

      return 0;
    }

    if (!append_atom_functions(&all, &total, &cap, atom_fns, n_atom_fns)) {
      basis_set_free_functions(atom_fns, n_atom_fns);
      basis_set_free_functions(all, total);
      *out = NULL;

      return 0;
    }

    free((void *)atom_fns);
  }

  *out = all;

  return total;
}
