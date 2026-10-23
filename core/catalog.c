#include "catalog.h"

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "qmc_version.h"

#include "catalog_data.inc"

#define HELP_WIDTH 78

static size_t module_total(void) {
  return sizeof qmc_catalog_modules / sizeof qmc_catalog_modules[0];
}

static size_t category_total(void) {
  return sizeof qmc_catalog_categories / sizeof qmc_catalog_categories[0];
}

static int ieq(const char *a, const char *b) {
  while (*a != '\0' && *b != '\0') {
    if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) {
      return 0;
    }

    a++;
    b++;
  }

  return *a == *b;
}

static int icontains(const char *hay, const char *needle) {
  size_t n = strlen(needle);
  if (n == 0) {
    return 1;
  }

  for (; *hay != '\0'; hay++) {
    size_t i = 0;
    while (i < n && hay[i] != '\0' &&
           tolower((unsigned char)hay[i]) ==
               tolower((unsigned char)needle[i])) {
      i++;
    }

    if (i == n) {
      return 1;
    }
  }

  return 0;
}

size_t qmc_module_count(void) { return module_total(); }

const qmc_module_t *qmc_module_at(size_t index) {
  return index < module_total() ? &qmc_catalog_modules[index] : NULL;
}

size_t qmc_category_count(void) { return category_total(); }

const char *qmc_category_at(size_t index) {
  return index < category_total() ? qmc_catalog_categories[index] : NULL;
}

const qmc_module_t *qmc_module_find(const char *name) {
  if (name == NULL) {
    return NULL;
  }

  const char *base = strrchr(name, '/');
  base = base != NULL ? base + 1 : name;
  size_t len = strlen(base);
  if (len > 2 && strcmp(base + len - 2, ".h") == 0) {
    len -= 2;
  }

  for (size_t i = 0; i < module_total(); i++) {
    const char *mn = qmc_catalog_modules[i].name;
    if (strlen(mn) == len && strncmp(mn, base, len) == 0) {
      return &qmc_catalog_modules[i];
    }
  }

  // case-insensitive second pass
  for (size_t i = 0; i < module_total(); i++) {
    const char *mn = qmc_catalog_modules[i].name;
    if (strlen(mn) == len) {
      char buf[64];
      if (len < sizeof buf) {
        memcpy(buf, base, len);
        buf[len] = '\0';
        if (ieq(mn, buf)) {
          return &qmc_catalog_modules[i];
        }
      }
    }
  }

  return NULL;
}

size_t qmc_module_search(const char *query, const qmc_module_t **out,
                         size_t max_out) {
  if (query == NULL || query[0] == '\0') {
    return 0;
  }

  size_t total = 0;
  for (size_t i = 0; i < module_total(); i++) {
    const qmc_module_t *m = &qmc_catalog_modules[i];
    if (icontains(m->name, query) || icontains(m->summary, query) ||
        icontains(m->category, query) || icontains(m->key_functions, query)) {
      if (out != NULL && total < max_out) {
        out[total] = m;
      }

      total++;
    }
  }

  return total;
}

// Help

// Prints text word-wrapped to HELP_WIDTH, continuation lines indented
static void print_wrapped(FILE *out, const char *first_prefix,
                          const char *indent, const char *text) {
  size_t col = strlen(first_prefix);
  size_t ind = strlen(indent);
  fputs(first_prefix, out);
  int at_line_start = 1;
  const char *p = text;
  while (*p != '\0') {
    while (*p == ' ') {
      p++;
    }

    const char *w = p;
    while (*p != '\0' && *p != ' ') {
      p++;
    }

    size_t wl = (size_t)(p - w);
    if (wl == 0) {
      break;
    }

    if (!at_line_start && col + 1 + wl > HELP_WIDTH) {
      fputc('\n', out);
      fputs(indent, out);
      col = ind;
      at_line_start = 1;
    }

    if (!at_line_start) {
      fputc(' ', out);
      col++;
    }

    fwrite(w, 1, wl, out);
    col += wl;
    at_line_start = 0;
  }

  fputc('\n', out);
}

static void print_module_detail(FILE *out, const qmc_module_t *m) {
  fprintf(out, "%s  [%s]\n", m->name, m->category);
  print_wrapped(out, "  ", "  ", m->summary);

  fprintf(out, "\n  #include <qmc/%s>", m->header);
  if (!m->in_qmc_h) {
    fputs("   (not in <qmc/qmc.h>: include it directly)", out);
  } else {
    fputs("   (or <qmc/qmc.h> for everything)", out);
  }

  fputc('\n', out);

  if (m->key_functions[0] != '\0') {
    print_wrapped(out, "  Start with: ", "              ", m->key_functions);
  }
  if (m->guide[0] != '\0') {
    fprintf(out, "  Guide:      %s/%.*s.html\n", QMC_DOCS_URL,
            (int)(strlen(m->guide) - 3), m->guide);
  }
  if (m->example[0] != '\0') {
    fprintf(out, "  Example:    %s/blob/main/%s\n", QMC_REPO_URL, m->example);
  }
  if (m->test[0] != '\0') {
    fprintf(
        out,
        "  Test:       %s/blob/main/%s  (every call, with expected results)\n",
        QMC_REPO_URL, m->test);
  }

  fprintf(out, "  API:        %s/api/%s.html  (all types and functions)\n",
          QMC_DOCS_URL, m->name);
  fprintf(out, "  Source:     %s/blob/main/%s\n", QMC_REPO_URL, m->header);
}

static void print_module_list(FILE *out, const qmc_module_t *const *mods,
                              size_t n) {
  for (size_t i = 0; i < n; i++) {
    char prefix[96];
    (void)snprintf(prefix, sizeof prefix, "  %-20s ", mods[i]->name);
    print_wrapped(out, prefix, "                       ", mods[i]->summary);
  }
}

static void print_overview(FILE *out) {
  fprintf(out, "QMC %s: quantum mechanics and quantum chemistry in C\n",
          QMC_VERSION_STRING);
  fprintf(
      out,
      "  %zu modules in %zu categories. #include <qmc/qmc.h> for everything,\n"
      "  or one header, e.g. #include <qmc/physics/hydrogen.h>.\n\n",
      module_total(), category_total());
  fputs("  qmc_help(stdout, \"<module>\")    how to use a module,  e.g. "
        "\"tight_binding\"\n"
        "  qmc_help(stdout, \"<category>\")  modules in a category, e.g. "
        "\"chemistry\"\n"
        "  qmc_help(stdout, \"<keyword>\")   search,                e.g. "
        "\"chern\", \"dmrg\"\n"
        "  from a shell:  qmc-help <topic>      full reference: " QMC_DOCS_URL
        "/api/\n",
        out);

  for (size_t c = 0; c < category_total(); c++) {
    fprintf(out, "\n%s\n", qmc_catalog_categories[c]);

    char line[HELP_WIDTH + 1];
    size_t len = 0;
    line[0] = '\0';
    for (size_t i = 0; i < module_total(); i++) {
      const qmc_module_t *m = &qmc_catalog_modules[i];
      if (strcmp(m->category, qmc_catalog_categories[c]) != 0) {
        continue;
      }

      size_t nl = strlen(m->name);
      if (len > 0 && 2 + len + 2 + nl > HELP_WIDTH) {
        fprintf(out, "  %s\n", line);
        len = 0;
        line[0] = '\0';
      }

      if (len > 0) {
        line[len++] = ' ';
        line[len++] = ' ';
      }

      memcpy(line + len, m->name, nl + 1);
      len += nl;
    }

    if (len > 0) {
      fprintf(out, "  %s\n", line);
    }
  }
}

size_t qmc_help(FILE *out, const char *topic) {
  if (out == NULL) {
    out = stdout;
  }

  if (topic == NULL || topic[0] == '\0') {
    print_overview(out);
    return module_total();
  }

  const qmc_module_t *m = qmc_module_find(topic);
  if (m != NULL) {
    print_module_detail(out, m);
    return 1;
  }

  // Category
  if (strlen(topic) >= 4) {
    size_t shown = 0;
    for (size_t c = 0; c < category_total(); c++) {
      if (!icontains(qmc_catalog_categories[c], topic)) {
        continue;
      }

      fprintf(out, "%s%s\n", shown > 0 ? "\n" : "", qmc_catalog_categories[c]);
      for (size_t i = 0; i < module_total(); i++) {
        if (strcmp(qmc_catalog_modules[i].category,
                   qmc_catalog_categories[c]) == 0) {
          const qmc_module_t *one = &qmc_catalog_modules[i];
          print_module_list(out, &one, 1);
          shown++;
        }
      }
    }

    if (shown > 0) {
      fputs("\nqmc_help(stdout, \"<module>\") shows how to use one of them.\n",
            out);
      return shown;
    }
  }

  // Keyword search
  enum { MAX_HITS = 16 };
  const qmc_module_t *hits[MAX_HITS];
  size_t total = qmc_module_search(topic, hits, MAX_HITS);
  if (total == 0) {
    fprintf(out,
            "No module, category or keyword matches \"%s\".\n"
            "Try qmc_help(stdout, NULL) for the full list.\n",
            topic);
    return 0;
  }

  size_t shown = 0;
  while (shown < total && shown < MAX_HITS) {
    shown++;
  }

  fprintf(out, "%zu module%s match \"%s\":\n", total, total == 1 ? "" : "s",
          topic);
  print_module_list(out, hits, shown);
  if (total > shown) {
    fprintf(out, "  ... and %zu more; use a longer keyword.\n", total - shown);
  }

  fputs("\nqmc_help(stdout, \"<module>\") shows how to use one of them.\n",
        out);

  return shown;
}
