/*
 * Test: module catalog and qmc_help()
 *
 * 1. Every entry is complete (name, header, category, summary), names are
 *    unique, and the categories listed are exactly those used by modules
 * 2. qmc_module_find: exact, header-path, case-insensitive and unknown names
 * 3. qmc_module_search: keyword hits, total count beyond max_out, empty query
 * 4. qmc_help output: overview, module detail (#include line, start-with list),
 *    category listing, keyword search and the no-match message
 * 5. Modules outside <qmc/qmc.h> (HDF5, NetCDF, plot) say so in their help
 */
#include <stdio.h>
#include <string.h>

#include "../core/catalog.h"
#include "../core/qmc_version.h"

static int failures = 0;

static void check(int cond, const char *what) {
  if (!cond) {
    fprintf(stderr, "FAIL: %s\n", what);
    failures++;
  }
}

/* Runs qmc_help into a buffer; returns the number of modules shown. */
static size_t help_to_buf(const char *topic, char *buf, size_t cap) {
  FILE *f = tmpfile();
  if (f == NULL) {
    check(0, "tmpfile available");
    buf[0] = '\0';
    return 0;
  }
  size_t shown = qmc_help(f, topic);
  rewind(f);
  size_t n = fread(buf, 1, cap - 1, f);
  buf[n] = '\0';
  fclose(f);
  return shown;
}

int main(void) {
  static char buf[32768];
  size_t n = qmc_module_count();
  check(n > 100, "catalog lists more than 100 modules");

  for (size_t i = 0; i < n; i++) {
    const qmc_module_t *m = qmc_module_at(i);
    check(m != NULL && m->name[0] != '\0' && m->header[0] != '\0' &&
              m->category[0] != '\0' && m->summary[0] != '\0',
          "every module has name, header, category and summary");
    check(strstr(m->header, ".h") != NULL, "header path ends in .h");
    int cat_known = 0;
    for (size_t c = 0; c < qmc_category_count(); c++) {
      if (strcmp(qmc_category_at(c), m->category) == 0) {
        cat_known = 1;
      }
    }
    check(cat_known, "module category is in the category list");
    for (size_t j = i + 1; j < n; j++) {
      check(strcmp(m->name, qmc_module_at(j)->name) != 0, "module names are unique");
    }
  }
  check(qmc_module_at(n) == NULL, "qmc_module_at past the end is NULL");
  check(qmc_category_at(qmc_category_count()) == NULL, "qmc_category_at past the end is NULL");

  /* find */
  const qmc_module_t *tb = qmc_module_find("tight_binding");
  check(tb != NULL && strcmp(tb->header, "physics/tight_binding.h") == 0,
        "find by name");
  check(qmc_module_find("physics/tight_binding.h") == tb, "find by header path");
  check(qmc_module_find("qmc/physics/tight_binding.h") == tb, "find by installed path");
  check(qmc_module_find("tight_binding.h") == tb, "find by file name");
  check(qmc_module_find("Tight_Binding") == tb, "find is case-insensitive");
  check(qmc_module_find("no_such_module") == NULL, "unknown name is NULL");
  check(qmc_module_find(NULL) == NULL, "NULL name is NULL");
  check(qmc_module_find("") == NULL, "empty name is NULL");

  /* search */
  const qmc_module_t *hits[4];
  size_t total = qmc_module_search("chern", hits, 4);
  int found_tb = 0;
  for (size_t i = 0; i < total && i < 4; i++) {
    found_tb |= hits[i] == tb;
  }
  check(total >= 1 && found_tb, "search 'chern' finds tight_binding");
  check(qmc_module_search("e", hits, 2) > 2, "search returns total beyond max_out");
  check(qmc_module_search("e", NULL, 0) > 2, "search with no output buffer still counts");
  check(qmc_module_search("", hits, 4) == 0, "empty query matches nothing");
  check(qmc_module_search(NULL, hits, 4) == 0, "NULL query matches nothing");
  check(qmc_module_search("zzzzqqqq", hits, 4) == 0, "nonsense query matches nothing");

  /* help */
  size_t shown = help_to_buf(NULL, buf, sizeof buf);
  check(shown == n, "overview shows every module");
  check(strstr(buf, "QMC " QMC_VERSION_STRING) != NULL, "overview names the version");
  check(strstr(buf, "tight_binding") != NULL && strstr(buf, "hydrogen") != NULL,
        "overview lists modules");

  shown = help_to_buf("tight_binding", buf, sizeof buf);
  check(shown == 1, "module help shows one module");
  check(strstr(buf, "#include <qmc/physics/tight_binding.h>") != NULL,
        "module help gives the #include line");
  check(strstr(buf, "Start with:") != NULL && strstr(buf, "tb_model_alloc") != NULL,
        "module help lists where to start");
  check(strstr(buf, "/api/tight_binding.html") != NULL, "module help points at the API page");

  shown = help_to_buf("chemistry", buf, sizeof buf);
  check(shown >= 5 && strstr(buf, "molecular_hf") != NULL && strstr(buf, "ccsd") != NULL,
        "category help lists its modules");

  shown = help_to_buf("chern", buf, sizeof buf);
  check(shown >= 1 && strstr(buf, "tight_binding") != NULL, "keyword help finds modules");

  shown = help_to_buf("zzzzqqqq", buf, sizeof buf);
  check(shown == 0 && strstr(buf, "No module") != NULL, "no-match message");

  help_to_buf("hdf5_writer", buf, sizeof buf);
  check(strstr(buf, "not in <qmc/qmc.h>") != NULL, "optional-dependency modules say so");
  help_to_buf("hydrogen", buf, sizeof buf);
  check(strstr(buf, "not in <qmc/qmc.h>") == NULL, "umbrella modules do not");

  if (failures != 0) {
    fprintf(stderr, "test_catalog: %d check(s) failed\n", failures);
    return 1;
  }
  printf("All test_catalog checks passed (%zu modules, %zu categories).\n", n,
         qmc_category_count());
  return 0;
}
