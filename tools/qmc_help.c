/* qmc-help: command-line front end to the module catalog
 *
 *   qmc-help                   overview of every category and module
 *   qmc-help tight_binding     how to use a module
 *   qmc-help chemistry         the modules in a category
 *   qmc-help chern dmrg        keyword search (several topics allowed)
 *   qmc-help --list            module names, one per line (for scripts/completion)
 *   qmc-help --headers         installed header paths, one per line
 *   qmc-help --version
 */
#include <stdio.h>
#include <string.h>

#include "catalog.h"
#include "qmc_version.h"

int main(int argc, char **argv) {
  if (argc >= 2 && strcmp(argv[1], "--version") == 0) {
    printf("qmc %s\n", qmc_version());

    return 0;
  }

  if (argc >= 2 && (strcmp(argv[1], "--list") == 0 || strcmp(argv[1], "--headers") == 0)) {
    int headers = strcmp(argv[1], "--headers") == 0;
    for (size_t i = 0; i < qmc_module_count(); i++) {
      const qmc_module_t *m = qmc_module_at(i);
      puts(headers ? m->header : m->name);
    }

    return 0;
  }

  if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
    puts("usage: qmc-help [topic ...]   (module, category or keyword; none = overview)\n"
         "       qmc-help --list | --headers | --version");

    return 0;
  }

  if (argc == 1) {
    (void)qmc_help(stdout, NULL);

    return 0;
  }

  int status = 0;
  for (int i = 1; i < argc; i++) {
    if (i > 1) {
      putchar('\n');
    }

    if (qmc_help(stdout, argv[i]) == 0) {
      status = 1;
    }
  }

  return status;
}
