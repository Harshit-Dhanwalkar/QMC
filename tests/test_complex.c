#include "../core/complex.h"
#include <stdio.h>

int main() {
  printf(" > Testing Complex...\n");

  complex_t a = c_new(1.0, 2.0);
  complex_t b = c_new(3.0, 4.0);
  complex_t c = c_add(a, b);
  if (c.re == 4.0 && c.im == 6.0) {
    printf("   complex test passed.\n");

    return 0;
  }

  printf("   complex test FAILED.\n");

  return 1;
}
