// Factorials and common utilities
#include "special.h"
#include <math.h>

double factorial(int n) {
  if (n < 0)
    return 0.0;
  double f = 1.0;
  for (int i = 2; i <= n; i++)
    f *= i;
  return f;
}

double log_factorial(int n) {
  if (n < 0)
    return -INFINITY;
  return lgamma(n + 1);
}
