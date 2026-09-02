#include "assert.hpp"

#include <cmath>
#include <cstdio>

int failures = 0;

void check(bool condition, const char* what) {
  if (!condition) {
    ++failures;
    printf("FAIL: %s\n", what);
  }
}

void built(bool ok) { check(ok, "the fixture store built"); }

void checkNear(double actual, double expected, double tolerance,
               const char* what) {
  if (!(std::fabs(actual - expected) <= tolerance)) {
    ++failures;
    printf("FAIL: %s (actual %.15g, expected %.15g)\n", what, actual, expected);
  }
}

uint64_t rngState = 0x9E3779B97F4A7C15ull;
double runif01() {
  rngState ^= rngState << 13;
  rngState ^= rngState >> 7;
  rngState ^= rngState << 17;
  return (double) (rngState >> 11) * 0x1.0p-53;
}
