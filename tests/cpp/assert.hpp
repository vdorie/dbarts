#ifndef TESTS_CPP_ASSERT_HPP
#define TESTS_CPP_ASSERT_HPP

// Zero-dependency micro-assert helpers, usable without pulling in any
// bartcore engine header (kept separate from common.hpp so TUs that only
// touch the low layers, e.g. data/tree, do not depend on the full stack).

#include <cstdint>
#include <cstdio>

extern int failures;

void check(bool condition, const char* what);
void checkNear(double actual, double expected, double tolerance,
               const char* what);

extern uint64_t rngState;
double runif01();

#endif  // TESTS_CPP_ASSERT_HPP
