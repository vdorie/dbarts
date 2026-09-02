#ifndef TESTS_CPP_ASSERT_HPP
#define TESTS_CPP_ASSERT_HPP

// Zero-dependency micro-assert helpers, usable without pulling in any
// bartcore engine header (kept separate from common.hpp so TUs that only
// touch the low layers, e.g. data/tree, do not depend on the full stack).

#include <cstdint>
#include <cstdio>

extern int failures;

void check(bool condition, const char* what);
// A fixture build the test expects to succeed. ColumnStore::build and
// buildTest are [[nodiscard]] so no product caller can drop a refusal; a
// fixture over valid data still says what it expects, since a store that
// stopped building half way would otherwise fail somewhere unrelated.
void built(bool ok);
void checkNear(double actual, double expected, double tolerance,
               const char* what);

extern uint64_t rngState;
double runif01();

#endif  // TESTS_CPP_ASSERT_HPP
