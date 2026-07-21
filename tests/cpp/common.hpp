#ifndef TESTS_CPP_COMMON_HPP
#define TESTS_CPP_COMMON_HPP

// Shared includes and fixtures for the component tests that need the full
// engine stack (data/tree/model/moves/chain/sampler/facade); test_data.cpp
// and test_tree.cpp stay off this header on purpose, so a touch to e.g.
// chain.hpp or sampler.hpp does not force them to recompile.

#include "assert.hpp"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <misc/partition.h>
#include <misc/simd.h>
#include <external/random.h>

#include <bartcore/bartcore.hpp>

using namespace bartcore;

// Structural round-trip gate for the state tests. With bitwise continuation
// dropped, a restored sampler must reconstruct the model - trees, leaf
// parameters, saved trees, latents, dart, rng - exactly, sigma to within the
// original-scale round trip. A flat split's payload (cut point or mask word)
// compares as its raw word; a leaf's compares to the last ulp, since a
// function leaf's value is a reporting mean whose sum order the canonical
// rebuild does not preserve.
bool sameFlatTrees(const std::vector<std::vector<FlatNode>>& a,
                   const std::vector<std::vector<FlatNode>>& b);
bool statesAgree(const SamplerStateData& a, const SamplerStateData& b);

// Gate (a): a state re-captured from the restored sampler reproduces the
// saved one, so restore reconstructs the model, not the accumulation history.
template <typename S>
static void checkStructuralRoundTrip(const SamplerStateData& saved,
                                     S& restored, const char* label) {
  SamplerStateData reState;
  restored.getState(reState);
  check(statesAgree(saved, reState), label);
}

// A burned-in sampler for mutation tests: strong signal in both columns so
// trees certainly split.
std::unique_ptr<ConstantLeafSampler> makeBurnedInSampler(
  std::vector<double>& x, std::vector<double>& y, size_t n, ext_rng* rng);
void makeMutationData(std::vector<double>& x, std::vector<double>& y,
                       size_t n);

// A logical matrix held both densely and as CSC arrays, for comparing the
// two build paths over identical values.
struct CscFixture {
  size_t n = 0, p = 0;
  std::vector<double> dense;   // column-major, zeros where nothing stored
  std::vector<int> pointers;   // p + 1
  std::vector<int> rows;
  std::vector<double> values;

  // fraction of rows stored per column; stored NaNs count as entries, the
  // Matrix convention for missing values
  void build(size_t n_, const std::vector<double>& nonzeroFractions,
             size_t numMissingPerColumn = 0) {
    n = n_;
    p = nonzeroFractions.size();
    dense.assign(n * p, 0.0);
    pointers.assign(p + 1, 0);
    rows.clear();
    values.clear();
    for (size_t j = 0; j < p; ++j) {
      size_t numMissing = 0;
      for (size_t i = 0; i < n; ++i) {
        if (runif01() >= nonzeroFractions[j]) continue;
        double value = 0.5 + runif01();
        if (numMissing < numMissingPerColumn) {
          value = std::nan("");
          ++numMissing;
        }
        dense[i + j * n] = value;
        rows.push_back(static_cast<int>(i));
        values.push_back(value);
      }
      pointers[j + 1] = static_cast<int>(rows.size());
    }
  }
};

// One entry point per translation unit, called from main() in original
// test order; each runs its area's tests, filtered by suite name there.
void runDataTests();
void runTreeTests(ext_rng* rng);
void runScanTests();
void runGrowTests(ext_rng* rng);
void runMovesTests(ext_rng* rng);
void runInteractionTests(ext_rng* rng);
void runModelTests(ext_rng* rng);
void runSamplerTests(ext_rng* rng);
void runStateTests(ext_rng* rng);
void runFuzzTests(int numSeeds);

#endif  // TESTS_CPP_COMMON_HPP
