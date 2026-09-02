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

// A bare-Tree probe feeds its index_t buffer straight into the misc.a suffstat
// kernels, which take misc_index_t*. If the two widths ever diverge, a probe
// that declares its indices as the wider type silently mis-strides the array
// and returns garbage (observed once with a size_t index buffer on small n).
// Guard the invariant here. Related gotcha for any new bare-kernel probe: call
// misc_simd_init() first (main.cpp does at startup) - misc_partition* are null
// function pointers until then.
static_assert(sizeof(index_t) == sizeof(misc_index_t),
              "tests/cpp index buffers feed the misc.a kernels; index_t and "
              "misc_index_t must be the same width");

// A storage-aware snapshot of a store's training codes: codeAt over every
// cell, laid out column-major over ALL numPredictors columns, which is what
// an index of j * numObservations + i means. train.codes packs the DENSELY
// stored columns only, and codeOffsets is a running cursor over those, so a
// snapshot taken off it cannot see a rank-stored column change and its
// per-column offsets are not j * numObservations on a mixed store.
inline std::vector<xint_t> storageDigest(const ColumnStore& data) {
  std::vector<xint_t> digest(data.numObservations * data.numPredictors);
  for (std::size_t j = 0; j < data.numPredictors; ++j)
    for (std::size_t i = 0; i < data.numObservations; ++i)
      digest[j * data.numObservations + i] =
        static_cast<xint_t>(data.codeAt(j, i));
  return digest;
}

// A canonical fingerprint of a tree's live structure alone (which nodes are
// split and on what), so a chain that MOVES can be told from one that only
// redraws its leaf parameters.
inline std::uint64_t treeStructureSignature(const Tree& tree) {
  std::vector<std::int32_t> subtree;
  tree.fillSubtree(0, subtree);
  std::uint64_t hash = 1469598103934665603ull;
  auto mix = [&hash](std::uint64_t value) {
    hash = (hash ^ value) * 1099511628211ull;
  };
  for (std::int32_t i : subtree) {
    const Node& node(tree.at(i));
    mix(node.isBottom() ? 0ull : 1ull);
    if (node.isBottom()) continue;
    mix(static_cast<std::uint64_t>(node.rule.variableIndex));
    mix(static_cast<std::uint64_t>(node.rule.splitIndex()));
  }
  return hash;
}

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

// A mixed dense + CSC predictor view for the engine builders: the fields a
// container-shaped fixture fills, gathered in one call.
inline PredictorSource mixedPredictorSource(
    size_t numRows, size_t numColumns, const double* denseValues,
    const int* pointers, const int* rows, const double* values,
    const std::int32_t* columnSources,
    const ColumnKind* columnTypes = nullptr,
    const std::uint32_t* categoryCounts = nullptr,
    const xint_t* referenceCodes = nullptr) {
  PredictorSource source;
  source.numRows = numRows;
  source.numColumns = numColumns;
  source.denseValues = denseValues;
  source.cscColumnPointers = pointers;
  source.cscRowIndices = rows;
  source.cscValues = values;
  source.columnSources = columnSources;
  source.columnTypes = columnTypes;
  source.categoryCounts = categoryCounts;
  source.referenceCodes = referenceCodes;
  return source;
}

// A logical matrix held both densely and as CSC arrays, for comparing the
// two build paths over identical values.
struct CscFixture {
  size_t n = 0, p = 0;
  std::vector<double> dense;   // column-major, zeros where nothing stored
  std::vector<int> pointers;   // p + 1
  std::vector<int> rows;
  std::vector<double> values;
  // the all-CSC column map (column j is CSC column j, the engine's ~j): the
  // spelling a bare sparse design takes through the one predictor view
  std::vector<std::int32_t> allCscSources;

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
    allCscSources.resize(p);
    for (size_t j = 0; j < p; ++j)
      allCscSources[j] = ~static_cast<std::int32_t>(j);
  }
};

// ext_printf is Rprintf (external/io.h), whose real implementation needs a
// live R session and segfaults without one, so this host defines the symbol
// itself and the executable's definition binds ahead of the framework's.
// Output is discarded unless a capture is armed, which is how the engine's
// info dumps become assertable here.
void beginPrintCapture(std::string& sink);
void endPrintCapture();

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
void runShapeTests(ext_rng* rng);
// no rng argument on purpose: the conformance fixtures own their generators
// and the suite restores the shared runif01 stream, so it neither shifts nor
// is shifted by any other suite's draws
void runFacadeTests();
void runStateTests(ext_rng* rng);
// no rng argument on purpose: the ensemble oracle owns its generator and
// restores the shared runif01 stream, so it neither shifts nor is shifted by
// any other suite's draws
void runEnsembleTests();
void runFuzzTests(int numSeeds);

#endif  // TESTS_CPP_COMMON_HPP
