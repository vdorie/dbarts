#ifndef TESTS_CPP_COMMON_HPP
#define TESTS_CPP_COMMON_HPP

// Shared includes and fixtures for the component tests that need the full
// engine stack (data/tree/model/moves/chain/sampler/facade); test_data.cpp
// and test_tree.cpp stay off this header on purpose, so a touch to e.g.
// chain.hpp or sampler.hpp does not force them to recompile.

#include "assert.hpp"
#include "test_peer.hpp"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
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

// Leaves of `tree` that hold rows and none of positive weight: legal under the
// membership rule, and what a weight-counting rule would never let a move or a
// prior draw produce.
inline std::size_t countWeightlessLeaves(const Tree& tree,
                                         const double* weights) {
  std::vector<std::int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  std::size_t weightless = 0;
  for (std::int32_t b : bottoms) {
    const Node& node(tree.at(b));
    bool anyWeight = false;
    for (std::size_t j = node.begin; j < node.end; ++j)
      anyWeight = anyWeight || weights[tree.indices[j]] > 0.0;
    weightless += node.numObservations() > 0 && !anyWeight ? 1 : 0;
  }
  return weightless;
}

// A generator at the same position as `rng`, for a reference arm.
inline ext_rng* cloneRng(const ext_rng* rng) {
  ext_rng* copy = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  std::vector<unsigned char> state(ext_rng_getSerializedStateLength(rng));
  ext_rng_writeSerializedState(rng, state.data());
  ext_rng_readSerializedState(copy, state.data());
  return copy;
}

// Advance both generators past the kernel under test and compare the streams
// they leave behind: identical tails prove identical consumption, which is
// what a skipped (rather than drawn-and-discarded) latent buys.
inline bool rngStreamsAgree(ext_rng* a, ext_rng* b, int numDraws = 32) {
  for (int j = 0; j < numDraws; ++j)
    if (ext_rng_simulateContinuousUniform(a) !=
        ext_rng_simulateContinuousUniform(b))
      return false;
  return true;
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

// Whether two reads of a sampler's state are one sampler: its chains,
// generators, latents and kept draws (statesAgree), its cut grid, its
// store's write position and the columns that can hold a missing value.
static inline bool samplerStatesAgree(const SamplerStateData& a,
                                      const SamplerStateData& b) {
  return statesAgree(a, b) && a.cutPoints == b.cutPoints &&
         a.cutMass == b.cutMass && a.currentSampleNum == b.currentSampleNum &&
         a.recordedDraws == b.recordedDraws &&
         a.missingColumns == b.missingColumns;
}
// Installs a state without force and then with it, and returns whether the
// two answers are one: the unforced call takes the state exactly when the
// forced install reports nothing altered, and then leaves the sampler the
// forced install does; it declines exactly when the forced install reports a
// repair; and a state one form refuses the other refuses. A call that does
// not install leaves the sampler as it was. Every flag is preset to the
// answer that would hide a call never writing it. installed and altered
// report the forced call.
template <typename S>
static bool verdictMatchesInstall(S& sampler, const SamplerStateData& state,
                                  const double* currentPredictors,
                                  bool& installed, bool& altered) {
  SamplerStateData before, unforced, forced;
  sampler.getState(before);
  bool notClean = false;
  bool took = sampler.setState(state, currentPredictors, nullptr, nullptr,
                               keepStoreCapacity, false, &notClean);
  sampler.getState(unforced);
  bool untouched = took || samplerStatesAgree(before, unforced);
  altered = !notClean;
  installed = sampler.setState(state, currentPredictors, nullptr, &altered);
  sampler.getState(forced);
  if (!installed)
    return !took && !notClean && !altered &&
           samplerStatesAgree(before, forced);
  if (took)
    return !notClean && !altered && samplerStatesAgree(unforced, forced);
  return untouched && notClean && altered;
}
// Installs a state both ways (verdictMatchesInstall) and holds the answer to
// the one named: installed as stored, or declined without force and repaired
// with it.
template <typename S>
static bool restoresWithStatus(S& sampler, const SamplerStateData& state,
                               const double* currentPredictors, bool expected) {
  bool installed = false, altered = !expected;
  return verdictMatchesInstall(sampler, state, currentPredictors, installed,
                               altered) &&
         installed && altered == expected;
}
template <typename S>
static bool restoresExactly(S& sampler, const SamplerStateData& state,
                            const double* currentPredictors = nullptr) {
  return restoresWithStatus(sampler, state, currentPredictors, false);
}
template <typename S>
static bool restoresAltered(S& sampler, const SamplerStateData& state,
                            const double* currentPredictors = nullptr) {
  return restoresWithStatus(sampler, state, currentPredictors, true);
}
// Whether an install without force declines a state as not clean and leaves
// the sampler as it was. The flag is preset down, so a refusal fails.
template <typename S>
static bool declinesUntouched(S& sampler, const SamplerStateData& state,
                              const double* currentPredictors = nullptr) {
  SamplerStateData before, after;
  sampler.getState(before);
  bool notClean = false;
  bool installed =
    sampler.setState(state, currentPredictors, nullptr, nullptr,
                     keepStoreCapacity, false, &notClean);
  sampler.getState(after);
  return !installed && notClean && samplerStatesAgree(before, after);
}
// Whether both forms of an install refuse a state, without a verdict, and
// leave the sampler as it was. The flag is preset up, so a decline fails.
template <typename S>
static bool refusesUntouched(S& sampler, const SamplerStateData& state,
                             const double* currentPredictors = nullptr) {
  SamplerStateData before, after;
  sampler.getState(before);
  bool notClean = true;
  bool refused =
    !sampler.setState(state, currentPredictors) &&
    !sampler.setState(state, currentPredictors, nullptr, nullptr,
                      keepStoreCapacity, false, &notClean) &&
    !notClean;
  sampler.getState(after);
  return refused && samplerStatesAgree(before, after);
}
// Whether an install without force takes a state, as one that needs no
// repair; the flag is preset up, so an install that never writes it fails.
template <typename S>
static bool installsClean(S& sampler, const SamplerStateData& state,
                          const double* currentPredictors = nullptr) {
  bool notClean = true;
  return sampler.setState(state, currentPredictors, nullptr, nullptr,
                          keepStoreCapacity, false, &notClean) &&
         !notClean;
}
// Flags the first ordinal rule among flat trees as sending missing values
// right; false when none splits on an ordinal column.
bool sendFirstOrdinalRuleMissingRight(std::vector<std::vector<FlatNode>>& trees);

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

/// Chain c's forest f totals through the sampler's own read.
template <typename S>
std::vector<double> forestTotals(const S& sampler, std::size_t c,
                                 std::size_t f = 0) {
  std::vector<double> out(sampler.numObservations());
  sampler.forestTotalFits(c, f, out.data());
  return out;
}

// The forest cache rule's measure. A forest's cached totalFits may differ from
// its tree fits summed in tree order by accumulated ADDITIVE rounding only: the
// sweep keeps the cache by difference updates, each rounding at the scale of
// the forest's working response at that sweep, so the gap is bounded by
// C eps max_s max_i |forestY_i(s)| sqrt(sweeps). Returns forest f's worst gap
// in units of eps scaleMax sqrt(sweeps), where scaleMax is the caller's
// running maximum of max_i |forestY_i| over the sweeps it has seen, updated
// here first: the gap keeps the rounding of a sweep whose multiplier was small
// and whose response was large after the multiplier grows again. forestY is
// recovered from the last sweep's residual as treeY + totalFits - the last
// tree's fits, which is the response those updates rounded against, and is
// floored at 1 so a forest with a tiny response is still held to an absolute
// floor.
template <typename C>
double forestCacheGapRatio(const C& chain, std::size_t f, std::size_t sweeps,
                           double& scaleMax) {
  const std::vector<double>& total = TestPeer::totalFitsInForest(chain, f);
  std::size_t n = total.size(), numTrees = chain.numTreesInForest(f);
  if (numTrees == 0) return 0.0;
  std::vector<double> fits(n * numTrees);
  TestPeer::forestTreeFits(chain, f, fits.data());
  const auto& resid = TestPeer::residual(chain, f);
  const double* last = fits.data() + (numTrees - 1) * n;
  double gap = 0.0;
  scaleMax = std::max(scaleMax, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    double gather = 0.0;
    for (std::size_t t = 0; t < numTrees; ++t) gather += fits[t * n + i];
    gap = std::max(gap, std::fabs(total[i] - gather));
    if (i < resid.size())
      scaleMax =
        std::max(scaleMax, std::fabs(static_cast<double>(resid[i]) + total[i] -
                                     last[i]));
  }
  double root = std::sqrt(static_cast<double>(sweeps > 0 ? sweeps : 1));
  return gap / (std::numeric_limits<double>::epsilon() * scaleMax * root);
}

// C in that bound, shared by the drift pin (testAmplitudeCacheDrift) and the
// fuzz invariant. Measured with no multiplicative leaf transform in the
// engine: the pin's worst ratio was 5.0 (probit, the 50 + 25 tree ensemble)
// and the fuzz's 5.0 over 1000 seeds; C is ten times that. The fuzz once read
// 184, when its running scale was sampled only after each op and so missed
// the working response of all but the last sweep of a multi-sweep
// grow-from-root; its grow op now runs a sweep at a time. The removed
// rescaling move multiplied the gap and reached 4.7e10 (probit) and 7.7e11
// (logistic) at the latent gate's shape within 5000 sweeps.
constexpr double forestCacheGapBound = 50.0;

// ext_printf is Rprintf (external/io.h), whose real implementation needs a
// live R session and segfaults without one, so this host defines the symbol
// itself and the executable's definition binds ahead of the framework's.
// Output is discarded unless a capture is armed, which is how the engine's
// info dumps become assertable here.
void beginPrintCapture(std::string& sink);
void endPrintCapture();

// Distribution checks shared by the suites (common.cpp): probabilities from
// log weights, Pearson's statistic against a fully specified distribution of
// probabilities over counts, and its chi-square upper tail.
std::vector<double> normalizedFromLogWeights(
  const std::vector<double>& logWeights);
double chiSquareStatistic(const std::vector<double>& counts,
                          const std::vector<double>& probabilities,
                          double numDraws);
double chiSquareUpperTail(double statistic, double df);

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
// the flat C API's per-draw struct and its registration; takes no rng
void runCapiTests();
void runStateTests(ext_rng* rng);
// no rng argument on purpose: the ensemble oracle owns its generator and
// restores the shared runif01 stream, so it neither shifts nor is shifted by
// any other suite's draws
void runEnsembleTests();
// The monotone leaf-order counter; own rng, restores the runif01 stream.
void runMonotoneTests();
void runFuzzTests(int numSeeds);

#endif  // TESTS_CPP_COMMON_HPP
