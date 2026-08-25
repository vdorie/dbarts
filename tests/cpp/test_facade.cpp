// Conformance gate for the type-erased boundary (bartcore/facade.hpp).
//
// SamplerBase is the one dispatch layer between the shipped flat C API and the
// engine: 59 pure virtuals, each forwarded by SamplerFacade<L> to a Sampler<L>
// method of the same name. Every other test in this suite drives Sampler<L>
// DIRECTLY, so a forwarder that drops an argument, inverts a flag, or reads the
// wrong forest/chain/slot is invisible to them - the engine twin it should have
// reached is tested, and the hop is not.
//
// This file calls every virtual THROUGH a SamplerBase reference and checks the
// answer against the concrete sampler's own state or against a direct call on
// its Sampler<L>. Three properties make a defective forwarder fail here:
//
//  - the fixtures make confusable quantities DIFFER. The saved-tree store is
//    partly filled (retained draws != capacity) and, where a row needs it,
//    wrapped (slot != draw index); samplers carry two chains, the combining
//    fixture two forests of unequal tree count and unequal amplitude width. A
//    forwarder that answers with the neighbouring quantity therefore answers
//    with a different value, not the same one.
//  - the table carries one row per virtual and the spy below records which
//    virtual each row reached, so a row wired to the wrong member fails on its
//    own coverage line rather than passing vacuously.
//  - the spy overrides EVERY pure virtual. A virtual added to SamplerBase makes
//    it abstract - a compile error - until it is named here, and the coverage
//    gate then fails until the table carries a row that calls it.
//
// The suite owns its generators and restores the shared runif01 stream, so it
// neither shifts nor is shifted by any other suite's draws.

#include "common.hpp"

namespace {

/// One enumerator per SamplerBase pure virtual, in declaration order.
enum class FacadeVirtual {
  shape, run, setOffset, setResponse, setWeights, weightsDigest,
  reapplyWeights, setSigma, setTestData, setTestOffset, setData, setPredictor,
  updatePredictor, setCutPoints, updatePredictorPerObservation,
  beginPredictorUpdate, currentSampleNum, savedSlotForDraw, savedTree,
  savedTreeSlopes, savedTreeMasks, flattenTree, predict, predictPerForest,
  predictVariance, getState, setState, installForests, sampleTreesFromPrior,
  sampleNodeParametersFromPrior, growFromRoot, setNumThreads, setNumThin,
  setVerbose, fitScale, setTreeStorage, setModel, sumOfSquaredResiduals,
  printTrees, rng, data, latents, sigma, dispersion, setForestBasis,
  setForestWeights, forestCalibration, setForestPriorScale, setActiveRows,
  setCounts, setCategoryOffset, setCategoryTestOffset, totalAmplitudes,
  numForestAmplitudes, amplitudes, forestTotalFits, fitsWithoutOffset,
  forestVariableCounts, numTreesInForest,
  count
};

constexpr std::size_t numFacadeVirtuals =
  static_cast<std::size_t>(FacadeVirtual::count);

/// The shared record of which virtual each row reached. Shared rather than
/// per-spy so a row that REBUILDS its fixture (setData replaces the store) does
/// not throw away its own evidence.
struct CallLog {
  std::vector<int> lastRow = std::vector<int>(numFacadeVirtuals, -1);
  int currentRow = -1;

  void mark(FacadeVirtual which) {
    lastRow[static_cast<std::size_t>(which)] = currentRow;
  }
  bool calledInRow(FacadeVirtual which, int row) const {
    return lastRow[static_cast<std::size_t>(which)] == row;
  }
  bool everCalled(FacadeVirtual which) const {
    return lastRow[static_cast<std::size_t>(which)] >= 0;
  }
};

/// A pass-through SamplerBase that records which virtual a row reached and
/// forwards the call unchanged to the facade under test. It exists for the
/// coverage gate, not to intercept anything: every argument and every return
/// value travels verbatim, so a row observes the facade's own behavior.
class SpySampler final : public SamplerBase {
public:
  SpySampler(SamplerBase& inner, CallLog& log) : inner_(inner), log_(log) {}

#define SPY_VOID(name, params, args)                                          \
  void name params override {                                                 \
    record(FacadeVirtual::name);                                              \
    inner_.name args;                                                         \
  }
#define SPY_RET(type, name, params, args)                                     \
  type name params override {                                                 \
    record(FacadeVirtual::name);                                              \
    return inner_.name args;                                                  \
  }

  SPY_RET(SamplerShape, shape, () const, ())
  SPY_RET(bool, run,
          (std::size_t b, std::size_t s, Results& r,
           const std::function<bool()>& p, const SweepCallback& c),
          (b, s, r, p, c))
  SPY_VOID(setOffset, (const double* o, bool u), (o, u))
  SPY_VOID(setResponse, (const double* y, bool u), (y, u))
  SPY_VOID(setWeights, (const double* w), (w))
  SPY_RET(std::uint64_t, weightsDigest, () const, ())
  SPY_VOID(reapplyWeights, (), ())
  SPY_VOID(setSigma, (double s), (s))
  SPY_RET(bool, setTestData, (const PredictorSource& s), (s))
  SPY_VOID(setTestOffset, (const double* o), (o))
  SPY_VOID(setData,
           (const double* x, const double* y, std::size_t n, const double* w,
            const double* o, const double* xt, std::size_t nt,
            const double* ot),
           (x, y, n, w, o, xt, nt, ot))
  SPY_RET(PredictorUpdateResult, setPredictor,
          (const PredictorSource& s, bool f, bool u), (s, f, u))
  SPY_RET(PredictorUpdateResult, updatePredictor,
          (const PredictorSource& s, const std::size_t* c, std::size_t nc,
           bool f, bool u),
          (s, c, nc, f, u))
  SPY_VOID(setCutPoints,
           (const double* const* nc, const std::uint32_t* n,
            const std::size_t* c, std::size_t ncol, const double* cur),
           (nc, n, c, ncol, cur))
  SPY_RET(bool, updatePredictorPerObservation,
          (const double* c, std::size_t j, bool* i), (c, j, i))
  SPY_RET(std::unique_ptr<PredictorUpdateSession>, beginPredictorUpdate,
          (const double* c, std::size_t j), (c, j))
  SPY_RET(std::size_t, currentSampleNum, () const, ())
  SPY_RET(std::size_t, savedSlotForDraw, (std::size_t d) const, (d))
  SPY_RET(const std::vector<FlatNode>&, savedTree,
          (std::size_t c, std::size_t s, std::size_t t, std::size_t f) const,
          (c, s, t, f))
  SPY_RET(const std::vector<double>&, savedTreeSlopes,
          (std::size_t c, std::size_t s, std::size_t t, std::size_t f) const,
          (c, s, t, f))
  SPY_RET(const std::vector<std::uint64_t>&, savedTreeMasks,
          (std::size_t c, std::size_t s, std::size_t t, std::size_t f) const,
          (c, s, t, f))
  SPY_VOID(flattenTree,
           (std::size_t c, std::size_t t, std::vector<FlatNode>& n,
            std::vector<std::uint32_t>& ct, std::vector<double>* sl,
            std::vector<std::uint64_t>* m, std::size_t f),
           (c, t, n, ct, sl, m, f))
  // The three replay virtuals are spelled out rather than macro-generated so
  // the spy can RECORD the numThreads it was handed. Recording it here proves
  // the argument survives the type-erased hop; the engine-side worker count
  // (bartcore::predictPartition) proves the facade then forwarded it rather
  // than passing a constant of its own.
  void predict(const PredictorSource& s, std::size_t nt, const double* co,
               std::size_t numThreads, double* o) override {
    record(FacadeVirtual::predict);
    predictThreads = numThreads;
    inner_.predict(s, nt, co, numThreads, o);
  }
  void predictPerForest(const PredictorSource& s, std::size_t nt,
                        std::size_t numThreads, double* o) override {
    record(FacadeVirtual::predictPerForest);
    predictPerForestThreads = numThreads;
    inner_.predictPerForest(s, nt, numThreads, o);
  }
  void predictVariance(const PredictorSource& s, std::size_t nt,
                       std::size_t numThreads, double* o) override {
    record(FacadeVirtual::predictVariance);
    predictVarianceThreads = numThreads;
    inner_.predictVariance(s, nt, numThreads, o);
  }

  std::size_t predictThreads = 0;
  std::size_t predictPerForestThreads = 0;
  std::size_t predictVarianceThreads = 0;
  SPY_VOID(getState, (SamplerStateData& s), (s))
  SPY_RET(bool, setState,
          (const SamplerStateData& s, const double* cp, bool* r), (s, cp, r))
  SPY_RET(WarmStartResult, installForests,
          (const SamplerStateData& d,
           const std::vector<std::pair<std::size_t, int>>& m),
          (d, m))
  SPY_VOID(sampleTreesFromPrior, (), ())
  SPY_VOID(sampleNodeParametersFromPrior, (), ())
  SPY_VOID(growFromRoot, (std::size_t s), (s))
  SPY_VOID(setNumThreads, (std::size_t t), (t))
  SPY_VOID(setNumThin, (std::size_t t), (t))
  SPY_VOID(setVerbose, (bool v, std::size_t e), (v, e))
  SPY_RET(double, fitScale, () const, ())
  SPY_VOID(setTreeStorage, (bool k, std::size_t n), (k, n))
  SPY_VOID(setModel, (const ModelParameters& m), (m))
  SPY_RET(double, sumOfSquaredResiduals, (std::size_t c), (c))
  SPY_VOID(printTrees,
           (const std::size_t* ci, std::size_t nci, const std::size_t* si,
            std::size_t nsi, const std::size_t* ti, std::size_t nti,
            std::size_t f, bool live),
           (ci, nci, si, nsi, ti, nti, f, live))
  SPY_RET(ext_rng*, rng, () const, ())
  SPY_RET(const ColumnStore&, data, () const, ())
  SPY_RET(const double*, latents, (std::size_t c) const, (c))
  SPY_RET(double, sigma, (std::size_t c) const, (c))
  SPY_RET(double, dispersion, (std::size_t c) const, (c))
  SPY_RET(bool, setForestBasis,
          (std::size_t f, const double* v, std::size_t n), (f, v, n))
  SPY_RET(bool, setForestWeights, (std::size_t f, const double* w), (f, w))
  SPY_RET(ForestCalibration, forestCalibration,
          (std::size_t c, std::size_t f) const, (c, f))
  SPY_RET(bool, setForestPriorScale, (std::size_t f, double s), (f, s))
  SPY_RET(bool, setActiveRows, (const double* a), (a))
  SPY_RET(bool, setCounts, (const int* c, const int* t), (c, t))
  SPY_RET(bool, setCategoryOffset, (const double* o), (o))
  SPY_RET(bool, setCategoryTestOffset, (const double* o), (o))
  SPY_RET(std::size_t, totalAmplitudes, () const, ())
  SPY_RET(std::size_t, numForestAmplitudes, (std::size_t f) const, (f))
  SPY_RET(bool, amplitudes, (std::size_t c, double* o) const, (c, o))
  SPY_VOID(forestTotalFits,
           (std::size_t c, std::size_t f, double* o) const, (c, f, o))
  SPY_RET(bool, fitsWithoutOffset, (std::size_t c, double* o), (c, o))
  SPY_VOID(forestVariableCounts,
           (std::size_t c, std::size_t f, std::uint32_t* o) const, (c, f, o))
  SPY_RET(std::size_t, numTreesInForest, (std::size_t f) const, (f))

#undef SPY_VOID
#undef SPY_RET

private:
  /// Marks the row in force; const because most of the boundary is.
  void record(FacadeVirtual which) const { log_.mark(which); }

  SamplerBase& inner_;
  CallLog& log_;
};

CallLog callLog;


/// One sampler held three ways: the concrete facade, the spy the rows call
/// through, and the typed impl they check the answer against.
template <typename L = ConstantGaussianLeaf>
struct FixtureT {
  std::unique_ptr<SamplerFacade<L>> facade;
  std::unique_ptr<SpySampler> spy;

  template <typename... Args>
  void build(Args&&... args) {
    facade = std::make_unique<SamplerFacade<L>>(std::forward<Args>(args)...);
    spy = std::make_unique<SpySampler>(*facade, callLog);
  }
  SamplerBase& base() { return *spy; }
  Sampler<L>& impl() { return facade->impl(); }
};

using Fixture = FixtureT<ConstantGaussianLeaf>;

/// The samplers the rows share. Each is built so that the quantities a
/// forwarder could confuse DIFFER: two chains, a partly filled saved-tree store
/// (retained draws 2, capacity 3), a combining fixture whose two forests carry
/// unequal tree counts (6 and 4) and unequal amplitude widths (1 and 2), and a
/// family fixture for every capability the boundary gates on.
struct Fixtures {
  static constexpr std::size_t n = 120, p = 2, nTest = 6, K = 3, capacity = 3;

  std::vector<double> x, y, y3, xTest, xTest2, offset, offset3, weights,
    yBinary, yCount, unitBasis, wideBasis, wideBasis2, categoryOffset,
    testCategoryOffset, testOffset, newColumn, replacement;
  std::vector<int> counts, trials, counts0;
  std::vector<double> xPooled, newColumn2;
  std::vector<ext_rng*> rngs;
  Fixture g, gt, d, b, m, v, l, nb;
  FixtureT<LinearGaussianLeaf> lin;
  std::size_t leafCovariate = 0;
  SamplerOptions gaussianOptions;

  ext_rng* newRng(std::uint_least32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  }

  Fixtures() {
    x.resize(n * p);
    y.resize(n);
    y3.resize(n);
    xTest.resize(nTest * p);
    xTest2.resize(4 * p);
    offset.resize(n);
    offset3.resize(n);
    weights.resize(n);
    yBinary.resize(n);
    yCount.resize(n);
    unitBasis.assign(n, 1.0);
    wideBasis.resize(2 * n);
    wideBasis2.assign(2 * n, 3.0);
    categoryOffset.assign(n * K, 0.0);
    testCategoryOffset.assign(nTest * K, 0.0);
    testOffset.assign(nTest, 1.0e6);
    newColumn.resize(n);
    newColumn2.resize(n);
    xPooled.resize(n * p);
    replacement.resize(n * p);
    counts.assign(n * K, 0);
    counts0.assign(n * K, 0);
    trials.assign(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = runif01();
      x[i + n] = runif01();
      y[i] = 2.0 * x[i] - x[i + n] + 0.3 * (runif01() - 0.5);
      y3[i] = 3.0 * y[i];
      offset[i] = 0.1 * static_cast<double>(i % 7);
      offset3[i] = 4.0 * std::sin(0.2 * static_cast<double>(i));
      weights[i] = 1.0 + static_cast<double>(i % 3);
      yBinary[i] = y[i] > 0.0 ? 1.0 : 0.0;
      yCount[i] = static_cast<double>(i % 4);
      newColumn[i] = runif01();
      newColumn2[i] = runif01();
      xPooled[i] = runif01();
      // 70 levels, so the column pools its masks and the saved-tree mask
      // channel exists at all
      xPooled[i + n] = static_cast<double>(i % 70);
      replacement[i] = runif01();
      replacement[i + n] = runif01();
      wideBasis[2 * i] = 1.0;
      wideBasis[2 * i + 1] = x[i] > 0.5 ? 1.0 : 0.0;
      counts[(i % K) * n + i] = 1;
      counts0[i] = 1;  // every observation in category 0
    }
    for (double& value : xTest) value = runif01();
    for (double& value : xTest2) value = runif01();

    gaussianOptions.numTrees = 8;
    gaussianOptions.numChains = 2;
    gaussianOptions.keepTrees = true;
    gaussianOptions.numSamplesToStore = capacity;

    buildGaussian(g, 51001u);
    buildGaussian(gt, 51002u);
    buildDisposable();
    buildCombining();
    buildMultinomial();
    buildVariance();
    buildLogistic();
    buildNegativeBinomial();
    buildVectorLeaf();
  }

  ~Fixtures() {
    for (ext_rng* r : rngs) ext_rng_destroy(r);
  }

  /// The main fixture: two chains, a saved-tree store of capacity 3 holding
  /// two recorded draws, so retained draws and capacity differ from the first
  /// row on.
  void buildGaussian(Fixture& fixture, std::uint_least32_t seed) {
    ext_rng* pair[2] = {newRng(seed), newRng(seed + 100u)};
    fixture.build(x.data(), y.data(), n, p, nullptr, nullptr,
                  ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542,
                  gaussianOptions, pair);
    fixture.impl().setTestPredictors(xTest.data(), nTest);
    Results results;
    fixture.impl().run(10, 2, results);
  }

  /// A single-chain sampler the destructive rows (the predictor conduit, the
  /// model and thread/thin/verbose settings, setData) work on, so the saved
  /// store the read rows depend on is never disturbed.
  void buildDisposable() {
    SamplerOptions options;
    options.numTrees = 4;
    ext_rng* one = newRng(51003u);
    d.build(x.data(), y.data(), n, p, nullptr, nullptr,
            ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options,
            &one);
    Results results;
    d.impl().run(5, 1, results);
  }

  void buildCombining() {
    SamplerOptions options;
    options.numChains = 2;
    options.keepTrees = true;
    options.numSamplesToStore = 2;
    AmplitudeSpec spec;
    spec.forests.resize(2);
    spec.forests[0].forest.numTrees = 6;
    spec.forests[0].basis = unitBasis.data();
    spec.forests[0].numBasisColumns = 1;
    spec.forests[0].nodeScaleFactor = 1.25;
    spec.forests[1].forest.numTrees = 4;
    spec.forests[1].basis = wideBasis.data();
    spec.forests[1].numBasisColumns = 2;
    spec.forests[1].nodeScaleFactor = 0.75;
    ext_rng* pair[2] = {newRng(51004u), newRng(51005u)};
    b.build(x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
            0.37804942330213542, options, spec, pair);
    Results results;
    b.impl().run(10, 2, results);
  }

  void buildMultinomial() {
    SamplerOptions options;
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = counts.data();
    spec.trials = trials.data();
    spec.forest.numTrees = 5;
    ext_rng* one = newRng(51006u);
    m.build(x.data(), n, p, options, spec, &one);
    m.impl().setTestPredictors(xTest.data(), nTest);
    Results results;
    m.impl().run(5, 1, results);
  }

  void buildVariance() {
    SamplerOptions options;
    options.numTrees = 6;
    options.numVarianceTrees = 6;
    options.keepTrees = true;
    options.numSamplesToStore = 2;
    ext_rng* one = newRng(51007u);
    v.build(x.data(), y.data(), n, p, nullptr, nullptr,
            ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options,
            &one);
    Results results;
    v.impl().run(10, 2, results);
  }

  void buildLogistic() {
    SamplerOptions options;
    options.numTrees = 6;
    options.numChains = 2;
    ext_rng* pair[2] = {newRng(51008u), newRng(51009u)};
    l.build(x.data(), yBinary.data(), n, p, weights.data(), nullptr,
            ResponseFamily::logistic, 1.0, 3.0, 0.37804942330213542, options,
            pair);
    Results results;
    l.impl().run(5, 1, results);
  }

  /// The one fixture whose saved store carries the two side channels: a vector
  /// leaf sizes the slopes, a pooled categorical column the masks. Without it
  /// those two readers have nothing to address.
  void buildVectorLeaf() {
    SamplerOptions options;
    options.numTrees = 4;
    options.keepTrees = true;
    options.numSamplesToStore = 2;
    options.leafCovariateColumns = &leafCovariate;
    options.numLeafCovariates = 1;
    std::vector<ColumnType> types = {ColumnType::ordinal,
                                     ColumnType::categorical};
    options.predictors.columnTypes = types.data();
    ext_rng* one = newRng(51011u);
    lin.build(xPooled.data(), y.data(), n, p, nullptr, nullptr,
              ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options,
              &one);
    Results results;
    lin.impl().run(5, 2, results);
  }

  void buildNegativeBinomial() {
    SamplerOptions options;
    options.numTrees = 6;
    options.dispersion = 2.5;  // positive: r is fixed there
    ext_rng* one = newRng(51010u);
    nb.build(x.data(), yCount.data(), n, p, nullptr, nullptr,
             ResponseFamily::nbinom, 1.0, 3.0, 0.37804942330213542, options,
             &one);
    Results results;
    nb.impl().run(5, 1, results);
  }
};

/// One virtual, the row that exercises it, and the name its failures carry.
struct Row {
  FacadeVirtual which;
  const char* name;
  void (*body)(Fixtures&);
};

/// Sum of a variable-count channel, the observable several rows read a forest's
/// split usage through.
std::size_t totalSplits(const std::vector<std::uint32_t>& counts) {
  std::size_t total = 0;
  for (std::uint32_t value : counts) total += value;
  return total;
}

/// Every tree of forest 0, chain 0, in one fingerprint: a structural row wants
/// "some tree moved", which one tree alone does not answer.
std::uint64_t forestSignature(Fixture& fixture) {
  std::uint64_t hash = 0;
  for (std::size_t t = 0; t < fixture.impl().numTrees(); ++t)
    hash = hash * 31u +
           treeStructureSignature(fixture.impl().chain(0).tree(t));
  return hash;
}

/// One tree's dump through whichever SamplerBase the caller hands over, so a
/// row can compare the boundary's rendering with the impl's own.
std::string dumpTrees(SamplerBase& sampler, std::size_t forestIndex) {
  const std::size_t chain = 0, slot = 0, tree = 0;
  std::string text;
  beginPrintCapture(text);
  sampler.printTrees(&chain, 1, &slot, 1, &tree, 1, forestIndex, false);
  endPrintCapture();
  return text;
}

/// The table. Rows run in order and each leaves the fixture it touched usable
/// by the rows after it; where a row needs a particular store state (a partly
/// filled store, a wrapped one) it says so and arranges it itself.
const Row rows[] = {
  {FacadeVirtual::shape, "shape", [](Fixtures& f) {
    // the five fields test_shape.cpp's oracle does not carry, each read where
    // it differs from the quantity a mis-fill would take it from
    SamplerShape s = f.g.base().shape();
    check(s.numSavedDraws == f.g.impl().filledSavedDraws() &&
            s.numSavedDraws != s.savedTreeCapacity,
          "facade shape: the retained draw count is not the store capacity");
    check(f.b.base().shape().numAmplitudes == f.b.impl().totalAmplitudes() &&
            f.b.base().shape().numAmplitudes == 3,
          "facade shape: the amplitude total is the combiner's");
    check(f.nb.base().shape().carriesDispersion &&
            !f.g.base().shape().carriesDispersion,
          "facade shape: the dispersion flag is the family's");
    check(f.b.base().shape().forestReportingIsDefined &&
            !f.g.base().shape().forestReportingIsDefined,
          "facade shape: the per-forest reporting flag is the coupling's");
    const ResidualPrior& prior = f.v.impl().varianceLeafPrior();
    ResidualPrior reported = f.v.base().shape().varianceLeafPrior;
    check(reported.sigmaEstimate == prior.sigmaEstimate &&
            reported.sigmaDf == prior.sigmaDf &&
            reported.sigmaRawScale == prior.sigmaRawScale &&
            prior.sigmaDf > 0.0,
          "facade shape: the variance leaf prior is the forest's own");
  }},
  {FacadeVirtual::run, "run", [](Fixtures& f) {
    std::size_t filled = f.g.impl().filledSavedDraws();
    std::size_t cursor = f.g.impl().currentSampleNum();
    std::size_t sweeps = 0;
    Results results;
    // the callback and the interrupt poll are arguments of this virtual too:
    // one fires per sweep, the other never asks to stop
    bool stoppedEarly = f.g.base().run(
      1, 1, results, []() { return false; },
      [&sweeps](std::size_t, std::size_t, bool) {
        ++sweeps;
        return false;
      });
    check(!stoppedEarly && sweeps == 2 * f.g.impl().numChains(),
          "facade run: every chain sweeps twice and nothing stops early");
    check(f.g.impl().filledSavedDraws() == filled + 1 &&
            f.g.impl().currentSampleNum() ==
              (cursor + 1) % Fixtures::capacity,
          "facade run: the recorded draw lands in the impl's own store");
  }},
  {FacadeVirtual::setOffset, "setOffset", [](Fixtures& f) {
    double scale = f.g.impl().fitScale();
    double residuals = f.g.impl().sumOfSquaredResiduals(0);
    f.g.base().setOffset(f.offset.data(), false);
    check(f.g.impl().fitScale() == scale &&
            f.g.impl().sumOfSquaredResiduals(0) != residuals,
          "facade setOffset: at updateScale false the offset moves the "
          "residual and not the transform");
    f.g.base().setOffset(f.offset3.data(), true);
    check(f.g.impl().fitScale() != scale,
          "facade setOffset: at updateScale true the transform re-anchors");
    f.g.base().setOffset(nullptr, true);
  }},
  {FacadeVirtual::setResponse, "setResponse", [](Fixtures& f) {
    double scale = f.g.impl().fitScale();
    double residuals = f.g.impl().sumOfSquaredResiduals(0);
    f.g.base().setResponse(f.y3.data(), false);
    check(f.g.impl().fitScale() == scale &&
            f.g.impl().sumOfSquaredResiduals(0) != residuals,
          "facade setResponse: at updateScale false the response moves the "
          "residual and not the transform");
    f.g.base().setResponse(f.y3.data(), true);
    check(f.g.impl().fitScale() != scale,
          "facade setResponse: at updateScale true the transform re-anchors");
    f.g.base().setResponse(f.y.data(), true);
  }},
  {FacadeVirtual::setWeights, "setWeights", [](Fixtures& f) {
    std::uint64_t digest = f.g.impl().weightsDigest();
    f.g.base().setWeights(f.weights.data());
    check(f.g.impl().weightsDigest() != digest,
          "facade setWeights: the impl carries the installed weights");
    f.g.base().setWeights(nullptr);
    check(f.g.impl().weightsDigest() == digest,
          "facade setWeights: clearing restores the unweighted digest");
  }},
  {FacadeVirtual::weightsDigest, "weightsDigest", [](Fixtures& f) {
    check(f.g.base().weightsDigest() == f.g.impl().weightsDigest() &&
            f.g.base().weightsDigest() != 0,
          "facade weightsDigest: the boundary reports the impl's digest");
  }},
  {FacadeVirtual::reapplyWeights, "reapplyWeights", [](Fixtures& f) {
    // a latent family re-derives its augmentation from the weights in force,
    // so the omega draws move
    std::vector<double> before(f.l.impl().latents(0),
                               f.l.impl().latents(0) + Fixtures::n);
    f.l.base().reapplyWeights();
    bool moved = false;
    for (std::size_t i = 0; i < Fixtures::n; ++i)
      moved |= f.l.impl().latents(0)[i] != before[i];
    check(moved, "facade reapplyWeights: the family redraws its latents");
  }},
  {FacadeVirtual::setSigma, "setSigma", [](Fixtures& f) {
    f.g.base().setSigma(2.5);
    checkNear(f.g.impl().sigma(0), 2.5, 1.0e-12,
              "facade setSigma: chain 0 carries the installed sigma");
    checkNear(f.g.impl().sigma(1), 2.5, 1.0e-12,
              "facade setSigma: chain 1 carries it too");
  }},
  {FacadeVirtual::setTestData, "setTestData", [](Fixtures& f) {
    check(f.g.base().setTestData(
            densePredictorSource(f.xTest2.data(), 4, Fixtures::p)) &&
            f.g.impl().numTestObservations() == 4,
          "facade setTestData: the impl's test store takes the new rows");
    f.g.base().setTestData(
      densePredictorSource(f.xTest.data(), Fixtures::nTest, Fixtures::p));
    check(f.g.impl().numTestObservations() == Fixtures::nTest,
          "facade setTestData: and takes them back");
  }},
  {FacadeVirtual::setTestOffset, "setTestOffset", [](Fixtures& f) {
    // 1e6 dominates any fit, so one recorded draw separates an installed
    // offset from a dropped one without a twin sampler
    std::vector<double> testFits(Fixtures::nTest * 2, 0.0);
    Results results;
    results.testFits = testFits.data();
    f.g.base().setTestOffset(f.testOffset.data());
    f.g.impl().run(0, 1, results);
    bool shifted = true;
    for (double value : testFits) shifted &= value > 1.0e5;
    check(shifted, "facade setTestOffset: the recorded test fits carry it");
    f.g.base().setTestOffset(nullptr);
  }},
  {FacadeVirtual::setData, "setData", [](Fixtures& f) {
    f.d.base().setData(f.x.data(), f.y.data(), Fixtures::n / 2, nullptr,
                       nullptr, nullptr, 0, nullptr);
    check(f.d.impl().numObservations() == Fixtures::n / 2 &&
            f.d.impl().data().numObservations == Fixtures::n / 2,
          "facade setData: the impl's store takes the replacement");
    f.buildDisposable();  // the rows after this one want the full store
  }},
  {FacadeVirtual::setPredictor, "setPredictor", [](Fixtures& f) {
    std::vector<xint_t> codes = f.d.impl().data().train.codes;
    check(f.d.base().setPredictor(
            densePredictorSource(f.replacement.data(), Fixtures::n,
                                 Fixtures::p),
            true, true) == PredictorUpdateResult::accepted,
          "facade setPredictor: the replacement is accepted");
    check(f.d.impl().data().train.codes != codes,
          "facade setPredictor: the impl's codes are the new ones");
  }},
  {FacadeVirtual::updatePredictor, "updatePredictor", [](Fixtures& f) {
    std::vector<xint_t> codes = f.d.impl().data().train.codes;
    const std::size_t columns[] = {1};
    check(f.d.base().updatePredictor(
            densePredictorSource(f.newColumn.data(), Fixtures::n, 1), columns,
            1, true, true) == PredictorUpdateResult::accepted,
          "facade updatePredictor: the named column is accepted");
    const std::vector<xint_t>& updated = f.d.impl().data().train.codes;
    bool columnZeroHeld = true, columnOneMoved = false;
    for (std::size_t i = 0; i < Fixtures::n; ++i) {
      columnZeroHeld &= updated[i] == codes[i];
      columnOneMoved |= updated[i + Fixtures::n] != codes[i + Fixtures::n];
    }
    check(columnZeroHeld && columnOneMoved,
          "facade updatePredictor: column 1 moves and column 0 does not");
  }},
  {FacadeVirtual::setCutPoints, "setCutPoints", [](Fixtures& f) {
    const double cuts[] = {0.2, 0.5, 0.8};
    const double* cutPointers[] = {cuts};
    const std::uint32_t numCuts[] = {3};
    const std::size_t columns[] = {0};
    f.d.base().setCutPoints(cutPointers, numCuts, columns, 1,
                            f.replacement.data());
    check(f.d.impl().data().numCuts[0] == 3 &&
            f.d.impl().data().cutPoints[0] ==
              std::vector<double>(cuts, cuts + 3),
          "facade setCutPoints: the impl's grid takes the named column");
  }},
  {FacadeVirtual::updatePredictorPerObservation, "perObservation",
   [](Fixtures& f) {
    std::vector<xint_t> codes = f.d.impl().data().train.codes;
    std::unique_ptr<bool[]> installed(new bool[Fixtures::n]);
    check(f.d.base().updatePredictorPerObservation(f.newColumn2.data(), 1,
                                                   installed.get()),
          "facade perObservation: the sweep finalizes");
    std::size_t numInstalled = 0;
    for (std::size_t i = 0; i < Fixtures::n; ++i)
      numInstalled += installed[i] ? 1 : 0;
    check(numInstalled > 0 && f.d.impl().data().train.codes != codes,
          "facade perObservation: the installed rows reach the impl's codes");
  }},
  {FacadeVirtual::beginPredictorUpdate, "beginPredictorUpdate",
   [](Fixtures& f) {
    std::unique_ptr<PredictorUpdateSession> session =
      f.d.base().beginPredictorUpdate(f.newColumn.data(), 1);
    check(session != nullptr,
          "facade beginPredictorUpdate: the boundary hands back a session");
    std::vector<xint_t> codes = f.d.impl().data().train.codes;
    std::size_t committed = 0;
    for (std::size_t i = 0; i < Fixtures::n; ++i)
      if (session->observationWouldRemainValid(i)) {
        session->commitObservation(i);
        ++committed;
      }
    check(session->finalize() && committed > 0 &&
            f.d.impl().data().train.codes != codes,
          "facade beginPredictorUpdate: the session drives the impl's store");
  }},
  {FacadeVirtual::currentSampleNum, "currentSampleNum", [](Fixtures& f) {
    Results results;
    // off zero, so a forwarder answering with a constant differs
    while (f.g.impl().currentSampleNum() == 0) f.g.impl().run(0, 1, results);
    check(f.g.base().currentSampleNum() == f.g.impl().currentSampleNum() &&
            f.g.base().currentSampleNum() != 0,
          "facade currentSampleNum: the boundary reports the write cursor");
  }},
  {FacadeVirtual::savedSlotForDraw, "savedSlotForDraw", [](Fixtures& f) {
    Results results;
    // a full store read from a nonzero cursor: the map is not the identity
    while (f.g.impl().currentSampleNum() == 0) f.g.impl().run(0, 1, results);
    check(f.g.impl().filledSavedDraws() == Fixtures::capacity,
          "facade savedSlotForDraw: the store is full for this read");
    bool agrees = true, identity = true;
    for (std::size_t i = 0; i < f.g.impl().filledSavedDraws(); ++i) {
      agrees &= f.g.base().savedSlotForDraw(i) == f.g.impl().savedSlotForDraw(i);
      identity &= f.g.base().savedSlotForDraw(i) == i;
    }
    check(agrees && !identity,
          "facade savedSlotForDraw: the ring map is the impl's, not identity");
  }},
  {FacadeVirtual::savedTree, "savedTree", [](Fixtures& f) {
    std::size_t slot = f.g.impl().savedSlotForDraw(0);
    check(&f.g.base().savedTree(1, slot, 0, 0) ==
            &f.g.impl().savedTree(1, slot, 0, 0),
          "facade savedTree: the named chain's slot is read");
    check(&f.b.base().savedTree(0, 0, 0, 1) ==
            &f.b.impl().savedTree(0, 0, 0, 1) &&
            &f.b.base().savedTree(0, 0, 0, 1) !=
              &f.b.base().savedTree(0, 0, 0, 0),
          "facade savedTree: the named forest is read, not forest 0");
  }},
  {FacadeVirtual::savedTreeSlopes, "savedTreeSlopes", [](Fixtures& f) {
    // only a vector-leaf store sizes this channel, so the vector-leaf fixture
    // is where a slot and tree index mean anything
    check(&f.lin.base().savedTreeSlopes(0, 1, 2, 0) ==
            &f.lin.impl().savedTreeSlopes(0, 1, 2, 0) &&
            &f.lin.base().savedTreeSlopes(0, 1, 2, 0) !=
              &f.lin.base().savedTreeSlopes(0, 0, 2, 0) &&
            &f.lin.base().savedTreeSlopes(0, 1, 2, 0) !=
              &f.lin.base().savedTreeSlopes(0, 1, 0, 0),
          "facade savedTreeSlopes: the named slot and tree are read");
  }},
  {FacadeVirtual::savedTreeMasks, "savedTreeMasks", [](Fixtures& f) {
    // and only a pooled categorical store sizes this one
    check(&f.lin.base().savedTreeMasks(0, 1, 2, 0) ==
            &f.lin.impl().savedTreeMasks(0, 1, 2, 0) &&
            &f.lin.base().savedTreeMasks(0, 1, 2, 0) !=
              &f.lin.base().savedTreeMasks(0, 0, 2, 0) &&
            &f.lin.base().savedTreeMasks(0, 1, 2, 0) !=
              &f.lin.base().savedTreeMasks(0, 1, 0, 0),
          "facade savedTreeMasks: the named slot and tree are read");
  }},
  {FacadeVirtual::flattenTree, "flattenTree", [](Fixtures& f) {
    std::vector<FlatNode> viaBase, viaImpl, forestZero;
    std::vector<std::uint32_t> countsBase, countsImpl, countsZero;
    f.b.base().flattenTree(0, 0, viaBase, countsBase, nullptr, nullptr, 1);
    f.b.impl().flattenTree(0, 0, viaImpl, countsImpl, nullptr, nullptr, 1);
    f.b.base().flattenTree(0, 0, forestZero, countsZero, nullptr, nullptr, 0);
    check(countsBase == countsImpl && viaBase.size() == viaImpl.size(),
          "facade flattenTree: the boundary flattens what the impl does");
    check(countsBase != countsZero || viaBase.size() != forestZero.size(),
          "facade flattenTree: the named forest is flattened, not forest 0");
  }},
  // The replay rows check two things at once: the fits agree with a direct
  // call on the impl, AND the thread count they were handed - a value chosen
  // to differ from every default and from the fixtures' own count - arrives at
  // the engine. A forwarder that drops the argument leaves the engine
  // resolving the sampler's own count instead, which the resolved-count check
  // catches; a boundary that never passed it on leaves the spy's record wrong.
  {FacadeVirtual::predict, "predict", [](Fixtures& f) {
    std::size_t slab = Fixtures::nTest * f.g.impl().filledSavedDraws() *
                       f.g.impl().numChains();
    std::vector<double> viaBase(slab, 0.0), viaImpl(slab, 1.0);
    PredictorSource source =
      densePredictorSource(f.xTest.data(), Fixtures::nTest, Fixtures::p);
    f.g.base().predict(source, Fixtures::nTest, nullptr, 5, viaBase.data());
    std::size_t resolved = predictPartition.resolvedThreads;
    f.g.impl().predict(source, Fixtures::nTest, nullptr, 5, viaImpl.data());
    check(viaBase == viaImpl,
          "facade predict: the boundary's fits are the impl's");
    check(f.g.spy->predictThreads == 5 && resolved == 5,
          "facade predict: the thread count crosses the boundary intact");
  }},
  {FacadeVirtual::predictPerForest, "predictPerForest", [](Fixtures& f) {
    std::size_t slab = Fixtures::nTest * f.b.impl().numForests() *
                       f.b.impl().filledSavedDraws() * f.b.impl().numChains();
    std::vector<double> viaBase(slab, 0.0), viaImpl(slab, 1.0);
    PredictorSource source =
      densePredictorSource(f.xTest.data(), Fixtures::nTest, Fixtures::p);
    f.b.base().predictPerForest(source, Fixtures::nTest, 5, viaBase.data());
    std::size_t resolved = predictPartition.resolvedThreads;
    f.b.impl().predictPerForest(source, Fixtures::nTest, 5, viaImpl.data());
    check(viaBase == viaImpl,
          "facade predictPerForest: the boundary's per-forest fits are the "
          "impl's");
    check(f.b.spy->predictPerForestThreads == 5 && resolved == 5,
          "facade predictPerForest: the thread count crosses the boundary "
          "intact");
  }},
  {FacadeVirtual::predictVariance, "predictVariance", [](Fixtures& f) {
    std::size_t slab = Fixtures::nTest * f.v.impl().filledSavedDraws();
    std::vector<double> viaBase(slab, 0.0), viaImpl(slab, 1.0);
    PredictorSource source =
      densePredictorSource(f.xTest.data(), Fixtures::nTest, Fixtures::p);
    f.v.base().predictVariance(source, Fixtures::nTest, 5, viaBase.data());
    std::size_t resolved = predictPartition.resolvedThreads;
    f.v.impl().predictVariance(source, Fixtures::nTest, 5, viaImpl.data());
    bool positive = true;
    for (double value : viaBase) positive &= value > 0.0;
    check(viaBase == viaImpl && positive,
          "facade predictVariance: the boundary's surface is the impl's");
    check(f.v.spy->predictVarianceThreads == 5 && resolved == 5,
          "facade predictVariance: the thread count crosses the boundary "
          "intact");
  }},
  {FacadeVirtual::getState, "getState", [](Fixtures& f) {
    SamplerStateData viaBase, viaImpl;
    f.g.base().getState(viaBase);
    f.g.impl().getState(viaImpl);
    check(statesAgree(viaBase, viaImpl),
          "facade getState: the boundary captures the impl's own state");
  }},
  {FacadeVirtual::setState, "setState", [](Fixtures& f) {
    SamplerStateData donor;
    f.g.impl().getState(donor);
    bool columnMaskRefused = true;
    check(f.gt.base().setState(donor, nullptr, &columnMaskRefused) &&
            !columnMaskRefused,
          "facade setState: the donor installs and the refusal flag is "
          "written");
    SamplerStateData restored;
    f.gt.impl().getState(restored);
    check(statesAgree(donor, restored),
          "facade setState: the twin now reproduces the donor's state");
  }},
  {FacadeVirtual::installForests, "installForests", [](Fixtures& f) {
    Results results;
    f.gt.impl().run(0, 2, results);  // move the twin off the donor first
    SamplerStateData donor;
    f.g.impl().getState(donor);
    std::uint64_t target = treeStructureSignature(f.g.impl().chain(1).tree(0));
    check(treeStructureSignature(f.gt.impl().chain(1).tree(0)) != target,
          "facade installForests: the twin's trees start elsewhere");
    std::vector<std::pair<std::size_t, int>> map = {{0, -1}, {1, -1}};
    check(f.gt.base().installForests(donor, map) == WarmStartResult::ok &&
            treeStructureSignature(f.gt.impl().chain(1).tree(0)) == target,
          "facade installForests: the donor's live trees seed the twin");
  }},
  {FacadeVirtual::sampleTreesFromPrior, "sampleTreesFromPrior",
   [](Fixtures& f) {
    std::uint64_t before = forestSignature(f.d);
    f.d.base().sampleTreesFromPrior();
    check(forestSignature(f.d) != before,
          "facade sampleTreesFromPrior: the impl's trees are redrawn");
  }},
  {FacadeVirtual::sampleNodeParametersFromPrior, "nodeParametersFromPrior",
   [](Fixtures& f) {
    std::uint64_t structure = forestSignature(f.d);
    std::vector<double> before = f.d.impl().chain(0).treeFits();
    f.d.base().sampleNodeParametersFromPrior();
    check(f.d.impl().chain(0).treeFits() != before &&
            forestSignature(f.d) == structure,
          "facade nodeParametersFromPrior: the leaf values move and the "
          "structure does not");
  }},
  {FacadeVirtual::growFromRoot, "growFromRoot", [](Fixtures& f) {
    std::uint64_t before = forestSignature(f.d);
    f.d.base().growFromRoot(2);
    check(forestSignature(f.d) != before,
          "facade growFromRoot: the impl's forest is grown");
  }},
  {FacadeVirtual::setNumThreads, "setNumThreads", [](Fixtures& f) {
    f.d.base().setNumThreads(2);
    check(f.d.impl().numThreads() == 2,
          "facade setNumThreads: the impl holds the thread count");
    f.d.base().setNumThreads(1);
  }},
  {FacadeVirtual::setNumThin, "setNumThin", [](Fixtures& f) {
    // the sweep hook counts iterations, which is what thinning multiplies
    std::size_t sweeps = 0;
    Results results;
    auto count = [&sweeps](std::size_t, std::size_t, bool) {
      ++sweeps;
      return false;
    };
    f.d.base().setNumThin(3);
    f.d.base().run(0, 1, results, {}, count);
    check(sweeps == 3, "facade setNumThin: one kept draw costs three sweeps");
    f.d.base().setNumThin(1);
    sweeps = 0;
    f.d.base().run(0, 1, results, {}, count);
    check(sweeps == 1, "facade setNumThin: and one again when it is reset");
  }},
  {FacadeVirtual::setVerbose, "setVerbose", [](Fixtures& f) {
    Results results;
    std::string quiet, loud;
    beginPrintCapture(quiet);
    f.d.impl().run(0, 1, results);
    endPrintCapture();
    f.d.base().setVerbose(true, 1);
    beginPrintCapture(loud);
    f.d.impl().run(0, 1, results);
    endPrintCapture();
    check(quiet.empty() && !loud.empty(),
          "facade setVerbose: the impl's progress reporting turns on");
    f.d.base().setVerbose(false, 100);
  }},
  {FacadeVirtual::fitScale, "fitScale", [](Fixtures& f) {
    check(f.g.base().fitScale() == f.g.impl().fitScale() &&
            f.g.base().fitScale() != 1.0,
          "facade fitScale: the boundary reports the impl's transform");
  }},
  {FacadeVirtual::setTreeStorage, "setTreeStorage", [](Fixtures& f) {
    check(f.d.impl().savedTreeCapacity() == 0,
          "facade setTreeStorage: the disposable fixture starts without a "
          "store");
    f.d.base().setTreeStorage(true, 5);
    check(f.d.impl().savedTreeCapacity() == 5 &&
            f.d.impl().filledSavedDraws() == 0,
          "facade setTreeStorage: the impl takes the named capacity");
  }},
  {FacadeVirtual::setModel, "setModel", [](Fixtures& f) {
    ModelParameters model;
    model.k = 5.0;
    check(f.d.impl().k(0) != 5.0,
          "facade setModel: the fixture does not already carry the model");
    f.d.base().setModel(model);
    check(f.d.impl().k(0) == 5.0,
          "facade setModel: the impl's chains take the replacement prior");
  }},
  {FacadeVirtual::sumOfSquaredResiduals, "sumOfSquaredResiduals",
   [](Fixtures& f) {
    check(f.g.base().sumOfSquaredResiduals(1) ==
            f.g.impl().sumOfSquaredResiduals(1) &&
            f.g.base().sumOfSquaredResiduals(0) !=
              f.g.base().sumOfSquaredResiduals(1),
          "facade sumOfSquaredResiduals: the named chain answers");
  }},
  {FacadeVirtual::printTrees, "printTrees", [](Fixtures& f) {
    std::string viaBase = dumpTrees(f.b.base(), 1);
    std::string viaImpl = dumpTrees(*f.b.facade, 1);
    check(!viaBase.empty() && viaBase == viaImpl,
          "facade printTrees: the boundary's dump is the facade's own");
    check(viaBase != dumpTrees(f.b.base(), 0),
          "facade printTrees: the named forest is dumped, not forest 0");
  }},
  {FacadeVirtual::rng, "rng", [](Fixtures& f) {
    check(f.g.base().rng() == f.g.impl().rng() && f.g.base().rng() != nullptr,
          "facade rng: the boundary hands back the impl's generator");
  }},
  {FacadeVirtual::data, "data", [](Fixtures& f) {
    check(&f.g.base().data() == &f.g.impl().data(),
          "facade data: the boundary hands back the impl's store");
  }},
  {FacadeVirtual::latents, "latents", [](Fixtures& f) {
    check(f.l.base().latents(1) == f.l.impl().latents(1) &&
            f.l.base().latents(0) != f.l.base().latents(1),
          "facade latents: the named chain's latents answer");
  }},
  {FacadeVirtual::sigma, "sigma", [](Fixtures& f) {
    check(f.g.base().sigma(1) == f.g.impl().sigma(1) &&
            f.g.base().sigma(0) != f.g.base().sigma(1),
          "facade sigma: the named chain's sigma answers");
  }},
  {FacadeVirtual::dispersion, "dispersion", [](Fixtures& f) {
    check(f.nb.base().dispersion(0) == f.nb.impl().dispersion(0) &&
            f.nb.base().dispersion(0) == 2.5 &&
            f.g.base().dispersion(0) == 0.0,
          "facade dispersion: the family that carries one reports it");
  }},
  {FacadeVirtual::setForestBasis, "setForestBasis", [](Fixtures& f) {
    double before = f.b.impl().forestCalibration(0, 0).basisRowNorm;
    double target = f.b.impl().forestCalibration(0, 1).basisRowNorm;
    check(!f.b.base().setForestBasis(2, f.wideBasis2.data(), 2) &&
            !f.g.base().setForestBasis(0, f.unitBasis.data(), 1),
          "facade setForestBasis: an absent forest and an amplitude-free "
          "sampler refuse");
    check(f.b.base().setForestBasis(1, f.wideBasis2.data(), 2),
          "facade setForestBasis: the named forest installs");
    check(f.b.impl().forestCalibration(0, 1).basisRowNorm != target &&
            f.b.impl().forestCalibration(0, 0).basisRowNorm == before,
          "facade setForestBasis: forest 1 re-anchors and forest 0 does not");
  }},
  {FacadeVirtual::setForestWeights, "setForestWeights", [](Fixtures& f) {
    std::vector<double> zeros(Fixtures::n, 0.0);
    std::vector<std::uint32_t> forestZero(Fixtures::p, 0),
      forestOne(Fixtures::p, 0);
    check(!f.b.base().setForestWeights(2, zeros.data()) &&
            !f.g.base().setForestWeights(0, zeros.data()),
          "facade setForestWeights: an absent forest and a coupling-free "
          "sampler refuse");
    check(f.b.base().setForestWeights(1, zeros.data()),
          "facade setForestWeights: the named forest installs");
    // with no positive-weight row a forest's prior draw is all bare roots,
    // which is the per-forest observable: forest 1 stops splitting, forest 0
    // does not
    f.b.impl().sampleTreesFromPrior();
    f.b.impl().forestVariableCounts(0, 0, forestZero.data());
    f.b.impl().forestVariableCounts(0, 1, forestOne.data());
    check(totalSplits(forestOne) == 0 && totalSplits(forestZero) > 0,
          "facade setForestWeights: the weight lands on forest 1, not on "
          "forest 0");
    f.b.base().setForestWeights(1, nullptr);
    f.b.impl().sampleTreesFromPrior();
    Results results;
    f.b.impl().run(0, 1, results);  // the prior draw leaves the fits stale
  }},
  {FacadeVirtual::forestCalibration, "forestCalibration", [](Fixtures& f) {
    ForestCalibration viaBase = f.b.base().forestCalibration(1, 1);
    ForestCalibration viaImpl = f.b.impl().forestCalibration(1, 1);
    check(viaBase.priorScale == viaImpl.priorScale &&
            viaBase.nodeScaleFactor == viaImpl.nodeScaleFactor,
          "facade forestCalibration: the named chain and forest answer");
    check(viaBase.nodeScaleFactor !=
            f.b.base().forestCalibration(1, 0).nodeScaleFactor,
          "facade forestCalibration: the two forests report their own");
  }},
  {FacadeVirtual::setForestPriorScale, "setForestPriorScale", [](Fixtures& f) {
    check(!f.g.base().setForestPriorScale(1, 0.7) &&
            !f.b.base().setForestPriorScale(0, 0.7),
          "facade setForestPriorScale: an absent forest and a combiner-owned "
          "calibration refuse");
    check(f.g.base().setForestPriorScale(0, 0.7),
          "facade setForestPriorScale: the single forest takes it");
    checkNear(f.g.impl().forestCalibration(0, 0).priorScale, 0.7, 1.0e-12,
              "facade setForestPriorScale: the impl reports the new scale");
  }},
  {FacadeVirtual::setActiveRows, "setActiveRows", [](Fixtures& f) {
    std::vector<double> mask(Fixtures::n, 1.0);
    mask[0] = 0.0;
    std::vector<double> refused(Fixtures::n, 0.5);
    check(!f.l.base().setActiveRows(refused.data()),
          "facade setActiveRows: a non-binary mask is refused");
    check(f.l.base().setActiveRows(mask.data()),
          "facade setActiveRows: a binary mask installs");
    // an inactive row is out of the model, so the sweep draws no latent for it
    // where every active row's moves
    double masked = f.l.impl().latents(0)[0], active = f.l.impl().latents(0)[1];
    Results results;
    f.l.impl().run(0, 1, results);
    check(f.l.impl().latents(0)[0] == masked &&
            f.l.impl().latents(0)[1] != active,
          "facade setActiveRows: the mask reaches the impl's own sweep");
    check(f.l.base().setActiveRows(nullptr),
          "facade setActiveRows: and clears");
  }},
  {FacadeVirtual::setCounts, "setCounts", [](Fixtures& f) {
    check(!f.g.base().setCounts(f.counts0.data(), f.trials.data()),
          "facade setCounts: a sampler owning no counts refuses");
    check(f.m.base().setCounts(f.counts0.data(), f.trials.data()),
          "facade setCounts: the counts-owning coupling takes them");
    std::vector<double> fits(Fixtures::n * Fixtures::K, 0.0);
    Results results;
    results.trainingFits = fits.data();
    f.m.impl().run(20, 1, results);
    double first = 0.0, second = 0.0;
    for (std::size_t i = 0; i < Fixtures::n; ++i) {
      first += fits[i];
      second += fits[i + Fixtures::n];
    }
    check(first > 2.0 * second,
          "facade setCounts: the replacement counts move the fitted softmax");
  }},
  {FacadeVirtual::setCategoryOffset, "setCategoryOffset", [](Fixtures& f) {
    check(!f.g.base().setCategoryOffset(f.categoryOffset.data()),
          "facade setCategoryOffset: a sampler owning no counts refuses");
    for (std::size_t i = 0; i < Fixtures::n; ++i)
      f.categoryOffset[i + 2 * Fixtures::n] = 1.0e3;  // category 2 saturates
    check(f.m.base().setCategoryOffset(f.categoryOffset.data()),
          "facade setCategoryOffset: the coupling takes it");
    std::vector<double> fits(Fixtures::n * Fixtures::K, 0.0);
    Results results;
    results.trainingFits = fits.data();
    f.m.impl().run(0, 1, results);
    bool saturated = true;
    for (std::size_t i = 0; i < Fixtures::n; ++i)
      saturated &= fits[i + 2 * Fixtures::n] > 0.99;
    check(saturated,
          "facade setCategoryOffset: the offset enters the reported softmax");
    f.m.base().setCategoryOffset(nullptr);
  }},
  {FacadeVirtual::setCategoryTestOffset, "setCategoryTestOffset",
   [](Fixtures& f) {
    check(!f.g.base().setCategoryTestOffset(f.testCategoryOffset.data()),
          "facade setCategoryTestOffset: a sampler owning no counts refuses");
    for (std::size_t i = 0; i < Fixtures::nTest; ++i)
      f.testCategoryOffset[i + Fixtures::nTest] = 1.0e3;  // category 1
    check(f.m.base().setCategoryTestOffset(f.testCategoryOffset.data()),
          "facade setCategoryTestOffset: the coupling takes it");
    std::vector<double> testFits(Fixtures::nTest * Fixtures::K, 0.0);
    Results results;
    results.testFits = testFits.data();
    f.m.impl().run(0, 1, results);
    bool saturated = true;
    for (std::size_t i = 0; i < Fixtures::nTest; ++i)
      saturated &= testFits[i + Fixtures::nTest] > 0.99;
    check(saturated,
          "facade setCategoryTestOffset: the offset enters the reported test "
          "blend");
    f.m.base().setCategoryTestOffset(nullptr);
  }},
  {FacadeVirtual::totalAmplitudes, "totalAmplitudes", [](Fixtures& f) {
    check(f.b.base().totalAmplitudes() == f.b.impl().totalAmplitudes() &&
            f.b.base().totalAmplitudes() == 3 &&
            f.g.base().totalAmplitudes() == 0,
          "facade totalAmplitudes: the coupling's ragged total answers");
  }},
  {FacadeVirtual::numForestAmplitudes, "numForestAmplitudes", [](Fixtures& f) {
    check(f.b.base().numForestAmplitudes(1) ==
            f.b.impl().numForestAmplitudes(1) &&
            f.b.base().numForestAmplitudes(0) == 1 &&
            f.b.base().numForestAmplitudes(1) == 2,
          "facade numForestAmplitudes: each forest's own width answers");
  }},
  {FacadeVirtual::amplitudes, "amplitudes", [](Fixtures& f) {
    std::size_t width = f.b.impl().totalAmplitudes();
    std::vector<double> viaBase(width, 0.0), viaImpl(width, 1.0),
      otherChain(width, 2.0);
    check(f.b.base().amplitudes(1, viaBase.data()) &&
            f.b.impl().amplitudes(1, viaImpl.data()) &&
            f.b.base().amplitudes(0, otherChain.data()),
          "facade amplitudes: the coupling writes the channel");
    check(viaBase == viaImpl && viaBase != otherChain,
          "facade amplitudes: the named chain's block is written");
  }},
  {FacadeVirtual::forestTotalFits, "forestTotalFits", [](Fixtures& f) {
    std::vector<double> viaBase(Fixtures::n, 0.0), viaImpl(Fixtures::n, 1.0),
      otherForest(Fixtures::n, 2.0);
    f.b.base().forestTotalFits(1, 1, viaBase.data());
    f.b.impl().forestTotalFits(1, 1, viaImpl.data());
    f.b.base().forestTotalFits(1, 0, otherForest.data());
    check(viaBase == viaImpl && viaBase != otherForest,
          "facade forestTotalFits: the named chain and forest are written");
  }},
  {FacadeVirtual::fitsWithoutOffset, "fitsWithoutOffset", [](Fixtures& f) {
    std::vector<double> viaBase(Fixtures::n, 0.0), viaImpl(Fixtures::n, 1.0);
    check(f.g.base().fitsWithoutOffset(1, viaBase.data()) &&
            f.g.impl().fitsWithoutOffset(1, viaImpl.data()) &&
            viaBase == viaImpl,
          "facade fitsWithoutOffset: the named chain's location is written");
    check(!f.m.base().fitsWithoutOffset(0, viaBase.data()),
          "facade fitsWithoutOffset: a multi-location coupling refuses");
  }},
  {FacadeVirtual::forestVariableCounts, "forestVariableCounts",
   [](Fixtures& f) {
    std::vector<std::uint32_t> viaBase(Fixtures::p, 0), viaImpl(Fixtures::p, 1),
      otherForest(Fixtures::p, 2);
    f.b.base().forestVariableCounts(1, 1, viaBase.data());
    f.b.impl().forestVariableCounts(1, 1, viaImpl.data());
    f.b.base().forestVariableCounts(1, 0, otherForest.data());
    check(viaBase == viaImpl && viaBase != otherForest,
          "facade forestVariableCounts: the named chain and forest answer");
  }},
  {FacadeVirtual::numTreesInForest, "numTreesInForest", [](Fixtures& f) {
    check(f.b.base().numTreesInForest(1) == f.b.impl().numTreesInForest(1) &&
            f.b.base().numTreesInForest(0) == 6 &&
            f.b.base().numTreesInForest(1) == 4,
          "facade numTreesInForest: each forest's own count answers");
  }},
};

static_assert(std::size(rows) == numFacadeVirtuals,
              "the conformance table carries one row per SamplerBase virtual");

}  // namespace

void runFacadeTests() {
  // own runif01 stream, restored on the way out, so the suite reads the same
  // under a filter as it does in the full run
  std::uint64_t savedRngState = rngState;
  rngState = 0xd1b54a32d192ed03ull;

  Fixtures fixtures;
  for (int index = 0; index < static_cast<int>(std::size(rows)); ++index) {
    const Row& row = rows[index];
    callLog.currentRow = index;
    row.body(fixtures);
    bool reached = callLog.calledInRow(row.which, index);
    char what[128];
    std::snprintf(what, sizeof what,
                  "facade conformance: row %s calls the virtual it names",
                  row.name);
    check(reached, what);
  }

  // and no virtual is left to a duplicate row: every enumerator was reached
  bool covered = true;
  for (std::size_t i = 0; i < numFacadeVirtuals; ++i)
    covered &= callLog.everCalled(static_cast<FacadeVirtual>(i));
  check(covered, "facade conformance: every virtual is called by a row");

  rngState = savedRngState;
  printf("ok: facade conformance, %zu virtuals through the base pointer\n",
         numFacadeVirtuals);
}
