// SamplerShape oracle: the type-erased shape must agree with the typed impl
// on every field, over every construction path a concrete facade type can be
// held by. A swapped or forgotten fill line shows up here rather than as a
// mis-sized buffer in the bridge.

#include "common.hpp"

namespace {

/// One field, named in the failure message so a mismatch identifies itself.
void checkField(bool matches, const char* path, const char* field) {
  char what[128];
  std::snprintf(what, sizeof what, "shape (%s): %s", path, field);
  check(matches, what);
}

/// The oracle. Every SamplerShape field has a Sampler<L> counterpart of the
/// same name, so the comparison is one macro per field; a new field that
/// forgets its typed accessor fails to compile here.
#define CHECK_SHAPE_FIELD(field) \
  checkField(shape.field == impl.field(), path, #field)

template <typename L, typename ResidT>
void checkShapeMatchesImpl(const SamplerFacade<L, ResidT>& facade,
                           const char* path) {
  SamplerShape shape = facade.shape();
  const Sampler<L, ResidT>& impl(facade.impl());

  CHECK_SHAPE_FIELD(numObservations);
  CHECK_SHAPE_FIELD(numPredictors);
  CHECK_SHAPE_FIELD(numTestObservations);
  CHECK_SHAPE_FIELD(numChains);
  CHECK_SHAPE_FIELD(numThreads);
  CHECK_SHAPE_FIELD(numTrees);
  CHECK_SHAPE_FIELD(numForests);
  CHECK_SHAPE_FIELD(numGroups);
  CHECK_SHAPE_FIELD(numLeafCovariates);
  CHECK_SHAPE_FIELD(leafCovariateColumns);
  CHECK_SHAPE_FIELD(numReportedLocations);
  CHECK_SHAPE_FIELD(numVariableCountForests);
  CHECK_SHAPE_FIELD(numCutpoints);
  CHECK_SHAPE_FIELD(savedTreeCapacity);
  CHECK_SHAPE_FIELD(family);
  CHECK_SHAPE_FIELD(hasVarianceForest);
  CHECK_SHAPE_FIELD(usesFunctionLeaves);
  CHECK_SHAPE_FIELD(kIsSampled);
  CHECK_SHAPE_FIELD(usesDart);
  CHECK_SHAPE_FIELD(supportsResponseMutation);
  CHECK_SHAPE_FIELD(supportsForestWeights);
  CHECK_SHAPE_FIELD(supportsCountsMutation);
  CHECK_SHAPE_FIELD(testFitsAreDefined);
}

#undef CHECK_SHAPE_FIELD

/// A fixture whose columns all carry signal, so trees split and the shape is
/// read off a sampler that has actually run.
struct ShapeFixture {
  static constexpr size_t n = 150, p = 4;
  std::vector<double> x, y, xTest, z;
  std::vector<std::uint32_t> groups;

  ShapeFixture() : x(n * p), y(n), xTest(20 * p), z(n), groups(n) {
    for (double& v : x) v = runif01();
    for (double& v : xTest) v = runif01();
    for (size_t i = 0; i < n; ++i) {
      z[i] = runif01() < 0.5 ? 1.0 : 0.0;
      y[i] = std::sin(3.0 * x[i]) + x[i + n] + z[i] + 0.2 * runif01();
      groups[i] = static_cast<std::uint32_t>(i % 5);
    }
  }
};

/// A short run before the read, so the state-dependent fields
/// (savedTreeCapacity in particular) are exercised past their initial values.
template <typename L>
void runBriefly(SamplerFacade<L>& facade, size_t numTestObservations) {
  SamplerShape shape = facade.shape();
  const size_t numSamples = 4;
  size_t slab = shape.numReportedLocations * shape.numChains * numSamples;
  std::vector<double> trainingFits(shape.numObservations * slab);
  std::vector<double> testFits(
    numTestObservations * slab + 1);  // + 1 keeps a zero-row buffer addressable
  std::vector<std::uint32_t> variableCounts(
    shape.numPredictors * shape.numVariableCountForests * numSamples *
    shape.numChains);
  std::vector<double> sigma(numSamples * shape.numChains);

  Results results;
  results.sigma = sigma.data();
  results.trainingFits = trainingFits.data();
  results.testFits = numTestObservations > 0 ? testFits.data() : nullptr;
  results.variableCounts = variableCounts.data();
  results.numReportedLocations = shape.numReportedLocations;
  results.numVariableCountForests = shape.numVariableCountForests;
  facade.run(5, numSamples, results);
}

// The plain constant-leaf gaussian sampler, in the configuration that moves
// the most fields at once: several chains, saved trees, dart, a sampled k,
// grouped intercepts, and a test set.
void testConstantGaussian(ShapeFixture& fixture) {
  const size_t numChains = 2, numTest = 20;
  ext_rng* rngs[numChains];
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], 8100 + static_cast<std::uint32_t>(c));
  }

  SamplerOptions options;
  options.numTrees = 12;
  options.numChains = numChains;
  options.keepTrees = true;
  options.numSamplesToStore = 4;
  options.useDart = true;
  options.updateK = true;
  options.groupIndices = fixture.groups.data();
  options.numGroups = 5;

  SamplerFacade<ConstantGaussianLeaf> facade(
    fixture.x.data(), fixture.y.data(), ShapeFixture::n, ShapeFixture::p,
    nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542,
    options, rngs);
  facade.setTestPredictors(fixture.xTest.data(), numTest);

  SamplerShape shape = facade.shape();
  check(shape.numForests == 1, "shape: constant gaussian is single-forest");
  check(shape.numReportedLocations == 1,
        "shape: constant gaussian reports one location");
  check(shape.testFitsAreDefined,
        "shape: constant gaussian test fits are defined");
  check(!shape.supportsForestWeights,
        "shape: a single-forest sampler admits no per-forest weight");
  check(!shape.supportsCountsMutation,
        "shape: a single-forest sampler owns no count response");
  check(shape.numGroups == 5, "shape: grouped intercept count");
  check(shape.usesDart && shape.kIsSampled, "shape: dart and sampled k");
  check(!shape.usesFunctionLeaves, "shape: constant leaf is not function-valued");
  checkShapeMatchesImpl(facade, "constant gaussian, before run");

  runBriefly(facade, numTest);
  check(facade.shape().savedTreeCapacity == 4,
        "shape: saved-tree capacity after a run");
  checkShapeMatchesImpl(facade, "constant gaussian, after run");

  for (size_t c = 0; c < numChains; ++c) ext_rng_destroy(rngs[c]);
  printf("ok: shape matches the typed impl, constant gaussian\n");
}

// Vector- and function-parameter leaves, which are the only paths where the
// leaf covariate designation and usesFunctionLeaves are nonzero/true.
void testLeafCovariateSamplers(ShapeFixture& fixture) {
  const size_t covariateColumns[] = {2, 3};

  SamplerOptions options;
  options.numTrees = 8;
  options.leafCovariateColumns = covariateColumns;
  options.numLeafCovariates = 2;

  {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 8110);
    SamplerFacade<LinearGaussianLeaf> facade(
      fixture.x.data(), fixture.y.data(), ShapeFixture::n, ShapeFixture::p,
      nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0,
      0.37804942330213542, options, &rng);
    SamplerShape shape = facade.shape();
    check(shape.numLeafCovariates == 2, "shape: linear leaf covariate count");
    check(shape.leafCovariateColumns != nullptr &&
            shape.leafCovariateColumns[0] == 2 &&
            shape.leafCovariateColumns[1] == 3,
          "shape: linear leaf covariate columns");
    check(!shape.usesFunctionLeaves,
          "shape: linear leaf is not function-valued");
    checkShapeMatchesImpl(facade, "linear leaf");
    ext_rng_destroy(rng);
  }

  {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 8111);
    SamplerOptions gpOptions = options;
    gpOptions.gpLeaves = true;
    SamplerFacade<GPGaussianLeaf> facade(
      fixture.x.data(), fixture.y.data(), ShapeFixture::n, ShapeFixture::p,
      nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0,
      0.37804942330213542, gpOptions, &rng);
    check(facade.shape().usesFunctionLeaves,
          "shape: gp leaf is function-valued");
    checkShapeMatchesImpl(facade, "gp leaf");
    ext_rng_destroy(rng);
  }
  printf("ok: shape matches the typed impl, linear and gp leaves\n");
}

// The response families that move a shape field of their own: ordinal carries
// K - 1 cutpoints, and a variance forest flips hasVarianceForest.
void testResponseSurfaces(ShapeFixture& fixture) {
  {
    const size_t K = 4;
    std::vector<double> yOrdinal(ShapeFixture::n);
    for (size_t i = 0; i < ShapeFixture::n; ++i)
      yOrdinal[i] = static_cast<double>(1 + i % K);

    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 8120);
    SamplerOptions options;
    options.numTrees = 8;
    options.numCategories = K;
    SamplerFacade<ConstantGaussianLeaf> facade(
      fixture.x.data(), yOrdinal.data(), ShapeFixture::n, ShapeFixture::p,
      nullptr, nullptr, ResponseFamily::ordinal, 1.0, 3.0,
      0.37804942330213542, options, &rng);
    SamplerShape shape = facade.shape();
    check(shape.numCutpoints == K - 1, "shape: ordinal cutpoint count");
    check(shape.family == ResponseFamily::ordinal, "shape: ordinal family");
    checkShapeMatchesImpl(facade, "ordinal");
    ext_rng_destroy(rng);
  }

  {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 8121);
    SamplerOptions options;
    options.numTrees = 8;
    options.numVarianceTrees = 5;
    SamplerFacade<ConstantGaussianLeaf> facade(
      fixture.x.data(), fixture.y.data(), ShapeFixture::n, ShapeFixture::p,
      nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0,
      0.37804942330213542, options, &rng);
    check(facade.shape().hasVarianceForest,
          "shape: heteroscedastic variance forest");
    checkShapeMatchesImpl(facade, "heteroscedastic");
    ext_rng_destroy(rng);
  }
  printf("ok: shape matches the typed impl, ordinal and heteroscedastic\n");
}

// BCF: two forests, and the one path whose recorded test fits are undefined.
void testBCF(ShapeFixture& fixture) {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 8130);

  SamplerOptions options;
  BCFSpec spec;
  spec.mu.numTrees = 10;
  spec.tau.numTrees = 5;
  spec.z = fixture.z.data();

  SamplerFacade<ConstantGaussianLeaf> facade(
    fixture.x.data(), fixture.y.data(), ShapeFixture::n, ShapeFixture::p,
    nullptr, nullptr, 1.0, 3.0, 0.37804942330213542, options, spec, &rng);

  SamplerShape shape = facade.shape();
  check(shape.numForests == 2, "shape: bcf carries two forests");
  check(!shape.testFitsAreDefined, "shape: bcf test fits are undefined");
  // The condition the R bridge's and the flat C API's shared test-surface
  // guard fires on, spelled here because the guard itself lives in the bridge
  // and cannot link into this build. Positive half.
  check(shape.numForests >= 2 && !shape.testFitsAreDefined,
        "shape: bcf meets the test-surface refusal condition");
  check(shape.numTrees == spec.mu.numTrees,
        "shape: bcf numTrees addresses the prognostic forest");
  check(shape.numReportedLocations == 1, "shape: bcf reports one location");
  check(shape.supportsForestWeights,
        "shape: bcf admits a per-forest weight");
  // the counts probe is a capability, not a forest count: bcf carries two
  // forests and a gaussian response, so it must NOT answer the counts channel
  check(!shape.supportsCountsMutation,
        "shape: bcf owns no count response despite two forests");
  checkShapeMatchesImpl(facade, "bcf, before run");

  runBriefly(facade, 0);
  checkShapeMatchesImpl(facade, "bcf, after run");

  ext_rng_destroy(rng);
  printf("ok: shape matches the typed impl, bcf\n");
}

// Multinomial: the only path that widens the reported-location and
// variable-count-forest channels past 1.
void testMultinomial(ShapeFixture& fixture) {
  const size_t K = 3;
  std::vector<int> counts(ShapeFixture::n * K, 0), trials(ShapeFixture::n, 1);
  for (size_t i = 0; i < ShapeFixture::n; ++i)
    counts[i + (i % K) * ShapeFixture::n] = 1;

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 8140);

  SamplerOptions options;
  options.numTrees = 10;
  MultinomialSpec spec;
  spec.numCategories = K;
  spec.counts = counts.data();
  spec.trials = trials.data();
  spec.forest.numTrees = 10;

  SamplerFacade<ConstantGaussianLeaf> facade(fixture.x.data(), ShapeFixture::n,
                                             ShapeFixture::p, options, spec,
                                             &rng);

  SamplerShape shape = facade.shape();
  check(shape.numForests == K, "shape: multinomial carries K forests");
  check(shape.numReportedLocations == K,
        "shape: multinomial reports K locations");
  check(shape.numVariableCountForests == K,
        "shape: multinomial counts variables over K forests");
  check(shape.testFitsAreDefined,
        "shape: multinomial test fits are defined");
  // The negative half of the condition asserted for bcf above: a multi-forest
  // model whose test blend IS defined passes the refusal, which is why the
  // guard tests testFitsAreDefined and not the forest count - a numForests
  // test alone would refuse this sampler's whole test surface.
  check(!(shape.numForests >= 2 && !shape.testFitsAreDefined),
        "shape: multinomial escapes the test-surface refusal condition");
  check(!shape.supportsForestWeights,
        "shape: multinomial admits no per-forest weight despite K forests");
  // the counts channel's capability probe, and the pair the bridge reads with
  // it: K comes off numReportedLocations rather than a field of its own
  check(shape.supportsCountsMutation,
        "shape: multinomial owns a replaceable count response");
  check(!shape.supportsResponseMutation,
        "shape: multinomial's counts channel is not the response conduit");
  check(shape.numReportedLocations == K,
        "shape: the counts channel reads K off numReportedLocations");
  checkShapeMatchesImpl(facade, "multinomial, before run");

  runBriefly(facade, 0);
  checkShapeMatchesImpl(facade, "multinomial, after run");

  ext_rng_destroy(rng);
  printf("ok: shape matches the typed impl, multinomial\n");
}

}  // namespace

void runShapeTests(ext_rng*) {
  // SamplerShape rides Rf_error's longjmp in every bridge entry point, so an
  // owning field would leak on each error path. facade.hpp asserts this too;
  // repeat it here so a test build catches a regression on its own.
  static_assert(std::is_trivially_destructible<SamplerShape>::value,
                "SamplerShape must not own storage");

  ShapeFixture fixture;
  testConstantGaussian(fixture);
  testLeafCovariateSamplers(fixture);
  testResponseSurfaces(fixture);
  testBCF(fixture);
  testMultinomial(fixture);
}
