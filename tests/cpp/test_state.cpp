#include "common.hpp"

static void testFlattenRoundTrip() {
  // one categorical column (codes 0..3) and one ordinal, a hand-built tree
  // with one rule of each kind
  const size_t n = 12;
  std::vector<double> x(n * 2), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 4);
  for (size_t i = 0; i < n; ++i) x[i + n] = static_cast<double>(i) / (n - 1.0);
  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};

  ColumnStore store;
  store.build(x.data(), n, 2, 10, false, types);

  std::vector<index_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rootRule;
  rootRule.variableIndex = 0;
  rootRule.setCategoryDirections(0x6);  // categories 1, 2 right
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  Rule leftRule;
  leftRule.variableIndex = 1;
  leftRule.setSplitIndex(4);
  tree.birth(store, tree.at(0).leftChild, leftRule, y.data(), nullptr);

  std::vector<double> params(tree.nodes.size(), 0.0);
  std::vector<int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  for (size_t k = 0; k < bottoms.size(); ++k)
    params[static_cast<size_t>(bottoms[k])] = static_cast<double>(k + 1);

  std::vector<FlatNode> flat;
  std::vector<std::uint32_t> counts;
  tree.flatten(store, params.data(), flat, &counts);

  check(flat.size() == 5 && counts.size() == 5, "flatten emits pre-order");
  check(flat[0].variable == 0 &&
          flatKindOf(flat[0]) == FlatKind::categoricalInline &&
          flat[0].mask == 0x6,
        "flatten tags an inline categorical mask");
  check(flat[1].variable == 1 &&
          flatKindOf(flat[1]) == FlatKind::ordinal &&
          flat[1].value == store.cutPoints[1][4],
        "flatten tags an ordinal cut");
  check(flat[2].variable == invalidVariable && flat[2].value == 1.0 &&
        flat[3].value == 2.0 && flat[4].value == 3.0,
        "flatten stores leaf parameters left-first");
  check(counts[0] == n &&
        counts[1] == static_cast<std::uint32_t>(
                       tree.at(tree.at(0).leftChild).numObservations()) &&
        counts[2] + counts[3] == counts[1],
        "flatten counts mirror the partitions");

  // replay against the raw predictors reproduces the live counts
  std::vector<std::uint32_t> replayed(flat.size());
  std::vector<size_t> replayIndices(n);
  for (size_t i = 0; i < n; ++i) replayIndices[i] = i;
  countFlatObservationsBelow(flat.data(), x.data(), n,
                             replayIndices.data(), 0, n, replayed.data());
  check(replayed == counts, "flat replay reproduces the live counts");

  // per-row prediction accumulates each row's leaf parameter
  std::vector<double> fits(n, 0.0);
  for (size_t i = 0; i < n; ++i) replayIndices[i] = i;
  addFlatPredictionsBelow(flat.data(), x.data(), n,
                          replayIndices.data(), 0, n, fits.data());
  bool fitsMatch = true;
  for (size_t i = 0; i < n; ++i) {
    int32_t leaf = tree.findBottomNodeForObservation(store, i);
    fitsMatch &= fits[i] == params[static_cast<size_t>(leaf)];
  }
  check(fitsMatch, "flat prediction routes rows to their leaves");

  // rebuild recovers the rules exactly
  std::vector<index_t> indices2(n);
  Tree tree2;
  tree2.initialize(indices2.data(), n);
  std::vector<double> params2;
  check(tree2.buildFromFlat(store, flat.data(), flat.size(), params2),
        "buildFromFlat accepts its own flatten");
  check(tree2.at(0).rule.variableIndex == 0 &&
        tree2.at(0).rule.categoryDirections() == 0x6,
        "buildFromFlat recovers the categorical mask");
  check(tree2.at(tree2.at(0).leftChild).rule.splitIndex() == 4,
        "buildFromFlat recovers the ordinal split index exactly");
  tree2.repartitionSubtree(store, 0);
  bool partitionsMatch = true;
  for (size_t i = 0; i < tree.nodes.size(); ++i)
    partitionsMatch &= tree2.at(static_cast<int32_t>(i)).numObservations() ==
                       tree.at(static_cast<int32_t>(i)).numObservations();
  check(partitionsMatch, "rebuilt tree partitions identically");
  bool paramsMatch = true;
  for (int32_t i : bottoms)
    paramsMatch &= params2[static_cast<size_t>(i)] ==
                   params[static_cast<size_t>(i)];
  check(paramsMatch, "buildFromFlat recovers leaf parameters");

  // malformed inputs: an ordinal value off the cut grid, a mask outside the
  // canonical gauge, and a truncated pre-order all refuse
  std::vector<FlatNode> bad(flat);
  bad[1].value += 1.0e-3;
  tree2.initialize(indices2.data(), n);
  check(!tree2.buildFromFlat(store, bad.data(), bad.size(), params2),
        "buildFromFlat rejects a value off the cut grid");
  bad = flat;
  bad[1].variable = 0;
  setFlatKind(bad[1], FlatKind::categoricalInline);
  bad[1].mask = 0x2;  // category 1 goes right, but 1 is unreachable here
  tree2.initialize(indices2.data(), n);
  check(!tree2.buildFromFlat(store, bad.data(), bad.size(), params2),
        "buildFromFlat rejects an out-of-gauge mask");
  tree2.initialize(indices2.data(), n);
  check(!tree2.buildFromFlat(store, flat.data(), flat.size() - 1, params2),
        "buildFromFlat rejects a truncated pre-order");
  check(flatTreeIsWellFormed(store, flat.data(), flat.size()) &&
        !flatTreeIsWellFormed(store, flat.data(), flat.size() - 1),
        "well-formedness check matches");

  printf("ok: flatten round trip\n");
}

static void testCategoricalFlattenBoundaries() {
  // the inline (<= 63) / pooled (>= 64) edge: for K at and around 63/64 a
  // flatten -> replay must reproduce the live routing, and a rebuild must
  // restore partitions identically
  const size_t Ks[] = {53, 54, 63, 64};
  for (size_t K : Ks) {
    const size_t n = 8 * K;
    std::vector<double> x(n), y(n, 0.0);
    for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % K);
    ColumnType types[] = {ColumnType::categorical};
    ColumnStore store;
    store.build(x.data(), n, 1, 10, false, types);
    check(store.numCuts[0] == K && store.columnIsPooled(0) == (K >= 64),
          "the pooling boundary is 64 categories");

    std::vector<index_t> indices(n);
    Tree tree;
    tree.initialize(indices.data(), n);
    Rule rule;
    rule.variableIndex = 0;
    // send a low category and the top one right, straddling a word top
    if (store.columnIsPooled(0)) {
      size_t offset =
        tree.allocateMask(maskWordsForCount(static_cast<std::uint32_t>(K)));
      std::uint64_t* words = tree.mutableMaskWordsFor(offset);
      maskSetBit(words, 2);
      maskSetBit(words, static_cast<std::uint32_t>(K - 1));
      rule.setMaskOffset(offset);
    } else {
      rule.setCategoryDirections((1ull << 2) | (1ull << (K - 1)));
    }
    tree.birth(store, 0, rule, y.data(), nullptr);

    std::vector<double> params(tree.nodes.size(), 0.0);
    std::vector<FlatNode> flat;
    std::vector<std::uint32_t> counts;
    std::vector<std::uint64_t> masks;
    tree.flatten(store, params.data(), flat, &counts, 1, nullptr, &masks);
    check(store.columnIsPooled(0) ==
            (flatKindOf(flat[0]) == FlatKind::categoricalPooled),
          "the tag matches the pooling tier");

    std::vector<std::uint32_t> replayed(flat.size());
    std::vector<size_t> replayIndices(n);
    for (size_t i = 0; i < n; ++i) replayIndices[i] = i;
    countFlatObservationsBelow(flat.data(), x.data(), n, replayIndices.data(),
                               0, n, replayed.data(), masks.data());
    check(replayed == counts, "flatten -> replay reproduces live routing");
    check(flatTreeIsWellFormed(store, flat.data(), flat.size(), masks.data(),
                               masks.size()),
          "boundary flat tree is well formed");

    Tree rebuilt;
    std::vector<index_t> rebuiltIndices(n);
    rebuilt.initialize(rebuiltIndices.data(), n);
    std::vector<double> rebuiltParams;
    check(rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams,
                                1, nullptr, masks.data(), masks.size()),
          "boundary rule rebuilds from flat");
    rebuilt.repartitionSubtree(store, 0);
    bool partsMatch = tree.nodes.size() == rebuilt.nodes.size();
    for (size_t i = 0; partsMatch && i < tree.nodes.size(); ++i)
      partsMatch &= rebuilt.at(static_cast<int32_t>(i)).numObservations() ==
                    tree.at(static_cast<int32_t>(i)).numObservations();
    check(partsMatch, "boundary rebuilt tree partitions identically");
  }
  printf("ok: categorical flatten boundaries\n");
}

static void testKeepTrees(ext_rng* rng) {
  const size_t n = 200, nTest = 20, numSamples = 4;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = runif01();

  SamplerOptions options;
  options.numTrees = 25;
  options.keepTrees = true;
  options.numSamplesToStore = numSamples;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  sampler.setTestPredictors(xTest.data(), nTest);

  Results empty;
  sampler.run(100, 0, empty);
  check(sampler.currentSampleNum() == 0, "burn-in does not advance the slot");

  std::vector<double> sigma(numSamples), testFits(nTest * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.testFits = testFits.data();
  sampler.run(0, numSamples, results);
  check(sampler.currentSampleNum() == 0, "a full run wraps the slot");

  // replaying the saved trees against the raw test rows must reproduce the
  // run's recorded test fits exactly: same parameters, same addition order
  std::vector<double> predicted(nTest * numSamples);
  sampler.predict(xTest.data(), nTest, predicted.data());
  check(predicted == testFits, "saved-tree predictions equal the run's test fits");

  // a partial run overwrites the oldest slots and leaves the rest
  std::vector<double> sigma2(2), testFits2(nTest * 2);
  Results results2;
  results2.sigma = sigma2.data();
  results2.testFits = testFits2.data();
  sampler.run(0, 2, results2);
  check(sampler.currentSampleNum() == 2, "a partial run advances the slot");

  std::vector<double> predicted2(nTest * numSamples);
  sampler.predict(xTest.data(), nTest, predicted2.data());
  bool overwritten = std::equal(testFits2.begin(), testFits2.end(),
                                predicted2.begin());
  bool preserved = std::equal(predicted.begin() + nTest * 2, predicted.end(),
                              predicted2.begin() + nTest * 2);
  check(overwritten, "new samples overwrite the oldest slots");
  check(preserved, "later slots survive a partial run");
  check(sampler.savedTreeCapacity() == numSamples,
        "keepTrees capacity comes from the options");

  printf("ok: keepTrees\n");
}

static void testPredictCurrentTrees(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  check(sampler.savedTreeCapacity() == 0, "keepTrees defaults off");

  // routing the training rows through the live trees agrees with the run's
  // recorded training fits (both are scale * sum-of-tree-fits + shift; only
  // the accumulation order differs)
  std::vector<double> trainingFits(n);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(0, 1, results);

  std::vector<double> predicted(n);
  sampler.predict(x.data(), n, predicted.data());
  bool match = true;
  for (size_t i = 0; i < n; ++i)
    match &= std::fabs(predicted[i] - trainingFits[i]) <= 1.0e-8;
  check(match, "live-tree predictions equal the recorded training fits");

  printf("ok: predict from current trees\n");
}

static void testStateRoundTripScaledOffset() {
  // setOffset(updateScale) moves the gaussian response transform after
  // creation; the state must carry it or a restored sampler mis-scales
  // every internal quantity
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> offset(n);
  for (size_t i = 0; i < n; ++i)
    offset[i] = 2.0 * std::sin(0.1 * (double) i);  // widens y - offset

  SamplerOptions options;
  options.numTrees = 25;

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 77) != 0 ||
      ext_rng_setSeed(rngB, 78) != 0) {
    check(false, "scaled state: rng creation");
    return;
  }

  ConstantLeafSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngA);
  Results empty;
  original.run(20, 0, empty);
  original.setOffset(offset.data(), true);
  original.run(20, 0, empty);

  SamplerStateData state;
  original.getState(state);
  check(state.chains[0].fitMax > state.chains[0].fitMin,
        "gaussian state captures the response transform");

  ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngB);
  // a host reinstalls the current offset but cannot reproduce the scale
  // trajectory that produced the state; the stored transform must win
  restored.setOffset(offset.data(), false);
  check(restored.setState(state, nullptr), "scaled state restores");

  // the moved transform round-trips: the restored model matches the source,
  // and its live-tree predictions land on the original scale, both before
  // either chain continues past the save point
  checkStructuralRoundTrip(state, restored,
                           "scaled restore reproduces the moved-scale model");

  std::vector<double> xTest(20 * 2);
  for (double& v : xTest) v = runif01();
  std::vector<double> predictionsA(20), predictionsB(20);
  original.predict(xTest.data(), 20, predictionsA.data());
  restored.predict(xTest.data(), 20, predictionsB.data());
  check(predictionsA == predictionsB,
        "scaled restore predicts on the original scale");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: state round-trip with a moved scale\n");
}

static void testStateRoundTrip() {
  // store the state, restore it into a fresh sampler, and gate the round trip
  // structurally (gate a) and by continued-vs-uninterrupted agreement (gate b)
  const size_t n = 200, numChains = 2;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t seed) {
    for (size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], seed + static_cast<std::uint32_t>(c));
    }
  };

  SamplerOptions options;
  options.numTrees = 25;
  options.numChains = numChains;
  options.keepTrees = true;
  options.numSamplesToStore = 3;

  std::vector<ext_rng*> rngs(numChains, nullptr);
  makeRngs(rngs, 1001);
  ConstantLeafSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(60, 2, empty);  // leave a nonzero slot position behind

  SamplerStateData state;
  original.getState(state);
  check(state.chains.size() == numChains && state.currentSampleNum == 2,
        "state captures every chain and the slot position");

  // different seeds on purpose: the serialized rng state must win
  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 9999);
  ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state, nullptr), "a stored state restores");

  checkStructuralRoundTrip(state, restored,
                           "restored state reproduces the saved model");

  // saved trees round-tripped: predictions from the copied buffer agree
  // exactly, before either chain continues past the save point
  std::vector<double> xTest(20 * 2);
  for (double& v : xTest) v = runif01();
  size_t capacity = original.savedTreeCapacity();
  std::vector<double> predictA(20 * capacity * numChains);
  std::vector<double> predictB(20 * capacity * numChains);
  original.predict(xTest.data(), 20, predictA.data());
  restored.predict(xTest.data(), 20, predictB.data());
  check(predictA == predictB, "saved trees ride along with the state");

  // gate (b): draws diverge in the last ulp and then chaotically, but a
  // continued chain and the uninterrupted original track the same posterior
  // well inside Monte Carlo error
  const size_t window = 300, total = window * numChains;
  std::vector<double> sigmaA(total), sigmaB(total);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsB.sigma = sigmaB.data();
  original.run(0, window, resultsA);
  restored.run(0, window, resultsB);
  double sumA = 0.0, sumB = 0.0, sumSqA = 0.0;
  for (size_t i = 0; i < total; ++i) {
    sumA += sigmaA[i];
    sumB += sigmaB[i];
    sumSqA += sigmaA[i] * sigmaA[i];
  }
  double meanA = sumA / total, meanB = sumB / total;
  double mcse = std::sqrt((sumSqA / total - meanA * meanA) / total);
  check(std::fabs(meanA - meanB) < 8.0 * mcse,
        "restored chain continues the same posterior");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: state round trip\n");
}

static void testStateRoundTripLatents(ext_rng* rng) {
  // logistic Polya-Gamma latents and dart state must ride along too
  const size_t n = 150;
  std::vector<double> x(n * 2), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = (x[i] + 0.25 * (runif01() - 0.5) > 0.5) ? 1.0 : 0.0;

  SamplerOptions options;
  options.numTrees = 10;
  options.nodeScale = 3.0;
  options.useDart = true;
  options.dart.updateDelay = 10;

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngA, 555);
  ConstantLeafSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::logistic, 1.0, 3.0, 1.0, options,
                          &rngA);
  Results empty;
  original.run(40, 0, empty);

  SamplerStateData state;
  original.getState(state);
  check(!state.chains[0].latents.empty(), "logistic state carries omega");
  check(!state.chains[0].dartProbabilities.empty(), "dart state is captured");

  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngB, 777);
  ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::logistic, 1.0, 3.0, 1.0, options,
                          &rngB);
  check(restored.setState(state, nullptr), "a binary+dart state restores");

  checkStructuralRoundTrip(state, restored,
                           "restored logistic + dart state reproduces the model");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  (void) rng;
  printf("ok: state round trip with latents\n");
}

static void testStateRoundTripStudentT(ext_rng* /*rng*/) {
  // Student-t continuous errors: the mixing precisions lambda (in the latents
  // block) and the residual df nu (its companion scalar) must ride the state,
  // restore exactly, and let a fresh sampler continue the same posterior; both
  // fixed-nu and estimated-nu modes. RNG-insulated per testLeafOfConsistency
  // (local generator + restored rngState): a stray runif01 would shift every
  // downstream snapshot.
  std::uint64_t savedRngState = rngState;
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  auto runOneMode = [&](double residualDf, bool estimated, const char* tag) {
    SamplerOptions options;
    options.numTrees = 25;
    options.residualDf = residualDf;

    ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngA, 20240717u);
    ext_rng_setSeed(rngB, 90909u);  // different: the serialized rng must win
    ConstantLeafSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                            ResponseFamily::gaussian, 1.0, 3.0,
                            0.37804942330213542, options, &rngA);
    Results empty;
    original.run(40, 0, empty);

    SamplerStateData state;
    original.getState(state);
    check(state.chains[0].latents.size() == n,
          "t state carries the mixing precisions");
    check(state.chains[0].residualDf > 0.0,
          "t state carries a positive residual df");
    if (!estimated)
      check(state.chains[0].residualDf == residualDf,
            "fixed-nu state records the exact df");

    ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                            ResponseFamily::gaussian, 1.0, 3.0,
                            0.37804942330213542, options, &rngB);
    check(restored.setState(state, nullptr), "a t state restores");

    // gate (a): the restored model reproduces the saved one bitwise, lambda
    // and nu included (statesAgree compares both)
    checkStructuralRoundTrip(state, restored,
                             "restored t state reproduces the model");
    SamplerStateData reState;
    restored.getState(reState);
    check(reState.chains[0].residualDf == state.chains[0].residualDf,
          "nu round-trips exactly");
    check(reState.chains[0].latents == state.chains[0].latents,
          "lambda round-trips exactly");

    // gate (b): draws diverge in the last ulp of the re-accumulated fits and
    // then chaotically, but the restored chain tracks the same sigma posterior
    // well inside Monte Carlo error
    const size_t window = 500;
    std::vector<double> sigmaA(window), sigmaB(window);
    Results resultsA, resultsB;
    resultsA.sigma = sigmaA.data();
    resultsB.sigma = sigmaB.data();
    original.run(0, window, resultsA);
    restored.run(0, window, resultsB);
    double sumA = 0.0, sumB = 0.0, sumSqA = 0.0;
    for (size_t i = 0; i < window; ++i) {
      sumA += sigmaA[i];
      sumB += sigmaB[i];
      sumSqA += sigmaA[i] * sigmaA[i];
    }
    double meanA = sumA / window, meanB = sumB / window;
    double mcse = std::sqrt((sumSqA / window - meanA * meanA) / window);
    check(std::fabs(meanA - meanB) < 8.0 * mcse,
          "restored t chain continues the same posterior");

    ext_rng_destroy(rngB);
    ext_rng_destroy(rngA);
    printf("ok: state round trip, student-t (%s)\n", tag);
  };

  runOneMode(4.0, false, "fixed nu");
  runOneMode(0.0, true, "estimated nu");

  // a t sampler refuses a state without the t blocks: a gaussian state carries
  // neither latents nor nu, so the shared latents check and the nu requirement
  // together reject it
  SamplerOptions gaussOptions;
  gaussOptions.numTrees = 25;
  ext_rng* rngG = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngG, 4242u);
  ConstantLeafSampler gaussian(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, gaussOptions, &rngG);
  Results emptyG;
  gaussian.run(20, 0, emptyG);
  SamplerStateData gaussState;
  gaussian.getState(gaussState);
  check(gaussState.chains[0].latents.empty(),
        "gaussian state carries no latents");
  check(std::isnan(gaussState.chains[0].residualDf),
        "gaussian state carries no residual df");

  SamplerOptions tOptions;
  tOptions.numTrees = 25;
  tOptions.residualDf = 4.0;
  ext_rng* rngT = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngT, 5353u);
  ConstantLeafSampler tSampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, tOptions, &rngT);
  check(!tSampler.setState(gaussState, nullptr),
        "a t sampler refuses a gaussian state");

  ext_rng_destroy(rngT);
  ext_rng_destroy(rngG);
  rngState = savedRngState;
  printf("ok: state round trip, student-t rejects a gaussian state\n");
}

static void testStateValidation(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  SamplerStateData state;
  sampler.getState(state);
  std::vector<double> treeFitsBefore(sampler.chain(0).treeFits());
  std::vector<std::vector<double>> cutsBefore(sampler.data().cutPoints);

  SamplerStateData bad(state);
  bad.chains.pop_back();
  check(!sampler.setState(bad, nullptr), "setState rejects a chain-count mismatch");

  bad = state;
  bool corrupted = false;
  for (std::vector<FlatNode>& tree : bad.chains[0].forests[0].trees) {
    if (tree.size() > 1) {
      tree[0].value += 1.0e-3;  // off the cut grid
      corrupted = true;
      break;
    }
  }
  check(corrupted && !sampler.setState(bad, nullptr),
        "setState rejects a split value off the cut grid");
  check(sampler.data().cutPoints == cutsBefore,
        "a rejected state leaves the cuts untouched");
  check(sampler.chain(0).treeFits() == treeFitsBefore,
        "a rejected state leaves the fits untouched");

  check(sampler.setState(state, nullptr), "the original state still restores");
  check(sampler.chain(0).treeFits() == treeFitsBefore,
        "restoring the current state is an identity");

  printf("ok: state validation\n");
}

// ---------------------------------------------------------------------------
// Interaction containment (docs/design/interaction-constraints.md,
// "Containment"): a state install or warm start must not admit a tree that
// violates a forest's interaction constraint - the availability predicate is
// not self-checking (splitVariableLogProbability never re-tests the node's own
// variable), so treeLogProbability would silently mis-score a donor grown
// unconstrained. Grow an UNCONSTRAINED donor whose pure-interaction signal
// forces order-2 paths, prove by hand that at least one of its trees violates
// max.order = 1, then assert both setState and installForests refuse it -
// while a same-constraint donor is accepted, so the gate is specific.
// ---------------------------------------------------------------------------
static void testInteractionContainment() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20260721);
  std::uint64_t savedRngState = rngState;
  rngState = 424242u;

  const size_t n = 300, p = 2, numTrees = 30;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();          // x0
    x[i + n] = runif01();      // x1
    // a pure AND-interaction: large only where BOTH exceed 0.5, unrepresentable
    // by a sum of single-variable steps, so a fitting tree MUST split on x0 and
    // x1 on one path (order 2) - impossible under max.order = 1 or forbid(x0,x1)
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = ((x[i] > 0.5 && x[i + n] > 0.5) ? 4.0 : -1.0) + 0.1 * z;
  }

  // borrow-lifetime bookkeeping: each Sampler holds its rng pointer, so the
  // rng objects must outlive their samplers
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  };
  auto makeSampler = [&](std::uint32_t seed, size_t maxOrder,
                         const size_t* pair) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.interactionMaxOrder = maxOrder;
    if (pair != nullptr) {
      options.interactionForbiddenPairs = pair;
      options.interactionNumForbiddenPairs = 1;
    }
    ext_rng* r = newRng(seed);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, options, &r);
  };

  // (1) the unconstrained donor, grown until its interaction is captured
  auto donor = makeSampler(555, 0, nullptr);
  Results empty;
  donor->run(200, 0, empty);
  SamplerStateData donorState;
  donor->getState(donorState);

  // by-hand proof the donor really is infeasible under max.order = 1: build
  // each donor tree against an independent store carrying a K = 1 constraint
  ColumnStore store;
  store.build(x.data(), n, p, 100);  // the sampler's default cut grid
  InteractionConstraint k1;
  k1.build(p, 1, nullptr, 0);
  std::vector<index_t> idx(n);
  std::vector<double> params;
  Tree scratch;
  size_t violators = 0;
  for (const std::vector<FlatNode>& flat : donorState.chains[0].forests[0].trees) {
    scratch.initialize(idx.data(), n);
    scratch.setInteractionConstraint(&k1);
    if (scratch.buildFromFlat(store, flat.data(), flat.size(), params) &&
        !scratch.interactionSubtreeIsValid(0))
      ++violators;
  }
  scratch.setInteractionConstraint(nullptr);
  check(violators > 0,
        "containment: the unconstrained donor holds an order-2 (infeasible) tree");

  // (2) a max.order = 1 target refuses the donor on both install paths
  auto k1Target = makeSampler(777, 1, nullptr);
  check(!k1Target->setState(donorState, nullptr),
        "containment: setState refuses an interaction-violating donor");
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(k1Target->installForests(donorState, liveMap) ==
          WarmStartResult::interactionMismatch,
        "containment: warm start refuses an interaction-violating donor");

  // (3) a forbid(x0, x1) target refuses it likewise (a distinct predicate)
  size_t pair[] = {0, 1};
  auto forbidTarget = makeSampler(888, 0, pair);
  check(!forbidTarget->setState(donorState, nullptr),
        "containment: setState refuses a forbidden-pair violation");
  check(forbidTarget->installForests(donorState, liveMap) ==
          WarmStartResult::interactionMismatch,
        "containment: warm start refuses a forbidden-pair violation");

  // (4) specificity: a same-constraint donor's trees are feasible and install
  auto k1Donor = makeSampler(999, 1, nullptr);
  k1Donor->run(200, 0, empty);
  SamplerStateData k1DonorState;
  k1Donor->getState(k1DonorState);
  auto k1Target2 = makeSampler(111, 1, nullptr);
  check(k1Target2->installForests(k1DonorState, liveMap) == WarmStartResult::ok,
        "containment: a same-constraint donor warm-starts cleanly");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: interaction containment (%zu donor violators)\n", violators);
}

// Column-mask containment (F1): like the interaction gate above, the warm-start
// install path had no feasibility gate for a forest's per-forest column mask
// (BCF moderators, a column-restricted variance forest). A donor tree splitting
// on a column OUTSIDE the target forest's mask was rebuilt and accepted, then
// splitVariableLogProbability scored it over an availability menu
// (collectAvailableVariables) that excludes the split - a silent wrong tree-prior
// density. Grow an UNRESTRICTED donor whose signal forces an x1 split, prove by
// hand its trees violate a column-0-only mask, then assert both setState and
// installForests refuse it while a same-mask donor is accepted (a specific gate).
// ---------------------------------------------------------------------------
static void testColumnMaskContainment() {
  std::uint64_t savedRngState = rngState;
  rngState = 313131u;

  const size_t n = 300, p = 2, numTrees = 30;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();          // x0
    x[i + n] = runif01();      // x1
    // signal lives ONLY in x1, so an unrestricted fitting tree must split on it -
    // a split a column-0-only forest forbids
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = 3.0 * x[i + n] + 0.1 * z;
  }

  // each Sampler holds its rng pointer, so the rng objects must outlive it
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  };
  // the allow list outlives construction (consumed into forest.columnMask there)
  const std::vector<size_t> allowZero = {0};
  auto makeSampler = [&](std::uint32_t seed, bool restrictToZero) {
    SamplerOptions options;
    options.numTrees = numTrees;
    if (restrictToZero) {
      options.forestColumns = allowZero.data();
      options.numForestColumns = allowZero.size();
    }
    ext_rng* r = newRng(seed);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, options, &r);
  };

  // (1) the unrestricted donor, grown until it captures the x1 signal
  auto donor = makeSampler(555, false);
  Results empty;
  donor->run(200, 0, empty);
  SamplerStateData donorState;
  donor->getState(donorState);

  // by-hand proof the donor is infeasible under a column-0-only mask: build each
  // donor tree against an independent store carrying that mask
  ColumnStore store;
  store.build(x.data(), n, p, 100);  // the sampler's default cut grid
  std::vector<std::uint8_t> maskZeroOnly(p, 0);
  maskZeroOnly[0] = 1;
  std::vector<index_t> idx(n);
  std::vector<double> params;
  Tree scratch;
  size_t violators = 0;
  for (const std::vector<FlatNode>& flat : donorState.chains[0].forests[0].trees) {
    scratch.initialize(idx.data(), n);
    scratch.setColumnMask(maskZeroOnly.data());
    if (scratch.buildFromFlat(store, flat.data(), flat.size(), params) &&
        !scratch.columnMaskSubtreeIsValid(0))
      ++violators;
  }
  scratch.setColumnMask(nullptr);
  check(violators > 0,
        "column mask: the unrestricted donor holds an out-of-mask tree");

  // (2) a column-0-only target refuses the donor on both install paths
  auto target = makeSampler(777, true);
  check(!target->setState(donorState, nullptr),
        "column mask: setState refuses an out-of-mask donor");
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(target->installForests(donorState, liveMap) ==
          WarmStartResult::columnMaskMismatch,
        "column mask: warm start refuses an out-of-mask donor");

  // (3) specificity: a same-mask donor's trees are feasible and install cleanly
  auto compliantDonor = makeSampler(999, true);
  compliantDonor->run(200, 0, empty);
  SamplerStateData compliantState;
  compliantDonor->getState(compliantState);
  auto target2 = makeSampler(111, true);
  check(target2->installForests(compliantState, liveMap) == WarmStartResult::ok,
        "column mask: a same-mask donor warm-starts cleanly");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: column-mask containment (%zu donor violators)\n", violators);
}

// Block-additive constraint (variant A, docs/design/interaction-constraints.md):
// each whole tree is confined to one declared group of predictors via the static
// per-tree column mask, so the ensemble is exactly f = sum_G f_G. Two groups
// {x0} / {x1} split the 30 trees evenly (15 each, deterministic contiguous
// assignment). Assert (a) every tree of a block-built forest splits only within
// its group, and (b) the shipped columnMask containment gate (F1) refuses a
// warm-start / setState donor whose trees split outside their block - the block
// masks lower onto the same per-tree columnMask_, so no second gate is needed.
// RNG-neutral (saves/restores rngState + a local seeded ext_rng).
// ---------------------------------------------------------------------------
static void testBlockAdditiveConfinement() {
  std::uint64_t savedRngState = rngState;
  rngState = 606060u;

  const size_t n = 300, p = 2, numTrees = 30, half = numTrees / 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();          // x0
    x[i + n] = runif01();      // x1
    // both columns carry signal, so group-0 trees split on x0 and group-1 trees
    // on x1 - confinement is then non-trivially exercised on both blocks
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = 2.0 * x[i] + 3.0 * x[i + n] + 0.1 * z;
  }

  // group 0 = {col 0}, group 1 = {col 1}; the trees split evenly and
  // contiguously, so tree t belongs to group (t < half ? 0 : 1)
  const std::vector<std::int32_t> blockOfColumn = {0, 1};
  const std::vector<size_t> blockTreeCounts = {half, numTrees - half};
  auto groupOfTree = [&](size_t t) { return t < half ? 0u : 1u; };
  // a one-column allow mask per group, for the by-hand feasibility checks
  std::vector<std::uint8_t> maskGroup[2] = {{1, 0}, {0, 1}};

  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  };
  auto makeSampler = [&](std::uint32_t seed, bool blocked) {
    SamplerOptions options;
    options.numTrees = numTrees;
    if (blocked) {
      options.numBlocks = blockTreeCounts.size();
      options.blockOfColumn = blockOfColumn.data();
      options.blockTreeCounts = blockTreeCounts.data();
    }
    ext_rng* r = newRng(seed);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, options, &r);
  };

  ColumnStore store;
  store.build(x.data(), n, p, 100);  // the sampler's default cut grid
  std::vector<index_t> idx(n);
  std::vector<double> params;
  Tree scratch;

  // (a) a block-built forest confines every tree to its group's columns
  auto blocked = makeSampler(555, true);
  Results empty;
  blocked->run(200, 0, empty);
  SamplerStateData blockedState;
  blocked->getState(blockedState);
  const auto& blockedTrees = blockedState.chains[0].forests[0].trees;
  bool allConfined = true;
  size_t splitters = 0;
  for (size_t t = 0; t < blockedTrees.size(); ++t) {
    scratch.initialize(idx.data(), n);
    scratch.setColumnMask(maskGroup[groupOfTree(t)].data());
    if (!scratch.buildFromFlat(store, blockedTrees[t].data(),
                               blockedTrees[t].size(), params))
      continue;
    if (!scratch.at(0).isBottom()) ++splitters;  // a tree that actually split
    if (!scratch.columnMaskSubtreeIsValid(0)) allConfined = false;
  }
  scratch.setColumnMask(nullptr);
  check(allConfined, "blocks: every tree splits only within its group");
  check(splitters > 0, "blocks: some trees actually split (non-trivial)");

  // (b) an UNRESTRICTED donor's trees split on both columns, so a donor tree in
  // a group-0 slot that split on x1 (or vice versa) violates its block
  auto donor = makeSampler(777, false);
  donor->run(200, 0, empty);
  SamplerStateData donorState;
  donor->getState(donorState);
  const auto& donorTrees = donorState.chains[0].forests[0].trees;
  size_t violators = 0;
  for (size_t t = 0; t < donorTrees.size(); ++t) {
    scratch.initialize(idx.data(), n);
    scratch.setColumnMask(maskGroup[groupOfTree(t)].data());
    if (scratch.buildFromFlat(store, donorTrees[t].data(), donorTrees[t].size(),
                              params) &&
        !scratch.columnMaskSubtreeIsValid(0))
      ++violators;
  }
  scratch.setColumnMask(nullptr);
  check(violators > 0, "blocks: the unrestricted donor holds out-of-block trees");

  // the shipped gate (F1) refuses the out-of-block donor on both install paths:
  // each live tree's columnMask_ is its block row, so rebuildLiveForest's
  // columnMaskSubtreeIsValid catches the violation. No second gate added.
  auto target = makeSampler(999, true);
  check(!target->setState(donorState, nullptr),
        "blocks: setState refuses an out-of-block donor");
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(target->installForests(donorState, liveMap) != WarmStartResult::ok,
        "blocks: warm start refuses an out-of-block donor");

  // (c) specificity: a same-blocks donor's trees are feasible and install cleanly
  auto compliantDonor = makeSampler(1234, true);
  compliantDonor->run(200, 0, empty);
  SamplerStateData compliantState;
  compliantDonor->getState(compliantState);
  auto target2 = makeSampler(4321, true);
  check(target2->installForests(compliantState, liveMap) == WarmStartResult::ok,
        "blocks: a same-blocks donor warm-starts cleanly");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: block-additive confinement (%zu splitters, %zu donor violators)\n",
         splitters, violators);
}

void runStateTests(ext_rng* rng) {
  testFlattenRoundTrip();
  testCategoricalFlattenBoundaries();
  testKeepTrees(rng);
  testPredictCurrentTrees(rng);
  testStateRoundTrip();
  testStateRoundTripScaledOffset();
  testStateRoundTripLatents(rng);
  testStateRoundTripStudentT(rng);
  testStateValidation(rng);
  testInteractionContainment();
  testColumnMaskContainment();
  testBlockAdditiveConfinement();
}
