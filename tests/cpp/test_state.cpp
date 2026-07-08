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

  std::vector<size_t> indices(n);
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
  std::vector<size_t> indices2(n);
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

    std::vector<size_t> indices(n);
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
    std::vector<size_t> rebuiltIndices(n);
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
  check(restored.setState(state), "scaled state restores");

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
  check(restored.setState(state), "a stored state restores");

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
  check(restored.setState(state), "a binary+dart state restores");

  checkStructuralRoundTrip(state, restored,
                           "restored logistic + dart state reproduces the model");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  (void) rng;
  printf("ok: state round trip with latents\n");
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
  check(!sampler.setState(bad), "setState rejects a chain-count mismatch");

  bad = state;
  bool corrupted = false;
  for (std::vector<FlatNode>& tree : bad.chains[0].forests[0].trees) {
    if (tree.size() > 1) {
      tree[0].value += 1.0e-3;  // off the cut grid
      corrupted = true;
      break;
    }
  }
  check(corrupted && !sampler.setState(bad),
        "setState rejects a split value off the cut grid");
  check(sampler.data().cutPoints == cutsBefore,
        "a rejected state leaves the cuts untouched");
  check(sampler.chain(0).treeFits() == treeFitsBefore,
        "a rejected state leaves the fits untouched");

  check(sampler.setState(state), "the original state still restores");
  check(sampler.chain(0).treeFits() == treeFitsBefore,
        "restoring the current state is an identity");

  printf("ok: state validation\n");
}

void runStateTests(ext_rng* rng) {
  testFlattenRoundTrip();
  testCategoricalFlattenBoundaries();
  testKeepTrees(rng);
  testPredictCurrentTrees(rng);
  testStateRoundTrip();
  testStateRoundTripScaledOffset();
  testStateRoundTripLatents(rng);
  testStateValidation(rng);
}
