#include "common.hpp"

static void testFlattenRoundTrip() {
  // one categorical column (codes 0..3) and one ordinal, a hand-built tree
  // with one rule of each kind
  const size_t n = 12;
  std::vector<double> x(n * 2), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 4);
  for (size_t i = 0; i < n; ++i) x[i + n] = static_cast<double>(i) / (n - 1.0);
  ColumnKind types[] = {ColumnKind::categorical, ColumnKind::numeric};

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
    ColumnKind types[] = {ColumnKind::categorical};
    ColumnStore store;
    store.build(x.data(), n, 1, 10, false, types);
    check(store.categoryCounts[0] == K &&
            store.columnIsPooled(0) == (K >= 64),
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
  sampler.predict(xTest.data(), nTest, 1, predicted.data());
  check(predicted == testFits, "saved-tree predictions equal the run's test fits");

  // a second recorded run overwrites the OLDEST slots, and the read walks the
  // ring from the write cursor: the draws that survive come first, in the
  // order they were drawn, and the new ones land at the tail
  std::vector<double> sigma2(2), testFits2(nTest * 2);
  Results results2;
  results2.sigma = sigma2.data();
  results2.testFits = testFits2.data();
  sampler.run(0, 2, results2);
  check(sampler.currentSampleNum() == 2, "a partial run advances the slot");
  check(sampler.filledSavedDraws() == numSamples,
        "a full store stays full across a partial run");

  std::vector<double> predicted2(nTest * numSamples);
  sampler.predict(xTest.data(), nTest, 1, predicted2.data());
  bool preserved = std::equal(predicted.begin() + nTest * 2, predicted.end(),
                              predicted2.begin());
  bool overwritten = std::equal(testFits2.begin(), testFits2.end(),
                                predicted2.begin() + nTest * 2);
  check(preserved, "the surviving draws lead, in the order they were drawn");
  check(overwritten, "the new draws land at the tail");
  check(sampler.savedTreeCapacity() == numSamples,
        "keepTrees capacity comes from the options");

  printf("ok: keepTrees\n");
}

// The saved-tree read map, at the two shapes testKeepTrees does not reach: a
// store still filling, and one a run wrapped past. Output draw i is slot
// (currentSampleNum + capacity - filled + i) % capacity over the filled most
// recent recorded draws, oldest first.
static void testSavedDrawOrder(ext_rng* rng) {
  const size_t n = 200, nTest = 20, capacity = 4;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = runif01();

  SamplerOptions options;
  options.numTrees = 25;
  options.keepTrees = true;
  options.numSamplesToStore = capacity;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  sampler.setTestPredictors(xTest.data(), nTest);

  check(sampler.filledSavedDraws() == 0, "a fresh store holds no draws");

  // partial fill: three draws in a store of four, read at the head
  std::vector<double> sigma1(3), testFits1(nTest * 3);
  Results results1;
  results1.sigma = sigma1.data();
  results1.testFits = testFits1.data();
  sampler.run(10, 3, results1);
  check(sampler.currentSampleNum() == 3 && sampler.filledSavedDraws() == 3,
        "a short run fills three of four slots");
  bool headMap = true;
  for (size_t i = 0; i < 3; ++i)
    headMap &= sampler.savedSlotForDraw(i) == i;
  check(headMap, "a partly filled store reads from slot 0");

  std::vector<double> predicted1(nTest * 3);
  sampler.predict(xTest.data(), nTest, 1, predicted1.data());
  check(predicted1 == testFits1,
        "a partly filled store replays its own draws, and only those");

  // wrapping: three more draws past the capacity boundary
  std::vector<double> sigma2(3), testFits2(nTest * 3);
  Results results2;
  results2.sigma = sigma2.data();
  results2.testFits = testFits2.data();
  sampler.run(0, 3, results2);
  check(sampler.currentSampleNum() == 2 && sampler.filledSavedDraws() == 4,
        "wrapping past capacity fills the store and moves the cursor");
  bool wrappedMap = sampler.savedSlotForDraw(0) == 2 &&
                    sampler.savedSlotForDraw(1) == 3 &&
                    sampler.savedSlotForDraw(2) == 0 &&
                    sampler.savedSlotForDraw(3) == 1;
  check(wrappedMap, "a wrapped store reads from the cursor, oldest first");

  std::vector<double> predicted2(nTest * 4);
  sampler.predict(xTest.data(), nTest, 1, predicted2.data());
  check(std::equal(testFits2.begin(), testFits2.end(),
                   predicted2.begin() + nTest),
        "the second run's draws are the three most recent, in order");
  check(std::equal(testFits1.begin() + nTest * 2, testFits1.end(),
                   predicted2.begin()),
        "the one surviving draw of the first run leads");

  // the store's draws belong to the fit that recorded them: resizing it drops
  // them, and so does a read of the emptied store
  sampler.setTreeStorage(true, capacity + 1);
  check(sampler.filledSavedDraws() == 0 && sampler.currentSampleNum() == 0,
        "resizing the store discards what it held");

  printf("ok: saved draw order\n");
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
  sampler.predict(x.data(), n, 1, predicted.data());
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
  original.predict(xTest.data(), 20, 1, predictionsA.data());
  restored.predict(xTest.data(), 20, 1, predictionsB.data());
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
  check(state.chains.size() == numChains && state.currentSampleNum == 2 &&
          state.recordedDraws == 2,
        "state captures every chain, the slot position and the draw count");

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
  original.predict(xTest.data(), 20, 1, predictA.data());
  restored.predict(xTest.data(), 20, 1, predictB.data());
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
    std::vector<double> sigmaA(window), sigmaB(window), nuA(window, -1.0);
    Results resultsA, resultsB;
    resultsA.sigma = sigmaA.data();
    // the recorded per-draw df, the run-channel twin of the state's nu: it is
    // written from settled state, so the last draw must be the nu the sampler
    // now holds and a fixed nu must repeat exactly
    resultsA.residualDf = nuA.data();
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

    bool dfRecorded = true;
    for (size_t i = 0; i < window; ++i)
      if (!(nuA[i] > 0.0) || (!estimated && nuA[i] != residualDf))
        dfRecorded = false;
    check(dfRecorded, "the df channel records a positive nu every draw");
    SamplerStateData postState;
    original.getState(postState);
    check(nuA[window - 1] == postState.chains[0].residualDf,
          "the last recorded df is the nu the sampler holds");

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

  // the df write is gated on the response, not on the caller: a gaussian
  // sampler handed the channel leaves it exactly as it found it
  std::vector<double> nuG(4, -7.0);
  Results gaussChannel;
  gaussChannel.residualDf = nuG.data();
  gaussian.run(0, 4, gaussChannel);
  bool dfUntouched = true;
  for (size_t i = 0; i < nuG.size(); ++i)
    if (nuG[i] != -7.0) dfUntouched = false;
  check(dfUntouched, "a gaussian sampler leaves the df channel untouched");

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

// Cross-grid warm start (docs/plans/warm-starts.md): a donor grown on a fine
// cut grid seeds a destination on a coarser one. The refusal is lifted; the
// donor's split indices remap onto the live grid and any the coarser grid
// starves collapse - setData's contract - so every installed tree stays
// occupied (no empty leaf a naive remap would leave), the install reports ok,
// and the sampler runs. Mirrors testSetDataQuantileShrink's collapse gate.
static void testCrossGridWarmStart() {
  std::uint64_t savedRngState = rngState;
  rngState = 909090u;

  // a single tree is forced to be deep: it must fit the whole box signal
  // below alone, splitting column 0 twice on one path (a sum of shallow trees
  // would spread the two thresholds across trees and never starve)
  const size_t n = 200, p = 2, numTrees = 1;

  // borrow-lifetime bookkeeping: each Sampler holds its rng pointer, so the
  // rng objects must outlive their samplers
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  };

  // donor: column 0 takes 10 discrete levels -> 9 quantile cuts. A box signal
  // (large only for x0 in [3, 6]) needs TWO column-0 thresholds on one path, so
  // the donor's trees split column 0 twice along a branch - splits a coarse
  // destination grid cannot both keep
  std::vector<double> xDonor(n * p), yDonor(n);
  for (size_t i = 0; i < n; ++i) {
    xDonor[i] = static_cast<double>(i % 10);
    xDonor[i + n] = runif01();
    yDonor[i] = 4.0 * ((xDonor[i] >= 3.0 && xDonor[i] <= 6.0) ? 1.0 : 0.0) +
                0.5 * xDonor[i + n] + 0.2 * (runif01() - 0.5);
  }
  SamplerOptions donorOptions;
  donorOptions.numTrees = numTrees;
  donorOptions.useQuantiles = true;
  ext_rng* donorRng = newRng(313);
  ConstantLeafSampler donor(xDonor.data(), yDonor.data(), n, p, nullptr,
                            nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                            0.37804942330213542, donorOptions, &donorRng);
  Results empty;
  donor.run(200, 0, empty);
  check(donor.data().numCuts[0] == 9, "cross-grid: donor holds a 9-cut grid");
  SamplerStateData donorState;
  donor.getState(donorState);

  // destination: column 0 takes just 2 levels -> a single cut, so any donor
  // path that split column 0 more than once starves and must collapse
  std::vector<double> xDest(n * p), yDest(n);
  for (size_t i = 0; i < n; ++i) {
    xDest[i] = static_cast<double>(i % 2);
    xDest[i + n] = runif01();
    yDest[i] = 2.0 * xDest[i] + 0.5 * xDest[i + n] + 0.2 * (runif01() - 0.5);
  }
  SamplerOptions destOptions;
  destOptions.numTrees = numTrees;
  destOptions.useQuantiles = true;
  ext_rng* destRng = newRng(717);
  ConstantLeafSampler dest(xDest.data(), yDest.data(), n, p, nullptr, nullptr,
                           ResponseFamily::gaussian, 1.0, 3.0,
                           0.37804942330213542, destOptions, &destRng);
  check(dest.data().numCuts[0] == 1,
        "cross-grid: destination holds a single-cut column-0 grid");

  // the refusal is lifted: a genuinely different grid now installs
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(dest.installForests(donorState, liveMap) == WarmStartResult::ok,
        "cross-grid: a different-grid donor warm-starts by remapping");

  // no starved split left an empty leaf; the remap only ever merges leaves, so
  // the live leaf count is below the donor's (a genuine collapse fired) and yet
  // real structure survived (more than one leaf remains)
  bool occupied = true;
  size_t donorLeaves = 0, liveLeaves = 0;
  std::vector<int32_t> bottoms;
  for (size_t t = 0; t < numTrees; ++t) {
    occupied &= dest.chain(0).tree(t).bottomNodesAreOccupied();
    for (const FlatNode& node : donorState.chains[0].forests[0].trees[t])
      if (node.variable == invalidVariable) ++donorLeaves;
    dest.chain(0).tree(t).fillBottom(0, bottoms);
    liveLeaves += bottoms.size();
  }
  check(occupied, "cross-grid: remap collapses starved splits, no empty leaves");
  check(liveLeaves < donorLeaves,
        "cross-grid: a starved split collapsed (fewer live leaves than donor)");
  check(liveLeaves > numTrees,
        "cross-grid: the donor's structure survived the remap");

  // the sampler runs from the warm-started state
  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  dest.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "cross-grid: sampler runs after a cross-grid warm start");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: cross-grid warm start (%zu donor -> %zu live leaves)\n",
         donorLeaves, liveLeaves);
}

// Warm start under a variance forest (docs/plans/variance-forest-mutation-
// routing.md, slice S5). installForests used to reassemble a state carrying no
// variance trees at all, so the destination adopted the donor's mean forest
// while keeping its own cold scale surface. The gate is STATE-level - the
// donor's variance trees ARE the destination's immediately after the install -
// and deliberately not behavioral: a behavioral probe of this shape was run
// twice during design and could not separate warm from cold, since one sweep
// of a variance forest on a strong scale signal already recovers the surface.
// Three refusals ride the same fixture: a rebuilt variance tree that leaves a
// bottom unoccupied (that tree would report a scale this data never
// supported), one that splits outside a restricted destination's variance
// columns, and their specificity arms.
static void testVarianceWarmStart() {
  std::uint64_t savedRngState = rngState;
  rngState = 515151u;

  const size_t n = 300, p = 2, numTrees = 20, numVarianceTrees = 5;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    // signal in the mean (x0) and in the scale (x1), so both forests split
    y[i] = 3.0 * x[i] + (x[i + n] < 0.5 ? 0.2 : 1.4) * z;
  }

  std::vector<ext_rng*> rngs;  // each Sampler holds its rng; outlive them all
  const std::vector<size_t> allowZero = {0};
  auto makeSampler = [&](std::uint32_t seed, std::uint32_t maxNumCuts,
                         bool restrictVarianceToZero) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numVarianceTrees = numVarianceTrees;
    options.maxNumCuts = maxNumCuts;
    if (restrictVarianceToZero) {
      options.varianceForestColumns = allowZero.data();
      options.numVarianceForestColumns = allowZero.size();
    }
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, options, &r);
  };
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  Results empty;

  auto donor = makeSampler(2024, 100, false);
  donor->run(60, 0, empty);
  SamplerStateData donorState;
  donor->getState(donorState);

  // (1) same grid: the donor's variance trees install verbatim
  auto dest = makeSampler(4048, 100, false);
  dest->run(60, 0, empty);
  SamplerStateData before;
  dest->getState(before);
  check(!sameFlatTrees(before.chains[0].varianceTrees,
                       donorState.chains[0].varianceTrees),
        "variance warm start: the destination's own surface differs first");
  check(dest->installForests(donorState, liveMap) == WarmStartResult::ok,
        "variance warm start: a heteroscedastic donor installs");
  SamplerStateData after;
  dest->getState(after);
  check(sameFlatTrees(after.chains[0].varianceTrees,
                      donorState.chains[0].varianceTrees),
        "variance warm start: the donor's variance trees are the "
        "destination's");

  // (2) cross grid: the donor's thresholds remap onto a coarser destination
  // grid, so the trees are NOT identical - what must hold is that every factor
  // survives positive and the flattened state is legal against the new grid
  auto coarse = makeSampler(6072, 8, false);
  check(coarse->data().cutPoints != donorState.cutPoints,
        "variance warm start: the coarse destination is on another grid");
  check(coarse->installForests(donorState, liveMap) == WarmStartResult::ok,
        "variance warm start: a cross-grid donor installs by remapping");
  const double* remapped = coarse->chain(0).varianceFits();
  bool positive = true;
  for (size_t i = 0; i < n; ++i)
    positive &= std::isfinite(remapped[i]) && remapped[i] > 0.0;
  check(positive, "variance warm start: every remapped factor stays positive");
  SamplerStateData remappedState;
  coarse->getState(remappedState);
  check(coarse->setState(remappedState, nullptr),
        "variance warm start: the remapped surface serializes legally");

  // (3) an install leaving a variance bottom unoccupied is refused: x0 <=
  // cut[2] then x0 > cut[8] is a region no row can reach
  SamplerStateData emptyBottom = donorState;
  const std::vector<double>& cuts(donorState.cutPoints[0]);
  check(cuts.size() > 8, "variance warm start: enough cuts for the nesting");
  std::vector<FlatNode>& stranded(emptyBottom.chains[0].varianceTrees[0]);
  stranded.assign(5, FlatNode());
  stranded[0].variable = 0;
  stranded[0].value = cuts[2];
  setFlatKind(stranded[0], FlatKind::ordinal);
  stranded[1].variable = 0;
  stranded[1].value = cuts[8];
  setFlatKind(stranded[1], FlatKind::ordinal);
  stranded[2].value = 1.3;
  stranded[3].value = 0.7;  // the unreachable bottom
  stranded[4].value = 1.1;
  auto strandTarget = makeSampler(8096, 100, false);
  check(strandTarget->installForests(emptyBottom, liveMap) ==
          WarmStartResult::varianceMismatch,
        "variance warm start: an unoccupied variance bottom is refused");

  // (4) a donor variance tree splitting outside a `variance = ~ x0`
  // destination's columns is refused, and a compliant one is not
  SamplerStateData outOfMask = donorState;
  std::vector<FlatNode>& onOne(outOfMask.chains[0].varianceTrees[0]);
  onOne.assign(3, FlatNode());
  onOne[0].variable = 1;
  onOne[0].value = donorState.cutPoints[1][10];
  setFlatKind(onOne[0], FlatKind::ordinal);
  onOne[1].value = 1.2;
  onOne[2].value = 0.8;
  auto restricted = makeSampler(1120, 100, true);
  check(restricted->installForests(outOfMask, liveMap) ==
          WarmStartResult::columnMaskMismatch,
        "variance warm start: an out-of-mask variance tree is refused");
  SamplerStateData compliant = outOfMask;
  for (std::vector<FlatNode>& tree : compliant.chains[0].varianceTrees) {
    tree.assign(1, FlatNode());
    tree[0].value = 1.1;
  }
  auto restricted2 = makeSampler(1344, 100, true);
  check(restricted2->installForests(compliant, liveMap) == WarmStartResult::ok,
        "variance warm start: an in-mask variance donor installs");

  // (5) setState is held to the rule by the SAME predicate, so the state one
  // entry refuses the other refuses too, and the refusal is named rather than
  // folded into "not consistent". It is taken before any chain is touched, so
  // the destination keeps the surface it had.
  restricted->run(30, 0, empty);
  SamplerStateData restrictedBefore;
  restricted->getState(restrictedBefore);
  bool columnMaskRefused = false;
  check(!restricted->setState(outOfMask, nullptr, &columnMaskRefused),
        "variance setState: an out-of-mask variance tree is refused");
  check(columnMaskRefused,
        "variance setState: the refusal is named as the column-mask one");
  SamplerStateData restrictedAfter;
  restricted->getState(restrictedAfter);
  check(sameFlatTrees(restrictedBefore.chains[0].varianceTrees,
                      restrictedAfter.chains[0].varianceTrees),
        "variance setState: a refused restore leaves the surface untouched");
  check(restricted->setState(compliant, nullptr, &columnMaskRefused) &&
          !columnMaskRefused,
        "variance setState: an in-mask variance state restores");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: variance-forest warm start\n");
}

// The per-forest leaf scale rides the state (docs/plans/multiforest-mutation-
// gaps.md item 3). BCF derives both forests' scales from the response's SHAPE,
// so a destination built on a differently shaped response constructs different
// ones and a restore must install the donor's. Forest 1's scale is unreadable
// through Chain::leaf() (forest 0 only), so gate both through a re-captured
// state; statesAgree compares the field, which is what keeps the fuzzer's
// OP_STATE round trip covering it.
static void testStateLeafScale(ext_rng* rng) {
  const size_t n = 300, p = 2;
  std::vector<double> x(n * p), z(n), y1(n), y2(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    y1[i] = std::sin(3.0 * x[i]) + x[i + n] + z[i] * (1.0 + x[i + n]);
  }
  // y2 keeps y1's RANGE and squeezes its interior: the BCF anchor is the sd of
  // the range-scaled response, so an affine rescaling would leave it identical
  // and only a shape change moves it
  double lo = *std::min_element(y1.begin(), y1.end());
  double hi = *std::max_element(y1.begin(), y1.end());
  double mid = 0.5 * (lo + hi);
  for (size_t i = 0; i < n; ++i) y2[i] = mid + 0.2 * (y1[i] - mid);
  y2[0] = lo;
  y2[n - 1] = hi;

  SamplerOptions options;
  AmplitudeSpec spec;
  spec.mu.numTrees = 20;
  spec.tau.numTrees = 10;
  spec.z = z.data();
  auto make = [&](const double* y) {
    return std::make_unique<Sampler<ConstantGaussianLeaf>>(
      x.data(), y, n, p, nullptr, nullptr, 1.0, 3.0, 0.37804942330213542,
      options, spec, &rng);
  };

  auto donor = make(y1.data());
  Results empty;
  donor->run(10, 2, empty);
  SamplerStateData donorState;
  donor->getState(donorState);
  const ForestStateData& d0 = donorState.chains[0].forests[0];
  const ForestStateData& d1 = donorState.chains[0].forests[1];
  check(d0.leafScale > 0.0 && d1.leafScale > 0.0,
        "leaf scale: every forest's scale is stored");
  check(d0.leafScale != d1.leafScale,
        "leaf scale: BCF's two forests calibrate differently");

  auto dest = make(y2.data());
  SamplerStateData destState;
  dest->getState(destState);
  // not a vacuous arm: the destination constructed its own, different scales
  check(destState.chains[0].forests[0].leafScale != d0.leafScale &&
          destState.chains[0].forests[1].leafScale != d1.leafScale,
        "leaf scale: a different-shape response constructs different scales");

  check(dest->setState(donorState, nullptr), "leaf scale: the donor restores");
  SamplerStateData reState;
  dest->getState(reState);
  check(reState.chains[0].forests[0].leafScale == d0.leafScale &&
          reState.chains[0].forests[1].leafScale == d1.leafScale,
        "leaf scale: both forests install the donor's scale");

  // a state stripped of the block (every pre-block state) restores exactly as
  // it did before: the destination keeps the scales it constructed
  SamplerStateData stripped = donorState;
  for (ForestStateData& fs : stripped.chains[0].forests) fs.leafScale = 0.0;
  auto old = make(y2.data());
  check(old->setState(stripped, nullptr),
        "leaf scale: a state without the block restores");
  SamplerStateData oldState;
  old->getState(oldState);
  check(oldState.chains[0].forests[0].leafScale ==
          destState.chains[0].forests[0].leafScale &&
          oldState.chains[0].forests[1].leafScale ==
          destState.chains[0].forests[1].leafScale,
        "leaf scale: an absent block leaves construction's scales");

  // installForests is the OTHER restore path, and it reassembles its own
  // ForestStateData from the donor's - so it needs the field copied through
  // (sampler.hpp) or its install arm would be dead
  auto warm = make(y2.data());
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(warm->installForests(donorState, liveMap) == WarmStartResult::ok,
        "leaf scale: the donor warm-starts");
  SamplerStateData warmState;
  warm->getState(warmState);
  check(warmState.chains[0].forests[0].leafScale == d0.leafScale &&
          warmState.chains[0].forests[1].leafScale == d1.leafScale,
        "leaf scale: a warm start adopts the donor's scale too");

  printf("ok: per-forest leaf scale rides the state\n");
}

// The variance forest's SAVED (keepTrees) trees ride the state: a re-created
// sampler must replay the recorded s^2(x) slot for slot rather than the
// multiplicative identity initializeSavedTrees left in the buffer. The live
// trees alone are not enough - predict addresses the saved buffer, never them.
static void testVarianceSavedTreeState() {
  std::uint64_t savedRngState = rngState;
  rngState = 828282u;

  const size_t n = 300, p = 2, numTrees = 20, numVarianceTrees = 5;
  const size_t numSamples = 4, nTest = 6;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    // signal in the mean (x0) and in the scale (x1), so both forests split
    y[i] = 3.0 * x[i] + (x[i + n] < 0.5 ? 0.2 : 1.4) * z;
  }
  std::vector<double> xTest(nTest * p);
  for (size_t i = 0; i < nTest; ++i) {
    xTest[i] = static_cast<double>(i) / (nTest - 1.0);
    xTest[i + nTest] = static_cast<double>(nTest - 1 - i) / (nTest - 1.0);
  }

  std::vector<ext_rng*> rngs;  // each Sampler holds its rng; outlive them all
  auto makeSampler = [&](std::uint32_t seed, size_t numVariance) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numVarianceTrees = numVariance;
    options.keepTrees = true;
    options.numSamplesToStore = numSamples;
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, options, &r);
  };

  std::vector<double> sigma(numSamples);
  Results results;
  results.sigma = sigma.data();

  auto donor = makeSampler(3131, numVarianceTrees);
  donor->run(60, numSamples, results);
  std::vector<double> before(nTest * numSamples);
  donor->predictVariance(xTest.data(), nTest, 1, before.data());
  bool positive = true;
  for (double v : before) positive = positive && v > 0.0;
  check(positive, "variance saved state: the donor's replay is positive");

  SamplerStateData state;
  donor->getState(state);
  check(state.chains[0].savedVarianceTrees.size() ==
          numSamples * numVarianceTrees,
        "variance saved state: the block is capacity x variance trees");

  auto dest = makeSampler(6262, numVarianceTrees);
  dest->run(60, numSamples, results);
  std::vector<double> own(nTest * numSamples);
  dest->predictVariance(xTest.data(), nTest, 1, own.data());
  check(own != before,
        "variance saved state: the destination's own surface differs first");
  check(dest->setState(state, nullptr),
        "variance saved state: a heteroscedastic state installs");
  std::vector<double> after(nTest * numSamples);
  dest->predictVariance(xTest.data(), nTest, 1, after.data());
  check(after == before,
        "variance saved state: the restored slots replay bitwise");
  checkStructuralRoundTrip(state, *dest,
                           "variance saved state: the re-captured state agrees");

  // refusals, each leaving the destination untouched
  SamplerStateData nonPositive = state;
  // the LAST record of a pre-order tree is always a leaf
  nonPositive.chains[0].savedVarianceTrees[0].back().value = 0.0;
  check(!dest->setState(nonPositive, nullptr),
        "variance saved state: a non-positive saved leaf is refused");
  SamplerStateData malformed = state;
  malformed.chains[0].savedVarianceTrees[1].clear();
  check(!dest->setState(malformed, nullptr),
        "variance saved state: a malformed saved tree is refused");
  SamplerStateData truncated = state;
  truncated.chains[0].savedVarianceTrees.pop_back();
  check(!dest->setState(truncated, nullptr),
        "variance saved state: a truncated saved block is refused");
  // an EMPTY block against a live capacity can only be a state written before
  // the channel existed; accepting it would restore the identity fill and
  // report a plausible constant s(x)
  SamplerStateData noBlock = state;
  noBlock.chains[0].savedVarianceTrees.clear();
  check(!dest->setState(noBlock, nullptr),
        "variance saved state: an empty block under a live capacity is "
        "refused");
  std::vector<double> unchanged(nTest * numSamples);
  dest->predictVariance(xTest.data(), nTest, 1, unchanged.data());
  check(unchanged == after,
        "variance saved state: a refused install changes nothing");

  // and a homoscedastic sampler refuses a state carrying the block
  auto plain = makeSampler(9393, 0);
  check(!plain->setState(state, nullptr),
        "variance saved state: a homoscedastic sampler refuses the block");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: variance-forest saved-tree state\n");
}

// A warm start from a SAVED sample takes that sample's own scale surface, not
// the donor's live one: the mean and variance saved buffers are index-aligned
// by construction (one slot base drives both, both index by the sample number),
// so the destination's live variance trees after the install are exactly the
// donor's saved slice. State-level for the reason testVarianceWarmStart's gate
// is. Four refusals ride the same fixture: a buffer that does not cover the
// named slot, an absent one, one whose STRIDE disagrees with the donor's
// variance tree count (the only way a slice can cross slot boundaries, and
// invisible to any check downstream, which sees a correctly sized vector), and
// a saved slot carrying a non-positive scale leaf.
static void testVarianceWarmStartSlot() {
  std::uint64_t savedRngState = rngState;
  rngState = 727272u;

  const size_t n = 300, p = 2, numTrees = 20, numVarianceTrees = 5;
  const size_t numSamples = 4;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = 3.0 * x[i] + (x[i + n] < 0.5 ? 0.2 : 1.4) * z;
  }

  std::vector<ext_rng*> rngs;
  auto makeSampler = [&](std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numVarianceTrees = numVarianceTrees;
    options.keepTrees = true;
    options.numSamplesToStore = numSamples;
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, options, &r);
  };
  std::vector<double> sigma(numSamples);
  Results results;
  results.sigma = sigma.data();

  auto donor = makeSampler(1717);
  donor->run(60, numSamples, results);
  SamplerStateData donorState;
  donor->getState(donorState);
  const std::vector<std::vector<FlatNode>>& savedVariance(
    donorState.chains[0].savedVarianceTrees);
  check(savedVariance.size() == numSamples * numVarianceTrees,
        "variance slot warm start: the donor's saved block is capacity x trees");
  auto slotTrees = [&](const SamplerStateData& state, size_t slot) {
    const std::vector<std::vector<FlatNode>>& block(
      state.chains[0].savedVarianceTrees);
    return std::vector<std::vector<FlatNode>>(
      block.begin() + slot * numVarianceTrees,
      block.begin() + (slot + 1) * numVarianceTrees);
  };

  // the first and last slots both install their own surface; the first is also
  // the discriminating one, since the last sweep's saved trees are the live
  // ones and would pass even under the old live-copy reassembly
  for (size_t slot : {size_t(0), numSamples - 1}) {
    auto dest = makeSampler(3434 + static_cast<std::uint32_t>(slot));
    dest->run(60, numSamples, results);
    SamplerStateData before;
    dest->getState(before);
    check(!sameFlatTrees(before.chains[0].varianceTrees,
                         slotTrees(donorState, slot)),
          "variance slot warm start: the destination's own surface differs "
          "first");
    std::vector<std::pair<size_t, int>> slotMap = {{0, static_cast<int>(slot)}};
    check(dest->installForests(donorState, slotMap) == WarmStartResult::ok,
          "variance slot warm start: a slot-sourced heteroscedastic donor "
          "installs");
    SamplerStateData after;
    dest->getState(after);
    check(sameFlatTrees(after.chains[0].varianceTrees,
                        slotTrees(donorState, slot)),
          "variance slot warm start: the named sample's scale surface is the "
          "destination's");
    if (slot == 0)
      check(!sameFlatTrees(after.chains[0].varianceTrees,
                           donorState.chains[0].varianceTrees),
            "variance slot warm start: and it is not the donor's live surface");
  }

  auto target = makeSampler(5656);
  target->run(60, numSamples, results);
  SamplerStateData untouched;
  target->getState(untouched);

  // a one-short buffer: the last slot is the one it strands, so name it
  SamplerStateData shortBuffer = donorState;
  shortBuffer.chains[0].savedVarianceTrees.pop_back();
  std::vector<std::pair<size_t, int>> lastMap = {
    {0, static_cast<int>(numSamples - 1)}};
  check(target->installForests(shortBuffer, lastMap) ==
          WarmStartResult::varianceSlotMismatch,
        "variance slot warm start: a short saved buffer is refused");

  SamplerStateData noBlock = donorState;
  noBlock.chains[0].savedVarianceTrees.clear();
  check(target->installForests(noBlock, lastMap) ==
          WarmStartResult::varianceSlotMismatch,
        "variance slot warm start: an absent saved buffer is refused");

  // stride: a live block shorter than the buffer's stride would slice across
  // two sweeps' trees and still hand installVarianceForest the right COUNT
  SamplerStateData spliced = donorState;
  spliced.chains[0].varianceTrees.resize(numVarianceTrees - 2);
  check(target->installForests(spliced, lastMap) ==
          WarmStartResult::varianceSlotMismatch,
        "variance slot warm start: a stride the live block contradicts is "
        "refused");

  // the positivity law the saved buffer is held to on the state side: a zero
  // scale leaf annihilates the product a rebuild forms
  SamplerStateData nonPositive = donorState;
  nonPositive.chains[0].savedVarianceTrees[0].back().value = 0.0;
  std::vector<std::pair<size_t, int>> firstMap = {{0, 0}};
  check(target->installForests(nonPositive, firstMap) ==
          WarmStartResult::varianceMismatch,
        "variance slot warm start: a non-positive saved scale leaf is refused");
  // and the refusal is SLOT-specific: the same state installs from any slot
  // whose own trees are intact
  check(target->installForests(nonPositive, lastMap) == WarmStartResult::ok,
        "variance slot warm start: another slot of the same state installs");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: variance-forest slot-sourced warm start\n");
}

static void testWeightsDigest() {
  // The two primitives the saved-state seam pairs to keep restored latents and
  // the weights in force from disagreeing: a digest that separates weight
  // vectors, and a repair that re-derives whatever is stated against them.
  // Engine level, so no digest is written or read here - a round trip at this
  // level carries none, which is why the byte-for-byte pins above are
  // untouched by the seam that does. RNG-insulated per testLeafOfConsistency.
  std::uint64_t savedRngState = rngState;
  rngState = 818181u;

  const size_t n = 90;
  std::vector<double> x(n * 2), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = (x[i] + 0.25 * (runif01() - 0.5) > 0.5) ? 1.0 : 0.0;

  // 1,4,4 against 2,2,5: equal in sum and in sum of squares, different in
  // bytes. A digest built from moments would call these the same weights.
  std::vector<double> wA(n), wB(n), ones(n, 1.0);
  double sumA = 0.0, sumB = 0.0, sumSqA = 0.0, sumSqB = 0.0;
  for (size_t i = 0; i < n; ++i) {
    wA[i] = i % 3 == 0 ? 1.0 : 4.0;
    wB[i] = i % 3 == 2 ? 5.0 : 2.0;
    sumA += wA[i];
    sumB += wB[i];
    sumSqA += wA[i] * wA[i];
    sumSqB += wB[i] * wB[i];
  }
  check(sumA == sumB && sumSqA == sumSqB,
        "the probe weight vectors share their moments");

  SamplerOptions options;
  options.numTrees = 10;
  options.nodeScale = 3.0;

  ext_rng* rngs[6];
  for (size_t i = 0; i < 6; ++i) {
    rngs[i] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[i], 8181);
  }

  ConstantLeafSampler weightedA(x.data(), y.data(), n, 2, wA.data(), nullptr,
                                ResponseFamily::logistic, 1.0, 3.0, 1.0,
                                options, &rngs[0]);
  ConstantLeafSampler weightedB(x.data(), y.data(), n, 2, wB.data(), nullptr,
                                ResponseFamily::logistic, 1.0, 3.0, 1.0,
                                options, &rngs[1]);
  check(weightedA.weightsDigest() != weightedB.weightsDigest(),
        "weight vectors sharing every moment still digest apart");

  ConstantLeafSampler unweighted(x.data(), y.data(), n, 2, nullptr, nullptr,
                                 ResponseFamily::logistic, 1.0, 3.0, 1.0,
                                 options, &rngs[2]);
  ConstantLeafSampler unitWeighted(x.data(), y.data(), n, 2, ones.data(),
                                   nullptr, ResponseFamily::logistic, 1.0,
                                   3.0, 1.0, options, &rngs[3]);
  check(unweighted.weightsDigest() == unitWeighted.weightsDigest(),
        "no weights and all ones are one sampler and so one digest");

  Results empty;
  weightedA.run(20, 0, empty);
  std::vector<double> beforeRepair(weightedA.latents(0),
                                   weightedA.latents(0) + n);
  weightedA.reapplyWeights();
  std::vector<double> afterRepair(weightedA.latents(0),
                                  weightedA.latents(0) + n);
  check(beforeRepair != afterRepair,
        "reapplyWeights redraws a logistic chain's Polya-Gamma latents");

  // a gaussian chain states nothing against its weights, so the repair draws
  // nothing and moves nothing: the chain that took it draws what the chain
  // that did not draws
  std::vector<double> yContinuous(n);
  for (size_t i = 0; i < n; ++i)
    yContinuous[i] = 2.0 * x[i] + 0.3 * (runif01() - 0.5);
  ConstantLeafSampler repaired(x.data(), yContinuous.data(), n, 2, wA.data(),
                               nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                               1.0, options, &rngs[4]);
  ConstantLeafSampler untouched(x.data(), yContinuous.data(), n, 2, wA.data(),
                                nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                                1.0, options, &rngs[5]);
  repaired.run(10, 0, empty);
  untouched.run(10, 0, empty);
  repaired.reapplyWeights();
  const size_t draws = 5;
  std::vector<double> sigmaRepaired(draws), sigmaUntouched(draws);
  Results resultsRepaired, resultsUntouched;
  resultsRepaired.sigma = sigmaRepaired.data();
  resultsUntouched.sigma = sigmaUntouched.data();
  repaired.run(0, draws, resultsRepaired);
  untouched.run(0, draws, resultsUntouched);
  check(sigmaRepaired == sigmaUntouched,
        "reapplyWeights is inert on a chain that reads weights as precisions");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: weights digest and repair\n");
}

// The saved-store readers' DESTINATION stride, on the one shape every other
// saved-tree test leaves out: more than one chain over a store that is not
// full. The readers loop over the retained draws and must stride the output by
// that same count - out is sized slab x filledSavedDraws x numChains - so a
// stride taken from the CAPACITY instead lands chain c >= 1 past the end of
// its slab, which is a heap write past the buffer for the last chain and a
// never-written hole for the first. Single-chain fixtures cannot see it:
// c * numDraws and c * capacity are both zero there.
//
// The oracle is the run's own recorded test fits, which the replay must
// reproduce draw for draw, plus a poisoned guard region past the buffer that a
// capacity-strided write falls into.
static void testMultiChainPartialFillPredict() {
  // own generators and own runif01 stream, restored on the way out, so the
  // suites after this one read the same draws with or without it
  std::uint64_t savedRngState = rngState;
  const size_t n = 200, nTest = 20, capacity = 4, numChains = 2, numDraws = 3;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = runif01();

  std::vector<ext_rng*> rngs(numChains);
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(rngs[c], 20260824u + static_cast<uint_least32_t>(c));
  }

  SamplerOptions options;
  options.numTrees = 25;
  options.numChains = numChains;
  options.keepTrees = true;
  options.numSamplesToStore = capacity;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, rngs.data());
  sampler.setTestPredictors(xTest.data(), nTest);

  std::vector<double> sigma(numDraws * numChains),
    testFits(nTest * numDraws * numChains);
  Results results;
  results.sigma = sigma.data();
  results.testFits = testFits.data();
  sampler.run(10, numDraws, results);
  check(sampler.filledSavedDraws() == numDraws &&
          sampler.savedTreeCapacity() == capacity,
        "the multi-chain store holds three of four slots");

  const double poison = -1.0e300;
  size_t slab = nTest, live = slab * numDraws * numChains;
  std::vector<double> out(live + 2 * slab, poison);
  sampler.predict(xTest.data(), nTest, 1, out.data());
  check(std::equal(testFits.begin(), testFits.end(), out.begin()),
        "every chain's slab replays that chain's own recorded test fits");
  bool guardHeld = true;
  for (size_t i = live; i < out.size(); ++i) guardHeld &= out[i] == poison;
  check(guardHeld, "no reader writes past the retained-draw destination");

  // the same for the per-forest and variance readers' shape argument: both
  // take the destination stride from the same pair, so pin that they agree
  // with predict's on this fixture rather than only on a full store
  check(sampler.savedSlotForDraw(numDraws - 1) == numDraws - 1,
        "a partly filled store still reads head-first");

  for (ext_rng* generator : rngs) ext_rng_destroy(generator);
  rngState = savedRngState;
  printf("ok: multi-chain partial-fill predict\n");
}

void runStateTests(ext_rng* rng) {
  testFlattenRoundTrip();
  testCategoricalFlattenBoundaries();
  testKeepTrees(rng);
  testSavedDrawOrder(rng);
  testPredictCurrentTrees(rng);
  testStateRoundTrip();
  testStateRoundTripScaledOffset();
  testStateRoundTripLatents(rng);
  testStateRoundTripStudentT(rng);
  testStateValidation(rng);
  testInteractionContainment();
  testBlockAdditiveConfinement();
  testCrossGridWarmStart();
  testVarianceWarmStart();
  testVarianceWarmStartSlot();
  testVarianceSavedTreeState();
  testStateLeafScale(rng);
  testWeightsDigest();
  testMultiChainPartialFillPredict();
}
