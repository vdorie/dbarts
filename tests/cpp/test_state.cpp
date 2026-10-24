#include "common.hpp"

#include <stdexcept>

static void testFlattenRoundTrip() {
  // one categorical column (codes 0..3) and one ordinal, a hand-built tree
  // with one rule of each kind
  const size_t n = 12;
  std::vector<double> x(n * 2), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 4);
  for (size_t i = 0; i < n; ++i) x[i + n] = static_cast<double>(i) / (n - 1.0);
  ColumnKind types[] = {ColumnKind::categorical, ColumnKind::numeric};

  ColumnStore store;
  built(store.build(x.data(), n, 2, 10, false, types));

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
    built(store.build(x.data(), n, 1, 10, false, types));
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
  // creation; the state records the transform it was read under. A host that
  // recorded the moved transform re-creates the sampler in it and the state
  // continues the chain; one that did not keeps its own, and the state goes
  // in as stored all the same, its numbers read against that transform
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> offset(n);
  for (size_t i = 0; i < n; ++i)
    offset[i] = 2.0 * std::sin(0.1 * (double) i);  // widens y - offset

  SamplerOptions options;
  options.numTrees = 25;
  // a drawn k, so the state carries one
  options.updateK = true;

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
  // trajectory that produced the state; its record of the transform can
  restored.setOffset(offset.data(), false);
  // the record moves the chains at once: no install does
  restored.setAnchor(state.chains[0].fitMin, state.chains[0].fitMax);
  {
    ForestCalibration source = original.forestCalibration(0, 0);
    ForestCalibration moved = restored.forestCalibration(0, 0);
    check(moved.responseScale == source.responseScale &&
            moved.responseShift == source.responseShift &&
            moved.priorScale == source.priorScale,
          "scaled state: the record puts the chains in the transform before "
          "any state goes in");
  }
  check(restoresExactly(restored, state), "scaled state restores as stored");
  // the state holds sigma on the internal scale, as the chain does
  check(state.chains[0].sigma * original.fitScale() == original.sigma(0) &&
          state.chains[0].sigma * restored.fitScale() == restored.sigma(0),
        "scaled state: sigma is stored and installed on the internal scale");

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

  // without the record the sampler keeps the transform its own creation
  // derived, and the state goes in as stored: no value is converted, so the
  // state read back is the installed one bit for bit under this sampler's
  // pair, and the function is the stored one read on this sampler's scale
  ext_rng* rngC = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngC, 79);
  ConstantLeafSampler converted(x.data(), y.data(), n, 2, nullptr, nullptr,
                                ResponseFamily::gaussian, 1.0, 3.0,
                                0.37804942330213542, options, &rngC);
  converted.setOffset(offset.data(), false);
  double ownMin, ownMax;
  converted.getAnchor(ownMin, ownMax);
  check(ownMin != state.chains[0].fitMin && ownMax != state.chains[0].fitMax,
        "state from another transform: the two transforms differ");
  check(restoresExactly(converted, state),
        "state from another transform: installs and reports no change");
  SamplerStateData readBack;
  converted.getState(readBack);
  double anchorMin, anchorMax;
  converted.getAnchor(anchorMin, anchorMax);
  check(anchorMin == ownMin && anchorMax == ownMax &&
          readBack.chains[0].fitMin == ownMin &&
          readBack.chains[0].fitMax == ownMax,
        "state from another transform: the sampler keeps its own transform");
  SamplerStateData relabeled(state);
  relabeled.chains[0].fitMin = ownMin;
  relabeled.chains[0].fitMax = ownMax;
  check(statesAgree(relabeled, readBack) &&
          !std::isnan(state.chains[0].forests[0].k) &&
          std::memcmp(&readBack.chains[0].sigma, &state.chains[0].sigma,
                      sizeof(double)) == 0 &&
          std::memcmp(&readBack.chains[0].forests[0].k,
                      &state.chains[0].forests[0].k, sizeof(double)) == 0,
        "state from another transform: every tree, leaf value, k and sigma "
        "reads back as stored, bit for bit");
  // the same internal fits under two transforms: equal once each is taken
  // back through its own
  std::vector<double> predictionsC(20);
  converted.predict(xTest.data(), 20, 1, predictionsC.data());
  ForestCalibration from = original.forestCalibration(0, 0);
  ForestCalibration to = converted.forestCalibration(0, 0);
  double worst = 0.0;
  for (size_t i = 0; i < 20; ++i)
    worst = std::max(
      worst,
      std::fabs((predictionsC[i] - to.responseShift) / to.responseScale -
                (predictionsA[i] - from.responseShift) / from.responseScale));
  check(worst < 1.0e-12 && predictionsC != predictionsA &&
          to.responseScale != from.responseScale,
        "state from another transform: the stored fits are read on this "
        "sampler's scale, not converted");
  check(converted.sigma(0) == state.chains[0].sigma * converted.fitScale() &&
          converted.sigma(0) != original.sigma(0),
        "state from another transform: sigma goes in on the internal scale");

  // an install writes the internal sigma it is handed, with no pass through
  // response units: pinned on a value the multiplier does not carry there
  // and back, within one transform and across two. Whether a value is
  // carried depends on its significand and the multiplier's, so the search
  // sweeps a whole binade.
  {
    double range = restored.fitScale();
    double sigma = state.chains[0].sigma;
    bool found = false;
    for (int step = 1; step < 8192 && !found; ++step) {
      sigma = state.chains[0].sigma * (1.0 + step / 4096.0);
      found = (sigma * range) / range != sigma;
    }
    check(found, "sigma round trip: a value the multiplier does not carry");
    SamplerStateData probe(state);
    probe.chains[0].sigma = sigma;
    SamplerStateData within, across;
    check(restoresExactly(restored, probe) && restoresExactly(converted, probe),
          "sigma round trip: the probing states install");
    restored.getState(within);
    converted.getState(across);
    check(std::memcmp(&within.chains[0].sigma, &sigma, sizeof(double)) == 0 &&
            restored.sigma(0) == sigma * range,
          "sigma round trip: an install within one transform leaves the "
          "internal sigma bit for bit");
    check(std::memcmp(&across.chains[0].sigma, &sigma, sizeof(double)) == 0 &&
            converted.sigma(0) == sigma * converted.fitScale(),
          "sigma round trip: an install across transforms leaves the "
          "internal sigma bit for bit");
  }

  // a warm start carries the donor's trees, k and sigma the same way: as
  // stored, read against the recipient's transform
  {
    ext_rng* rngW = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngW, 82);
    ConstantLeafSampler warm(x.data(), y.data(), n, 2, nullptr, nullptr,
                             ResponseFamily::gaussian, 1.0, 3.0,
                             0.37804942330213542, options, &rngW);
    warm.setOffset(offset.data(), false);
    const std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
    check(warm.installForests(state, liveMap, x.data()) == WarmStartResult::ok,
          "warm start from another transform: installs");
    SamplerStateData warmed;
    warm.getState(warmed);
    double warmMin, warmMax;
    warm.getAnchor(warmMin, warmMax);
    check(warmMin == ownMin && warmMax == ownMax &&
            warmed.chains[0].fitMin == ownMin &&
            warmed.chains[0].fitMax == ownMax,
          "warm start from another transform: the sampler keeps its own "
          "transform");
    check(sameFlatTrees(warmed.chains[0].forests[0].trees,
                        state.chains[0].forests[0].trees) &&
            std::memcmp(&warmed.chains[0].sigma, &state.chains[0].sigma,
                        sizeof(double)) == 0 &&
            std::memcmp(&warmed.chains[0].forests[0].k,
                        &state.chains[0].forests[0].k, sizeof(double)) == 0,
          "warm start from another transform: trees, leaf values, k and "
          "sigma go in as stored, bit for bit");
    ext_rng_destroy(rngW);
  }

  // unequal pairs naming the same units: a constant response's (0, 0), the
  // window centered on 0, and a response spanning exactly (-0.5, 0.5) both
  // take multiplier 1 and shift 0, and the install is the stored chain as
  // under any other pair
  std::vector<double> yConstant(n, 0.0), yUnit(n);
  for (size_t i = 0; i < n; ++i)
    yUnit[i] = static_cast<double>(i % 5) / 4.0 - 0.5;
  ext_rng* rngD = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngE = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngD, 80);
  ext_rng_setSeed(rngE, 81);
  ConstantLeafSampler constant(x.data(), yConstant.data(), n, 2, nullptr,
                               nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                               0.37804942330213542, options, &rngD);
  ConstantLeafSampler unit(x.data(), yUnit.data(), n, 2, nullptr, nullptr,
                           ResponseFamily::gaussian, 1.0, 3.0,
                           0.37804942330213542, options, &rngE);
  constant.run(20, 0, empty);
  SamplerStateData constantState, unitState;
  constant.getState(constantState);
  unit.getAnchor(ownMin, ownMax);
  bool pairsDiffer = constantState.chains[0].fitMin == 0.0 &&
    constantState.chains[0].fitMax == 0.0 && ownMin == -0.5 && ownMax == 0.5;
  bool sameUnits = restoresExactly(unit, constantState);
  unit.getState(unitState);
  const auto& storedTrees = constantState.chains[0].forests[0].trees;
  const auto& installedTrees = unitState.chains[0].forests[0].trees;
  bool bitwise = storedTrees.size() == installedTrees.size();
  for (size_t t = 0; bitwise && t < storedTrees.size(); ++t) {
    bitwise = storedTrees[t].size() == installedTrees[t].size();
    for (size_t i = 0; bitwise && i < storedTrees[t].size(); ++i)
      bitwise = storedTrees[t][i].variable == installedTrees[t][i].variable &&
        std::memcmp(&storedTrees[t][i].value, &installedTrees[t][i].value,
                    sizeof(double)) == 0;
  }
  check(pairsDiffer && sameUnits && bitwise,
        "state from another transform: another pair naming the same units "
        "installs as stored, bit for bit");

  ext_rng_destroy(rngE);
  ext_rng_destroy(rngD);
  ext_rng_destroy(rngC);
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: state round-trip with a moved scale (as stored %.2e)\n", worst);
}

static void testSetAnchorCarriesSigmaAndVarianceCalibration() {
  // a re-created sampler is moved to the recorded mapping before any state
  // goes in; on every chain the held sigma, the variance forest's calibration
  // and its surface are those of the recorded mapping, not the creation
  // mapping's
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  const double heldSigma = 0.7;
  const double sigmaDf = 3.0, rawScale = 0.37804942330213542;

  ext_rng* rngsA[2];
  ext_rng* rngsB[2];
  for (size_t c = 0; c < 2; ++c) {
    rngsA[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngsA[c], 81 + static_cast<unsigned>(c));
    rngsB[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngsB[c], 91 + static_cast<unsigned>(c));
  }
  SamplerOptions held;
  held.numTrees = 20;
  held.numChains = 2;
  held.sigmaIsFixed = true;
  ConstantLeafSampler fixedSigma(x.data(), y.data(), n, 2, nullptr, nullptr,
                                 ResponseFamily::gaussian, heldSigma, sigmaDf,
                                 rawScale, held, rngsA);
  double min, max;
  fixedSigma.getAnchor(min, max);
  fixedSigma.setAnchor(min - 1.0, max + 2.0);
  check(std::fabs(fixedSigma.sigma(0) - heldSigma) <= 1.0e-12 * heldSigma &&
          std::fabs(fixedSigma.sigma(1) - heldSigma) <= 1.0e-12 * heldSigma,
        "setAnchor: a held sigma keeps its original-scale value on every "
        "chain");

  SamplerOptions withVariance;
  withVariance.numTrees = 20;
  withVariance.numChains = 2;
  withVariance.numVarianceTrees = 5;
  ConstantLeafSampler variance(x.data(), y.data(), n, 2, nullptr, nullptr,
                               ResponseFamily::gaussian, 1.0, sigmaDf, rawScale,
                               withVariance, rngsB);
  double scaleBefore[2];
  std::vector<double> surfaceBefore[2];
  for (size_t c = 0; c < 2; ++c) {
    scaleBefore[c] = TestPeer::varianceLeaf(variance.chain(c)).scale;
    surfaceBefore[c].resize(n);
    check(variance.currentVarianceFits(c, false, surfaceBefore[c].data()),
          "setAnchor: the variance surface reads before the move");
  }
  variance.getAnchor(min, max);
  variance.setAnchor(min - 1.0, max + 2.0);
  for (size_t c = 0; c < 2; ++c) {
    double working = 1.0 / variance.chain(c).sigmaScale();
    ConstantVarianceLeaf expected = ConstantVarianceLeaf::calibrated(
      sigmaDf, working * working * rawScale, withVariance.numVarianceTrees);
    const ConstantVarianceLeaf& leaf =
      TestPeer::varianceLeaf(variance.chain(c));
    check(leaf.scale != scaleBefore[c] &&
            leaf.degreesOfFreedom == expected.degreesOfFreedom &&
            leaf.scale == expected.scale,
          "setAnchor: the variance forest's calibration is the recorded "
          "mapping's on every chain");
    // the surface is an original-scale quantity, which the move keeps
    std::vector<double> surfaceAfter(n);
    check(variance.currentVarianceFits(c, false, surfaceAfter.data()),
          "setAnchor: the variance surface reads after the move");
    double worst = 0.0;
    for (size_t i = 0; i < n; ++i)
      worst = std::max(worst, std::fabs(surfaceAfter[i] / surfaceBefore[c][i] -
                                        1.0));
    check(worst < 1.0e-13,
          "setAnchor: the variance surface keeps its original-scale value");
  }
  for (size_t c = 0; c < 2; ++c) {
    ext_rng_destroy(rngsB[c]);
    ext_rng_destroy(rngsA[c]);
  }
  printf("ok: setAnchor carries a held sigma and the variance forest\n");
}

static void testInstallLeavesHeldSigma() {
  // a state that carries a drawn sigma goes into a chain that holds sigma
  // fixed without moving it
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  const double rawScale = 0.37804942330213542;
  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngA, 83);
  ext_rng_setSeed(rngB, 84);
  SamplerOptions drawn;
  drawn.numTrees = 20;
  ConstantLeafSampler donor(x.data(), y.data(), n, 2, nullptr, nullptr,
                            ResponseFamily::gaussian, 1.0, 3.0, rawScale, drawn,
                            &rngA);
  Results empty;
  donor.run(20, 0, empty);
  SamplerStateData state;
  donor.getState(state);
  check(!std::isnan(state.chains[0].sigma) &&
          std::fabs(donor.sigma(0) - 0.7) > 1.0e-3,
        "held sigma install: the donor's state carries a drawn sigma");

  SamplerOptions held;
  held.numTrees = 20;
  held.sigmaIsFixed = true;
  ConstantLeafSampler recipient(x.data(), y.data(), n, 2, nullptr, nullptr,
                                ResponseFamily::gaussian, 0.7, 3.0, rawScale,
                                held, &rngB);
  double before = recipient.sigma(0);
  bool altered = false;
  check(recipient.setState(state, nullptr, nullptr, &altered),
        "held sigma install: the state installs");
  check(recipient.sigma(0) == before,
        "held sigma install: a held sigma is left alone");
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: a state install leaves a held sigma alone\n");
}

// A state stored under another response shift installs on every leaf model,
// those with no mean term to carry a shift included: a gp leaf's saved draws
// and forests coupled through amplitudes go in as stored, like any other.
static void testShiftedStateInstalls() {
  std::uint64_t savedRngState = rngState;
  rngState = 424242u;
  const size_t n = 150, p = 2, numKept = 3;
  std::vector<double> x(n * p), y(n), shifted(n), z(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    y[i] = std::sin(3.0 * x[i]) + x[i + n] + z[i] * (0.5 + x[i + n]) +
           0.2 * (runif01() - 0.5);
    shifted[i] = y[i] + 5.0;
  }
  const double rawScale = 0.37804942330213542;
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, seed);
    rngs.push_back(rng);
    return rng;
  };
  auto checkInstalls = [&](auto& donor, auto& recipient, const char* what) {
    Results empty;
    donor.run(10, numKept, empty);
    SamplerStateData state, readBack;
    donor.getState(state);
    double min, max;
    recipient.getAnchor(min, max);
    std::string label(what);
    check(min != state.chains[0].fitMin && max != state.chains[0].fitMax &&
            std::fabs((max - min) -
                      (state.chains[0].fitMax - state.chains[0].fitMin)) <
              1.0e-12,
          (label + ": the recipient's transform is the donor's, shifted")
            .c_str());
    check(restoresExactly(recipient, state),
          (label + ": a state under another shift installs").c_str());
    recipient.getState(readBack);
    SamplerStateData relabeled(state);
    relabeled.chains[0].fitMin = min;
    relabeled.chains[0].fitMax = max;
    check(statesAgree(relabeled, readBack) &&
            !state.chains[0].forests[0].savedTrees.empty(),
          (label + ": the trees, fits and saved draws read back as stored")
            .c_str());
  };
  {
    const size_t covariates[] = {1};
    SamplerOptions options;
    options.numTrees = 8;
    options.gpLeaves = true;
    options.leafCovariateColumns = covariates;
    options.numLeafCovariates = 1;
    options.keepTrees = true;
    options.numSamplesToStore = numKept;
    ext_rng* rngA = newRng(4301u);
    ext_rng* rngB = newRng(4302u);
    Sampler<GPGaussianLeaf> donor(x.data(), y.data(), n, p, nullptr, nullptr,
                                  ResponseFamily::gaussian, 1.0, 3.0, rawScale,
                                  options, &rngA);
    options.leafCovariateColumns = covariates;
    Sampler<GPGaussianLeaf> recipient(x.data(), shifted.data(), n, p, nullptr,
                                      nullptr, ResponseFamily::gaussian, 1.0,
                                      3.0, rawScale, options, &rngB);
    checkInstalls(donor, recipient, "shifted state, gp leaf");
  }
  {
    SamplerOptions options;
    options.keepTrees = true;
    options.numSamplesToStore = numKept;
    AmplitudeSpec spec;
    spec.mu.numTrees = 10;
    spec.tau.numTrees = 6;
    spec.z = z.data();
    ext_rng* rngA = newRng(4303u);
    ext_rng* rngB = newRng(4304u);
    ConstantLeafSampler donor(x.data(), y.data(), n, p, nullptr, nullptr, 1.0,
                              3.0, rawScale, options, spec, &rngA);
    ConstantLeafSampler recipient(x.data(), shifted.data(), n, p, nullptr,
                                  nullptr, 1.0, 3.0, rawScale, options, spec,
                                  &rngB);
    checkInstalls(donor, recipient, "shifted state, amplitude forests");
  }
  for (ext_rng* rng : rngs) ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: a state under another response shift installs as stored\n");
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
    // a fixed nu is model and the state holds no block for it
    check(estimated ? state.chains[0].residualDf > 0.0
                    : std::isnan(state.chains[0].residualDf),
          "a drawn-nu t state carries a positive df, a fixed-nu one none");

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
    check(std::isnan(state.chains[0].residualDf) ||
            reState.chains[0].residualDf == state.chains[0].residualDf,
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
    check(nuA[window - 1] ==
            (estimated ? postState.chains[0].residualDf : residualDf),
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
  std::vector<double> treeFitsBefore(TestPeer::treeFits(sampler.chain(0)));
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
  check(TestPeer::treeFits(sampler.chain(0)) == treeFitsBefore,
        "a rejected state leaves the fits untouched");

  check(sampler.setState(state, nullptr), "the original state still restores");
  check(TestPeer::treeFits(sampler.chain(0)) == treeFitsBefore,
        "restoring the current state is an identity");

  printf("ok: state validation\n");
}

/// A latent block of precisions - Student-t scales, the Polya-Gamma variates
/// of a logistic or negative-binomial sampler - is refused when one value is
/// not positive and finite, and the refusal leaves the sampler holding what it
/// held; the smallest positive value installs. A probit sampler's latents are
/// responses and take any sign. Draws no value from the shared generator.
static void testStateLatentFloor() {
  const size_t n = 60;
  std::vector<double> x(n * 2), real(n), binary(n), count(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>((i * 37) % n) / static_cast<double>(n);
    x[n + i] = static_cast<double>((i * 11) % n) / static_cast<double>(n);
    real[i] = x[i] + 0.2 * std::sin(static_cast<double>(i));
    binary[i] = (x[i] > 0.5) != (i % 7 == 0) ? 1.0 : 0.0;
    count[i] = std::floor(4.0 * x[i]) + static_cast<double>(i % 3);
  }
  struct Case {
    const char* name;
    ResponseFamily family;
    const double* y;
    double residualDf, shape;
    bool precisions;
  };
  const double absent = std::numeric_limits<double>::quiet_NaN();
  const Case cases[] = {
    {"student-t", ResponseFamily::gaussian, real.data(), 4.0, absent, true},
    {"student-t, df drawn", ResponseFamily::gaussian, real.data(), 0.0, absent,
     true},
    {"logistic", ResponseFamily::logistic, binary.data(), absent, absent, true},
    {"nbinom", ResponseFamily::nbinom, count.data(), absent, 3.0, true},
    {"nbinom, shape drawn", ResponseFamily::nbinom, count.data(), absent, -1.0,
     true},
    {"probit", ResponseFamily::probit, binary.data(), absent, absent, false}};
  const double refused[] = {0.0, -1.0, absent,
                            std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity()};
  for (const Case& c : cases) {
    std::string name(c.name);
    SamplerOptions options;
    options.numTrees = 5;
    options.residualDf = c.residualDf;
    options.shape = c.shape;
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 20261007u);
    ConstantLeafSampler sampler(x.data(), c.y, n, 2, nullptr, nullptr,
                                c.family, 1.0, 3.0, 0.37804942330213542,
                                options, &rng);
    Results empty;
    sampler.run(10, 0, empty);
    SamplerStateData state;
    sampler.getState(state);
    check(state.chains[0].latents.size() == n,
          (name + ": the state carries a latent per row").c_str());
    if (state.chains[0].latents.size() != n) {
      ext_rng_destroy(rng);
      continue;
    }

    if (!c.precisions) {
      SamplerStateData negative(state);
      negative.chains[0].latents[n / 2] = -2.5;
      check(sampler.setState(negative, nullptr),
            (name + ": a negative latent installs").c_str());
      SamplerStateData held;
      sampler.getState(held);
      check(held.chains[0].latents == negative.chains[0].latents,
            (name + ": and is the latent the sampler holds").c_str());
      ext_rng_destroy(rng);
      continue;
    }

    bool positive = true;
    for (double value : state.chains[0].latents)
      if (!(value > 0.0) || !std::isfinite(value)) positive = false;
    check(positive, (name + ": its own precisions are positive").c_str());

    // every row is asked, the first and the last included
    const size_t rows[] = {0, n / 2, n - 1};
    for (size_t row : rows) {
      for (double value : refused) {
        SamplerStateData bad(state);
        bad.chains[0].latents[row] = value;
        check(!sampler.setState(bad, nullptr),
              (name + ": a precision that is not positive and finite is "
                      "refused").c_str());
        SamplerStateData after;
        sampler.getState(after);
        check(statesAgree(state, after),
              (name + ": and the refusal leaves the state as it was").c_str());
      }
    }
    check(sampler.setState(state, nullptr),
          (name + ": its own state installs after the refusals").c_str());

    SamplerStateData tiny(state);
    tiny.chains[0].latents[n / 2] = 1.0e-300;
    check(sampler.setState(tiny, nullptr),
          (name + ": the smallest positive precision installs").c_str());
    SamplerStateData held;
    sampler.getState(held);
    check(held.chains[0].latents == tiny.chains[0].latents,
          (name + ": and is the precision the sampler holds").c_str());
    ext_rng_destroy(rng);
  }
  printf("ok: state latent floor\n");
}

// A column whose range is one value (constant, or narrower than its uniform
// spacing) gets a grid of equal cuts, and an infinite value must not stretch
// the uniform grid: the store's own grids restore, while NaN or an over-long
// grid is refused.
static void testDegenerateGridRestores(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  for (size_t i = 0; i < n; ++i) x[i + n] = 0.0;  // constant second column
  x[0] = std::numeric_limits<double>::infinity();
  std::unique_ptr<ConstantLeafSampler> samplerPtr =
    makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  const std::vector<double>& cuts0 = sampler.data().cutPoints[0];
  const std::vector<double>& cuts1 = sampler.data().cutPoints[1];
  check(std::isfinite(cuts0.front()) && std::isfinite(cuts0.back()) &&
          cuts0.back() < 1.0,
        "the uniform grid spans the finite values only");
  check(cuts1.size() == 1 && cuts1[0] == 0.0 &&
          sampler.data().numCuts[1] == 1,
        "a constant column's uniform grid is one point at its value");

  SamplerStateData state;
  sampler.getState(state);
  check(restoresExactly(sampler, state),
        "a state carrying a one-point grid restores as stored");

  SamplerStateData bad(state);
  bad.cutPoints[0][1] = std::numeric_limits<double>::quiet_NaN();
  check(!sampler.setState(bad, nullptr), "setState rejects a NaN cut");
  bad = state;
  bad.cutPoints[0][0] = std::numeric_limits<double>::quiet_NaN();
  check(!sampler.setState(bad, nullptr), "setState rejects a leading NaN cut");
  bad = state;
  bad.cutPoints[0].assign(maxNumCutsRepresentable + 1u, 0.5);
  check(!sampler.setState(bad, nullptr),
        "setState rejects a grid past the representable count");
  check(sampler.setState(state, nullptr), "the original state still restores");

  const double repeated[] = {0.0, 0.5, 0.5, 1.0};
  check(cutGridIsValid(cuts1.data(), cuts1.size()) &&
          !cutGridIsValid(repeated, 4) && cutGridIsValid(repeated, 2),
        "a grid is valid only where it strictly increases");
  printf("ok: degenerate grid restores\n");
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
  built(store.build(x.data(), n, p, 100));  // the sampler's default cut grid
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

  // whether every live tree of a sampler satisfies a constraint, judged on
  // the independent store
  auto liveTreesSatisfy = [&](ConstantLeafSampler& s,
                              const InteractionConstraint& constraint) {
    SamplerStateData live;
    s.getState(live);
    bool satisfied = true;
    for (const std::vector<FlatNode>& flat : live.chains[0].forests[0].trees) {
      scratch.initialize(idx.data(), n);
      scratch.setInteractionConstraint(&constraint);
      satisfied &= scratch.buildFromFlat(store, flat.data(), flat.size(),
                                         params) &&
                   scratch.interactionSubtreeIsValid(0);
    }
    scratch.setInteractionConstraint(nullptr);
    return satisfied;
  };

  // (2) a max.order = 1 target declines the donor's state without force and
  // collapses its violating splits with it, and a warm start does the same
  auto k1Target = makeSampler(777, 1, nullptr);
  check(declinesUntouched(*k1Target, donorState),
        "containment: setState without force declines an "
        "interaction-violating state, the sampler untouched");
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  SamplerStateData cutBack, repaired;
  check(warmStartsAltered(*k1Target, donorState, liveMap, x.data()) &&
          liveTreesSatisfy(*k1Target, k1),
        "containment: a warm start declines an interaction-violating donor "
        "and, forced, cuts every tree back into the constraint");
  k1Target->getState(cutBack);

  check(restoresAltered(*k1Target, donorState) &&
          liveTreesSatisfy(*k1Target, k1),
        "containment: forced, every tree is collapsed into the constraint "
        "and the install reports it");
  k1Target->getState(repaired);
  check(sameFlatTrees(cutBack.chains[0].forests[0].trees,
                      repaired.chains[0].forests[0].trees),
        "containment: the warm start's trees are the forced state install's");
  check(installsClean(*k1Target, repaired),
        "containment: the repaired sampler's own state is clean");

  // (3) a forbid(x0, x1) target does likewise (a distinct predicate)
  size_t pair[] = {0, 1};
  auto forbidTarget = makeSampler(888, 0, pair);
  InteractionConstraint forbid01;
  forbid01.build(p, 0, pair, 1);
  check(declinesUntouched(*forbidTarget, donorState),
        "containment: setState without force declines a forbidden-pair "
        "violation");
  check(warmStartsAltered(*forbidTarget, donorState, liveMap, x.data()) &&
          liveTreesSatisfy(*forbidTarget, forbid01),
        "containment: a warm start declines a forbidden-pair violation and, "
        "forced, leaves no forbidden pair on a path");
  check(restoresAltered(*forbidTarget, donorState) &&
          liveTreesSatisfy(*forbidTarget, forbid01),
        "containment: forced, no tree keeps a forbidden pair on a path");

  // (4) specificity: a same-constraint donor's trees are feasible and install
  auto k1Donor = makeSampler(999, 1, nullptr);
  k1Donor->run(200, 0, empty);
  SamplerStateData k1DonorState;
  k1Donor->getState(k1DonorState);
  auto k1Target2 = makeSampler(111, 1, nullptr);
  check(warmStartsClean(*k1Target2, k1DonorState, liveMap, x.data()),
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
  built(store.build(x.data(), n, p, 100));  // the sampler's default cut grid
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

  // each live tree's columnMask_ is its block row, so the one column-mask
  // rule holds a state to its blocks: without force the out-of-block state is
  // declined, with it each tree is collapsed into its own group; a warm start
  // is declined too and, forced, keeps by rank the trees each block holds
  // within its columns
  auto target = makeSampler(999, true);
  check(declinesUntouched(*target, donorState),
        "blocks: setState without force declines an out-of-block state");
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  bool warmStarts = warmStartsAltered(*target, donorState, liveMap, x.data());
  SamplerStateData cutBack, collapsedState;
  target->getState(cutBack);
  bool installs = restoresAltered(*target, donorState);
  target->getState(collapsedState);
  auto confinedTrees = [&](const std::vector<std::vector<FlatNode>>& trees) {
    bool confined = true;
    for (size_t t = 0; t < trees.size(); ++t) {
      scratch.initialize(idx.data(), n);
      scratch.setColumnMask(maskGroup[groupOfTree(t)].data());
      confined &= scratch.buildFromFlat(store, trees[t].data(), trees[t].size(),
                                        params) &&
                  scratch.columnMaskSubtreeIsValid(0);
    }
    scratch.setColumnMask(nullptr);
    return confined;
  };
  // the warm start, forced, chooses block by block (Sampler::installForests)
  check(warmStarts && confinedTrees(cutBack.chains[0].forests[0].trees),
        "blocks: a warm start declines an out-of-block donor and, forced, "
        "takes trees each within its block");
  check(installs && confinedTrees(collapsedState.chains[0].forests[0].trees),
        "blocks: forced, every tree is collapsed into its own group");

  // (c) specificity: a same-blocks donor's trees are feasible and install cleanly
  auto compliantDonor = makeSampler(1234, true);
  compliantDonor->run(200, 0, empty);
  SamplerStateData compliantState;
  compliantDonor->getState(compliantState);
  auto target2 = makeSampler(4321, true);
  check(warmStartsClean(*target2, compliantState, liveMap, x.data()),
        "blocks: a same-blocks donor warm-starts cleanly");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: block-additive confinement (%zu splitters, %zu donor violators)\n",
         splitters, violators);
}

// ---------------------------------------------------------------------------
// A single forest's column list: SamplerOptions::forestColumns confines the
// one mean forest, and MultinomialForestSpec::columns every category forest,
// to the listed columns. The response leans hardest on x2, the excluded
// column, so an unrestricted chain splits on it in every state read. Asserted:
// no split outside the list over 500 sweeps on the constant, linear and
// monotone leaves, under blocks (each tree within its group as well) and on
// every category forest; the chain is bit for bit the chain on the listed
// columns alone, under DART too, whose excluded columns report probability 0;
// an out-of-list donor is refused by both install entries; and an empty list,
// or one naming every column, is bit for bit the chain without one.
// RNG-neutral (saves/restores rngState; local seeded generators).
// ---------------------------------------------------------------------------
template <typename S>
static size_t countSplitsOutside(S& sampler, const std::uint8_t* allowed,
                                 size_t numBlocks = 0, size_t* total = nullptr) {
  SamplerStateData state;
  sampler.getState(state);
  size_t outside = 0;
  for (const ForestStateData& forest : state.chains[0].forests)
    for (size_t t = 0; t < forest.trees.size(); ++t)
      for (const FlatNode& node : forest.trees[t]) {
        if (node.variable == invalidVariable) continue;
        if (total != nullptr) ++*total;
        // under blocks, tree t of the two equal groups may split on its own
        // group's column alone: column 0 for the first half, 1 for the second
        bool ok = numBlocks == 0
          ? allowed[node.variable] != 0
          : static_cast<size_t>(node.variable) ==
              (t < forest.trees.size() / 2 ? 0u : 1u);
        outside += ok ? 0 : 1;
      }
  return outside;
}

static void testSingleForestColumnRestriction() {
  std::uint64_t savedRngState = rngState;
  rngState = 717171u;

  const size_t n = 300, p = 3, numTrees = 20, K = 3;
  std::vector<double> x(n * p), y(n);
  std::vector<int> counts(n * K, 0), trials(n, 1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) x[i + j * n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = x[i] + 2.0 * x[i + n] + 4.0 * x[i + 2 * n] + 0.1 * z;
    counts[(x[i + 2 * n] < 0.34 ? 0 : (x[i + 2 * n] < 0.67 ? 1 : 2)) * n + i] = 1;
  }
  const std::vector<size_t> allowed = {0, 1}, everyColumn = {0, 1, 2};
  const std::uint8_t allowedMask[] = {1, 1, 0};
  const std::vector<std::int32_t> blockOfColumn = {0, 1, -1};
  const std::vector<size_t> blockTreeCounts = {numTrees / 2, numTrees / 2};
  const std::vector<size_t> covariates = {0};
  const std::int8_t directions[] = {1, 0, 0};
  const double rawScale = 0.37804942330213542;

  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  };
  auto restrictTo = [](SamplerOptions& options, const std::vector<size_t>& list,
                       size_t count) {
    options.forestColumns = list.data();
    options.numForestColumns = count;
  };
  auto make = [&]<typename L>(SamplerOptions options, std::uint32_t seed,
                              size_t numColumns = 3) {
    options.numTrees = numTrees;
    ext_rng* r = newRng(seed);
    return std::make_unique<Sampler<L>>(
      x.data(), y.data(), n, numColumns, nullptr, nullptr,
      ResponseFamily::gaussian, 1.0, 3.0, rawScale, options, &r);
  };
  Results empty;
  // every sweep's live trees are read, so a split made and undone is seen
  auto sweepsOutside = [&](auto& sampler, size_t numBlocks, size_t& total) {
    size_t outside = 0;
    for (size_t sweep = 0; sweep < 500; ++sweep) {
      sampler.run(1, 0, empty);
      outside += countSplitsOutside(sampler, allowedMask, numBlocks, &total);
    }
    return outside;
  };

  SamplerOptions plain, restricted, blocked, linear, monotone;
  restrictTo(restricted, allowed, allowed.size());
  blocked = restricted;
  blocked.numBlocks = blockTreeCounts.size();
  blocked.blockOfColumn = blockOfColumn.data();
  blocked.blockTreeCounts = blockTreeCounts.data();
  linear = restricted;
  linear.leafCovariateColumns = covariates.data();
  linear.numLeafCovariates = covariates.size();
  monotone = restricted;
  monotone.monotoneDirections = directions;

  size_t total = 0, donorTotal = 0;
  auto constant = make.operator()<ConstantGaussianLeaf>(restricted, 101);
  check(sweepsOutside(*constant, 0, total) == 0 && total > 0,
        "forest columns: a restricted constant-leaf chain stays in its list");
  auto donor = make.operator()<ConstantGaussianLeaf>(plain, 101);
  check(sweepsOutside(*donor, 0, donorTotal) > 0,
        "forest columns: the unrestricted chain splits on the excluded column");
  total = 0;
  auto inBlocks = make.operator()<ConstantGaussianLeaf>(blocked, 102);
  check(sweepsOutside(*inBlocks, 2, total) == 0 && total > 0,
        "forest columns: under blocks each tree stays in its group");
  total = 0;
  auto linearLeaf = make.operator()<LinearGaussianLeaf>(linear, 103);
  check(sweepsOutside(*linearLeaf, 0, total) == 0 && total > 0,
        "forest columns: a restricted linear-leaf chain stays in its list");
  total = 0;
  auto monotoneLeaf =
    make.operator()<MonotoneConstantGaussianLeaf>(monotone, 104);
  check(sweepsOutside(*monotoneLeaf, 0, total) == 0 && total > 0,
        "forest columns: a restricted monotone chain stays in its list");

  auto makeMultinomial = [&](bool withList, std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = numTrees;
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = counts.data();
    spec.trials = trials.data();
    spec.forest.numTrees = numTrees;
    if (withList) {
      spec.forest.columns = allowed.data();
      spec.forest.numColumns = allowed.size();
    }
    ext_rng* r = newRng(seed);
    return std::make_unique<Sampler<ConstantGaussianLeaf>>(x.data(), n, p,
                                                           options, spec, &r);
  };
  total = donorTotal = 0;
  auto categories = makeMultinomial(true, 105);
  check(sweepsOutside(*categories, 0, total) == 0 && total > 0,
        "forest columns: every category forest stays in its list");
  auto freeCategories = makeMultinomial(false, 105);
  check(sweepsOutside(*freeCategories, 0, donorTotal) > 0,
        "forest columns: unrestricted category forests split outside it");

  // the restricted chain is the chain on the listed columns alone: the list is
  // the store's first two columns, so the two states compare node for node
  auto identicalForests = [](const SamplerStateData& a,
                             const SamplerStateData& b) {
    const auto& ta = a.chains[0].forests[0].trees;
    const auto& tb = b.chains[0].forests[0].trees;
    if (ta.size() != tb.size() || a.chains[0].sigma != b.chains[0].sigma)
      return false;
    for (size_t t = 0; t < ta.size(); ++t) {
      if (ta[t].size() != tb[t].size()) return false;
      for (size_t i = 0; i < ta[t].size(); ++i)
        if (ta[t][i].variable != tb[t][i].variable ||
            ta[t][i].mask != tb[t][i].mask || ta[t][i].flags != tb[t][i].flags)
          return false;
    }
    return true;
  };
  SamplerOptions dart, dartRestricted, dartEveryColumn, emptyList;
  dart.useDart = true;
  dartRestricted = dart;
  restrictTo(dartRestricted, allowed, allowed.size());
  dartEveryColumn = dart;
  restrictTo(dartEveryColumn, everyColumn, everyColumn.size());
  restrictTo(emptyList, allowed, 0);
  auto stateAfter = [&]<typename L>(const SamplerOptions& options,
                                    size_t numColumns) {
    auto sampler = make.operator()<L>(options, 106, numColumns);
    sampler->run(100, 0, empty);
    SamplerStateData state;
    sampler->getState(state);
    return state;
  };
  using Constant = ConstantGaussianLeaf;
  SamplerStateData onList = stateAfter.operator()<Constant>(restricted, 3);
  SamplerStateData onSubMatrix = stateAfter.operator()<Constant>(plain, 2);
  check(identicalForests(onList, onSubMatrix),
        "forest columns: the restricted chain is the chain on its columns");
  SamplerStateData dartOnList =
    stateAfter.operator()<Constant>(dartRestricted, 3);
  SamplerStateData dartOnSubMatrix = stateAfter.operator()<Constant>(dart, 2);
  const std::vector<double>& probabilities(dartOnList.chains[0].dartProbabilities);
  check(identicalForests(dartOnList, dartOnSubMatrix) &&
          probabilities.size() == p &&
          probabilities[0] == dartOnSubMatrix.chains[0].dartProbabilities[0] &&
          probabilities[1] == dartOnSubMatrix.chains[0].dartProbabilities[1] &&
          dartOnList.chains[0].dartAlpha == dartOnSubMatrix.chains[0].dartAlpha,
        "forest columns: restricted DART is DART on its columns");
  check(probabilities[2] == 0.0 && probabilities[0] > 0.0 &&
          probabilities[1] > 0.0,
        "forest columns: DART reports 0 for an excluded column");
  SamplerStateData unlisted = stateAfter.operator()<Constant>(plain, 3);
  SamplerStateData dartUnlisted = stateAfter.operator()<Constant>(dart, 3);
  check(statesAgree(unlisted, stateAfter.operator()<Constant>(emptyList, 3)) &&
          identicalForests(unlisted,
                           stateAfter.operator()<Constant>(emptyList, 3)),
        "forest columns: an empty list is the chain without one");
  SamplerStateData dartListed =
    stateAfter.operator()<Constant>(dartEveryColumn, 3);
  check(identicalForests(dartUnlisted, dartListed) &&
          dartUnlisted.chains[0].dartProbabilities ==
            dartListed.chains[0].dartProbabilities,
        "forest columns: DART over every column is DART without a list");

  // a state is installed as it is given, a probability on an excluded column
  // included; the next Dirichlet draw puts that column back at 0
  auto dartChain = make.operator()<Constant>(dartRestricted, 107);
  dartChain->run(20, 0, empty);
  SamplerStateData carried;
  dartChain->getState(carried);
  carried.chains[0].dartProbabilities = {0.25, 0.25, 0.5};
  bool carriedInstalls = dartChain->setState(carried, nullptr);
  SamplerStateData asInstalled, afterDraw;
  dartChain->getState(asInstalled);
  dartChain->run(1, 0, empty);
  dartChain->getState(afterDraw);
  const std::vector<double>& redrawn(afterDraw.chains[0].dartProbabilities);
  check(carriedInstalls &&
          asInstalled.chains[0].dartProbabilities[2] == 0.5 &&
          redrawn[2] == 0.0 && redrawn[0] > 0.0 && redrawn[1] > 0.0 &&
          std::fabs(redrawn[0] + redrawn[1] - 1.0) < 1e-12,
        "forest columns: a DART draw zeroes an excluded column a state carried");
  // a warm-start donor's DART weights, one a predictor, are its shape
  {
    const std::vector<std::pair<size_t, int>> live = {{0, -1}};
    bool refused = true;
    for (size_t length : {size_t{1}, p - 1, p + 1, size_t{100000}}) {
      SamplerStateData misfit = carried;
      misfit.chains[0].dartProbabilities.assign(length, 1.0 / length);
      WarmStartResult result;
      bool altered;
      refused = refused &&
        warmVerdictMatchesInstall(*dartChain, misfit, live, x.data(), result,
                                  altered) &&
        result == WarmStartResult::shapeMismatch;
    }
    check(refused, "forest columns: a warm-start donor whose DART weights "
                   "are not one a predictor is refused in both forms");
  }

  // an out-of-list donor is declined by a warm start without force, as its
  // state is, and the chain's own state is not
  SamplerStateData donorState, ownState;
  donor->getState(donorState);
  constant->getState(ownState);
  check(declinesUntouched(*constant, donorState),
        "forest columns: setState without force declines an out-of-list "
        "state");
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(constant->installForests(donorState, liveMap, x.data(), false) ==
          WarmStartResult::declined,
        "forest columns: a warm start without force declines an out-of-list "
        "donor");
  SamplerStateData afterRefusals;
  constant->getState(afterRefusals);
  check(statesAgree(ownState, afterRefusals) &&
          installsClean(*constant, ownState),
        "forest columns: a decline leaves the chain, whose own state restores");
  SamplerStateData freeCategoryState;
  freeCategories->getState(freeCategoryState);
  check(declinesUntouched(*categories, freeCategoryState),
        "forest columns: category forests decline an out-of-list state");
  // forced, both collapse every split outside the list
  auto offList = [&](const SamplerStateData& state) {
    for (const ForestStateData& forest : state.chains[0].forests)
      for (const std::vector<FlatNode>& tree : forest.trees)
        for (const FlatNode& node : tree)
          if (node.variable == 2) return true;
    return false;
  };
  SamplerStateData constantRepaired, categoriesRepaired;
  bool constantInstalls = restoresAltered(*constant, donorState);
  bool categoriesInstall = restoresAltered(*categories, freeCategoryState);
  constant->getState(constantRepaired);
  categories->getState(categoriesRepaired);
  check(offList(donorState) && constantInstalls && !offList(constantRepaired),
        "forest columns: forced, no split outside the list is left");
  SamplerStateData constantWarm;
  bool constantWarmStarts =
    warmStartsAltered(*constant, donorState, liveMap, x.data());
  constant->getState(constantWarm);
  check(constantWarmStarts &&
          sameFlatTrees(constantWarm.chains[0].forests[0].trees,
                        constantRepaired.chains[0].forests[0].trees),
        "forest columns: a forced warm start cuts the donor back as the "
        "forced state install does");
  check(offList(freeCategoryState) && categoriesInstall &&
          !offList(categoriesRepaired),
        "forest columns: forced, no category forest keeps one either");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: single-forest column restriction\n");
}

// Cross-grid warm start: a donor grown on a fine cut grid seeds a destination
// on a coarser one, by force, since no such donor stands as stored. The
// donor's trees are built and routed on its own points, a leaf none of the
// destination's rows reaches merging there, and the splits then move onto
// the destination's grid, one left no point merging at the mean of the
// leaves beneath by the rows each held: what a forced state install gives.
// Every installed tree stays occupied, the store's grid and codes are as
// they were, and the sampler runs.
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

  // a twin takes the donor as a state, forced and keeping its own grid
  ext_rng* twinRng = newRng(717);
  ConstantLeafSampler twin(xDest.data(), yDest.data(), n, p, nullptr, nullptr,
                           ResponseFamily::gaussian, 1.0, 3.0,
                           0.37804942330213542, destOptions, &twinRng);
  bool twinMoved = false;
  check(twin.setState(donorState, xDest.data(), nullptr, &twinMoved,
                      keepStoreCapacity, true) && twinMoved,
        "cross-grid: the twin takes the donor's state across the grids");

  // declined as it stands, installed by force, the store's grid untouched
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  const std::vector<std::vector<double>> ownCutPoints(dest.data().cutPoints);
  const std::vector<xint_t> ownCodes(dest.data().train.codes);
  check(warmStartsAltered(dest, donorState, liveMap, xDest.data()),
        "cross-grid: a different-grid donor is declined and, forced, "
        "installs across the grids");
  check(dest.data().cutPoints == ownCutPoints &&
          dest.data().train.codes == ownCodes,
        "cross-grid: the destination's cut points and codes are as they were");
  SamplerStateData warm, moved;
  dest.getState(warm);
  twin.getState(moved);
  check(sameFlatTrees(warm.chains[0].forests[0].trees,
                      moved.chains[0].forests[0].trees),
        "cross-grid: the trees are the forced state install's");

  // no leaf is left empty and no split outside its interval, and the live
  // leaf count is below the donor's (a genuine collapse fired)
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
  check(occupied, "cross-grid: the install leaves no empty leaf");
  bool inside = true;
  for (size_t t = 0; t < numTrees; ++t) {
    const Tree& tree = dest.chain(0).tree(t);
    std::vector<int32_t> subtree;
    tree.fillSubtree(0, subtree);
    for (int32_t i : subtree)
      inside &= tree.at(i).isBottom() ||
                !tree.splitIsOutsideInterval(dest.data(), i);
  }
  check(inside,
        "cross-grid: no split sits outside its interval on the shorter grid");
  check(liveLeaves < donorLeaves,
        "cross-grid: a starved split collapsed (fewer live leaves than donor)");

  // a split left no point merges by rows. Column 1 of a coarse destination
  // holds one cut; the hand tree splits it at two of the donor's points, a
  // below b, both nearest that cut. The root moves there and the split
  // beneath it, with no point left, becomes one leaf at the mean of its two
  // leaves weighted by the destination's rows each held before the move
  {
    SamplerOptions coarseOptions;
    coarseOptions.numTrees = numTrees;
    coarseOptions.maxNumCuts = 1;
    ext_rng* coarseRng = newRng(919);
    ConstantLeafSampler coarse(xDest.data(), yDest.data(), n, p, nullptr,
                               nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                               0.37804942330213542, coarseOptions, &coarseRng);
    const std::vector<double>& points(donorState.cutPoints[1]);
    double a = points[(2 * points.size()) / 5], b = points[points.size() / 2];
    SamplerStateData hand = donorState;
    std::vector<FlatNode>& tree(hand.chains[0].forests[0].trees[0]);
    tree.assign(5, FlatNode());
    tree[0].variable = tree[2].variable = 1;
    tree[0].value = a;
    tree[2].value = b;
    setFlatKind(tree[0], FlatKind::ordinal);
    setFlatKind(tree[2], FlatKind::ordinal);
    tree[1].value = 0.10;
    tree[3].value = -0.20;
    tree[4].value = 0.05;
    size_t numMiddle = 0, numRight = 0;
    for (size_t i = 0; i < n; ++i) {
      numMiddle += xDest[i + n] > a && xDest[i + n] <= b;
      numRight += xDest[i + n] > b;
    }
    double merged = (numMiddle * -0.20 + numRight * 0.05) /
                    static_cast<double>(numMiddle + numRight);
    bool installs = coarse.data().numCuts[1] == 1 && a < b &&
      warmStartsAltered(coarse, hand, liveMap, xDest.data());
    SamplerStateData after;
    coarse.getState(after);
    const std::vector<FlatNode>& got(after.chains[0].forests[0].trees[0]);
    check(installs && numMiddle > 0 && numRight > 0 && got.size() == 3 &&
            got[0].variable == 1 &&
            got[0].value == coarse.data().cutPoints[1][0] &&
            got[1].value == 0.10 && std::fabs(got[2].value - merged) < 1e-12 &&
            std::fabs(merged - -0.075) > 1e-3,
          "cross-grid: a split left no point merges at the mean of its "
          "leaves by the rows each held before the move");
  }

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

// A warm start between samplers that standardize a leaf covariate
// differently: the donor's linear coefficients go in as stored, bit for bit,
// numbers per unit of the recipient's own scaled covariate, from the donor's
// live trees or a kept draw, on the donor's grid or another, so the seeded
// forest is not the donor's function where the standardizations differ; the
// recipient's center, scale and lengthscale are never moved, and the donor's
// are not read, so a donor recording none, or a record a hand has broken,
// installs the same. A gp leaf's per-row fits hold no kernel: copied on the
// same grid, cold on another, as before. Two chains, several trees, the
// edited tree and chain last.
static void testWarmStartStandardization() {
  std::uint64_t savedRngState = rngState;
  rngState = 577215u;
  const size_t n = 150, p = 2, numChains = 2, numTrees = 8, numKept = 3;
  std::vector<double> x(n * p), xStretched(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 2.0 * runif01() - 1.0;
    xStretched[i] = x[i];
    xStretched[i + n] = 3.0 * x[i + n] + 5.0;
    y[i] = 2.0 * x[i] + (x[i] > 0.5 ? x[i + n] : 0.0) +
           0.2 * (runif01() - 0.5);
  }

  std::vector<ext_rng*> rngs;
  rngs.reserve(32 * numChains);
  auto newRngs = [&](std::uint32_t seed) {
    for (size_t c = 0; c < numChains; ++c) {
      rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
      ext_rng_setSeed(rngs.back(), seed + static_cast<std::uint32_t>(c));
    }
    return rngs.data() + rngs.size() - numChains;
  };
  const size_t covariates[] = {1};
  auto make = [&](const std::vector<double>& xs, bool gp, bool keepTrees,
                  std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numChains = numChains;
    options.maxNumCuts = 12;
    options.keepTrees = keepTrees;
    options.numSamplesToStore = keepTrees ? numKept : 0;
    options.leafCovariateColumns = covariates;
    options.numLeafCovariates = 1;
    options.gpLeaves = gp;
    return createSampler(xs.data(), y.data(), n, p, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, newRngs(seed));
  };
  // a sampler standardized on the stretched covariate that holds the donor's
  // rows: an in-place update keeps the creation-time constants
  auto makeOffStandard = [&](bool gp, std::uint32_t seed) {
    std::unique_ptr<SamplerBase> sampler = make(xStretched, gp, false, seed);
    check(sampler->setPredictor(x.data(), true, true) ==
            PredictorUpdateResult::accepted,
          "warm standardization: the recipient takes the donor's rows");
    return sampler;
  };
  // every donor cut point and the midpoints between them
  auto refineGrid = [&](SamplerBase& sampler,
                        const std::vector<std::vector<double>>& grid) {
    std::vector<std::vector<double>> fine(p);
    const double* pointers[p];
    std::uint32_t counts[p];
    size_t columns[p];
    for (size_t j = 0; j < p; ++j) {
      for (size_t k = 0; k < grid[j].size(); ++k) {
        if (k > 0) fine[j].push_back(0.5 * (grid[j][k - 1] + grid[j][k]));
        fine[j].push_back(grid[j][k]);
      }
      pointers[j] = fine[j].data();
      counts[j] = static_cast<std::uint32_t>(fine[j].size());
      columns[j] = j;
    }
    sampler.setCutPoints(pointers, counts, columns, p, x.data());
  };
  auto liveFits = [&](SamplerBase& sampler) {
    std::vector<double> fits(n * numChains);
    for (size_t c = 0; c < numChains; ++c)
      sampler.fitsWithoutOffset(c, fits.data() + c * n);
    return fits;
  };
  auto worstGap = [](const std::vector<double>& a,
                     const std::vector<double>& b) {
    double worst = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
      worst = std::max(worst, std::fabs(a[i] - b[i]));
    return worst;
  };
  auto sameConstants = [](const SamplerStateData& a,
                          const SamplerStateData& b) {
    for (size_t c = 0; c < a.chains.size(); ++c) {
      const ForestStateData& fa = a.chains[c].forests[0];
      const ForestStateData& fb = b.chains[c].forests[0];
      if (fa.leafCovariateCenters != fb.leafCovariateCenters ||
          fa.leafCovariateScales != fb.leafCovariateScales ||
          fa.leafLengthscales != fb.leafLengthscales)
        return false;
    }
    return true;
  };
  // the live trees' values bit for bit, a zero's sign included
  auto sameLiveValues = [](const SamplerStateData& a,
                           const SamplerStateData& b) {
    for (size_t c = 0; c < a.chains.size(); ++c) {
      const ForestStateData& fa = a.chains[c].forests[0];
      const ForestStateData& fb = b.chains[c].forests[0];
      if (fa.trees.size() != fb.trees.size() ||
          fa.treeParams.size() != fb.treeParams.size())
        return false;
      for (size_t t = 0; t < fa.trees.size(); ++t) {
        if (fa.trees[t].size() != fb.trees[t].size() ||
            fa.treeParams[t].size() != fb.treeParams[t].size())
          return false;
        for (size_t i = 0; i < fa.trees[t].size(); ++i)
          if (fa.trees[t][i].variable != fb.trees[t][i].variable ||
              std::memcmp(&fa.trees[t][i].value, &fb.trees[t][i].value,
                          sizeof(double)) != 0)
            return false;
        if (!fa.treeParams[t].empty() &&
            std::memcmp(fa.treeParams[t].data(), fb.treeParams[t].data(),
                        fa.treeParams[t].size() * sizeof(double)) != 0)
          return false;
      }
    }
    return true;
  };
  const std::vector<std::pair<size_t, int>> liveMap = {{0, -1}, {1, -1}};
  Results none;

  // --- linear leaves
  std::unique_ptr<SamplerBase> donor = make(x, false, true, 100);
  donor->run(80, numKept, none);
  SamplerStateData donorState;
  donor->getState(donorState);
  std::vector<double> donorFits = liveFits(*donor);

  // the donor's grid, another standardization: clean, the values as stored
  {
    std::unique_ptr<SamplerBase> recipient = makeOffStandard(false, 200);
    SamplerStateData before, after;
    recipient->getState(before);
    check(recipient->data().cutPoints == donorState.cutPoints &&
            !sameConstants(before, donorState),
          "warm standardization: the donor's grid, another standardization");
    // the store is not read: the donor's live trees seed both chains
    check(warmStartsClean(*recipient, donorState, liveMap, x.data()),
          "warm standardization: a donor under another standardization "
          "installs as it stands");
    recipient->getState(after);
    check(sameLiveValues(after, donorState),
          "warm standardization: the coefficients go in as stored");
    check(sameConstants(before, after),
          "warm standardization: the recipient keeps its standardization");
    check(worstGap(liveFits(*recipient), donorFits) > 1e-3,
          "warm standardization: so the seeded fit is not the donor's "
          "function on the recipient's rows");
    std::vector<double> fits(n * numChains * 2);
    Results results;
    results.trainingFits = fits.data();
    recipient->run(0, 2, results);
    bool finite = true;
    for (double fit : fits) finite = finite && std::isfinite(fit);
    check(finite, "warm standardization: the seeded chains draw");
  }

  // a kept draw as the seed, chains crossed: each chain's live trees are the
  // draw it was given, as stored
  {
    std::unique_ptr<SamplerBase> recipient = makeOffStandard(false, 300);
    const std::vector<std::pair<size_t, int>> slotMap = {
      {1, static_cast<int>(donor->savedSlotForDraw(numKept - 1))},
      {0, static_cast<int>(donor->savedSlotForDraw(0))}};
    check(recipient->installForests(donorState, slotMap, x.data()) ==
            WarmStartResult::ok,
          "warm standardization: a kept draw installs");
    SamplerStateData expected = donorState, after;
    for (size_t c = 0; c < numChains; ++c) {
      const ForestStateData& from =
        donorState.chains[slotMap[c].first].forests[0];
      ForestStateData& to = expected.chains[c].forests[0];
      size_t base = static_cast<size_t>(slotMap[c].second) * numTrees;
      to.trees.assign(from.savedTrees.begin() + base,
                      from.savedTrees.begin() + base + numTrees);
      to.treeParams.assign(from.savedTreeParams.begin() + base,
                           from.savedTreeParams.begin() + base + numTrees);
    }
    recipient->getState(after);
    check(sameLiveValues(after, expected),
          "warm standardization: a kept draw seeds its coefficients as "
          "stored");
  }

  // a grid that holds every donor split point and more: the values as
  // stored, and the recipient's standardization is not re-derived
  {
    std::unique_ptr<SamplerBase> recipient = makeOffStandard(false, 400);
    refineGrid(*recipient, donorState.cutPoints);
    SamplerStateData before, after;
    recipient->getState(before);
    check(recipient->data().cutPoints != donorState.cutPoints,
          "warm standardization: the recipient is on another grid");
    check(recipient->installForests(donorState, liveMap, x.data()) ==
            WarmStartResult::ok,
          "warm standardization: a donor on another grid installs");
    recipient->getState(after);
    check(sameLiveValues(after, donorState),
          "warm standardization: across grids the coefficients go in as "
          "stored");
    check(sameConstants(before, after),
          "warm standardization: across grids the recipient's standardization "
          "is bit for bit what it was");
  }

  // equal standardizations: no arithmetic, down to the sign of a zero in the
  // last chain's last tree, which adding a slope's zero term would flip
  {
    SamplerStateData edited = donorState;
    ForestStateData& lastForest = edited.chains[numChains - 1].forests[0];
    for (FlatNode& node : lastForest.trees.back())
      if (flatKindOf(node) == FlatKind::leaf) {
        node.value = -0.0;
        break;
      }
    lastForest.treeParams.back()[0] = 0.25;
    std::unique_ptr<SamplerBase> recipient = make(x, false, false, 500);
    SamplerStateData before, after;
    recipient->getState(before);
    check(sameConstants(before, donorState),
          "warm standardization: the twin standardizes as the donor does");
    check(recipient->installForests(edited, liveMap, x.data()) ==
            WarmStartResult::ok,
          "warm standardization: an equal-standardization donor installs");
    recipient->getState(after);
    check(sameLiveValues(after, edited),
          "warm standardization: between equal standardizations the installed "
          "values are identical to the donor's");
  }

  // the donor's record is not read: one recording none, and one a hand has
  // broken, install as stored with nothing of the recipient's moved
  {
    SamplerStateData unrecorded = donorState, zeroScale = donorState,
      tooLong = donorState, halfRecord = donorState;
    for (ChainStateData& chain : unrecorded.chains) {
      chain.forests[0].leafCovariateCenters.clear();
      chain.forests[0].leafCovariateScales.clear();
    }
    const size_t last = numChains - 1;
    zeroScale.chains[last].forests[0].leafCovariateScales[0] = 0.0;
    tooLong.chains[last].forests[0].leafCovariateCenters.push_back(0.0);
    tooLong.chains[last].forests[0].leafCovariateScales.push_back(1.0);
    halfRecord.chains[last].forests[0].leafCovariateScales.clear();
    std::uint32_t seed = 600;
    for (const SamplerStateData* given :
         {&unrecorded, &zeroScale, &tooLong, &halfRecord}) {
      std::unique_ptr<SamplerBase> recipient = makeOffStandard(false, seed++);
      SamplerStateData before, after;
      recipient->getState(before);
      check(recipient->installForests(*given, liveMap, x.data()) ==
              WarmStartResult::ok,
            "warm standardization: a donor's record is not read");
      recipient->getState(after);
      check(sameLiveValues(after, donorState) && sameConstants(before, after),
            "warm standardization: and its values arrive as stored");
    }
  }

  // --- gp leaves: fits copied on the donor's grid, cold on another, the
  // kernel's constants the recipient's own on both
  {
    std::unique_ptr<SamplerBase> gpDonor = make(x, true, false, 800);
    gpDonor->run(60, 0, none);
    SamplerStateData gpDonorState;
    gpDonor->getState(gpDonorState);

    std::unique_ptr<SamplerBase> recipient = makeOffStandard(true, 900);
    SamplerStateData before, after;
    recipient->getState(before);
    check(!sameConstants(before, gpDonorState),
          "warm standardization, gp: the recipient's constants differ");
    check(warmStartsClean(*recipient, gpDonorState, liveMap, x.data()),
          "warm standardization, gp: a same-grid donor installs as it "
          "stands");
    recipient->getState(after);
    check(worstGap(liveFits(*recipient), liveFits(*gpDonor)) < 1e-12,
          "warm standardization, gp: the donor's fits are copied");
    check(sameConstants(before, after),
          "warm standardization, gp: the kernel's constants are unchanged on "
          "the donor's grid");

    std::unique_ptr<SamplerBase> other = makeOffStandard(true, 1000);
    refineGrid(*other, gpDonorState.cutPoints);
    other->getState(before);
    WarmStartReport coldStart;
    check(warmStartsAltered(*other, gpDonorState, liveMap, x.data(), nullptr,
                            &coldStart) &&
            coldStart.fitsCold && !coldStart.mean[0].changed() &&
            !coldStart.mean[1].changed(),
          "warm standardization, gp: a donor recorded on other points is "
          "declined though every split stands and, forced, installs with "
          "its fits reported cold and no tree changed");
    other->getState(after);
    bool cold = true;
    std::vector<double> totals(n);
    for (size_t c = 0; c < numChains; ++c) {
      other->forestTotalFits(c, 0, totals.data());
      for (double total : totals) cold = cold && total == 0.0;
    }
    check(cold, "warm standardization, gp: fits start at zero on another grid");
    check(sameConstants(before, after),
          "warm standardization, gp: the kernel's constants are unchanged on "
          "another grid");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: warm start across leaf standardizations\n");
}

// A leaf covariate without spread, across the calls that restate live
// coefficients. Such a column sits at its center under the placeholder scale
// and reads zero on every training row, so the intercept alone is the fit:
// between two of them a live block does not move; off one or onto one the
// formula runs over the placeholder 1 and the leaf is the function it was,
// taking its value at the new constant onto one. Whether
// a column is one is read off the rows and constants in force, so a column
// that was given values, moved or made constant after creation converts as
// what it is then. Each case asserts the fitted function - the live fit on
// rows the call kept and on rows it appended, kept draws at fresh points - to
// 1e-10 relative to max(1, |value|), and that the next draws stay within the
// range of a twin made the same way and left alone, widened by twice its
// width either way. The appended rows sit inside the routing column's range
// and the response's, so its grid, every route and the response transform
// stand. Two chains, several trees.
static void testNoSpreadLiveConversion() {
  std::uint64_t savedRngState = rngState;
  rngState = 141421u;
  const size_t n = 150, n2 = 190, p = 2, numChains = 2, numTrees = 8,
    numKept = 3;
  const double tolerance = 1e-10;
  // column 0 routes; z holds the leaf covariate's real values
  std::vector<double> route(n2), z(n2), y(n2);
  double routeMin = 1.0, routeMax = 0.0, yMin = 0.0, yMax = 0.0;
  for (size_t i = 0; i < n; ++i) {
    route[i] = runif01();
    z[i] = 2.0 * runif01() - 1.0;
    y[i] = 2.0 * route[i] + (route[i] > 0.5 ? z[i] : -z[i]) +
           0.2 * (runif01() - 0.5);
    routeMin = std::min(routeMin, route[i]);
    routeMax = std::max(routeMax, route[i]);
    yMin = i == 0 ? y[i] : std::min(yMin, y[i]);
    yMax = i == 0 ? y[i] : std::max(yMax, y[i]);
  }
  for (size_t i = n; i < n2; ++i) {
    route[i] = routeMin + (routeMax - routeMin) * (0.05 + 0.9 * runif01());
    z[i] = 3.0 * runif01() - 1.5;
    y[i] = yMin + (yMax - yMin) * (0.05 + 0.9 * runif01());
  }
  // n2 rows, column-major; a sampler made on fewer reads a prefix of each
  // column through its own copy
  auto design = [&](size_t rows, double constant, double slope) {
    std::vector<double> x(rows * p);
    for (size_t i = 0; i < rows; ++i) {
      x[i] = route[i];
      x[i + rows] = constant + slope * z[i];
    }
    return x;
  };
  // rows [0, n) at the constant and the appended ones spread about it
  auto grown = [&](double constant) {
    std::vector<double> x = design(n2, constant, 0.0);
    for (size_t i = n; i < n2; ++i) x[i + n2] += z[i];
    return x;
  };

  std::vector<ext_rng*> rngs;
  rngs.reserve(64 * numChains);
  auto newRngs = [&](std::uint32_t seed) {
    for (size_t c = 0; c < numChains; ++c) {
      rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
      ext_rng_setSeed(rngs.back(), seed + static_cast<std::uint32_t>(c));
    }
    return rngs.data() + rngs.size() - numChains;
  };
  const size_t covariates[] = {1};
  // a derivation over a column of one value is refused, so a constant
  // column's center comes from a supplied record without spread, as a host's
  // carried one would: the value the placeholder centers at
  const double noScale = std::numeric_limits<double>::quiet_NaN();
  auto constantOf = [&](const std::vector<double>& x, double* value) {
    size_t rows = x.size() / p;
    *value = x[rows];
    for (size_t i = 1; i < rows; ++i)
      if (x[i + rows] != *value) return false;
    return true;
  };
  auto make = [&](const std::vector<double>& x, std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numChains = numChains;
    options.keepTrees = true;
    options.numSamplesToStore = numKept;
    options.leafCovariateColumns = covariates;
    options.numLeafCovariates = 1;
    double constant;
    if (constantOf(x, &constant))
      options.suppliedStandardization =
        SuppliedStandardization{&constant, &noScale};
    return createSampler(x.data(), y.data(), x.size() / p, p, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, newRngs(seed));
  };
  Results none;
  auto liveFits = [&](SamplerBase& sampler) {
    size_t rows = sampler.data().numObservations;
    std::vector<double> fits(rows * numChains);
    for (size_t c = 0; c < numChains; ++c)
      sampler.fitsWithoutOffset(c, fits.data() + c * rows);
    return fits;
  };
  auto keptDraws = [&](SamplerBase& sampler, const std::vector<double>& x) {
    std::vector<double> kept(x.size() / p * numKept * numChains);
    sampler.predict(x.data(), x.size() / p, 1, kept.data());
    return kept;
  };
  // the newest kept draw of each chain at the rows of x: the live function
  // while the last sweep run is the last one kept
  auto newestDraw = [&](SamplerBase& sampler, const std::vector<double>& x) {
    size_t rows = x.size() / p;
    std::vector<double> kept = keptDraws(sampler, x), newest(rows * numChains);
    for (size_t c = 0; c < numChains; ++c)
      std::copy_n(kept.begin() + rows * (numKept - 1 + numKept * c), rows,
                  newest.begin() + c * rows);
    return newest;
  };
  // rows [0, count) of each chain's block of rowsA in a against the same of b
  auto worstGap = [&](const std::vector<double>& a, size_t rowsA,
                      const std::vector<double>& b, size_t rowsB,
                      size_t first, size_t count) {
    double worst = 0.0;
    for (size_t c = 0; c < a.size() / rowsA; ++c)
      for (size_t i = first; i < first + count; ++i) {
        double u = a[i + c * rowsA], v = b[i + c * rowsB];
        worst = std::max(worst, std::fabs(u - v) / std::max(1.0, std::fabs(v)));
      }
    return worst;
  };
  auto standardization = [](SamplerBase& sampler, double& center,
                            double& scale) {
    SamplerStateData state;
    sampler.getState(state);
    center = state.chains.back().forests[0].leafCovariateCenters[0];
    scale = state.chains.back().forests[0].leafCovariateScales[0];
  };
  // the largest live slope, in size, over every chain
  auto largestSlope = [](SamplerBase& sampler) {
    SamplerStateData state;
    sampler.getState(state);
    double largest = 0.0;
    for (const ChainStateData& chain : state.chains)
      for (const std::vector<double>& params : chain.forests[0].treeParams)
        for (double slope : params)
          largest = std::max(largest, std::fabs(slope));
    return largest;
  };
  // three more draws of each: the converted sampler's stay within the twin's
  // range widened by twice its width on either side
  auto drawsLikeTwin = [&](SamplerBase& sampler, SamplerBase& twin) {
    double range[2][2];
    SamplerBase* samplers[] = {&twin, &sampler};
    for (size_t k = 0; k < 2; ++k) {
      std::vector<double> fits(samplers[k]->data().numObservations *
                               numChains * 3);
      Results results;
      results.trainingFits = fits.data();
      samplers[k]->run(0, 3, results);
      range[k][0] = *std::min_element(fits.begin(), fits.end());
      range[k][1] = *std::max_element(fits.begin(), fits.end());
    }
    double width = 2.0 * (range[0][1] - range[0][0]);
    return std::isfinite(range[1][0]) && std::isfinite(range[1][1]) &&
           range[1][0] >= range[0][0] - width &&
           range[1][1] <= range[0][1] + width;
  };
  auto replace = [&](SamplerBase& sampler, const std::vector<double>& x) {
    double anchorMin, anchorMax, newMin, newMax;
    sampler.getAnchor(anchorMin, anchorMax);
    std::vector<double> grid(sampler.data().cutPoints[0]);
    double constant;
    DataReplacement replacement;
    if (constantOf(x, &constant)) {
      replacement.standardization = DerivedValueChoice::supply;
      replacement.suppliedStandardization =
        SuppliedStandardization{&constant, &noScale};
    }
    bool installed = sampler.setData(x.data(), y.data(), x.size() / p, nullptr,
                                     nullptr, nullptr, 0, nullptr, replacement);
    sampler.getAnchor(newMin, newMax);
    check(installed && newMin == anchorMin && newMax == anchorMax &&
            sampler.data().cutPoints[0] == grid,
          "no-spread conversion: the replacement installs, the routing grid "
          "and the response transform standing");
  };
  const std::vector<std::pair<size_t, int>> liveMap = {{0, -1}, {1, -1}};
  const double nan = std::numeric_limits<double>::quiet_NaN();
  double center = nan, scale = nan;

  // --- neither side has spread: one constant replaced by another, at values
  // whose mean is exact (1000) and whose mean rounds (0.1)
  for (double constant : {1000.0, 0.1}) {
    const std::vector<double> x = design(n, constant, 0.0),
      xOther = design(n2, 2.0 * constant, 0.0),
      xDonee = design(n, 2.0 * constant, 0.0),
      xFresh = design(n, constant, 1.0);
    auto build = [&] {
      std::unique_ptr<SamplerBase> sampler = make(x, 10);
      sampler->run(60, numKept, none);
      return sampler;
    };
    std::unique_ptr<SamplerBase> sampler = build(), twin = build(),
      donor = build();
    standardization(*sampler, center, scale);
    check(center == constant && std::isnan(scale),
          "no-spread conversion: a constant column is centered at its value "
          "and stored without a scale");
    check(largestSlope(*sampler) > 1e-3,
          "no-spread conversion: the uninformed slopes are drawn, not zero");
    SamplerStateData before, after, donorState;
    sampler->getState(before);
    const std::vector<double> liveBefore = liveFits(*sampler),
      keptBefore = keptDraws(*sampler, xFresh);
    replace(*sampler, xOther);
    sampler->getState(after);
    check(worstGap(liveFits(*sampler), n2, liveBefore, n, 0, n) < tolerance,
          "no-spread conversion: between two constants the live fit does not "
          "move");
    bool sameBlocks = true;
    for (size_t c = 0; c < numChains; ++c) {
      const ForestStateData& a = before.chains[c].forests[0];
      const ForestStateData& b = after.chains[c].forests[0];
      sameBlocks = sameBlocks && a.treeParams == b.treeParams &&
        a.trees.size() == b.trees.size();
      for (size_t t = 0; sameBlocks && t < a.trees.size(); ++t)
        for (size_t i = 0; sameBlocks && i < a.trees[t].size(); ++i)
          sameBlocks = a.trees[t][i].value == b.trees[t][i].value;
    }
    check(sameBlocks, "no-spread conversion: between two constants neither "
                      "the intercepts nor the slopes move");
    check(worstGap(keptDraws(*sampler, xFresh), n, keptBefore, n, 0, n) <
            tolerance,
          "no-spread conversion: between two constants the kept draws stay "
          "the functions they were");
    check(drawsLikeTwin(*sampler, *twin),
          "no-spread conversion: between two constants the next draws stay "
          "in the twin's range");

    // the same pair through a warm start, which installs the coefficients
    // as stored: each side reads the column at its own center, so the fit is
    // the intercepts' on both
    donor->getState(donorState);
    std::unique_ptr<SamplerBase> recipient = make(xDonee, 20);
    check(warmStartsClean(*recipient, donorState, liveMap, xDonee.data()),
          "no-spread conversion: a donor at another constant installs as it "
          "stands");
    check(worstGap(liveFits(*recipient), n, liveFits(*donor), n, 0, n) <
            tolerance,
          "no-spread conversion: between two constants the seeded fit is the "
          "donor's");
    check(drawsLikeTwin(*recipient, *donor),
          "no-spread conversion: and the seeded chains draw in the donor's "
          "range");
  }

  // --- a placeholder of zeros given values after creation: the column has
  // spread from then on, its slopes are informed, and they convert by the
  // formula
  {
    const std::vector<double> xZero = design(n, 0.0, 0.0),
      xReal = design(n, 0.0, 1.0), xMore = design(n2, 0.0, 1.0);
    auto build = [&] {
      std::unique_ptr<SamplerBase> sampler = make(xZero, 30);
      check(sampler->setPredictor(xReal.data(), true, false) ==
              PredictorUpdateResult::accepted,
            "no-spread conversion: the placeholder takes values");
      sampler->run(60, numKept, none);
      return sampler;
    };
    std::unique_ptr<SamplerBase> sampler = build(), twin = build(),
      donor = build();
    standardization(*sampler, center, scale);
    check(center == 0.0 && scale == 1.0,
          "no-spread conversion: a placeholder column given values is stored "
          "with the 1 it is divided by");
    const std::vector<double> liveBefore = liveFits(*sampler),
      functionBefore = newestDraw(*sampler, xMore);
    check(worstGap(functionBefore, n2, liveBefore, n, 0, n) < tolerance,
          "no-spread conversion: the newest kept draw is the live function");
    replace(*sampler, xMore);
    standardization(*sampler, center, scale);
    check(center != 0.0 && scale != 1.0 && largestSlope(*sampler) > 0.1,
          "no-spread conversion: the standardization is re-derived and the "
          "informed slopes are kept");
    const std::vector<double> liveAfter = liveFits(*sampler);
    check(worstGap(liveAfter, n2, liveBefore, n, 0, n) < tolerance,
          "no-spread conversion: a placeholder given values keeps its live "
          "fit on the old rows");
    check(worstGap(liveAfter, n2, functionBefore, n2, n, n2 - n) < tolerance,
          "no-spread conversion: and is the same function on the appended "
          "rows");
    check(drawsLikeTwin(*sampler, *twin),
          "no-spread conversion: a placeholder given values draws on in the "
          "twin's range");

    // as a donor into a sampler made on the values, whose standardization
    // is other: the slopes go in as stored, under the recipient's own
    SamplerStateData donorState, before, seeded;
    donor->getState(donorState);
    std::unique_ptr<SamplerBase> recipient = make(xReal, 40);
    recipient->getState(before);
    check(recipient->installForests(donorState, liveMap, xReal.data()) ==
            WarmStartResult::ok,
          "no-spread conversion: a donor whose placeholder took values "
          "installs");
    recipient->getState(seeded);
    bool asStored = largestSlope(*recipient) > 0.1;
    for (size_t c = 0; c < numChains; ++c) {
      const ForestStateData& got = seeded.chains[c].forests[0];
      asStored = asStored &&
        got.treeParams == donorState.chains[c].forests[0].treeParams &&
        got.leafCovariateCenters ==
          before.chains[c].forests[0].leafCovariateCenters &&
        got.leafCovariateScales ==
          before.chains[c].forests[0].leafCovariateScales;
    }
    check(asStored && worstGap(liveFits(*recipient), n, liveFits(*donor), n,
                               0, n) > 1e-3,
          "no-spread conversion: its slopes go in as stored under the "
          "recipient's standardization, so the seeded fit is not the "
          "donor's");
  }

  // --- made on one constant and moved to another: every row reads the same
  // value off the center, so the column has spread and the function is kept
  // where the next data spreads the covariate about the new constant
  {
    const std::vector<double> x = design(n, 1000.0, 0.0),
      xMoved = design(n, 2000.0, 0.0), xSpread = design(n2, 2000.0, 1.0);
    std::unique_ptr<SamplerBase> sampler = make(x, 50);
    check(sampler->setPredictor(xMoved.data(), true, false) ==
            PredictorUpdateResult::accepted,
          "no-spread conversion: the constant moves");
    sampler->run(60, numKept, none);
    standardization(*sampler, center, scale);
    check(center == 1000.0 && scale == 1.0,
          "no-spread conversion: a constant column moved off its center is "
          "stored with the 1 it is divided by");
    const std::vector<double> functionBefore = newestDraw(*sampler, xSpread);
    replace(*sampler, xSpread);
    check(worstGap(liveFits(*sampler), n2, functionBefore, n2, 0, n2) <
            tolerance,
          "no-spread conversion: a constant moved off its center keeps its "
          "function on every row of the next data");
  }

  // --- values made constant after creation. Off the center the column goes
  // on reading a value that is not zero under its real scale; at the center
  // exactly it reads zero under a real scale, which is not the placeholder.
  // Either way it has spread until a replacement re-derives the pair
  {
    const std::vector<double> xReal = design(n, 0.0, 1.0),
      xHeld = design(n, 0.25, 0.0), xHeldMore = design(n2, 0.25, 0.0),
      xMore = design(n2, 0.0, 1.0);
    std::unique_ptr<SamplerBase> sampler = make(xReal, 60);
    double realCenter, realScale;
    standardization(*sampler, realCenter, realScale);
    check(sampler->setPredictor(xHeld.data(), true, false) ==
            PredictorUpdateResult::accepted,
          "no-spread conversion: the values are made constant");
    sampler->run(60, numKept, none);
    standardization(*sampler, center, scale);
    check(center == realCenter && scale == realScale && scale != 1.0,
          "no-spread conversion: values made constant keep the center and "
          "scale they were given");
    std::vector<double> functionBefore = newestDraw(*sampler, xHeldMore);
    replace(*sampler, xHeldMore);
    standardization(*sampler, center, scale);
    check(center == 0.25 && std::isnan(scale) && largestSlope(*sampler) > 0.0,
          "no-spread conversion: onto a constant column the slopes are "
          "carried by the formula");
    check(worstGap(liveFits(*sampler), n2, functionBefore, n2, 0, n2) <
            tolerance,
          "no-spread conversion: and each leaf takes its function's value at "
          "the constant");

    const std::vector<double> xAtCenter = design(n, realCenter, 0.0);
    std::unique_ptr<SamplerBase> centered = make(xReal, 70);
    check(centered->setPredictor(xAtCenter.data(), true, false) ==
            PredictorUpdateResult::accepted,
          "no-spread conversion: the values are made the center");
    centered->run(60, numKept, none);
    standardization(*centered, center, scale);
    check(center == realCenter && scale == realScale,
          "no-spread conversion: a column held at its center under a real "
          "scale keeps that scale in the state");
    functionBefore = newestDraw(*centered, xMore);
    replace(*centered, xMore);
    check(worstGap(liveFits(*centered), n2, functionBefore, n2, 0, n2) <
            tolerance,
          "no-spread conversion: a column held at its center under a real "
          "scale converts by the formula");
  }

  // --- rows appended that give a constant column spread, at a constant whose
  // mean rounds: the slopes convert by the formula over the placeholder 1, so
  // each leaf is the function it was, on the old rows and on the appended ones
  {
    const std::vector<double> x = design(n, 0.1, 0.0), xGrown = grown(0.1),
      xFresh = design(n, 0.1, 1.0);
    auto build = [&] {
      std::unique_ptr<SamplerBase> sampler = make(x, 80);
      sampler->run(60, numKept, none);
      return sampler;
    };
    std::unique_ptr<SamplerBase> sampler = build(), twin = build();
    const std::vector<double> liveBefore = liveFits(*sampler),
      functionBefore = newestDraw(*sampler, xGrown),
      keptBefore = keptDraws(*sampler, xFresh);
    replace(*sampler, xGrown);
    standardization(*sampler, center, scale);
    check(std::isfinite(scale) && scale > 0.1 && largestSlope(*sampler) > 0.0,
          "no-spread conversion: off a constant column the slopes are "
          "carried by the formula");
    const std::vector<double> liveAfter = liveFits(*sampler);
    check(worstGap(liveAfter, n2, liveBefore, n, 0, n) < tolerance &&
            worstGap(liveAfter, n2, functionBefore, n2, n, n2 - n) <
              tolerance,
          "no-spread conversion: off a constant column each leaf is the "
          "function it was, on old rows and appended ones");
    check(worstGap(keptDraws(*sampler, xFresh), n, keptBefore, n, 0, n) <
            tolerance,
          "no-spread conversion: off a constant column the kept draws stay "
          "the functions they were");
    check(drawsLikeTwin(*sampler, *twin),
          "no-spread conversion: off a constant column the next draws stay "
          "in the twin's range");
  }

  // --- a store whose kept draws neither fill it nor start at its first
  // slot. No run leaves one; a state whose cursor is not its draw count
  // installs one, and a rewrite of the kept draws must find them in the
  // slots a replay reads
  {
    const std::vector<double> xReal = design(n, 0.0, 1.0),
      xStretched = design(n, 5.0, 3.0), xFresh = design(n, 0.7, 2.0);
    std::unique_ptr<SamplerBase> sampler = make(xReal, 90),
      shifted = make(xReal, 91);
    sampler->run(40, numKept - 1, none);
    SamplerStateData state;
    sampler->getState(state);
    for (ChainStateData& chain : state.chains) {
      ForestStateData& fs = chain.forests[0];
      std::rotate(fs.savedTrees.begin(), fs.savedTrees.end() - numTrees,
                  fs.savedTrees.end());
      std::rotate(fs.savedTreeParams.begin(),
                  fs.savedTreeParams.end() - numTrees,
                  fs.savedTreeParams.end());
    }
    state.currentSampleNum = 0;
    check(shifted->setState(state, nullptr) &&
            shifted->shape().numSavedDraws == numKept - 1 &&
            shifted->savedSlotForDraw(0) == 1,
          "no-spread conversion: a state installs kept draws from the "
          "second slot on");
    const std::vector<double> keptBefore = keptDraws(*shifted, xFresh);
    check(keptBefore == keptDraws(*sampler, xFresh),
          "no-spread conversion: and they replay as the draws they are");
    check(shifted->setData(xStretched.data(), y.data(), n, nullptr, nullptr,
                           nullptr, 0, nullptr),
          "no-spread conversion: the stretched covariate installs");
    check(worstGap(keptDraws(*shifted, xFresh), n, keptBefore, n, 0, n) <
            tolerance,
          "no-spread conversion: kept draws that start past the first slot "
          "are the ones rewritten");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: live conversion across a leaf covariate without spread\n");
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
  check(dest->installForests(donorState, liveMap, x.data()) ==
          WarmStartResult::ok,
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
  check(coarse->installForests(donorState, liveMap, x.data()) ==
          WarmStartResult::ok,
        "variance warm start: a cross-grid donor installs by remapping");
  const double* remapped = TestPeer::varianceFits(coarse->chain(0));
  bool positive = true;
  for (size_t i = 0; i < n; ++i)
    positive &= std::isfinite(remapped[i]) && remapped[i] > 0.0;
  check(positive, "variance warm start: every remapped factor stays positive");
  SamplerStateData remappedState;
  coarse->getState(remappedState);
  check(coarse->setState(remappedState, nullptr),
        "variance warm start: the remapped surface serializes legally");

  // (2b) a variance threshold off the grid alone: the destination's own
  // state as the donor, recorded on its points and one more on column 1,
  // where one variance tree splits. The mean trees stand, the donor is
  // declined, and forced that one tree is the only one moved
  {
    auto own = makeSampler(5060, 100, false);
    own->run(20, 0, empty);
    SamplerStateData offGrid;
    own->getState(offGrid);
    std::vector<double>& points(offGrid.cutPoints[1]);
    double between = 0.5 * (points[40] + points[41]);
    points.insert(points.begin() + 41, between);
    std::vector<FlatNode>& tree(offGrid.chains[0].varianceTrees[0]);
    tree.assign(3, FlatNode());
    tree[0].variable = 1;
    tree[0].value = between;
    setFlatKind(tree[0], FlatKind::ordinal);
    tree[1].value = 1.2;
    tree[2].value = 0.8;
    WarmStartReport report;
    check(warmStartsAltered(*own, offGrid, liveMap, x.data(), nullptr,
                            &report) &&
            !report.mean[0].changed() && report.variance[0].moved == 1 &&
            report.variance[0].cutBack == 0 &&
            report.variance[0].taken == numVarianceTrees,
          "variance warm start: a variance threshold off the grid alone "
          "declines the donor and, forced, is the one tree moved");
  }

  // (3) a variance bottom no row reaches merges at install: x0 <= cut[2] then
  // x0 > cut[8] is a region no row can reach
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
  check(strandTarget->installForests(emptyBottom, liveMap, x.data()) ==
          WarmStartResult::ok,
        "variance warm start: an unoccupied variance bottom installs");
  SamplerStateData strandedAfter;
  strandTarget->getState(strandedAfter);
  check(strandedAfter.chains[0].varianceTrees[0].size() == 3 &&
          strandTarget->chain(0).varianceTree(0).bottomNodesAreOccupied(),
        "variance warm start: and merges into its sibling");

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
  {
    // the variance split is the only change: 1.2 over the rows at or below
    // the cut and 0.8 over the rest merge at their geometric mean by rows
    SamplerStateData meanOnly = outOfMask, merged;
    for (std::vector<FlatNode>& tree : meanOnly.chains[0].varianceTrees) {
      tree.assign(1, FlatNode());
      tree[0].value = 1.1;
    }
    meanOnly.chains[0].varianceTrees[0] = onOne;
    size_t numLeft = 0;
    for (size_t i = 0; i < n; ++i) numLeft += x[i + n] <= onOne[0].value;
    double expected = std::exp((numLeft * std::log(1.2) +
                                (n - numLeft) * std::log(0.8)) / n);
    bool cutBack = warmStartsAltered(*restricted, meanOnly, liveMap, x.data());
    restricted->getState(merged);
    const std::vector<FlatNode>& leaf(merged.chains[0].varianceTrees[0]);
    check(cutBack && numLeft > 0 && numLeft < n && leaf.size() == 1 &&
            std::fabs(leaf[0].value - expected) < 1e-12,
          "variance warm start: an out-of-mask variance tree is declined "
          "and, forced, cut back to one leaf at the geometric mean by rows");
  }
  SamplerStateData compliant = outOfMask;
  for (std::vector<FlatNode>& tree : compliant.chains[0].varianceTrees) {
    tree.assign(1, FlatNode());
    tree[0].value = 1.1;
  }
  auto restricted2 = makeSampler(1344, 100, true);
  check(restricted2->installForests(compliant, liveMap, x.data()) ==
          WarmStartResult::ok,
        "variance warm start: an in-mask variance donor installs");

  // (5) setState holds a state's variance trees to the same mask: without
  // force the state is declined before any chain is touched, so the
  // destination keeps the surface it had
  restricted->run(30, 0, empty);
  check(declinesUntouched(*restricted, outOfMask),
        "variance setState: an out-of-mask variance tree is declined without "
        "force, the surface untouched");
  check(installsClean(*restricted, compliant),
        "variance setState: an in-mask variance state restores");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: variance-forest warm start\n");
}

// The per-forest leaf scale is model: no install carries it. BCF derives both
// forests' scales from the response's SHAPE, so a destination built on a
// differently shaped response constructs different ones, and both a restore
// and a warm start leave them as constructed. y2 keeps y1's range, so the two
// samplers share their units and nothing is converted.
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
  auto scales = [](const Sampler<ConstantGaussianLeaf>& sampler) {
    return std::vector<double>{sampler.forestCalibration(0, 0).priorScale,
                               sampler.forestCalibration(0, 1).priorScale};
  };

  auto donor = make(y1.data());
  Results empty;
  donor->run(10, 2, empty);
  SamplerStateData donorState;
  donor->getState(donorState);

  auto dest = make(y2.data());
  std::vector<double> constructed = scales(*dest);
  // not a vacuous arm: the destination constructed its own, different scales
  check(constructed[0] != scales(*donor)[0] &&
          constructed[1] != scales(*donor)[1],
        "leaf scale: a different-shape response constructs different scales");
  check(dest->setState(donorState, nullptr), "leaf scale: the donor restores");
  check(scales(*dest) == constructed,
        "leaf scale: a restore leaves the destination's scales");

  auto warm = make(y2.data());
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(warm->installForests(donorState, liveMap, x.data()) ==
          WarmStartResult::ok,
        "leaf scale: the donor warm-starts");
  check(scales(*warm) == constructed,
        "leaf scale: a warm start leaves the destination's scales");

  printf("ok: per-forest leaf scale stays the sampler's\n");
}

// The variance forest's own prior draw (docs/design/aft-status-setter.md
// slice 3). The two forest prior-draw entries are mean-only by contract, so
// before this entry a heteroscedastic chain had no path to a prior-drawn
// s(x) at all. Three things are gated here: that the draw leaves LIVE state -
// the surface is the product of the drawn factors, every bottom is occupied
// and the state restores - that it is FOREST-LOCAL, leaving the mean forest,
// sigma and the leaf calibration where it found them, and that both halves
// really are the priors: the structure varies draw to draw and the leaf factor
// matches the calibrated inverse-chi-squared in a moment whose standard error
// is known exactly.
static void testVarianceForestPriorDraw() {
  std::uint64_t savedRngState = rngState;
  rngState = 818181u;

  const size_t n = 300, p = 2, numTrees = 20, numVarianceTrees = 5;
  const double sigmaDf = 3.0, sigmaRawScale = 0.37804942330213542;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = 3.0 * x[i] + (x[i + n] < 0.5 ? 0.2 : 1.4) * z;
  }

  std::vector<ext_rng*> rngs;
  auto makeSampler = [&](std::uint32_t seed, size_t m, double varianceBase) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numVarianceTrees = m;
    options.varianceBase = varianceBase;
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      sigmaDf, sigmaRawScale, options, &r);
  };

  Results empty;
  auto sampler = makeSampler(3131, numVarianceTrees, 0.95);
  sampler->run(60, 0, empty);
  SamplerStateData before;
  sampler->getState(before);
  double sigmaBefore = sampler->chain(0).sigma();
  std::vector<double> surfaceBefore(
    TestPeer::varianceFits(sampler->chain(0)),
    TestPeer::varianceFits(sampler->chain(0)) + n);

  sampler->sampleVarianceForestFromPrior();

  const auto& chain = sampler->chain(0);
  const double* factors = TestPeer::varianceFactors(chain);
  const double* surface = TestPeer::varianceFits(chain);
  bool surfaceMoved = false, productHolds = true, positive = true;
  for (size_t i = 0; i < n; ++i) {
    double product = 1.0;
    for (size_t j = 0; j < numVarianceTrees; ++j) product *= factors[j * n + i];
    if (std::fabs(product - surface[i]) > 1.0e-12 * std::fabs(product))
      productHolds = false;
    if (!(surface[i] > 0.0)) positive = false;
    if (surface[i] != surfaceBefore[i]) surfaceMoved = true;
  }
  check(surfaceMoved, "variance prior draw: the surface moves");
  check(productHolds,
        "variance prior draw: s^2(x) is the product of the drawn factors");
  check(positive, "variance prior draw: every drawn scale is positive");
  for (size_t j = 0; j < numVarianceTrees; ++j)
    check(chain.varianceTree(j).bottomNodesAreOccupied(),
          "variance prior draw: the empty-leaf veto holds on every drawn tree");

  SamplerStateData after;
  sampler->getState(after);
  check(sameFlatTrees(after.chains[0].forests[0].trees,
                      before.chains[0].forests[0].trees),
        "variance prior draw: the mean forest is untouched");
  check(sampler->chain(0).sigma() == sigmaBefore,
        "variance prior draw: sigma is untouched");
  auto restored = makeSampler(4141, numVarianceTrees, 0.95);
  check(restored->setState(after, nullptr),
        "variance prior draw: the drawn state is live state and restores");

  // the tree half really runs: repeated draws do not agree on a node count,
  // and at least one carries a split. A leaf-only entry would report the same
  // count every time.
  size_t firstCount = chain.varianceTree(0).nodes.size();
  bool countVaries = false, anySplit = false;
  for (int rep = 0; rep < 40; ++rep) {
    sampler->sampleVarianceForestFromPrior();
    for (size_t j = 0; j < numVarianceTrees; ++j) {
      size_t nodes = chain.varianceTree(j).nodes.size();
      if (nodes != firstCount) countVaries = true;
      if (nodes > 1) anySplit = true;
    }
  }
  check(anySplit, "variance prior draw: the structure draw grows trees");
  check(countVaries, "variance prior draw: the structure varies by draw");

  // the leaf half is ConstantVarianceLeaf's own prior. One tree and a tree
  // prior that never grows leaves a bare root, so each draw IS one leaf factor
  // h ~ chi^-2(nu, lambda^2) at the m = 1 calibration (nu = sigmaDf, lambda^2 =
  // initialVariance * rawScale). Score the reciprocal, whose distribution is
  // exactly chisq(nu) / (nu lambda^2): mean 1 / lambda^2 and variance 2 / (nu
  // lambda^4), so the Monte Carlo error of the mean is known in closed form
  // and the bound below is five standard errors rather than a guess.
  auto flat = makeSampler(5151, 1, 0.0);
  // the seeded surface before any draw IS the initial variance the
  // calibration is stated against, on the WORKING scale the chain holds it in,
  // so the target below assumes nothing about the response transform
  const double initialVariance = TestPeer::varianceFits(flat->chain(0))[0];
  const double leafScale = initialVariance * sigmaRawScale;
  const int numDraws = 4000;
  double sum = 0.0, sumSquares = 0.0;
  for (int rep = 0; rep < numDraws; ++rep) {
    flat->sampleVarianceForestFromPrior();
    check(flat->chain(0).varianceTree(0).nodes.size() == 1,
          "variance prior draw: a zero-growth prior leaves a bare root");
    // u = lambda^2 / h is chisq(nu) / nu: mean 1, variance 2 / nu, both free
    // of the scale, so the two moments below pin lambda^2 and nu separately
    double u = leafScale / TestPeer::varianceFits(flat->chain(0))[0];
    sum += u;
    sumSquares += u * u;
  }
  double mean = sum / numDraws;
  double variance = (sumSquares - numDraws * mean * mean) / (numDraws - 1);
  double standardError = std::sqrt(2.0 / (sigmaDf * numDraws));
  check(std::fabs(mean - 1.0) < 5.0 * standardError,
        "variance prior draw: the leaf factor's scale is the calibrated one");
  check(std::fabs(2.0 / variance - sigmaDf) < 0.6,
        "variance prior draw: the leaf factor's degrees of freedom are the "
        "calibrated ones");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: variance-forest prior draw (%d leaf draws, mean %.4f, "
         "df %.3f)\n", numDraws, mean, 2.0 / variance);
}

// The CURRENT-state read of the variance surface: the same quantity a recorded
// sweep's variance and varianceTest channels carry, on the original scale and
// with no run. The test arm rebuilds before it reports, which is what makes the
// read correct at a state no recorded sweep produced - here a test-predictor
// swap, whose new rows the stored test product knows nothing about.
static void testCurrentVarianceRead() {
  std::uint64_t savedRngState = rngState;
  rngState = 838383u;

  const size_t n = 300, p = 2, numTrees = 20, numVarianceTrees = 5;
  const size_t numSamples = 4, nTest = 6, testBase = 11, testBase2 = 97;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = 3.0 * x[i] + (x[i + n] < 0.5 ? 0.2 : 1.4) * z;
  }
  // the test rows ARE training rows, so the two reads must agree entry for
  // entry: one surface, addressed two ways
  auto testRowsFrom = [&](size_t base) {
    std::vector<double> block(nTest * p);
    for (size_t i = 0; i < nTest; ++i) {
      block[i] = x[base + i];
      block[i + nTest] = x[base + i + n];
    }
    return block;
  };
  std::vector<double> xTest = testRowsFrom(testBase);
  std::vector<double> xTest2 = testRowsFrom(testBase2);

  std::vector<ext_rng*> rngs;
  auto makeSampler = [&](std::uint32_t seed, size_t numVariance) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numVarianceTrees = numVariance;
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, options, &r);
  };

  auto sampler = makeSampler(3131, numVarianceTrees);
  sampler->setTestPredictors(xTest.data(), nTest);
  std::vector<double> sigma(numSamples), recordedTrain(n * numSamples),
    recordedTest(nTest * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.varianceFits = recordedTrain.data();
  results.varianceTestFits = recordedTest.data();
  sampler->run(40, numSamples, results);

  std::vector<double> train(n), test(nTest);
  check(sampler->currentVarianceFits(0, false, train.data()),
        "current variance: a heteroscedastic sampler answers the train read");
  check(sampler->currentVarianceFits(0, true, test.data()),
        "current variance: and the test read with test rows installed");
  bool trainAgrees = true, testAgrees = true, rowsAgree = true;
  const double* lastTrain = recordedTrain.data() + (numSamples - 1) * n;
  const double* lastTest = recordedTest.data() + (numSamples - 1) * nTest;
  for (size_t i = 0; i < n; ++i)
    if (train[i] != lastTrain[i]) trainAgrees = false;
  for (size_t i = 0; i < nTest; ++i) {
    if (test[i] != lastTest[i]) testAgrees = false;
    if (test[i] != train[testBase + i]) rowsAgree = false;
  }
  check(trainAgrees,
        "current variance: the train read is the recorded channel bitwise");
  check(testAgrees,
        "current variance: the test read is the recorded channel bitwise");
  check(rowsAgree,
        "current variance: a test row that is a training row reads the same");

  // the swap: the stored test product was built for the OLD rows, so a read
  // that did not rebuild would report them
  sampler->setTestPredictors(xTest2.data(), nTest);
  check(sampler->currentVarianceFits(0, true, test.data()),
        "current variance: the test read answers after a predictor swap");
  bool swapped = true;
  for (size_t i = 0; i < nTest; ++i)
    if (test[i] != train[testBase2 + i]) swapped = false;
  check(swapped, "current variance: the test read follows the new test rows");

  // and the two refusals, which are the two states the channels report
  // nothing in
  auto plain = makeSampler(4141, 0);
  check(!plain->currentVarianceFits(0, false, train.data()),
        "current variance: a homoscedastic sampler refuses");
  auto noTest = makeSampler(5151, numVarianceTrees);
  check(!noTest->currentVarianceFits(0, true, test.data()),
        "current variance: a test read with no test rows refuses");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: current variance-surface read\n");
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
  nonPositive.chains[0].savedVarianceTrees[0].back().value =
    std::numeric_limits<double>::infinity();
  check(!dest->setState(nonPositive, nullptr),
        "variance saved state: an infinite saved leaf is refused");
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
    check(dest->installForests(donorState, slotMap, x.data()) ==
            WarmStartResult::ok,
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
  check(target->installForests(shortBuffer, lastMap, x.data()) ==
          WarmStartResult::varianceSlotMismatch,
        "variance slot warm start: a short saved buffer is refused");

  SamplerStateData noBlock = donorState;
  noBlock.chains[0].savedVarianceTrees.clear();
  check(target->installForests(noBlock, lastMap, x.data()) ==
          WarmStartResult::varianceSlotMismatch,
        "variance slot warm start: an absent saved buffer is refused");

  // stride: a live block shorter than the buffer's stride would slice across
  // two sweeps' trees and still hand installVarianceForest the right COUNT
  SamplerStateData spliced = donorState;
  spliced.chains[0].varianceTrees.resize(numVarianceTrees - 2);
  check(target->installForests(spliced, lastMap, x.data()) ==
          WarmStartResult::varianceSlotMismatch,
        "variance slot warm start: a stride the live block contradicts is "
        "refused");

  // the positivity rule the saved buffer is held to on the state side: a zero
  // scale leaf annihilates the product a rebuild forms
  SamplerStateData nonPositive = donorState;
  nonPositive.chains[0].savedVarianceTrees[0].back().value = 0.0;
  std::vector<std::pair<size_t, int>> firstMap = {{0, 0}};
  check(target->installForests(nonPositive, firstMap, x.data()) ==
          WarmStartResult::varianceMismatch,
        "variance slot warm start: a non-positive saved scale leaf is refused");
  SamplerStateData infinite = donorState;
  infinite.chains[0].savedVarianceTrees[0].back().value =
    std::numeric_limits<double>::infinity();
  check(target->installForests(infinite, firstMap, x.data()) ==
          WarmStartResult::varianceMismatch,
        "variance slot warm start: an infinite saved scale leaf is refused");
  // and the refusal is SLOT-specific: the same state installs from any slot
  // whose own trees are intact
  check(target->installForests(nonPositive, lastMap, x.data()) ==
          WarmStartResult::ok,
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
        "reapplyWeights is inert on a gaussian chain");

  // the engine installs a state stored under other weights as it stands; the
  // repair above is the host's, and moves latents, not the chain
  SamplerStateData weightedState;
  weightedA.getState(weightedState);
  check(restoresExactly(weightedB, weightedState),
        "a state stored under other weights installs as stored");

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

/// Whether every live tree of every chain, mean and variance, routes a row to
/// each bottom node.
template <typename L>
static bool liveTreesAreOccupied(Sampler<L>& sampler) {
  for (size_t c = 0; c < sampler.numChains(); ++c) {
    auto& chain = sampler.chain(c);
    for (size_t f = 0; f < chain.numForests(); ++f)
      for (size_t t = 0; t < chain.numTreesInForest(f); ++t)
        if (!chain.treeInForest(f, t).bottomNodesAreOccupied()) return false;
    if (chain.hasVarianceForest())
      for (size_t j = 0; j < chain.numVarianceTrees(); ++j)
        if (!chain.varianceTree(j).bottomNodesAreOccupied()) return false;
  }
  return true;
}

/// The live trees and leaf parameters of two states, every chain and forest.
static bool sameLiveTrees(const SamplerStateData& a,
                          const SamplerStateData& b) {
  if (a.chains.size() != b.chains.size()) return false;
  for (size_t c = 0; c < a.chains.size(); ++c) {
    const ChainStateData& x(a.chains[c]);
    const ChainStateData& y(b.chains[c]);
    if (x.forests.size() != y.forests.size()) return false;
    for (size_t f = 0; f < x.forests.size(); ++f)
      if (!sameFlatTrees(x.forests[f].trees, y.forests[f].trees) ||
          x.forests[f].treeParams != y.forests[f].treeParams)
        return false;
    if (!sameFlatTrees(x.varianceTrees, y.varianceTrees)) return false;
  }
  return true;
}

/// A state stored before a forced setPredictor routes no row to some of its
/// bottom nodes. Restoring it, or warm-starting from it on the same grid,
/// merges those nodes exactly as the forced update merged them, leaves every
/// live tree occupied, and the sampler then runs and restores itself.
template <typename L>
static void checkStaleStateMerges(Sampler<L>& sampler,
                                  const std::vector<double>& xNew,
                                  const char* label) {
  Results empty;
  sampler.run(40, 0, empty);
  SamplerStateData stale;
  sampler.getState(stale);
  sampler.setPredictor(xNew.data(), true, false);
  SamplerStateData forced;
  sampler.getState(forced);
  bool merged = !sameLiveTrees(stale, forced);

  bool restores = restoresAltered(sampler, stale);
  SamplerStateData restored;
  sampler.getState(restored);
  bool restoredAsForced = statesAgree(forced, restored);
  bool restoredOccupied = liveTreesAreOccupied(sampler);

  std::vector<std::pair<size_t, int>> liveMap;
  for (size_t c = 0; c < sampler.numChains(); ++c) liveMap.push_back({c, -1});
  bool warmStarts = sampler.installForests(stale, liveMap, xNew.data()) ==
    WarmStartResult::ok;
  SamplerStateData warm;
  sampler.getState(warm);
  bool warmAsForced = sameLiveTrees(forced, warm);
  bool warmOccupied = liveTreesAreOccupied(sampler);

  sampler.run(5, 0, empty);
  SamplerStateData after;
  sampler.getState(after);
  bool continues = restoresExactly(sampler, after);
  // one rule sent the way its column cannot route: dropped, nothing merged
  SamplerStateData flagged(after);
  bool dropReported = sendFirstOrdinalRuleMissingRight(
                        flagged.chains[0].forests[0].trees) &&
    restoresAltered(sampler, flagged);

  char line[160];
  snprintf(line, sizeof line, "%s: a dropped direction is reported", label);
  check(dropReported, line);
  snprintf(line, sizeof line, "%s: the forced update merged a stale leaf",
           label);
  check(merged, line);
  snprintf(line, sizeof line, "%s: the stale state restores, reported altered",
           label);
  check(restores, line);
  snprintf(line, sizeof line, "%s: the restore merges as the forced update",
           label);
  check(restoredAsForced, line);
  snprintf(line, sizeof line, "%s: no restored tree has an empty leaf", label);
  check(restoredOccupied, line);
  snprintf(line, sizeof line, "%s: a same-grid warm start installs", label);
  check(warmStarts, line);
  snprintf(line, sizeof line, "%s: the warm start merges as the forced update",
           label);
  check(warmAsForced, line);
  snprintf(line, sizeof line, "%s: no warm-started tree has an empty leaf",
           label);
  check(warmOccupied, line);
  snprintf(line, sizeof line,
           "%s: the sampler runs and restores itself as stored", label);
  check(continues, line);
}

static void testStaleStateMerge() {
  // The data come from a pinned generator state, not from wherever the tests
  // before this one left it: a forced update merges a leaf only where a chain
  // nested two splits on x0, which 40 sweeps leave on some draws of the data
  // and not on others. Each kind's "merged a stale leaf" check asserts that
  // premise. The state is left advanced, not restored, so the tests after
  // this one read the same stream whichever suites ran first.
  rngState = 303243635367295568ull;
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), y(n), z(n), xNew;
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    y[i] = 4.0 * (x[i] > 0.5 ? 1.0 : 0.0) + x[i + n] +
      0.2 * (runif01() - 0.5);
    z[i] = i % 2 == 0 ? 1.0 : 0.0;
  }
  // the extremes stay, so the grid would not move; the interior collapses
  xNew = x;
  double lo = *std::min_element(x.begin(), x.begin() + n);
  double hi = *std::max_element(x.begin(), x.begin() + n);
  for (size_t i = 0; i < n; ++i)
    if (xNew[i] > lo && xNew[i] < hi) xNew[i] = 0.5;

  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  };
  const double rawScale = 0.37804942330213542;

  {
    SamplerOptions options;
    options.numTrees = 10;
    options.numVarianceTrees = 4;
    ext_rng* r = newRng(711u);
    ConstantLeafSampler s(x.data(), y.data(), n, p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0, rawScale,
                          options, &r);
    checkStaleStateMerges(s, xNew, "constant and variance leaves");
  }
  std::vector<size_t> covariates = {0};
  {
    SamplerOptions options;
    options.numTrees = 10;
    options.leafCovariateColumns = covariates.data();
    options.numLeafCovariates = 1;
    ext_rng* r = newRng(712u);
    Sampler<LinearGaussianLeaf> s(x.data(), y.data(), n, p, nullptr, nullptr,
                                  ResponseFamily::gaussian, 1.0, 3.0, rawScale,
                                  options, &r);
    checkStaleStateMerges(s, xNew, "linear leaf");
  }
  {
    SamplerOptions options;
    options.numTrees = 10;
    options.gpLeaves = true;
    options.leafCovariateColumns = covariates.data();
    options.numLeafCovariates = 1;
    ext_rng* r = newRng(713u);
    Sampler<GPGaussianLeaf> s(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0, rawScale,
                              options, &r);
    checkStaleStateMerges(s, xNew, "gp leaf");
  }
  {
    std::int8_t directions[] = {1, 0};
    SamplerOptions options;
    options.numTrees = 10;
    options.monotoneDirections = directions;
    ext_rng* r = newRng(714u);
    Sampler<MonotoneConstantGaussianLeaf> s(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, rawScale, options, &r);
    checkStaleStateMerges(s, xNew, "monotone leaf");
  }
  {
    SamplerOptions options;
    AmplitudeSpec spec;
    spec.mu.numTrees = 10;
    spec.mu.base = 0.95;
    spec.mu.power = 2.0;
    spec.tau.numTrees = 6;
    spec.tau.base = 0.25;
    spec.tau.power = 3.0;
    spec.z = z.data();
    ext_rng* r = newRng(715u);
    ConstantLeafSampler s(x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
                          rawScale, options, spec, &r);
    checkStaleStateMerges(s, xNew, "two-forest sampler");
  }
  for (ext_rng* r : rngs) ext_rng_destroy(r);
  printf("ok: a stale state merges empty leaves on install\n");
}

/// A state or a warm-start donor in which one column's grid holds a point
/// twice is refused by the engine's own check, whatever reader let it by,
/// and the sampler it was offered to draws what its untouched twin draws.
static void testRepeatedGridRefused() {
  std::uint64_t savedRngState = rngState;
  rngState = 60606u;
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = 4.0 * (x[i] > 0.5 ? 1.0 : 0.0) + x[i + n] + 0.2 * (runif01() - 0.5);
  ext_rng* rngs[2];
  std::unique_ptr<ConstantLeafSampler> samplers[2];
  for (size_t k = 0; k < 2; ++k) {
    rngs[k] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(rngs[k], 2718u);
    SamplerOptions options;
    options.numTrees = 10;
    samplers[k] = std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &rngs[k]);
    Results empty;
    samplers[k]->run(30, 0, empty);
  }
  ConstantLeafSampler& sampler = *samplers[0];
  SamplerStateData own, repeated;
  sampler.getState(own);
  repeated = own;
  std::vector<double>& cuts = repeated.cutPoints[1];
  cuts.insert(cuts.begin() + 7, cuts[7]);  // one point given twice
  check(!sampler.setState(repeated, nullptr),
        "a state whose grid repeats a point is refused");
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(sampler.installForests(repeated, liveMap, x.data()) !=
          WarmStartResult::ok,
        "a warm-start donor whose grid repeats a point is refused");
  check(sampler.data().cutPoints == own.cutPoints,
        "the refusals leave the sampler's grids as they were");

  std::vector<double> draws[2];
  for (size_t k = 0; k < 2; ++k) {
    draws[k].resize(5);
    Results results;
    results.sigma = draws[k].data();
    samplers[k]->run(0, 5, results);
  }
  check(draws[0] == draws[1],
        "after the refusals the sampler draws what its twin draws");
  check(restoresExactly(sampler, own) &&
          sampler.installForests(own, liveMap, x.data()) == WarmStartResult::ok,
        "the state and the donor with distinct grids are taken");
  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a grid that repeats a point is refused\n");
}

/// A weighted grid's weights ride the state and every grid rollback: a state
/// installs the weights it carries and one carrying none leaves the grids
/// unweighted; a refused state, a refused refresh and a warm start's scoped
/// donor grid each leave them as they were; setCutPoints keeps them for the
/// grid the column holds and drops them for another. A continuation after
/// the round trip draws what the twin draws.
static void testWeightedGridState() {
  std::uint64_t savedRngState = rngState;
  rngState = 61616u;
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>((i % 6) * (i % 6));  // widths 1, 3, ..., 9
    x[i + n] = runif01();
    y[i] = 0.3 * x[i] + x[i + n] + 0.2 * (runif01() - 0.5);
  }
  ext_rng* rngs[2];
  std::unique_ptr<ConstantLeafSampler> samplers[2];
  for (size_t k = 0; k < 2; ++k) {
    rngs[k] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(rngs[k], 3141u);
    SamplerOptions options;
    options.numTrees = 10;
    samplers[k] = std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &rngs[k]);
    Results empty;
    samplers[k]->run(30, 0, empty);
  }
  ConstantLeafSampler& sampler = *samplers[0];
  const std::vector<std::vector<double>> weights = sampler.data().cutMass;
  SamplerStateData own;
  sampler.getState(own);
  check(sampler.data().cutsWeighted(0) && !sampler.data().cutsWeighted(1) &&
          own.cutMass == weights,
        "a state carries the grid's weights");

  // a state without weights goes onto the store's weights, which stay, as
  // its own goes onto them
  SamplerStateData bare(own);
  bare.cutMass.clear();
  bool gridDiffers = false;
  check(sampler.setState(bare, nullptr, nullptr, nullptr, keepStoreCapacity,
                         true, nullptr, &gridDiffers) &&
          gridDiffers && sampler.data().cutMass == weights,
        "a state carrying no weights goes onto the store's, which stay");
  check(restoresExactly(sampler, own) && sampler.data().cutMass == weights,
        "a state carrying the store's weights installs as stored");

  // refusals leave the weights: an invalid weight, a state refused after its
  // grid went in, a refused refresh, and a cross-grid donor's scoped grid
  SamplerStateData badMass(own), badTrees(bare);
  badMass.cutMass[0][2] = badMass.cutMass[0][1];
  badTrees.chains[0].forests[0].trees.pop_back();
  check(!sampler.setState(badMass, nullptr) &&
          !sampler.setState(badTrees, nullptr) &&
          sampler.data().cutMass == weights,
        "a refused state leaves the weights as they were");
  std::vector<double> constant(n, 0.0);
  const size_t column = 0;
  check(sampler.updatePredictor(constant.data(), &column, 1, false, true) ==
            PredictorUpdateResult::rolledBack &&
          sampler.data().cutMass == weights,
        "a refresh rolled back puts the weights back with the grid");
  std::vector<double> whole(x);
  std::fill(whole.begin(), whole.begin() + n, 0.0);
  check(sampler.setPredictor(whole.data(), false, true) ==
            PredictorUpdateResult::rolledBack &&
          sampler.data().cutMass == weights,
        "and so does a whole-matrix refresh");
  SamplerStateData donor(own);
  donor.cutPoints[1].push_back(donor.cutPoints[1].back() + 1.0);
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  WarmStartResult warm = sampler.installForests(donor, liveMap, x.data());
  check(warm == WarmStartResult::ok &&
          sampler.data().cutMass == weights,
        "a cross-grid warm start leaves the live weights");
  check(restoresExactly(sampler, own) && restoresExactly(*samplers[1], own),
        "the state goes back in, and into the twin");

  // the continuation after the round trips draws what the twin, restored
  // from the state alone, draws
  std::vector<double> draws[2];
  for (size_t k = 0; k < 2; ++k) {
    draws[k].resize(5);
    Results results;
    results.sigma = draws[k].data();
    samplers[k]->run(0, 5, results);
  }
  // to rounding: a rebuild routes rows in an order its history sets
  double gap = 0.0;
  for (size_t d = 0; d < draws[0].size(); ++d)
    gap = std::max(gap, std::fabs(draws[0][d] - draws[1][d]));
  check(gap < 1e-10,
        "a sampler whose weights went out and back draws what its twin draws");

  // setCutPoints: the held grid keeps its weights, another drops them
  std::vector<double> held(sampler.data().cutPoints[0]);
  const double* grids[] = {held.data()};
  std::uint32_t counts[] = {static_cast<std::uint32_t>(held.size())};
  sampler.setCutPoints(grids, counts, &column, 1, x.data());
  check(sampler.data().cutMass == weights,
        "setCutPoints handed the held grid keeps its weights");
  held.back() += 1.0;
  sampler.setCutPoints(grids, counts, &column, 1, x.data());
  check(!sampler.data().cutsWeighted(0),
        "setCutPoints handed another grid leaves it unweighted");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a grid's weights ride states and rollbacks\n");
}

/// Whether any live tree of any chain, mean or variance, holds an ordinal
/// split outside the interval its ancestors leave.
template <typename L>
static bool liveTreesHoldSplitOutsideInterval(Sampler<L>& sampler) {
  auto holds = [&](const Tree& tree) {
    std::vector<int32_t> subtree;
    tree.fillSubtree(0, subtree);
    for (int32_t i : subtree)
      if (!tree.at(i).isBottom() &&
          tree.splitIsOutsideInterval(sampler.data(), i))
        return true;
    return false;
  };
  for (size_t c = 0; c < sampler.numChains(); ++c) {
    auto& chain = sampler.chain(c);
    for (size_t f = 0; f < chain.numForests(); ++f)
      for (size_t t = 0; t < chain.numTreesInForest(f); ++t)
        if (holds(chain.treeInForest(f, t))) return true;
    if (chain.hasVarianceForest())
      for (size_t j = 0; j < chain.numVarianceTrees(); ++j)
        if (holds(chain.varianceTree(j))) return true;
  }
  return false;
}

/// A state no sampler writes: a tree whose root and left child split one
/// column at one value, beside missing values that keep both sides of the
/// child occupied. It installs with the child merged and the install
/// reported as altered, in the mean forest and in the variance forest each
/// alone; the same nodes on distinct values install as stored.
static void testStackedSplitsMerge() {
  std::uint64_t savedRngState = rngState;
  rngState = 90210u;
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = 4.0 * (x[i] > 0.5 ? 1.0 : 0.0) + x[i + n] + 0.2 * (runif01() - 0.5);
  for (size_t i = 0; i < n; i += 4) x[i] = std::nan("");

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 4242u);
  SamplerOptions options;
  options.numTrees = 10;
  options.numVarianceTrees = 4;
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(30, 0, empty);
  SamplerStateData own;
  sampler.getState(own);
  check(restoresExactly(sampler, own) &&
          !liveTreesHoldSplitOutsideInterval(sampler),
        "stacked splits: the sampler's own state installs as stored");

  const std::vector<double>& cuts = own.cutPoints[0];
  size_t middle = cuts.size() / 2;
  auto split = [&](double value, bool missingGoesRight) {
    FlatNode node;
    node.variable = 0;
    setFlatKind(node, FlatKind::ordinal);
    node.value = value;
    if (missingGoesRight) node.flags |= flatMissingGoesRight;
    return node;
  };
  auto leaf = [](double value) {
    FlatNode node;
    node.value = value;
    return node;
  };
  auto treeOn = [&](double childValue, double scale) {
    return std::vector<FlatNode>{split(cuts[middle], false),
                                 split(childValue, true), leaf(0.5 * scale),
                                 leaf(1.5 * scale), leaf(2.0 * scale)};
  };

  SamplerStateData ordered(own);
  ordered.chains[0].forests[0].trees[0] = treeOn(cuts[middle / 2], 0.01);
  ordered.chains[0].varianceTrees[0] = treeOn(cuts[middle / 2], 1.0);
  check(restoresExactly(sampler, ordered),
        "stacked splits: distinct values in order install as stored");

  SamplerStateData meanStacked(ordered), varianceStacked(ordered);
  meanStacked.chains[0].forests[0].trees[0] = treeOn(cuts[middle], 0.01);
  varianceStacked.chains[0].varianceTrees[0] = treeOn(cuts[middle], 1.0);
  SamplerStateData restored;
  bool meanAltered = restoresAltered(sampler, meanStacked);
  sampler.getState(restored);
  check(meanAltered && restored.chains[0].forests[0].trees[0].size() == 3 &&
          restored.chains[0].varianceTrees[0].size() == 5 &&
          !liveTreesHoldSplitOutsideInterval(sampler) &&
          liveTreesAreOccupied(sampler),
        "stacked splits: the mean tree's child is merged, reported altered");
  check(restoresExactly(sampler, restored),
        "stacked splits: the merged state installs as stored");
  bool varianceAltered = restoresAltered(sampler, varianceStacked);
  sampler.getState(restored);
  check(varianceAltered && restored.chains[0].varianceTrees[0].size() == 3 &&
          restored.chains[0].forests[0].trees[0].size() == 5 &&
          !liveTreesHoldSplitOutsideInterval(sampler) &&
          liveTreesAreOccupied(sampler),
        "stacked splits: the variance tree's child is merged, reported altered");

  // a same-grid warm start merges it too
  std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  check(sampler.installForests(meanStacked, liveMap, x.data()) ==
          WarmStartResult::ok &&
          !liveTreesHoldSplitOutsideInterval(sampler) &&
          liveTreesAreOccupied(sampler),
        "stacked splits: a warm start from the state merges the child");

  sampler.run(5, 0, empty);
  SamplerStateData after;
  sampler.getState(after);
  bool finite = true;
  for (const std::vector<FlatNode>& tree : after.chains[0].forests[0].trees)
    for (const FlatNode& node : tree) finite = finite && std::isfinite(node.value);
  check(finite && restoresExactly(sampler, after),
        "stacked splits: the sampler runs on and restores itself as stored");

  ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: a state stacking two splits on one value installs merged\n");
}

/// Sampler::setState's altered flag, one cause at a time on states that
/// differ from an exact one in that cause alone: a merge in the variance
/// trees only and in the mean trees only, a direction dropped from a mean
/// and from a variance tree with nothing merged, and one inexact chain of
/// two. A state refused after its scratch builds dropped a direction reports
/// nothing.
static void testRestoreStatus() {
  std::uint64_t savedRngState = rngState;
  rngState = 31337u;
  const size_t n = 200, p = 2, numChains = 2;
  std::vector<double> x(n * p), y(n), xNew;
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = 4.0 * (x[i] > 0.5 ? 1.0 : 0.0) + x[i + n] +
      (x[i] > 0.5 ? 2.0 : 0.2) * (runif01() - 0.5);
  xNew = x;
  double lo = *std::min_element(x.begin(), x.begin() + n);
  double hi = *std::max_element(x.begin(), x.begin() + n);
  for (size_t i = 0; i < n; ++i)
    if (xNew[i] > lo && xNew[i] < hi) xNew[i] = 0.5;

  std::vector<ext_rng*> rngs(numChains);
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(rngs[c], 731u + static_cast<std::uint_least32_t>(c));
  }
  SamplerOptions options;
  options.numTrees = 10;
  options.numVarianceTrees = 6;
  options.numChains = numChains;
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, rngs.data());
  Results empty;
  sampler.run(60, 0, empty);
  SamplerStateData own, forced, restored;
  sampler.getState(own);
  check(restoresExactly(sampler, own),
        "restore status: a two-chain sampler's own state installs as stored");

  // one rule of one tree of the second chain: no row moves, nothing merges
  SamplerStateData meanFlagged(own), varianceFlagged(own);
  bool flagged = sendFirstOrdinalRuleMissingRight(
                   meanFlagged.chains[1].forests[0].trees) &&
    sendFirstOrdinalRuleMissingRight(varianceFlagged.chains[1].varianceTrees);
  check(flagged, "restore status: a mean and a variance tree split");
  bool meanDropped = restoresAltered(sampler, meanFlagged);
  sampler.getState(restored);
  check(meanDropped && statesAgree(own, restored) &&
          restoresExactly(sampler, own),
        "restore status: a direction dropped from one chain's mean tree");
  bool varianceDropped = restoresAltered(sampler, varianceFlagged);
  sampler.getState(restored);
  check(varianceDropped && statesAgree(own, restored) &&
          restoresExactly(sampler, own),
        "restore status: a direction dropped from one chain's variance tree");

  // the forced update's state with the stored trees of one kind put back
  sampler.setPredictor(xNew.data(), true, false);
  sampler.getState(forced);
  SamplerStateData meanStale(forced), varianceStale(forced);
  bool meanMerged = false, varianceMerged = false;
  for (size_t c = 0; c < numChains; ++c) {
    meanStale.chains[c].forests[0].trees = own.chains[c].forests[0].trees;
    varianceStale.chains[c].varianceTrees = own.chains[c].varianceTrees;
    meanMerged = meanMerged ||
      !sameFlatTrees(own.chains[c].forests[0].trees,
                     forced.chains[c].forests[0].trees);
    varianceMerged = varianceMerged ||
      !sameFlatTrees(own.chains[c].varianceTrees,
                     forced.chains[c].varianceTrees);
  }
  check(meanMerged && varianceMerged,
        "restore status: the forced update merged trees of each kind");
  check(restoresAltered(sampler, varianceStale) &&
          liveTreesAreOccupied(sampler) && restoresExactly(sampler, forced),
        "restore status: a merge in the variance trees alone");
  check(restoresAltered(sampler, meanStale) && liveTreesAreOccupied(sampler) &&
          restoresExactly(sampler, forced),
        "restore status: a merge in the mean trees alone");
  // refused for its last tree, after the first chain's merges and a dropped
  // direction passed validation
  SamplerStateData bad(meanStale);
  FlatNode childless;
  childless.variable = 0;
  setFlatKind(childless, FlatKind::ordinal);
  childless.value = bad.cutPoints[0][0];
  bad.chains[1].forests[0].trees.back().assign(1, childless);
  bool altered = sendFirstOrdinalRuleMissingRight(
    bad.chains[0].forests[0].trees);
  bool notClean = true;
  check(!sampler.setState(bad, nullptr, nullptr, &altered) && !altered,
        "restore status: a refused state reports nothing");
  check(!sampler.setState(bad, nullptr, nullptr, nullptr, keepStoreCapacity,
                          false, &notClean) && !notClean,
        "restore status: a refusal is not a verdict, whatever its other "
        "trees need");

  // declined across cut grids: a state on a coarser grid, drawn over the rows
  // before the update, holds bottoms the collapsed rows do not reach. The
  // grid, the codes and the chains are left; forced, the state's splits are
  // moved onto the store's grid, which stays
  {
    std::vector<ext_rng*> coarseRngs(numChains);
    for (size_t c = 0; c < numChains; ++c) {
      coarseRngs[c] =
        ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
      ext_rng_setSeed(coarseRngs[c], 911u + static_cast<std::uint_least32_t>(c));
    }
    SamplerOptions coarseOptions(options);
    coarseOptions.maxNumCuts = 20;
    ConstantLeafSampler coarse(x.data(), y.data(), n, p, nullptr, nullptr,
                               ResponseFamily::gaussian, 1.0, 3.0,
                               0.37804942330213542, coarseOptions,
                               coarseRngs.data());
    coarse.run(60, 0, empty);
    SamplerStateData crossGrid;
    coarse.getState(crossGrid);
    std::vector<xint_t> codes = storageDigest(sampler.data());
    check(crossGrid.cutPoints != forced.cutPoints &&
            declinesUntouched(sampler, crossGrid, xNew.data()) &&
            sampler.data().cutPoints == forced.cutPoints &&
            storageDigest(sampler.data()) == codes,
          "restore status: a state declined on another grid leaves the grid, "
          "the codes and the chains");
    check(restoresAltered(sampler, crossGrid, xNew.data()) &&
            sampler.data().cutPoints == forced.cutPoints &&
            storageDigest(sampler.data()) == codes &&
            liveTreesAreOccupied(sampler),
          "restore status: forced, the state goes onto the store's grid");
    for (ext_rng* r : coarseRngs) ext_rng_destroy(r);
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a state install reports whether it altered the state\n");
}

/// A flat tree from a sampler whose column could hold a missing value builds
/// against a store whose column never has: the direction is dropped, a rule
/// that sends every level one way is refused, and a mask malformed in the
/// donor's own gauge is still refused. Against a store whose column has held
/// one and is complete again, the rule that splits the missing value from
/// every level builds with its direction.
static void testStaleMissingDirectionBuild() {
  const size_t n = 12;
  std::vector<double> x(n * 2);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 4);
  for (size_t i = 0; i < n; ++i) x[i + n] = static_cast<double>(i) / (n - 1.0);
  ColumnKind types[] = {ColumnKind::categorical, ColumnKind::numeric};
  ColumnStore store;
  built(store.build(x.data(), n, 2, 10, false, types));

  std::vector<index_t> indices(n);
  std::vector<double> params;
  Tree tree;
  bool dropped = false;
  auto buildsOn = [&](const ColumnStore& onto,
                      const std::vector<FlatNode>& flat) {
    tree.initialize(indices.data(), n);
    dropped = false;
    return tree.buildFromFlat(onto, flat.data(), flat.size(), params, 1,
                              nullptr, nullptr, 0, &dropped);
  };
  auto builds = [&](const std::vector<FlatNode>& flat) {
    return buildsOn(store, flat);
  };
  auto rule = [](int32_t variable, FlatKind kind, bool missingRight) {
    FlatNode node;
    node.variable = variable;
    setFlatKind(node, kind);
    if (missingRight) node.flags |= flatMissingGoesRight;
    return node;
  };
  auto levels = [&](std::uint64_t mask, bool missingRight) {
    FlatNode node = rule(0, FlatKind::categoricalInline, missingRight);
    node.mask = mask;
    return node;
  };
  FlatNode leaf, cut = rule(1, FlatKind::ordinal, true);
  cut.value = store.cutPoints[1][4];

  check(builds({levels(0x6, true), cut, leaf, leaf, leaf}) &&
          tree.at(0).rule.categoryDirections() == 0x6 &&
          !tree.at(tree.at(0).leftChild).rule.missingGoesRight() &&
          tree.at(tree.at(0).leftChild).rule.splitIndex() == 4,
        "stale direction: an ordinal and a categorical rule build without it");
  bool bothDropped = dropped;
  FlatNode plain = cut;
  plain.flags &= static_cast<std::uint8_t>(~flatMissingGoesRight);
  check(bothDropped &&
          builds({levels(0x6, false), plain, leaf, leaf, leaf}) && !dropped &&
          builds({levels(0x6, false), cut, leaf, leaf, leaf}) && dropped &&
          builds({levels(0x6, true), plain, leaf, leaf, leaf}) && dropped,
        "stale direction: a build reports a drop of either kind, and none "
        "where the tree carries none");
  check(!builds({levels(0x0, true), leaf, leaf}) &&
          !builds({levels(0xf, false), leaf, leaf}) &&
          !builds({levels(0x6, false), leaf, levels(0x6, false), leaf, leaf}) &&
          !builds({levels(0x6, true), leaf, levels(0x0, true), leaf, leaf}),
        "stale direction: on a column never missing, a rule sending every "
        "level that reaches it one way is refused");
  // the same column once it has held a missing value, complete again
  std::vector<double> xOnceMissing(x);
  xOnceMissing[0] = std::nan("");
  ColumnStore seen;
  built(seen.build(xOnceMissing.data(), n, 2, 10, false, types));
  seen.setPredictors(x.data(), false);
  bool complete = seen.hasMissing[0] != 0 && seen.hasMissing[1] == 0;
  for (size_t i = 0; i < n; ++i)
    complete = complete && seen.codeAt(0, i) == store.codeAt(0, i);
  check(complete && buildsOn(seen, {levels(0x0, true), leaf, leaf}) &&
          tree.at(0).rule.categoryDirections() == Rule::missingDirectionBit &&
          !dropped && buildsOn(seen, {levels(0xf, false), leaf, leaf}) &&
          tree.at(0).rule.categoryDirections() == 0xf && !dropped &&
          buildsOn(seen, {levels(0x6, true), cut, leaf, leaf, leaf}) &&
          tree.at(0).rule.categoryDirections() ==
            (0x6 | Rule::missingDirectionBit) &&
          dropped,
        "stale direction: on a column that has held one, a rule splitting "
        "missing from every level builds, its direction kept");
  check(builds({levels(0x6, true), leaf, levels(0x2, true), leaf, leaf}) &&
          !builds({levels(0x6, false), leaf, levels(0x2, true), leaf, leaf}),
        "stale direction: it builds only where its ancestors pass missing");
  check(!builds({levels(0xf, true), leaf, leaf}) &&
          !builds({levels(0x0, false), leaf, leaf}),
        "stale direction: a rule sending everything one way is refused");
  check(!builds({levels(0x16, true), leaf, leaf}),
        "stale direction: a level past the column's count is refused");

  // a pooled mask keeps the bit in its words, as a data mutation leaves it
  const std::uint32_t K = 70;
  std::vector<double> wide(8 * K);
  for (size_t i = 0; i < wide.size(); ++i) wide[i] = static_cast<double>(i % K);
  ColumnStore pooled;
  built(pooled.build(wide.data(), wide.size(), 1, 10, false, types));
  std::vector<std::uint64_t> words(maskWordsForCount(K), 0), wordsAfter;
  maskSetBit(words.data(), 2);
  std::vector<FlatNode> flat = {rule(0, FlatKind::categoricalPooled, true),
                                leaf, leaf}, flatAfter;
  flat[0].numMaskWords = static_cast<std::uint32_t>(words.size());
  std::vector<index_t> wideIndices(wide.size());
  auto buildsWide = [&]() {
    tree.initialize(wideIndices.data(), wide.size());
    dropped = false;
    return tree.buildFromFlat(pooled, flat.data(), flat.size(), params, 1,
                              nullptr, words.data(), words.size(), &dropped);
  };
  bool wideBuilt = pooled.columnIsPooled(0) && !pooled.hasMissing[0] &&
    buildsWide() && tree.ruleMissingGoesRight(pooled, tree.at(0).rule);
  check(wideBuilt && !dropped,
        "stale direction: a pooled rule builds and keeps its word, unreported");
  // a refused build leaves the tree half-built
  if (wideBuilt)
    tree.flatten(pooled, params.data(), flatAfter, nullptr, 1, nullptr,
                 &wordsAfter);
  check(wideBuilt && flatAfter[0].flags == flat[0].flags && wordsAfter == words,
        "stale direction: and flattens to the state it was built from");
  maskSetBit(words.data(), K + 1);
  check(!buildsWide(),
        "stale direction: a pooled bit past the missing position is refused");
  // every level left with missing right, and every level right
  words.assign(words.size(), 0);
  bool nothingRight = !buildsWide();
  for (std::uint32_t level = 0; level < K; ++level)
    maskSetBit(words.data(), level);
  flat[0].flags &= static_cast<std::uint8_t>(~flatMissingGoesRight);
  check(nothingRight && !buildsWide(),
        "stale direction: on a pooled column never missing, a rule sending "
        "every level one way is refused");
  flat[0].flags |= flatMissingGoesRight;
  // levels 2 and 3 right at the root, level 2 and missing right beneath it
  FlatNode below = flat[0];
  below.maskOffset = words.size();
  flat = {flat[0], leaf, below, leaf, leaf};
  words.assign(2 * words.size(), 0);
  maskSetBit(words.data(), 2);
  maskSetBit(words.data(), 3);
  maskSetBit(words.data() + words.size() / 2, 2);
  bool passedDown = buildsWide();
  flat[0].flags &= static_cast<std::uint8_t>(~flatMissingGoesRight);
  check(passedDown && !buildsWide(),
        "stale direction: a pooled rule builds only where its ancestors pass "
        "missing");
  printf("ok: a stale missing direction builds\n");
}

/// Rules flagged missing-right among flat trees.
static size_t countMissingRight(const std::vector<std::vector<FlatNode>>& trees) {
  size_t count = 0;
  for (const std::vector<FlatNode>& tree : trees)
    for (const FlatNode& node : tree)
      count += (node.flags & flatMissingGoesRight) != 0 ? 1u : 0u;
  return count;
}

/// A state stored while columns held missing values, restored after a forced
/// setPredictor filled them. The columns stay able to hold one, so the forced
/// update keeps every direction and the install reproduces its trees; the
/// saved draws keep theirs; a sampler created over the filled rows with the
/// first's flags, as a copy or a reload is, takes the state to the same
/// trees; one created over them without the flags refuses it or installs it
/// with no inline direction left; a same-grid warm start takes it too; and
/// the sampler then runs and restores itself. make(x, seed, seen) builds the
/// sampler over a predictor matrix, `seen` being the flags its creation
/// adopts or null; `seed` is one whose draws leave every forest carrying a
/// missing-right rule, which the first report checks. The two-forest arm
/// takes its own: a two-forest stream moved when a leaf of only control rows
/// became legal in the treatment forest, and the shared seed's draws then
/// left that forest no missing-right rule.
template <typename Make>
static void checkStaleDirectionRestores(Make make, const std::vector<double>& x,
                                        const std::vector<double>& xFilled,
                                        bool inlineOnly, const char* label,
                                        std::uint32_t seed = 721u) {
  auto sampler = make(x, seed, nullptr);
  Results empty;
  sampler->run(60, 2, empty);
  SamplerStateData stale, forced, restored, other, warm, after, taken;
  sampler->getState(stale);
  const ChainStateData& staleChain(stale.chains[0]);
  bool carries = !sampler->hasVarianceForest() ||
    countMissingRight(staleChain.varianceTrees) > 0;
  for (const ForestStateData& fs : staleChain.forests)
    carries = carries && countMissingRight(fs.trees) > 0 &&
      countMissingRight(fs.savedTrees) > 0;

  std::vector<std::uint8_t> flags(sampler->data().hasMissing);
  bool filled = sampler->setPredictor(xFilled.data(), true, false) ==
    PredictorUpdateResult::accepted;
  bool flagsKept = sampler->data().hasMissing == flags;
  sampler->getState(forced);
  bool restores = restoresAltered(*sampler, stale);
  sampler->getState(restored);
  size_t staleRight = countMissingRight(staleChain.varianceTrees);
  size_t right = countMissingRight(restored.chains[0].varianceTrees);
  bool savedKept = true;
  for (size_t f = 0; f < staleChain.forests.size(); ++f) {
    staleRight += countMissingRight(staleChain.forests[f].trees);
    right += countMissingRight(restored.chains[0].forests[f].trees);
    savedKept = savedKept &&
      sameFlatTrees(staleChain.forests[f].savedTrees,
                    restored.chains[0].forests[f].savedTrees);
  }
  bool occupied = liveTreesAreOccupied(*sampler);

  // created over the filled rows: with the flags, and without them
  auto recipient = make(xFilled, seed + 1u, flags.data());
  bool flagsAdopted = recipient->data().hasMissing == flags;
  bool otherRestores = restoresAltered(*recipient, stale, xFilled.data());
  recipient->getState(other);
  auto stranger = make(xFilled, seed + 2u, nullptr);
  bool neverMissing = true;
  for (std::uint8_t flag : stranger->data().hasMissing)
    neverMissing = neverMissing && flag == 0;
  bool strangerAltered = false;
  bool strangerDeclines = declinesUntouched(*stranger, stale, xFilled.data());
  bool strangerTook = stranger->setState(stale, xFilled.data(), nullptr,
                                         &strangerAltered);
  stranger->getState(taken);
  size_t strangerRight = countMissingRight(taken.chains[0].varianceTrees);
  for (const ForestStateData& fs : taken.chains[0].forests)
    strangerRight += countMissingRight(fs.trees);

  bool warmStarts = sampler->installForests(stale, {{0, -1}}, x.data()) ==
    WarmStartResult::ok;
  sampler->getState(warm);
  sampler->run(5, 0, empty);
  recipient->run(5, 0, empty);
  sampler->getState(after);
  bool continues = restoresExactly(*sampler, after);

  char line[160];
  auto report = [&](bool ok, const char* what) {
    snprintf(line, sizeof line, "stale direction, %s: %s", label, what);
    check(ok, line);
  };
  report(carries, "every forest's live and saved trees send missing right");
  report(filled && flagsKept && restores,
         "the columns stay flagged once filled and the stale state restores, "
         "reported altered");
  report(statesAgree(forced, restored),
         "the restore reproduces the forced update");
  report(right > 0 && right <= staleRight,
         "the restored trees keep their directions");
  report(occupied && savedKept,
         "no leaf is left empty and the saved draws keep their directions");
  report(flagsAdopted && otherRestores && sameLiveTrees(forced, other),
         "a sampler created over the filled rows with the flags restores it "
         "to the same trees");
  report(neverMissing &&
           (!strangerTook ||
            (strangerDeclines && strangerAltered &&
             (inlineOnly ? strangerRight == 0 : strangerRight < right))),
         "a sampler created without them refuses it or, declining it "
         "unforced, drops the inline directions when forced");
  report(warmStarts && sameLiveTrees(forced, warm),
         "a same-grid warm start installs the same trees");
  report(continues, "the sampler runs and restores itself as stored");
}

static void testStaleMissingDirectionRestores() {
  // a private data stream, so the fixture is the same under a suite filter
  std::uint64_t savedRngState = rngState;
  rngState = 2718u;
  const size_t n = 300, p = 3;
  const std::uint32_t K = 70;
  // an ordinal, an inline categorical and a pooled categorical column, each
  // missing in rows whose mean, treatment effect or spread differs
  std::vector<double> x(n * p), y(n), z(n), xFilled;
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = static_cast<double>(i % 4);
    x[i + 2 * n] = static_cast<double>(i % K);
    z[i] = i % 2 == 0 ? 1.0 : 0.0;
  }
  xFilled = x;
  for (size_t i = 0; i < n; ++i) {
    bool gone[] = {i % 5 == 0, i % 7 == 0, i % 3 == 0};
    y[i] = x[i] + (gone[0] ? 2.0 + 2.0 * z[i] : 0.0) +
      (gone[1] ? -2.0 - 2.0 * z[i] : 0.0) + (gone[2] ? 1.5 : 0.0) +
      (gone[0] || gone[1] ? 1.5 : 0.1) * (runif01() - 0.5);
    for (size_t j = 0; j < p; ++j)
      if (gone[j] && i >= K) x[i + j * n] = std::nan("");
  }
  ColumnKind kinds[] = {ColumnKind::numeric, ColumnKind::categorical,
                        ColumnKind::categorical};
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr));
    ext_rng_setSeed(rngs.back(), seed);
    return &rngs.back();
  };
  const double rawScale = 0.37804942330213542;
  SamplerOptions options;
  options.numTrees = 20;
  options.predictors.columnTypes = kinds;
  options.keepTrees = true;
  options.numSamplesToStore = 2;
  auto single = [&](size_t numColumns) {
    return [&, numColumns](const std::vector<double>& data,
                           std::uint32_t seed, const std::uint8_t* seen) {
      SamplerOptions created(options);
      created.missingSeen = seen;
      return std::make_unique<ConstantLeafSampler>(
        data.data(), y.data(), n, numColumns, nullptr, nullptr,
        ResponseFamily::gaussian, 1.0, 3.0, rawScale, created, newRng(seed));
    };
  };
  checkStaleDirectionRestores(single(p), x, xFilled, false, "pooled column");
  options.numVarianceTrees = 10;
  checkStaleDirectionRestores(single(2), x, xFilled, true,
                              "mean and variance trees");
  options.numVarianceTrees = 0;
  AmplitudeSpec spec;
  spec.mu.numTrees = 20;
  spec.mu.base = 0.95;
  spec.mu.power = 2.0;
  spec.tau.numTrees = 12;
  spec.tau.base = 0.25;
  spec.tau.power = 3.0;
  spec.z = z.data();
  checkStaleDirectionRestores(
    [&](const std::vector<double>& data, std::uint32_t seed,
        const std::uint8_t* seen) {
      SamplerOptions created(options);
      created.missingSeen = seen;
      return std::make_unique<ConstantLeafSampler>(
        data.data(), y.data(), n, 2, nullptr, nullptr, 1.0, 3.0, rawScale,
        created, spec, newRng(seed));
    },
    x, xFilled, true, "two-forest sampler", 724u);
  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a filled column's directions are kept on install\n");
}

/// A forced install collapses a tree at the first split on each path from the
/// root that its forest's column list or interaction constraint forbids: the
/// splits above it stay, everything beneath it becomes one leaf at the
/// row-weighted mean of the leaves beneath, geometric in the variance forest,
/// and the state's other trees go in as stored. Without force the same state
/// is declined and the sampler left as it was; on a sampler with no
/// restriction it is clean. The mean tree, in pre-order, is x0 at the root,
/// x2 on its left over a leaf and a split on x1, a leaf on its right; the
/// variance tree x0 over a split on x1 and a leaf.
static void testForcedConstraintRepair() {
  std::uint64_t savedRngState = rngState;
  rngState = 424243u;
  const size_t n = 240, p = 3, numTrees = 4, numVarianceTrees = 3;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) x[i + j * n] = runif01();
    y[i] = x[i] + x[i + n] + 0.1 * (runif01() - 0.5);
  }
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr));
    ext_rng_setSeed(rngs.back(), seed);
    return &rngs.back();
  };
  const std::vector<size_t> zeroAndOne = {0, 1}, zero = {0};
  const size_t pair01[] = {0, 1};
  enum Restriction { unrestricted, columnList, maxOrder, forbiddenPair,
                     varianceColumns };
  auto make = [&](std::uint32_t seed, Restriction restriction) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numVarianceTrees = numVarianceTrees;
    if (restriction == columnList) {
      options.forestColumns = zeroAndOne.data();
      options.numForestColumns = zeroAndOne.size();
    } else if (restriction == maxOrder) {
      options.interactionMaxOrder = 1;
    } else if (restriction == forbiddenPair) {
      options.interactionForbiddenPairs = pair01;
      options.interactionNumForbiddenPairs = 1;
    } else if (restriction == varianceColumns) {
      options.varianceForestColumns = zero.data();
      options.numVarianceForestColumns = zero.size();
    }
    auto sampler = std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, newRng(seed));
    Results empty;
    sampler->run(20, 0, empty);
    return sampler;
  };
  auto plain = make(11, unrestricted);
  const ColumnStore& store(plain->data());
  auto split = [&](size_t j) {
    FlatNode node;
    node.variable = static_cast<int32_t>(j);
    node.value = store.cutPoints[j][store.numCuts[j] / 2];
    setFlatKind(node, FlatKind::ordinal);
    return node;
  };
  auto leaf = [](double value) {
    FlatNode node;
    node.value = value;
    return node;
  };
  const double v1 = 0.03, v2 = -0.05, v3 = 0.11, v4 = -0.02;
  const double h1 = 0.7, h2 = 1.9, h3 = 1.1;
  const std::vector<FlatNode> meanTree = {split(0), split(2), leaf(v1),
                                          split(1), leaf(v2), leaf(v3),
                                          leaf(v4)};
  const std::vector<FlatNode> varianceTree = {split(0), split(1), leaf(h1),
                                              leaf(h2), leaf(h3)};
  // the rows each node of a tree holds, in pre-order, by the engine's own
  // routing
  auto rowCounts = [&](const std::vector<FlatNode>& flat) {
    std::vector<index_t> indices(n);
    std::vector<double> params;
    std::vector<FlatNode> nodes;
    std::vector<std::uint32_t> counts;
    Tree scratch;
    scratch.initialize(indices.data(), n);
    built(scratch.buildFromFlat(store, flat.data(), flat.size(), params));
    scratch.repartitionSubtree(store, 0);
    scratch.flatten(store, params.data(), nodes, &counts);
    return counts;
  };
  const std::vector<std::uint32_t> meanCounts = rowCounts(meanTree),
    varianceCounts = rowCounts(varianceTree);
  check(meanCounts[2] > 0 && meanCounts[4] > 0 && meanCounts[5] > 0 &&
          meanCounts[6] > 0 && meanCounts[2] != meanCounts[4] &&
          varianceCounts[2] > 0 && varianceCounts[3] > 0 &&
          varianceCounts[4] > 0 && varianceCounts[2] != varianceCounts[3],
        "forced repair: every leaf of the hand trees holds rows, unequally");
  auto weighted = [](std::initializer_list<std::pair<double, double>> terms) {
    double total = 0.0, weight = 0.0;
    for (const std::pair<double, double>& term : terms) {
      total += term.first * term.second;
      weight += term.first;
    }
    return total / weight;
  };
  double n1 = meanCounts[2], n2 = meanCounts[4], n3 = meanCounts[5];
  const std::vector<FlatNode> belowRoot = {
    split(0), leaf(weighted({{n1, v1}, {n2, v2}, {n3, v3}})), leaf(v4)};
  const std::vector<FlatNode> belowSecond = {
    split(0), split(2), leaf(v1), leaf(weighted({{n2, v2}, {n3, v3}})),
    leaf(v4)};
  const std::vector<FlatNode> varianceBelowRoot = {
    split(0),
    leaf(std::exp(weighted({{varianceCounts[2], std::log(h1)},
                            {varianceCounts[3], std::log(h2)}}))),
    leaf(h3)};
  auto sameTree = [](const std::vector<FlatNode>& got,
                     const std::vector<FlatNode>& want) {
    if (got.size() != want.size()) return false;
    for (size_t i = 0; i < got.size(); ++i) {
      if (got[i].variable != want[i].variable) return false;
      bool isLeaf = want[i].variable == invalidVariable;
      if (isLeaf ? std::fabs(got[i].value - want[i].value) > 1e-12
                 : got[i].value != want[i].value)
        return false;
    }
    return true;
  };

  // the hand trees stand as they are where nothing forbids them
  {
    SamplerStateData state;
    plain->getState(state);
    state.chains[0].forests[0].trees[0] = meanTree;
    state.chains[0].varianceTrees[0] = varianceTree;
    SamplerStateData installed;
    bool exact = restoresExactly(*plain, state);
    plain->getState(installed);
    check(exact && sameTree(installed.chains[0].forests[0].trees[0],
                            meanTree) &&
            sameTree(installed.chains[0].varianceTrees[0], varianceTree),
          "forced repair: an unrestricted sampler takes the trees as stored");
  }

  struct Case {
    Restriction restriction;
    bool inVariance;
    const std::vector<FlatNode>* expected;
    const char* label;
  };
  const Case cases[] = {
    {columnList, false, &belowRoot, "a split off the forest's column list"},
    {maxOrder, false, &belowRoot, "a second variable under max order 1"},
    {forbiddenPair, false, &belowSecond, "a forbidden pair on one path"},
    {varianceColumns, true, &varianceBelowRoot,
     "a variance split off the variance columns"},
  };
  std::uint32_t seed = 21;
  for (const Case& item : cases) {
    auto sampler = make(seed++, item.restriction);
    SamplerStateData state, repaired;
    sampler->getState(state);
    if (item.inVariance)
      state.chains[0].varianceTrees[0] = varianceTree;
    else
      state.chains[0].forests[0].trees[0] = meanTree;
    char line[160];
    snprintf(line, sizeof line,
             "forced repair, %s: declined without force, the sampler "
             "untouched", item.label);
    check(declinesUntouched(*sampler, state), line);
    bool installs = restoresAltered(*sampler, state);
    sampler->getState(repaired);
    const std::vector<std::vector<FlatNode>>& got(
      item.inVariance ? repaired.chains[0].varianceTrees
                      : repaired.chains[0].forests[0].trees);
    const std::vector<std::vector<FlatNode>>& stored(
      item.inVariance ? state.chains[0].varianceTrees
                      : state.chains[0].forests[0].trees);
    snprintf(line, sizeof line,
             "forced repair, %s: collapsed at the first offending split, the "
             "split above it kept, the leaf at the weighted mean", item.label);
    check(installs && sameTree(got[0], *item.expected), line);
    bool othersStored = true;
    for (size_t t = 1; t < got.size(); ++t)
      othersStored = othersStored && sameTree(got[t], stored[t]);
    const std::vector<std::vector<FlatNode>>& otherGot(
      item.inVariance ? repaired.chains[0].forests[0].trees
                      : repaired.chains[0].varianceTrees);
    const std::vector<std::vector<FlatNode>>& otherStored(
      item.inVariance ? state.chains[0].forests[0].trees
                      : state.chains[0].varianceTrees);
    snprintf(line, sizeof line,
             "forced repair, %s: every other tree goes in as stored",
             item.label);
    check(othersStored && sameFlatTrees(otherGot, otherStored), line);
    snprintf(line, sizeof line,
             "forced repair, %s: the repaired sampler's own state is clean "
             "and it runs on", item.label);
    Results empty;
    bool clean = restoresExactly(*sampler, repaired);
    sampler->run(5, 0, empty);
    check(clean && liveTreesAreOccupied(*sampler), line);
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a forced install collapses at the first forbidden split\n");
}

/// A state installs whatever store its kept draws came from. A store of the
/// sampler's own size is copied slot for slot, write position included; from
/// any other the newest recorded draws that fit go to the first slots, oldest
/// first, the capacity unchanged, and every such install is clean. Saved-tree
/// blocks that name different draw counts stay refused.
static void testStateStoreSizes() {
  std::uint64_t savedRngState = rngState;
  rngState = 868687u;
  const size_t n = 200, nTest = 12, p = 2, numTrees = 6;
  std::vector<double> x, y, xTest(nTest * p);
  makeMutationData(x, y, n);
  for (double& v : xTest) v = runif01();
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint32_t seed) {
    rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr));
    ext_rng_setSeed(rngs.back(), seed);
    return rngs.back();
  };
  const size_t covariates[] = {1};
  Results none;
  // a sampler keeping `capacity` draws that has recorded `draws` of them
  auto make = [&]<typename L>(size_t capacity, size_t draws,
                              std::uint32_t seed, size_t numChains = 1,
                              size_t numVarianceTrees = 0) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.numChains = numChains;
    options.numVarianceTrees = numVarianceTrees;
    options.keepTrees = capacity > 0;
    options.numSamplesToStore = capacity;
    if constexpr (L::hasVectorParams) {
      options.leafCovariateColumns = covariates;
      options.numLeafCovariates = 1;
    }
    std::vector<ext_rng*> chainRngs;
    for (size_t c = 0; c < numChains; ++c)
      chainRngs.push_back(newRng(seed + static_cast<std::uint32_t>(c)));
    auto sampler = std::make_unique<Sampler<L>>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, chainRngs.data());
    sampler->run(10, draws, none);
    return sampler;
  };
  // the kept draws' predictions, chain-major, oldest draw first
  auto kept = [&](auto& sampler) {
    std::vector<double> out(nTest * sampler.filledSavedDraws() *
                            sampler.numChains());
    if (!out.empty()) sampler.predict(xTest.data(), nTest, 1, out.data());
    return out;
  };
  auto keptVariance = [&](auto& sampler) {
    std::vector<double> out(nTest * sampler.filledSavedDraws());
    if (!out.empty())
      sampler.predictVariance(xTest.data(), nTest, 1, out.data());
    return out;
  };
  // the last `count` draws of one chain's slab
  auto newest = [&](const std::vector<double>& draws, size_t count) {
    return std::vector<double>(draws.end() - nTest * count, draws.end());
  };
  using Constant = ConstantGaussianLeaf;

  // a wrapped ring of ten: thirteen draws leave the cursor at slot 3
  auto source = make.operator()<Constant>(10, 13, 100);
  SamplerStateData ten;
  source->getState(ten);
  const std::vector<double> tenDraws = kept(*source);
  check(ten.currentSampleNum == 3 && ten.recordedDraws == 10,
        "store sizes: the source ring has wrapped");

  // ten into four: the newest four, in order, at the store's start
  {
    auto four = make.operator()<Constant>(4, 3, 110);
    check(installsClean(*four, ten) && four->savedTreeCapacity() == 4 &&
            four->filledSavedDraws() == 4 && four->currentSampleNum() == 0 &&
            four->savedSlotForDraw(0) == 0,
          "store sizes: ten into four is clean and fills the store");
    check(kept(*four) == newest(tenDraws, 4),
          "store sizes: ten into four keeps the newest four, in order");
    SamplerStateData reinstalled;
    four->getState(reinstalled);
    check(restoresExactly(*four, reinstalled) &&
            kept(*four) == newest(tenDraws, 4),
          "store sizes: the smaller store's own state installs as stored");
    four->run(0, 1, none);
    std::vector<double> after = kept(*four);
    check(four->filledSavedDraws() == 4 &&
            std::equal(after.begin(), after.end() - nTest,
                       tenDraws.end() - nTest * 3),
          "store sizes: the next draw displaces the oldest of the four");
  }
  // a ring of four, wrapped, into ten: all four at the start, then it fills
  {
    auto small = make.operator()<Constant>(4, 6, 120);
    SamplerStateData four;
    small->getState(four);
    const std::vector<double> fourDraws = kept(*small);
    auto large = make.operator()<Constant>(10, 7, 121);
    check(four.currentSampleNum == 2 && installsClean(*large, four) &&
            large->savedTreeCapacity() == 10 &&
            large->filledSavedDraws() == 4 && large->currentSampleNum() == 4,
          "store sizes: four into ten is clean and holds four");
    check(kept(*large) == fourDraws,
          "store sizes: four into ten keeps a wrapped ring's draws in order");
    large->run(0, 3, none);
    std::vector<double> after = kept(*large);
    check(large->filledSavedDraws() == 7 &&
            std::equal(fourDraws.begin(), fourDraws.end(), after.begin()),
          "store sizes: three more sweeps make seven, the four leading");

    // four into none, and none into four
    auto bare = make.operator()<Constant>(0, 0, 122);
    check(installsClean(*bare, four) && bare->savedTreeCapacity() == 0 &&
            bare->filledSavedDraws() == 0,
          "store sizes: a sampler keeping no draws takes none");
    SamplerStateData bareState;
    bare->getState(bareState);
    check(installsClean(*small, bareState) && small->savedTreeCapacity() == 4 &&
            small->filledSavedDraws() == 0 && small->currentSampleNum() == 0,
          "store sizes: a state with no kept draws leaves the store empty");
  }
  // a store still filling: three recorded draws into four, and into two
  {
    auto filling = make.operator()<Constant>(10, 3, 130);
    SamplerStateData three;
    filling->getState(three);
    const std::vector<double> threeDraws = kept(*filling);
    auto four = make.operator()<Constant>(4, 4, 131);
    auto two = make.operator()<Constant>(2, 2, 132);
    check(installsClean(*four, three) && four->filledSavedDraws() == 3 &&
            four->currentSampleNum() == 3 && kept(*four) == threeDraws,
          "store sizes: three recorded draws into four hold three");
    check(installsClean(*two, three) && two->filledSavedDraws() == 2 &&
            two->currentSampleNum() == 0 &&
            kept(*two) == newest(threeDraws, 2),
          "store sizes: and into two, the newest two");
  }
  // equal capacity: the ring is the state's, slot for slot
  {
    auto twin = make.operator()<Constant>(10, 4, 140);
    SamplerStateData copied;
    bool exact = restoresExactly(*twin, ten);
    twin->getState(copied);
    check(exact && samplerStatesAgree(ten, copied) &&
            twin->currentSampleNum() == 3 && twin->savedSlotForDraw(0) == 3,
          "store sizes: an equal store continues the ring slot for slot");
    source->run(0, 5, none);
    twin->run(0, 5, none);
    SamplerStateData sourceAfter, twinAfter;
    source->getState(sourceAfter);
    twin->getState(twinAfter);
    check(samplerStatesAgree(sourceAfter, twinAfter),
          "store sizes: five sweeps on, the two store the same state");
  }
  // a capacity adopted at the install: the state's own continues the ring,
  // another takes the newest that fit
  {
    auto adopting = make.operator()<Constant>(4, 2, 150);
    bool notClean = true;
    check(adopting->setState(ten, nullptr, nullptr, nullptr, 10, false,
                             &notClean) &&
            !notClean && adopting->savedTreeCapacity() == 10 &&
            adopting->currentSampleNum() == 3 && kept(*adopting) == tenDraws,
          "store sizes: adopting the state's capacity continues its ring");
    auto six = make.operator()<Constant>(4, 2, 151);
    check(six->setState(ten, nullptr, nullptr, nullptr, 6) &&
            six->savedTreeCapacity() == 6 && six->filledSavedDraws() == 6 &&
            six->currentSampleNum() == 0 && kept(*six) == newest(tenDraws, 6),
          "store sizes: adopting another capacity keeps the newest that fit");
  }
  // two chains, each chain's newest draws
  {
    auto wide = make.operator()<Constant>(10, 13, 160, 2);
    SamplerStateData state;
    wide->getState(state);
    const std::vector<double> wideDraws = kept(*wide);
    auto narrow = make.operator()<Constant>(4, 1, 162, 2);
    bool clean = installsClean(*narrow, state);
    const std::vector<double> narrowDraws = kept(*narrow);
    bool newestKept = narrowDraws.size() == 2 * nTest * 4;
    for (size_t c = 0; newestKept && c < 2; ++c)
      newestKept = std::equal(
        narrowDraws.begin() + c * nTest * 4,
        narrowDraws.begin() + (c + 1) * nTest * 4,
        wideDraws.begin() + c * nTest * 10 + nTest * 6);
    check(clean && newestKept,
          "store sizes: each of two chains keeps its own newest four");
  }
  // a linear leaf's saved slopes travel with their trees
  {
    using Linear = LinearGaussianLeaf;
    auto wide = make.operator()<Linear>(10, 13, 170);
    SamplerStateData state;
    wide->getState(state);
    auto narrow = make.operator()<Linear>(4, 2, 171);
    check(installsClean(*narrow, state) &&
            kept(*narrow) == newest(kept(*wide), 4),
          "store sizes: a linear leaf's newest four predict as the source's");
  }
  // the variance forest's saved trees travel with the mean forest's, and a
  // block naming another draw count is refused in both forms
  {
    auto wide = make.operator()<Constant>(10, 13, 180, 1, 4);
    SamplerStateData state;
    wide->getState(state);
    auto narrow = make.operator()<Constant>(4, 2, 181, 1, 4);
    check(installsClean(*narrow, state) &&
            kept(*narrow) == newest(kept(*wide), 4) &&
            keptVariance(*narrow) == newest(keptVariance(*wide), 4),
          "store sizes: the variance draws kept are the mean draws' own");
    auto bare = make.operator()<Constant>(0, 0, 182, 1, 4);
    SamplerStateData bareState;
    bare->getState(bareState);
    check(installsClean(*bare, state) && bare->filledSavedDraws() == 0 &&
            installsClean(*narrow, bareState) &&
            narrow->filledSavedDraws() == 0,
          "store sizes: a heteroscedastic sampler takes none, and a state "
          "with none");
    SamplerStateData stripped(state), ragged(state), before, after;
    stripped.chains[0].savedVarianceTrees.clear();
    ragged.chains[0].forests[0].savedTrees.pop_back();
    narrow->getState(before);
    bool notClean = true;
    bool strippedRefused =
      !narrow->setState(stripped, nullptr) &&
      !narrow->setState(stripped, nullptr, nullptr, nullptr,
                        keepStoreCapacity, false, &notClean) &&
      !notClean;
    notClean = true;
    bool raggedRefused =
      !narrow->setState(ragged, nullptr) &&
      !narrow->setState(ragged, nullptr, nullptr, nullptr, keepStoreCapacity,
                        false, &notClean) &&
      !notClean;
    narrow->getState(after);
    check(strippedRefused && raggedRefused && samplerStatesAgree(before, after),
          "store sizes: blocks naming different draw counts are refused, "
          "forced or not, the sampler untouched");
  }
  // one capacity is asked of a state across its chains and its mean forests:
  // a block holding another whole number of draws is refused in both forms
  {
    auto wide = make.operator()<Constant>(10, 13, 200, 2);
    SamplerStateData state;
    wide->getState(state);
    state.chains[1].forests[0].savedTrees.resize(4 * numTrees);
    check(refusesUntouched(*wide, state),
          "store sizes: a second chain holding another number of draws is "
          "refused, forced or not, the sampler untouched");

    std::vector<double> z(n);
    for (size_t i = 0; i < n; ++i) z[i] = i % 2 == 0 ? 1.0 : 0.0;
    SamplerOptions options;
    options.keepTrees = true;
    options.numSamplesToStore = 10;
    AmplitudeSpec spec;
    spec.mu.numTrees = numTrees;
    spec.tau.numTrees = 4;
    spec.z = z.data();
    ext_rng* rng = newRng(210);
    Sampler<Constant> twoForests(x.data(), y.data(), n, p, nullptr, nullptr,
                                 1.0, 3.0, 0.37804942330213542, options, spec,
                                 &rng);
    twoForests.run(10, 13, none);
    SamplerStateData own;
    twoForests.getState(own);
    SamplerStateData shorter(own);
    shorter.chains[0].forests[1].savedTrees.resize(4 * 4);
    check(own.chains[0].forests[1].savedTrees.size() == 10 * 4 &&
            restoresExactly(twoForests, own) &&
            refusesUntouched(twoForests, shorter),
          "store sizes: a second forest holding another number of draws is "
          "refused, forced or not, the sampler untouched");
  }
  // a pooled factor's saved rules keep their words in a side channel, which
  // travels with the trees it belongs to, mean and variance
  {
    const size_t numLevels = 70, numVarianceTrees = 4;
    std::vector<double> xWide(x), yWide(n), xTestWide(xTest);
    for (size_t i = 0; i < n; ++i) {
      xWide[n + i] = static_cast<double>(i % numLevels);
      yWide[i] = (xWide[i] > 0.5 ? 2.0 : 0.0) +
        (i % numLevels < 30 ? 1.5 : -1.5) +
        (i % numLevels % 2 == 0 ? 2.0 : 0.1) * (runif01() - 0.5);
    }
    for (size_t i = 0; i < nTest; ++i)
      xTestWide[nTest + i] = static_cast<double>((7 * i) % numLevels);
    static const ColumnKind kinds[p] = {ColumnKind::numeric,
                                        ColumnKind::categorical};
    auto makeWide = [&](size_t capacity, size_t draws, std::uint32_t seed) {
      SamplerOptions options;
      options.numTrees = numTrees;
      options.numVarianceTrees = numVarianceTrees;
      options.keepTrees = true;
      options.numSamplesToStore = capacity;
      options.predictors.columnTypes = kinds;
      ext_rng* rng = newRng(seed);
      auto sampler = std::make_unique<Sampler<Constant>>(
        xWide.data(), yWide.data(), n, p, nullptr, nullptr,
        ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options,
        &rng);
      sampler->run(40, draws, none);
      return sampler;
    };
    auto keptWide = [&](auto& sampler, bool variance) {
      std::vector<double> out(nTest * sampler.filledSavedDraws());
      if (variance)
        sampler.predictVariance(xTestWide.data(), nTest, 1, out.data());
      else
        sampler.predict(xTestWide.data(), nTest, 1, out.data());
      return out;
    };
    // whether draws [first, first + count) of `got`, oldest first from slot
    // 0, hold the mask words of the draws at `slots` of `source`; numWords
    // counts them
    auto sameMasks = [&](const SamplerStateData& got,
                         const SamplerStateData& source,
                         const std::vector<size_t>& slots, size_t numWords[2]) {
      const ChainStateData& a(got.chains[0]);
      const ChainStateData& b(source.chains[0]);
      bool same = true;
      numWords[0] = numWords[1] = 0;
      for (size_t d = 0; d < slots.size(); ++d) {
        for (size_t t = 0; t < numTrees; ++t) {
          const std::vector<std::uint64_t>& words(
            a.forests[0].savedTreeMasks[d * numTrees + t]);
          same = same &&
            words == b.forests[0].savedTreeMasks[slots[d] * numTrees + t];
          numWords[0] += words.size();
        }
        for (size_t j = 0; j < numVarianceTrees; ++j) {
          const std::vector<std::uint64_t>& words(
            a.savedVarianceTreeMasks[d * numVarianceTrees + j]);
          same = same &&
            words == b.savedVarianceTreeMasks[slots[d] * numVarianceTrees + j];
          numWords[1] += words.size();
        }
      }
      return same;
    };
    size_t numWords[2];
    auto wide = makeWide(10, 13, 220);
    SamplerStateData ten, got;
    wide->getState(ten);
    auto narrow = makeWide(4, 2, 221);
    bool clean = installsClean(*narrow, ten);
    narrow->getState(got);
    // the masks first: a rule replayed against another tree's words reads
    // past them
    bool masksKept = sameMasks(got, ten, {9, 0, 1, 2}, numWords);
    check(clean && ten.currentSampleNum == 3 && masksKept && numWords[0] > 0 &&
            numWords[1] > 0,
          "store sizes: a pooled factor's newest four keep their mask words, "
          "mean and variance");
    if (masksKept)
      check(keptWide(*narrow, false) == newest(keptWide(*wide, false), 4) &&
              keptWide(*narrow, true) == newest(keptWide(*wide, true), 4),
            "store sizes: and predict as the source's newest four");

    auto small = makeWide(4, 6, 222);
    SamplerStateData four;
    small->getState(four);
    auto large = makeWide(10, 7, 223);
    clean = installsClean(*large, four);
    large->getState(got);
    masksKept = sameMasks(got, four, {2, 3, 0, 1}, numWords);
    check(clean && four.currentSampleNum == 2 && masksKept &&
            numWords[0] > 0 && numWords[1] > 0,
          "store sizes: a pooled factor's four into ten keep their mask "
          "words, mean and variance");
    if (masksKept)
      check(keptWide(*large, false) == keptWide(*small, false) &&
              keptWide(*large, true) == keptWide(*small, true),
            "store sizes: and predict as the source's four");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a state installs the newest kept draws that fit\n");
}

/// The first ordinal split among flat trees, null where none splits on a
/// threshold.
static FlatNode* firstOrdinalSplit(std::vector<std::vector<FlatNode>>& trees) {
  for (std::vector<FlatNode>& tree : trees)
    for (FlatNode& node : tree)
      if (node.variable != invalidVariable &&
          flatKindOf(node) == FlatKind::ordinal)
        return &node;
  return nullptr;
}

/// Chain::liveThresholdsAreOnGrid reads the live mean and variance trees of a
/// chain state against a grid by value: a grid without one threshold fails
/// for either kind of tree alone, a finer grid holding every threshold
/// passes, and leaves, subset rules and kept draws are not read.
static void testLiveThresholdsOnGrid() {
  std::uint64_t savedRngState = rngState;
  rngState = 424242u;
  const size_t n = 200, p = 3;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    x[i + 2 * n] = static_cast<double>(i % 4);
    y[i] = 4.0 * (x[i] > 0.5 ? 1.0 : 0.0) + x[i + n] + (i % 4 == 2 ? 1.0 : 0.0) +
      (x[i] > 0.5 ? 2.0 : 0.2) * (runif01() - 0.5);
  }
  ColumnKind types[] = {ColumnKind::numeric, ColumnKind::numeric,
                        ColumnKind::categorical};
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 4242u);
  SamplerOptions options;
  options.numTrees = 10;
  options.numVarianceTrees = 6;
  options.keepTrees = true;
  options.numSamplesToStore = 3;
  options.predictors.columnTypes = types;
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(60, 3, empty);
  SamplerStateData own;
  sampler.getState(own);
  const auto& chain = sampler.chain(0);
  ChainStateData& state = own.chains[0];
  FlatNode* meanSplit = firstOrdinalSplit(state.forests[0].trees);
  FlatNode* varianceSplit = firstOrdinalSplit(state.varianceTrees);
  bool subsetRule = false;
  for (const std::vector<FlatNode>& tree : state.forests[0].trees)
    for (const FlatNode& node : tree) subsetRule = subsetRule || node.variable == 2;
  check(meanSplit != nullptr && varianceSplit != nullptr && subsetRule &&
          own.cutPoints[2].empty() && !state.forests[0].savedTrees.empty(),
        "thresholds on grid: the fixture splits by threshold in both kinds of "
        "tree, by subset, and keeps draws");
  check(chain.liveThresholdsAreOnGrid(state, own.cutPoints),
        "thresholds on grid: a state's thresholds are points of its own grid");

  // one point taken out, or moved by an ulp, fails exactly the trees on it
  auto without = [&](const FlatNode& split, bool nudge) {
    std::vector<std::vector<double>> grid(own.cutPoints);
    std::vector<double>& cuts = grid[static_cast<size_t>(split.variable)];
    auto point = std::find(cuts.begin(), cuts.end(), split.value);
    if (nudge)
      *point = std::nextafter(*point, std::numeric_limits<double>::infinity());
    else
      cuts.erase(point);
    return grid;
  };
  check(!chain.liveThresholdsAreOnGrid(state, without(*meanSplit, false)) &&
          !chain.liveThresholdsAreOnGrid(state, without(*meanSplit, true)),
        "thresholds on grid: a mean threshold the grid lacks, or holds an ulp "
        "away, is off it");
  // a variance threshold no mean tree holds: moved in the state to a value
  // its grid gains, so the mean trees stay on both grids
  double added = 0.5 * (own.cutPoints[0][0] + own.cutPoints[0][1]);
  SamplerStateData moved(own);
  FlatNode* movedSplit = firstOrdinalSplit(moved.chains[0].varianceTrees);
  movedSplit->variable = 0;
  movedSplit->value = added;
  std::vector<std::vector<double>> finer(own.cutPoints);
  finer[0].insert(finer[0].begin() + 1, added);
  check(chain.liveThresholdsAreOnGrid(moved.chains[0], finer) &&
          !chain.liveThresholdsAreOnGrid(moved.chains[0], own.cutPoints) &&
          chain.liveThresholdsAreOnGrid(state, finer),
        "thresholds on grid: a variance tree's threshold alone is read, and a "
        "finer grid holds every stored one");

  // not read: a kept draw's threshold, a leaf's value and a subset rule's mask
  SamplerStateData kept(own);
  FlatNode* keptSplit = firstOrdinalSplit(kept.chains[0].forests[0].savedTrees);
  FlatNode* keptVarianceSplit =
    firstOrdinalSplit(kept.chains[0].savedVarianceTrees);
  check(keptSplit != nullptr && keptVarianceSplit != nullptr,
        "thresholds on grid: the kept draws split by threshold");
  keptSplit->value = added;
  keptVarianceSplit->value = added;
  check(chain.liveThresholdsAreOnGrid(kept.chains[0], own.cutPoints),
        "thresholds on grid: kept draws are held to no grid");
  // an ordinal record on a column whose entry is empty holds no point of it
  SamplerStateData mistagged(own);
  firstOrdinalSplit(mistagged.chains[0].forests[0].trees)->variable = 2;
  check(!chain.liveThresholdsAreOnGrid(mistagged.chains[0], own.cutPoints),
        "thresholds on grid: a threshold on a column with no points is off "
        "its grid");

  ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: a chain state's live thresholds are read against a grid\n");
}

/// One Sampler::setState call that keeps the store's grid, its answers
/// preset to the ones that would hide a call never writing them.
struct LiveInstall {
  bool installed, notClean, altered, gridDiffers;
  template <typename S>
  LiveInstall(S& sampler, const SamplerStateData& state, const double* x,
              bool force)
    : notClean(!force), altered(!force), gridDiffers(false) {
    installed = sampler.setState(state, x, nullptr, &altered,
                                 keepStoreCapacity, force, &notClean,
                                 &gridDiffers);
  }
  bool clean() const { return installed && !notClean && !altered; }
  bool declined() const { return !installed && notClean; }
  bool refused() const { return !installed && !notClean; }
};

/// Whether both forms of an install that keeps the store's grid refuse a
/// state and leave the sampler, its grid, weights and codes as they were.
template <typename S>
static bool refusesLiveUntouched(S& sampler, const SamplerStateData& state,
                                 const double* x) {
  SamplerStateData before, after;
  sampler.getState(before);
  std::vector<xint_t> codes = storageDigest(sampler.data());
  bool refused = LiveInstall(sampler, state, x, false).refused() &&
    LiveInstall(sampler, state, x, true).refused();
  sampler.getState(after);
  return refused && samplerStatesAgree(before, after) &&
    storageDigest(sampler.data()) == codes;
}

/// Sampler::setState onto the store's grid: its grid, weights and codes are
/// never written. A state on the store's own grid installs as stored, as into
/// a twin; one whose thresholds the grid holds installs as stored
/// whatever grid or weights it records; one with a threshold the grid lacks,
/// in a mean or a variance tree of any chain, is declined without force and
/// with it leaves what installing it on its own grid and then calling
/// setCutPoints by value leaves, kept draws untouched. A threshold on neither
/// grid and another grid on an ordered factor are refused.
static void testLiveGridInstall() {
  std::uint64_t savedRngState = rngState;
  rngState = 271828u;
  const size_t n = 200, p = 2, numChains = 2;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = 4.0 * (x[i] > 0.5 ? 1.0 : 0.0) + x[i + n] +
      (x[i] > 0.5 ? 2.0 : 0.2) * (runif01() - 0.5);
  SamplerOptions options;
  options.numTrees = 10;
  options.numVarianceTrees = 6;
  options.numChains = numChains;
  options.keepTrees = true;
  options.numSamplesToStore = 3;
  std::vector<ext_rng*> rngs;
  auto make = [&]() {
    size_t first = rngs.size();
    for (size_t c = 0; c < numChains; ++c) {
      rngs.push_back(
        ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr));
      ext_rng_setSeed(rngs.back(), 577u + static_cast<std::uint_least32_t>(c));
    }
    auto sampler = std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, rngs.data() + first);
    Results empty;
    sampler->run(60, 3, empty);
    return sampler;
  };
  rngs.reserve(4 * numChains);
  auto live = make(), twin = make();
  ConstantLeafSampler& sampler = *live;
  SamplerStateData own, got, other;
  sampler.getState(own);
  auto sigmaDraws = [](ConstantLeafSampler& from) {
    std::vector<double> draws(5 * numChains);
    Results results;
    results.sigma = draws.data();
    from.run(0, 5, results);
    return draws;
  };
  auto onGrid = [](ConstantLeafSampler& of, const SamplerStateData& state) {
    for (size_t c = 0; c < state.chains.size(); ++c)
      if (!of.chain(c).liveThresholdsAreOnGrid(state.chains[c],
                                               of.data().cutPoints))
        return false;
    return true;
  };

  // its own grid: what the twin's install does, and no difference reported
  Results none;
  sampler.run(5, 0, none);
  twin->run(5, 0, none);
  LiveInstall same(sampler, own, nullptr, false);
  sampler.getState(got);
  check(same.clean() && !same.gridDiffers && samplerStatesAgree(own, got) &&
          restoresExactly(*twin, own),
        "live grid: a state on the store's grid installs as stored");
  check(sigmaDraws(sampler) == sigmaDraws(*twin),
        "live grid: and continues as the twin continues");
  check(LiveInstall(sampler, own, nullptr, false).clean() &&
          restoresExactly(*twin, own),
        "live grid: both go back to the stored state");

  // three of column 0's points: the stored thresholds on it are mostly gone
  const std::vector<double> fine(own.cutPoints[0]);
  check(fine.size() > 80, "live grid: the fixture's first column holds a grid");
  const std::vector<double> coarse = {fine[24], fine[49], fine[74]};
  const size_t column = 0;
  auto setGrid = [&](ConstantLeafSampler& on, const std::vector<double>& grid,
                     SplitPlacement placement) {
    const double* grids[] = {grid.data()};
    std::uint32_t counts[] = {static_cast<std::uint32_t>(grid.size())};
    on.setCutPoints(grids, counts, &column, 1, x.data(), placement);
  };
  setGrid(sampler, coarse, SplitPlacement::byPosition);
  SamplerStateData onCoarse;
  sampler.getState(onCoarse);
  std::vector<xint_t> coarseCodes = storageDigest(sampler.data());
  check(!onGrid(sampler, own) && onGrid(sampler, onCoarse),
        "live grid: the stored state holds thresholds the coarse grid lacks");
  LiveInstall unforced(sampler, own, x.data(), false);
  sampler.getState(got);
  check(unforced.declined() && samplerStatesAgree(onCoarse, got) &&
          storageDigest(sampler.data()) == coarseCodes,
        "live grid: a threshold off the grid declines, nothing written");

  LiveInstall forced(sampler, own, x.data(), true);
  sampler.getState(got);
  check(forced.installed && forced.altered && forced.gridDiffers &&
          got.cutPoints == onCoarse.cutPoints &&
          got.cutMass == onCoarse.cutMass &&
          storageDigest(sampler.data()) == coarseCodes,
        "live grid: forced, the store's grid, weights and codes are as they "
        "were");
  bool keptAsStored = true;
  for (size_t c = 0; c < numChains; ++c)
    keptAsStored = keptAsStored &&
      sameFlatTrees(got.chains[c].forests[0].savedTrees,
                    own.chains[c].forests[0].savedTrees) &&
      sameFlatTrees(got.chains[c].savedVarianceTrees,
                    own.chains[c].savedVarianceTrees);
  check(onGrid(sampler, got) && liveTreesAreOccupied(sampler) &&
          !sameLiveTrees(got, own) && keptAsStored &&
          got.currentSampleNum == own.currentSampleNum &&
          got.recordedDraws == own.recordedDraws,
        "live grid: forced, the trees are moved onto the grid and the kept "
        "draws are not");
  // the twin: the state on its own grid, then the grid replaced by value
  check(restoresExactly(*twin, own), "live grid: the twin takes the state");
  setGrid(*twin, coarse, SplitPlacement::byValue);
  twin->getState(other);
  check(samplerStatesAgree(got, other) &&
          storageDigest(twin->data()) == coarseCodes,
        "live grid: a forced install across grids leaves what the install "
        "on the state's grid followed by setCutPoints by value leaves");
  check(sigmaDraws(sampler) == sigmaDraws(*twin),
        "live grid: and the two continue alike");

  // a grid holding every stored point: clean, kept, the trees as stored
  std::vector<double> finer;
  for (size_t k = 0; k + 1 < fine.size(); ++k) {
    finer.push_back(fine[k]);
    finer.push_back(0.5 * (fine[k] + fine[k + 1]));
  }
  finer.push_back(fine.back());
  setGrid(sampler, finer, SplitPlacement::byValue);
  std::vector<xint_t> finerCodes = storageDigest(sampler.data());
  LiveInstall onto(sampler, own, x.data(), false);
  sampler.getState(got);
  check(onto.clean() && onto.gridDiffers && got.cutPoints[0] == finer &&
          storageDigest(sampler.data()) == finerCodes && statesAgree(own, got),
        "live grid: a state whose thresholds the grid holds installs as "
        "stored on it");
  LiveInstall ontoForced(sampler, own, x.data(), true);
  sampler.getState(got);
  check(ontoForced.installed && !ontoForced.altered &&
          ontoForced.gridDiffers && statesAgree(own, got),
        "live grid: and forced, with nothing altered");

  // one tree of the second chain: a threshold of the recorded grid the
  // store's lacks, in a mean tree and in a variance tree each alone
  setGrid(sampler, coarse, SplitPlacement::byValue);
  SamplerStateData base;
  sampler.getState(base);
  base.cutPoints = own.cutPoints;
  check(LiveInstall(sampler, base, x.data(), false).clean(),
        "live grid: thresholds on the coarse grid, recorded on the fine one, "
        "install as stored");
  SamplerStateData meanOff(base), varianceOff(base);
  FlatNode* meanSplit = firstOrdinalSplit(meanOff.chains[1].forests[0].trees);
  FlatNode* varianceSplit = firstOrdinalSplit(varianceOff.chains[1].varianceTrees);
  check(meanSplit != nullptr && varianceSplit != nullptr,
        "live grid: the second chain splits in both kinds of tree");
  meanSplit->variable = 0;
  meanSplit->value = fine[10];
  varianceSplit->variable = 0;
  varianceSplit->value = fine[60];
  for (const SamplerStateData* state : {&meanOff, &varianceOff}) {
    sampler.getState(other);
    LiveInstall declined(sampler, *state, x.data(), false);
    sampler.getState(got);
    bool untouched = declined.declined() && samplerStatesAgree(other, got);
    LiveInstall moved(sampler, *state, x.data(), true);
    sampler.getState(got);
    check(untouched && moved.installed && moved.altered && onGrid(sampler, got) &&
            got.cutPoints[0] == coarse && liveTreesAreOccupied(sampler),
          state == &meanOff
            ? "live grid: one mean threshold of the second chain declines "
              "the state, and forced it is moved"
            : "live grid: one variance threshold of the second chain "
              "declines the state, and forced it is moved");
  }

  // a threshold that is a point of neither grid, whatever grid is recorded
  SamplerStateData neither(base), neitherEqual(base);
  firstOrdinalSplit(neither.chains[1].forests[0].trees)->variable = 0;
  firstOrdinalSplit(neither.chains[1].forests[0].trees)->value =
    0.5 * (fine[10] + fine[11]);
  neitherEqual.chains = neither.chains;
  neitherEqual.cutPoints[0] = coarse;
  check(refusesLiveUntouched(sampler, neither, x.data()) &&
          refusesLiveUntouched(sampler, neitherEqual, x.data()),
        "live grid: a threshold on neither grid is refused under the store's "
        "recorded grid and under another");
  // and one the record lacks though the store's grid holds it
  SamplerStateData recordLacks(base);
  recordLacks.cutPoints[0] = {fine[5], fine[49], fine[90]};
  bool usesOthers = false;
  for (size_t c = 0; c < numChains; ++c)
    usesOthers = usesOthers ||
      !sampler.chain(c).liveThresholdsAreOnGrid(recordLacks.chains[c],
                                                recordLacks.cutPoints);
  check(usesOthers && refusesLiveUntouched(sampler, recordLacks, x.data()),
        "live grid: a threshold the recorded grid lacks is refused though "
        "the store's holds it");

  // an install step that throws, standing on the recorded grid: the store's
  // grid, weights and codes are back when the throw leaves
  {
    const std::vector<std::vector<double>> points = sampler.data().cutPoints;
    const std::vector<std::vector<double>> mass = sampler.data().cutMass;
    std::vector<xint_t> codes = storageDigest(sampler.data());
    bool thrown = false, onRecorded = false;
    try {
      TestPeer::installAcrossGrids(sampler, own.cutPoints, x.data(),
                                   [&]() -> bool {
        onRecorded = sampler.data().cutPoints[0] == fine &&
          storageDigest(sampler.data()) != codes;
        throw std::runtime_error("install step");
      });
    } catch (const std::runtime_error&) {
      thrown = true;
    }
    check(thrown && onRecorded && points[0] == coarse &&
            sampler.data().cutPoints == points &&
            sampler.data().cutMass == mass &&
            storageDigest(sampler.data()) == codes,
          "live grid: a throw out of the install step puts the store's grid, "
          "weights and codes back");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);

  // weights are the store's: a state recording none, or others, on equal
  // points installs as stored and writes none
  {
    std::vector<double> xw(n * p);
    for (size_t i = 0; i < n; ++i) {
      xw[i] = static_cast<double>((i % 6) * (i % 6));
      xw[i + n] = x[i + n];
    }
    ext_rng* weightedRngs[2];
    std::unique_ptr<ConstantLeafSampler> weighted[2];
    for (size_t k = 0; k < 2; ++k) {
      weightedRngs[k] =
        ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
      ext_rng_setSeed(weightedRngs[k], 3141u);
      SamplerOptions weightedOptions;
      weightedOptions.numTrees = 10;
      weighted[k] = std::make_unique<ConstantLeafSampler>(
        xw.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
        1.0, 3.0, 0.37804942330213542, weightedOptions, &weightedRngs[k]);
      Results empty;
      weighted[k]->run(30, 0, empty);
    }
    const std::vector<std::vector<double>> weights = weighted[0]->data().cutMass;
    SamplerStateData withWeights, bare, after;
    weighted[0]->getState(withWeights);
    bare = withWeights;
    bare.cutMass.clear();
    LiveInstall bareUnforced(*weighted[0], bare, xw.data(), false);
    LiveInstall bareForced(*weighted[0], bare, xw.data(), true);
    weighted[0]->getState(after);
    check(weighted[0]->data().cutsWeighted(0) && bareUnforced.clean() &&
            bareUnforced.gridDiffers && bareForced.installed &&
            !bareForced.altered && bareForced.gridDiffers &&
            weighted[0]->data().cutMass == weights &&
            samplerStatesAgree(withWeights, after),
          "live grid: a state recording no weights leaves a weighted grid "
          "its weights");
    // the twin made unweighted by another grid and then its own, which
    // leaves the points and not their weights
    {
      std::vector<double> held(weighted[1]->data().cutPoints[0]);
      std::vector<double> moved(held);
      moved.back() += 1.0;
      const size_t column = 0;
      const double* movedGrid[] = {moved.data()};
      const double* heldGrid[] = {held.data()};
      std::uint32_t counts[] = {static_cast<std::uint32_t>(held.size())};
      weighted[1]->setCutPoints(movedGrid, counts, &column, 1, xw.data());
      weighted[1]->setCutPoints(heldGrid, counts, &column, 1, xw.data());
    }
    check(restoresExactly(*weighted[1], bare) &&
            !weighted[1]->data().cutsWeighted(0),
          "live grid: the twin's grid is unweighted");
    LiveInstall heavyUnforced(*weighted[1], withWeights, xw.data(), false);
    LiveInstall heavyForced(*weighted[1], withWeights, xw.data(), true);
    weighted[1]->getState(after);
    check(heavyUnforced.clean() && heavyUnforced.gridDiffers &&
            heavyForced.installed && !heavyForced.altered &&
            heavyForced.gridDiffers && !weighted[1]->data().cutsWeighted(0) &&
            samplerStatesAgree(bare, after),
          "live grid: a state recording weights leaves an unweighted grid "
          "unweighted");
    for (ext_rng* r : weightedRngs) ext_rng_destroy(r);
  }

  // an ordered factor's grid is its level count's: a state recording another
  // on that column is another column's, on thresholds both grids hold
  {
    std::vector<double> xo(n * p);
    for (size_t i = 0; i < n; ++i) {
      xo[i] = static_cast<double>(i % 4);
      xo[i + n] = x[i + n];
    }
    std::vector<double> yo(n);
    for (size_t i = 0; i < n; ++i)
      yo[i] = (i % 4 >= 2 ? 3.0 : 0.0) + xo[i + n] + 0.2 * (runif01() - 0.5);
    ColumnKind kinds[] = {ColumnKind::orderedFactor, ColumnKind::numeric};
    ext_rng* orderedRng =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(orderedRng, 1618u);
    SamplerOptions orderedOptions;
    orderedOptions.numTrees = 10;
    orderedOptions.predictors.columnTypes = kinds;
    ConstantLeafSampler ordered(xo.data(), yo.data(), n, p, nullptr, nullptr,
                                ResponseFamily::gaussian, 1.0, 3.0,
                                0.37804942330213542, orderedOptions,
                                &orderedRng);
    Results empty;
    ordered.run(30, 0, empty);
    SamplerStateData levels, fiveLevels;
    ordered.getState(levels);
    fiveLevels = levels;
    fiveLevels.cutPoints[0].push_back(3.5);
    check(levels.cutPoints[0] == std::vector<double>({0.5, 1.5, 2.5}) &&
            firstOrdinalSplit(levels.chains[0].forests[0].trees) != nullptr &&
            LiveInstall(ordered, levels, xo.data(), false).clean(),
          "live grid: an ordered factor's state installs on its own grid");
    check(refusesLiveUntouched(ordered, fiveLevels, xo.data()),
          "live grid: another grid recorded on an ordered factor is refused "
          "by an install that keeps the store's");
    ext_rng_destroy(orderedRng);
  }

  rngState = savedRngState;
  printf("ok: a state installs onto the store's cut grid\n");
}

// The two forms of a warm start on hand-built donors, each held to one answer
// both ways (warmVerdictMatchesInstall): what the report counts, the verdict
// taken after the sides are drawn, a donor that does not build refused in
// both forms, and the unforced answer against the forced install over seeded
// pairs.
static void testWarmStartForms() {
  std::uint64_t savedRngState = rngState;
  rngState = 271828u;
  const size_t n = 200, p = 2, numTrees = 4;
  const double missing = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(n * p), xOther(n * p), xHoled, y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    xOther[i] = runif01();
    xOther[i + n] = runif01();
    y[i] = 2.0 * x[i] + std::sin(4.0 * x[i + n]) + 0.2 * (runif01() - 0.5);
  }
  // column 1 with a missing value in every tenth row and one value far above
  // the rest, in the row of the largest column 0
  xHoled = x;
  size_t top = 0;
  for (size_t i = 0; i < n; ++i) {
    if (i % 10 == 3) xHoled[i + n] = missing;
    if (x[i] > x[top]) top = i;
  }
  xHoled[top + n] = 2.0;

  std::vector<ext_rng*> rngs;
  auto make = [&](const std::vector<double>& xs, std::uint32_t seed,
                  size_t maxOrder = 0) {
    SamplerOptions options;
    options.numTrees = numTrees;
    options.interactionMaxOrder = maxOrder;
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return std::make_unique<ConstantLeafSampler>(
      xs.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &r);
  };
  auto split = [](int32_t variable, double value, bool missingRight = false) {
    FlatNode node;
    node.variable = variable;
    node.value = value;
    setFlatKind(node, FlatKind::ordinal);
    if (missingRight) node.flags |= flatMissingGoesRight;
    return node;
  };
  auto leaf = [](double value) {
    FlatNode node;
    node.value = value;
    return node;
  };
  const std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  Results none;

  // --- the report, a tree counted once under the first name that holds. The
  // donor is recorded on other points of column 1: a tree as stored; one
  // moved, which also loses a side; one cut back, whose splits are off the
  // grid too; and one that only loses its side
  auto recipient = make(x, 11);
  SamplerStateData hand;
  recipient->getState(hand);
  const std::vector<double> own0(hand.cutPoints[0]), own1(hand.cutPoints[1]);
  std::vector<double>& other1(hand.cutPoints[1]);
  double half = 0.5 * (own1[1] - own1[0]);
  other1.pop_back();
  for (double& point : other1) point += half;
  std::vector<std::vector<FlatNode>>& trees(hand.chains[0].forests[0].trees);
  trees[0] = {split(0, own0[40]), leaf(0.1), leaf(-0.1)};
  trees[1] = {split(1, other1[30], true), leaf(0.2), leaf(-0.2)};
  trees[2] = {split(1, other1[20]), split(1, other1[60]), leaf(0.3),
              leaf(0.5), leaf(-0.3)};
  trees[3] = {split(0, own0[50], true), leaf(0.05), leaf(-0.05)};
  {
    WarmStartReport report;
    bool installs =
      warmStartsAltered(*recipient, hand, liveMap, x.data(), nullptr, &report);
    SamplerStateData after;
    recipient->getState(after);
    const WarmStartTreeCounts& counts(report.mean[0]);
    check(installs && counts.offered == 4 && counts.taken == 4 &&
            counts.moved == 1 && counts.cutBack == 1 &&
            counts.sideDropped == 1 &&
            !report.variance[0].changed() && !report.fitsCold &&
            !report.slopesLeftBehind &&
            sameFlatTrees({after.chains[0].forests[0].trees[0]}, {trees[0]}) &&
            after.chains[0].forests[0].trees[2].size() == 3,
          "warm start forms: the report counts each changed tree once, moved "
          "ahead of a side dropped and cut back ahead of moved");
  }
  // a donor of linear leaves into constant ones leaves its slopes behind, a
  // change of its own with no tree counted
  {
    SamplerStateData clean;
    recipient->getState(clean);
    SamplerStateData sloped = clean;
    for (ForestStateData& forest : sloped.chains[0].forests)
      for (const std::vector<FlatNode>& tree : forest.trees)
        forest.treeParams.emplace_back((tree.size() + 1) / 2, 0.25);
    DonorLeafShape linearLeaves;
    linearLeaves.kind = DonorLeafShape::Kind::linear;
    linearLeaves.numCovariates = 1;
    WarmStartReport report;
    check(warmStartsClean(*recipient, clean, liveMap, x.data()) &&
            warmStartsAltered(*recipient, sloped, liveMap, x.data(),
                              &linearLeaves, &report) &&
            report.slopesLeftBehind && !report.mean[0].changed(),
          "warm start forms: a clean linear donor into constant leaves is "
          "declined and, forced, reports its slopes left behind alone");
  }

  // --- a donor that does not build is refused in both forms, never
  // declined: a threshold that is no point of the grid the donor records,
  // beside one off the store's grid, and where every threshold is a point of
  // the store's
  {
    SamplerStateData offOwn = hand, onStore = hand;
    offOwn.chains[0].forests[0].trees[1][0].value += 0.25 * half;
    for (size_t t = 1; t < numTrees; ++t)
      onStore.chains[0].forests[0].trees[t] = {leaf(0.0)};
    onStore.chains[0].forests[0].trees[1] =
      {split(1, own1[30]), leaf(0.2), leaf(-0.2)};
    WarmStartResult result;
    bool altered;
    bool refused = true;
    for (const SamplerStateData* bad : {&offOwn, &onStore})
      refused = refused &&
        warmVerdictMatchesInstall(*recipient, *bad, liveMap, x.data(), result,
                                  altered) &&
        result == WarmStartResult::rebuildFailed;
    check(refused, "warm start forms: a threshold off the donor's own "
                   "recorded grid is the rebuild refusal in both forms, the "
                   "sampler untouched");
  }

  // --- the verdict follows the draw. The donor could hold no missing value
  // in column 1 and the recipient holds some. The lower split sits above
  // every value of column 1 its branch holds, so its right leaf is reached by
  // missing rows alone: sent right the split stands and the donor is clean,
  // sent left the leaf is empty and merges. Both are met over the seeds, and
  // the unforced answer is the forced install's each time
  {
    size_t stood = 0, merged = 0;
    bool agree = true;
    for (std::uint32_t seed = 1; seed <= 16; ++seed) {
      auto holed = make(xHoled, 100 + seed);
      SamplerStateData donor;
      holed->getState(donor);
      donor.missingColumns.assign(p, 0);
      std::vector<std::vector<FlatNode>>& given(
        donor.chains[0].forests[0].trees);
      for (std::vector<FlatNode>& tree : given) tree = {leaf(0.0)};
      given[0] = {split(0, donor.cutPoints[0][50]),
                  split(1, donor.cutPoints[1].back()), leaf(0.1), leaf(0.3),
                  leaf(-0.2)};
      WarmStartResult result;
      bool altered;
      agree = agree &&
        warmVerdictMatchesInstall(*holed, donor, liveMap, xHoled.data(),
                                  result, altered) &&
        result == WarmStartResult::ok;
      SamplerStateData after;
      holed->getState(after);
      const std::vector<FlatNode>& got(after.chains[0].forests[0].trees[0]);
      bool stands = got.size() == 5 &&
        (got[1].flags & flatMissingGoesRight) != 0;
      agree = agree && (stands ? !altered : altered && got.size() == 3);
      ++(stands ? stood : merged);
    }
    check(agree && stood > 0 && merged > 0,
          "warm start forms: the unforced answer is taken after the sides "
          "are drawn, clean where the drawn side keeps the split");
  }

  // --- the unforced answer is the forced install's over seeded pairs: the
  // donor's own rows, other rows, rows moved under the donor's grid, a limit
  // on interactions, and sides to draw
  {
    size_t taken = 0, changed = 0;
    bool agree = true;
    for (std::uint32_t seed = 0; seed < 40; ++seed) {
      bool holes = seed % 5 == 4;
      auto donor = make(holes ? xHoled : x, 200 + seed);
      donor->run(30, 0, none);
      SamplerStateData state;
      donor->getState(state);
      if (holes) state.missingColumns.assign(p, 0);
      const std::vector<double>& rows(
        seed % 5 == 1 || seed % 5 == 2 ? xOther : holes ? xHoled : x);
      auto pair = make(seed % 5 == 2 ? x : rows, 300 + seed,
                       seed % 5 == 3 ? 1 : 0);
      if (seed % 5 == 2) pair->setPredictor(xOther.data(), true, false);
      WarmStartResult result;
      bool altered;
      agree = agree &&
        warmVerdictMatchesInstall(*pair, state, liveMap, rows.data(), result,
                                  altered) &&
        result == WarmStartResult::ok;
      ++(altered ? changed : taken);
    }
    check(agree && taken >= 8 && changed >= 8,
          "warm start forms: over 40 seeded pairs the unforced call declines "
          "exactly when the forced install is altered, and the report says "
          "so");
    printf("ok: warm start forms (%zu pairs taken as stored, %zu changed)\n",
           taken, changed);
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
}

// A warm start across kinds of leaf (Sampler::installForests): each tree's
// block is held to the donor's own leaves, gp leaves go only to gp leaves on
// as many rows, a constant donor gives linear leaves slopes of 0 and is
// clean, and a linear donor on another number of covariates leaves its
// slopes behind, each leaf at the intercept it stores.
static void testWarmStartLeafKinds() {
  std::uint64_t savedRngState = rngState;
  rngState = 271828u;
  const size_t n = 80, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    y[i] = 2.0 * x[i] + x[i + n] * x[i] + 0.2 * (runif01() - 0.5);
  }
  std::vector<ext_rng*> rngs;
  const size_t one[] = {1}, two[] = {0, 1};
  auto make = [&](const size_t* covariates, size_t numCovariates, bool gp,
                  std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = 6;
    options.numChains = 1;
    options.maxNumCuts = 12;
    options.leafCovariateColumns = covariates;
    options.numLeafCovariates = numCovariates;
    options.gpLeaves = gp;
    rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
    ext_rng_setSeed(rngs.back(), seed);
    std::unique_ptr<SamplerBase> sampler = createSampler(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &rngs.back());
    Results none;
    sampler->run(40, 0, none);
    return sampler;
  };
  const std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  using Kind = DonorLeafShape::Kind;
  auto shapeOf = [](Kind kind, size_t covariates, size_t rows) {
    DonorLeafShape shape;
    shape.kind = kind;
    shape.numCovariates = covariates;
    shape.numObservations = rows;
    return shape;
  };
  // an error in both forms, the sampler as it was
  auto refused = [&](SamplerBase& sampler, const SamplerStateData& donor,
                     const DonorLeafShape& shape, WarmStartResult expected) {
    WarmStartResult result;
    bool altered;
    return warmVerdictMatchesInstall(sampler, donor, liveMap, x.data(),
                                     result, altered, &shape) &&
           result == expected;
  };

  std::unique_ptr<SamplerBase> linear = make(one, 1, false, 11);
  std::unique_ptr<SamplerBase> gp = make(one, 1, true, 12);
  SamplerStateData linearState, gpState;
  linear->getState(linearState);
  gp->getState(gpState);
  SamplerStateData longSlopes = linearState, shortFits = gpState;
  longSlopes.chains[0].forests[0].treeParams[2].push_back(0.0);
  shortFits.chains[0].forests[0].treeParams[1].pop_back();
  check(refused(*linear, longSlopes, shapeOf(Kind::linear, 1, 0),
                WarmStartResult::leafBlockMismatch) &&
          refused(*gp, shortFits, shapeOf(Kind::gp, 0, n),
                  WarmStartResult::leafBlockMismatch),
        "warm start kinds: a slope or fit block of another size than the "
        "donor's leaves give is the malformed-donor error");
  check(refused(*linear, gpState, shapeOf(Kind::gp, 0, n),
                WarmStartResult::leafKindMismatch) &&
          refused(*gp, linearState, shapeOf(Kind::linear, 1, 0),
                  WarmStartResult::leafKindMismatch) &&
          refused(*gp, gpState, shapeOf(Kind::gp, 0, n + 1),
                  WarmStartResult::leafRowsMismatch),
        "warm start kinds: gp leaves on one side, or on other rows, are an "
        "error in both forms");

  // a constant donor: slopes of 0, the trees as they stand, clean
  SamplerStateData constantState = linearState;
  for (ForestStateData& forest : constantState.chains[0].forests) {
    forest.treeParams.clear();
    forest.leafCovariateCenters.clear();
    forest.leafCovariateScales.clear();
  }
  DonorLeafShape constant = shapeOf(Kind::constant, 0, 0);
  WarmStartResult result;
  bool altered;
  WarmStartReport report;
  std::unique_ptr<SamplerBase> intoLinear = make(one, 1, false, 13);
  bool clean = warmVerdictMatchesInstall(*intoLinear, constantState, liveMap,
                                         x.data(), result, altered, &constant,
                                         &report) &&
               result == WarmStartResult::ok && !altered;
  SamplerStateData after;
  intoLinear->getState(after);
  bool zeros = true;
  for (const std::vector<double>& slopes : after.chains[0].forests[0].treeParams)
    for (double slope : slopes) zeros = zeros && slope == 0.0;
  check(clean && zeros && report.slopesStartedAtZero &&
          !report.slopesLeftBehind &&
          sameFlatTrees(after.chains[0].forests[0].trees,
                        linearState.chains[0].forests[0].trees),
        "warm start kinds: a constant donor gives linear leaves slopes of 0, "
        "taken as it stands");

  // a linear donor on two covariates into one: each leaf its stored
  // intercept, slopes of 0, declined and forced
  std::unique_ptr<SamplerBase> wide = make(two, 2, false, 14);
  SamplerStateData wideState;
  wide->getState(wideState);
  std::unique_ptr<SamplerBase> narrow = make(one, 1, false, 15);
  DonorLeafShape twoSlopes = shapeOf(Kind::linear, 2, 0);
  check(warmStartsAltered(*narrow, wideState, liveMap, x.data(), &twoSlopes,
                          &report) &&
          report.slopesLeftBehind && !report.mean[0].changed(),
        "warm start kinds: a linear donor on another number of covariates "
        "is declined and, forced, leaves its slopes behind alone");
  narrow->getState(after);
  zeros = true;
  for (const std::vector<double>& slopes : after.chains[0].forests[0].treeParams)
    for (double slope : slopes) zeros = zeros && slope == 0.0;
  check(zeros && sameFlatTrees(after.chains[0].forests[0].trees,
                               wideState.chains[0].forests[0].trees),
        "warm start kinds: and each leaf is the donor's stored intercept, "
        "not restated");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: warm start across kinds of leaf\n");
}

// A warm start from a donor of another tree count. The sides a donor grown
// without missing values owes are drawn once, over every tree of the sample a
// chain takes, before any round: a recipient of fewer trees ends with its
// generator where a twin holding all of the donor's trees ends, and the trees
// it keeps carry the twin's sides. A donor of fewer trees is clean, the rest
// single leaves; one of more is declined and, forced, altered.
static void testWarmStartTreeCounts() {
  std::uint64_t savedRngState = rngState;
  rngState = 161803u;
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), xHoled, y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    y[i] = 2.0 * x[i] + std::sin(4.0 * x[i + n]) + 0.2 * (runif01() - 0.5);
  }
  xHoled = x;
  for (size_t i = 3; i < n; i += 10)
    xHoled[i + n] = std::numeric_limits<double>::quiet_NaN();

  std::vector<ext_rng*> rngs;
  auto make = [&](const std::vector<double>& xs, std::uint32_t seed,
                  size_t numTrees) {
    SamplerOptions options;
    options.numTrees = numTrees;
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return std::make_unique<ConstantLeafSampler>(
      xs.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &r);
  };
  const std::vector<std::pair<size_t, int>> liveMap = {{0, -1}};
  Results none;

  // six trees on column 1, which the donor could hold no missing value in and
  // the recipient does, at distinct points and values
  auto source = make(xHoled, 5, 6);
  SamplerStateData donor;
  source->getState(donor);
  donor.missingColumns.assign(p, 0);
  const std::vector<double>& points(donor.cutPoints[1]);
  for (size_t t = 0; t < 6; ++t) {
    FlatNode split, low, high;
    split.variable = 1;
    split.value = points[20 + 10 * t];
    setFlatKind(split, FlatKind::ordinal);
    low.value = 0.1 * static_cast<double>(t + 1);
    high.value = -low.value;
    donor.chains[0].forests[0].trees[t] = {split, low, high};
  }

  auto recipient = make(xHoled, 41, 3);
  auto twin = make(xHoled, 42, 6);
  SamplerStateData before, twinState;
  recipient->getState(before);
  twin->getState(twinState);
  twinState.chains[0].rngState = before.chains[0].rngState;
  check(twin->setState(twinState, xHoled.data()) &&
          warmStartsClean(*twin, donor, liveMap, xHoled.data()),
        "warm start counts: a recipient of as many trees takes the donor as "
        "it stands");
  WarmStartReport report;
  check(warmStartsAltered(*recipient, donor, liveMap, xHoled.data(), nullptr,
                          &report) &&
          report.mean[0].offered == 6 && report.mean[0].taken == 3 &&
          report.mean[0].singleLeaves == 0 && report.mean[0].moved == 0 &&
          report.mean[0].cutBack == 0 && report.mean[0].sideDropped == 0,
        "warm start counts: a recipient of fewer trees declines the donor "
        "and, forced, takes 3 of its 6");
  SamplerStateData got, all;
  recipient->getState(got);
  twin->getState(all);
  bool sides = got.chains[0].rngState == all.chains[0].rngState;
  size_t right = 0;
  for (const std::vector<FlatNode>& kept : got.chains[0].forests[0].trees) {
    bool found = false;
    for (const std::vector<FlatNode>& drawn : all.chains[0].forests[0].trees)
      found = found || sameFlatTrees({kept}, {drawn});
    sides = sides && found && kept.size() == 3;
    if ((kept[0].flags & flatMissingGoesRight) != 0) ++right;
  }
  check(sides && right > 0 && right < 3,
        "warm start counts: the sides are drawn once over the donor sample, "
        "the generator past one pass of coins and the trees kept holding "
        "that pass's sides");

  // three trees into five: clean, the rest single leaves at 0
  auto three = make(x, 6, 3);
  three->run(20, 0, none);
  SamplerStateData small;
  three->getState(small);
  auto five = make(x, 43, 5);
  WarmStartResult result;
  bool altered;
  check(warmVerdictMatchesInstall(*five, small, liveMap, x.data(), result,
                                  altered, nullptr, &report) &&
          result == WarmStartResult::ok && !altered &&
          report.mean[0].offered == 3 && report.mean[0].taken == 3 &&
          report.mean[0].singleLeaves == 2,
        "warm start counts: a recipient of more trees takes the donor as it "
        "stands, the rest single leaves");
  SamplerStateData padded;
  five->getState(padded);
  std::vector<std::vector<FlatNode>> expected(small.chains[0].forests[0].trees);
  expected.resize(5, std::vector<FlatNode>(1));
  check(sameFlatTrees(padded.chains[0].forests[0].trees, expected),
        "warm start counts: the donor's trees in its order, then leaves at 0");

  // the unforced answer is the forced install's over seeded pairs of other
  // counts, other rows and sides to draw
  {
    size_t taken = 0, changed = 0;
    bool agree = true;
    for (std::uint32_t seed = 0; seed < 24; ++seed) {
      bool holes = seed % 3 == 2;
      auto pairDonor = make(holes ? xHoled : x, 500 + seed, 2 + seed % 5);
      pairDonor->run(20, 0, none);
      SamplerStateData state;
      pairDonor->getState(state);
      if (holes) state.missingColumns.assign(p, 0);
      auto pair = make(holes ? xHoled : x, 600 + seed, 2 + (seed / 5) % 5);
      agree = agree &&
        warmVerdictMatchesInstall(*pair, state, liveMap,
                                  (holes ? xHoled : x).data(), result,
                                  altered) &&
        result == WarmStartResult::ok &&
        (2 + seed % 5 <= 2 + (seed / 5) % 5 || altered);
      ++(altered ? changed : taken);
    }
    check(agree && taken >= 4 && changed >= 4,
          "warm start counts: over 24 seeded pairs of other counts the "
          "unforced call declines exactly when the forced install is altered");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: warm start from another tree count\n");
}

// A leaf standardization a host supplies at construction, as a sampler
// re-created from its host's record takes it: the build holds the supplied
// centers and scales, a scale not a number for a column without spread, and a
// gp leaf the supplied lengthscales where its model gives none, and the reader
// returns them; supplied equal to a derivation the sampler is the derived one
// bit for bit, draws included; a sampler over a pre-built store derives its
// own. An install keeps the standardization and puts the slopes in as
// stored; kept gp draws made under another center, scale or lengthscale are
// refused, live fits alone installed.
static void testSuppliedStandardization() {
  std::uint64_t savedRngState = rngState;
  rngState = 1732050u;
  const size_t n = 140, p = 2, numChains = 2, numKept = 3;
  std::vector<double> x(n * p), xOther(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 2.0 * runif01() - 1.0;
    xOther[i] = x[i];
    xOther[i + n] = 4.0 * runif01() + 3.0;
    y[i] = std::sin(4.0 * x[i]) + x[i + n] + 0.2 * (runif01() - 0.5);
  }
  std::vector<ext_rng*> rngs;
  auto newRngs = [&](std::uint32_t seed) {
    for (size_t c = 0; c < numChains; ++c) {
      rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
      ext_rng_setSeed(rngs.back(), seed + static_cast<std::uint32_t>(c));
    }
    return rngs.data() + rngs.size() - numChains;
  };
  const size_t covariates[] = {1};
  const double fixedLengthscale[] = {0.5};
  auto optionsFor = [&](bool gp, bool keepTrees,
                        const SuppliedStandardization& supplied,
                        const double* lengthscales = nullptr) {
    SamplerOptions options;
    options.numTrees = 6;
    options.numChains = numChains;
    options.keepTrees = keepTrees;
    options.numSamplesToStore = keepTrees ? numKept : 0;
    options.leafCovariateColumns = covariates;
    options.numLeafCovariates = 1;
    options.gpLeaves = gp;
    options.gpLengthscales = lengthscales;
    options.suppliedStandardization = supplied;
    return options;
  };
  auto make = [&](const std::vector<double>& xs, const SamplerOptions& options,
                  std::uint32_t seed) {
    return createSampler(xs.data(), y.data(), n, p, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, newRngs(seed));
  };
  auto draws = [&](SamplerBase& sampler) {
    std::vector<double> fits(n * numChains * 4);
    Results results;
    results.trainingFits = fits.data();
    sampler.run(2, 4, results);
    return fits;
  };
  auto bits = [](const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() &&
           (a.empty() ||
            std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  };
  Results none;

  for (bool gp : {false, true}) {
    const char* kind = gp ? "gp" : "linear";
    std::unique_ptr<SamplerBase> derived = make(x, optionsFor(gp, true, {}), 10);
    DerivedRecords records;
    derived->derivedRecords(records);
    SamplerStateData derivedState;
    derived->getState(derivedState);
    const ForestStateData& blocks = derivedState.chains[0].forests[0];
    check(bits(records.leafCenters, blocks.leafCovariateCenters) &&
            bits(records.leafScales, blocks.leafCovariateScales) &&
            (gp ? bits(records.leafLengthscales, blocks.leafLengthscales)
                : records.leafLengthscales.empty()),
          "supplied standardization: the reader gives the standardization a "
          "build holds");

    // the derived values handed back: the same sampler, draws included
    SuppliedStandardization same{
      records.leafCenters.data(), records.leafScales.data(),
      gp ? records.leafLengthscales.data() : nullptr};
    std::unique_ptr<SamplerBase> again = make(x, optionsFor(gp, true, same), 10);
    derived = make(x, optionsFor(gp, true, {}), 10);
    SamplerStateData againState;
    again->getState(againState);
    derived->getState(derivedState);
    check(statesAgree(againState, derivedState) &&
            draws(*again) == draws(*derived),
          gp ? "supplied standardization, gp: equal to a derivation, the "
               "sampler and its draws are the derived ones"
             : "supplied standardization, linear: equal to a derivation, the "
               "sampler and its draws are the derived ones");

    // other values, held whatever the rows say; a column without spread
    // divides by 1
    const double centers[] = {0.25}, scales[] = {2.5},
      lengthscales[] = {0.75};
    std::unique_ptr<SamplerBase> other = make(
      xOther, optionsFor(gp, false, {centers, scales, lengthscales}), 20);
    DerivedRecords held;
    other->derivedRecords(held);
    check(held.leafCenters == std::vector<double>{0.25} &&
            held.leafScales == std::vector<double>{2.5} &&
            (gp ? held.leafLengthscales == std::vector<double>{0.75}
                : held.leafLengthscales.empty()),
          "supplied standardization: other values are held");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double flatScales[] = {nan};
    std::vector<double> xFlat(x);
    for (size_t i = 0; i < n; ++i) xFlat[i + n] = 0.25;
    std::unique_ptr<SamplerBase> flat =
      make(xFlat, optionsFor(gp, false, {centers, flatScales, nullptr}), 30);
    flat->derivedRecords(held);
    check(held.leafCenters == std::vector<double>{0.25} &&
            held.leafScales.size() == 1 && std::isnan(held.leafScales[0]),
          "supplied standardization: a scale without spread divides by 1 and "
          "is recorded so");
    // a column the record does not give, its center not a number, derives;
    // a lengthscale not a number derives over the supplied constants
    const double unknown[] = {nan};
    std::unique_ptr<SamplerBase> underived =
      make(x, optionsFor(gp, false, {unknown, unknown, unknown}), 50);
    underived->derivedRecords(held);
    bool derivesColumn = bits(held.leafCenters, records.leafCenters) &&
      bits(held.leafScales, records.leafScales) &&
      bits(held.leafLengthscales, records.leafLengthscales);
    DerivedRecords overSupplied;
    make(x, optionsFor(gp, false, {centers, scales, nullptr}), 60)
      ->derivedRecords(overSupplied);
    make(x, optionsFor(gp, false, {centers, scales, unknown}), 60)
      ->derivedRecords(held);
    check(derivesColumn &&
            bits(held.leafLengthscales, overSupplied.leafLengthscales),
          "supplied standardization: a value not a number derives its column");
    (void) kind;
  }

  // a lengthscale the model supplies is model: a record's is not read, and
  // the reader returns none
  {
    const double centers[] = {0.25}, scales[] = {2.5}, lengthscales[] = {3.0};
    std::unique_ptr<SamplerBase> fixed = make(
      x, optionsFor(true, false, {centers, scales, lengthscales},
                    fixedLengthscale), 40);
    DerivedRecords held;
    fixed->derivedRecords(held);
    SamplerStateData state;
    fixed->getState(state);
    check(held.leafLengthscales.empty() &&
            state.chains[0].forests[0].leafLengthscales ==
              std::vector<double>{0.5} &&
            held.leafCenters == std::vector<double>{0.25},
          "supplied standardization: a model's lengthscale is the one in "
          "force and no record's");
  }

  // a sampler over a pre-built store derives its constants over its own
  // rows, whatever the options supply: here every row, a full data's
  {
    ColumnStore parent;
    check(parent.build(x.data(), n, p, 100, false, nullptr, covariates, 1),
          "supplied standardization: the parent store builds");
    std::vector<size_t> rows(n);
    for (size_t i = 0; i < n; ++i) rows[i] = i;
    ColumnStore store;
    store.buildFromParent(parent, rows.data(), n, nullptr, 0, covariates, 1);
    const double centers[] = {0.25}, scales[] = {2.5};
    std::unique_ptr<SamplerBase> view = createSamplerOverStore(
      std::move(store), y.data(), nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542,
      optionsFor(false, false, {centers, scales, nullptr}), newRngs(50));
    std::unique_ptr<SamplerBase> derived =
      make(x, optionsFor(false, false, {}), 51);
    DerivedRecords fromView, fromData;
    view->derivedRecords(fromView);
    derived->derivedRecords(fromData);
    check(fromView.leafCenters == fromData.leafCenters &&
            fromView.leafScales == fromData.leafScales,
          "supplied standardization: a view derives its own");
  }

  // an install keeps the sampler's standardization and puts slopes in as
  // stored: the kept draws, which no install moves, read the same numbers
  {
    std::unique_ptr<SamplerBase> source = make(x, optionsFor(false, true, {}),
                                               60);
    source->run(30, numKept, none);
    SamplerStateData stored;
    source->getState(stored);
    std::unique_ptr<SamplerBase> target =
      make(xOther, optionsFor(false, true, {}), 61);
    SamplerStateData before, after;
    target->getState(before);
    // column 1's grid is another, so splits on it move by value
    bool installed = target->setState(stored, xOther.data());
    target->getState(after);
    bool slopesAsStored = installed;
    for (size_t c = 0; c < numChains; ++c) {
      const ForestStateData& got = after.chains[c].forests[0];
      const ForestStateData& was = before.chains[c].forests[0];
      slopesAsStored = slopesAsStored &&
        got.leafCovariateCenters == was.leafCovariateCenters &&
        got.leafCovariateScales == was.leafCovariateScales &&
        got.leafCovariateCenters !=
          stored.chains[c].forests[0].leafCovariateCenters &&
        got.savedTreeParams == stored.chains[c].forests[0].savedTreeParams;
    }
    check(slopesAsStored,
          "supplied standardization: an install keeps the standardization "
          "and puts the slopes in as stored");
  }

  // kept gp draws made under another standardization or lengthscale are
  // refused, nothing touched; live fits alone install
  {
    std::unique_ptr<SamplerBase> source = make(x, optionsFor(true, true, {}),
                                               70);
    source->run(10, numKept, none);
    SamplerStateData kept, live;
    source->getState(kept);
    std::unique_ptr<SamplerBase> unkept = make(x, optionsFor(true, false, {}),
                                               71);
    unkept->run(10, 0, none);
    unkept->getState(live);
    DerivedRecords records;
    source->derivedRecords(records);
    std::vector<double> longer(records.leafLengthscales);
    longer[0] *= 2.0;
    std::vector<double> shifted(records.leafCenters);
    shifted[0] += 0.125;
    SuppliedStandardization otherLengthscale{
      records.leafCenters.data(), records.leafScales.data(), longer.data()};
    SuppliedStandardization otherCenter{
      shifted.data(), records.leafScales.data(),
      records.leafLengthscales.data()};
    for (const SuppliedStandardization& supplied :
         {otherLengthscale, otherCenter}) {
      std::unique_ptr<SamplerBase> target =
        make(x, optionsFor(true, true, supplied), 72);
      SamplerStateData before, after;
      target->getState(before);
      bool refused = true;
      check(!target->setState(kept, x.data(), &refused) && refused,
            "supplied standardization, gp: kept draws under another "
            "standardization or lengthscale are refused");
      target->getState(after);
      check(statesAgree(before, after),
            "supplied standardization, gp: a refusal touches nothing");
      check(target->setState(live, x.data(), &refused) && !refused,
            "supplied standardization, gp: live fits alone install");
    }
    std::unique_ptr<SamplerBase> twin = make(x, optionsFor(true, true, {}), 73);
    bool refused = true;
    check(twin->setState(kept, x.data(), &refused) && !refused,
          "supplied standardization, gp: kept draws under the sampler's own "
          "install");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: supplied leaf standardization\n");
}

/// A creation build handed cut grids (SamplerOptions::suppliedCutGrid) holds
/// them in place of its derivations: a numeric column's points and weights as
/// given, a column given none derived, and an ordered factor's grid its
/// levels' whatever its entry. The grids a build derives, handed to another,
/// give that build's codes and draws through the single-forest, combining and
/// multinomial constructors; another grid gives the codes setCutPoints to it
/// gives; derivedRecords reads back what the build holds.
static void testSuppliedCutGrid() {
  std::uint64_t savedRngState = rngState;
  rngState = 1414213u;
  const size_t n = 150, p = 3, K = 3;
  std::vector<double> x(n * p), y(n), z(n);
  std::vector<int> counts(n * K, 0), trials(n, 1);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = static_cast<double>((i % 6) * (i % 6));  // a weighted grid
    x[i + 2 * n] = static_cast<double>(i % 4);           // an ordered factor
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    y[i] = std::sin(4.0 * x[i]) + 0.1 * x[i + n] + z[i] +
           0.2 * (runif01() - 0.5);
    counts[(i % K) * n + i] = 1;
  }
  ColumnKind kinds[] = {ColumnKind::numeric, ColumnKind::numeric,
                        ColumnKind::orderedFactor};
  std::vector<ext_rng*> rngs;
  auto newRng = [&](std::uint_least32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(r, seed);
    rngs.push_back(r);
    return r;
  };
  auto optionsWith = [&](const SuppliedCutGrid& supplied) {
    SamplerOptions options;
    options.numTrees = 10;
    options.predictors.columnTypes = kinds;
    options.suppliedCutGrid = supplied;
    return options;
  };
  using Built = std::unique_ptr<Sampler<ConstantGaussianLeaf>>;
  auto single = [&](const SuppliedCutGrid& supplied) {
    ext_rng* r = newRng(2718u);
    return Built(new Sampler<ConstantGaussianLeaf>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, optionsWith(supplied), &r));
  };
  auto combining = [&](const SuppliedCutGrid& supplied) {
    AmplitudeSpec spec;
    spec.mu.numTrees = 8;
    spec.tau.numTrees = 4;
    spec.z = z.data();
    ext_rng* r = newRng(3141u);
    return Built(new Sampler<ConstantGaussianLeaf>(
      x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
      0.37804942330213542, optionsWith(supplied), spec, &r));
  };
  auto multinomial = [&](const SuppliedCutGrid& supplied) {
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = counts.data();
    spec.trials = trials.data();
    spec.forest.numTrees = 5;
    ext_rng* r = newRng(1618u);
    return Built(new Sampler<ConstantGaussianLeaf>(
      x.data(), n, p, optionsWith(supplied), spec, &r));
  };
  auto draws = [](Sampler<ConstantGaussianLeaf>& sampler) {
    std::vector<double> train(n * 5 * K);
    Results results;
    results.trainingFits = train.data();
    sampler.run(5, 5, results);
    return train;
  };

  Built derived = single({});
  DerivedRecords records;
  derived->derivedRecords(records);
  check(records.cutPoints == derived->data().cutPoints &&
          records.cutMass == derived->data().cutMass &&
          derived->data().cutsWeighted(1) && !derived->data().cutsWeighted(0) &&
          records.cutPoints[2] == std::vector<double>({0.5, 1.5, 2.5}),
        "supplied grid: the reader gives the grids and weights a build holds");

  // the derived grids handed back: the same store and the same draws
  SuppliedCutGrid same{&records.cutPoints, &records.cutMass};
  bool sameAll = true;
  for (int kind = 0; kind < 3; ++kind) {
    Built plain = kind == 0 ? single({}) : kind == 1 ? combining({})
                                                     : multinomial({});
    Built given = kind == 0 ? single(same) : kind == 1 ? combining(same)
                                                       : multinomial(same);
    sameAll = sameAll && given->data().cutPoints == plain->data().cutPoints &&
              given->data().cutMass == plain->data().cutMass &&
              storageDigest(given->data()) == storageDigest(plain->data()) &&
              draws(*given) == draws(*plain);
  }
  check(sameAll,
        "supplied grid: a build's own grids give its codes and draws through "
        "all three constructors");

  // another grid on the first column, none on the second, a foreign one on
  // the ordered factor
  std::vector<std::vector<double>> points(records.cutPoints);
  std::vector<std::vector<double>> mass(records.cutMass);
  points[0] = {0.25, 0.5, 0.75};
  mass[0].clear();
  points[1].clear();
  mass[1].clear();
  points[2] = {9.0};
  SuppliedCutGrid other{&points, &mass};
  bool otherAll = true;
  for (int kind = 0; kind < 3; ++kind) {
    Built given = kind == 0 ? single(other) : kind == 1 ? combining(other)
                                                        : multinomial(other);
    Built twin = kind == 0 ? single({}) : kind == 1 ? combining({})
                                                    : multinomial({});
    const double* grids[] = {points[0].data()};
    std::uint32_t numCuts[] = {3u};
    const size_t column = 0;
    twin->setCutPoints(grids, numCuts, &column, 1, x.data());
    otherAll = otherAll && given->data().cutPoints[0] == points[0] &&
               !given->data().cutsWeighted(0) &&
               given->data().cutPoints[1] == records.cutPoints[1] &&
               given->data().cutMass[1] == records.cutMass[1] &&
               given->data().cutPoints[2] == records.cutPoints[2] &&
               storageDigest(given->data()) == storageDigest(twin->data());
  }
  check(otherAll,
        "supplied grid: a numeric column takes the grid given, one given "
        "none derives, and an ordered factor keeps its levels'");

  // weights given with a grid are the store's
  std::vector<std::vector<double>> weighted(records.cutPoints);
  std::vector<std::vector<double>> weights(records.cutMass);
  weights[0].resize(weighted[0].size() + 1);
  for (size_t k = 0; k < weights[0].size(); ++k)
    weights[0][k] = static_cast<double>(k * k);
  Built heavy = single(SuppliedCutGrid{&weighted, &weights});
  check(heavy->data().cutMass[0] == weights[0] &&
          heavy->data().cutPoints[0] == records.cutPoints[0],
        "supplied grid: weights given with a grid weigh its points");

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a creation build takes the grids it is handed\n");
}

/// A linear or gp leaf derives its covariates' center, scale and default
/// lengthscale over the rows in the likelihood: membership by positive weight
/// and the mask, every row where none is in. Too little of a covariate over
/// those rows refuses a derivation (the stated-lengthscale gp taking the
/// placeholder) and keeps it at a re-derivation.
static void testLeafStandardizationRows() {
  std::uint64_t savedRngState = rngState;
  rngState = 2236067u;
  const size_t n = 1200, p = 2, numChains = 2;
  std::vector<double> x(n * p), y(n), binary(n), count(n), w(n), mask(n);
  for (size_t i = 0; i < n; ++i) {
    bool in = i % 5 != 2;
    x[i] = runif01();
    // the rows out lie far off, so the rows in's numbers are not every row's
    x[i + n] = in ? 2.0 * runif01() - 1.0 : 40.0 + 10.0 * runif01();
    y[i] = std::sin(4.0 * x[i]) + 0.2 * (runif01() - 0.5);
    binary[i] = y[i] > 0.0 ? 1.0 : 0.0;
    count[i] = std::floor(3.0 * x[i]) + static_cast<double>(i % 2);
    w[i] = in ? 1.0 : 0.0;
    mask[i] = w[i];
  }
  // a row in missing it, among those the heuristic's stride takes
  x[n + 3] = std::numeric_limits<double>::quiet_NaN();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<std::uint8_t> rowsIn(n);
  for (size_t i = 0; i < n; ++i) rowsIn[i] = w[i] > 0.0 ? 1 : 0;
  // the moments over the rows in, by hand
  double total = 0.0, squares = 0.0;
  size_t observed = 0;
  for (size_t i = 0; i < n; ++i)
    if (rowsIn[i] != 0 && !std::isnan(x[i + n])) {
      total += x[i + n];
      ++observed;
    }
  double center = total / static_cast<double>(observed);
  for (size_t i = 0; i < n; ++i)
    if (rowsIn[i] != 0 && !std::isnan(x[i + n]))
      squares += (x[i + n] - center) * (x[i + n] - center);
  double scale = std::sqrt(squares / static_cast<double>(observed - 1));
  // the heuristic over the rows in: a stride over those rows, a missing row
  // counted at 0
  std::vector<double> sampled;
  {
    size_t m = 0;
    for (size_t i = 0; i < n; ++i) m += rowsIn[i];
    size_t stride = 1 + (m - 1) / 512;
    for (size_t i = 0, k = 0; i < n; ++i) {
      if (rowsIn[i] == 0) continue;
      if (k++ % stride != 0) continue;
      sampled.push_back(std::isnan(x[i + n]) ? 0.0
                                             : (x[i + n] - center) / scale);
    }
  }
  std::vector<double> distances;
  for (size_t a = 0; a < sampled.size(); ++a)
    for (size_t b = a + 1; b < sampled.size(); ++b)
      distances.push_back(std::fabs(sampled[a] - sampled[b]));
  std::sort(distances.begin(), distances.end());
  size_t numDistances = distances.size();
  double heuristic = numDistances % 2 == 1
    ? distances[numDistances / 2]
    : 0.5 * (distances[numDistances / 2 - 1] + distances[numDistances / 2]);
  auto near = [](double a, double b) {
    return std::fabs(a - b) <= 1e-12 * std::max(1.0, std::fabs(b));
  };

  // the moments: over a row set the literals; a null or all-ones set
  // today's bits; too little answered for one row, one value, no observed
  // value, and every row at one value
  {
    double mean, sd, allMean, allSd, onesMean, onesSd;
    bool tooLittle = standardizationMomentsForColumn(x.data() + n, n, &mean,
                                                     &sd, rowsIn.data());
    check(!tooLittle && near(mean, center) && near(sd, scale),
          "leaf rows: the moments over the rows in");
    std::vector<std::uint8_t> ones(n, 1);
    standardizationMomentsForColumn(x.data() + n, n, &allMean, &allSd);
    standardizationMomentsForColumn(x.data() + n, n, &onesMean, &onesSd,
                                    ones.data());
    check(std::memcmp(&allMean, &onesMean, sizeof(double)) == 0 &&
            std::memcmp(&allSd, &onesSd, sizeof(double)) == 0,
          "leaf rows: an all-ones row set reads today's bits");
    std::vector<std::uint8_t> one(n, 0);
    one[0] = 1;
    std::vector<double> atTwo(x.begin() + n, x.end());
    for (size_t i = 0; i < n; ++i)
      if (rowsIn[i] != 0) atTwo[i] = 2.0;
    std::vector<double> missing(atTwo);
    for (size_t i = 0; i < n; ++i)
      if (rowsIn[i] != 0) missing[i] = nan;
    std::vector<double> flat(n, 0.5);
    check(standardizationMomentsForColumn(x.data() + n, n, &mean, &sd,
                                          one.data()) &&
            standardizationMomentsForColumn(atTwo.data(), n, &mean, &sd,
                                            rowsIn.data()) &&
            mean == 2.0 && sd == 1.0 &&
            standardizationMomentsForColumn(missing.data(), n, &mean, &sd,
                                            rowsIn.data()) &&
            standardizationMomentsForColumn(flat.data(), n, &mean, &sd) &&
            !standardizationMomentsForColumn(atTwo.data(), n, &mean, &sd),
          "leaf rows: one row, one value, none observed and a constant "
          "column hold too little, the rows out at other values not");
  }

  std::vector<ext_rng*> rngs;
  auto newRngs = [&](std::uint32_t seed) {
    for (size_t c = 0; c < numChains; ++c) {
      rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
      ext_rng_setSeed(rngs.back(), seed + static_cast<std::uint32_t>(c));
    }
    return rngs.data() + rngs.size() - numChains;
  };
  rngs.reserve(64);
  const size_t covariates[] = {1};
  const size_t bothCovariates[] = {0, 1};
  auto optionsFor = [&](bool gp, const double* active = nullptr,
                        SuppliedStandardization supplied = {},
                        bool both = false) {
    SamplerOptions options;
    options.numTrees = 6;
    options.numChains = numChains;
    options.leafCovariateColumns = both ? bothCovariates : covariates;
    options.numLeafCovariates = both ? 2 : 1;
    options.gpLeaves = gp;
    options.gpMaxLeafSize = 60;
    options.activeRows = active;
    options.suppliedStandardization = supplied;
    options.shape = 3.0;
    return options;
  };
  auto make = [&](const std::vector<double>& xs, const double* ys,
                  const double* weights, ResponseFamily family,
                  const SamplerOptions& options, std::uint32_t seed) {
    return createSampler(xs.data(), ys, n, p, weights, nullptr, family, 1.0,
                         3.0, 0.37804942330213542, options, newRngs(seed));
  };
  auto draws = [&](SamplerBase& sampler, size_t sweeps) {
    std::vector<double> fits(n * numChains * sweeps);
    Results results;
    results.trainingFits = fits.data();
    sampler.run(0, sweeps, results);
    return fits;
  };
  auto bits = [](const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() &&
           (a.empty() ||
            std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  };

  // at creation, by weight under gaussian and by the creation mask under
  // probit and nbinom: the literals, the gp heuristic among them, and the
  // draws of a sampler handed the same numbers as a record, bit for bit
  struct Family {
    const char* name;
    ResponseFamily family;
    const double* y;
    bool byMask;
  };
  const Family families[] = {
    {"gaussian", ResponseFamily::gaussian, y.data(), false},
    {"probit", ResponseFamily::probit, binary.data(), true},
    {"nbinom", ResponseFamily::nbinom, count.data(), true}};
  for (const Family& f : families) {
    for (bool gp : {false, true}) {
      std::string name = std::string("leaf rows, ") + f.name +
                         (gp ? ", gp" : ", linear");
      const double* weights = f.byMask ? nullptr : w.data();
      std::unique_ptr<SamplerBase> derived =
        make(x, f.y, weights, f.family,
             optionsFor(gp, f.byMask ? mask.data() : nullptr), 10);
      DerivedRecords records;
      derived->derivedRecords(records);
      check(records.leafCenters.size() == 1 &&
              near(records.leafCenters[0], center) &&
              near(records.leafScales[0], scale) &&
              (!gp || near(records.leafLengthscales[0], heuristic)) &&
              records.leafFromAllRows == std::vector<std::uint8_t>{0},
            (name + ": the numbers over the rows in").c_str());
      SuppliedStandardization same{
        records.leafCenters.data(), records.leafScales.data(),
        gp ? records.leafLengthscales.data() : nullptr};
      std::unique_ptr<SamplerBase> supplied =
        make(x, f.y, weights, f.family, optionsFor(gp, nullptr, same), 10);
      if (f.byMask) {
        check(derived->setActiveRows(mask.data()) &&
                supplied->setActiveRows(mask.data()),
              (name + ": the mask installs").c_str());
      }
      check(bits(draws(*derived, 20), draws(*supplied, 20)),
            (name + ": draws a record's twin draws, for 20 sweeps").c_str());
    }
  }

  // no row in: every row read, and said; too little refused, naming the
  // column, the gp leaf with stated lengthscales taking the placeholder
  {
    std::vector<double> none(n, 0.0);
    std::unique_ptr<SamplerBase> sampler =
      make(x, y.data(), none.data(), ResponseFamily::gaussian,
           optionsFor(false), 20);
    DerivedRecords records;
    sampler->derivedRecords(records);
    double mean, sd;
    standardizationMomentsForColumn(x.data() + n, n, &mean, &sd);
    check(records.leafCenters[0] == mean && records.leafScales[0] == sd &&
            records.leafFromAllRows == std::vector<std::uint8_t>{1},
          "leaf rows: no row in reads every row and says so");

    std::vector<double> atTwo(x);
    for (size_t i = 0; i < n; ++i)
      if (rowsIn[i] != 0) atTwo[i + n] = 2.0;
    for (bool gp : {false, true}) {
      bool refused = false;
      try {
        make(atTwo, y.data(), w.data(), ResponseFamily::gaussian,
             optionsFor(gp), 21);
      } catch (const LeafStandardizationError& error) {
        refused = error.covariate == 0 && error.column == 1 &&
                  error.numRowsIn == 960 && error.numObservations == n &&
                  error.observed;
      }
      check(refused, gp ? "leaf rows: gp refused, rows in at one value"
                        : "leaf rows: linear refused, rows in at one value");
    }
    const double stated[] = {0.8};
    SamplerOptions options = optionsFor(true);
    options.gpLengthscales = stated;
    std::unique_ptr<SamplerBase> placeholder =
      make(atTwo, y.data(), w.data(), ResponseFamily::gaussian, options, 22);
    placeholder->derivedRecords(records);
    check(records.leafCenters[0] == 2.0 && records.leafScales[0] == 1.0 &&
            records.leafLengthscales.empty(),
          "leaf rows: a stated gp lengthscale takes the placeholder");
  }

  // a re-derivation over the rows in now, the weights and the mask in force:
  // a covariate with too little kept and reported, the other derived, and a
  // covariate not marked held
  {
    std::unique_ptr<SamplerBase> sampler =
      make(x, y.data(), nullptr, ResponseFamily::gaussian,
           optionsFor(false, nullptr, {}, true), 30);
    DerivedRecords before, after;
    sampler->derivedRecords(before);
    sampler->setWeights(w.data());
    const std::uint8_t onlyCovariate[] = {0, 1};
    std::vector<std::uint8_t> kept =
      sampler->rederiveStandardization(onlyCovariate);
    sampler->derivedRecords(after);
    check(kept == std::vector<std::uint8_t>{0, 0} &&
            after.leafCenters[0] == before.leafCenters[0] &&
            after.leafScales[0] == before.leafScales[0] &&
            near(after.leafCenters[1], center) &&
            near(after.leafScales[1], scale),
          "leaf rows: the marked covariate over the rows in, another held");
    // the mask in force joins the weights: the rows both keep
    std::vector<double> mask2(n);
    std::vector<std::uint8_t> both(n);
    for (size_t i = 0; i < n; ++i) {
      mask2[i] = i % 2 == 0 ? 1.0 : 0.0;
      both[i] = rowsIn[i] != 0 && mask2[i] != 0.0 ? 1 : 0;
    }
    check(sampler->setActiveRows(mask2.data()),
          "leaf rows: a mask installs on a gaussian sampler");
    // a mask the chains refuse leaves the copy the one they hold
    std::vector<double> refusedMask(n, 1.0);
    refusedMask[1] = 0.5;
    check(!sampler->setActiveRows(refusedMask.data()),
          "leaf rows: a mask off 0 and 1 is refused");
    sampler->rederiveStandardization(onlyCovariate);
    sampler->derivedRecords(after);
    double mean, sd;
    standardizationMomentsForColumn(x.data() + n, n, &mean, &sd, both.data());
    check(after.leafCenters[1] == mean && after.leafScales[1] == sd,
          "leaf rows: a re-derivation reads the mask in force");
    // a covariate at one value over the rows in is kept and reported
    std::vector<double> atTwo(x.begin() + n, x.end());
    for (size_t i = 0; i < n; ++i)
      if (both[i] != 0) atTwo[i] = 2.0;
    const size_t column[] = {1};
    check(sampler->updatePredictor(atTwo.data(), column, 1, true, false) ==
            PredictorUpdateResult::accepted,
          "leaf rows: the column update installs");
    DerivedRecords held;
    sampler->derivedRecords(held);
    kept = sampler->rederiveStandardization();
    sampler->derivedRecords(after);
    check(kept == std::vector<std::uint8_t>{0, 1} &&
            after.leafCenters[1] == held.leafCenters[1] &&
            after.leafScales[1] == held.leafScales[1],
          "leaf rows: too little kept at a re-derivation, and reported");
    // no row in: every covariate kept
    std::vector<double> none(n, 0.0);
    sampler->setWeights(none.data());
    kept = sampler->rederiveStandardization();
    check(kept == std::vector<std::uint8_t>{1, 1},
          "leaf rows: no row in keeps every covariate");
  }

  // a whole-data replacement derives over its own rows in, the mask it is
  // handed with them; a refused one leaves the sampler as an untouched
  // twin's, bit for bit
  {
    std::unique_ptr<SamplerBase> sampler =
      make(x, binary.data(), nullptr, ResponseFamily::probit,
           optionsFor(false, mask.data()), 40);
    check(sampler->setActiveRows(mask.data()),
          "leaf rows: the creation mask installs");
    const size_t m = n / 2;
    std::vector<double> xHalf(m * p), yHalf(m), maskHalf(m);
    std::vector<std::uint8_t> halfIn(m);
    for (size_t i = 0; i < m; ++i) {
      xHalf[i] = x[i];
      xHalf[i + m] = x[i + n] + 1.0;
      yHalf[i] = binary[i];
      maskHalf[i] = i % 4 == 0 ? 0.0 : 1.0;
      halfIn[i] = maskHalf[i] != 0.0 ? 1 : 0;
    }
    DataReplacement replacement;
    replacement.activeRows = maskHalf.data();
    check(sampler->setData(xHalf.data(), yHalf.data(), m, nullptr, nullptr,
                           nullptr, 0, nullptr, replacement),
          "leaf rows: a replacement of another row count installs");
    DerivedRecords records;
    sampler->derivedRecords(records);
    double mean, sd;
    standardizationMomentsForColumn(xHalf.data() + m, m, &mean, &sd,
                                    halfIn.data());
    check(records.leafCenters[0] == mean && records.leafScales[0] == sd,
          "leaf rows: a replacement derives over the rows its mask keeps");
    // the replacement dropped the mask, so a re-derivation reads its rows
    // and no row past them
    sampler->rederiveStandardization();
    sampler->derivedRecords(records);
    standardizationMomentsForColumn(xHalf.data() + m, m, &mean, &sd);
    check(records.leafCenters[0] == mean && records.leafScales[0] == sd,
          "leaf rows: the mask goes with the replacement");

    std::unique_ptr<SamplerBase> twin =
      make(x, binary.data(), nullptr, ResponseFamily::probit,
           optionsFor(false, mask.data()), 41);
    std::unique_ptr<SamplerBase> refused =
      make(x, binary.data(), nullptr, ResponseFamily::probit,
           optionsFor(false, mask.data()), 41);
    twin->setActiveRows(mask.data());
    refused->setActiveRows(mask.data());
    std::vector<double> xFlat(xHalf);
    for (size_t i = 0; i < m; ++i)
      if (halfIn[i] != 0) xFlat[i + m] = 3.0;
    bool thrown = false;
    try {
      (void) refused->setData(xFlat.data(), yHalf.data(), m, nullptr, nullptr,
                              nullptr, 0, nullptr, replacement);
    } catch (const LeafStandardizationError& error) {
      thrown = error.numRowsIn == 3 * m / 4 && error.numObservations == m;
    }
    SamplerStateData twinState, refusedState;
    twin->getState(twinState);
    refused->getState(refusedState);
    check(thrown && refused->data().numObservations == n &&
            bits(draws(*twin, 5), draws(*refused, 5)),
          "leaf rows: a refused replacement leaves the sampler untouched");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a leaf derives over the rows in the likelihood\n");
}

/// The setters of rows re-derive every leaf covariate through
/// rederiveStandardization with no column marked: after setWeights or
/// setActiveRows that call leaves the sampler bit for bit where a twin's
/// predictor update of every leaf column to the values in force, followed by
/// a re-derivation of those columns, leaves it (a gp leaf's fits, which the
/// update rebuilds in another order, rebuilt on both). A covariate the rows
/// in hold at one value is kept and reported, its numbers unmoved, while the
/// other re-derives.
static void testRederiveStandardization() {
  std::uint64_t savedRngState = rngState;
  rngState = 1732050u;
  const size_t n = 300, p = 3, numChains = 2;
  std::vector<double> x(n * p), y(n), binary(n), w(n), mask(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    // the rows out by weight lie far off
    x[i + 2 * n] = i % 4 == 0 ? 10.0 + runif01() : 2.0 * runif01() - 1.0;
    y[i] = std::sin(4.0 * x[i]) + 0.3 * x[i + 2 * n] +
           0.2 * (runif01() - 0.5);
    binary[i] = y[i] > 0.5 ? 1.0 : 0.0;
    w[i] = i % 4 == 0 ? 0.0 : 1.0;
    mask[i] = i % 3 == 0 ? 0.0 : 1.0;
  }
  // the leaf columns' values in force, column-major, as a predictor update
  // takes them
  const size_t leafColumns[] = {0, 2};
  const std::uint8_t leafMarks[] = {1, 0, 1};
  std::vector<double> leafValues(x.begin(), x.begin() + n);
  leafValues.insert(leafValues.end(), x.begin() + 2 * n, x.end());

  std::vector<ext_rng*> rngs;
  auto make = [&](bool gp, bool byMask, std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = 6;
    options.numChains = numChains;
    options.leafCovariateColumns = leafColumns;
    options.numLeafCovariates = 2;
    options.gpLeaves = gp;
    options.gpMaxLeafSize = 60;
    // a linear leaf's saved slopes are converted with the live ones
    options.keepTrees = !gp;
    options.numSamplesToStore = gp ? 0 : 5;
    ext_rng* chainRngs[numChains];
    for (size_t c = 0; c < numChains; ++c) {
      rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
      ext_rng_setSeed(rngs.back(), seed + static_cast<std::uint32_t>(c));
      chainRngs[c] = rngs.back();
    }
    std::unique_ptr<SamplerBase> sampler = createSampler(
      x.data(), byMask ? binary.data() : y.data(), n, p, nullptr, nullptr,
      byMask ? ResponseFamily::probit : ResponseFamily::gaussian, 1.0, 3.0,
      0.37804942330213542, options, chainRngs);
    Results results;
    sampler->run(5, gp ? 0 : 5, results);
    return sampler;
  };
  auto draws = [&](SamplerBase& sampler) {
    std::vector<double> fits(n * numChains * 5);
    Results results;
    results.trainingFits = fits.data();
    sampler.run(0, 5, results);
    return fits;
  };
  auto same = [](const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() &&
           (a.empty() ||
            std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  };
  auto sameRecords = [&](const DerivedRecords& a, const DerivedRecords& b) {
    return same(a.leafCenters, b.leafCenters) &&
           same(a.leafScales, b.leafScales) &&
           same(a.leafLengthscales, b.leafLengthscales);
  };
  auto setRows = [&](SamplerBase& sampler, bool byMask) {
    if (byMask) return sampler.setActiveRows(mask.data());
    sampler.setWeights(w.data());
    return true;
  };
  auto refresh = [&](SamplerBase& sampler) {
    return sampler.updatePredictor(leafValues.data(), leafColumns, 2, true,
                                   false) == PredictorUpdateResult::accepted;
  };

  for (bool byMask : {false, true}) {
    for (bool gp : {false, true}) {
      std::string name = std::string("row setters, ") +
                         (byMask ? "probit mask" : "gaussian weights") +
                         (gp ? ", gp" : ", linear");
      std::unique_ptr<SamplerBase> sampler = make(gp, byMask, 40);
      std::unique_ptr<SamplerBase> twin = make(gp, byMask, 40);
      DerivedRecords before, after, twinAfter;
      sampler->derivedRecords(before);
      bool installed = (!gp || refresh(*sampler)) && setRows(*sampler, byMask);
      std::vector<std::uint8_t> kept = sampler->rederiveStandardization();
      installed = installed && setRows(*twin, byMask) && refresh(*twin);
      std::vector<std::uint8_t> twinKept =
        twin->rederiveStandardization(leafMarks);
      sampler->derivedRecords(after);
      twin->derivedRecords(twinAfter);
      SamplerStateData state, twinState;
      sampler->getState(state);
      twin->getState(twinState);
      check(installed && kept == std::vector<std::uint8_t>{0, 0} &&
              twinKept == kept && sameRecords(after, twinAfter) &&
              !same(after.leafCenters, before.leafCenters) &&
              statesAgree(state, twinState) &&
              same(draws(*sampler), draws(*twin)),
            (name + ": every covariate as a twin's update and "
                    "re-derivation of the leaf columns")
              .c_str());
    }
  }

  // a keep: the rows in hold the first covariate at one value, so it keeps
  // its numbers while the second re-derives over those rows
  {
    std::unique_ptr<SamplerBase> sampler = make(false, false, 50);
    std::vector<double> atHalf(x.begin(), x.begin() + n);
    for (size_t i = 0; i < n; ++i)
      if (w[i] > 0.0) atHalf[i] = 0.5;
    const size_t column[] = {0};
    bool installed = sampler->updatePredictor(atHalf.data(), column, 1, true,
                                              false) ==
                     PredictorUpdateResult::accepted;
    DerivedRecords held, after;
    sampler->derivedRecords(held);
    sampler->setWeights(w.data());
    std::vector<std::uint8_t> kept = sampler->rederiveStandardization();
    sampler->derivedRecords(after);
    std::vector<std::uint8_t> rowsIn(n);
    for (size_t i = 0; i < n; ++i) rowsIn[i] = w[i] > 0.0 ? 1 : 0;
    double mean, sd;
    standardizationMomentsForColumn(x.data() + 2 * n, n, &mean, &sd,
                                    rowsIn.data());
    check(installed && kept == std::vector<std::uint8_t>{1, 0} &&
            after.leafCenters[0] == held.leafCenters[0] &&
            after.leafScales[0] == held.leafScales[0] &&
            after.leafCenters[1] == mean && after.leafScales[1] == sd,
          "row setters: a covariate the rows in hold at one value kept and "
          "reported, the other derived over them");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: the setters of rows re-derive every leaf covariate\n");
}

/// A whole-data replacement keeps, takes or derives the grid, the leaf
/// standardization and the transform as it is told (DataReplacement), and
/// rederiveStandardization takes a linear leaf's standardization over the
/// values in force, its coefficients rewritten into it.
static void testDataReplacement() {
  std::uint64_t savedRngState = rngState;
  rngState = 1618033u;
  const size_t n = 120, p = 2, more = 40, numKept = 3, nProbe = 30;
  std::vector<double> x(n * p), y(n), xEdit(n * p), yEdit(n),
    xMore((n + more) * p), yMore(n + more), xProbe(nProbe * p);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 2.0 * runif01() - 1.0;
    y[i] = std::sin(4.0 * x[i]) + x[i + n] + 0.2 * (runif01() - 0.5);
    xEdit[i] = 2.0 * x[i];
    xEdit[i + n] = 3.0 * x[i + n] + 5.0;
    yEdit[i] = 3.0 * y[i] + 10.0;
  }
  for (size_t i = 0; i < n + more; ++i) {
    xMore[i] = runif01();
    xMore[i + n + more] = 2.0 * runif01() - 1.0;
    yMore[i] = std::sin(4.0 * xMore[i]) + xMore[i + n + more];
  }
  for (size_t i = 0; i < nProbe * p; ++i) xProbe[i] = runif01();
  std::vector<ext_rng*> rngs;
  const size_t covariates[] = {1};
  auto make = [&](std::uint32_t seed, bool linear = true) {
    SamplerOptions options;
    options.numTrees = 6;
    options.keepTrees = true;
    options.numSamplesToStore = numKept;
    if (linear) {
      options.leafCovariateColumns = covariates;
      options.numLeafCovariates = 1;
    }
    rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
    ext_rng_setSeed(rngs.back(), seed);
    std::unique_ptr<SamplerBase> sampler = createSampler(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &rngs.back());
    Results none;
    sampler->run(10, numKept, none);
    return sampler;
  };
  auto kept = [&](SamplerBase& sampler) {
    std::vector<double> out(nProbe * numKept);
    sampler.predict(xProbe.data(), nProbe, 1, out.data());
    return out;
  };
  auto worstGap = [](const std::vector<double>& a,
                     const std::vector<double>& b) {
    double worst = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
      worst = std::max(worst, std::fabs(a[i] - b[i]) /
                                std::max(1.0, std::fabs(a[i])));
    return worst;
  };
  auto standardizationOf = [](const DerivedRecords& r) {
    return std::make_pair(r.leafCenters, r.leafScales);
  };
  using Choice = DerivedValueChoice;

  // keep: every value held, the new values quantized on the grid in force,
  // and nothing rewritten
  {
    std::unique_ptr<SamplerBase> sampler = make(4u);
    DerivedRecords before, after;
    sampler->derivedRecords(before);
    std::vector<double> keptBefore = kept(*sampler);
    double sigmaBefore = sampler->sigma(0);
    DataReplacement keep;
    keep.grid = keep.standardization = keep.range = Choice::keep;
    check(sampler->setData(xEdit.data(), yEdit.data(), n, nullptr, nullptr,
                           nullptr, 0, nullptr, keep),
          "data replacement, keep: the replacement is taken");
    sampler->derivedRecords(after);
    bool codesOnGrid = true;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < p; ++j)
        codesOnGrid = codesOnGrid && sampler->data().codeAt(j, i) ==
                                       sampler->data().codeFor(j, xEdit[i + j * n]);
    check(after.cutPoints == before.cutPoints &&
            after.cutMass == before.cutMass &&
            after.rangeMin == before.rangeMin &&
            after.rangeMax == before.rangeMax &&
            standardizationOf(after) == standardizationOf(before) &&
            sampler->sigma(0) == sigmaBefore && codesOnGrid,
          "data replacement, keep: the grid, the range, the standardization "
          "and sigma are held, the new values quantized on the grid");
    check(kept(*sampler) == keptBefore,
          "data replacement, keep: the kept draws are not rewritten");
    // nor sigma, however many times the data are replaced under it
    bool sigmaHeld = true;
    for (int k = 0; k < 50 && sigmaHeld; ++k)
      sigmaHeld = sampler->setData(k % 2 == 0 ? x.data() : xEdit.data(),
                                   k % 2 == 0 ? y.data() : yEdit.data(), n,
                                   nullptr, nullptr, nullptr, 0, nullptr,
                                   keep) &&
                  sampler->sigma(0) == sigmaBefore;
    check(sigmaHeld, "data replacement, keep: sigma is held bit for bit");

    // more rows: the covariates gathered for them under the constants in
    // force
    check(sampler->setData(xMore.data(), yMore.data(), n + more, nullptr,
                           nullptr, nullptr, 0, nullptr, keep),
          "data replacement, keep: a longer replacement is taken");
    sampler->derivedRecords(after);
    std::vector<double> fits((n + more) * 2);
    Results results;
    results.trainingFits = fits.data();
    sampler->run(0, 2, results);
    bool finite = true;
    for (double fit : fits) finite = finite && std::isfinite(fit);
    check(after.cutPoints == before.cutPoints &&
            standardizationOf(after) == standardizationOf(before) &&
            sampler->data().numObservations == n + more && finite,
          "data replacement, keep: 40 rows more are read on the grid and "
          "the standardization in force, and the fits are finite");
  }

  // supply: the host's grid, range and standardization, the splits on a
  // changed grid moved by value and the kept draws rewritten into the rest
  {
    std::unique_ptr<SamplerBase> sampler = make(5u);
    DerivedRecords before, after;
    sampler->derivedRecords(before);
    std::vector<double> keptBefore = kept(*sampler);
    std::vector<std::vector<double>> points(before.cutPoints),
      mass(before.cutMass);
    points[0] = {0.25, 0.5, 0.75};
    mass[0].clear();
    const double centers[] = {0.5}, scales[] = {2.0};
    DataReplacement supply;
    supply.grid = supply.standardization = supply.range = Choice::supply;
    supply.suppliedGrid = SuppliedCutGrid{&points, &mass};
    supply.suppliedStandardization = SuppliedStandardization{centers, scales};
    supply.rangeMin = before.rangeMin - 1.0;
    supply.rangeMax = before.rangeMax + 2.0;
    check(sampler->setData(x.data(), y.data(), n, nullptr, nullptr, nullptr, 0,
                           nullptr, supply),
          "data replacement, supply: the replacement is taken");
    sampler->derivedRecords(after);
    check(after.cutPoints[0] == points[0] && after.cutMass[0].empty() &&
            after.cutPoints[1] == before.cutPoints[1] &&
            after.rangeMin == supply.rangeMin &&
            after.rangeMax == supply.rangeMax &&
            after.leafCenters == std::vector<double>{0.5} &&
            after.leafScales == std::vector<double>{2.0},
          "data replacement, supply: the grid, range and standardization are "
          "the ones handed");
    check(worstGap(kept(*sampler), keptBefore) < 1e-12,
          "data replacement, supply: the kept draws predict what they did");

    // the values in force, handed back, are a keep
    std::vector<double> keptNow = kept(*sampler);
    double sigmaNow = sampler->sigma(0);
    DerivedRecords now;
    sampler->derivedRecords(now);
    DataReplacement back;
    back.grid = back.standardization = back.range = Choice::supply;
    back.suppliedGrid = SuppliedCutGrid{&now.cutPoints, &now.cutMass};
    back.suppliedStandardization = SuppliedStandardization{
      now.leafCenters.data(), now.leafScales.data()};
    back.rangeMin = now.rangeMin;
    back.rangeMax = now.rangeMax;
    check(sampler->setData(x.data(), y.data(), n, nullptr, nullptr, nullptr, 0,
                           nullptr, back),
          "data replacement, supply: the values in force are taken");
    sampler->derivedRecords(after);
    check(after.cutPoints == now.cutPoints && after.cutMass == now.cutMass &&
            kept(*sampler) == keptNow && sampler->sigma(0) == sigmaNow,
          "data replacement, supply: the values in force move nothing");
  }

  // derive, placed by position: a grid whose count is held keeps every
  // split on its position
  {
    std::unique_ptr<SamplerBase> sampler = make(6u, false);
    SamplerStateData before, after;
    sampler->getState(before);
    std::vector<double> oldGrid = sampler->data().cutPoints[0];
    DataReplacement position;
    position.placement = SplitPlacement::byPosition;
    check(sampler->setData(xEdit.data(), y.data(), n, nullptr, nullptr,
                           nullptr, 0, nullptr, position),
          "data replacement, position: the replacement is taken");
    sampler->getState(after);
    const std::vector<double>& newGrid = sampler->data().cutPoints[0];
    bool positionsHeld = newGrid.size() == oldGrid.size();
    size_t numSplits = 0;
    const std::vector<std::vector<FlatNode>>& treesBefore =
      before.chains[0].forests[0].trees;
    const std::vector<std::vector<FlatNode>>& treesAfter =
      after.chains[0].forests[0].trees;
    for (size_t t = 0; t < treesBefore.size() && positionsHeld; ++t) {
      if (treesBefore[t].size() != treesAfter[t].size()) continue;
      for (size_t k = 0; k < treesBefore[t].size() && positionsHeld; ++k) {
        if (treesBefore[t][k].variable != 0) continue;
        size_t index = static_cast<size_t>(
          std::find(oldGrid.begin(), oldGrid.end(), treesBefore[t][k].value) -
          oldGrid.begin());
        positionsHeld = index < newGrid.size() &&
                        treesAfter[t][k].value == newGrid[index];
        ++numSplits;
      }
    }
    check(positionsHeld && numSplits > 0,
          "data replacement, position: a doubled column's splits keep their "
          "positions on its grid");
  }

  // rederiveStandardization: the derivation over the values in force, the
  // kept draws and the live fit the functions they were
  {
    std::unique_ptr<SamplerBase> sampler = make(7u),
      twin = make(7u);
    std::vector<double> keptBefore = kept(*sampler);
    const size_t column = 1;
    std::vector<double> moved(xEdit.begin() + n, xEdit.end());
    for (SamplerBase* s : {sampler.get(), twin.get()})
      check(s->updatePredictor(moved.data(), &column, 1, true, false) ==
              PredictorUpdateResult::accepted,
            "rederive standardization: the column is replaced");
    DerivedRecords held, derived;
    sampler->derivedRecords(held);
    // a mark on a column the leaf does not take derives nothing
    const std::uint8_t onlyOther[] = {1, 0};
    sampler->rederiveStandardization(onlyOther);
    sampler->derivedRecords(derived);
    check(standardizationOf(derived) == standardizationOf(held),
          "rederive standardization: a covariate not replaced is held");
    sampler->rederiveStandardization();
    sampler->derivedRecords(derived);
    double mean, sd;
    standardizationMomentsForColumn(moved.data(), n, &mean, &sd);
    std::vector<double> liveSampler(n), liveTwin(n);
    sampler->fitsWithoutOffset(0, liveSampler.data());
    twin->fitsWithoutOffset(0, liveTwin.data());
    check(derived.leafCenters == std::vector<double>{mean} &&
            derived.leafScales == std::vector<double>{sd} &&
            held.leafCenters != derived.leafCenters,
          "rederive standardization: the reader gives the derivation over "
          "the values in force");
    check(worstGap(kept(*sampler), keptBefore) < 1e-12 &&
            worstGap(liveSampler, liveTwin) < 1e-12,
          "rederive standardization: the kept draws and the live fit are "
          "the functions they were");
  }

  // keep on an aft sampler, which R refuses the whole-data conduit and a
  // header-only host does not: the transform of the log times is held
  {
    std::vector<double> status(n, 1.0);
    for (size_t i = 0; i < n; i += 3) status[i] = 0.0;
    SamplerOptions options;
    options.numTrees = 6;
    options.keepTrees = true;
    options.numSamplesToStore = numKept;
    options.survivalStatus = status.data();
    rngs.push_back(ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL));
    ext_rng_setSeed(rngs.back(), 8u);
    std::unique_ptr<SamplerBase> sampler = createSampler(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::aft, 1.0,
      3.0, 0.37804942330213542, options, &rngs.back());
    Results none;
    sampler->run(10, numKept, none);
    DerivedRecords before, after;
    sampler->derivedRecords(before);
    std::vector<double> keptBefore = kept(*sampler);
    DataReplacement keep;
    keep.range = Choice::keep;
    check(sampler->setData(x.data(), yEdit.data(), n, nullptr, nullptr,
                           nullptr, 0, nullptr, keep),
          "data replacement, aft keep: the replacement is taken");
    sampler->derivedRecords(after);
    check(before.hasRange && after.rangeMin == before.rangeMin &&
            after.rangeMax == before.rangeMax &&
            kept(*sampler) == keptBefore,
          "data replacement, aft keep: the transform is held and the kept "
          "draws are not rewritten");
  }

  for (ext_rng* r : rngs) ext_rng_destroy(r);
  rngState = savedRngState;
  printf("ok: a whole-data replacement derives, keeps or takes each value\n");
}

void runStateTests(ext_rng* rng) {
  testFlattenRoundTrip();
  testCategoricalFlattenBoundaries();
  testKeepTrees(rng);
  testSavedDrawOrder(rng);
  testPredictCurrentTrees(rng);
  testStateRoundTrip();
  testStateRoundTripScaledOffset();
  testSetAnchorCarriesSigmaAndVarianceCalibration();
  testInstallLeavesHeldSigma();
  testShiftedStateInstalls();
  testStateRoundTripLatents(rng);
  testStateRoundTripStudentT(rng);
  testStateValidation(rng);
  testStateLatentFloor();
  testDegenerateGridRestores(rng);
  testInteractionContainment();
  testBlockAdditiveConfinement();
  testSingleForestColumnRestriction();
  testCrossGridWarmStart();
  testWarmStartForms();
  testWarmStartTreeCounts();
  testWarmStartLeafKinds();
  testWarmStartStandardization();
  testNoSpreadLiveConversion();
  testVarianceWarmStart();
  testVarianceWarmStartSlot();
  testStaleStateMerge();
  testRestoreStatus();
  testForcedConstraintRepair();
  testStateStoreSizes();
  testStackedSplitsMerge();
  testRepeatedGridRefused();
  testWeightedGridState();
  testLiveThresholdsOnGrid();
  testLiveGridInstall();
  testSuppliedCutGrid();
  testSuppliedStandardization();
  testDataReplacement();
  testLeafStandardizationRows();
  testRederiveStandardization();
  testStaleMissingDirectionBuild();
  testStaleMissingDirectionRestores();
  testVarianceForestPriorDraw();
  testCurrentVarianceRead();
  testVarianceSavedTreeState();
  testStateLeafScale(rng);
  testWeightsDigest();
  testMultiChainPartialFillPredict();
}
