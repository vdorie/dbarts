// Needs the data and tree layers (cut points route through Tree/Rule), but
// not model/moves/chain/sampler/facade, so a touch there does not force
// this TU to recompile.
#include "assert.hpp"

#include <cstring>
#include <limits>

#include <bartcore/data.hpp>
#include <bartcore/tree.hpp>

using namespace bartcore;

static void testColumnStoreCodes() {
  const size_t n = 500, p = 3;
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();

  ColumnStore store;
  store.build(x.data(), n, p, 100);

  // reference: linear scan against uniformly spaced cuts
  for (size_t j = 0; j < p; ++j) {
    const double* column = x.data() + j * n;
    double xMin = column[0], xMax = column[0];
    for (size_t i = 1; i < n; ++i) {
      xMin = std::min(xMin, column[i]);
      xMax = std::max(xMax, column[i]);
    }
    double increment = (xMax - xMin) / 101.0;
    for (size_t i = 0; i < n; ++i) {
      uint32_t k = 0;
      while (k < 100 && column[i] > xMin + (double) (k + 1) * increment) ++k;
      if (store.train.codes[i + j * n] != (xint_t) k) {
        check(false, "column store code mismatch");
        return;
      }
    }
  }
  check(true, "");
  printf("ok: column store codes\n");
}

static void testColumnStoreView() {
  const size_t n = 120, p = 3;
  std::vector<double> x(n * p);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    x[i + 2 * n] = (double) (i % 4);  // categorical codes 0..3
  }
  // rows 0 and 1 carry column 0's extremes; the subset below excludes them,
  // so a store built over only the subset would bin that column differently
  x[0] = 0.0;
  x[1] = 1.0;

  std::vector<ColumnType> types = {ColumnType::ordinal, ColumnType::ordinal,
                                   ColumnType::categorical};
  ColumnStore parent;
  parent.build(x.data(), n, p, 25, false, types.data());

  std::vector<size_t> rows, testRows;
  for (size_t i = 2; i < n; i += 2) rows.push_back(i);
  testRows.push_back(0);  // the extremes land in the test rows
  testRows.push_back(1);
  for (size_t i = 3; i < n; i += 2) testRows.push_back(i);

  ColumnStore view;
  view.buildFromParent(parent, rows.data(), rows.size(), testRows.data(),
                       testRows.size());

  check(view.isView, "a row-subset view records its provenance");
  check(!view.hasRequantizeSource() && !view.acceptsNewRawPredictors(),
        "a view retains no re-quantize source and refuses raw mutation");
  check(parent.hasRequantizeSource() && parent.acceptsNewRawPredictors(),
        "a dense-built parent keeps both capabilities");
  check(view.numObservations == rows.size() &&
        view.numTestObservations == testRows.size() &&
        view.numPredictors == p, "view dimensions");
  check(view.types == parent.types && view.numCuts == parent.numCuts &&
        view.cutPoints == parent.cutPoints &&
        view.maxNumCuts == parent.maxNumCuts,
        "view copies the parent cut grid");

  bool codesMatch = true;
  for (size_t j = 0; j < p && codesMatch; ++j)
    for (size_t i = 0; i < rows.size() && codesMatch; ++i)
      codesMatch =
        view.train.codes[i + j * rows.size()] == parent.train.codes[rows[i] + j * n];
  check(codesMatch, "view gathers subset codes");

  bool testCodesMatch = true;
  for (size_t i = 0; i < testRows.size() && testCodesMatch; ++i)
    for (size_t j = 0; j < p && testCodesMatch; ++j)
      testCodesMatch =
        view.testCodeAt(j, i) == parent.train.codes[testRows[i] + j * n];
  check(testCodesMatch, "view gathers test codes from parent rows");

  // demonstrate the property matters: a store built over the subset's raw
  // values bins column 0 onto a different grid than the view keeps
  std::vector<double> subsetX(rows.size() * p);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < rows.size(); ++i)
      subsetX[i + j * rows.size()] = x[rows[i] + j * n];
  ColumnStore rebuilt;
  rebuilt.build(subsetX.data(), rows.size(), p, 25, false, types.data());
  bool anyDiffer = false;
  for (size_t i = 0; i < rows.size() && !anyDiffer; ++i)
    anyDiffer = rebuilt.train.codes[i] != view.train.codes[i];
  check(anyDiffer, "subset-built store bins differently than the view");

  printf("ok: column store view\n");
}

// A view may span a subset of the parent's columns; view-local column j binds
// to parent column columns[j]. The absent (all-columns) list must reproduce
// the default view byte-for-byte, and a subset must bin each of its columns
// identically to the mapped parent column, gathering leaf raw through the map.
static void testColumnStoreColumnSubset() {
  // leave the shared draw stream where downstream tests expect it
  uint64_t savedRngState = rngState;
  const size_t n = 120, p = 3;
  std::vector<double> x(n * p);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();  // the gathered leaf covariate
    x[i + 2 * n] = (double) (i % 4);  // categorical codes 0..3
  }
  x[0] = 0.0;  // column 0 extremes land in the test rows below, so a store
  x[1] = 1.0;  // built over only the training subset would bin it differently

  std::vector<ColumnType> types = {ColumnType::ordinal, ColumnType::ordinal,
                                   ColumnType::categorical};
  size_t covariate[] = {1};
  ColumnStore parent;
  parent.build(x.data(), n, p, 25, false, types.data(), covariate, 1);

  std::vector<size_t> rows, testRows;
  for (size_t i = 2; i < n; i += 2) rows.push_back(i);
  testRows.push_back(0);
  testRows.push_back(1);
  for (size_t i = 3; i < n; i += 2) testRows.push_back(i);

  // the identity list reproduces the default view byte-for-byte, leaf raw too
  size_t viewCovariate[] = {1};
  std::vector<size_t> allColumns = {0, 1, 2};
  ColumnStore viewDefault, viewFull;
  viewDefault.buildFromParent(parent, rows.data(), rows.size(),
                              testRows.data(), testRows.size(), viewCovariate,
                              1);
  viewFull.buildFromParent(parent, rows.data(), rows.size(), testRows.data(),
                           testRows.size(), viewCovariate, 1,
                           allColumns.data(), allColumns.size());
  check(viewFull.numPredictors == p && viewFull.train.codes == viewDefault.train.codes &&
          viewFull.cutPoints == viewDefault.cutPoints &&
          viewFull.numCuts == viewDefault.numCuts &&
          viewFull.types == viewDefault.types &&
          viewFull.maxNumCuts == viewDefault.maxNumCuts &&
          viewFull.test.codes == viewDefault.test.codes &&
          viewFull.gatheredRawColumns == viewDefault.gatheredRawColumns &&
          viewFull.gatheredRawValues == viewDefault.gatheredRawValues &&
          viewFull.gatheredMeans == viewDefault.gatheredMeans &&
          viewFull.gatheredSds == viewDefault.gatheredSds,
        "an all-columns view matches the default view path byte-for-byte");

  // a column subset spans only its columns, each binned identically to the
  // mapped parent column
  std::vector<size_t> subset = {2, 0};
  ColumnStore view;
  view.buildFromParent(parent, rows.data(), rows.size(), testRows.data(),
                       testRows.size(), nullptr, 0, subset.data(),
                       subset.size());
  check(view.numPredictors == subset.size(), "subset view spans its columns");
  bool gridMatches = true;
  for (size_t j = 0; j < subset.size(); ++j)
    gridMatches = gridMatches && view.types[j] == parent.types[subset[j]] &&
      view.numCuts[j] == parent.numCuts[subset[j]] &&
      view.cutPoints[j] == parent.cutPoints[subset[j]] &&
      view.maxNumCuts[j] == parent.maxNumCuts[subset[j]];
  check(gridMatches, "subset view copies each mapped parent column's cut grid");

  bool codesMatch = true;
  for (size_t j = 0; j < subset.size() && codesMatch; ++j)
    for (size_t i = 0; i < rows.size() && codesMatch; ++i)
      codesMatch =
        view.train.codes[i + j * rows.size()] == parent.codeAt(subset[j], rows[i]);
  for (size_t i = 0; i < testRows.size() && codesMatch; ++i)
    for (size_t j = 0; j < subset.size() && codesMatch; ++j)
      codesMatch = view.testCodeAt(j, i) == parent.codeAt(subset[j], testRows[i]);
  check(codesMatch, "subset view bins its columns identically to the parent");

  // a leaf covariate over a subset gathers through the same map: view-local
  // column 0 here is parent column 1, the covariate the parent owns raw for
  std::vector<size_t> gatherSubset = {1, 0};
  size_t gatherLocal[] = {0};
  ColumnStore gatherView;
  gatherView.buildFromParent(parent, rows.data(), rows.size(), testRows.data(),
                             testRows.size(), gatherLocal, 1,
                             gatherSubset.data(), gatherSubset.size());
  const double* raw = gatherView.rawColumn(0);
  bool gatherMatches = gatherView.gatheredRawColumns.size() == 1 &&
    gatherView.gatheredRawColumns[0] == 0 && raw != nullptr;
  for (size_t i = 0; i < rows.size() && gatherMatches; ++i)
    gatherMatches = raw[i] == x[rows[i] + n];
  double mean, sd, parentMean, parentSd;
  standardizationMomentsForColumn(x.data() + n, n, &parentMean, &parentSd);
  gatherMatches = gatherMatches &&
    gatherView.suppliedStandardization(0, &mean, &sd) && mean == parentMean &&
    sd == parentSd;
  check(gatherMatches,
        "subset view gathers a leaf covariate through the column map");

  rngState = savedRngState;
  printf("ok: column store column subset\n");
}

// The engine keeps no predictor matrix: build owns bitwise copies of the
// designated (leaf-covariate) columns so rawColumn serves them after the build
// borrow releases, and buildTest owns the test values for rawTestColumn.
static void testColumnStoreLeafGather() {
  const size_t n = 150, p = 4;
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();

  size_t gather[] = {1, 3};
  ColumnStore store;
  store.build(x.data(), n, p, 100, false, nullptr, gather, 2);

  bool gatheredMatch = true;
  const double* c1 = store.rawColumn(1);
  const double* c3 = store.rawColumn(3);
  if (c1 == nullptr || c3 == nullptr) {
    gatheredMatch = false;
  } else {
    for (size_t i = 0; i < n; ++i) {
      gatheredMatch &= c1[i] == x[i + 1 * n];
      gatheredMatch &= c3[i] == x[i + 3 * n];
    }
  }
  check(gatheredMatch, "leaf gather owns bitwise copies of designated columns");
  check(store.rawColumn(0) == nullptr && store.rawColumn(2) == nullptr,
        "undesignated columns own no raw");
  check(c1 != x.data() + n, "gathered raw is owned, not the source borrow");

  // the source going away must not disturb the owned copy
  std::fill(x.begin(), x.end(), -99.0);
  bool survives = true;
  for (size_t i = 0; i < n; ++i) survives &= store.rawColumn(3)[i] != -99.0;
  check(survives, "gathered raw survives the source borrow releasing");

  // a re-quantize refreshes the gathered copy from the new values
  std::vector<double> newCol1(n);
  for (double& v : newCol1) v = 2.0 + runif01();
  size_t col1 = 1;
  store.setColumns(newCol1.data(), &col1, 1, false);
  bool refreshed = true;
  const double* g = store.rawColumn(1);
  for (size_t i = 0; i < n; ++i) refreshed &= g[i] == newCol1[i];
  check(refreshed, "the gathered copy refreshes on re-quantize");

  // buildTest owns the test values; rawTestColumn serves the owned copy
  const size_t numTest = 20;
  std::vector<double> xTest(numTest * p);
  for (double& v : xTest) v = runif01();
  store.buildTest(xTest.data(), numTest);
  bool testOwned = store.rawTestColumn(1) != nullptr;
  const double* t1 = store.rawTestColumn(1);
  for (size_t i = 0; i < numTest && testOwned; ++i)
    testOwned &= t1[i] == xTest[i + 1 * numTest];
  check(testOwned, "buildTest owns test values for rawTestColumn");

  // the storage-aware test accessor reproduces the retired row-major code at
  // every (var, i): codeFor against the owned test value, dense-backed
  bool testCodesMatch = true;
  for (size_t j = 0; j < p && testCodesMatch; ++j)
    for (size_t i = 0; i < numTest && testCodesMatch; ++i)
      testCodesMatch =
        store.testCodeAt(j, i) == store.codeFor(j, xTest[i + j * numTest]);
  check(testCodesMatch, "testCodeAt equals the dense row-major test code");

  printf("ok: column store leaf gather\n");
}

static void testColumnStoreMutation() {
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();

  ColumnStore store;
  store.build(x.data(), n, p, 100);
  std::vector<xint_t> originalCodes(store.train.codes);
  std::vector<double> originalCuts0(store.cutPoints[0]);

  // column overwrite with cut refresh: column 0 codes untouched, column 1
  // re-quantized against cuts spanning the new range
  std::vector<double> newColumn(n);
  for (double& v : newColumn) v = 2.0 + 3.0 * runif01();
  size_t columnIndex = 1;
  store.setColumns(newColumn.data(), &columnIndex, 1, true);

  bool codesMatch = true;
  for (size_t i = 0; i < n; ++i) {
    // column 0's codes untouched; the store owns codes and never writes the
    // new values back into the caller's matrix (no write-through)
    codesMatch &= store.train.codes[i] == originalCodes[i];
    codesMatch &= store.train.codes[i + n] == store.codeFor(1, newColumn[i]);
  }
  check(codesMatch, "setColumns re-quantizes only the target column");
  check(store.cutPoints[1].front() > 2.0 && store.cutPoints[1].back() < 5.0,
        "setColumns recomputes cuts over the new range");
  check(store.cutPoints[0] == originalCuts0, "setColumns leaves other cuts");

  // single-cell overwrite against existing cuts
  xint_t before = store.train.codes[7];
  store.setCell(7, 0, x[8]);
  check(store.train.codes[7] == originalCodes[8], "setCell re-quantizes one cell");
  check(before == originalCodes[7], "");  // silence unused warning

  // whole-matrix replacement without cut refresh: quantized on old cuts, codes
  // owned (the raw source is a call argument, retained nowhere)
  std::vector<double> x2(n * p);
  for (double& v : x2) v = runif01();
  store.setPredictors(x2.data(), false);
  codesMatch = true;
  for (size_t i = 0; i < n; ++i)
    codesMatch &= store.train.codes[i] == store.codeFor(0, x2[i]);
  check(codesMatch, "setPredictors re-quantizes against existing cuts");

  printf("ok: column store mutation\n");
}

// codeFor's ordinal arm binary-searches the cuts; it must return the exact code
// the old linear scan did for every relation to the grid, ties and NA included.
// Deterministic inputs and a restored rngState keep the shared stream (and so
// the downstream snapshots) untouched.
static void testCodeForOrdinalBoundaries() {
  uint64_t savedRngState = rngState;
  const size_t n = 40;
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i);  // 0..39
  ColumnStore store;
  store.build(x.data(), n, 1, 8);  // ordinal, eight uniform cuts

  const std::vector<double>& cuts = store.cutPoints[0];
  const uint32_t numCuts = store.numCuts[0];

  // the pre-lower_bound codeFor transcribed: NA reserved, else the count of
  // cuts strictly below value
  auto oldCodeFor = [&](double value) -> xint_t {
    if (value != value) return naCode;
    uint32_t k = 0;
    while (k < numCuts && value > cuts[k]) ++k;
    return static_cast<xint_t>(k);
  };

  double probes[] = {
    cuts.front() - 1.0,          // below all cuts
    cuts[numCuts / 2],           // equal to an interior cut (tie)
    0.5 * (cuts[0] + cuts[1]),   // strictly between two cuts
    cuts.back() + 1.0,           // above all cuts
    std::nan(""),                // missing
  };
  bool match = true;
  for (double v : probes) match &= store.codeFor(0, v) == oldCodeFor(v);
  check(match, "codeFor ordinal binary search matches the linear scan");
  check(store.codeFor(0, std::nan("")) == naCode, "codeFor keeps the NA code");

  rngState = savedRngState;
  printf("ok: codeFor ordinal boundaries\n");
}

// setCutPoints shrinking an ordinal grid must not leave a split indexing past
// the new cuts: a rule sending its column's missing values right keeps both
// children occupied after the shrink, so an empty-child collapse alone would
// spare it and flatten would then read past cutPoints[j].
static void testSetCutPointsOrphan() {
  const size_t n = 12;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < 8; ++i) x[i] = static_cast<double>(i + 1);  // 1..8
  for (size_t i = 8; i < n; ++i) x[i] = std::nan("");

  ColumnStore store;
  store.build(x.data(), n, 1, 7, true);  // quantile cuts {1.5, ..., 7.5}
  check(store.hasMissing[0], "column carries missing values");

  std::vector<index_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rule;  rule.variableIndex = 0;  rule.setSplitIndex(6);
  rule.setMissingGoesRight(true);
  tree.birth(store, 0, rule, y.data(), nullptr);
  check(tree.at(tree.at(0).leftChild).numObservations() > 0 &&
        tree.at(tree.at(0).leftChild + 1).numObservations() > 0,
        "both children occupied before the shrink");

  // shrink below the split index; missing still routes right, so neither child
  // empties and the empty-child collapse alone would leave the split orphaned
  double newCuts[] = {2.0, 5.0};
  store.setCutPointsForColumn(0, newCuts, 2, x.data());
  check(store.hasMissing[0], "missing survives the re-quantize");

  std::vector<double> params(tree.nodes.size(), 0.0);
  tree.dropStaleMissingDirections(store);
  tree.repartitionSubtree(store, 0);
  tree.collapseEmptyNodes(store, nullptr, params);

  bool inRange = true;
  for (const Node& node : tree.nodes)
    if (!node.isBottom() &&
        store.types[static_cast<size_t>(node.rule.variableIndex)] !=
          ColumnType::categorical)
      inRange &= node.rule.splitIndex() <
                 static_cast<int32_t>(store.numCuts[static_cast<size_t>(
                   node.rule.variableIndex)]);
  check(inRange, "no split indexes past the shrunken grid");

  printf("ok: setCutPoints orphan\n");
}

static void testQuantileCutPoints() {
  const size_t n = 200;
  // column 0: continuous with more uniques than cuts; column 1: 10 discrete
  // levels, so quantile mode induces exactly 9 midpoint cuts
  std::vector<double> x(n * 2);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = static_cast<double>(i % 10);
  }

  ColumnStore store;
  store.build(x.data(), n, 2, 20, true);

  check(store.numCuts[0] == 20, "quantile cuts capped at maxNumCuts");
  check(store.numCuts[1] == 9, "few uniques induce numUnique - 1 cuts");
  bool discreteCutsMatch = true;
  for (std::uint32_t k = 0; k < 9; ++k)
    discreteCutsMatch &= store.cutPoints[1][k] == static_cast<double>(k) + 0.5;
  check(discreteCutsMatch, "discrete quantile cuts are unique-value midpoints");
  bool discreteCodesMatch = true;
  for (size_t i = 0; i < n; ++i)
    discreteCodesMatch &= store.train.codes[i + n] == static_cast<xint_t>(i % 10);
  check(discreteCodesMatch, "discrete quantile codes are value ranks");

  // continuous column: reference the thinning directly
  std::vector<double> sorted(x.begin(), x.begin() + n);
  std::sort(sorted.begin(), sorted.end());
  size_t step = n / 20, offset = step / 2;
  bool continuousCutsMatch = true;
  for (std::uint32_t k = 0; k < 20; ++k) {
    size_t index = std::min(static_cast<size_t>(k) * step + offset, n - 2);
    continuousCutsMatch &=
      store.cutPoints[0][k] == 0.5 * (sorted[index] + sorted[index + 1]);
  }
  check(continuousCutsMatch, "continuous quantile cuts thin sorted uniques");

  // refresh feasibility: fewer uniques than existing cuts is invalid
  std::vector<double> coarse(n);
  for (size_t i = 0; i < n; ++i) coarse[i] = static_cast<double>(i % 4);
  check(!store.cutsWouldRemainValid(1, coarse.data()),
        "coarser column fails the quantile feasibility check");
  std::vector<double> finer(n);
  for (size_t i = 0; i < n; ++i) finer[i] = static_cast<double>(i % 25);
  check(store.cutsWouldRemainValid(1, finer.data()),
        "finer column passes the quantile feasibility check");

  printf("ok: quantile cut points\n");
}

static void testMapOldCutPointsOntoNew() {
  // quantile grids over 1..8 give cuts {1.5, ..., 7.5}; hand-check the remap
  const size_t n = 8;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i + 1);

  ColumnStore store;
  store.build(x.data(), n, 1, 7, true);
  std::vector<std::vector<double>> oldCuts(store.cutPoints);

  // root splits at index 3 (cut 4.5), its left child at index 1 (cut 2.5)
  std::vector<index_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rootRule;  rootRule.variableIndex = 0;  rootRule.setSplitIndex(3);
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  Rule leftRule;  leftRule.variableIndex = 0;  leftRule.setSplitIndex(1);
  tree.birth(store, tree.at(0).leftChild, leftRule, y.data(), nullptr);
  int32_t leftChild = tree.at(0).leftChild;

  // new values 2..16 by twos: cuts {3, 5, ..., 15}; 4.5 is nearest 5 (index
  // 1), and 2.5 clamps into the left child's shrunken interval [0, 1)
  std::vector<double> x2(n);
  for (size_t i = 0; i < n; ++i) x2[i] = 2.0 * static_cast<double>(i + 1);
  store.setData(x2.data(), n);

  std::vector<double> params(tree.nodes.size(), 0.0);
  tree.mapOldCutPointsOntoNew(store, oldCuts, params);
  check(tree.at(0).rule.splitIndex() == 1, "root split remaps to nearest cut");
  check(tree.at(leftChild).rule.splitIndex() == 0,
        "child split clamps into the ancestor-constrained interval");

  // shift the grid entirely above the old cuts: the root clamps to index 0,
  // leaving the left child no interval, so it collapses with plain-mean param
  std::vector<double> x3(n);
  for (size_t i = 0; i < n; ++i) x3[i] = 20.0 + 2.0 * static_cast<double>(i);
  oldCuts = store.cutPoints;
  int32_t oldRootIndex = tree.at(0).rule.splitIndex();
  check(oldRootIndex == 1, "");  // silence unused in release
  store.setData(x3.data(), n);

  params.assign(tree.nodes.size(), 0.0);
  std::vector<int32_t> bottoms;
  tree.fillBottom(leftChild, bottoms);
  for (size_t k = 0; k < bottoms.size(); ++k)
    params[static_cast<size_t>(bottoms[k])] = static_cast<double>(k + 1);
  tree.mapOldCutPointsOntoNew(store, oldCuts, params);
  check(tree.at(0).rule.splitIndex() == 0, "root split clamps to the low end");
  check(tree.at(leftChild).isBottom(), "interval-starved subtree collapses");
  checkNear(params[static_cast<size_t>(leftChild)], 1.5, 1e-15,
            "collapsed subtree takes the plain mean of its leaf parameters");

  printf("ok: mapOldCutPointsOntoNew\n");
}

static void testMapOldCutPointsStarvedWeightedMerge() {
  // an interval-starved collapse merges its occupied leaves by effective
  // observation count, not the plain mean: split unevenly so the two leaves
  // carry unequal weight, then starve the subtree and check the merged param
  const size_t n = 8;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i + 1);

  ColumnStore store;
  store.build(x.data(), n, 1, 7, true);
  std::vector<std::vector<double>> oldCuts(store.cutPoints);

  std::vector<index_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rootRule;  rootRule.variableIndex = 0;  rootRule.setSplitIndex(5);
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  Rule leftRule;  leftRule.variableIndex = 0;  leftRule.setSplitIndex(0);
  tree.birth(store, tree.at(0).leftChild, leftRule, y.data(), nullptr);
  int32_t leftChild = tree.at(0).leftChild;

  std::vector<int32_t> bottoms;
  tree.fillBottom(leftChild, bottoms);
  check(bottoms.size() == 2, "starved subtree has two occupied leaves");
  std::vector<double> params(tree.nodes.size(), 0.0);
  double num = 0.0, den = 0.0;
  for (size_t k = 0; k < bottoms.size(); ++k) {
    double p = static_cast<double>(2 * k);
    params[static_cast<size_t>(bottoms[k])] = p;
    double w = tree.at(bottoms[k]).sumWeights;
    num += w * p;
    den += w;
  }
  double weightedMean = num / den;
  double plainMean = (params[static_cast<size_t>(bottoms[0])] +
                      params[static_cast<size_t>(bottoms[1])]) / 2.0;
  check(tree.at(bottoms[0]).sumWeights != tree.at(bottoms[1]).sumWeights,
        "leaves carry unequal weight");
  check(weightedMean != plainMean, "weighted and plain merges differ");

  // shift the grid entirely above the old cuts: the root clamps to index 0,
  // leaving the left child no interval, so its subtree collapses
  std::vector<double> x2(n);
  for (size_t i = 0; i < n; ++i) x2[i] = 20.0 + 2.0 * static_cast<double>(i);
  store.setData(x2.data(), n);

  tree.mapOldCutPointsOntoNew(store, oldCuts, params);
  check(tree.at(leftChild).isBottom(), "interval-starved subtree collapses");
  checkNear(params[static_cast<size_t>(leftChild)], weightedMean, 1e-15,
            "collapsed subtree takes the effective-count-weighted mean");

  printf("ok: mapOldCutPointsStarvedWeightedMerge\n");
}

static void testMissingIngestion() {
  const size_t n = 100;
  double na = std::nan("");
  // column 0: ordinal with NAs; column 1: categorical with NAs; column 2:
  // complete ordinal
  std::vector<double> x(n * 3);
  for (size_t i = 0; i < n; ++i) {
    x[i] = i % 10 == 0 ? na : runif01();
    x[i + n] = i % 7 == 0 ? na : static_cast<double>(i % 4);
    x[i + 2 * n] = runif01();
  }

  ColumnType types[] = {ColumnType::ordinal, ColumnType::categorical,
                        ColumnType::ordinal};
  ColumnStore store;
  store.build(x.data(), n, 3, 10, false, types);

  check(store.hasMissing[0] == 1 && store.hasMissing[1] == 1 &&
          store.hasMissing[2] == 0,
        "hasMissing marks exactly the columns with NAs");
  check(store.numCuts[1] == 4,
        "categorical categories are counted over observed values");

  bool codesRight = true;
  for (size_t i = 0; i < n; ++i) {
    codesRight &= (store.train.codes[i] == naCode) == (i % 10 == 0);
    codesRight &=
      (store.train.codes[i + n] == static_cast<xint_t>(naCategory)) == (i % 7 == 0);
  }
  check(codesRight, "missing values take the reserved codes");

  double loObserved = 2.0, hiObserved = -1.0;
  for (size_t i = 0; i < n; ++i) {
    if (i % 10 == 0) continue;
    if (x[i] < loObserved) loObserved = x[i];
    if (x[i] > hiObserved) hiObserved = x[i];
  }
  check(store.cutPoints[0].front() > loObserved &&
          store.cutPoints[0].back() < hiObserved,
        "uniform cuts span the observed range only");

  ColumnStore quantileStore;
  quantileStore.build(x.data(), n, 3, 10, true, types);
  check(quantileStore.numCuts[0] == 10 && quantileStore.hasMissing[0] == 1,
        "quantile grids skip NaN and keep the requested count");

  check(store.categoricalValueIsValid(1, na),
        "a missing categorical value is representable");
  printf("ok: missing ingestion\n");
}

// The bridge's dense-container ingestion assembles a transient contiguous
// block from per-column sources - numeric columns copied through, 1-based
// factor codes coerced to zero-based doubles, missing codes to NaN -
// instead of quantizing a retained cbound matrix. A store built from the
// assembled block must match a store built from the cbind reference
// bitwise, codes and cut grid both.
static void testTransientBlockAssembly() {
  const size_t n = 300, p = 3;
  std::vector<double> numeric1(n), numeric2(n);
  std::vector<int> factorCodes(n);  // 1-based, as an R factor holds them
  const int naCode = -2147483647 - 1;  // INT_MIN, R's NA_INTEGER
  for (size_t i = 0; i < n; ++i) {
    numeric1[i] = runif01();
    numeric2[i] = 10.0 * runif01() - 5.0;
    factorCodes[i] = 1 + (int) (i % 4);
  }
  factorCodes[17] = naCode;

  // the cbind reference: what makeCategoricalModelMatrix's retained matrix
  // held - as.double(numeric), as.double(as.integer(factor) - 1)
  std::vector<double> reference(n * p);
  for (size_t i = 0; i < n; ++i) {
    reference[i] = numeric1[i];
    reference[i + n] = factorCodes[i] == naCode
      ? std::numeric_limits<double>::quiet_NaN()
      : (double) (factorCodes[i] - 1);
    reference[i + 2 * n] = numeric2[i];
  }

  // the bridge assembly: per-column copies with the factor coercion
  std::vector<double> assembled(n * p);
  std::memcpy(assembled.data(), numeric1.data(), n * sizeof(double));
  for (size_t i = 0; i < n; ++i)
    assembled[i + n] = factorCodes[i] == naCode
      ? std::numeric_limits<double>::quiet_NaN()
      : (double) (factorCodes[i] - 1);
  std::memcpy(assembled.data() + 2 * n, numeric2.data(),
              n * sizeof(double));

  std::vector<ColumnType> types = {ColumnType::ordinal,
                                   ColumnType::categorical,
                                   ColumnType::ordinal};
  ColumnStore fromReference, fromAssembled;
  fromReference.build(reference.data(), n, p, 100, false, types.data());
  fromAssembled.build(assembled.data(), n, p, 100, false, types.data());

  check(fromAssembled.train.codes == fromReference.train.codes,
        "container-assembled block quantizes to the cbind reference codes");
  check(fromAssembled.cutPoints == fromReference.cutPoints,
        "container-assembled block builds the cbind reference cut grid");
  printf("ok: transient block assembly\n");
}

// The bitwise device on the test side: a mixed dense + CSC test set, built
// against a training cut grid through a mapped test view, bins IDENTICALLY to
// a dense test matrix of the same values. Both storage tiers appear - a rank
// column (nonzero fraction at or below the threshold) and a densified one - so
// the check covers testCodeAt on each, the storage-aware test descent through a
// grown tree, and one recorded test fit (the leaf parameter each row lands on).
static void testSparseTestColumnStore() {
  uint64_t savedRngState = rngState;  // leave the shared draw stream in place
  const size_t nTrain = 400, numTest = 200, p = 4;

  // training is a plain dense ordinal matrix; its cut grid is what the test
  // codes quantize against, shared by identity
  std::vector<double> xTrain(nTrain * p);
  for (double& v : xTrain) v = runif01();
  ColumnStore store;
  store.build(xTrain.data(), nTrain, p, 25);

  // test columns: 0 and 3 dense-backed, 1 a rank-tier CSC column (~8% stored),
  // 2 a densified-tier CSC column (~60% stored)
  std::vector<double> dense0(numTest), dense3(numTest);
  for (size_t i = 0; i < numTest; ++i) {
    dense0[i] = runif01();
    dense3[i] = runif01();
  }
  std::vector<int> pointers(3, 0), rows;
  std::vector<double> values;
  const double fractions[] = {0.08, 0.6};
  for (int csc = 0; csc < 2; ++csc) {
    for (size_t i = 0; i < numTest; ++i)
      if (runif01() < fractions[csc]) {
        rows.push_back(static_cast<int>(i));
        values.push_back(0.3 + runif01());
      }
    pointers[static_cast<size_t>(csc) + 1] = static_cast<int>(rows.size());
  }

  // the dense-matrix reference: CSC columns densified with 0.0 in implicit rows
  std::vector<double> denseTest(numTest * p, 0.0);
  for (size_t i = 0; i < numTest; ++i) {
    denseTest[i + 0 * numTest] = dense0[i];
    denseTest[i + 3 * numTest] = dense3[i];
  }
  for (int csc = 0; csc < 2; ++csc)
    for (int k = pointers[static_cast<size_t>(csc)];
         k < pointers[static_cast<size_t>(csc) + 1]; ++k)
      denseTest[static_cast<size_t>(rows[static_cast<size_t>(k)]) +
                static_cast<size_t>(csc + 1) * numTest] =
        values[static_cast<size_t>(k)];

  // grow a tree over the training grid, splitting on a dense-backed column and
  // both CSC-backed columns so the descent visits every test storage kind
  std::vector<index_t> indices(nTrain);
  std::vector<double> yTrain(nTrain, 0.0);
  Tree tree;
  tree.initialize(indices.data(), nTrain);
  Rule root;  root.variableIndex = 0;  root.setSplitIndex(12);
  tree.birth(store, 0, root, yTrain.data(), nullptr);
  Rule left;  left.variableIndex = 1;  left.setSplitIndex(5);
  tree.birth(store, tree.at(0).leftChild, left, yTrain.data(), nullptr);
  Rule right;  right.variableIndex = 2;  right.setSplitIndex(5);
  tree.birth(store, tree.at(0).leftChild + 1, right, yTrain.data(), nullptr);

  std::vector<double> params(tree.nodes.size());
  for (size_t node = 0; node < params.size(); ++node)
    params[node] = static_cast<double>(node) + 0.5;

  // dense test first: snapshot its codes, leaf indices, and recorded fits
  store.buildTest(denseTest.data(), numTest);
  std::vector<xint_t> denseCodes(p * numTest);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < numTest; ++i)
      denseCodes[j * numTest + i] = store.testCodeAt(j, i);
  std::vector<int32_t> denseLeaf(numTest);
  std::vector<double> denseFit(numTest);
  for (size_t i = 0; i < numTest; ++i) {
    denseLeaf[i] = tree.findBottomNodeForRow(store, i);
    denseFit[i] = params[static_cast<size_t>(denseLeaf[i])];
  }

  // the same values through the mixed build: dense sources 0 and 1 carry test
  // columns 0 and 3, CSC sources 0 and 1 carry columns 1 and 2
  std::vector<double> denseBlock(numTest * 2);
  std::memcpy(denseBlock.data(), dense0.data(), numTest * sizeof(double));
  std::memcpy(denseBlock.data() + numTest, dense3.data(),
              numTest * sizeof(double));
  std::vector<std::int32_t> columnSources = {0, ~0, ~1, 1};
  bartcore::PredictorSource testSource;
  testSource.numRows = numTest;
  testSource.numColumns = p;
  testSource.denseValues = denseBlock.data();
  testSource.cscColumnPointers = pointers.data();
  testSource.cscRowIndices = rows.data();
  testSource.cscValues = values.data();
  testSource.columnSources = columnSources.data();
  store.buildTest(testSource);

  check(store.testColumnIsSparse(1) && !store.testColumnIsSparse(0) &&
        !store.testColumnIsSparse(2) && !store.testColumnIsSparse(3),
        "the density threshold splits the test storage tiers");

  bool codesMatch = true;
  for (size_t j = 0; j < p && codesMatch; ++j)
    for (size_t i = 0; i < numTest && codesMatch; ++i)
      codesMatch = store.testCodeAt(j, i) == denseCodes[j * numTest + i];
  check(codesMatch, "mixed test codes match the dense test matrix at every cell");

  check(store.rawTestColumn(0) != nullptr && store.rawTestColumn(3) != nullptr &&
        store.rawTestColumn(1) == nullptr && store.rawTestColumn(2) == nullptr,
        "dense-backed test columns serve raw; CSC-backed refuse it");

  bool leavesMatch = true, fitsMatch = true;
  for (size_t i = 0; i < numTest; ++i) {
    int32_t leaf = tree.findBottomNodeForRow(store, i);
    leavesMatch &= leaf == denseLeaf[i];
    fitsMatch &= params[static_cast<size_t>(leaf)] == denseFit[i];
  }
  check(leavesMatch, "test descent lands on the dense path's leaf every row");
  check(fitsMatch, "recorded test fits match the dense path every row");

  rngState = savedRngState;
  printf("ok: sparse test column store\n");
}

// The CSC-categorical code contract on the test side: a sparse categorical test
// column's implicit rows carry the reference level's ACTUAL level-order code
// (never a numeric zero), so its zeroCode equals that code and testCodeAt
// reproduces the densified factor at every row. The reference is a middle level
// so its code is not zero.
static void testSparseCategoricalTestColumnStore() {
  uint64_t savedRngState = rngState;
  const size_t nTrain = 300, numTest = 200;
  const std::uint32_t K = 6;
  const xint_t reference = static_cast<xint_t>(K / 2);

  // training plants every level so the categorical count is K
  std::vector<double> xTrain(nTrain);
  for (size_t i = 0; i < nTrain; ++i)
    xTrain[i] = static_cast<double>(i % K);
  ColumnType type = ColumnType::categorical;
  ColumnStore store;
  store.build(xTrain.data(), nTrain, 1, 100, false, &type);
  check(store.numCuts[0] == K, "training counts the categorical levels");

  // test: the reference level is implicit, other levels stored (~10%), so the
  // column lands rank-tier
  std::vector<double> denseTest(numTest);
  std::vector<int> pointers(2, 0), rows;
  std::vector<double> values;
  for (size_t i = 0; i < numTest; ++i) {
    xint_t code = reference;
    if (runif01() < 0.1) {
      code = static_cast<xint_t>(runif01() * static_cast<double>(K - 1));
      if (code >= reference) ++code;  // uniform over the non-reference levels
      rows.push_back(static_cast<int>(i));
      values.push_back(static_cast<double>(code));
    }
    denseTest[i] = static_cast<double>(code);
  }
  pointers[1] = static_cast<int>(rows.size());

  store.buildTest(denseTest.data(), numTest);
  std::vector<xint_t> denseCodes(numTest);
  for (size_t i = 0; i < numTest; ++i) denseCodes[i] = store.testCodeAt(0, i);

  std::int32_t source = ~0;
  bartcore::PredictorSource testSource;
  testSource.numRows = numTest;
  testSource.numColumns = 1;
  testSource.cscColumnPointers = pointers.data();
  testSource.cscRowIndices = rows.data();
  testSource.cscValues = values.data();
  testSource.columnSources = &source;
  testSource.referenceCodes = &reference;
  store.buildTest(testSource);

  check(store.testColumnIsSparse(0), "the sparse categorical test column is rank-tier");
  check(reference != 0, "the reference level's code is not a numeric zero");
  check(store.testSparseColumn(0).zeroCode == reference,
        "the test zero code carries the reference level's own code");

  bool codesMatch = true;
  for (size_t i = 0; i < numTest && codesMatch; ++i)
    codesMatch = store.testCodeAt(0, i) == denseCodes[i];
  check(codesMatch, "sparse categorical test codes match the densified factor");

  rngState = savedRngState;
  printf("ok: sparse categorical test column store\n");
}

// One predictor view drives every builder, so the collapsed
// ColumnStore::build/buildTest must produce the SAME store from any equivalent
// spelling of a source. Four equivalences over identical values: a dense build
// through an explicit null-map view; a dense build through an identity-MAPPED
// view (which additionally copies and owns its block); a dense TEST build
// through an identity-mapped view (the buildTest/buildTestMixed collapse); and
// an all-negative-map (bare CSC) build against the DENSE-EQUIVALENT build on
// both storage tiers. A categorical column carries a declared level count
// wider than its observed codes throughout, so the declared-K channel rides
// each spelling too.
static void testPredictorViewEquivalence() {
  uint64_t savedRngState = rngState;
  const size_t n = 260, numTest = 90, p = 4;
  const std::uint32_t declaredK = 5;  // the observed codes reach 2

  std::vector<double> x(n * p), xTest(numTest * p);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 3.0 * runif01() - 1.0;
    x[i + 2 * n] = static_cast<double>(i % 3);
    x[i + 3 * n] = runif01() < 0.3 ? 0.0 : runif01();
  }
  for (size_t i = 0; i < numTest; ++i) {
    xTest[i] = runif01();
    xTest[i + numTest] = 3.0 * runif01() - 1.0;
    xTest[i + 2 * numTest] = static_cast<double>(i % 3);
    xTest[i + 3 * numTest] = runif01();
  }
  const ColumnType types[p] = { ColumnType::ordinal, ColumnType::ordinal,
                                ColumnType::categorical, ColumnType::ordinal };
  const std::uint32_t counts[p] = { 0, 0, declaredK, 0 };
  std::vector<std::uint32_t> maxCuts(p, 20);
  std::vector<std::int32_t> identity(p);
  for (size_t j = 0; j < p; ++j) identity[j] = static_cast<std::int32_t>(j);

  auto gridsAgree = [](const ColumnStore& a, const ColumnStore& b) {
    return a.types == b.types && a.numCuts == b.numCuts &&
           a.cutPoints == b.cutPoints && a.maxNumCuts == b.maxNumCuts &&
           a.hasMissing == b.hasMissing;
  };
  auto trainCodesAgree = [](const ColumnStore& a, const ColumnStore& b) {
    for (size_t j = 0; j < a.numPredictors; ++j)
      for (size_t i = 0; i < a.numObservations; ++i)
        if (a.codeAt(j, i) != b.codeAt(j, i)) return false;
    return true;
  };
  auto testCodesAgree = [](const ColumnStore& a, const ColumnStore& b) {
    for (size_t j = 0; j < a.numPredictors; ++j)
      for (size_t i = 0; i < a.numTestObservations; ++i)
        if (a.testCodeAt(j, i) != b.testCodeAt(j, i)) return false;
    return true;
  };

  for (bool useQuantiles : { false, true }) {
    ColumnStore dense;
    dense.build(x.data(), n, p, maxCuts.data(), useQuantiles, types, nullptr, 0,
                counts);
    check(dense.numCuts[2] == declaredK,
          "the declared level count reaches the dense build");

    PredictorSource nullMap;
    nullMap.numRows = n;
    nullMap.numColumns = p;
    nullMap.denseValues = x.data();
    nullMap.columnTypes = types;
    nullMap.categoryCounts = counts;
    ColumnStore viaView;
    viaView.build(nullMap, maxCuts.data(), 0, useQuantiles);
    check(gridsAgree(dense, viaView) && trainCodesAgree(dense, viaView) &&
          !viaView.builtFromCsc && viaView.ownedDenseValues.empty(),
          "a null-map view builds the dense store, retaining no raw");

    PredictorSource mapped = nullMap;
    mapped.columnSources = identity.data();
    ColumnStore viaMap;
    viaMap.build(mapped, maxCuts.data(), 0, useQuantiles);
    check(gridsAgree(dense, viaMap) && trainCodesAgree(dense, viaMap),
          "an identity-mapped view bins exactly as the dense build");
    check(viaMap.ownedDenseValues.size() == n * p &&
          viaMap.rawColumn(1) == viaMap.ownedDenseValues.data() + n,
          "an identity-mapped build owns its dense block");

    dense.buildTest(xTest.data(), numTest);
    PredictorSource testMap;
    testMap.numRows = numTest;
    testMap.numColumns = p;
    testMap.denseValues = xTest.data();
    testMap.columnSources = identity.data();
    ColumnStore denseTestViaMap;
    denseTestViaMap.build(nullMap, maxCuts.data(), 0, useQuantiles);
    denseTestViaMap.buildTest(testMap);
    check(testCodesAgree(dense, denseTestViaMap),
          "an identity-mapped test view bins exactly as the dense test build");
    bool rawAgrees = true;
    for (size_t j = 0; j < p && rawAgrees; ++j)
      for (size_t i = 0; i < numTest; ++i)
        rawAgrees &= dense.rawTestColumn(j)[i] ==
                     denseTestViaMap.rawTestColumn(j)[i];
    check(rawAgrees, "both test spellings serve the same owned raw");
  }

  // the all-negative map: two CSC columns, one per storage tier, the densified
  // one categorical with a reference level past zero and a declared count
  // wider than its observed codes
  const size_t q = 2;
  const std::uint32_t cscDeclaredK = 7;  // the observed codes reach 4
  const xint_t reference = 2;
  std::vector<double> cscDense(n * q, 0.0);
  std::vector<int> pointers(q + 1, 0), rows;
  std::vector<double> values;
  for (size_t i = 0; i < n; ++i)
    if (runif01() < 0.08) {
      cscDense[i] = 0.5 + runif01();
      rows.push_back(static_cast<int>(i));
      values.push_back(cscDense[i]);
    }
  pointers[1] = static_cast<int>(rows.size());
  for (size_t i = 0; i < n; ++i) {
    xint_t code = reference;
    if (runif01() < 0.6) {
      code = static_cast<xint_t>(runif01() * 4.0);
      if (code >= reference) ++code;  // uniform over the non-reference levels
      rows.push_back(static_cast<int>(i));
      values.push_back(static_cast<double>(code));
    }
    cscDense[i + n] = static_cast<double>(code);
  }
  pointers[2] = static_cast<int>(rows.size());

  const ColumnType cscTypes[q] = { ColumnType::ordinal,
                                   ColumnType::categorical };
  const std::uint32_t cscCounts[q] = { 0, cscDeclaredK };
  const xint_t cscReferences[q] = { 0, reference };
  std::vector<std::int32_t> allCsc(q);
  for (size_t j = 0; j < q; ++j) allCsc[j] = ~static_cast<std::int32_t>(j);

  for (bool useQuantiles : { false, true }) {
    PredictorSource cscView;
    cscView.numRows = n;
    cscView.numColumns = q;
    cscView.cscColumnPointers = pointers.data();
    cscView.cscRowIndices = rows.data();
    cscView.cscValues = values.data();
    cscView.columnSources = allCsc.data();
    cscView.columnTypes = cscTypes;
    cscView.categoryCounts = cscCounts;
    cscView.referenceCodes = cscReferences;
    ColumnStore fromCsc;
    fromCsc.build(cscView, nullptr, 20, useQuantiles);
    ColumnStore fromDense;
    fromDense.build(cscDense.data(), n, q, 20, useQuantiles, cscTypes, nullptr,
                    0, cscCounts);
    check(fromCsc.columnIsSparse(0) && !fromCsc.columnIsSparse(1),
          "the all-negative map splits the two storage tiers");
    check(fromCsc.numCuts[1] == cscDeclaredK,
          "the declared level count reaches the CSC-backed categorical column");
    check(gridsAgree(fromDense, fromCsc) && trainCodesAgree(fromDense, fromCsc),
          "an all-negative-map view bins as the dense-equivalent build");
  }

  rngState = savedRngState;
  printf("ok: predictor view builder equivalence\n");
}

void runDataTests() {
  testColumnStoreCodes();
  testColumnStoreView();
  testColumnStoreColumnSubset();
  testColumnStoreLeafGather();
  testColumnStoreMutation();
  testCodeForOrdinalBoundaries();
  testSetCutPointsOrphan();
  testQuantileCutPoints();
  testMapOldCutPointsOntoNew();
  testMapOldCutPointsStarvedWeightedMerge();
  testMissingIngestion();
  testTransientBlockAssembly();
  testSparseTestColumnStore();
  testSparseCategoricalTestColumnStore();
  testPredictorViewEquivalence();
}
