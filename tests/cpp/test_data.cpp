// Needs the data and tree layers (cut points route through Tree/Rule), but
// not model/moves/chain/sampler/facade, so a touch there does not force
// this TU to recompile.
#include "assert.hpp"

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
      if (store.codes[i + j * n] != (xint_t) k) {
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

  check(view.isView, "a row-subset view is flagged as one");
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
        view.codes[i + j * rows.size()] == parent.codes[rows[i] + j * n];
  check(codesMatch, "view gathers subset codes");

  bool testCodesMatch = true;
  for (size_t i = 0; i < testRows.size() && testCodesMatch; ++i)
    for (size_t j = 0; j < p && testCodesMatch; ++j)
      testCodesMatch =
        view.testCodes[i * p + j] == parent.codes[testRows[i] + j * n];
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
    anyDiffer = rebuilt.codes[i] != view.codes[i];
  check(anyDiffer, "subset-built store bins differently than the view");

  printf("ok: column store view\n");
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

  printf("ok: column store leaf gather\n");
}

static void testColumnStoreMutation() {
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();

  ColumnStore store;
  store.build(x.data(), n, p, 100);
  std::vector<xint_t> originalCodes(store.codes);
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
    codesMatch &= store.codes[i] == originalCodes[i];
    codesMatch &= store.codes[i + n] == store.codeFor(1, newColumn[i]);
  }
  check(codesMatch, "setColumns re-quantizes only the target column");
  check(store.cutPoints[1].front() > 2.0 && store.cutPoints[1].back() < 5.0,
        "setColumns recomputes cuts over the new range");
  check(store.cutPoints[0] == originalCuts0, "setColumns leaves other cuts");

  // single-cell overwrite against existing cuts
  xint_t before = store.codes[7];
  store.setCell(7, 0, x[8]);
  check(store.codes[7] == originalCodes[8], "setCell re-quantizes one cell");
  check(before == originalCodes[7], "");  // silence unused warning

  // whole-matrix replacement without cut refresh: quantized on old cuts, codes
  // owned (the raw source is a call argument, retained nowhere)
  std::vector<double> x2(n * p);
  for (double& v : x2) v = runif01();
  store.setPredictors(x2.data(), false);
  codesMatch = true;
  for (size_t i = 0; i < n; ++i)
    codesMatch &= store.codes[i] == store.codeFor(0, x2[i]);
  check(codesMatch, "setPredictors re-quantizes against existing cuts");

  printf("ok: column store mutation\n");
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

  std::vector<size_t> indices(n);
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
    discreteCodesMatch &= store.codes[i + n] == static_cast<xint_t>(i % 10);
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
  std::vector<size_t> indices(n);
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

  std::vector<size_t> indices(n);
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
    codesRight &= (store.codes[i] == naCode) == (i % 10 == 0);
    codesRight &=
      (store.codes[i + n] == static_cast<xint_t>(naCategory)) == (i % 7 == 0);
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

void runDataTests() {
  testColumnStoreCodes();
  testColumnStoreView();
  testColumnStoreLeafGather();
  testColumnStoreMutation();
  testSetCutPointsOrphan();
  testQuantileCutPoints();
  testMapOldCutPointsOntoNew();
  testMapOldCutPointsStarvedWeightedMerge();
  testMissingIngestion();
}
