// Needs the data/tree/model layers (CGMTreePrior) but not moves/chain/
// sampler/facade, so a touch there does not force this TU to recompile.
#include "assert.hpp"

#include <external/random.h>

#include <bartcore/data.hpp>
#include <bartcore/tree.hpp>
#include <bartcore/model.hpp>

using namespace bartcore;

// Build x with a known partition structure and verify splitting mechanics.
static void testTreeMechanics() {
  const size_t n = 100;
  std::vector<double> x(n);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = (double) i / (double) (n - 1);
    y[i] = x[i] < 0.5 ? 1.0 : 3.0;
  }

  ColumnStore store;
  store.build(x.data(), n, 1, 100);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  checkNear(tree.at(0).sumWeightedResponse / tree.at(0).sumWeights, 2.0, 1e-12,
            "root average");

  // split at the median cut
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(49);
  tree.birth(store, 0, rule, y.data(), nullptr);

  const Node& left(tree.at(tree.at(0).leftChild));
  const Node& right(tree.at(tree.at(0).leftChild + 1));
  check(left.numObservations() + right.numObservations() == n,
        "birth partitions all observations");

  // verify against a direct count with the store's own codes
  size_t expectedLeft = 0;
  for (size_t i = 0; i < n; ++i)
    if (store.train.codes[i] <= 49) ++expectedLeft;
  check(left.numObservations() == expectedLeft, "partition count matches codes");

  double leftSum = 0.0, rightSum = 0.0;
  for (size_t i = 0; i < n; ++i)
    (store.train.codes[i] <= 49 ? leftSum : rightSum) += y[i];
  checkNear(left.sumWeightedResponse, leftSum, 1e-12, "left child sum wz");
  checkNear(right.sumWeightedResponse, rightSum, 1e-12, "right child sum wz");

  // split interval of the left child is bounded by the parent's rule
  int32_t lo, hi;
  tree.splitInterval(store, tree.at(0).leftChild, 0, &lo, &hi);
  check(lo == 0 && hi == 48, "left child split interval");
  tree.splitInterval(store, tree.at(0).leftChild + 1, 0, &lo, &hi);
  check(lo == 50 && hi == 99, "right child split interval");

  // orphanChildren sums the children's sufficient statistics
  double mergedSumWZ = left.sumWeightedResponse + right.sumWeightedResponse;
  tree.orphanChildren(0);
  checkNear(tree.at(0).sumWeightedResponse, mergedSumWZ, 1e-12,
            "orphan merge sum wz");
  check(tree.at(0).isBottom(), "orphan clears children");

  printf("ok: tree mechanics\n");
}

static void testTreePriorMath() {
  const size_t n = 64;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = (double) i / (double) (n - 1);

  ColumnStore store;
  store.build(x.data(), n, 1, 100);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);

  CGMTreePrior prior;  // base 0.95, power 2

  checkNear(prior.growthProbability(tree, store, 0), 0.95, 1e-15,
            "root growth probability");

  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(49);
  tree.birth(store, 0, rule, y.data(), nullptr);
  int32_t left = tree.at(0).leftChild;

  checkNear(prior.growthProbability(tree, store, left), 0.95 / 4.0, 1e-15,
            "depth-1 growth probability");

  // single-split tree: log p = log(g0) + log(1/p) + log(1/numCuts)
  //                          + 2 log(1 - g1)
  double expected = std::log(0.95) + std::log(1.0) - std::log(100.0) +
                    2.0 * std::log(1.0 - 0.95 / 4.0);
  checkNear(prior.treeLogProbability(tree, store), expected, 1e-12,
            "tree log probability");

  printf("ok: tree prior math\n");
}

static void testCategoricalMechanics() {
  const size_t n = 120;
  // column 0: categorical with 4 levels; column 1: continuous
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 4);
    x[i + n] = runif01();
    y[i] = static_cast<double>(i % 4);
  }

  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
  ColumnStore store;
  store.build(x.data(), n, 2, 10, false, types);

  check(store.numCuts[0] == 4, "categorical column counts its categories");
  check(store.cutPoints[0].empty(), "categorical column keeps no cut points");
  bool codesMatch = true;
  for (size_t i = 0; i < n; ++i)
    codesMatch &= store.train.codes[i] == static_cast<xint_t>(i % 4);
  check(codesMatch, "categorical codes are the values");
  check(store.categoricalValueIsValid(0, 3.0) &&
          !store.categoricalValueIsValid(0, 4.0) &&
          !store.categoricalValueIsValid(0, 1.5),
        "categorical value validity");

  // rule sending {1, 3} right partitions by bit test
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rule;
  rule.variableIndex = 0;
  rule.setCategoryDirections((1u << 1) | (1u << 3));
  tree.birth(store, 0, rule, y.data(), nullptr);

  const Node& left(tree.at(tree.at(0).leftChild));
  const Node& right(tree.at(tree.at(0).leftChild + 1));
  check(left.numObservations() == n / 2 && right.numObservations() == n / 2,
        "categorical partition splits by mask");
  bool sidesMatch = true;
  for (size_t i = left.begin; i < left.end; ++i) {
    size_t category = tree.indices[i] % 4;
    sidesMatch &= category == 0 || category == 2;
  }
  check(sidesMatch, "left side holds only left-bound categories");
  checkNear(left.sumWeightedResponse / left.sumWeights, (0.0 + 2.0) / 2.0, 1e-12,
            "left average by category");

  // reachability filters through nested rules
  check(tree.reachableCategories(store, 0, 0) == 0xfu,
        "root reaches all categories");
  check(tree.reachableCategories(store, tree.at(0).leftChild, 0) ==
          ((1u << 0) | (1u << 2)),
        "left child reaches left-bound categories");
  check(tree.reachableCategories(store, tree.at(0).leftChild + 1, 0) ==
          ((1u << 1) | (1u << 3)),
        "right child reaches right-bound categories");
  check(tree.variableAvailable(store, tree.at(0).leftChild, 0),
        "two reachable categories keep the variable available");

  // a second split below exhausts one side
  Rule childRule;
  childRule.variableIndex = 0;
  childRule.setCategoryDirections(1u << 2);
  tree.birth(store, tree.at(0).leftChild, childRule, y.data(), nullptr);
  int32_t grandLeft = tree.at(tree.at(0).leftChild).leftChild;
  check(tree.reachableCategories(store, grandLeft, 0) == (1u << 0),
        "grandchild reaches a single category");
  check(!tree.variableAvailable(store, grandLeft, 0),
        "one reachable category exhausts the variable");

  // routing agrees with the partitions, including a code override
  bool routesMatch = true;
  for (size_t i = 0; i < n; i += 7) {
    int32_t leaf = tree.findBottomNodeForObservation(store, i);
    size_t inLeaf = 0;
    for (size_t k = tree.at(leaf).begin; k < tree.at(leaf).end; ++k)
      inLeaf += tree.indices[k] == i ? 1 : 0;
    routesMatch &= inLeaf == 1;
  }
  check(routesMatch, "categorical routing matches partitions");
  int32_t overridden = tree.findBottomNodeForObservation(
    store, 0, 0, static_cast<xint_t>(3));
  check(overridden == tree.at(0).leftChild + 1,
        "code override routes to the right side");

  printf("ok: categorical mechanics\n");
}

static void testCategoricalPriorMath(ext_rng* rng) {
  const size_t n = 120;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 5);

  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);

  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);

  CGMTreePrior prior;

  // R = 5 reachable: 2^5 - 2 = 30 valid assignments, drawn uniformly with
  // unreachable bits (none here) zero and neither side empty
  const int numDraws = 60000;
  std::vector<int> patternCounts(1u << 5, 0);
  for (int i = 0; i < numDraws; ++i) {
    Rule rule = prior.drawRuleForVariable(tree, store, rng, 0, 0);
    check(rule.categoryDirections() > 0 && rule.categoryDirections() < 31u + 1u &&
            rule.categoryDirections() != 31u,
          "categorical draw leaves neither side empty");
    ++patternCounts[rule.categoryDirections()];
  }
  bool uniform = patternCounts[0] == 0 && patternCounts[31] == 0;
  double expected = static_cast<double>(numDraws) / 30.0;
  for (std::uint32_t pattern = 1; pattern < 31; ++pattern)
    uniform = uniform &&
      std::fabs(patternCounts[pattern] - expected) < 5.0 * std::sqrt(expected);
  check(uniform, "categorical rule draw is uniform over valid assignments");

  Rule rootRule;
  rootRule.variableIndex = 0;
  rootRule.setCategoryDirections((1u << 3) | (1u << 4));
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  checkNear(prior.ruleForVariableLogProbability(tree, store, 0),
            -std::log(30.0), 1e-13, "root categorical rule probability");

  // left child reaches {0, 1, 2}: R = 3 gives 2^3 - 2 = 6
  Rule childRule;
  childRule.variableIndex = 0;
  childRule.setCategoryDirections(1u << 1);
  tree.birth(store, tree.at(0).leftChild, childRule, y.data(), nullptr);
  checkNear(prior.ruleForVariableLogProbability(tree, store, tree.at(0).leftChild),
            -std::log(6.0), 1e-13, "child categorical rule probability");

  // full tree log probability: root splits (g0), children of root: left
  // splits (g1), right leaf (1 - g1'); right of root reaches {3,4} so g1'
  // uses availability; grandchildren reach single categories, forced leaves
  double g0 = prior.growthProbability(tree, store, 0);
  double g1 = prior.growthProbability(tree, store, tree.at(0).leftChild);
  double gRight = prior.growthProbability(tree, store, tree.at(0).leftChild + 1);
  int32_t leftChild = tree.at(0).leftChild;
  double gGrandLeft = prior.growthProbability(tree, store, tree.at(leftChild).leftChild);
  double gGrandRight = prior.growthProbability(tree, store, tree.at(leftChild).leftChild + 1);
  // the left child's rule sends category 1 right: its left grandchild
  // reaches {0, 2} and can grow, the right reaches only {1} and cannot
  check(gGrandLeft != 0.0, "two-category node can grow");
  check(gGrandRight == 0.0, "single-category node cannot grow");
  double expectedLogProbability =
    std::log(g0) - std::log(30.0) +
    std::log(g1) - std::log(6.0) +
    std::log(1.0 - gRight) +
    std::log(1.0 - gGrandLeft) + std::log(1.0 - gGrandRight);
  checkNear(prior.treeLogProbability(tree, store), expectedLogProbability,
            1e-12, "categorical tree log probability");

  printf("ok: categorical prior math\n");
}

static void testPooledMaskMechanics(ext_rng* rng) {
  // columns past 63 categories store per-tree pooled mask words; exercise
  // the partition kernel, reachability, draws, the pool's compaction, and
  // the flattened side channel (docs/design/pooled-masks.md)
  const size_t K = 70;  // ceil(71 / 64) = 2 words
  const size_t n = 10 * K;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    size_t category = (i * 37) % K;
    x[i] = static_cast<double>(category);
    y[i] = (category % 3 == 0 ? 2.0 : 0.0) + 0.3 * (runif01() - 0.5);
  }

  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);
  check(store.numCuts[0] == K, "pooled column counts its categories");
  check(store.columnIsPooled(0) && store.hasPooledCategorical,
        "pooled tier predicates");
  check(maskWordsForCount(K) == 2, "70 categories need two words");
  check(missingCategoryCode(K) == 70 &&
          store.codeFor(0, std::nan("")) == 70,
        "a pooled column's missing code is K");

  // a hand-built rule whose direction bits straddle both words
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  size_t offset = tree.allocateMask(2);
  {
    std::uint64_t* words = tree.mutableMaskWordsFor(offset);
    maskSetBit(words, 1);
    maskSetBit(words, 40);
    maskSetBit(words, 64);
    maskSetBit(words, 69);
  }
  Rule rule;
  rule.variableIndex = 0;
  rule.setMaskOffset(offset);
  tree.birth(store, 0, rule, y.data(), nullptr);
  const Node& left(tree.at(tree.at(0).leftChild));
  const Node& right(tree.at(tree.at(0).leftChild + 1));
  size_t expectedRight = 0;
  for (size_t i = 0; i < n; ++i) {
    size_t category = (i * 37) % K;
    if (category == 1 || category == 40 || category == 64 || category == 69)
      ++expectedRight;
  }
  check(right.numObservations() == expectedRight &&
          left.numObservations() == n - expectedRight,
        "pooled mask partitions across words");
  bool sidesMatch = true;
  for (size_t i = right.begin; i < right.end; ++i) {
    size_t category = (tree.indices[i] * 37) % K;
    sidesMatch &= category == 1 || category == 40 || category == 64 ||
                  category == 69;
  }
  check(sidesMatch, "right side holds only right-bound categories");

  std::uint64_t reachable[2];
  tree.reachableCategoriesWide(store, tree.at(0).leftChild + 1, 0, reachable);
  check(maskEquals(reachable, tree.maskWordsFor(tree.at(0).rule), 2),
        "right child reaches the mask's categories");
  tree.reachableCategoriesWide(store, tree.at(0).leftChild, 0, reachable);
  check(maskPopcount(reachable, 2) == K - 4 && !maskTestBit(reachable, 40) &&
          maskTestBit(reachable, 0) && maskTestBit(reachable, 68),
        "left child reaches the complement");

  // each pooled direction bit is marginally 1/2 under the uniform prior on
  // assignments that leave neither side empty
  Tree drawTree;
  drawTree.initialize(indices.data(), n);
  CGMTreePrior prior;
  const int numDraws = 20000;
  std::vector<int> bitCounts(K, 0);
  bool neverEmpty = true;
  for (int d = 0; d < numDraws; ++d) {
    size_t mark = drawTree.maskPoolMark();
    Rule drawn = prior.drawRuleForVariable(drawTree, store, rng, 0, 0);
    const std::uint64_t* words = drawTree.maskWordsFor(drawn);
    size_t numRight = maskPopcount(words, 2);
    neverEmpty &= numRight > 0 && numRight < K;
    for (size_t bit = 0; bit < K; ++bit)
      bitCounts[bit] += maskTestBit(words, static_cast<std::uint32_t>(bit))
        ? 1 : 0;
    drawTree.truncateMaskPool(mark);
  }
  check(neverEmpty, "pooled draws leave neither side empty");
  bool marginsMatch = true;
  double tolerance = 5.0 * std::sqrt(0.25 * numDraws);
  for (size_t bit = 0; bit < K; ++bit)
    marginsMatch &= std::fabs(bitCounts[bit] - 0.5 * numDraws) < tolerance;
  check(marginsMatch, "pooled draw direction bits are marginally uniform");
  checkNear(prior.ruleForVariableLogProbability(tree, store, 0),
            -(70.0 * std::log(2.0) + std::log1p(-std::pow(2.0, -69.0))),
            1e-13, "pooled rule probability uses the closed form");

  // compaction: strand garbage past the high-water mark, then verify the
  // live rule's words survive at a fresh offset and still partition
  std::vector<std::uint64_t> liveMask(tree.maskWordsFor(tree.at(0).rule),
                                      tree.maskWordsFor(tree.at(0).rule) + 2);
  for (int g = 0; g < 200; ++g) tree.allocateMask(2);
  tree.compactMaskPoolIfNeeded(store);
  check(tree.maskPool.size() == 2, "compaction keeps only live masks");
  check(maskEquals(tree.maskWordsFor(tree.at(0).rule), liveMask.data(), 2),
        "compaction preserves mask contents");
  tree.repartitionSubtree(store, 0);
  check(tree.at(tree.at(0).leftChild + 1).numObservations() == expectedRight,
        "compacted rules still partition");

  // the flattened side channel: offsets sequential in pre-order, category
  // bits only, gauge re-checked on rebuild
  std::vector<FlatNode> flat;
  std::vector<std::uint64_t> flatMasks;
  std::vector<double> paramByNode(tree.nodes.size(), 0.0);
  tree.flatten(store, paramByNode.data(), flat, nullptr, 1, nullptr,
               &flatMasks);
  check(flat.size() == 3 && flatMasks.size() == 2 &&
          flatKindOf(flat[0]) == FlatKind::categoricalPooled &&
          flat[0].maskOffset == 0 && flat[0].numMaskWords == 2,
        "pooled flatten emits the side channel");
  check(maskEquals(flatMasks.data(), liveMask.data(), 2),
        "flattened words match the live mask");
  check(flatTreeIsWellFormed(store, flat.data(), flat.size(),
                             flatMasks.data(), flatMasks.size()),
        "pooled flat tree is well formed");
  check(!flatTreeIsWellFormed(store, flat.data(), flat.size()),
        "a missing mask channel is rejected");

  std::vector<size_t> rebuiltIndices(n);
  std::vector<double> rebuiltParams;
  Tree rebuilt;
  rebuilt.initialize(rebuiltIndices.data(), n);
  check(rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams,
                              1, nullptr, flatMasks.data(), flatMasks.size()),
        "pooled tree rebuilds from flat");
  check(maskEquals(rebuilt.maskWordsFor(rebuilt.at(0).rule),
                   liveMask.data(), 2),
        "pooled mask round-trips exactly");

  flat[0].maskOffset = 1;  // offsets must be the running cursor
  rebuilt.initialize(rebuiltIndices.data(), n);
  check(!rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams,
                               1, nullptr, flatMasks.data(),
                               flatMasks.size()),
        "a non-sequential mask offset is rejected");
  flat[0].maskOffset = 0;
  std::vector<std::uint64_t> badMasks(flatMasks);
  maskSetBit(badMasks.data(), 70);  // the missing position must stay clear
  rebuilt.initialize(rebuiltIndices.data(), n);
  check(!rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams,
                               1, nullptr, badMasks.data(), badMasks.size()),
        "mask bits past the categories are rejected");

  // missing values compose: the NA pseudo-category is bit K
  std::vector<double> xMissing(x);
  xMissing[0] = std::nan("");
  xMissing[7] = std::nan("");
  ColumnStore storeMissing;
  storeMissing.build(xMissing.data(), n, 1, 10, false, types);
  check(storeMissing.numCuts[0] == K && storeMissing.hasMissing[0] == 1 &&
          storeMissing.train.codes[0] == 70,
        "pooled NA takes code K");
  Tree missingTree;
  missingTree.initialize(indices.data(), n);
  missingTree.reachableCategoriesWide(storeMissing, 0, 0, reachable);
  check(maskPopcount(reachable, 2) == K + 1 && maskTestBit(reachable, 70),
        "the missing pseudo-category enters the reachable set");

  printf("ok: pooled mask mechanics\n");
}

static void testMissingMechanics() {
  const size_t n = 90;
  double na = std::nan("");
  std::vector<double> x(n * 2), y(n);
  size_t numMissingOrdinal = 0, numMissingCategorical = 0;
  for (size_t i = 0; i < n; ++i) {
    if (i % 3 == 0) { x[i] = na; ++numMissingOrdinal; }
    else x[i] = runif01();
    if (i % 5 == 0) { x[i + n] = na; ++numMissingCategorical; }
    else x[i + n] = static_cast<double>(i % 4);
    y[i] = runif01();
  }
  ColumnType types[] = {ColumnType::ordinal, ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 2, 10, false, types);

  // an ordinal rule routes the reserved code by its missing direction
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(4);
  rule.setMissingGoesRight(true);
  check(rule.splitIndex() == 4 && rule.missingGoesRight(),
        "rule accessors pack the split index and missing direction");
  tree.birth(store, 0, rule, y.data(), nullptr);

  const Node& right(tree.at(tree.at(0).leftChild + 1));
  size_t missingOnRight = 0;
  for (size_t i = right.begin; i < right.end; ++i)
    if (store.train.codes[tree.indices[i]] == naCode) ++missingOnRight;
  const Node& left(tree.at(tree.at(0).leftChild));
  size_t missingOnLeft = 0;
  for (size_t i = left.begin; i < left.end; ++i)
    if (store.train.codes[tree.indices[i]] == naCode) ++missingOnLeft;
  check(missingOnRight == numMissingOrdinal && missingOnLeft == 0,
        "missing ordinal rows all route by the rule's direction");

  std::vector<size_t> indices2(n);
  Tree treeLeft;
  treeLeft.initialize(indices2.data(), n);
  Rule ruleLeft;
  ruleLeft.variableIndex = 0;
  ruleLeft.setSplitIndex(4);
  treeLeft.birth(store, 0, ruleLeft, y.data(), nullptr);
  const Node& leftLeft(treeLeft.at(treeLeft.at(0).leftChild));
  size_t missingOnLeft2 = 0;
  for (size_t i = leftLeft.begin; i < leftLeft.end; ++i)
    if (store.train.codes[treeLeft.indices[i]] == naCode) ++missingOnLeft2;
  check(missingOnLeft2 == numMissingOrdinal,
        "the canonical zero direction sends missing rows left");

  // a categorical rule treats missing as one more category: the reachable
  // mask seeds it and a rule can isolate it
  std::vector<size_t> indices3(n);
  Tree catTree;
  catTree.initialize(indices3.data(), n);
  check(catTree.reachableCategories(store, 0, 1) ==
          (0xfull | Rule::missingDirectionBit),
        "the reachable mask includes the missing category");
  Rule catRule;
  catRule.variableIndex = 1;
  catRule.setCategoryDirections(Rule::missingDirectionBit);
  catTree.birth(store, 0, catRule, y.data(), nullptr);
  const Node& catRight(catTree.at(catTree.at(0).leftChild + 1));
  check(catRight.numObservations() == numMissingCategorical,
        "a mask can send exactly the missing rows down one side");
  check(catTree.reachableCategories(store, catTree.at(0).leftChild, 1) ==
          0xfull &&
        catTree.reachableCategories(store, catTree.at(0).leftChild + 1, 1) ==
          Rule::missingDirectionBit,
        "children inherit the filtered missing category");

  // the flat format moves the missing direction out of the value
  std::vector<double> params(tree.nodes.size(), 0.0);
  std::vector<FlatNode> flat;
  tree.flatten(store, params.data(), flat);
  check((flat[0].flags & flatMissingGoesRight) != 0 &&
          flatKindOf(flat[0]) == FlatKind::ordinal &&
          flat[0].value == store.cutPoints[0][4],
        "an ordinal flat node carries its missing direction in flags");

  std::vector<double> catParams(catTree.nodes.size(), 0.0);
  std::vector<FlatNode> catFlat;
  catTree.flatten(store, catParams.data(), catFlat);
  check(flatKindOf(catFlat[0]) == FlatKind::categoricalInline &&
          catFlat[0].mask == 0 &&
          (catFlat[0].flags & flatMissingGoesRight) != 0,
        "a missing-only mask flattens to an empty mask plus the flag");
  check(flatTreeIsWellFormed(store, catFlat.data(), catFlat.size()),
        "the missing-only mask is well formed");

  std::vector<size_t> indices4(n);
  Tree rebuilt;
  rebuilt.initialize(indices4.data(), n);
  std::vector<double> rebuiltParams;
  check(rebuilt.buildFromFlat(store, catFlat.data(), catFlat.size(),
                              rebuiltParams),
        "a flat tree with a missing direction rebuilds");
  check(rebuilt.at(0).rule.categoryDirections() == Rule::missingDirectionBit,
        "the rebuilt rule recovers the missing direction");

  // raw-value replay routes NaN by the flag
  std::vector<std::uint32_t> counts(flat.size());
  std::vector<size_t> replayIndices(n);
  for (size_t i = 0; i < n; ++i) replayIndices[i] = i;
  countFlatObservationsBelow(flat.data(), x.data(), n,
                             replayIndices.data(), 0, n, counts.data());
  size_t expectedRight = 0;
  for (size_t i = 0; i < n; ++i)
    if (isNA(x[i]) || x[i] > store.cutPoints[0][4]) ++expectedRight;
  check(counts[2] == expectedRight,
        "replay against raw values routes NaN by the stored direction");

  FlatNode bad = catFlat[0];
  bad.flags |= 0x8u;  // a bit outside the missing flag and the kind field
  std::vector<FlatNode> badTree(catFlat);
  badTree[0] = bad;
  check(!flatTreeIsWellFormed(store, badTree.data(), badTree.size()),
        "unknown flag bits are malformed");
  std::vector<FlatNode> badLeaf(catFlat);
  badLeaf[1].flags = flatMissingGoesRight;
  check(!flatTreeIsWellFormed(store, badLeaf.data(), badLeaf.size()),
        "a flagged leaf is malformed");
  printf("ok: missing mechanics\n");
}

void runTreeTests(ext_rng* rng) {
  testTreeMechanics();
  testTreePriorMath();
  testCategoricalMechanics();
  testCategoricalPriorMath(rng);
  testPooledMaskMechanics(rng);
  testMissingMechanics();
}
