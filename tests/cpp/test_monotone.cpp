#include "common.hpp"

#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <set>

// The monotone leaf-order counter: log Z_T, the birth/death ratio counted on
// the finer tree's side, the position laws and the linear-extension draw,
// against brute-force permutation counts through monotoneTreeIsFeasible, the
// geometry the counter must reproduce. Takes no shared rng and restores the
// runif01 stream.

namespace {

struct TestTree {
  std::vector<index_t> index;
  Tree tree;
  const ColumnStore* store;
  std::vector<double> y;

  explicit TestTree(const ColumnStore& s)
      : index(s.numObservations), store(&s), y(s.numObservations, 0.0) {
    tree.initialize(index.data(), s.numObservations);
    tree.computeLeafStats(0, y.data(), nullptr);
  }
  std::int32_t split(std::int32_t node, std::int32_t var, std::int32_t cut) {
    Rule rule;
    rule.variableIndex = var;
    rule.setSplitIndex(cut);
    tree.birth(*store, node, rule, y.data(), nullptr);
    return tree.at(node).leftChild;
  }
  // a chain along var over [from, to] cuts, starting at leaf node; returns
  // the last (highest-code) leaf
  std::int32_t chain(std::int32_t node, std::int32_t var, std::int32_t from,
                     std::int32_t to) {
    for (std::int32_t s = from; s <= to; ++s) node = split(node, var, s) + 1;
    return node;
  }
  void growRandom(int numBirths) {
    for (int step = 0; step < numBirths; ++step) {
      std::vector<std::int32_t> leaves;
      tree.fillBottom(0, leaves);
      std::int32_t leaf =
          leaves[static_cast<size_t>(runif01() * leaves.size())];
      int var = static_cast<int>(runif01() * store->numPredictors);
      if (store->splitsBySubset(static_cast<size_t>(var))) {
        splitSubset(leaf, var);
        continue;
      }
      std::int32_t left, right;
      tree.splitInterval(*store, leaf, var, &left, &right);
      if (right < left) continue;
      Rule rule;
      rule.variableIndex = var;
      rule.setSplitIndex(
          left + static_cast<std::int32_t>(runif01() * (right - left + 1)));
      if (store->hasMissing[static_cast<size_t>(var)])
        rule.setMissingGoesRight(runif01() < 0.5);
      tree.birth(*store, leaf, rule, y.data(), nullptr);
    }
  }
  // a random nontrivial split of the positions reaching the leaf, missing
  // included when the column has missing values, as the rule draw makes one
  void splitSubset(std::int32_t leaf, int var) {
    size_t j = static_cast<size_t>(var);
    size_t numWords = maskWordsForCount(store->categoryCounts[j]);
    std::vector<std::uint64_t> reachable(numWords);
    if (numWords > 1)
      tree.reachableCategoriesWide(*store, leaf, var, reachable.data());
    else
      reachable[0] = tree.reachableCategories(*store, leaf, var);
    size_t count = maskPopcount(reachable.data(), numWords);
    if (count < 2) return;
    std::vector<std::uint64_t> right(numWords);
    do {
      std::fill(right.begin(), right.end(), 0);
      for (std::uint32_t b = 0; b < 64 * numWords; ++b)
        if (maskTestBit(reachable.data(), b) && runif01() < 0.5)
          maskSetBit(right.data(), b);
    } while (maskPopcount(right.data(), numWords) % count == 0);
    Rule rule;
    rule.variableIndex = var;
    if (numWords > 1) {
      size_t offset = tree.allocateMask(numWords);
      std::copy(right.begin(), right.end(), tree.mutableMaskWordsFor(offset));
      rule.setMaskOffset(offset);
    } else {
      rule.setCategoryDirections(right[0]);
    }
    tree.birth(*store, leaf, rule, y.data(), nullptr);
  }
};

void makeStore(ColumnStore& store, size_t p, std::uint32_t cuts, size_t n) {
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();
  built(store.build(x.data(), n, p, cuts));
}

// x0 numeric; x1 numeric, 30% missing; x2 a 4-level factor, 20% missing; x3
// a 70-level factor (pooled), 10% missing
void makeMixedStore(ColumnStore& store, size_t n) {
  const double na = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(n * 4);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[n + i] = runif01() < 0.3 ? na : runif01();
    x[2 * n + i] = runif01() < 0.2 ? na : std::floor(4.0 * runif01());
    x[3 * n + i] = runif01() < 0.1 ? na : static_cast<double>(i % 70);
  }
  const ColumnKind types[4] = {ColumnKind::numeric, ColumnKind::numeric,
                               ColumnKind::categorical,
                               ColumnKind::categorical};
  built(store.build(x.data(), n, 4, 6, false, types));
}

// The order the constraint requires, from points, independently of the
// engine's geometry: every combination of the tree's split variables' values
// (each code or level, and missing where training had one, as predict takes
// it) routed as prediction routes it. Two points differing only in a constrained
// predictor, by one observed code, require their leaves ordered. Returns the
// (lower, upper) leaf pairs; withMissing false leaves missing values out.
using LeafPairs = std::set<std::pair<std::int32_t, std::int32_t>>;
LeafPairs pointOrder(const Tree& tree, const ColumnStore& store,
                     const std::int8_t* dir, bool withMissing = true) {
  std::vector<std::int32_t> internal, axes;
  tree.fillNotBottom(0, internal);
  for (std::int32_t node : internal) {
    std::int32_t v = tree.at(node).rule.variableIndex;
    if (std::find(axes.begin(), axes.end(), v) == axes.end()) axes.push_back(v);
  }
  // per axis its values, missing last, and the stride of its digit
  std::vector<std::vector<xint_t>> values(axes.size());
  std::vector<size_t> stride(axes.size());
  size_t numPoints = 1;
  for (size_t a = 0; a < axes.size(); ++a) {
    size_t j = static_cast<size_t>(axes[a]);
    bool subset = store.splitsBySubset(j);
    std::uint32_t count = subset ? store.categoryCounts[j] : store.numCuts[j] + 1;
    for (std::uint32_t c = 0; c < count; ++c)
      values[a].push_back(static_cast<xint_t>(c));
    if (withMissing && store.hasMissing[j])
      values[a].push_back(subset ? missingCategoryCode(store.categoryCounts[j])
                                 : naCode);
    stride[a] = numPoints;
    numPoints *= values[a].size();
  }
  std::vector<std::int32_t> leafOf(numPoints);
  std::vector<xint_t> code(store.numPredictors, 0);
  for (size_t point = 0; point < numPoints; ++point) {
    for (size_t a = 0; a < axes.size(); ++a)
      code[static_cast<size_t>(axes[a])] =
          values[a][(point / stride[a]) % values[a].size()];
    std::int32_t node = 0;
    while (!tree.at(node).isBottom()) {
      const Rule& rule = tree.at(node).rule;
      bool right = tree.ruleSendsRight(
          store, rule, code[static_cast<size_t>(rule.variableIndex)]);
      node = tree.at(node).leftChild + (right ? 1 : 0);
    }
    leafOf[point] = node;
  }
  LeafPairs out;
  for (size_t a = 0; a < axes.size(); ++a) {
    std::int8_t d = dir[axes[a]];
    if (d == 0 || store.splitsBySubset(static_cast<size_t>(axes[a]))) continue;
    // observed codes c and c + 1: the missing digit, if any, is last
    size_t last = values[a].size() -
                  (withMissing && store.hasMissing[static_cast<size_t>(axes[a])]
                       ? 2
                       : 1);
    for (size_t point = 0; point < numPoints; ++point) {
      if ((point / stride[a]) % values[a].size() >= last) continue;
      std::int32_t low = leafOf[point], high = leafOf[point + stride[a]];
      if (low == high) continue;
      out.insert(d > 0 ? std::make_pair(low, high) : std::make_pair(high, low));
    }
  }
  return out;
}

size_t numLeaves(const Tree& tree) {
  std::vector<std::int32_t> leaves;
  tree.fillBottom(0, leaves);
  return leaves.size();
}

// Every ordering of the leaves by distinct values, kept when the geometry
// admits it: the whole tree's extensions, and among them the share where c2's
// rank among the leaves in `group` is c1's plus one (restricted to a union of
// components, a uniform extension of the tree is one of the union).
struct Brute {
  double extensions = 0.0, adjacent = 0.0;
};
Brute bruteForce(const Tree& tree, const ColumnStore& store,
                 const std::int8_t* dir, std::int32_t c1 = invalidNode,
                 std::int32_t c2 = invalidNode,
                 const std::vector<std::int32_t>& group = {}) {
  std::vector<std::int32_t> leaves;
  tree.fillBottom(0, leaves);
  std::vector<int> perm(leaves.size());
  for (size_t i = 0; i < perm.size(); ++i) perm[i] = static_cast<int>(i);
  std::vector<double> mu(tree.nodes.size(), 0.0);
  Brute out;
  do {
    for (size_t i = 0; i < leaves.size(); ++i) mu[leaves[i]] = perm[i];
    if (!monotoneTreeIsFeasible(tree, store, dir, mu.data(), 0.0)) continue;
    out.extensions += 1.0;
    if (c1 == invalidNode) continue;
    int r1 = 0, r2 = 0;
    for (std::int32_t g : group) {
      r1 += mu[g] < mu[c1];
      r2 += mu[g] < mu[c2];
    }
    out.adjacent += r2 == r1 + 1;
  } while (std::next_permutation(perm.begin(), perm.end()));
  return out;
}

double logFactorial(size_t n) {
  return std::lgamma(static_cast<double>(n) + 1.0);
}
double logChoose(size_t n, size_t k) {
  return logFactorial(n) - logFactorial(k) - logFactorial(n - k);
}

// log Z_T0 - log Z_T* by two whole-tree counts, T0 the pair merged
double directLogRatio(const Tree& tree, const ColumnStore& store,
                      const std::int8_t* dir, std::int32_t parent,
                      MonotoneCountScratch& s) {
  double logStar = monotoneLogNormalizer(tree, store, dir, s);
  Tree merged = tree;
  merged.orphanChildren(parent);
  return monotoneLogNormalizer(merged, store, dir, s) - logStar;
}

// c1 is the child lower in the order: the higher-code one on a decreasing
// axis
void pairOf(const Tree& tree, const std::int8_t* dir, std::int32_t parent,
            std::int32_t& c1, std::int32_t& c2) {
  std::int32_t left = tree.at(parent).leftChild;
  bool flip = dir[tree.at(parent).rule.variableIndex] < 0;
  c1 = flip ? left + 1 : left;
  c2 = flip ? left : left + 1;
}

// the leaves in the components of T* holding the pair
std::vector<std::int32_t> unionOf(const MonotoneLeafOrder& order,
                                  std::int32_t c1, std::int32_t c2) {
  std::vector<std::int32_t> out;
  for (size_t k : {order.componentOf[order.positionOf(c1)],
                   order.componentOf[order.positionOf(c2)]})
    for (std::int32_t leaf : order.components[k].leaves)
      if (std::find(out.begin(), out.end(), leaf) == out.end())
        out.push_back(leaf);
  return out;
}

// Whether a sequence of local labels is a linear extension of c.
bool isExtension(const MonotoneOrderComponent& c, const size_t* seq) {
  std::vector<size_t> rank(c.size);
  for (size_t r = 0; r < c.size; ++r) rank[seq[r]] = r;
  for (size_t x = 0; x < c.size; ++x)
    for (size_t y = 0; y < c.size; ++y)
      if ((c.predecessorsOf(x)[y / 64] >> (y % 64)) & 1u)
        if (rank[y] > rank[x]) return false;
  return true;
}

double seconds(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}

}  // namespace

// The geometry against the point oracle on random trees over a numeric store
// and the mixed one (factor splits, pooled ones, and missing values routed
// either way, in free and constrained predictors): the builder's relation is
// exactly the required pairs, each leaf's bounds read them, and the
// feasibility check agrees on random and on feasible leaf values.
static void testMonotoneGeometryPoints() {
  ColumnStore numeric, mixed;
  makeStore(numeric, 3, 24, 600);
  makeMixedStore(mixed, 700);
  const std::int8_t numericDirs[3][3] = {{1, 0, 0}, {1, -1, 0}, {1, -1, 1}};
  const std::int8_t mixedDirs[4][4] = {
      {1, 0, 0, 0}, {-1, 0, 0, 0}, {1, -1, 0, 0}, {0, 1, 0, 0}};
  MonotoneLeafOrder order;
  MonotoneNeighborScratch scratch;
  int numTrees = 0, numPairs = 0, factorPairs = 0, missingPairs = 0;
  bool relationOk = true, boundsOk = true, feasibleOk = true;
  for (int which = 0; which < 2; ++which) {
    const ColumnStore& store = which == 0 ? numeric : mixed;
    for (int trial = 0; trial < (which == 0 ? 200 : 600); ++trial) {
      const std::int8_t* dir =
          which == 0 ? numericDirs[trial % 3] : mixedDirs[trial % 4];
      TestTree t(store);
      t.growRandom(4 + trial % 8);
      LeafPairs required = pointOrder(t.tree, store, dir);
      buildMonotoneLeafOrder(t.tree, store, dir, order);
      LeafPairs relation;
      for (const MonotoneOrderComponent& c : order.components)
        for (size_t x = 0; x < c.size; ++x)
          for (size_t y = 0; y < c.size; ++y)
            if ((c.predecessorsOf(x)[y / 64] >> (y % 64)) & 1u)
              relation.insert({c.leaves[y], c.leaves[x]});
      relationOk = relationOk && relation == required;
      ++numTrees;
      numPairs += static_cast<int>(required.size());
      if (which == 1) {
        LeafPairs observed = pointOrder(t.tree, store, dir, false);
        for (const auto& pair : required)
          missingPairs += observed.count(pair) == 0;
        std::vector<std::int32_t> internal;
        t.tree.fillNotBottom(0, internal);
        for (std::int32_t node : internal)
          if (store.splitsBySubset(static_cast<size_t>(
                  t.tree.at(node).rule.variableIndex))) {
            factorPairs += static_cast<int>(required.size());
            break;
          }
      }

      std::vector<std::int32_t> leaves;
      t.tree.fillBottom(0, leaves);
      std::vector<double> mu(t.tree.nodes.size(), 0.0);
      for (std::int32_t leaf : leaves) mu[leaf] = runif01();
      for (std::int32_t k : leaves) {
        double a = -HUGE_VAL, b = HUGE_VAL, ea, eb;
        bool c = false, ec;
        for (const auto& pair : required) {
          if (pair.second == k) a = std::max(a, mu[pair.first]);
          if (pair.first == k) b = std::min(b, mu[pair.second]);
          c = c || pair.first == k || pair.second == k;
        }
        monotoneNeighborBounds(t.tree, store, dir, leaves, k, mu.data(),
                               nullptr, 0, scratch, &ea, &eb, &ec);
        boundsOk = boundsOk && ea == a && eb == b && ec == c;
      }
      bool feasible = true;
      for (const auto& pair : required)
        feasible = feasible && mu[pair.first] <= mu[pair.second];
      feasibleOk = feasibleOk && monotoneTreeIsFeasible(t.tree, store, dir,
                                                        mu.data(), 0.0) ==
                                     feasible;
      // each leaf at its longest chain of required predecessors is feasible
      for (std::int32_t leaf : leaves) mu[leaf] = 0.0;
      for (size_t pass = 0; pass < leaves.size(); ++pass)
        for (const auto& pair : required)
          mu[pair.second] = std::max(mu[pair.second], mu[pair.first] + 1.0);
      feasibleOk = feasibleOk &&
                   monotoneTreeIsFeasible(t.tree, store, dir, mu.data(), 0.0);
    }
  }
  // hand-built: a relation only the pooled factor's missing value carries.
  // x0 cut; below, x3's right side holds levels 0-34 and missing; above,
  // levels 35-69 and missing. Three pairs: two through levels, one through
  // missing alone.
  {
    const std::int8_t dir[4] = {1, 0, 0, 0};
    TestTree t(mixed);
    std::int32_t low = t.split(0, 0, 3);
    size_t numWords = maskWordsForCount(mixed.categoryCounts[3]);
    for (std::int32_t half : {low, low + 1}) {
      size_t offset = t.tree.allocateMask(numWords);
      std::uint64_t* words = t.tree.mutableMaskWordsFor(offset);
      for (std::uint32_t c = 0; c < 70; ++c)
        if ((c < 35) == (half == low)) maskSetBit(words, c);
      maskSetBit(words, static_cast<std::uint32_t>(missingCategoryCode(70)));
      Rule rule;
      rule.variableIndex = 3;
      rule.setMaskOffset(offset);
      t.tree.birth(mixed, half, rule, t.y.data(), nullptr);
    }
    LeafPairs required = pointOrder(t.tree, mixed, dir);
    buildMonotoneLeafOrder(t.tree, mixed, dir, order);
    LeafPairs relation;
    for (const MonotoneOrderComponent& c : order.components)
      for (size_t x = 0; x < c.size; ++x)
        for (size_t y = 0; y < c.size; ++y)
          if ((c.predecessorsOf(x)[y / 64] >> (y % 64)) & 1u)
            relation.insert({c.leaves[y], c.leaves[x]});
    check(required.size() == 3 && relation == required,
          "monotone geometry: a pooled factor's missing value relates leaves");
  }
  // a factor without missing values: the order does not depend on which side
  // of a level split is called right (one partition, labelled two ways)
  {
    ColumnStore plain;
    std::vector<double> x(2 * 400);
    for (size_t i = 0; i < 400; ++i) {
      x[i] = runif01();
      x[400 + i] = static_cast<double>(i % 4);
    }
    const ColumnKind types[2] = {ColumnKind::numeric, ColumnKind::categorical};
    built(plain.build(x.data(), 400, 2, 6, false, types));
    const std::int8_t dir[2] = {1, 0};
    MonotoneCountScratch count;
    double logZ[2];
    size_t numEdges[2];
    for (int label = 0; label < 2; ++label) {
      TestTree t(plain);
      std::int32_t low = t.split(0, 0, 3);
      for (std::int32_t half : {low, low + 1}) {
        Rule rule;
        rule.variableIndex = 1;
        rule.setCategoryDirections(half == low || label == 0 ? 0x3u : 0xcu);
        t.tree.birth(plain, half, rule, t.y.data(), nullptr);
      }
      logZ[label] = monotoneLogNormalizer(t.tree, plain, dir, count);
      numEdges[label] = pointOrder(t.tree, plain, dir).size();
    }
    check(logZ[0] == logZ[1] && numEdges[0] == 2 && numEdges[1] == 2,
          "monotone geometry: a level split's labels leave the order alone");
  }
  check(relationOk, "monotone geometry: relation equals the point oracle");
  check(boundsOk, "monotone geometry: bounds read the required pairs");
  check(feasibleOk, "monotone geometry: feasibility agrees with the oracle");
  check(factorPairs > 1000 && missingPairs > 100,
        "monotone geometry: factor and missing-value pairs covered");
  printf("ok: monotone geometry vs point oracle (%d trees, %d pairs, %d in "
         "factor-split trees, %d through missing values)\n",
         numTrees, numPairs, factorPairs, missingPairs);
}

// Hand-built orders with closed forms, each also against brute force when
// small: chains over 64 and 65 leaves (one and two words), a star, an N, an
// N plus a chain, and the 2x2x2 grid over three mixed-direction axes (the
// Boolean lattice B3, e = 48).
static void testMonotoneCountHandBuilt() {
  ColumnStore store;
  makeStore(store, 3, 80, 400);
  MonotoneCountScratch s;
  const std::int8_t inc[3] = {1, 0, 1}, dec[3] = {-1, 0, -1},
                    mixed[3] = {1, -1, 1};

  for (size_t length : {64u, 65u}) {
    for (const std::int8_t* dir : {inc, dec}) {
      TestTree t(store);
      t.chain(0, 0, 0, static_cast<std::int32_t>(length) - 2);
      check(numLeaves(t.tree) == length, "monotone count: chain built");
      checkNear(monotoneLogNormalizer(t.tree, store, dir, s),
                -logFactorial(length), 1e-9, "monotone count: chain e = 1");
    }
  }

  {  // a minimum below 12 strips cut on the free axis: e = 12!
    TestTree t(store);
    std::int32_t right = t.split(0, 0, 0) + 1;
    t.chain(right, 1, 0, 10);
    checkNear(monotoneLogNormalizer(t.tree, store, inc, s),
              logFactorial(12) - logFactorial(13), 1e-12,
              "monotone count: star e = 12!");
    const MonotoneOrderComponent& star = s.order.components[0];
    check(star.size == 13, "monotone count: star is one component");
    monotoneLogExtensions(star, s, false);
    const MonotoneDownSetLayer& whole = s.layers[star.size % 2];
    check(std::ldexp(whole.forward[0],
                     static_cast<int>(whole.forwardExponent)) == 479001600.0,
          "monotone count: star count exact");
  }

  {  // N: x1 cut, the halves cut on x2 at different codes; then a chain on top
    for (const std::int8_t* dir : {inc, dec}) {
      TestTree t(store);
      std::int32_t left = t.split(0, 0, 40);
      t.split(left, 1, 30);
      t.split(left + 1, 1, 50);
      checkNear(monotoneLogNormalizer(t.tree, store, dir, s),
                std::log(5.0 / 24.0), 1e-12, "monotone count: N e = 5");
      checkNear(bruteForce(t.tree, store, dir).extensions, 5.0, 0.0,
                "monotone count: N brute force");
      t.split(t.tree.at(left + 1).leftChild + 1, 0, 60);
      Brute brute = bruteForce(t.tree, store, dir);
      checkNear(monotoneLogNormalizer(t.tree, store, dir, s),
                std::log(brute.extensions) - logFactorial(5), 1e-12,
                "monotone count: N plus chain vs brute force");
    }
  }

  {  // 2x2x2 over three axes, directions mixed: B3
    TestTree t(store);
    std::int32_t a = t.split(0, 0, 40);
    for (std::int32_t half : {a, a + 1}) {
      std::int32_t b = t.split(half, 1, 40);
      for (std::int32_t quarter : {b, b + 1}) t.split(quarter, 2, 40);
    }
    checkNear(monotoneLogNormalizer(t.tree, store, mixed, s),
              std::log(48.0) - logFactorial(8), 1e-12,
              "monotone count: 2x2x2 e = 48");
    checkNear(bruteForce(t.tree, store, mixed).extensions, 48.0, 0.0,
              "monotone count: 2x2x2 brute force");
  }
  printf("ok: monotone count, hand-built orders\n");
}

using Directions = std::vector<std::vector<std::int8_t>>;

// 200 random trees over 1-3 constrained axes with mixed directions: log Z_T
// against the brute-force e / L!.
static void testMonotoneCountRandom(const ColumnStore& store,
                                    const Directions& dirs,
                                    const char* label) {
  MonotoneCountScratch s;
  double worst = 0.0;
  size_t largest = 0;
  for (int trial = 0; trial < 200; ++trial) {
    TestTree t(store);
    t.growRandom(2 + trial % 5);
    const std::int8_t* dir = dirs[trial % dirs.size()].data();
    double logZ = monotoneLogNormalizer(t.tree, store, dir, s);
    size_t L = numLeaves(t.tree);
    double expected =
        std::log(bruteForce(t.tree, store, dir).extensions) - logFactorial(L);
    worst = std::max(worst, std::fabs(logZ - expected));
    for (const MonotoneOrderComponent& c : s.order.components)
      largest = std::max(largest, c.size);
  }
  checkNear(worst, 0.0, 1e-12, "monotone count: random trees vs brute force");
  check(largest >= 5, "monotone count: random trees reach 5-leaf components");
  printf("ok: monotone count, 200 random trees%s (worst |log Z err| %.2g)\n",
         label, worst);
}

// The ratio on T*'s side against two direct whole-tree counts over every
// death (both children leaves) of random trees, one component and two, and
// against brute-force theta on the smaller ones.
static void testMonotoneRatioRandom(const ColumnStore& store,
                                    const Directions& dirs,
                                    const char* label) {
  MonotoneCountScratch s, direct;
  int one = 0, two = 0, decreasingPairs = 0, bruteChecked = 0;
  double worst = 0.0, worstTheta = 0.0;
  for (int trial = 0; trial < 200; ++trial) {
    TestTree t(store);
    t.growRandom(3 + trial % 9);
    const std::int8_t* dir = dirs[trial % dirs.size()].data();
    std::vector<std::int32_t> internal;
    t.tree.fillNotBottom(0, internal);
    for (std::int32_t parent : internal) {
      if (!t.tree.childrenAreBottom(parent)) continue;
      double ratio = monotoneLogNormalizerRatio(t.tree, store, dir, parent, s);
      std::int32_t c1, c2;
      pairOf(t.tree, dir, parent, c1, c2);
      bool same = s.order.componentOf[s.order.positionOf(c1)] ==
                  s.order.componentOf[s.order.positionOf(c2)];
      std::vector<std::int32_t> group = unionOf(s.order, c1, c2);
      (same ? one : two)++;
      decreasingPairs += dir[t.tree.at(parent).rule.variableIndex] < 0;
      worst =
          std::max(worst, std::fabs(ratio - directLogRatio(t.tree, store, dir,
                                                           parent, direct)));
      if (numLeaves(t.tree) <= 7) {
        Brute brute = bruteForce(t.tree, store, dir, c1, c2, group);
        double theta = std::exp(ratio) / static_cast<double>(group.size());
        worstTheta = std::max(
            worstTheta, std::fabs(theta - brute.adjacent / brute.extensions));
        ++bruteChecked;
      }
    }
  }
  checkNear(worst, 0.0, 1e-12, "monotone ratio: T* side vs direct counts");
  checkNear(worstTheta, 0.0, 1e-12, "monotone ratio: theta vs brute force");
  check(one > 100 && two > 20 && decreasingPairs > 20 && bruteChecked > 100,
        "monotone ratio: one- and two-component and decreasing pairs covered");
  printf("ok: monotone ratio%s, %d one-component and %d two-component moves "
         "(%d decreasing, %d brute-forced)\n",
         label, one, two, decreasingPairs, bruteChecked);
}

// The plan's named examples, x1 constrained (both directions) and x2 free.
static void testMonotoneRatioNamed() {
  ColumnStore store;
  makeStore(store, 2, 40, 300);
  MonotoneCountScratch s;
  for (std::int8_t d : {1, -1}) {
    const std::int8_t dir[2] = {d, 0};
    {  // x1 cut twice into a < b < c, b cut on x2: one component
      TestTree t(store);
      std::int32_t b = t.split(t.split(0, 0, 10) + 1, 0, 20);
      std::int32_t pair = t.split(b, 1, 20);
      (void)pair;
      double ratio = monotoneLogNormalizerRatio(t.tree, store, dir, b, s);
      checkNear(std::exp(ratio), 2.0, 1e-12, "monotone ratio: middle of three");
      std::int32_t c1, c2;
      pairOf(t.tree, dir, b, c1, c2);
      Brute brute =
          bruteForce(t.tree, store, dir, c1, c2, unionOf(s.order, c1, c2));
      checkNear(brute.adjacent / brute.extensions, 0.5, 1e-15,
                "monotone ratio: middle of three brute theta");
    }
    {  // x1 cut once, both halves cut on x2 at one code: a sibling pair in two
      TestTree t(store);
      std::int32_t left = t.split(0, 0, 20);
      t.split(left, 1, 20);
      t.split(left + 1, 1, 20);
      for (std::int32_t half : {left, left + 1}) {
        double ratio = monotoneLogNormalizerRatio(t.tree, store, dir, half, s);
        checkNear(std::exp(ratio), 4.0 / 3.0, 1e-12, "monotone ratio: 2x2");
        std::int32_t c1, c2;
        pairOf(t.tree, dir, half, c1, c2);
        check(s.order.componentOf[s.order.positionOf(c1)] !=
                  s.order.componentOf[s.order.positionOf(c2)],
              "monotone ratio: 2x2 pair in two components");
        Brute brute =
            bruteForce(t.tree, store, dir, c1, c2, unionOf(s.order, c1, c2));
        checkNear(brute.adjacent / brute.extensions, 1.0 / 3.0, 1e-15,
                  "monotone ratio: 2x2 brute theta");
      }
    }
    {  // a single x1 cut: theta 1 with c1 lower in the order, 0 by code
      TestTree t(store);
      std::int32_t left = t.split(0, 0, 20);
      double ratio = monotoneLogNormalizerRatio(t.tree, store, dir, 0, s);
      checkNear(std::exp(ratio), 2.0, 1e-12, "monotone ratio: single cut");
      std::int32_t c1, c2;
      pairOf(t.tree, dir, 0, c1, c2);
      std::vector<std::int32_t> both = {left, left + 1};
      checkNear(bruteForce(t.tree, store, dir, c1, c2, both).adjacent, 1.0, 0.0,
                "monotone ratio: single cut, c2 follows c1 by order");
      checkNear(bruteForce(t.tree, store, dir, left, left + 1, both).adjacent,
                d > 0 ? 1.0 : 0.0, 0.0,
                "monotone ratio: single cut, by code only when increasing");
    }
  }
  printf("ok: monotone ratio, named examples\n");
}

// Position laws against enumerated extensions, every element of the
// components of random trees.
static void testMonotonePositionLaw() {
  ColumnStore store;
  makeStore(store, 3, 12, 300);
  MonotoneCountScratch s;
  const std::int8_t dir[3] = {1, -1, 0};
  double worst = 0.0;
  int checked = 0;
  for (int trial = 0; trial < 60; ++trial) {
    TestTree t(store);
    t.growRandom(4 + trial % 6);
    buildMonotoneLeafOrder(t.tree, store, dir, s.order);
    for (const MonotoneOrderComponent& c : s.order.components) {
      if (c.size < 2 || c.size > 8) continue;
      std::vector<size_t> perm(c.size);
      for (size_t i = 0; i < c.size; ++i) perm[i] = i;
      std::vector<double> counts(c.size * c.size, 0.0);
      double total = 0.0;
      do {
        if (!isExtension(c, perm.data())) continue;
        total += 1.0;
        for (size_t r = 0; r < c.size; ++r) counts[perm[r] * c.size + r] += 1.0;
      } while (std::next_permutation(perm.begin(), perm.end()));
      double logE = monotoneLogExtensions(c, s, true);
      checkNear(logE, std::log(total), 1e-12, "monotone position law: count");
      monotoneBackwardCounts(c, s);
      for (size_t x = 0; x < c.size; ++x) {
        monotoneLogPositionLaw(c, s, x, logE, s.firstLaw);
        for (size_t r = 0; r < c.size; ++r)
          worst = std::max(worst, std::fabs(std::exp(s.firstLaw[r]) -
                                            counts[x * c.size + r] / total));
        ++checked;
      }
    }
  }
  checkNear(worst, 0.0, 1e-13, "monotone position law vs enumeration");
  check(checked > 100, "monotone position law: elements checked");
  printf("ok: monotone position law (%d elements)\n", checked);
}

// The linear-extension draw is uniform: chi-square over the extensions of a
// 5-leaf N plus a chain; every draw of a two-word component is an extension;
// and the exact prior leaf draw lands in the cone with the rejection
// sampler's per-leaf means.
static void testMonotoneExtensionDraw() {
  ColumnStore store;
  makeStore(store, 2, 80, 400);
  MonotoneCountScratch s;
  const std::int8_t dir[2] = {1, 0};
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20260930u);

  // x2 above 70 is one isolated leaf; below it, the N plus a chain
  TestTree t(store);
  std::int32_t left = t.split(t.split(0, 1, 70), 0, 40);
  t.split(left, 1, 30);
  t.split(left + 1, 1, 50);
  t.split(t.tree.at(left + 1).leftChild + 1, 0, 60);
  buildMonotoneLeafOrder(t.tree, store, dir, s.order);
  check(s.order.components.size() == 2 && s.order.components[0].size == 5,
        "monotone draw: N plus chain, and an isolated leaf");
  const MonotoneOrderComponent c = s.order.components[0];
  double e = std::exp(monotoneLogExtensions(c, s, true));
  const int numDraws = 40000;
  std::map<std::vector<size_t>, int> seen;
  std::vector<size_t> seq(c.size);
  bool valid = true;
  for (int i = 0; i < numDraws; ++i) {
    monotoneDrawExtension(rng, c, s, seq.data());
    valid = valid && isExtension(c, seq.data());
    ++seen[seq];
  }
  check(valid, "monotone draw: every draw an extension");
  e = std::round(e);
  checkNear(static_cast<double>(seen.size()), e, 0.0,
            "monotone draw: every extension drawn");
  double expected = numDraws / e, chiSquare = 0.0;
  for (const auto& kv : seen)
    chiSquare += (kv.second - expected) * (kv.second - expected) / expected;
  // Wilson-Hilferty upper 1e-4 point
  double df = e - 1.0, h = 2.0 / (9.0 * df);
  double critical = df * std::pow(1.0 - h + 3.719 * std::sqrt(h), 3.0);
  check(chiSquare < critical, "monotone draw: uniform over extensions");
  printf(
      "ok: monotone extension draw uniform (e %.0f, chi-square %.1f < %.1f)\n",
      e, chiSquare, critical);

  // two chains of 32 on a shared minimum: 65 elements, two words
  TestTree wide(store);
  std::int32_t right = wide.split(0, 0, 0) + 1;
  std::int32_t halves = wide.split(right, 1, 40);
  wide.chain(halves, 0, 1, 31);
  wide.chain(halves + 1, 0, 1, 31);
  buildMonotoneLeafOrder(wide.tree, store, dir, s.order);
  const MonotoneOrderComponent w = s.order.components[0];
  check(w.size == 65 && w.words == 2, "monotone draw: two-word component");
  checkNear(monotoneLogExtensions(w, s, true), logChoose(64, 32), 1e-9,
            "monotone draw: two chains on a minimum count");
  seq.resize(w.size);
  valid = true;
  for (int i = 0; i < 200; ++i) {
    monotoneDrawExtension(rng, w, s, seq.data());
    valid = valid && isExtension(w, seq.data());
  }
  check(valid, "monotone draw: two-word draws are extensions");

  // the exact prior leaf draw against rejection on the same tree
  std::vector<std::int32_t> leaves;
  t.tree.fillBottom(0, leaves);
  buildMonotoneLeafOrder(t.tree, store, dir, s.order);
  std::vector<bool> related(t.tree.nodes.size(), false);
  for (const MonotoneOrderComponent& comp : s.order.components)
    for (std::int32_t leaf : comp.leaves) related[leaf] = comp.size > 1;
  const double constrainedSd = 1.3, freeSd = 0.7;
  std::vector<double> mu(t.tree.nodes.size(), 0.0), sumExact(mu.size(), 0.0),
      sumReject(mu.size(), 0.0), sqExact(mu.size(), 0.0),
      sqReject(mu.size(), 0.0);
  bool feasible = true;
  const int priorDraws = 20000;
  for (int i = 0; i < priorDraws; ++i) {
    monotoneDrawPriorLeaves(rng, t.tree, store, dir, constrainedSd, freeSd, s,
                            mu.data());
    feasible =
        feasible && monotoneTreeIsFeasible(t.tree, store, dir, mu.data());
    for (std::int32_t leaf : leaves) {
      sumExact[leaf] += mu[leaf];
      sqExact[leaf] += mu[leaf] * mu[leaf];
    }
    do {
      for (std::int32_t leaf : leaves)
        mu[leaf] = (related[leaf] ? constrainedSd : freeSd) *
                   ext_rng_simulateStandardNormal(rng);
    } while (!monotoneTreeIsFeasible(t.tree, store, dir, mu.data()));
    for (std::int32_t leaf : leaves) {
      sumReject[leaf] += mu[leaf];
      sqReject[leaf] += mu[leaf] * mu[leaf];
    }
  }
  check(feasible, "monotone prior draw: every draw in the cone");
  double worstZ = 0.0;
  for (std::int32_t leaf : leaves) {
    double m1 = sumExact[leaf] / priorDraws, m2 = sumReject[leaf] / priorDraws;
    double v1 = sqExact[leaf] / priorDraws - m1 * m1,
           v2 = sqReject[leaf] / priorDraws - m2 * m2;
    worstZ = std::max(worstZ,
                      std::fabs(m1 - m2) / std::sqrt((v1 + v2) / priorDraws));
  }
  check(worstZ < 4.5, "monotone prior draw: means match rejection");
  printf("ok: monotone exact prior draw (%zu leaves, worst mean |z| %.2f)\n",
         leaves.size(), worstZ);
  ext_rng_destroy(rng);
}

// Past the double range: two 515-leaf chains on a shared minimum (e =
// C(1030, 515) > 1.8e308), and a two-component pair atop two 515-leaf chains
// (theta = ab / ((a+b)(a+b-1)), with C(a+b, a) past the range).
static void testMonotoneCountScale() {
  ColumnStore store;
  makeStore(store, 2, 600, 3000);
  check(store.numCuts[0] >= 520, "monotone scale: enough cuts");
  MonotoneCountScratch s, direct;
  const std::int8_t dir[2] = {1, 0};

  TestTree t(store);
  std::int32_t halves = t.split(t.split(0, 0, 0) + 1, 1, 300);
  t.chain(halves, 0, 1, 514);
  t.chain(halves + 1, 0, 1, 514);
  auto t0 = std::chrono::steady_clock::now();
  double logZ = monotoneLogNormalizer(t.tree, store, dir, s);
  double elapsed = seconds(t0);
  check(s.order.components.size() == 1 && s.order.components[0].size == 1031,
        "monotone scale: one 1031-leaf component");
  checkNear(logZ, logChoose(1030, 515) - logFactorial(1031), 1e-9,
            "monotone scale: log e of two chains on a minimum");
  check(logChoose(1030, 515) > std::log(DBL_MAX),
        "monotone scale: past doubles");

  TestTree u(store);
  std::int32_t left = u.split(0, 0, 515);
  std::int32_t lower = u.split(left, 1, 300);
  u.chain(lower, 0, 0, 512);
  u.chain(lower + 1, 0, 0, 512);
  u.split(left + 1, 1, 300);
  t0 = std::chrono::steady_clock::now();
  double ratio = monotoneLogNormalizerRatio(u.tree, store, dir, left + 1, s);
  double ratioElapsed = seconds(t0);
  double a = 515.0, b = 515.0;
  checkNear(ratio, std::log(a * b / (a + b - 1.0)), 1e-9,
            "monotone scale: two-component ratio closed form");
  checkNear(ratio, directLogRatio(u.tree, store, dir, left + 1, direct), 1e-9,
            "monotone scale: two-component ratio vs direct counts");
  check(logChoose(1030, 515) > std::log(DBL_MAX),
        "monotone scale: the binomial is past doubles");
  printf(
      "ok: monotone scale (1031-leaf count %.3f s, 1030-leaf ratio %.3f s)\n",
      elapsed, ratioElapsed);
}

// A factor column that gains missing values gains a value on its axis, which
// can relate leaves the order did not: x1 constrained, cut once, each half cut
// on the factor f with the partition {a, b} | {c, d} labelled the other way on
// the other half. Without missing values the order is two pairs and these
// leaf values are feasible; once f has one, it goes left at both f rules and
// relates the two left leaves, whose values then fall along x1. An accepted
// predictor update, whole-matrix or row by row, reseeds the tree to all-zero.
static void testMonotoneMissingArrives() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20261001u);
  const size_t n = 400;
  std::vector<double> x(2 * n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 4);
    x[n + i] = runif01();
    y[i] = x[n + i];
  }
  const std::int8_t dir[2] = {0, 1};
  const ColumnKind types[2] = {ColumnKind::categorical, ColumnKind::numeric};
  SamplerOptions options;
  options.numTrees = 1;
  options.birthOrDeathProbability = 1.0;
  options.swapProbability = 0.0;
  options.changeProbability = 0.0;
  options.monotoneDirections = dir;
  options.predictors.columnTypes = types;

  std::vector<double> withMissing(x);
  withMissing[5] = std::numeric_limits<double>::quiet_NaN();
  for (int path = 0; path < 2; ++path) {
    Sampler<MonotoneConstantGaussianLeaf> sampler(
        x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
        1.0, 3.0, 0.37804942330213542, options, &rng);
    const ColumnStore& store(sampler.data());
    const std::vector<double>& cuts(store.cutPoints[1]);
    double cut = *std::lower_bound(cuts.begin(), cuts.end(), 0.5);
    std::vector<FlatNode> flat(7, FlatNode());
    const double values[4] = {0.05, -0.05, -0.04, 0.06};
    const std::uint64_t masks[2] = {0xcu, 0x3u};
    flat[0].variable = 1;
    flat[0].value = cut;
    setFlatKind(flat[0], FlatKind::ordinal);
    for (int half = 0; half < 2; ++half) {
      FlatNode& rule = flat[1 + 3 * half];
      rule.variable = 0;
      rule.mask = masks[half];
      setFlatKind(rule, FlatKind::categoricalInline);
      flat[2 + 3 * half].value = values[2 * half];
      flat[3 + 3 * half].value = values[2 * half + 1];
    }
    SamplerStateData state;
    sampler.getState(state);
    state.chains[0].forests[0].trees[0] = flat;
    check(sampler.setState(state, nullptr),
          "monotone missing arrives: the state installs");

    auto liveValues = [&](std::vector<double>& out) {
      std::vector<FlatNode> live;
      std::vector<std::uint32_t> counts;
      sampler.flattenTree(0, 0, live, counts);
      out.clear();
      for (const FlatNode& node : live)
        if (node.variable == invalidVariable) out.push_back(node.value);
    };
    std::vector<double> before, after;
    liveValues(before);
    check(before.size() == 4 && before[0] == values[0] &&
              before[2] == values[2],
          "monotone missing arrives: feasible without missing values");
    if (path == 0) {
      check(sampler.setPredictor(withMissing.data(), false, false) ==
                PredictorUpdateResult::accepted,
            "monotone missing arrives: the update is accepted");
    } else {
      std::unique_ptr<bool[]> installed(new bool[n]);
      check(sampler.updatePredictorPerObservation(withMissing.data(), 0,
                                                  installed.get()) &&
                installed[5],
            "monotone missing arrives: the row update installs");
    }
    check(sampler.data().hasMissing[0] != 0,
          "monotone missing arrives: the column has a missing value");
    liveValues(after);
    bool zero = after.size() == 4;
    for (double v : after) zero = zero && v == 0.0;
    check(zero, path == 0
                    ? "monotone missing arrives: setPredictor reseeds"
                    : "monotone missing arrives: a row update reseeds");
  }
  ext_rng_destroy(rng);
  printf("ok: monotone order gains the missing value an update brings\n");
}

void runMonotoneTests() {
  std::uint64_t saved = rngState;
  rngState = 7071u;
  testMonotoneCountHandBuilt();
  {
    ColumnStore store;
    makeStore(store, 3, 12, 300);
    testMonotoneCountRandom(
        store, {{1, 0, 0}, {-1, 1, 0}, {1, -1, 1}, {-1, -1, -1}, {0, -1, 1}},
        "");
  }
  {
    ColumnStore store;
    makeStore(store, 3, 12, 300);
    testMonotoneRatioRandom(
        store, {{1, 0, 0}, {-1, 0, 1}, {1, -1, 0}, {-1, 0, 0}, {0, 1, -1}},
        "");
  }
  testMonotoneRatioNamed();
  testMonotonePositionLaw();
  testMonotoneExtensionDraw();
  testMonotoneCountScale();
  testMonotoneGeometryPoints();
  testMonotoneMissingArrives();
  {  // factor splits and missing values, in free and constrained predictors
    ColumnStore store;
    makeMixedStore(store, 700);
    const Directions dirs = {
        {1, 0, 0, 0}, {-1, 1, 0, 0}, {1, -1, 0, 0}, {0, 1, 0, 0}};
    testMonotoneCountRandom(store, dirs, " with factors and missing values");
    testMonotoneRatioRandom(store, dirs, " with factors and missing values");
  }
  rngState = saved;
}
