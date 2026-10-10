#include "common.hpp"

#include <chrono>
#include <functional>
#include <future>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <thread>

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
  // variances too: an isolated leaf drawn at the constrained sd keeps its
  // mean at zero and moves only its spread. The z is against the normal-
  // theory sd of a sample variance, 2 v^2 / N per arm.
  double worstZ = 0.0, worstVarianceZ = 0.0;
  for (std::int32_t leaf : leaves) {
    double m1 = sumExact[leaf] / priorDraws, m2 = sumReject[leaf] / priorDraws;
    double v1 = sqExact[leaf] / priorDraws - m1 * m1,
           v2 = sqReject[leaf] / priorDraws - m2 * m2;
    worstZ = std::max(worstZ,
                      std::fabs(m1 - m2) / std::sqrt((v1 + v2) / priorDraws));
    worstVarianceZ =
        std::max(worstVarianceZ,
                 std::fabs(v1 - v2) /
                     std::sqrt(2.0 * (v1 * v1 + v2 * v2) / priorDraws));
  }
  check(worstZ < 4.5, "monotone prior draw: means match rejection");
  check(worstVarianceZ < 4.5,
        "monotone prior draw: variances match rejection, the isolated leaf's "
        "at the free sd");
  printf("ok: monotone exact prior draw (%zu leaves, worst mean |z| %.2f, "
         "variance |z| %.2f)\n",
         leaves.size(), worstZ, worstVarianceZ);
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

namespace {

using MonotoneFacade = SamplerFacade<MonotoneConstantGaussianLeaf>;

// The data of the missing-value tests is f, a 4-level factor with no missing
// value, x1 uniform and y rising in x1.
const size_t arrivalRows = 400;
const ColumnKind arrivalTypes[2] = {ColumnKind::categorical,
                                    ColumnKind::numeric};

// The hand tree: x1, constrained, cut once, and each half cut on f with the
// partition {a, b} | {c, d} labelled the other way on the other half. Leaves
// in pre-order: (low; a, b), (low; c, d), (high; c, d), (high; a, b). Under
// an increasing constraint and without missing values the order is two
// pairs, each (low; S) below (high; S). Once f has a missing value each f
// rule draws the side it goes to, the low one first. Sent left at both it
// relates the two left leaves, (low; a, b) below (high; c, d), and sent right
// at both the two right ones, (low; c, d) below (high; a, b); sent one way
// each it joins two leaves already related. So values can be out of order
// under one of the two like draws and under no other, and a chain that is to
// break is seeded so that its draw is a like one and given the values that
// draw breaks. A decreasing constraint mirrors the order, and the values.
const double breaksSentLeft[4] = {0.05, -0.05, -0.04, 0.06};
const double breaksSentRight[4] = {-0.05, 0.05, 0.06, -0.04};
const double holdsEitherWay[4] = {-0.05, -0.06, 0.04, 0.06};

// Which sampler: per chain whether its hand tree's values break once f has a
// missing value, how many trees a chain has and which of them is the hand
// tree, the others being single leaves, and the direction on x1.
struct ArrivalShape {
  std::vector<bool> breaks;
  size_t numTrees = 1, at = 0;
  std::int8_t direction = 1;
};

// A sampler of that shape over the data above, chain c seeded by c or, where
// it is to break, by the first seed after it whose two coins agree.
// scansFirst says chain 0's generator draws a scan order before the coins: a
// row-by-row update of this sampler, or a joint one it leads.
struct Arrival {
  ArrivalShape shape;
  std::int8_t directions[2];
  std::vector<ext_rng*> rngs;
  std::vector<const double*> values;
  std::unique_ptr<MonotoneFacade> facade;
  size_t numChains;
  double cut = 0.0;

  // the two coins chain c's generator gives next, as a first missing value
  // will draw them; the generator is left where it was
  void nextCoins(size_t c, bool scansFirst, bool coins[2]) {
    ext_rng* peek = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    std::vector<unsigned char> image(ext_rng_getSerializedStateLength(rngs[c]));
    ext_rng_writeSerializedState(rngs[c], image.data());
    ext_rng_readSerializedState(peek, image.data());
    if (c == 0 && scansFirst) {
      std::vector<size_t> order(arrivalRows);
      ext_rng_drawPermutation(peek, order.data(), arrivalRows);
    }
    for (int k = 0; k < 2; ++k)
      coins[k] = ext_rng_simulateBernoulli(peek, 0.5) == 1;
    ext_rng_destroy(peek);
  }

  Arrival(const std::vector<double>& x, const std::vector<double>& y,
          const ArrivalShape& s, bool scansFirst = false)
      : shape(s), directions{0, s.direction}, numChains(s.breaks.size()) {
    for (size_t c = 0; c < numChains; ++c) {
      rngs.push_back(
          ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr));
      ext_rng_setSeed(rngs.back(), 20261001u + static_cast<std::uint32_t>(c));
    }
    SamplerOptions options;
    options.numTrees = shape.numTrees;
    options.numChains = numChains;
    options.numThreads = 1;
    options.birthOrDeathProbability = 1.0;
    options.swapProbability = 0.0;
    options.changeProbability = 0.0;
    options.monotoneDirections = directions;
    options.predictors.columnTypes = arrivalTypes;
    facade = std::make_unique<MonotoneFacade>(
        x.data(), y.data(), arrivalRows, size_t(2), nullptr, nullptr,
        ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options,
        rngs.data());
    const std::vector<double>& cuts(sampler().data().cutPoints[1]);
    cut = *std::lower_bound(cuts.begin(), cuts.end(), 0.5);
    // nothing between here and the update draws from a chain's generator
    values.assign(numChains, holdsEitherWay);
    for (size_t c = 0; c < numChains; ++c) {
      if (!shape.breaks[c]) continue;
      bool coins[2];
      nextCoins(c, scansFirst, coins);
      for (std::uint32_t attempt = 1; coins[0] != coins[1]; ++attempt) {
        ext_rng_setSeed(rngs[c], 20261001u + static_cast<std::uint32_t>(c) +
                                     1000u * attempt);
        nextCoins(c, scansFirst, coins);
      }
      values[c] = coins[0] ? breaksSentRight : breaksSentLeft;
    }
    SamplerStateData state;
    sampler().getState(state);
    for (size_t c = 0; c < numChains; ++c) {
      std::vector<FlatNode> flat(7, FlatNode());
      const std::uint64_t masks[2] = {0xcu, 0x3u};
      flat[0].variable = 1;
      flat[0].value = cut;
      setFlatKind(flat[0], FlatKind::ordinal);
      for (int half = 0; half < 2; ++half) {
        FlatNode& rule = flat[1 + 3 * half];
        rule.variable = 0;
        rule.mask = masks[half];
        setFlatKind(rule, FlatKind::categoricalInline);
        flat[2 + 3 * half].value = installed(c, 2 * half);
        flat[3 + 3 * half].value = installed(c, 2 * half + 1);
      }
      state.chains[c].forests[0].trees[shape.at] = flat;
    }
    check(sampler().setState(state, nullptr),
          "monotone missing arrives: the state installs");
  }
  ~Arrival() {
    for (ext_rng* rng : rngs) ext_rng_destroy(rng);
  }

  Sampler<MonotoneConstantGaussianLeaf>& sampler() { return facade->impl(); }
  bool factorHasMissing() { return sampler().data().hasMissing[0] != 0; }
  // the value installed at leaf k of chain c's hand tree
  double installed(size_t c, size_t k) const {
    return shape.direction * values[c][k];
  }
  // whether the leaves of chain c's hand tree are the installed values, or
  // all zero
  bool leavesAre(size_t c, bool zero) {
    std::vector<FlatNode> live;
    std::vector<std::uint32_t> counts;
    sampler().flattenTree(c, shape.at, live, counts);
    size_t numLeaves = 0;
    bool same = true;
    for (const FlatNode& node : live)
      if (node.variable == invalidVariable) {
        same = same && numLeaves < 4 &&
               node.value == (zero ? 0.0 : installed(c, numLeaves));
        ++numLeaves;
      }
    return same && numLeaves == 4;
  }
  bool leavesKept(size_t c) { return leavesAre(c, false); }
  bool leavesZero(size_t c) { return leavesAre(c, true); }
  void sweep(size_t numSweeps) {
    Results none;
    sampler().run(numSweeps, 0, none);
  }
};

// What of b differs from a among the things a refused update leaves alone,
// or null: the codes, has-missing flags and cut grids, and per chain the
// generator, every tree with its leaf values and row-to-leaf map, the cached
// fits and sigma. At tolerance 0 all of it bit for bit; a positive tolerance
// is for the leaf values, fits and sigma alone. The order a leaf's rows are
// held in is left out: a refusal re-routes twice and may leave another. Each
// partition is held to its own chain's map instead.
const char* arrivalDifference(Arrival& a, Arrival& b, double tolerance) {
  auto same = [tolerance](double u, double v) {
    return tolerance > 0.0 ? std::fabs(u - v) <= tolerance
                           : std::memcmp(&u, &v, sizeof(double)) == 0;
  };
  auto& sa = a.sampler();
  auto& sb = b.sampler();
  const ColumnStore& da(sa.data());
  const ColumnStore& db(sb.data());
  size_t n = da.numObservations;
  if (storageDigest(da) != storageDigest(db)) return "codes";
  if (da.hasMissing != db.hasMissing) return "has-missing flags";
  if (da.cutPoints != db.cutPoints) return "cut points";
  for (size_t c = 0; c < a.numChains; ++c) {
    std::vector<unsigned char> ga(ext_rng_getSerializedStateLength(a.rngs[c]));
    std::vector<unsigned char> gb(ga.size());
    ext_rng_writeSerializedState(a.rngs[c], ga.data());
    ext_rng_writeSerializedState(b.rngs[c], gb.data());
    if (ga != gb) return "generator";
    const std::vector<double>& ta(TestPeer::totalFitsInForest(sa.chain(c), 0));
    const std::vector<double>& tb(TestPeer::totalFitsInForest(sb.chain(c), 0));
    for (size_t i = 0; i < n; ++i)
      if (!same(ta[i], tb[i])) return "fits";
    if (!same(sa.sigma(c), sb.sigma(c))) return "sigma";
    for (size_t t = 0; t < a.shape.numTrees; ++t) {
      std::vector<FlatNode> fa, fb;
      std::vector<std::uint32_t> counts;
      sa.flattenTree(c, t, fa, counts);
      sb.flattenTree(c, t, fb, counts);
      if (fa.size() != fb.size()) return "tree";
      for (size_t i = 0; i < fa.size(); ++i) {
        bool leaf = fa[i].variable == invalidVariable;
        if (fa[i].variable != fb[i].variable || fa[i].flags != fb[i].flags ||
            (!leaf && fa[i].mask != fb[i].mask))
          return "tree";
        if (leaf && !same(fa[i].value, fb[i].value)) return "leaf values";
      }
      const std::uint32_t* la = TestPeer::leafOf(sa.chain(c), t);
      const std::uint32_t* lb = TestPeer::leafOf(sb.chain(c), t);
      if (std::memcmp(la, lb, n * sizeof(std::uint32_t)) != 0)
        return "row-to-leaf map";
      for (int which = 0; which < 2; ++which) {
        const Tree& tree = (which == 0 ? sa : sb).chain(c).tree(t);
        const std::uint32_t* leaf = which == 0 ? la : lb;
        std::vector<std::int32_t> bottoms;
        tree.fillBottom(0, bottoms);
        size_t numRows = 0;
        for (std::int32_t bottom : bottoms) {
          const Node& node(tree.at(bottom));
          numRows += node.numObservations();
          for (size_t m = node.begin; m < node.end; ++m)
            if (leaf[tree.indices[m]] != static_cast<std::uint32_t>(bottom))
              return "partition against its map";
        }
        if (numRows != n) return "partition against its map";
      }
    }
  }
  return nullptr;
}

void checkAt(bool condition, const std::string& label, const char* what) {
  std::string message = "monotone missing arrives: " + label + ": " + what;
  check(condition, message.c_str());
}

// a and b are twins bit for bit, and three sweeps of each later draw the same
// trees with values, fits and sigma equal to rounding. Not bit for bit then:
// a sweep sums over a leaf's rows in the order the leaf holds them, which any
// re-route can change, a refused update's included.
void checkTwins(Arrival& a, Arrival& b, const std::string& label) {
  for (int pass = 0; pass < 2; ++pass) {
    if (pass == 1) {
      a.sweep(3);
      b.sweep(3);
    }
    const char* difference = arrivalDifference(a, b, pass == 0 ? 0.0 : 1e-12);
    std::string what = pass == 0 ? "the sampler is its twin"
                                 : "three sweeps later it draws as its twin";
    if (difference != nullptr)
      what += std::string(", differs in ") + difference;
    checkAt(difference == nullptr, label, what.c_str());
  }
}

size_t countTrue(const bool* flags, size_t n) {
  size_t count = 0;
  for (size_t i = 0; i < n; ++i) count += flags[i] ? 1 : 0;
  return count;
}

}  // namespace

// A factor column's first missing value on a sampler with a monotone
// constraint. Unforced, an update that brings it while a tree's leaf values
// would then be out of order, under the directions the value's arrival
// draws, is refused and the sampler is an untouched twin's, its generators
// and directions included: whole matrix and by column, with and without a
// cut refresh, and row by row, where each row bringing the value is refused
// and the others install. Where the values stay in order it is accepted and
// they are kept.
// Forced, and through setData, it is taken and such a tree is set to zero.
// The tree is a chain's only one, or one of three and not the first, in every
// chain or the second alone, under either direction.
static void testMonotoneMissingArrives() {
  const double na = std::numeric_limits<double>::quiet_NaN();
  const size_t n = arrivalRows;
  std::vector<double> x(2 * n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 4);
    x[n + i] = runif01();
    y[i] = x[n + i];
  }
  const size_t factorColumn = 0;
  const size_t bothColumns[2] = {0, 1};

  // the row given the missing value is one it moves to another leaf, x1 above
  // the cut and level a or b, so a fit rebuilt from the new partition shows;
  // `second` is another such row and `moved` a row whose level changes
  size_t row = n, second = n, moved = n;
  {
    Arrival probe(x, y, {{true}});
    check(probe.leavesKept(0) && !probe.factorHasMissing(),
          "monotone missing arrives: in order without missing values");
    for (size_t i = 0; i < n; ++i) {
      bool moves = x[n + i] > probe.cut && x[i] < 2.0;
      if (moves && row == n) row = i;
      else if (moves && second == n) second = i;
      else if (!moves && moved == n) moved = i;
    }
  }
  check(row < n && second < n && moved < n,
        "monotone missing arrives: the fixture has its three rows");
  if (row == n || second == n || moved == n) return;

  // the matrix with the missing value; the same with x1 rescaled, which keeps
  // every code under a cut refresh and moves every cut point
  std::vector<double> withMissing(x), rescaled(x);
  withMissing[row] = na;
  rescaled[row] = na;
  for (size_t i = 0; i < n; ++i) rescaled[n + i] = 0.25 + 0.5 * x[n + i];
  // the factor column with two missing values and one row moved to another
  // level, and the column a call that refuses the two leaves behind
  std::vector<double> twoMissing(x.begin(), x.begin() + n);
  twoMissing[moved] = std::fmod(twoMissing[moved] + 1.0, 4.0);
  std::vector<double> twoRefused(twoMissing);
  twoMissing[row] = na;
  twoMissing[second] = na;
  // a matrix whose factor holds a missing value from the start, elsewhere
  std::vector<double> alreadyMissing(x);
  alreadyMissing[moved] = na;
  std::vector<double> anotherMissing(alreadyMissing);
  anotherMissing[row] = na;

  // one unforced update of `whole` or, row by row, of the factor `column`:
  // whole matrix and by column, each without and with a cut refresh over the
  // rescaled matrix, then row by row
  std::unique_ptr<bool[]> installed(new bool[n]), twinInstalled(new bool[n]);
  static const char* const unforcedForms[5] = {
      ", whole matrix", ", by column", ", whole matrix with a cut refresh",
      ", by column with a cut refresh", ", row by row"};
  auto update = [&](Arrival& s, int form, bool* valid, const double* whole,
                    const double* column) {
    auto& sampler = s.sampler();
    PredictorUpdateResult result = PredictorUpdateResult::accepted;
    *valid = true;
    if (form == 0)
      result = sampler.setPredictor(whole, false, false);
    else if (form == 1)
      result = sampler.updatePredictor(whole, &factorColumn, 1, false, false);
    else if (form == 2)
      result = sampler.setPredictor(rescaled.data(), false, true);
    else if (form == 3)
      result =
          sampler.updatePredictor(rescaled.data(), bothColumns, 2, false, true);
    else
      *valid =
          sampler.updatePredictorPerObservation(column, 0, installed.get());
    return result;
  };

  // out of order once f has a missing value in every chain or the second
  // alone, and in order either way; three trees with the hand tree second or
  // last where a check of a chain's first tree alone would pass
  struct Layout {
    const char* name;
    ArrivalShape shape;
  };
  const Layout layouts[5] = {
      {"one chain", {{true}}},
      {"two chains of three trees", {{true, true}, 3, 1}},
      {"the second chain alone, decreasing", {{false, true}, 3, 2, -1}},
      {"in order, one chain", {{false}}},
      {"in order, two chains, decreasing", {{false, false}, 3, 1, -1}}};
  int numRefused = 0, numAccepted = 0, numForced = 0;
  for (const Layout& layout : layouts) {
    const ArrivalShape& shape = layout.shape;
    size_t numChains = shape.breaks.size();
    bool breaks = std::count(shape.breaks.begin(), shape.breaks.end(), true);
    for (int form = 0; form < 5; ++form) {
      std::string label = std::string(layout.name) + unforcedForms[form];
      Arrival s(x, y, shape, form == 4), twin(x, y, shape, form == 4);
      std::vector<double> cutsBefore(s.sampler().data().cutPoints[1]);
      bool valid;
      PredictorUpdateResult result =
          update(s, form, &valid, withMissing.data(), twoMissing.data());
      checkAt(valid, label, "the session ends valid");
      if (!breaks) {
        checkAt(result == PredictorUpdateResult::accepted &&
                    (form < 4 || countTrue(installed.get(), n) == n),
                label, "accepted");
        bool kept = s.factorHasMissing();
        for (size_t c = 0; c < numChains; ++c) kept = kept && s.leavesKept(c);
        checkAt(kept, label, "the values are kept and the flag is set");
        // the refresh the refused forms roll back is one that moves the grid
        if (form == 2 || form == 3)
          checkAt(s.sampler().data().cutPoints[1] != cutsBefore, label,
                  "the cut points moved");
        ++numAccepted;
        continue;
      }
      if (form < 4) {
        checkAt(result == PredictorUpdateResult::rolledBack, label, "refused");
      } else {
        // the twin is given the column the call left
        checkAt(!installed[row] && !installed[second], label,
                "each row bringing the missing value is refused");
        checkAt(countTrue(installed.get(), n) == n - 2, label,
                "every other row installs");
        checkAt(twin.sampler().updatePredictorPerObservation(
                    twoRefused.data(), 0, twinInstalled.get()) &&
                    countTrue(twinInstalled.get(), n) == n,
                label, "the twin takes the column without the two");
      }
      checkAt(!s.factorHasMissing(), label, "the factor's flag is clear");
      checkTwins(s, twin, label);
      ++numRefused;
    }

    // forced whole and by column, and setData: taken, a tree left out of
    // order set to zero and a tree in order kept
    for (int form = 0; form < 3; ++form) {
      static const char* const forms[3] = {", forced whole matrix",
                                           ", forced by column", ", setData"};
      std::string label = std::string(layout.name) + forms[form];
      Arrival s(x, y, shape);
      auto& sampler = s.sampler();
      bool taken;
      if (form == 0)
        taken = sampler.setPredictor(withMissing.data(), true, false) ==
                PredictorUpdateResult::accepted;
      else if (form == 1)
        taken = sampler.updatePredictor(withMissing.data(), &factorColumn, 1,
                                        true, false) ==
                PredictorUpdateResult::accepted;
      else
        taken = sampler.setData(withMissing.data(), y.data(), n, nullptr,
                                nullptr, nullptr, 0);
      checkAt(taken && s.factorHasMissing(), label, "taken");
      bool reset = true;
      for (size_t c = 0; c < numChains; ++c)
        reset = reset && (shape.breaks[c] ? s.leavesZero(c) : s.leavesKept(c));
      checkAt(reset, label, "a tree left out of order is zero, another kept");
      ++numForced;
    }
  }

  // two samplers swept jointly, in either order: a row one of them refuses is
  // declined in both, and each ends a twin swept with the column without it
  for (int first = 0; first < 2; ++first) {
    std::string label = first == 0 ? "jointly, the sampler in order first"
                                   : "jointly, the sampler in order second";
    const size_t columns[2] = {0, 0};
    // the sampler that leads the sweep draws its scan order
    Arrival held(x, y, {{false}}), broken(x, y, {{true}}, first == 1);
    Arrival heldTwin(x, y, {{false}}), brokenTwin(x, y, {{true}}, first == 1);
    SamplerBase* samplers[2] = {held.facade.get(), broken.facade.get()};
    SamplerBase* twins[2] = {heldTwin.facade.get(), brokenTwin.facade.get()};
    if (first == 1) {
      std::swap(samplers[0], samplers[1]);
      std::swap(twins[0], twins[1]);
    }
    checkAt(updatePredictorPerObservationJointly(
                samplers, 2, twoMissing.data(), columns, installed.get()),
            label, "both sessions end valid");
    checkAt(!installed[row] && !installed[second] &&
                countTrue(installed.get(), n) == n - 2,
            label, "the two rows are declined, the others install");
    checkAt(!held.factorHasMissing() && !broken.factorHasMissing(), label,
            "neither flag is set");
    checkAt(updatePredictorPerObservationJointly(
                twins, 2, twoRefused.data(), columns, twinInstalled.get()),
            label, "the twins take the column without the two");
    checkTwins(held, heldTwin, label + ", the one in order");
    checkTwins(broken, brokenTwin, label + ", the other");
    ++numRefused;
  }

  // a sampler with kept draws, swept jointly behind and ahead of one that
  // refuses the missing value. Leading, its session finds the row valid and
  // draws for its live and kept rules before the other declines, and takes
  // the draw back: it ends as a twin swept with the column without the two
  // rows, its kept draws and generators included
  {
    using PlainFacade = SamplerFacade<ConstantGaussianLeaf>;
    std::vector<double> yPlain(n);
    for (size_t i = 0; i < n; ++i)
      yPlain[i] = (x[i] < 2.0 ? -1.0 : 1.0) + x[n + i];
    struct Plain {
      std::vector<ext_rng*> rngs;
      std::unique_ptr<PlainFacade> facade;
      Plain(const std::vector<double>& x, const std::vector<double>& y) {
        for (std::uint32_t c = 0; c < 2; ++c) {
          rngs.push_back(
              ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr));
          ext_rng_setSeed(rngs.back(), 20261010u + c);
        }
        SamplerOptions options;
        options.numTrees = 10;
        options.numChains = 2;
        options.numThreads = 1;
        options.keepTrees = true;
        options.numSamplesToStore = 4;
        options.predictors.columnTypes = arrivalTypes;
        facade = std::make_unique<PlainFacade>(
            x.data(), y.data(), arrivalRows, size_t(2), nullptr, nullptr,
            ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options,
            rngs.data());
        Results none;
        facade->impl().run(60, 4, none);
      }
      ~Plain() {
        facade.reset();
        for (ext_rng* rng : rngs) ext_rng_destroy(rng);
      }
      SamplerStateData state() {
        SamplerStateData result;
        facade->impl().getState(result);
        return result;
      }
      bool flagged() { return facade->impl().data().hasMissing[0] != 0; }
    };
    auto numRight = [](const SamplerStateData& state) {
      size_t count = 0;
      for (const ChainStateData& chain : state.chains)
        for (const std::vector<FlatNode>& tree : chain.forests[0].savedTrees)
          for (const FlatNode& node : tree)
            count += (node.flags & flatMissingGoesRight) != 0 ? 1u : 0u;
      return count;
    };
    // alone it takes the two rows, and its kept rules draw
    {
      Plain alone(x, yPlain);
      bool valid = alone.facade->impl().updatePredictorPerObservation(
          twoMissing.data(), 0, installed.get());
      checkAt(valid && installed[row] && installed[second] && alone.flagged() &&
                  numRight(alone.state()) > 0,
              "a sampler with kept draws",
              "alone it takes the rows and its kept rules draw");
    }
    for (int first = 0; first < 2; ++first) {
      std::string label = first == 0
                              ? "a sampler with kept draws leading the sweep"
                              : "a sampler with kept draws led in the sweep";
      const size_t columns[2] = {0, 0};
      Plain plain(x, yPlain), plainTwin(x, yPlain);
      Arrival broken(x, y, {{true}}, first == 1);
      Arrival brokenTwin(x, y, {{true}}, first == 1);
      SamplerBase* samplers[2] = {plain.facade.get(), broken.facade.get()};
      SamplerBase* twins[2] = {plainTwin.facade.get(),
                               brokenTwin.facade.get()};
      if (first == 1) {
        std::swap(samplers[0], samplers[1]);
        std::swap(twins[0], twins[1]);
      }
      bool valid = updatePredictorPerObservationJointly(
          samplers, 2, twoMissing.data(), columns, installed.get());
      checkAt(valid && !installed[row] && !installed[second] &&
                  countTrue(installed.get(), n) == n - 2,
              label, "the two rows are declined, the others install");
      checkAt(updatePredictorPerObservationJointly(
                  twins, 2, twoRefused.data(), columns, twinInstalled.get()) &&
                  countTrue(twinInstalled.get(), n) == n,
              label, "the twins take the column without the two");
      SamplerStateData held = plain.state(), heldTwin = plainTwin.state();
      bool generators = true;
      for (size_t c = 0; c < held.chains.size(); ++c)
        generators = generators &&
                     held.chains[c].rngState == heldTwin.chains[c].rngState;
      checkAt(!plain.flagged() && !broken.factorHasMissing() && generators &&
                  numRight(held) == 0 && statesAgree(held, heldTwin),
              label,
              "its flag is clear and its trees, kept draws and generators are "
              "its twin's");
      checkTwins(broken, brokenTwin, label + ", the other");
      ++numRefused;
    }
  }

  // a column that already holds a missing value takes another
  for (int form : {0, 1, 4}) {
    std::string label =
        std::string("a second missing value") + unforcedForms[form];
    Arrival s(alreadyMissing, y, {{false, false}});
    checkAt(s.factorHasMissing(), label, "the column holds one from the start");
    bool valid;
    bool accepted = update(s, form, &valid, anotherMissing.data(),
                           anotherMissing.data()) ==
                        PredictorUpdateResult::accepted &&
                    valid && (form < 4 || countTrue(installed.get(), n) == n);
    bool kept = s.sampler().data().codeAt(0, row) == missingCategoryCode(4) &&
                s.leavesKept(0) && s.leavesKept(1);
    checkAt(accepted && kept, label, "accepted, the values kept");
    ++numAccepted;
  }
  printf("ok: monotone order and the missing value an update brings (%d "
         "refused against a twin, %d accepted, %d forced)\n",
         numRefused, numAccepted, numForced);
}

// Which has-missing flag can relate leaves the order did not, on random trees
// through the engine's geometry. Over a store with no missing value no rule
// carries a direction, so a first missing value goes left at every rule on
// its column: it changes no relation on a numeric or an ordered column, and
// on an unordered factor's adds some and removes none. Read with the flag
// down, a store whose rules hold directions relates fewer leaves and no
// others.
static void testMonotoneMissingRelates() {
  // x0, x1 constrained; x2 numeric, x3 ordered, x4 a 4-level factor and x5 a
  // 70-level one (pooled), all free
  const size_t n = 700, p = 6;
  const std::int8_t dir[p] = {1, -1, 0, 0, 0, 0};
  const ColumnKind types[p] = {
      ColumnKind::numeric,       ColumnKind::numeric,
      ColumnKind::numeric,       ColumnKind::orderedFactor,
      ColumnKind::categorical,   ColumnKind::categorical};
  const double levels[p] = {0.0, 0.0, 0.0, 5.0, 4.0, 70.0};
  const int numTrees = 4000;
  // per column, the related pairs the flag change adds, removes and reverses
  struct Tally {
    long added = 0, removed = 0, reversed = 0;
  };
  Tally gained[p], lost[p];
  long numRelated = 0;
  bool flagsAsBuilt = true;
  std::vector<int> before, after;
  auto relations = [&](const Tree& tree, const ColumnStore& store,
                       std::vector<int>& out) {
    std::vector<std::int32_t> leaves;
    tree.fillBottom(0, leaves);
    MonotoneLeafGeometry geometry;
    geometry.build(tree, store, leaves);
    out.assign(leaves.size() * leaves.size(), 0);
    for (size_t l = 0; l < leaves.size(); ++l)
      for (size_t r = l + 1; r < leaves.size(); ++r)
        out[l * leaves.size() + r] = geometry.relation(l, r, dir);
  };
  auto compare = [&](Tally& tally) {
    for (size_t k = 0; k < before.size(); ++k) {
      if (before[k] == after[k]) continue;
      if (before[k] == 0) ++tally.added;
      else if (after[k] == 0) ++tally.removed;
      else ++tally.reversed;
    }
  };
  for (int withMissing = 0; withMissing < 2; ++withMissing) {
    const double na = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> x(n * p);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < p; ++j) {
        double u = runif01();
        x[j * n + i] = j == 5 ? static_cast<double>(i % 70)
                       : levels[j] > 0.0 ? std::floor(levels[j] * u) : u;
        if (withMissing && runif01() < 0.15) x[j * n + i] = na;
      }
    ColumnStore store;
    built(store.build(x.data(), n, p, 6, false, types));
    for (size_t j = 0; j < p; ++j)
      flagsAsBuilt = flagsAsBuilt && (store.hasMissing[j] != 0) == withMissing;
    for (int trial = 0; trial < numTrees; ++trial) {
      TestTree t(store);
      t.growRandom(2 + trial % 9);
      relations(t.tree, store, before);
      for (int relation : before) numRelated += relation != 0;
      for (size_t j = 0; j < p; ++j) {
        Tree changed = t.tree;
        store.hasMissing[j] = withMissing ? 0 : 1;
        relations(changed, store, after);
        compare(withMissing ? lost[j] : gained[j]);
        store.hasMissing[j] = withMissing ? 1 : 0;
      }
    }
  }
  bool thresholdsUnchanged = true, noneLost = true, noneGained = true;
  long numRemoved = 0;
  for (size_t j = 0; j < p; ++j) {
    thresholdsUnchanged =
        thresholdsUnchanged && (j > 3 || gained[j].added == 0);
    noneLost = noneLost && gained[j].removed == 0 && gained[j].reversed == 0;
    noneGained = noneGained && lost[j].added == 0 && lost[j].reversed == 0;
    numRemoved += lost[j].removed;
  }
  check(flagsAsBuilt && numRelated > 20000,
        "monotone missing relates: the stores and trees are as described");
  check(thresholdsUnchanged && noneLost && gained[4].added > 0,
        "monotone missing relates: a first missing value adds relations "
        "through an unordered factor alone, and removes and reverses none");
  check(noneGained && numRemoved > 0,
        "monotone missing relates: read with the flag down, relations are "
        "removed and none added or reversed");
  printf("ok: monotone order and a has-missing flag (%d trees each way; a "
         "first missing value adds %ld pairs through the 4-level factor, "
         "%ld through the pooled one, none through the others; the flag "
         "down removes %ld)\n",
         numTrees, gained[4].added, gained[5].added, numRemoved);
}

namespace {

// The constant leaf's integrated marginal and posterior at prior sd tau, as
// the engine drops the response's raw sum of squares (it cancels in a move).
double refBase(double sw, double swy, double sig2, double tau) {
  double pp = 1.0 / (tau * tau), prec = pp + sw / sig2;
  return 0.5 * std::log(pp / prec) + 0.5 * (swy / sig2) * (swy / sig2) / prec;
}
void refPost(double sw, double swy, double sig2, double tau, double& m,
             double& sd) {
  double prec = 1.0 / (tau * tau) + sw / sig2;
  m = (swy / sig2) / prec;
  sd = std::sqrt(1.0 / prec);
}
double upperMass(double a, double m, double sd) {
  return 0.5 * std::erfc((a - m) / (sd * std::sqrt(2.0)));
}
// P(a <= X1 <= X2) for independent normals, composite Simpson over X2
double orderedPairMass(double a, double m1, double s1, double m2, double s2) {
  double lo = a, hi = std::max(a, m2 + 14.0 * s2);
  lo = std::max(lo, m2 - 14.0 * s2);
  if (!(hi > lo)) return 0.0;
  const int n = 40000;
  double h = (hi - lo) / n, sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    double x = lo + i * h;
    double f = std::exp(-0.5 * ((x - m2) / s2) * ((x - m2) / s2)) /
               (s2 * std::sqrt(2.0 * std::numbers::pi)) *
               (upperMass(a, m1, s1) - upperMass(x, m1, s1));
    sum += f * (i == 0 || i == n ? 1.0 : (i % 2 ? 4.0 : 2.0));
  }
  return sum * h / 3.0;
}

}  // namespace

// The move's acceptance, RNG-free in value: a death whose split the data
// strongly favours, so its ratio is far below 1 and the free bound cannot
// decide it, against the closed form of the corrected statement. T* is
// A < B1 < B2 (x1 cut twice; constrained) or A < {B1, B2} (B cut on the free
// x2), mu_A frozen, the touched leaves integrated. The "leaf" prior adds
// log(Z_T* / Z_T0): -log 3 and -log 3/2; "joint" nothing, and never counts.
// The old engine's value, the touched marginals divided by their prior cone
// mass d given mu_A, differs from both.
static void testMonotoneMoveClosedForm() {
  const size_t n = 400;
  std::vector<double> x(2 * n), y(n), weights(n, 1.0);
  for (size_t i = 0; i < n; ++i) {
    x[i] = (static_cast<double>(i) + 0.5) / n;
    x[n + i] = runif01();
  }
  ColumnStore store;
  built(store.build(x.data(), n, 2, 40));
  const double sig2 = 0.01, scale = 0.5, k = 2.0, muA = 0.3;
  const double c = std::sqrt(std::numbers::pi / (std::numbers::pi - 1.0));
  const double tauC = c * scale / k;
  const double base = 0.95, power = 2.0;
  CGMTreePrior prior;
  prior.base = base;
  prior.power = power;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  int checked = 0;
  for (int design = 0; design < 2; ++design) {
    bool freeSplit = design == 1;
    for (size_t i = 0; i < n; ++i) {
      double x1 = x[i], x2 = x[n + i];
      double level = x1 <= 0.5 ? 0.3 : (freeSplit ? (x2 > 0.5 ? 0.6 : 0.4)
                                                  : (x1 <= 0.75 ? 0.4 : 0.6));
      y[i] = level + 0.1 * (runif01() - 0.5);
    }
    for (MonotonePrior which : {MonotonePrior::leaf, MonotonePrior::joint}) {
      MonotoneConstantGaussianLeaf leaf;
      leaf.scale = scale;
      leaf.data = &store;
      leaf.directions = {1, 0};
      leaf.cInflation = c;
      leaf.prior = which;
      std::vector<index_t> idx(n);
      Tree tree;
      std::int32_t a = 0, b = 0, b1 = 0, b2 = 0;
      auto build = [&]() {
        tree.initialize(idx.data(), n);
        tree.computeLeafStats(0, y.data(), weights.data());
        Rule rule;
        rule.variableIndex = 0;
        rule.setSplitIndex(19);
        tree.birth(store, 0, rule, y.data(), weights.data());
        a = tree.at(0).leftChild;
        b = a + 1;
        rule.variableIndex = freeSplit ? 1 : 0;
        rule.setSplitIndex(freeSplit ? 19 : 29);
        tree.birth(store, b, rule, y.data(), weights.data());
        b1 = tree.at(b).leftChild;
        b2 = b1 + 1;
      };
      MoveScratch scratch;
      std::vector<double> mu;
      double alpha = -1.0;
      bool counted = false, died = false;
      for (std::uint32_t seed = 1; seed < 200 && !died; ++seed) {
        build();
        mu.assign(tree.nodes.size(), 0.0);
        mu[a] = muA;
        MoveContext ctx{store, prior, 1.0, 0.0, 0.0, 0.0, 0.5,
                        weights.data(), k, scratch};
        ctx.leafParams = mu.data();
        ext_rng_setSeed(rng, seed);
        leaf.count.downSets = 777777;
        bool stepTaken = false, wasBirth = true;
        alpha = birthOrDeathMove(ctx, leaf, rng, tree, y.data(),
                                 std::sqrt(sig2), &stepTaken, &wasBirth);
        died = !wasBirth;
        counted = leaf.count.downSets != 777777;
      }
      check(died, "monotone closed form: a death step");
      // T* again, for its statistics
      build();
      const Node& nB1 = tree.at(b1);
      const Node& nB2 = tree.at(b2);
      double swB = nB1.sumWeights + nB2.sumWeights;
      double swyB = nB1.sumWeightedResponse + nB2.sumWeightedResponse;
      double m, sd, m1, s1, m2, s2;
      refPost(swB, swyB, sig2, tauC, m, sd);
      refPost(nB1.sumWeights, nB1.sumWeightedResponse, sig2, tauC, m1, s1);
      refPost(nB2.sumWeights, nB2.sumWeightedResponse, sig2, tauC, m2, s2);
      double logI0 = refBase(swB, swyB, sig2, tauC) + std::log(upperMass(muA, m, sd));
      double logIStar =
        refBase(nB1.sumWeights, nB1.sumWeightedResponse, sig2, tauC) +
        refBase(nB2.sumWeights, nB2.sumWeightedResponse, sig2, tauC) +
        std::log(freeSplit ? upperMass(muA, m1, s1) * upperMass(muA, m2, s2)
                           : orderedPairMass(muA, m1, s1, m2, s2));
      double gB = base / std::pow(2.0, power), gC = base / std::pow(3.0, power);
      double logPrior =
        std::log(1.0 - gB) - std::log(gB) - 2.0 * std::log(1.0 - gC);
      // reverse birth: probability 1/2, one of two birthable leaves; the
      // death: probability 1/2, the one nog node
      double logTransition = std::log(0.5 / 2.0) - std::log(0.5 / 1.0);
      double joint = logPrior + logTransition + logI0 - logIStar;
      double logZStarOverZ0 = freeSplit ? std::log(2.0 / 3.0) : -std::log(3.0);
      double closed =
        joint + (which == MonotonePrior::leaf ? logZStarOverZ0 : 0.0);
      double dStar = freeSplit ? upperMass(muA, 0.0, tauC) * upperMass(muA, 0.0, tauC)
                               : orderedPairMass(muA, 0.0, tauC, 0.0, tauC);
      double old = joint + std::log(dStar) - std::log(upperMass(muA, 0.0, tauC));
      check(closed < -20.0, "monotone closed form: the death is unfavoured");
      checkNear(std::log(alpha), closed, 1e-6,
                "monotone closed form: the move's log alpha");
      check(std::fabs(closed - old) > 0.25,
            "monotone closed form: differs from the d-divided value");
      check(counted == (which == MonotonePrior::leaf),
            "monotone closed form: only the leaf prior counts");
      ++checked;
    }
  }
  ext_rng_destroy(rng);
  printf("ok: monotone move log alpha against the closed form (%d moves)\n",
         checked);
}

namespace {

// log P(a <= Z <= b) for a standard normal Z, from tails in logs; the cone
// references' own arithmetic, -Inf where rounding leaves no mass
double refLogMass(double a, double b) {
  if (!(b > a)) return -HUGE_VAL;
  double near, far;
  if (a > 0.0) {
    near = Rf_pnorm5(a, 0.0, 1.0, 0, 1);
    far = std::isinf(b) ? -HUGE_VAL : Rf_pnorm5(b, 0.0, 1.0, 0, 1);
  } else {
    near = std::isinf(b) ? 0.0 : Rf_pnorm5(b, 0.0, 1.0, 1, 1);
    far = std::isinf(a) ? -HUGE_VAL : Rf_pnorm5(a, 0.0, 1.0, 1, 1);
  }
  return near + std::log1p(-std::exp(std::min(far - near, 0.0)));
}

// log P(aL <= L <= min(bL, U), lo <= U <= hi) for independent
// L ~ N(mL, sL^2), U ~ N(mR, sR^2), by brute force: the log integrand on a
// 1e5-point grid over a range wide enough to hold it, the region within 60
// nats of that grid's peak, then composite Simpson on 2e5 panels of
// exp(log integrand - peak) over that region, split at bL so the kink sits on
// a node
double refLogCone(double lo, double hi, double aL, double bL, double mL,
                  double sL, double mR, double sR) {
  auto logF = [&](double x) {
    double z = (x - mR) / sR;
    return -0.5 * z * z - std::log(sR) -
           0.5 * std::log(2.0 * std::numbers::pi) +
           refLogMass((aL - mL) / sL, (std::min(bL, x) - mL) / sL);
  };
  double span = 100.0 * std::max(sL, sR);
  double from = std::isfinite(lo) ? lo : std::min(mL, mR) - span;
  double to = std::isfinite(hi) ? hi : std::max(std::max(mL, mR), from) + span;
  const int coarse = 100000;
  double h = (to - from) / coarse, peak = -HUGE_VAL;
  std::vector<double> f(coarse + 1);
  for (int i = 0; i <= coarse; ++i) {
    f[static_cast<size_t>(i)] = logF(from + i * h);
    peak = std::max(peak, f[static_cast<size_t>(i)]);
  }
  int first = coarse, last = 0;
  for (int i = 0; i <= coarse; ++i)
    if (f[static_cast<size_t>(i)] > peak - 60.0) {
      first = std::min(first, i);
      last = std::max(last, i);
    }
  double l = from + std::max(first - 2, 0) * h;
  double r = from + std::min(last + 2, coarse) * h;
  auto simpson = [&](double a, double b, int n) {
    double w = (b - a) / n, sum = 0.0;
    for (int i = 0; i <= n; ++i) {
      double v = std::exp(logF(a + i * w) - peak);
      sum += v * (i == 0 || i == n ? 1.0 : (i % 2 ? 4.0 : 2.0));
    }
    return sum * w / 3.0;
  };
  double total = bL > l && bL < r
    ? simpson(l, bL, 100000) + simpson(bL, r, 100000)
    : simpson(l, r, 200000);
  return peak + std::log(total);
}

}  // namespace

// The constrained pair's log cone probability against closed forms and a
// brute-force log-space reference, deep in the tail: data that run against
// the constraint by up to 200 sd, children whose posterior sds differ by up
// to 50 times, and a frozen lower bound up to 20 sd above the lower leaf's
// mean. A linear-space quadrature with an absolute tolerance is off by tens
// of nats here and returns zero, which the move takes for an infeasible cone.
static void testMonotoneConeDeepTail() {
  using Leaf = MonotoneConstantGaussianLeaf;
  const double gaps[] = {0, 2, 5, 6.5, 7, 7.5, 8, 10, 15, 30, 37, 40, 60, 200};
  const double ratios[] = {0.05, 0.2, 1, 5, 20, 50};
  double worst = 0.0;
  bool finite = true;
  for (double gap : gaps)
    for (double ratio : ratios) {
      double sL = 0.3, sR = ratio * sL, mL = 0.7;
      double mR = mL - gap * std::sqrt(sL * sL + sR * sR);
      double exact = Rf_pnorm5(-gap, 0.0, 1.0, 1, 1);
      double scale = std::max(1.0, std::fabs(exact));
      // a finite upper bound far past the integrand forces the integral
      double hiR = std::max(mL, mR) + 80.0 * std::max(sL, sR);
      double bounded =
        Leaf::logConeProbability(-HUGE_VAL, hiR, -HUGE_VAL, HUGE_VAL, mL, sL,
                                 mR, sR);
      Leaf::PairUpperLogDensity logDensity{mR, sR, -HUGE_VAL, HUGE_VAL, mL,
                                           sL};
      double whole = monotoneLogIntegrateLogConcave(
                       logDensity, -HUGE_VAL, HUGE_VAL, mR, std::min(sR, sL),
                       std::numeric_limits<double>::quiet_NaN()) -
                     std::log(sR) - 0.5 * std::log(2.0 * std::numbers::pi);
      double closed = Leaf::logConeProbability(-HUGE_VAL, HUGE_VAL, -HUGE_VAL,
                                               HUGE_VAL, mL, sL, mR, sR);
      finite = finite && std::isfinite(bounded) && std::isfinite(whole);
      worst = std::max(worst, std::fabs(bounded - exact) / scale);
      worst = std::max(worst, std::fabs(whole - exact) / scale);
      worst = std::max(worst, std::fabs(closed - exact) / scale);
    }
  // children whose posterior sds differ by 1e3 to 1e10, either way round, and
  // a captured pair 655 apart sitting 1.97 sd against the order
  struct Wide {
    double gap, sL, sR;
  };
  const Wide wides[] = {
    {-1.97, 0.00177, 1.159}, {0.0, 1e-3, 1.0}, {-1.0, 1e-3, 1.0},
    {3.0, 1e-3, 1.0}, {-8.0, 1e-6, 1.0}, {-2.0, 1e-6, 1.0},
    {-30.0, 1e-6, 1.0}, {-2.0, 1e-10, 1.0}, {1.0, 1e-10, 1.0},
    {-2.0, 1.0, 1e-3}, {-8.0, 1.0, 1e-6}, {-2.0, 1.0, 1e-10},
  };
  for (const Wide& w : wides) {
    double mL = 0.4, mR = mL + w.gap * std::sqrt(w.sL * w.sL + w.sR * w.sR);
    double exact = Rf_pnorm5(w.gap, 0.0, 1.0, 1, 1);
    double hiR = std::max(mL, mR) + 80.0 * std::max(w.sL, w.sR);
    double got = Leaf::logConeProbability(-HUGE_VAL, hiR, -HUGE_VAL, HUGE_VAL,
                                          mL, w.sL, mR, w.sR);
    finite = finite && std::isfinite(got);
    worst = std::max(worst, std::fabs(got - exact) / std::max(1.0, std::fabs(exact)));
  }
  check(finite, "monotone cone: finite at every gap");
  check(worst <= 1e-10, "monotone cone: the integral matches the closed form "
                        "to 1e-10 (relative past 1 nat)");

  // bounded cones: {zA, bL - aL, lowR - aL, hiR - aL, mR - aL, sR}, the lower
  // leaf standard normal and aL = zA
  struct Case {
    double zA, bL, lowR, hiR, mR, sR;
  };
  const Case cases[] = {
    {6.0, 3.0, 0.0, HUGE_VAL, -5.0, 1.0},
    {6.0, HUGE_VAL, 0.5, 10.0, 1.0, 0.2},
    {8.7, 3.0, 0.0, 10.0, -5.0, 5.0},
    {8.7, HUGE_VAL, 0.0, HUGE_VAL, -20.0, 1.0},
    {8.7, 0.5, 0.2, 4.0, 2.0, 0.05},
    {20.0, 3.0, 0.0, HUGE_VAL, -5.0, 0.2},
    {20.0, HUGE_VAL, 0.1, 10.0, -30.0, 3.0},
    {20.0, 0.01, 0.0, 1.0, 0.0, 1.0},
  };
  double worstBounded = 0.0;
  for (const Case& cs : cases) {
    double aL = cs.zA, bL = aL + cs.bL, lowR = aL + cs.lowR;
    double hiR = aL + cs.hiR, mR = aL + cs.mR;
    double got =
      Leaf::logConeProbability(lowR, hiR, aL, bL, 0.0, 1.0, mR, cs.sR);
    double ref = refLogCone(lowR, hiR, aL, bL, 0.0, 1.0, mR, cs.sR);
    finite = finite && std::isfinite(got);
    worstBounded =
      std::max(worstBounded, std::fabs(got - ref) / std::max(1.0, std::fabs(ref)));
  }
  // a cone captured from a fit whose frozen bound sat 8.7 sd above the lower
  // leaf: the lower-tail difference rounded the integrand to zero everywhere
  {
    double lowR = -0.054755, aL = -0.054755, mL = -0.129222, sL = 0.008561;
    double mR = -0.059441, sR = 0.019915;
    double got = Leaf::logConeProbability(lowR, HUGE_VAL, aL, HUGE_VAL, mL, sL,
                                          mR, sR);
    double ref = refLogCone(lowR, HUGE_VAL, aL, HUGE_VAL, mL, sL, mR, sR);
    check(std::isfinite(got), "monotone cone: a frozen bound far above the "
                              "lower leaf is not a false zero");
    checkNear(got, -41.875, 1e-2, "monotone cone: the captured frozen-bound "
                                  "cone");
    worstBounded = std::max(worstBounded, std::fabs(got - ref) / std::fabs(ref));
  }
  check(finite, "monotone cone: finite with bounds set");
  check(worstBounded <= 1e-8,
        "monotone cone: bounded cones match the log-space reference");
  check(Leaf::logConeProbability(1.0, 1.0, 0.0, HUGE_VAL, 0.0, 1.0, 0.0,
                                 1.0) == -HUGE_VAL &&
          Leaf::logConeProbability(0.0, 2.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0) ==
            -HUGE_VAL,
        "monotone cone: -Inf on an empty cone");

  // logStandardNormalMass on intervals one ulp wide: pnorm is not monotone to
  // the ulp, so a difference of tails can come out positive
  bool notNaN = true;
  double worstUlp = 0.0;
  for (int i = 0; i < 100000; ++i) {
    double lo = -8.0 + 16.0 * (i + 0.5) / 100000;
    double hi = std::nextafter(lo, HUGE_VAL);
    double got = logStandardNormalMass(lo, hi);
    if (!std::isfinite(got)) notNaN = false;
    double ref = std::log(hi - lo) - 0.5 * lo * lo -
                 0.5 * std::log(2.0 * std::numbers::pi);
    worstUlp = std::max(worstUlp, std::fabs(got - ref) / std::fabs(ref));
  }
  check(notNaN, "logStandardNormalMass: finite on an interval one ulp wide");
  check(worstUlp < 1e-12, "logStandardNormalMass: one ulp is width times "
                          "density");
  // far out, where the log tails are so large that a narrow interval's
  // difference of them is rounding: against the asymptotic series
  // log Q(x) = -x^2 / 2 - log x - log sqrt(2 pi) + log(1 - x^-2 + 3 x^-4 -
  // 15 x^-6), differenced analytically
  bool finiteFar = true;
  double worstFar = 0.0;
  auto series = [](double x) {
    double r = 1.0 / (x * x);
    return std::log1p(-r + 3.0 * r * r - 15.0 * r * r * r);
  };
  for (double lo : {1e3, 3.5e4, 3.5e5, 1e6, 1e7})
    for (double k : {0.5, 1.01, 2.0, 10.0, 1e2, 1e4, 1e6}) {
      double width = k * 1e-5 / (1.0 + lo), hi = lo + width;
      if (!(hi > lo)) continue;  // narrower than an ulp of lo
      width = hi - lo;
      double got = logStandardNormalMass(lo, hi);
      double d = -0.5 * width * (lo + hi) - std::log1p(width / lo) +
                 series(hi) - series(lo);
      double ref = Rf_pnorm5(lo, 0.0, 1.0, 0, 1) + std::log(-std::expm1(d));
      finiteFar = finiteFar && std::isfinite(got);
      worstFar = std::max(worstFar, std::fabs(got - ref) / std::fabs(ref));
    }
  check(finiteFar, "logStandardNormalMass: finite on narrow far-tail "
                   "intervals");
  check(worstFar < 1e-12, "logStandardNormalMass: narrow far-tail intervals "
                          "match the asymptotic series");
  printf("ok: monotone cone deep tail (worst %.1e closed form, %.1e bounded)\n",
         worst, worstBounded);
}

// The same scores through the move's seam, logLikelihoodForBranchWithParams,
// on hand-built trees whose data run against an increasing constraint: a
// two-leaf birth with unequal children at gaps 8 to 45 sd (closed form), a
// pair under a frozen neighbor 9 and 25 sd above the lower child (log-space
// reference), and one merged leaf under a frozen bound 45 sd above it.
static void testMonotoneSeamDeepTail() {
  const size_t n = 400;
  std::vector<double> x(n), y(n), weights(n, 1.0);
  for (size_t i = 0; i < n; ++i) x[i] = (static_cast<double>(i) + 0.5) / n;
  ColumnStore store;
  built(store.build(x.data(), n, 1, 40));
  const double sig2 = 0.01, scale = 0.5, k = 2.0;
  const double c = std::sqrt(std::numbers::pi / (std::numbers::pi - 1.0));
  const double tauC = c * scale / k;
  MonotoneConstantGaussianLeaf leaf;
  leaf.scale = scale;
  leaf.data = &store;
  leaf.directions = {1};
  leaf.cInflation = c;
  std::vector<index_t> idx(n);
  Tree tree;
  auto split = [&](std::int32_t node, std::int32_t cut) {
    Rule rule;
    rule.variableIndex = 0;
    rule.setSplitIndex(cut);
    tree.birth(store, node, rule, y.data(), weights.data());
    return tree.at(node).leftChild;
  };
  double worst = 0.0;
  bool finite = true;
  int checked = 0;

  // one split at x = 0.1: 40 rows below, 360 above sitting delta lower
  for (double delta : {0.25, 0.6, 1.4}) {
    for (size_t i = 0; i < n; ++i) y[i] = x[i] <= 0.1 ? 0.5 : 0.5 - delta;
    tree.initialize(idx.data(), n);
    tree.computeLeafStats(0, y.data(), weights.data());
    std::int32_t lower = split(0, 3), upper = lower + 1;
    const Node& nL = tree.at(lower);
    const Node& nR = tree.at(upper);
    double mL, sL, mR, sR;
    refPost(nL.sumWeights, nL.sumWeightedResponse, sig2, tauC, mL, sL);
    refPost(nR.sumWeights, nR.sumWeightedResponse, sig2, tauC, mR, sR);
    double ref = refBase(nL.sumWeights, nL.sumWeightedResponse, sig2, tauC) +
                 refBase(nR.sumWeights, nR.sumWeightedResponse, sig2, tauC) +
                 Rf_pnorm5((mR - mL) / std::sqrt(sL * sL + sR * sR), 0.0, 1.0,
                           1, 1);
    double got = leaf.logLikelihoodForBranchWithParams(tree, 0, y.data(),
                                                       weights.data(), k, sig2,
                                                       nullptr);
    finite = finite && std::isfinite(got);
    worst = std::max(worst, std::fabs(got - ref) / std::max(1.0, std::fabs(ref)));
    ++checked;
  }

  // A | B1 | B2 at x = 0.5 and 0.75, mu_A frozen z sd above B1's mean
  for (double z : {9.0, 25.0}) {
    for (size_t i = 0; i < n; ++i)
      y[i] = x[i] <= 0.5 ? 0.0 : (x[i] <= 0.75 ? -0.2 : 0.1);
    tree.initialize(idx.data(), n);
    tree.computeLeafStats(0, y.data(), weights.data());
    std::int32_t a = split(0, 19), b = a + 1;
    std::int32_t b1 = split(b, 29), b2 = b1 + 1;
    const Node& n1 = tree.at(b1);
    const Node& n2 = tree.at(b2);
    double m1, s1, m2, s2;
    refPost(n1.sumWeights, n1.sumWeightedResponse, sig2, tauC, m1, s1);
    refPost(n2.sumWeights, n2.sumWeightedResponse, sig2, tauC, m2, s2);
    std::vector<double> mu(tree.nodes.size(), 0.0);
    mu[a] = m1 + z * s1;
    double ref = refBase(n1.sumWeights, n1.sumWeightedResponse, sig2, tauC) +
                 refBase(n2.sumWeights, n2.sumWeightedResponse, sig2, tauC) +
                 refLogCone(mu[a], HUGE_VAL, mu[a], HUGE_VAL, m1, s1, m2, s2);
    double got = leaf.logLikelihoodForBranchWithParams(
      tree, b, y.data(), weights.data(), k, sig2, mu.data());
    finite = finite && std::isfinite(got);
    worst = std::max(worst, std::fabs(got - ref) / std::max(1.0, std::fabs(ref)));
    ++checked;

    // the merged leaf B under the same frozen bound, 45 sd above its mean
    tree.initialize(idx.data(), n);
    tree.computeLeafStats(0, y.data(), weights.data());
    a = split(0, 19);
    b = a + 1;
    const Node& nB = tree.at(b);
    double m, s;
    refPost(nB.sumWeights, nB.sumWeightedResponse, sig2, tauC, m, s);
    mu.assign(tree.nodes.size(), 0.0);
    mu[a] = m + 45.0 * s;
    ref = refBase(nB.sumWeights, nB.sumWeightedResponse, sig2, tauC) +
          Rf_pnorm5(45.0, 0.0, 1.0, 0, 1);
    got = leaf.logLikelihoodForBranchWithParams(tree, b, y.data(),
                                                weights.data(), k, sig2,
                                                mu.data());
    finite = finite && std::isfinite(got);
    worst = std::max(worst, std::fabs(got - ref) / std::max(1.0, std::fabs(ref)));
    ++checked;
  }
  check(finite, "monotone seam: contrary data and far frozen bounds score "
                "finite");
  check(worst <= 1e-8, "monotone seam: deep-tail scores match their "
                       "references");
  printf("ok: monotone seam deep tail (%d scores, worst %.1e)\n", checked,
         worst);
}

// The free bound decides a move only where the count would decide it the
// same way, for the same u: over random trees, every nog node, births and
// deaths, and a spread of the rest of the ratio.
static void testMonotoneFreeBound() {
  ColumnStore store;
  makeStore(store, 3, 12, 300);
  const std::int8_t dir[3] = {1, -1, 0};
  MonotoneConstantGaussianLeaf leaf;
  leaf.data = &store;
  leaf.directions.assign(dir, dir + 3);
  MonotoneCountScratch s;
  int moves = 0, mismatches = 0, decidedBirth = 0, decidedDeath = 0,
      countedBoth = 0;
  for (int trial = 0; trial < 150; ++trial) {
    TestTree t(store);
    t.growRandom(2 + trial % 8);
    std::vector<std::int32_t> internal;
    t.tree.fillNotBottom(0, internal);
    for (std::int32_t parent : internal) {
      if (!t.tree.childrenAreBottom(parent)) continue;
      double logRatio = monotoneLogNormalizerRatio(t.tree, store, dir, parent, s);
      for (int draw = 0; draw < 20; ++draw) {
        double r1 = std::exp(6.0 * (runif01() - 0.5)), u = runif01();
        for (bool birth : {true, false}) {
          double logBound = leaf.prepareLogNormalizerRatio(t.tree, parent);
          double ratio = r1, logNormalizer;
          bool accept = decideNormalizedMove(leaf, u, logBound, birth, &ratio,
                                             &logNormalizer);
          bool reference =
            u < r1 * std::exp(birth ? logRatio : -logRatio);
          mismatches += accept != reference;
          if (std::isnan(logNormalizer)) {
            (birth ? decidedBirth : decidedDeath)++;
          } else {
            ++countedBoth;
            if (std::fabs(logNormalizer - logRatio) > 1e-12) ++mismatches;
          }
          ++moves;
        }
      }
    }
  }
  check(mismatches == 0, "monotone free bound: same decision as the count");
  check(decidedBirth > 100 && decidedDeath > 100 && countedBoth > 100,
        "monotone free bound: births and deaths decided both ways");
  printf("ok: monotone free bound (%d decisions, %d births and %d deaths "
         "without a count)\n",
         moves, decidedBirth, decidedDeath);
}

// The "joint" prior's structure draw: one predictor on three values (two
// cuts), one tree, so the trees are the root, a split at either cut, and each
// of those with its other side split too. Their CGM probabilities times Z_T
// (1, 1/2, 1/6 for 1, 2, 3 chained leaves) are the joint law; the "leaf" prior
// draws CGM's own. Chi-square on 20000 draws each, and the other law refused.
static void testMonotoneJointTreePrior() {
  const size_t n = 90;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 3);
  const double base = 0.95, power = 2.0;
  double g1 = base / std::pow(2.0, power);
  // ROOT, cut 0.5 with 2 / 3 leaves, cut 1.5 with 2 / 3 leaves
  double cgm[5] = {1.0 - base, 0.5 * base * (1.0 - g1), 0.5 * base * g1,
                   0.5 * base * (1.0 - g1), 0.5 * base * g1};
  double z[5] = {1.0, 0.5, 1.0 / 6.0, 0.5, 1.0 / 6.0};
  std::vector<std::int8_t> dir = {1};
  std::vector<FlatNode> flat;
  std::vector<std::uint32_t> counts;
  double stat[2][2];
  for (int which = 0; which < 2; ++which) {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 8100 + which);
    SamplerOptions options;
    options.numTrees = 1;
    options.maxNumCuts = 2;
    options.base = base;
    options.power = power;
    options.birthOrDeathProbability = 1.0;
    options.swapProbability = 0.0;
    options.changeProbability = 0.0;
    options.monotoneDirections = dir.data();
    options.monotonePrior = static_cast<std::uint8_t>(which);
    Sampler<MonotoneConstantGaussianLeaf> sampler(
      x.data(), y.data(), n, 1, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &rng);
    double observed[5] = {0, 0, 0, 0, 0};
    const int draws = 20000;
    for (int d = 0; d < draws; ++d) {
      sampler.sampleTreesFromPrior();
      sampler.flattenTree(0, 0, flat, counts);
      if (flat[0].variable == invalidVariable) {
        observed[0] += 1.0;
        continue;
      }
      int leaves = 0;
      for (const FlatNode& node : flat) leaves += node.variable == invalidVariable;
      observed[(flat[0].value < 1.0 ? 1 : 3) + (leaves == 3)] += 1.0;
    }
    for (int law = 0; law < 2; ++law) {
      double total = 0.0, expected[5];
      for (int i = 0; i < 5; ++i) total += expected[i] = cgm[i] * (law ? z[i] : 1.0);
      double chi2 = 0.0;
      for (int i = 0; i < 5; ++i) {
        double e = draws * expected[i] / total;
        chi2 += (observed[i] - e) * (observed[i] - e) / e;
      }
      stat[which][law] = chi2;
    }
    ext_rng_destroy(rng);
  }
  // 4 degrees of freedom: 23.5 is p = 1e-4
  check(stat[0][0] < 23.5 && stat[0][1] > 200.0,
        "monotone leaf prior: trees follow the CGM prior");
  check(stat[1][1] < 23.5 && stat[1][0] > 200.0,
        "monotone joint prior: trees follow p_CGM(T) Z_T");
  printf("ok: monotone joint tree prior (chi-square %.1f leaf, %.1f joint)\n",
         stat[0][0], stat[1][1]);
}

// The joint prior's accept step keeps a tree exactly when the iid leaves it
// drew lie in the cone, read here off the point oracle's required pairs. The
// tree law alone cannot see a reversed test: negating iid draws of a common
// sd maps the cone onto its reverse, so both are accepted at the rate Z_T.
// Own rng; restores the runif01 stream the random trees consume.
static void testMonotoneJointAcceptIsConeTest() {
  std::uint64_t saved = rngState;
  ColumnStore store;
  makeStore(store, 3, 8, 200);
  const std::int8_t dir[3] = {1, -1, 0};
  MonotoneConstantGaussianLeaf leaf;
  leaf.scale = 1.0;
  leaf.data = &store;
  leaf.directions.assign(dir, dir + 3);
  leaf.cInflation = std::sqrt(std::numbers::pi / (std::numbers::pi - 1.0));
  leaf.prior = MonotonePrior::joint;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20261001u);
  int mismatches = 0, accepted = 0, rejected = 0;
  for (int trial = 0; trial < 60; ++trial) {
    TestTree t(store);
    t.growRandom(1 + trial % 5);
    LeafPairs pairs = pointOrder(t.tree, store, dir);
    for (int draw = 0; draw < 20; ++draw) {
      bool accept = leaf.jointPriorAccepts(rng, t.tree, 2.0);
      bool inCone = true;
      for (const auto& pair : pairs)
        inCone = inCone &&
                 leaf.jointDraw[pair.first] <= leaf.jointDraw[pair.second];
      mismatches += accept != inCone;
      (accept ? accepted : rejected)++;
    }
  }
  check(mismatches == 0,
        "monotone joint prior: a tree is kept exactly when its draws lie in "
        "the cone");
  check(accepted > 100 && rejected > 100,
        "monotone joint prior: draws kept and refused both");
  ext_rng_destroy(rng);
  rngState = saved;
  printf("ok: monotone joint accept is the cone test (%d kept, %d refused)\n",
         accepted, rejected);
}

namespace {
using MonotoneSampler = Sampler<MonotoneConstantGaussianLeaf>;

/// x1 constrained increasing, x2 free, y rising in x1 with an x2 bump, so
/// a leaf-prior fit keeps counting orders.
std::unique_ptr<MonotoneSampler> makeCountingSampler(
    std::vector<double>& x, std::vector<double>& y, const std::int8_t* dir,
    std::size_t numChains, std::size_t numThreads, ext_rng** rngs) {
  const size_t n = 300;
  x.resize(2 * n);
  y.resize(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[n + i] = runif01();
    y[i] = 2.0 * x[i] + std::sin(6.0 * x[n + i]) + 0.1 * (runif01() - 0.5);
  }
  SamplerOptions options;
  options.numTrees = 4;
  options.numChains = numChains;
  options.numThreads = numThreads;
  options.birthOrDeathProbability = 1.0;
  options.swapProbability = 0.0;
  options.changeProbability = 0.0;
  options.monotoneDirections = dir;
  return std::make_unique<MonotoneSampler>(
      x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, rngs);
}

/// Chain c's derived state against a from-scratch rebuild of its trees:
/// every row's leaf map entry is the leaf its tree routes it to, totalFits
/// is the tree-order gather, and every tree lies in the cone. Exact after a
/// rebuild; a sweep's difference updates round, which tolerance admits.
bool chainMatchesItsTrees(MonotoneSampler& sampler, std::size_t c,
                          const std::int8_t* dir, double tolerance = 0.0) {
  auto& chain = sampler.chain(c);
  const ColumnStore& data = sampler.data();
  std::size_t n = data.numObservations, numTrees = chain.numTrees();
  std::vector<double> gather(n, 0.0);
  bool ok = true;
  for (std::size_t t = 0; t < numTrees; ++t) {
    const Tree& tree = chain.tree(t);
    const std::uint32_t* leaf = TestPeer::leafOf(chain, t);
    std::vector<double>& mu = TestPeer::muByTree(chain, t);
    std::vector<std::int32_t> bottoms;
    tree.fillBottom(0, bottoms);
    for (std::int32_t b : bottoms)
      for (std::size_t m = tree.at(b).begin; m < tree.at(b).end; ++m)
        ok = ok && leaf[tree.indices[m]] == static_cast<std::uint32_t>(b);
    ok = ok && mu.size() >= tree.nodes.size() &&
         monotoneTreeIsFeasible(tree, data, dir, mu.data());
    for (std::size_t i = 0; i < n; ++i) gather[i] += mu[leaf[i]];
  }
  const std::vector<double>& total = TestPeer::totalFitsInForest(chain, 0);
  for (std::size_t i = 0; i < n; ++i)
    ok = ok && std::fabs(total[i] - gather[i]) <= tolerance;
  return ok;
}

std::vector<std::vector<FlatNode>> flattenChain(MonotoneSampler& sampler,
                                                std::size_t c) {
  std::vector<std::vector<FlatNode>> trees(sampler.chain(c).numTrees());
  std::vector<std::uint32_t> counts;
  for (std::size_t t = 0; t < trees.size(); ++t)
    sampler.flattenTree(c, t, trees[t], counts);
  return trees;
}

bool sameTree(const std::vector<FlatNode>& a, const std::vector<FlatNode>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a[i].variable != b[i].variable || a[i].mask != b[i].mask ||
        a[i].flags != b[i].flags)
      return false;
  return true;
}

/// Trees [0, t*) changed and [t*, T) did not, t* < T: a sweep stopped at
/// tree t*, whose move was put back.
bool stoppedPartWay(const std::vector<std::vector<FlatNode>>& before,
                    const std::vector<std::vector<FlatNode>>& after) {
  std::size_t t = 0;
  while (t < before.size() && !sameTree(before[t], after[t])) ++t;
  if (t == before.size()) return false;
  for (std::size_t u = t; u < before.size(); ++u)
    if (!sameTree(before[u], after[u])) return false;
  return true;
}
}  // namespace

/// Step 15: a cancel inside a count, inline and on another thread; an
/// allocation failure inside one, inline and rethrown from a worker after
/// the join; and the slow-count tally.
static void testMonotoneCountInterrupt() {
  MonotoneCountHooks& hooks = monotoneCountHooks();
  const std::int8_t dir[2] = {1, 0};
  std::vector<double> x, y;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20261015u);
  auto sampler = makeCountingSampler(x, y, dir, 1, 1, &rng);
  Results none;
  sampler->run(200, 0, none);

  // a cancel on the second poll: the first is the sweep's own, so a true
  // return means a count was stopped
  hooks.pollInterval.store(1);
  for (int onThread = 0; onThread < 2; ++onThread) {
    bool stopped = false;
    for (int attempt = 0; attempt < 200 && !stopped; ++attempt) {
      auto before = flattenChain(*sampler, 0);
      int calls = 0;
      std::function<bool()> cancel = [&calls]() { return ++calls >= 2; };
      Results empty;
      auto body = [&]() {
        stopped = sampler->chain(0).run(1, 0, empty, nullptr, 0, &cancel);
      };
      if (onThread) {
        std::thread worker(body);
        worker.join();
      } else {
        body();
      }
      if (!stopped) continue;
      auto after = flattenChain(*sampler, 0);
      const char* where = onThread ? "on a worker" : "inline";
      check(calls == 2, "monotone count cancel: stopped at the count's poll");
      check(stoppedPartWay(before, after),
            onThread ? "monotone count cancel on a worker: the moved tree is "
                       "back at T0, later trees untouched"
                     : "monotone count cancel inline: the moved tree is back "
                       "at T0, later trees untouched");
      check(chainMatchesItsTrees(*sampler, 0, dir),
            onThread ? "monotone count cancel on a worker: the fits equal a "
                       "rebuild from the trees"
                     : "monotone count cancel inline: the fits equal a "
                       "rebuild from the trees");
      (void) where;
    }
    check(stopped, "monotone count cancel: some sweep's count was stopped");
    sampler->run(5, 0, none);
    check(chainMatchesItsTrees(*sampler, 0, dir, 1e-9),
          "monotone count cancel: the next run is valid");
  }
  hooks.pollInterval.store(std::size_t(1) << 16);

  // an allocation failure in a count, inline
  hooks.failNextCount.store(true);
  bool threw = false;
  for (int attempt = 0; attempt < 200 && !threw; ++attempt) {
    try {
      sampler->run(1, 0, none);
    } catch (const std::bad_alloc& e) {
      threw = std::strstr(e.what(), "prior = \"joint\"") != nullptr;
    }
  }
  hooks.failNextCount.store(false);
  check(threw, "monotone count allocation failure: inline, the run throws a "
               "bad_alloc naming the remedies");
  check(chainMatchesItsTrees(*sampler, 0, dir),
        "monotone count allocation failure: the fits equal a rebuild");
  sampler->run(5, 0, none);
  check(chainMatchesItsTrees(*sampler, 0, dir, 1e-9),
        "monotone count allocation failure: the next run is valid");

  // the tally: every count over a negative threshold, and none over a huge
  // one, since each run starts it afresh
  hooks.slowSeconds.store(-1.0);
  sampler->run(20, 0, none);
  SlowCountTally tally = sampler->slowCountTally();
  check(tally.slowCounts > 0 && tally.slowestLeaves >= 2 &&
            tally.slowestDownSets > 0 && tally.slowestSeconds >= 0.0,
        "monotone slow counts: a lowered threshold tallies counts and sizes");
  hooks.slowSeconds.store(1e9);
  sampler->run(20, 0, none);
  check(sampler->slowCountTally().slowCounts == 0,
        "monotone slow counts: the tally resets at the next run");
  hooks.slowSeconds.store(1.0);
  ext_rng_destroy(rng);

  // workers: one chain's allocation failure stops the other and is rethrown
  // after the join. Run on a thread of its own with a deadline, so a join
  // that never returns fails the suite rather than hanging it.
  ext_rng* rngs[2];
  for (int c = 0; c < 2; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(rngs[c], 20261016u + c);
  }
  auto workers = makeCountingSampler(x, y, dir, 2, 2, rngs);
  workers->run(100, 0, none);
  hooks.failNextCount.store(true);
  std::promise<int> outcome;
  std::future<int> done = outcome.get_future();
  std::thread driver([&]() {
    int result = 0;
    try {
      for (int attempt = 0; attempt < 200 && result == 0; ++attempt)
        workers->run(1, 0, none);
    } catch (const std::bad_alloc&) {
      result = 1;
    } catch (...) {
      result = 2;
    }
    outcome.set_value(result);
  });
  if (done.wait_for(std::chrono::seconds(60)) != std::future_status::ready) {
    std::printf("FAIL: monotone count allocation failure on a worker: the "
                "run did not return (a chain's exception skipped the join's "
                "count)\n");
    std::fflush(stdout);
    std::_Exit(1);
  }
  driver.join();
  hooks.failNextCount.store(false);
  check(done.get() == 1, "monotone count allocation failure on a worker: "
                         "rethrown after the join");
  for (std::size_t c = 0; c < 2; ++c)
    check(chainMatchesItsTrees(*workers, c, dir, 1e-9),
          "monotone count allocation failure on a worker: every chain's "
          "fits match its trees");
  workers->run(5, 0, none);
  for (std::size_t c = 0; c < 2; ++c)
    check(chainMatchesItsTrees(*workers, c, dir, 1e-9),
          "monotone count allocation failure on a worker: the next run is "
          "valid");
  for (int c = 0; c < 2; ++c) ext_rng_destroy(rngs[c]);
  printf("ok: monotone count cancel, allocation failure and slow-count "
         "tally\n");
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
  testMonotoneMissingRelates();
  testMonotoneMoveClosedForm();
  testMonotoneConeDeepTail();
  testMonotoneSeamDeepTail();
  testMonotoneFreeBound();
  testMonotoneJointTreePrior();
  testMonotoneJointAcceptIsConeTest();
  testMonotoneCountInterrupt();
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
