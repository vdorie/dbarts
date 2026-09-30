#include "common.hpp"

#include <chrono>
#include <map>

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
      std::int32_t left, right;
      tree.splitInterval(*store, leaf, var, &left, &right);
      if (right < left) continue;
      split(leaf, var,
            left + static_cast<std::int32_t>(runif01() * (right - left + 1)));
    }
  }
};

void makeStore(ColumnStore& store, size_t p, std::uint32_t cuts, size_t n) {
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();
  built(store.build(x.data(), n, p, cuts));
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

// 200 random trees over 1-3 constrained axes with mixed directions: log Z_T
// against the brute-force e / L!.
static void testMonotoneCountRandom() {
  ColumnStore store;
  makeStore(store, 3, 12, 300);
  MonotoneCountScratch s;
  const std::int8_t dirs[5][3] = {
      {1, 0, 0}, {-1, 1, 0}, {1, -1, 1}, {-1, -1, -1}, {0, -1, 1}};
  double worst = 0.0;
  size_t largest = 0;
  for (int trial = 0; trial < 200; ++trial) {
    TestTree t(store);
    t.growRandom(2 + trial % 5);
    const std::int8_t* dir = dirs[trial % 5];
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
  printf("ok: monotone count, 200 random trees (worst |log Z err| %.2g)\n",
         worst);
}

// The ratio on T*'s side against two direct whole-tree counts over every
// death (both children leaves) of random trees, one component and two, and
// against brute-force theta on the smaller ones.
static void testMonotoneRatioRandom() {
  ColumnStore store;
  makeStore(store, 3, 12, 300);
  MonotoneCountScratch s, direct;
  const std::int8_t dirs[5][3] = {
      {1, 0, 0}, {-1, 0, 1}, {1, -1, 0}, {-1, 0, 0}, {0, 1, -1}};
  int one = 0, two = 0, decreasingPairs = 0, bruteChecked = 0;
  double worst = 0.0, worstTheta = 0.0;
  for (int trial = 0; trial < 200; ++trial) {
    TestTree t(store);
    t.growRandom(3 + trial % 9);
    const std::int8_t* dir = dirs[trial % 5];
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
  printf("ok: monotone ratio, %d one-component and %d two-component moves "
         "(%d decreasing, %d brute-forced)\n",
         one, two, decreasingPairs, bruteChecked);
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

void runMonotoneTests() {
  std::uint64_t saved = rngState;
  rngState = 7071u;
  testMonotoneCountHandBuilt();
  testMonotoneCountRandom();
  testMonotoneRatioRandom();
  testMonotoneRatioNamed();
  testMonotonePositionLaw();
  testMonotoneExtensionDraw();
  testMonotoneCountScale();
  rngState = saved;
}
