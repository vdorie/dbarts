// Interaction-constraint tests (docs/design/interaction-constraints.md,
// deliverable B): the availability-predicate oracle, the swap sibling-strand
// de-risk toy, a change-strand feasibility invariant, and an enumerable
// exact-posterior structure match. Needs the data/tree/model/moves layers but
// not chain/sampler/facade, so a touch there does not force a recompile.
#include "assert.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

#include <external/random.h>

#include <bartcore/data.hpp>
#include <bartcore/tree.hpp>
#include <bartcore/model.hpp>
#include <bartcore/moves.hpp>

using namespace bartcore;

namespace {

// A p-bitset the reference math builds by hand, independent of the engine.
struct RefSet {
  std::vector<int> bits;
  explicit RefSet(size_t p) : bits(p, 0) {}
  void set(size_t j) { bits[j] = 1; }
  bool has(size_t j) const { return bits[j] != 0; }
  size_t count() const {
    size_t c = 0;
    for (int b : bits) c += static_cast<size_t>(b);
    return c;
  }
};

// Distinct split variables strictly above nodeIndex, by an independent walk.
RefSet referenceAncestors(const Tree& tree, int32_t nodeIndex, size_t p) {
  RefSet a(p);
  int32_t current = nodeIndex;
  while (tree.at(current).parent != invalidNode) {
    current = tree.at(current).parent;
    a.set(static_cast<size_t>(tree.at(current).rule.variableIndex));
  }
  return a;
}

// birth helper: split nodeIndex on an ordinal variable at a cut index.
void birthOrdinal(Tree& tree, const ColumnStore& store, int32_t nodeIndex,
                  int32_t variableIndex, int32_t cut, const double* y,
                  const double* weights) {
  Rule rule;
  rule.variableIndex = variableIndex;
  rule.setSplitIndex(cut);
  tree.birth(store, nodeIndex, rule, y, weights);
}

// ---------------------------------------------------------------------------
// (1) availability oracle: the constrained collectAvailableVariables and
// per-variable variableAvailable must equal an independent brute-force
// reference (cut-availability AND the interaction predicate) at every node.
// ---------------------------------------------------------------------------
void testAvailabilityOracle() {
  const size_t n = 200, p = 5;
  std::vector<double> x(n * p);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < n; ++i)
      x[i + j * n] = static_cast<double>((i * 7 + j * 3) % 11) / 11.0;
  std::vector<double> y(n, 0.0);
  for (size_t i = 0; i < n; ++i) y[i] = static_cast<double>(i % 3);

  ColumnStore store;
  store.build(x.data(), n, p, 10);  // ~10 cuts/column: cuts do not exhaust

  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);

  // a fixed, multi-level tree touching several variables on some paths
  birthOrdinal(tree, store, 0, 0, 4, y.data(), nullptr);        // root: x0
  int32_t l = tree.at(0).leftChild, r = l + 1;
  birthOrdinal(tree, store, l, 2, 3, y.data(), nullptr);        // left: x2
  birthOrdinal(tree, store, r, 1, 5, y.data(), nullptr);        // right: x1
  int32_t ll = tree.at(l).leftChild;
  birthOrdinal(tree, store, ll, 4, 2, y.data(), nullptr);       // left-left: x4

  std::vector<int32_t> nodesToCheck;
  tree.fillSubtree(0, nodesToCheck);

  // three constraint flavors: order-only, forbidden-only, and both
  struct Flavor { size_t maxOrder; std::vector<size_t> pairs; };
  Flavor flavors[] = {
    {3, {}},
    {0, {1, 2, 0, 4}},
    {2, {2, 4, 1, 3}},
  };

  bool ok = true;
  std::vector<std::uint8_t> avail(p);
  for (const Flavor& f : flavors) {
    InteractionConstraint constraint;
    constraint.build(p, f.maxOrder, f.pairs.empty() ? nullptr : f.pairs.data(),
                     f.pairs.size() / 2);

    for (int32_t node : nodesToCheck) {
      // cut-availability from the UNCONSTRAINED engine (constraint lifted)
      tree.setInteractionConstraint(nullptr);
      std::vector<int> cutAvail(p);
      for (size_t j = 0; j < p; ++j)
        cutAvail[j] =
          tree.variableAvailable(store, node, static_cast<int32_t>(j)) ? 1 : 0;

      // reference interaction predicate applied by hand
      RefSet ancestors = referenceAncestors(tree, node, p);
      size_t order = ancestors.count();
      std::vector<int> ref(p);
      for (size_t j = 0; j < p; ++j) {
        bool allowed = cutAvail[j] != 0;
        if (allowed && f.maxOrder > 0 && f.maxOrder < p)
          if (!ancestors.has(j) && order >= f.maxOrder) allowed = false;
        if (allowed)
          for (size_t k = 0; k < f.pairs.size(); k += 2) {
            size_t a = f.pairs[k], b = f.pairs[k + 1];
            if ((j == a && ancestors.has(b)) || (j == b && ancestors.has(a)))
              allowed = false;
          }
        ref[j] = allowed ? 1 : 0;
      }

      // the two constrained engine paths must both match the reference
      tree.setInteractionConstraint(&constraint);
      size_t count = tree.collectAvailableVariables(store, node, avail.data());
      size_t refCount = 0;
      for (size_t j = 0; j < p; ++j) {
        refCount += static_cast<size_t>(ref[j]);
        if (avail[j] != (ref[j] ? 1 : 0)) ok = false;
        bool one =
          tree.variableAvailable(store, node, static_cast<int32_t>(j));
        if (one != (ref[j] != 0)) ok = false;
      }
      if (count != refCount) ok = false;
    }
  }
  tree.setInteractionConstraint(nullptr);
  check(ok, "availability oracle: collect + per-variable match brute force");
  printf("ok: interaction availability oracle\n");
}

// ---------------------------------------------------------------------------
// (2) DE-RISK: the swap sibling-strand toy. Forbid (x1, x2); root=x0, left=x1,
// right=x2 (0-based; the memo's forbid (x2,x3), root x1, children x2/x3). A
// swap that lifts x1 above x2 co-occurs the forbidden pair on the x2 sibling
// path - with NEITHER swapped variable equal to the stranded x2, so the
// per-variable ruleIsValid checks miss it and only the whole-subtree walk
// catches it. Every swap from this tree must be rejected.
// ---------------------------------------------------------------------------
void testSwapSiblingStrand(ext_rng* rng) {
  const size_t n = 240, p = 3;
  std::vector<double> x(n * p), y(n), weights(n, 1.0);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 2);          // x0
    x[i + n] = static_cast<double>((i / 2) % 2);  // x1
    x[i + 2 * n] = static_cast<double>((i / 4) % 2);  // x2
    y[i] = static_cast<double>((i % 2) + (i / 4) % 2);
  }
  ColumnStore store;
  store.build(x.data(), n, p, 1);  // 2 values/column -> 1 cut each

  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), weights.data());

  // root x0; left child x1, right child x2, each with two leaves below
  birthOrdinal(tree, store, 0, 0, 0, y.data(), weights.data());
  int32_t left = tree.at(0).leftChild, right = left + 1;
  birthOrdinal(tree, store, left, 1, 0, y.data(), weights.data());
  birthOrdinal(tree, store, right, 2, 0, y.data(), weights.data());

  size_t pair[] = {1, 2};
  InteractionConstraint constraint;
  constraint.build(p, 0, pair, 1);  // forbid (x1, x2), no order cap
  tree.setInteractionConstraint(&constraint);

  // the original tree is feasible; both candidate swaps strand a sibling
  check(tree.interactionSubtreeIsValid(0),
        "swap-strand: original tree is interaction-feasible");

  // direct check: lifting x1 to the root strands the x2 sibling
  Rule rootRule = tree.at(0).rule, leftRule = tree.at(left).rule;
  tree.at(0).rule = leftRule;   // root <- x1
  tree.at(left).rule = rootRule;  // left <- x0
  check(!tree.interactionSubtreeIsValid(0),
        "swap-strand: x1 at root strands the x2 sibling (rejected)");
  tree.at(0).rule = rootRule;
  tree.at(left).rule = leftRule;

  // integration: the only swappable node is the root, and BOTH child choices
  // strand a sibling, so every swap must be a no-op that leaves the tree fixed
  CGMTreePrior prior;
  prior.base = 0.95;
  prior.power = 2.0;
  ConstantGaussianLeaf leaf{0.5};
  MoveScratch scratch;
  MoveContext ctx{store, prior, 0.5, 1.0, 0.0, 0.5, weights.data(), 2.0, scratch};

  bool everTaken = false, stayedValid = true;
  for (int iter = 0; iter < 500; ++iter) {
    bool stepTaken = false;
    swapMove(ctx, leaf, rng, tree, y.data(), 0.5, &stepTaken);
    everTaken |= stepTaken;
    stayedValid &= tree.interactionSubtreeIsValid(0);
    // structure must be unchanged: root x0, children x1 / x2
    stayedValid &= tree.at(0).rule.variableIndex == 0 &&
                   tree.at(left).rule.variableIndex == 1 &&
                   tree.at(right).rule.variableIndex == 2;
  }
  check(!everTaken, "swap-strand: every stranding swap rejected");
  check(stayedValid, "swap-strand: tree stays feasible and unchanged");

  tree.setInteractionConstraint(nullptr);
  printf("ok: swap sibling-strand de-risk\n");
}

// ---------------------------------------------------------------------------
// (3) change-strand feasibility invariant. Forbid (x0, x2); a change at an
// ancestor that redraws it to x2 would strand a descendant x0. Assert the
// direct rejection and that a long move chain never leaves an infeasible tree.
// ---------------------------------------------------------------------------
void testChangeStrandInvariant(ext_rng* rng) {
  const size_t n = 300, p = 3;
  std::vector<double> x(n * p), y(n), weights(n, 1.0);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 4) / 4.0;          // x0, 3 cuts
    x[i + n] = static_cast<double>((i / 4) % 4) / 4.0;  // x1
    x[i + 2 * n] = static_cast<double>((i / 16) % 4) / 4.0;  // x2
    y[i] = static_cast<double>(i % 4) + 0.3 * static_cast<double>((i / 4) % 4);
  }
  ColumnStore store;
  store.build(x.data(), n, p, 3);

  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), weights.data());

  // root x1, left child x0 (depth 2): changing the root to x2 strands x0
  birthOrdinal(tree, store, 0, 1, 1, y.data(), weights.data());
  int32_t left = tree.at(0).leftChild;
  birthOrdinal(tree, store, left, 0, 1, y.data(), weights.data());

  size_t pair[] = {0, 2};
  InteractionConstraint constraint;
  constraint.build(p, 0, pair, 1);  // forbid (x0, x2)
  tree.setInteractionConstraint(&constraint);

  check(tree.interactionSubtreeIsValid(0),
        "change-strand: original tree feasible");
  Rule saved = tree.at(0).rule;
  Rule rule2;
  rule2.variableIndex = 2;
  rule2.setSplitIndex(1);
  tree.at(0).rule = rule2;  // root <- x2
  check(!tree.interactionSubtreeIsValid(0),
        "change-strand: x2 at root strands the x0 descendant");
  tree.at(0).rule = saved;

  // feasibility invariant: drive the full move set; no accepted move may leave
  // an infeasible tree
  CGMTreePrior prior;
  prior.base = 0.95;
  prior.power = 2.0;
  ConstantGaussianLeaf leaf{0.5};
  MoveScratch scratch;
  MoveContext ctx{store, prior, 0.5, 0.1, 0.4, 0.5, weights.data(), 2.0, scratch};

  bool alwaysValid = true;
  for (int iter = 0; iter < 20000; ++iter) {
    bool stepTaken = false;
    StepType stepType;
    metropolisJumpForTree(ctx, leaf, rng, tree, y.data(), 0.5, &stepTaken,
                          &stepType);
    alwaysValid &= tree.interactionSubtreeIsValid(0);
  }
  check(alwaysValid, "change-strand: feasibility invariant after every move");

  tree.setInteractionConstraint(nullptr);
  printf("ok: change-strand feasibility invariant\n");
}

// classify a tree over 2 one-cut predictors into one of the 3 legal
// structures under max.order = 1: 0 = root leaf, 1 = split x0, 2 = split x1;
// -1 = anything else (an inert constraint would produce depth-2 / mixed trees).
int classifyStructure(const Tree& tree) {
  if (tree.at(0).isBottom()) return 0;
  int32_t v = tree.at(0).rule.variableIndex;
  int32_t l = tree.at(0).leftChild;
  if (!tree.at(l).isBottom() || !tree.at(l + 1).isBottom()) return -1;
  if (v == 0) return 1;
  if (v == 1) return 2;
  return -1;
}

// ---------------------------------------------------------------------------
// (4) exact-posterior gate. p = 2, one cut each, max.order = 1: the ONLY legal
// structures are {root leaf, split x0, split x1}. The constrained tree prior
// is pi = (1 - g0, g0/2, g0/2); the marginal is the ordinary conjugate
// Gaussian one (leaf model untouched). The birth/death/change kernel's
// structure-visit frequencies must match pi * m(y) / Z within MC error, and
// no illegal (depth-2 / mixed) structure may ever appear.
// ---------------------------------------------------------------------------
void testExactPosterior(ext_rng* rng) {
  const size_t n = 40, p = 2;
  // unit weights via the null-weights path (mathematically w = 1)
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 2);          // x0: alternating 0/1
    x[i + n] = static_cast<double>((i / 5) % 2);  // x1: blocks of five
    // signal on the leaf-prior scale (sd scale/k = 0.25) so all three
    // structures compete: a moderate x0 effect, a weak x1 effect, and fixed
    // deterministic within-leaf scatter
    y[i] = 0.14 * (2.0 * x[i] - 1.0) + 0.06 * (2.0 * x[i + n] - 1.0) +
           0.04 * (static_cast<double>(i % 7) - 3.0);
  }
  ColumnStore store;
  store.build(x.data(), n, p, 1);
  check(store.numCuts[0] == 1 && store.numCuts[1] == 1,
        "exact posterior: one cut per predictor");

  const double base = 0.6, power = 2.0, k = 2.0, sigma = 0.5, scale = 0.5;
  const double sigmaSq = sigma * sigma;

  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);

  ConstantGaussianLeaf leaf{scale};

  // marginal m(y | T) for each structure, from the conjugate leaf likelihood
  auto leafIL = [&](int32_t node) {
    return leaf.logIntegratedLikelihood(k, sigmaSq, tree.at(node).sumWeights,
                                        tree.at(node).sumWeightedResponse);
  };
  double logMarginal[3];
  logMarginal[0] = leafIL(0);  // root leaf
  for (int v = 0; v < 2; ++v) {
    birthOrdinal(tree, store, 0, v, 0, y.data(), nullptr);
    int32_t lc = tree.at(0).leftChild;
    logMarginal[1 + v] = leafIL(lc) + leafIL(lc + 1);
    tree.orphanChildren(0);
    tree.computeLeafStats(0, y.data(), nullptr);
  }

  // constrained tree prior: pi = (1 - g0, g0/2, g0/2), g0 = base (depth 0)
  double logPrior[3] = {std::log(1.0 - base), std::log(base / 2.0),
                        std::log(base / 2.0)};
  double logTarget[3], maxLog = -1e300;
  for (int s = 0; s < 3; ++s) {
    logTarget[s] = logPrior[s] + logMarginal[s];
    if (logTarget[s] > maxLog) maxLog = logTarget[s];
  }
  double target[3], norm = 0.0;
  for (int s = 0; s < 3; ++s) {
    target[s] = std::exp(logTarget[s] - maxLog);
    norm += target[s];
  }
  for (int s = 0; s < 3; ++s) target[s] /= norm;

  // run the constant-leaf structure kernel under the constraint and tally
  InteractionConstraint constraint;
  constraint.build(p, 1, nullptr, 0);  // max.order = 1
  tree.setInteractionConstraint(&constraint);

  CGMTreePrior prior;
  prior.base = base;
  prior.power = power;
  MoveScratch scratch;
  MoveContext ctx{store, prior, 0.5, 0.1, 0.4, 0.5, nullptr, k, scratch};

  const int numBurn = 5000, numDraws = 1500000;
  for (int iter = 0; iter < numBurn; ++iter) {
    bool stepTaken = false;
    StepType stepType;
    metropolisJumpForTree(ctx, leaf, rng, tree, y.data(), sigma, &stepTaken,
                          &stepType);
  }
  long counts[3] = {0, 0, 0};
  bool alwaysLegal = true;
  for (int iter = 0; iter < numDraws; ++iter) {
    bool stepTaken = false;
    StepType stepType;
    metropolisJumpForTree(ctx, leaf, rng, tree, y.data(), sigma, &stepTaken,
                          &stepType);
    int s = classifyStructure(tree);
    if (s < 0) { alwaysLegal = false; break; }
    ++counts[s];
  }
  check(alwaysLegal, "exact posterior: only the 3 legal structures ever visited");

  bool matches = true;
  for (int s = 0; s < 3; ++s) {
    double freq = static_cast<double>(counts[s]) / static_cast<double>(numDraws);
    if (std::fabs(freq - target[s]) > 0.01) matches = false;
  }
  check(matches, "exact posterior: structure-visit frequencies match pi * m / Z");
  printf("ok: interaction exact posterior "
         "(target %.3f/%.3f/%.3f, freq %.3f/%.3f/%.3f)\n",
         target[0], target[1], target[2],
         static_cast<double>(counts[0]) / numDraws,
         static_cast<double>(counts[1]) / numDraws,
         static_cast<double>(counts[2]) / numDraws);
  tree.setInteractionConstraint(nullptr);
}

}  // namespace

void runInteractionTests(ext_rng* rng) {
  testAvailabilityOracle();
  testSwapSiblingStrand(rng);
  testChangeStrandInvariant(rng);
  testExactPosterior(rng);
}
