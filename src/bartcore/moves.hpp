#ifndef BARTCORE_MOVES_HPP
#define BARTCORE_MOVES_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#ifdef BARTCORE_MOVE_CENSUS
#include <cstdio>
#include <cstdlib>
#endif

#include <external/random.h>

#include "data.hpp"
#include "model.hpp"
#include "tree.hpp"

#ifdef BARTCORE_MOVE_CENSUS
#include "scan.hpp"
#endif

// The conjugate Metropolis-Hastings tree moves: birth/death, change, swap and
// perturb proposals with their acceptance ratios.

namespace bartcore {

/// Per-chain scratch for the moves; reused across calls to avoid allocation.
struct MoveScratch {
  std::vector<int32_t> nodeScratch;
  std::vector<double> probabilityScratch;
  Tree::SubtreeSnapshot snapshot;
  // pooled categorical masks: the changed node's reachable set, the pattern
  // draw, and a depth-indexed arena for the validity walk's left-side masks
  std::vector<std::uint64_t> reachableWords;
  std::vector<std::uint64_t> patternWords;
  std::vector<std::uint64_t> maskArena;
};

struct MoveContext {
  const ColumnStore& data;
  const CGMTreePrior& treePrior;
  double birthOrDeathProbability;
  double swapProbability;
  double perturbProbability;
  double birthProbability;
  const double* weights;
  double k;
  MoveScratch& scratch;
  // the tree's current leaf-parameter vector M, handed to the scoring of a
  // ParamScoringLeafModel (monotone reads frozen neighbor mu here); null and
  // unread for the conjugate leaves, which integrate every leaf out.
  const double* leafParams = nullptr;
};

/// A branch's score: the veto rank the moves compare FIRST, then the
/// log-likelihood of the leaves that rank admits.
///
/// rank is the worst of the branch's leaves under Tree::leafVetoRank (2 a
/// member-empty leaf, 1 a leaf of only zero-weight rows, 0 all admissible),
/// and logLikelihood sums only the rank-0 leaves. A vetoed leaf enters no
/// likelihood term of its forest, so its marginal is not part of the branch's
/// score; the conjugate leaves return exactly 0.0 there anyway, but a leaf
/// model with a prior normalization (linear, GP) does not, and paying its
/// factorization for a leaf that has nothing to estimate is both wrong and
/// wasteful. Skipping makes an all-vetoed branch score prior x transition
/// exactly.
struct BranchScore {
  int rank;
  double logLikelihood;
};

template <MoveScorableLeafModel L, typename ResidT = double>
BranchScore logLikelihoodForBranch(const MoveContext& ctx, const L& leaf,
                                   Tree& tree, int32_t branchIndex,
                                   const ResidT* y, double sigma) {
  std::vector<int32_t>& bottoms(tree.bottomScratch);
  bottoms.clear();
  tree.fillBottom(branchIndex, bottoms);

  // The rank is leaf-model independent - it reads membership and weight only -
  // so it is taken here for every leaf model, the branch-owning ones included.
  int rank = 0;
  double result = 0.0;
  for (int32_t i : bottoms) {
    int leafRank = tree.leafVetoRank(i, ctx.weights);
    if (leafRank > rank) rank = leafRank;
    // a leaf whose score reads M owns the branch marginal outright (the
    // constrained joint over the touched leaves given frozen neighbors); the
    // conjugate leaves take the per-leaf marginal sum.
    if constexpr (!ParamScoringLeafModel<L>)
      if (leafRank == 0)
        result += leaf.logIntegratedLikelihoodForNode(tree, y, ctx.weights,
                                                      ctx.k, sigma * sigma, i);
  }

  if constexpr (ParamScoringLeafModel<L>)
    return {rank, leaf.logLikelihoodForBranchWithParams(tree, branchIndex, y,
                                                        ctx.weights, ctx.k,
                                                        sigma * sigma,
                                                        ctx.leafParams)};
  return {rank, result};
}

/// Resolves a (current, proposal) score pair into the two log-likelihoods the
/// move's acceptance expression consumes, applying the veto lexicographically:
/// the branch of worse rank takes -HUGE_VAL and the comparison is decided
/// there, ranks equal it is today's arithmetic on the finite parts.
///
/// This keeps empty leaves out of the chain state at every scale - a valid
/// branch's log-likelihood is
/// unbounded below, so a finite penalty would be out-penalized by a big node
/// or a small sigma and the empty leaf accepted - while leaving a chain whose
/// CURRENT state is vetoed a law to move under. Weights do not ride the tree,
/// so any weight or mask install can strand a leaf that was fine when it was
/// grown; comparing two vetoed branches by their penalties alone gave NaN,
/// which rejects every proposal and freezes the forest permanently. Ranked,
/// such a chain mixes under prior x transition at constant likelihood and any
/// move clearing the veto is accepted outright, so the promise that a forest
/// with nothing to fit "sits at its prior" holds. Rank 2 keeps the membership
/// law absolute: a weight-vetoed state can still never install a member-empty
/// leaf, which state export and the predictor surface both require.
inline void resolveVetoRank(const BranchScore& current,
                            const BranchScore& proposal,
                            double* currentLogLikelihood,
                            double* proposalLogLikelihood) {
  *currentLogLikelihood =
    current.rank > proposal.rank ? -HUGE_VAL : current.logLikelihood;
  *proposalLogLikelihood =
    proposal.rank > current.rank ? -HUGE_VAL : proposal.logLikelihood;
  // -HUGE_VAL on both sides is the leaf model's own FEASIBILITY sentinel (the
  // monotone empty cone), not the veto; the difference would be NaN, which
  // rejects by comparison but reports an acceptance probability of 1. Reject
  // explicitly instead.
  if (*currentLogLikelihood == -HUGE_VAL && *proposalLogLikelihood == -HUGE_VAL)
    *currentLogLikelihood = 0.0;
}

#ifdef BARTCORE_MOVE_CENSUS
// ===========================================================================
// Stage 0 move census: SCAFFOLDING, not part of the sampler.
//
// Compiled only under -DBARTCORE_MOVE_CENSUS. With the macro unset every hook
// expands to nothing and its arguments go unevaluated, so the ordinary build
// carries no code and no cost. Nothing here draws and the probe restores every
// byte it touches, so both builds walk the same RNG stream.
//
// One comma-separated line per structural proposal is appended to the file
// named by BARTCORE_MOVE_CENSUS_FILE; with the variable unset nothing is
// written. Six record kinds, told apart by the first field:
//
//   p,sweep,forest,tree,move,noop,accepted,nodeDepth,treeDepth,interior,
//     logLikelihood,logPrior,logCorrection
//   d,sweep,forest,tree,nodeDepth,displacement,logRatio
//   g,sweep,forest,tree,node,nodeDepth,isNog,interior,nog,scanned,
//     jointCandidates,jointEntropy,jointIncumbent,jointMaximum,jointRank,
//     cutCandidates,cutEntropy,cutIncumbent,cutMaximum,cutRank
//   x,sweep,forest,tree,candidates,entropy,pickWeight,pickRank,maxWeight
//   r,sweep,forest,tree,node,current,target,accepted
//   t,sweep,forest,tree,leaves,interior,nog
//
// A 'p' record's three log terms are that move's own acceptance expression:
// the veto-resolved log-likelihood difference, the log prior ratio (birth and
// death the growth factors, change and perturb the subtree strictly below the
// node, swap the swapped subtree), and the surviving proposal-density ratio
// (birth and death the transition ratio, change and perturb
// logProposalCorrection, swap 0). All three
// are NA when the proposal never reached a score (noop = 1: pi(T') = 0, an
// unsatisfiable rule draw, or no eligible node), where nodeDepth is -1 if it
// had no target node. treeDepth and interior are the shape the proposal saw,
// which for an accepted birth or death is not the shape it left.
//
// A 'd' record prices the same-variable cut move at a fixed schedule of signed
// displacements of the change proposal's target node, clipped to that node's
// descendant-valid interval and deduplicated after clipping. Same variable, so
// the correction of the move being priced is identically 1 and the log ratio
// is the subtree-below prior difference plus the resolved likelihood
// difference.
//
// A 'g' record rides the same change proposal and prices the CLOSED
// neighbourhood a collapsed rule draw at a nog node would have: every
// (available ordinal variable, admissible cut) pair, weighted by the cut
// scan's collapsed marginal times the prior factors that survive - the node's
// own rule prior 1/|SI| (its split-variable factor is uniform over the
// available set and cancels, and both children being leaves the good set IS
// the ancestor interval, so no proposal count enters) and the two
// log(1 - growth(child)) terms that are exactly changeMove's below-node prior.
// jointRank is the incumbent's position by weight, 1 being the largest. The
// cut* fields repeat the summary restricted to the incumbent variable.
// scanned is 0 - and the ten summary fields NA - when the node is not a nog
// node, when any variable available there is categorical, or when the leaf
// model carries no scalar (sum w, sum wz) marginal for the scan to score.
//
// An 'x' record prices informed death: at every death proposal the nog nodes
// are weighted by exp of the merged-leaf marginal ratio, the merged leaf's
// statistic formed as the two children's sum (computeLeafStats re-accumulates
// over a node's index span instead, which needs a pass this probe does not
// take), and the realized uniform pick is located in that distribution.
//
// An 'r' record carries the perturb proposal's node and its signed
// displacement, so a run of same-direction accepted displacements at one node
// can be counted offline. Node ids are arena slots a released pair can reuse,
// so a run is only as long as the tree's own identity holds.
//
// A 't' record is per tree per sweep, written after the tree's move settles.
//
// The location, the shape and the probe's snapshot are per-thread singletons,
// so a threaded run is correct but writes every chain's records to the one
// file interleaved: run the census one chain at a time.
inline void findGoodOrdinalRules(const MoveContext& ctx, const Tree& tree,
                                 int32_t nodeIndex, int32_t variableIndex,
                                 int32_t* lower, int32_t* upper);

namespace census {

inline std::FILE* stream() {
  static std::FILE* file = []() {
    const char* path = std::getenv("BARTCORE_MOVE_CENSUS_FILE");
    return path != nullptr ? std::fopen(path, "a") : nullptr;
  }();
  return file;
}

/// Where the current proposal sits and the shape it saw. The sweep loop sets
/// the location; the hooks stash the shape before the move disturbs it.
struct State {
  long sweep = -1;
  int forest = -1;
  int tree = -1;
  int nodeDepth = -1;
  int treeDepth = 0;
  int interior = 0;
};

inline State& state() {
  static thread_local State s;
  return s;
}

inline void setLocation(long sweep, int forest, int tree) {
  State& s = state();
  s.sweep = sweep;
  s.forest = forest;
  s.tree = tree;
}

inline void walk(const Tree& tree, int32_t i, int depth, State& s) {
  if (depth > s.treeDepth) s.treeDepth = depth;
  if (tree.at(i).isBottom()) return;
  ++s.interior;
  walk(tree, tree.at(i).leftChild, depth + 1, s);
  walk(tree, tree.at(i).leftChild + 1, depth + 1, s);
}

inline void shape(const Tree& tree, int32_t node) {
  State& s = state();
  s.treeDepth = 0;
  s.interior = 0;
  s.nodeDepth = node == invalidNode ? -1 : static_cast<int>(tree.depthOf(node));
  walk(tree, 0, 0, s);
}

/// R-readable numerics: NA for a term the proposal never had, Inf/-Inf for the
/// veto's sentinel, %.17g otherwise.
inline const char* number(double x, char* buffer, std::size_t size) {
  if (std::isnan(x)) return "NA";
  if (std::isinf(x)) return x > 0.0 ? "Inf" : "-Inf";
  std::snprintf(buffer, size, "%.17g", x);
  return buffer;
}

inline void proposal(const char* move, bool noop, bool accepted,
                     double logLikelihood, double logPrior,
                     double logCorrection) {
  std::FILE* file = stream();
  if (file == nullptr) return;
  const State& s = state();
  char b[3][32];
  std::fprintf(file, "p,%ld,%d,%d,%s,%d,%d,%d,%d,%d,%s,%s,%s\n", s.sweep,
               s.forest, s.tree, move, noop ? 1 : 0, accepted ? 1 : 0,
               s.nodeDepth, s.treeDepth, s.interior,
               number(logLikelihood, b[0], sizeof(b[0])),
               number(logPrior, b[1], sizeof(b[1])),
               number(logCorrection, b[2], sizeof(b[2])));
}

/// A proposal that never reached a score.
inline void noop(const char* move, const Tree& tree, int32_t node) {
  shape(tree, node);
  double na = std::nan("");
  proposal(move, true, false, na, na, na);
}

template <MoveScorableLeafModel L, typename ResidT>
void cutProbe(const MoveContext& ctx, const L& leaf, Tree& tree, int32_t node,
              const ResidT* y, double sigma) {
  std::FILE* file = stream();
  if (file == nullptr) return;
  Rule rule = tree.at(node).rule;
  if (ctx.data.splitsBySubset(static_cast<std::size_t>(rule.variableIndex)))
    return;  // a first cut move is ordinal-only

  int32_t lower, upper;
  findGoodOrdinalRules(ctx, tree, node, rule.variableIndex, &lower, &upper);
  if (upper < lower) return;
  int32_t current = rule.splitIndex();
  int32_t leftChild = tree.at(node).leftChild;

  const CGMTreePrior& prior(ctx.treePrior);
  double belowX = prior.treeLogProbability(tree, ctx.data, leftChild) +
                  prior.treeLogProbability(tree, ctx.data, leftChild + 1);
  BranchScore xScore = logLikelihoodForBranch(ctx, leaf, tree, node, y, sigma);

  static thread_local Tree::SubtreeSnapshot snapshot;
  const State& s = state();
  char b[32];
  int32_t last = 0;
  for (int32_t step : {-8, -4, -2, -1, 1, 2, 4, 8}) {
    int32_t target = std::clamp(current + step, lower, upper);
    if (target - current == 0 || target - current == last) continue;
    last = target - current;

    tree.snapshotSubtree(node, snapshot);
    tree.at(node).rule.setSplitIndex(target);
    tree.refreshSubtree(ctx.data, node, y, ctx.weights);
    double belowY = prior.treeLogProbability(tree, ctx.data, leftChild) +
                    prior.treeLogProbability(tree, ctx.data, leftChild + 1);
    BranchScore yScore =
      logLikelihoodForBranch(ctx, leaf, tree, node, y, sigma);
    double xLogL, yLogL;
    resolveVetoRank(xScore, yScore, &xLogL, &yLogL);
    tree.restoreSubtree(snapshot);

    std::fprintf(file, "d,%ld,%d,%d,%d,%d,%s\n", s.sweep, s.forest, s.tree,
                 static_cast<int>(tree.depthOf(node)), last,
                 number((belowY - belowX) + (yLogL - xLogL), b, sizeof(b)));
  }
}

/// Leaf models the neighbourhood probes can enumerate: the cut scan's scalar
/// marginal over a (sum w, sum wz) pair, which neither a vector-parameter leaf
/// nor the scale leaf carries.
template <typename L>
concept ScannableLeafModel =
  ScalarLeafModel<L> && requires(const L leaf, double d) {
    { leaf.logIntegratedLikelihood(d, d, d, d) } -> std::same_as<double>;
  };

/// A discrete neighbourhood, summarized: how many candidates carry finite
/// weight, the entropy of the normalized weights in nats, the incumbent's
/// share and its rank by weight (1 the largest), and the largest share.
struct Neighbourhood {
  double candidates = 0.0;
  double entropy = 0.0;
  double incumbent = 0.0;
  double maximum = 0.0;
  double rank = 0.0;
};

/// Normalize log weights and summarize them. A candidate whose weight is the
/// scan's occupancy sentinel counts for nothing: the empty-leaf veto scores it
/// -inf, so it is in the neighbourhood at weight zero. Entries are restricted
/// to those tagged `tag` when that is non-negative, which is how the
/// single-variable neighbourhood reuses the joint enumeration.
inline Neighbourhood summarizeNeighbourhood(const std::vector<double>& logWeight,
                                            const std::vector<int32_t>& tags,
                                            int32_t tag,
                                            std::size_t incumbentIndex) {
  Neighbourhood out;
  double largest = -HUGE_VAL;
  for (std::size_t i = 0; i < logWeight.size(); ++i) {
    if (tag >= 0 && tags[i] != tag) continue;
    if (!std::isfinite(logWeight[i])) continue;
    out.candidates += 1.0;
    if (logWeight[i] > largest) largest = logWeight[i];
  }
  if (out.candidates == 0.0) return out;

  double total = 0.0;
  for (std::size_t i = 0; i < logWeight.size(); ++i) {
    if (tag >= 0 && tags[i] != tag) continue;
    if (!std::isfinite(logWeight[i])) continue;
    total += std::exp(logWeight[i] - largest);
  }
  double incumbentWeight = std::isfinite(logWeight[incumbentIndex])
    ? std::exp(logWeight[incumbentIndex] - largest) / total
    : 0.0;
  double entropy = 0.0;
  double rank = 1.0;
  for (std::size_t i = 0; i < logWeight.size(); ++i) {
    if (tag >= 0 && tags[i] != tag) continue;
    if (!std::isfinite(logWeight[i])) continue;
    double p = std::exp(logWeight[i] - largest) / total;
    if (p > 0.0) entropy -= p * std::log(p);
    if (p > incumbentWeight) rank += 1.0;
  }
  out.entropy = entropy;
  out.incumbent = incumbentWeight;
  out.maximum = std::exp(0.0) / total;
  out.rank = rank;
  return out;
}

inline void countShape(const Tree& tree, int32_t i, int* leaves, int* interior,
                       int* nog) {
  if (tree.at(i).isBottom()) {
    ++*leaves;
    return;
  }
  ++*interior;
  if (tree.childrenAreBottom(i)) ++*nog;
  countShape(tree, tree.at(i).leftChild, leaves, interior, nog);
  countShape(tree, tree.at(i).leftChild + 1, leaves, interior, nog);
}

/// The closed rule neighbourhood at the change proposal's target node, priced
/// but not drawn from. Enumerates (available ordinal variable, admissible cut)
/// and weights each by the cut scan's collapsed marginal times the surviving
/// prior factors; the node's rule is written and restored, and no membership,
/// leaf statistic or draw is touched.
template <MoveScorableLeafModel L, typename ResidT>
void nogProbe(const MoveContext& ctx, const L& leaf, Tree& tree, int32_t node,
              const ResidT* y, double sigma) {
  std::FILE* file = stream();
  if (file == nullptr) return;

  int leaves = 0, interior = 0, nog = 0;
  countShape(tree, 0, &leaves, &interior, &nog);
  bool isNog = tree.childrenAreBottom(node);
  bool scanned = false;
  Neighbourhood joint, cut;

  if constexpr (ScannableLeafModel<L>) {
    const ColumnStore& data(ctx.data);
    static thread_local std::vector<std::uint8_t> available;
    static thread_local std::vector<ConstantLeafScanBin> bins;
    static thread_local std::vector<double> scanScores;
    static thread_local std::vector<double> logWeight;
    static thread_local std::vector<int32_t> tags;

    available.resize(data.numPredictors);
    std::size_t numAvailable =
      tree.collectAvailableVariables(data, node, available.data());
    bool allOrdinal = isNog && numAvailable > 0;
    for (std::size_t j = 0; allOrdinal && j < data.numPredictors; ++j)
      if (available[j] != 0 && data.splitsBySubset(j)) allOrdinal = false;

    if (allOrdinal) {
      const Node& target(tree.at(node));
      const index_t* members = tree.indices + target.begin;
      std::size_t numMembers = target.numObservations();
      int32_t leftChild = target.leftChild;
      const Rule incumbent = target.rule;
      std::size_t incumbentIndex = 0;
      bool foundIncumbent = false;
      logWeight.clear();
      tags.clear();

      for (std::size_t j = 0; j < data.numPredictors; ++j) {
        if (available[j] == 0) continue;
        int32_t low, high;
        tree.splitInterval(data, node, static_cast<int32_t>(j), &low, &high);
        if (high < low) continue;
        std::size_t numCuts = static_cast<std::size_t>(data.numCuts[j]);
        scanScores.assign(2 * numCuts, 0.0);
        std::size_t written =
          scanOrdinalCuts(data, j, members, numMembers, y, ctx.weights, leaf,
                          ctx.k, sigma * sigma, bins, scanScores.data());
        // the doubled layout means the node routes missing rows, so the rule's
        // missing direction is part of the candidate and the rule prior widens
        // by the same factor two the candidate count does
        bool doubled = numCuts > 0 && written == 2 * numCuts;
        double logRulePrior = -std::log(static_cast<double>(high - low + 1)) -
                              (doubled ? std::log(2.0) : 0.0);
        int directions = doubled ? 2 : 1;
        for (int32_t c = low; c <= high; ++c) {
          for (int direction = 0; direction < directions; ++direction) {
            Rule& rule(tree.at(node).rule);
            rule.variableIndex = static_cast<int32_t>(j);
            rule.setSplitIndex(c);
            if (doubled) rule.setMissingGoesRight(direction == 1);
            // both children are leaves, so the below-node prior is exactly
            // changeMove's two log(1 - growth) terms
            double below =
              std::log(1.0 -
                       ctx.treePrior.growthProbability(tree, data, leftChild)) +
              std::log(1.0 - ctx.treePrior.growthProbability(tree, data,
                                                             leftChild + 1));
            if (!foundIncumbent &&
                static_cast<int32_t>(j) == incumbent.variableIndex &&
                c == incumbent.splitIndex() &&
                (!doubled || (direction == 1) == incumbent.missingGoesRight())) {
              incumbentIndex = logWeight.size();
              foundIncumbent = true;
            }
            std::size_t entry =
              doubled ? 2 * static_cast<std::size_t>(c) +
                          static_cast<std::size_t>(direction)
                      : static_cast<std::size_t>(c);
            logWeight.push_back(scanScores[entry] + logRulePrior + below);
            tags.push_back(static_cast<int32_t>(j));
          }
        }
      }
      tree.at(node).rule = incumbent;

      if (foundIncumbent) {
        joint = summarizeNeighbourhood(logWeight, tags, -1, incumbentIndex);
        cut = summarizeNeighbourhood(logWeight, tags, incumbent.variableIndex,
                                     incumbentIndex);
        scanned = true;
      }
    }
  }

  const State& s = state();
  double na = std::nan("");
  char b[10][32];
  std::fprintf(file,
               "g,%ld,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
               s.sweep, s.forest, s.tree, node,
               static_cast<int>(tree.depthOf(node)), isNog ? 1 : 0, interior,
               nog, scanned ? 1 : 0,
               number(scanned ? joint.candidates : na, b[0], sizeof(b[0])),
               number(scanned ? joint.entropy : na, b[1], sizeof(b[1])),
               number(scanned ? joint.incumbent : na, b[2], sizeof(b[2])),
               number(scanned ? joint.maximum : na, b[3], sizeof(b[3])),
               number(scanned ? joint.rank : na, b[4], sizeof(b[4])),
               number(scanned ? cut.candidates : na, b[5], sizeof(b[5])),
               number(scanned ? cut.entropy : na, b[6], sizeof(b[6])),
               number(scanned ? cut.incumbent : na, b[7], sizeof(b[7])),
               number(scanned ? cut.maximum : na, b[8], sizeof(b[8])),
               number(scanned ? cut.rank : na, b[9], sizeof(b[9])));
}

/// Informed death, priced but not drawn from: the nog nodes weighted by exp of
/// the merged-leaf marginal ratio, with the uniform pick located in that
/// distribution. Pure arithmetic on the cached leaf statistics - the merged
/// leaf's pair is the two children's sum, which is what makes the weight
/// vector free.
template <MoveScorableLeafModel L>
void deathProbe(const MoveContext& ctx, const L& leaf, const Tree& tree,
                int32_t picked, double sigma) {
  std::FILE* file = stream();
  if (file == nullptr) return;

  bool scored = false;
  Neighbourhood out;
  if constexpr (ScannableLeafModel<L>) {
    static thread_local std::vector<int32_t> nogNodes;
    static thread_local std::vector<double> logWeight;
    static thread_local std::vector<int32_t> tags;
    nogNodes.clear();
    tree.fillNoGrand(0, nogNodes);
    logWeight.clear();
    tags.clear();
    std::size_t incumbentIndex = 0;
    double residualVariance = sigma * sigma;
    for (int32_t v : nogNodes) {
      const Node& left(tree.at(tree.at(v).leftChild));
      const Node& right(tree.at(tree.at(v).leftChild + 1));
      double merged = leaf.logIntegratedLikelihood(
        ctx.k, residualVariance, left.sumWeights + right.sumWeights,
        left.sumWeightedResponse + right.sumWeightedResponse);
      double split =
        leaf.logIntegratedLikelihood(ctx.k, residualVariance, left.sumWeights,
                                     left.sumWeightedResponse) +
        leaf.logIntegratedLikelihood(ctx.k, residualVariance, right.sumWeights,
                                     right.sumWeightedResponse);
      if (v == picked) incumbentIndex = logWeight.size();
      logWeight.push_back(merged - split);
      tags.push_back(-1);
    }
    if (!logWeight.empty()) {
      out = summarizeNeighbourhood(logWeight, tags, -1, incumbentIndex);
      scored = true;
    }
  }

  const State& s = state();
  double na = std::nan("");
  char b[5][32];
  std::fprintf(file, "x,%ld,%d,%d,%s,%s,%s,%s,%s\n", s.sweep, s.forest, s.tree,
               number(scored ? out.candidates : na, b[0], sizeof(b[0])),
               number(scored ? out.entropy : na, b[1], sizeof(b[1])),
               number(scored ? out.incumbent : na, b[2], sizeof(b[2])),
               number(scored ? out.rank : na, b[3], sizeof(b[3])),
               number(scored ? out.maximum : na, b[4], sizeof(b[4])));
}

/// The perturb proposal's node and signed displacement, so consecutive
/// same-direction accepted displacements at one node can be counted offline.
inline void perturbProbe(int32_t node, int32_t current, int32_t target,
                         bool accepted) {
  std::FILE* file = stream();
  if (file == nullptr) return;
  const State& s = state();
  std::fprintf(file, "r,%ld,%d,%d,%d,%d,%d,%d\n", s.sweep, s.forest, s.tree,
               node, current, target, accepted ? 1 : 0);
}

/// One tree's settled shape, written once per tree per sweep.
inline void treeShape(const Tree& tree) {
  std::FILE* file = stream();
  if (file == nullptr) return;
  int leaves = 0, interior = 0, nog = 0;
  countShape(tree, 0, &leaves, &interior, &nog);
  const State& s = state();
  std::fprintf(file, "t,%ld,%d,%d,%d,%d,%d\n", s.sweep, s.forest, s.tree,
               leaves, interior, nog);
}

}  // namespace census

#define BARTCORE_CENSUS_NOOP(move, tree, node) census::noop(move, tree, node)
#define BARTCORE_CENSUS_SHAPE(tree, node) census::shape(tree, node)
#define BARTCORE_CENSUS_PROPOSAL(...) census::proposal(__VA_ARGS__)
#define BARTCORE_CENSUS_CUTS(...) census::cutProbe(__VA_ARGS__)
#define BARTCORE_CENSUS_LOCATION(...) census::setLocation(__VA_ARGS__)
#define BARTCORE_CENSUS_NOG(...) census::nogProbe(__VA_ARGS__)
#define BARTCORE_CENSUS_DEATHS(...) census::deathProbe(__VA_ARGS__)
#define BARTCORE_CENSUS_PERTURB(...) census::perturbProbe(__VA_ARGS__)
#define BARTCORE_CENSUS_TREE(tree) census::treeShape(tree)
#else
#define BARTCORE_CENSUS_NOOP(move, tree, node) ((void)0)
#define BARTCORE_CENSUS_SHAPE(tree, node) ((void)0)
#define BARTCORE_CENSUS_PROPOSAL(...) ((void)0)
#define BARTCORE_CENSUS_CUTS(...) ((void)0)
#define BARTCORE_CENSUS_LOCATION(...) ((void)0)
#define BARTCORE_CENSUS_NOG(...) ((void)0)
#define BARTCORE_CENSUS_DEATHS(...) ((void)0)
#define BARTCORE_CENSUS_PERTURB(...) ((void)0)
#define BARTCORE_CENSUS_TREE(tree) ((void)0)
#endif  // BARTCORE_MOVE_CENSUS

inline double probabilityOfBirthStep(const MoveContext& ctx, const Tree& tree,
                                     bool birthableNodeExists) {
  if (!birthableNodeExists) return 0.0;
  if (tree.hasSingleNode()) return 1.0;
  return ctx.birthProbability;
}

inline bool birthableNodeExists(const MoveContext& ctx, Tree& tree) {
  std::vector<int32_t>& bottoms(tree.bottomScratch);
  bottoms.clear();
  tree.fillBottom(0, bottoms);
  for (int32_t i : bottoms)
    if (tree.hasAnyAvailableVariable(ctx.data, i)) return true;
  return false;
}

inline double probabilityOfSelectingNodeForDeath(Tree& tree,
                                                 std::vector<int32_t>& scratch) {
  scratch.clear();
  tree.fillNoGrand(0, scratch);
  if (scratch.empty()) return 0.0;
  return 1.0 / static_cast<double>(scratch.size());
}

inline double probabilityOfSelectingNodeForBirth(const MoveContext& ctx,
                                                 Tree& tree) {
  if (tree.hasSingleNode()) return 1.0;

  std::vector<int32_t>& bottoms(tree.bottomScratch);
  bottoms.clear();
  tree.fillBottom(0, bottoms);

  double totalProbability = 0.0;
  for (int32_t i : bottoms)
    totalProbability += tree.hasAnyAvailableVariable(ctx.data, i) ? 1.0 : 0.0;

  if (totalProbability <= 0.0) return 0.0;
  return 1.0 / totalProbability;
}

inline int32_t drawBirthableNode(const MoveContext& ctx, ext_rng* rng, Tree& tree,
                                 double* nodeSelectionProbability) {
  if (tree.hasSingleNode()) {
    *nodeSelectionProbability = 1.0;
    return 0;
  }

  std::vector<int32_t>& bottoms(ctx.scratch.nodeScratch);
  bottoms.clear();
  tree.fillBottom(0, bottoms);

  std::vector<double>& probabilities(ctx.scratch.probabilityScratch);
  probabilities.resize(bottoms.size());
  double totalProbability = 0.0;
  for (size_t i = 0; i < bottoms.size(); ++i) {
    probabilities[i] =
      tree.hasAnyAvailableVariable(ctx.data, bottoms[i]) ? 1.0 : 0.0;
    totalProbability += probabilities[i];
  }

  if (totalProbability <= 0.0) {
    *nodeSelectionProbability = 0.0;
    return invalidNode;
  }

  misc_scalarMultiplyVectorInPlace(probabilities.data(), probabilities.size(),
                                   1.0 / totalProbability);
  size_t index = ext_rng_drawFromDiscreteDistribution(rng, probabilities.data(),
                                                      probabilities.size());
  *nodeSelectionProbability = probabilities[index];
  return bottoms[index];
}

template <MoveScorableLeafModel L, typename ResidT = double>
double birthOrDeathMove(const MoveContext& ctx, const L& leaf, ext_rng* rng,
                        Tree& tree, const ResidT* y, double sigma,
                        bool* stepTaken, bool* stepWasBirth,
                        int32_t* changedNode = nullptr) {
  // A root-only tree whose lone leaf admits no split variable can neither birth
  // (no rule to draw) nor die (no children); its move is a no-op this sweep.
  // The single-node branch below would otherwise force a birth and draw a rule
  // for invalidVariable. A movable tree never reaches here, so RNG is untouched.
  if (tree.hasSingleNode() && !birthableNodeExists(ctx, tree)) {
    *stepTaken = false;
    *stepWasBirth = false;
    BARTCORE_CENSUS_NOOP("death", tree, invalidNode);
    return 0.0;
  }

  double ratio;

  double transitionProbabilityOfSelectingNodeForBirth;
  int32_t nodeToChange =
    drawBirthableNode(ctx, rng, tree, &transitionProbabilityOfSelectingNodeForBirth);

  double transitionProbabilityOfBirthStep =
    probabilityOfBirthStep(ctx, tree, nodeToChange != invalidNode);

  if (ext_rng_simulateBernoulli(rng, transitionProbabilityOfBirthStep) == 1) {
    *stepWasBirth = true;
    BARTCORE_CENSUS_SHAPE(tree, nodeToChange);

    double parentPriorGrowthProbability =
      ctx.treePrior.growthProbability(tree, ctx.data, nodeToChange);
    double oldPriorProbability = 1.0 - parentPriorGrowthProbability;
    BranchScore oldScore =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

    // The proposal rule is drawn from the prior, so its density cancels with
    // the prior's in the acceptance ratio.
    Node oldNode = tree.at(nodeToChange);
    size_t maskPoolMark = tree.maskPoolMark();
    Rule newRule = ctx.treePrior.drawRuleAndVariable(tree, ctx.data, rng, nodeToChange);
    tree.birth(ctx.data, nodeToChange, newRule, y, ctx.weights);

    double leftPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild);
    double rightPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild + 1);
    double newPriorProbability = parentPriorGrowthProbability *
                                 (1.0 - leftPriorGrowthProbability) *
                                 (1.0 - rightPriorGrowthProbability);

    BranchScore newScore =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

    double transitionProbabilityOfDeathStep =
      1.0 - probabilityOfBirthStep(ctx, tree, birthableNodeExists(ctx, tree));
    double transitionProbabilityOfSelectingNodeForDeath =
      probabilityOfSelectingNodeForDeath(tree, ctx.scratch.nodeScratch);

    double priorRatio = newPriorProbability / oldPriorProbability;
    double transitionRatio =
      (transitionProbabilityOfDeathStep * transitionProbabilityOfSelectingNodeForDeath) /
      (transitionProbabilityOfBirthStep * transitionProbabilityOfSelectingNodeForBirth);
    double oldLogLikelihood, newLogLikelihood;
    resolveVetoRank(oldScore, newScore, &oldLogLikelihood, &newLogLikelihood);
    double likelihoodRatio = std::exp(newLogLikelihood - oldLogLikelihood);

    ratio = priorRatio * likelihoodRatio * transitionRatio;

    if (ext_rng_simulateContinuousUniform(rng) < ratio) {
      *stepTaken = true;
      if (changedNode != nullptr) *changedNode = nodeToChange;
    } else {
      // Reference behavior: the index segment stays permuted; only structure
      // and cached leaf stats are restored. A rejected pooled draw is the
      // last pool allocation, so the mark reclaims it.
      tree.undoBirth(nodeToChange);
      tree.truncateMaskPool(maskPoolMark);
      tree.at(nodeToChange).sumWeights = oldNode.sumWeights;
      tree.at(nodeToChange).sumWeightedResponse = oldNode.sumWeightedResponse;
      *stepTaken = false;
    }
    BARTCORE_CENSUS_PROPOSAL("birth", false, *stepTaken,
                             newLogLikelihood - oldLogLikelihood,
                             std::log(priorRatio), std::log(transitionRatio));
  } else {
    *stepWasBirth = false;

    double transitionProbabilityOfDeathStep = 1.0 - transitionProbabilityOfBirthStep;

    double transitionProbabilityOfSelectingNodeForDeath;
    std::vector<int32_t>& noGrand(ctx.scratch.nodeScratch);
    noGrand.clear();
    tree.fillNoGrand(0, noGrand);
    size_t index =
      ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, noGrand.size());
    transitionProbabilityOfSelectingNodeForDeath =
      1.0 / static_cast<double>(noGrand.size());
    nodeToChange = noGrand[index];
    BARTCORE_CENSUS_SHAPE(tree, nodeToChange);
    BARTCORE_CENSUS_DEATHS(ctx, leaf, tree, nodeToChange, sigma);

    double parentPriorGrowthProbability =
      ctx.treePrior.growthProbability(tree, ctx.data, nodeToChange);
    double leftPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild);
    double rightPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild + 1);
    BranchScore oldScore =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

    Node oldNode = tree.at(nodeToChange);
    tree.orphanChildren(nodeToChange);

    BranchScore newScore =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);
    double transitionProbabilityOfBirthStepReverse =
      probabilityOfBirthStep(ctx, tree, true);
    double reverseTransitionProbabilityOfSelectingNodeForBirth =
      probabilityOfSelectingNodeForBirth(ctx, tree);

    double oldPriorProbability = parentPriorGrowthProbability *
                                 (1.0 - leftPriorGrowthProbability) *
                                 (1.0 - rightPriorGrowthProbability);
    double newPriorProbability = 1.0 - parentPriorGrowthProbability;

    double priorRatio = newPriorProbability / oldPriorProbability;
    double transitionRatio =
      (transitionProbabilityOfBirthStepReverse * reverseTransitionProbabilityOfSelectingNodeForBirth) /
      (transitionProbabilityOfDeathStep * transitionProbabilityOfSelectingNodeForDeath);
    double oldLogLikelihood, newLogLikelihood;
    resolveVetoRank(oldScore, newScore, &oldLogLikelihood, &newLogLikelihood);
    double likelihoodRatio = std::exp(newLogLikelihood - oldLogLikelihood);

    ratio = priorRatio * likelihoodRatio * transitionRatio;

    if (ext_rng_simulateContinuousUniform(rng) < ratio) {
      tree.releasePair(oldNode.leftChild);
      *stepTaken = true;
      if (changedNode != nullptr) *changedNode = nodeToChange;
    } else {
      tree.at(nodeToChange) = oldNode;  // reattaches children unchanged
      *stepTaken = false;
    }
    BARTCORE_CENSUS_PROPOSAL("death", false, *stepTaken,
                             newLogLikelihood - oldLogLikelihood,
                             std::log(priorRatio), std::log(transitionRatio));
  }

  return ratio < 1.0 ? ratio : 1.0;
}

/// Ordinal cut range at nodeIndex that keeps every descendant split on the
/// same variable logically satisfiable (findGoodOrdinalRules).
inline void findGoodOrdinalRules(const MoveContext& ctx, const Tree& tree,
                                 int32_t nodeIndex, int32_t variableIndex,
                                 int32_t* lower, int32_t* upper) {
  int32_t leftIndex, rightIndex;
  tree.splitInterval(ctx.data, nodeIndex, variableIndex, &leftIndex, &rightIndex);

  // min/max split index used for the variable in each child's subtree
  struct Walker {
    static void minMax(const Tree& tree, int32_t i, int32_t variableIndex,
                       int32_t* min, int32_t* max) {
      const Node& node(tree.at(i));
      if (node.isBottom()) return;
      if (node.rule.variableIndex == variableIndex) {
        if (node.rule.splitIndex() < *min) *min = node.rule.splitIndex();
        if (node.rule.splitIndex() > *max) *max = node.rule.splitIndex();
      }
      minMax(tree, node.leftChild, variableIndex, min, max);
      minMax(tree, node.leftChild + 1, variableIndex, min, max);
    }
  };

  int32_t leftMin = rightIndex + 1, leftMaxOut = leftIndex - 1;
  int32_t rightMinOut = rightIndex + 1, rightMax = leftIndex - 1;
  Walker::minMax(tree, tree.at(nodeIndex).leftChild, variableIndex, &leftMin,
                 &leftMaxOut);
  Walker::minMax(tree, tree.at(nodeIndex).leftChild + 1, variableIndex,
                 &rightMinOut, &rightMax);
  int32_t leftMax = leftMaxOut;
  int32_t rightMin = rightMinOut;

  *lower = std::max(leftIndex, leftMax + 1);
  *upper = std::min(rightIndex, rightMin - 1);
}

inline bool categoricalSubtreeIsValid(const Tree& tree, int32_t nodeIndex,
                                      int32_t variableIndex,
                                      std::uint64_t reachable);
inline bool categoricalSubtreeIsValidWide(const Tree& tree, int32_t nodeIndex,
                                          int32_t variableIndex,
                                          const std::uint64_t* reachable,
                                          size_t numWords,
                                          std::uint64_t* arena, size_t depth);

/// Draw a categorical assignment straight from the node prior (the
/// propose-from-prior mechanism): a single unrestricted gauge-pattern draw
/// over the
/// reachable set, with no descendant-validity rejection loop. Returns true and
/// fills newRule when the draw leaves every same-variable descendant
/// satisfiable; returns false (pi(T') = 0, an automatic no-op) otherwise. The
/// prior density 1/(2^R - 2) cancels the node's rule prior exactly, so the
/// acceptance keeps only the subtree-below and likelihood ratios.
inline bool drawCategoricalRuleFromPrior(const MoveContext& ctx, ext_rng* rng,
                                         Tree& tree, int32_t nodeToChange,
                                         int32_t newVariableIndex, Rule& newRule,
                                         size_t maskPoolMark) {
  int32_t leftChild = tree.at(nodeToChange).leftChild;
  return tree.withReachableMask(
    ctx.data, nodeToChange, newVariableIndex, ctx.scratch.reachableWords,
    nullptr,
    [&](const std::uint64_t* reachable, size_t numWords) {
      size_t numReachable = maskPopcount(reachable, numWords);
      ctx.scratch.patternWords.resize(numWords);
      ctx.scratch.maskArena.resize((tree.nodes.size() + 1) * numWords);
      size_t offset = tree.allocateMask(numWords);
      CGMTreePrior::drawCategoryPatternWide(
        rng, numReachable, ctx.scratch.patternWords.data(), numWords);
      std::uint64_t* directions = tree.mutableMaskWordsFor(offset);
      CGMTreePrior::categoryDirectionsForPatternWide(
        reachable, ctx.scratch.patternWords.data(), directions, numWords);
      std::uint64_t* leftReachable = ctx.scratch.maskArena.data();
      maskAndNot(reachable, directions, leftReachable, numWords);
      if (!categoricalSubtreeIsValidWide(tree, leftChild, newVariableIndex,
                                         leftReachable, numWords,
                                         ctx.scratch.maskArena.data(), 1) ||
          !categoricalSubtreeIsValidWide(tree, leftChild + 1, newVariableIndex,
                                         directions, numWords,
                                         ctx.scratch.maskArena.data(), 1)) {
        tree.truncateMaskPool(maskPoolMark);
        return false;
      }
      newRule.setMaskOffset(offset);
      return true;
    },
    [&](std::uint64_t reachable) {
      int numReachable = std::popcount(reachable);
      std::uint64_t pattern =
        CGMTreePrior::drawCategoryPattern(rng, numReachable);
      std::uint64_t directions =
        CGMTreePrior::categoryDirectionsForPattern(reachable, pattern);
      if (!categoricalSubtreeIsValid(tree, leftChild, newVariableIndex,
                                     reachable & ~directions) ||
          !categoricalSubtreeIsValid(tree, leftChild + 1, newVariableIndex,
                                     directions))
        return false;  // narrow columns allocate no pool words to reclaim
      newRule.setCategoryDirections(directions);
      return true;
    });
}

/// Change-move proposal kernel: redraw the split variable and rule at an
/// internal node, keeping the subtree below in place. The acceptance satisfies
/// detailed balance,
///   alpha = exp( B(T') - B(T) + yLogL - xLogL + logProposalCorrection ),
/// where B is the tree-prior log-probability of the subtree STRICTLY BELOW the
/// changed node; every prior factor at or above the node cancels between T and
/// T', as in birth/death. The correction is the surviving proposal-density
/// ratio q(T|T')/q(T'|T) and composes PER SIDE, because the forward density
/// uses the new variable v''s mechanism and the reverse the old variable v's:
///   logProposalCorrection =
///       (v' ordinal ? log|Valid_T(v')| - log|SI(v')| : 0)
///     + (v  ordinal ? log|SI(v)| - log|Valid_T'(v)| : 0).
/// An ordinal side draws uniformly over the descendant-valid good set while
/// the node's rule prior normalizes over the ancestor-only interval, leaving
/// the counted ratio (|SI| the interval size, |Valid| the good-set count; the
/// variable prior and missing coin cancel within the side). A categorical side
/// draws straight from the node prior (drawCategoricalRuleFromPrior), whose
/// density cancels its side's prior factor exactly and contributes nothing; a
/// forward draw that strands a descendant split is an automatic no-op
/// (pi(T') = 0). findGoodOrdinalRules and splitInterval both ignore the node's
/// OWN rule, so the reverse ordinal count, re-enumerated on the current tree,
/// equals its value on the changed tree and always contains the old rule
/// (>= 1, never zero); a same-variable redraw gives correction 1. Omitting the
/// correction is the CGM-lineage defect this repairs.
template <MoveScorableLeafModel L, typename ResidT = double>
double changeMove(const MoveContext& ctx, const L& leaf, ext_rng* rng, Tree& tree,
                  const ResidT* y, double sigma, bool* stepTaken,
                  int32_t* changedNode = nullptr) {
  *stepTaken = false;

  std::vector<int32_t>& notBottom(ctx.scratch.nodeScratch);
  notBottom.clear();
  tree.fillNotBottom(0, notBottom);
  if (notBottom.empty()) {
    BARTCORE_CENSUS_NOOP("change", tree, invalidNode);
    return -1.0;
  }

  size_t nodeNumber =
    ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, notBottom.size());
  int32_t nodeToChange = notBottom[nodeNumber];
  BARTCORE_CENSUS_CUTS(ctx, leaf, tree, nodeToChange, y, sigma);
  BARTCORE_CENSUS_NOG(ctx, leaf, tree, nodeToChange, y, sigma);
  BARTCORE_CENSUS_SHAPE(tree, nodeToChange);

  int32_t newVariableIndex =
    ctx.treePrior.drawSplitVariable(tree, ctx.data, rng, nodeToChange);

  Rule newRule;
  newRule.variableIndex = newVariableIndex;
  // covers every exit: a pooled proposal's words are the last allocation,
  // so aborts and MH rejections truncate back; acceptance keeps them (the
  // replaced rule's words become garbage the chain compacts later)
  size_t maskPoolMark = tree.maskPoolMark();

  const bool newIsCategorical =
    ctx.data.splitsBySubset(static_cast<size_t>(newVariableIndex));
  int32_t oldVariableIndex = tree.at(nodeToChange).rule.variableIndex;
  const bool oldIsCategorical =
    ctx.data.splitsBySubset(static_cast<size_t>(oldVariableIndex));
  int32_t forwardValid = 0, forwardInterval = 0;  // new-side ordinal counts
  int32_t reverseValid = 0, reverseInterval = 0;  // old-side ordinal counts

  if (newIsCategorical) {
    // counting descendant-valid gauge patterns is exponential for wide masks,
    // so a categorical proposal draws from the prior; the density cancels the
    // node's rule prior and its side contributes no correction
    if (!drawCategoricalRuleFromPrior(ctx, rng, tree, nodeToChange,
                                      newVariableIndex, newRule,
                                      maskPoolMark)) {
      BARTCORE_CENSUS_NOOP("change", tree, nodeToChange);
      return -1.0;  // pi(T') = 0: an unsatisfiable prior draw is a no-op
    }
  } else {
    int32_t left, right;
    tree.splitInterval(ctx.data, nodeToChange, newVariableIndex, &left, &right);
    int32_t lower, upper;
    findGoodOrdinalRules(ctx, tree, nodeToChange, newVariableIndex, &lower, &upper);
    if (upper - lower + 1 <= 0) {
      BARTCORE_CENSUS_NOOP("change", tree, nodeToChange);
      return -1.0;
    }

    newRule.setSplitIndex(static_cast<int32_t>(
      ext_rng_simulateIntegerUniformInRange(rng, lower, upper + 1)));
    // like the birth draw: the missing direction is a fresh symmetric coin
    // whenever the column can route a missing value
    if (ctx.data.hasMissing[static_cast<size_t>(newVariableIndex)])
      newRule.setMissingGoesRight(ext_rng_simulateBernoulli(rng, 0.5) == 1);

    forwardValid = upper - lower + 1;
    forwardInterval = right - left + 1;
  }

  if (!oldIsCategorical) {
    // the reverse count re-enumerates the OLD variable on the current tree;
    // splitInterval and findGoodOrdinalRules both ignore the node's own rule,
    // so it always contains the old rule and never vanishes. Never run these
    // ordinal counters on a categorical column - a categorical reverse side
    // is a prior draw whose density cancels, contributing nothing.
    int32_t oldLeft, oldRight, oldLower, oldUpper;
    tree.splitInterval(ctx.data, nodeToChange, oldVariableIndex, &oldLeft,
                       &oldRight);
    findGoodOrdinalRules(ctx, tree, nodeToChange, oldVariableIndex, &oldLower,
                         &oldUpper);
    reverseValid = oldUpper - oldLower + 1;
    reverseInterval = oldRight - oldLeft + 1;
  }

  // interaction constraint: the variable was drawn feasibly against
  // nodeToChange's ANCESTORS, but the unchanged subtree below may now strand a
  // descendant (a co-occurrence or order break the drawn variable introduces).
  // treeLogProbability delegates to the availability primitives, which cannot
  // self-detect a node whose OWN variable is barred, so score the -1.0 no-op
  // (pi(T') = 0) directly, exactly as the unsatisfiable categorical draw does.
  if (tree.hasInteractionConstraint()) {
    Rule savedRule = tree.at(nodeToChange).rule;
    tree.at(nodeToChange).rule = newRule;
    bool valid = tree.interactionSubtreeIsValid(nodeToChange);
    tree.at(nodeToChange).rule = savedRule;
    if (!valid) {
      tree.truncateMaskPool(maskPoolMark);
      BARTCORE_CENSUS_NOOP("change", tree, nodeToChange);
      return -1.0;
    }
  }

  double logProposalCorrection = 0.0;
  if (!newIsCategorical && !oldIsCategorical) {
    logProposalCorrection =
      std::log(static_cast<double>(reverseInterval)) -
      std::log(static_cast<double>(forwardInterval)) +
      std::log(static_cast<double>(forwardValid)) -
      std::log(static_cast<double>(reverseValid));
  } else if (!newIsCategorical) {
    logProposalCorrection =
      std::log(static_cast<double>(forwardValid)) -
      std::log(static_cast<double>(forwardInterval));
  } else if (!oldIsCategorical) {
    logProposalCorrection =
      std::log(static_cast<double>(reverseInterval)) -
      std::log(static_cast<double>(reverseValid));
  }

  // the node's own split-variable and rule-prior factors cancel against the
  // proposal (or are carried by logProposalCorrection), so the pi ratio
  // reduces to the subtree strictly below the changed node
  int32_t leftChild = tree.at(nodeToChange).leftChild;
  double belowX =
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild) +
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild + 1);
  BranchScore xScore =
    logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

  tree.snapshotSubtree(nodeToChange, ctx.scratch.snapshot);

  tree.at(nodeToChange).rule = newRule;
  tree.refreshSubtree(ctx.data, nodeToChange, y, ctx.weights);

  double belowY =
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild) +
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild + 1);
  BranchScore yScore =
    logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

  // the veto gates the one exp the acceptance has always taken: the resolved
  // pair reproduces its exact argument wherever the ranks differ or agree at
  // 0, and supplies the finite difference where both branches are vetoed
  double xLogL, yLogL;
  resolveVetoRank(xScore, yScore, &xLogL, &yLogL);
  double alpha =
    std::exp((belowY - belowX) + (yLogL - xLogL) + logProposalCorrection);
  alpha = alpha > 1.0 ? 1.0 : alpha;

  if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
    *stepTaken = true;
    if (changedNode != nullptr) *changedNode = nodeToChange;
  } else {
    tree.restoreSubtree(ctx.scratch.snapshot);
    tree.truncateMaskPool(maskPoolMark);
  }
  BARTCORE_CENSUS_PROPOSAL("change", false, *stepTaken, yLogL - xLogL,
                           belowY - belowX, logProposalCorrection);
  return alpha;
}

/// ordinalRuleIsValid: every descendant split on variableIndex must fall
/// inside the interval implied by its ancestors.
inline bool ordinalRuleIsValid(const Tree& tree, int32_t nodeIndex,
                               int32_t variableIndex, int32_t leftIndex,
                               int32_t rightIndex) {
  const Node& node(tree.at(nodeIndex));
  if (node.isBottom()) return true;

  if (node.rule.variableIndex == variableIndex) {
    int32_t splitIndex = node.rule.splitIndex();
    if (splitIndex < leftIndex || splitIndex > rightIndex) return false;
    return ordinalRuleIsValid(tree, node.leftChild, variableIndex, leftIndex,
                              splitIndex - 1) &&
           ordinalRuleIsValid(tree, node.leftChild + 1, variableIndex,
                              splitIndex + 1, rightIndex);
  }

  return ordinalRuleIsValid(tree, node.leftChild, variableIndex, leftIndex,
                            rightIndex) &&
         ordinalRuleIsValid(tree, node.leftChild + 1, variableIndex, leftIndex,
                            rightIndex);
}

/// categoricalSubtreeIsValid: every split on variableIndex in the subtree
/// must stay in the canonical gauge (its directions confined to the
/// categories reaching it) and keep at least one reachable category on each
/// side; reachable is the mask entering the subtree.
inline bool categoricalSubtreeIsValid(const Tree& tree, int32_t nodeIndex,
                                      int32_t variableIndex,
                                      std::uint64_t reachable) {
  const Node& node(tree.at(nodeIndex));
  if (node.isBottom()) return true;

  if (node.rule.variableIndex == variableIndex) {
    std::uint64_t directions = node.rule.categoryDirections();
    if ((directions & ~reachable) != 0 || directions == 0 ||
        directions == reachable)
      return false;
    return categoricalSubtreeIsValid(tree, node.leftChild, variableIndex,
                                     reachable & ~directions) &&
           categoricalSubtreeIsValid(tree, node.leftChild + 1, variableIndex,
                                     directions);
  }

  return categoricalSubtreeIsValid(tree, node.leftChild, variableIndex,
                                   reachable) &&
         categoricalSubtreeIsValid(tree, node.leftChild + 1, variableIndex,
                                   reachable);
}

/// The pooled-column analogue: reachable spans numWords words, the left
/// branch's filtered set lives in the arena's slot at depth (deeper matching
/// rules use higher slots, so a frame's mask survives its subtree), and the
/// right branch reads the rule's own immutable pool words.
inline bool categoricalSubtreeIsValidWide(const Tree& tree, int32_t nodeIndex,
                                          int32_t variableIndex,
                                          const std::uint64_t* reachable,
                                          size_t numWords,
                                          std::uint64_t* arena, size_t depth) {
  const Node& node(tree.at(nodeIndex));
  if (node.isBottom()) return true;

  if (node.rule.variableIndex == variableIndex) {
    const std::uint64_t* directions = tree.maskWordsFor(node.rule);
    if (!maskIsSubsetOf(directions, reachable, numWords) ||
        maskIsZero(directions, numWords) ||
        maskEquals(directions, reachable, numWords))
      return false;
    std::uint64_t* leftReachable = arena + depth * numWords;
    maskAndNot(reachable, directions, leftReachable, numWords);
    return categoricalSubtreeIsValidWide(tree, node.leftChild, variableIndex,
                                         leftReachable, numWords, arena,
                                         depth + 1) &&
           categoricalSubtreeIsValidWide(tree, node.leftChild + 1,
                                         variableIndex, directions, numWords,
                                         arena, depth + 1);
  }

  return categoricalSubtreeIsValidWide(tree, node.leftChild, variableIndex,
                                       reachable, numWords, arena, depth) &&
         categoricalSubtreeIsValidWide(tree, node.leftChild + 1, variableIndex,
                                       reachable, numWords, arena, depth);
}

inline bool ruleIsValid(const MoveContext& ctx, const Tree& tree, int32_t nodeIndex,
                        int32_t variableIndex) {
  if (ctx.data.splitsBySubset(static_cast<size_t>(variableIndex))) {
    return tree.withReachableMask(
      ctx.data, nodeIndex, variableIndex, ctx.scratch.reachableWords, nullptr,
      [&](const std::uint64_t* reachable, size_t numWords) {
        ctx.scratch.maskArena.resize((tree.nodes.size() + 1) * numWords);
        return categoricalSubtreeIsValidWide(tree, nodeIndex, variableIndex,
                                             reachable, numWords,
                                             ctx.scratch.maskArena.data(), 0);
      },
      [&](std::uint64_t reachable) {
        return categoricalSubtreeIsValid(tree, nodeIndex, variableIndex,
                                         reachable);
      });
  }

  int32_t leftIndex, rightIndex;
  tree.splitInterval(ctx.data, nodeIndex, variableIndex, &leftIndex, &rightIndex);
  return ordinalRuleIsValid(tree, nodeIndex, variableIndex, leftIndex, rightIndex);
}

template <MoveScorableLeafModel L, typename ResidT = double>
double swapMove(const MoveContext& ctx, const L& leaf, ext_rng* rng, Tree& tree,
                const ResidT* y, double sigma, bool* stepTaken,
                int32_t* changedNode = nullptr) {
  *stepTaken = false;

  std::vector<int32_t>& swappable(ctx.scratch.nodeScratch);
  swappable.clear();
  tree.fillSwappable(0, swappable);
  if (swappable.empty()) {
    BARTCORE_CENSUS_NOOP("swap", tree, invalidNode);
    return -1.0;
  }

  size_t nodeNumber =
    ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, swappable.size());
  int32_t parent = swappable[nodeNumber];
  BARTCORE_CENSUS_SHAPE(tree, parent);
  int32_t leftChild = tree.at(parent).leftChild;
  int32_t rightChild = leftChild + 1;

  bool leftHasRule = !tree.at(leftChild).isBottom();
  bool rightHasRule = !tree.at(rightChild).isBottom();
  bool childrenHaveSameRule = leftHasRule && rightHasRule &&
    tree.rulesAreEqual(ctx.data, tree.at(leftChild).rule,
                       tree.at(rightChild).rule);

  double alpha;

  // The swap gives the parent a child's rule and each swapped child the
  // parent's. When the two children share a rule both are swapped; otherwise
  // one is picked (a fair coin when both carry a rule, else the only one that
  // does). Either way every swapped child's original rule is childRule, so the
  // two cases differ only in which children move - captured in swapChildren.
  Rule parentRule = tree.at(parent).rule;
  Rule childRule;
  int32_t swapChildren[2];
  int numSwapChildren;
  if (!childrenHaveSameRule) {
    int32_t child;
    if (leftHasRule && rightHasRule) {
      child = ext_rng_simulateBernoulli(rng, 0.5) == 1 ? leftChild : rightChild;
    } else {
      child = leftHasRule ? leftChild : rightChild;
    }
    childRule = tree.at(child).rule;
    swapChildren[0] = child;
    numSwapChildren = 1;
  } else {
    childRule = tree.at(leftChild).rule;
    swapChildren[0] = leftChild;
    swapChildren[1] = rightChild;
    numSwapChildren = 2;
  }

  auto applySwap = [&]() {
    tree.at(parent).rule = childRule;
    for (int i = 0; i < numSwapChildren; ++i)
      tree.at(swapChildren[i]).rule = parentRule;
  };
  auto undoSwap = [&]() {
    tree.at(parent).rule = parentRule;
    for (int i = 0; i < numSwapChildren; ++i)
      tree.at(swapChildren[i]).rule = childRule;
  };

  // test the swap for logical consistency before scoring it
  applySwap();
  bool swapIsSensible = ruleIsValid(ctx, tree, parent, childRule.variableIndex);
  if (childRule.variableIndex != parentRule.variableIndex && swapIsSensible)
    swapIsSensible = ruleIsValid(ctx, tree, parent, parentRule.variableIndex);
  // interaction is a WHOLE-subtree, all-variables property the per-variable
  // ruleIsValid checks above cannot see (the swap sibling-strand break): a
  // swap that lifts x2 above x3 co-occurs a forbidden pair with neither
  // swapped variable equal to x3. Score it the -1.0 no-op (pi(T') = 0).
  if (swapIsSensible) swapIsSensible = tree.interactionSubtreeIsValid(parent);
  undoSwap();

  if (!swapIsSensible) {
    BARTCORE_CENSUS_NOOP("swap", tree, parent);
    return -1.0;
  }

  // as in changeMove, prior terms outside the swapped subtree cancel
  double xLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data, parent);
  BranchScore xScore =
    logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

  tree.snapshotSubtree(parent, ctx.scratch.snapshot);
  applySwap();
  tree.refreshSubtree(ctx.data, parent, y, ctx.weights);

  double yLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data, parent);
  BranchScore yScore =
    logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

  // as in changeMove, the veto gates the existing single exp
  double xLogL, yLogL;
  resolveVetoRank(xScore, yScore, &xLogL, &yLogL);
  alpha = std::exp(yLogPi + yLogL - xLogPi - xLogL);
  alpha = alpha > 1.0 ? 1.0 : alpha;

  if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
    *stepTaken = true;
    if (changedNode != nullptr) *changedNode = parent;
  } else {
    tree.restoreSubtree(ctx.scratch.snapshot);
  }
  BARTCORE_CENSUS_PROPOSAL("swap", false, *stepTaken, yLogL - xLogL,
                           yLogPi - xLogPi, 0.0);

  return alpha;
}

/// Window half-width for the perturb move, in GRID POSITIONS: a displacement
/// is drawn from the cuts within perturbWidth of the current one. A
/// compile-time constant rather than a knob - acceptance falls off steeply
/// with the displacement and no caller can set it from evidence - and a width
/// arm therefore needs a private build.
inline constexpr int32_t perturbWidth = 1;

/// Perturb-move proposal kernel: displace one interior node's ordinal cut by
/// at most perturbWidth grid positions, keeping its split VARIABLE and the
/// whole tree's shape. The acceptance is changeMove's with the node's own
/// prior factors cancelling exactly rather than against a proposal density:
///   alpha = exp( B(T') - B(T) + yLogL - xLogL + logProposalCorrection ),
/// where B is the tree-prior log-probability of the subtree STRICTLY BELOW the
/// node. splitVariableLogProbability reads ancestors only and the rule prior
/// normalizes over the ancestor-constrained interval, both of which an
/// unchanged variable leaves fixed, so changeMove's per-side |Valid|/|SI|
/// machinery collapses to the window ratio alone:
///   W(c) = { j in [lo, hi] : 0 < |j - c| <= w },  |W(c)| = min(hi, c + w) -
///                                                          max(lo, c - w),
///   logProposalCorrection = log|W(c)| - log|W(c')|.
/// The correction is exact and needs no re-enumeration on T': splitInterval
/// and findGoodOrdinalRules both ignore the node's OWN rule and read only
/// ancestors and descendants, neither of which a displacement touches, so
/// [lo, hi] is identical on T and T' and the reverse count is taken on the
/// unmodified tree. Equalling or crossing an ancestor's or a descendant's cut
/// on the same variable is impossible rather than handled, [lo, hi] being set
/// one index inside both. hi == lo leaves an empty window and a no-op.
///
/// Eligible nodes are the interior nodes on an ORDINAL column, and the filter
/// is not an optimization: the selected set must be a tree function a
/// displacement cannot move, or its reciprocal stops cancelling between T and
/// T'. For the same reason no width filter may skip a node - a degenerate
/// interval is a no-op INSIDE the kernel. A categorical rule has no cut to
/// displace and is never proposed, so the move is inert on an all-categorical
/// design.
///
/// Three of changeMove's checks drop with the variable held: no mask pool (an
/// ordinal rule allocates no words), no interaction walk (co-occurrence and
/// order are properties of the split VARIABLES), and no stranding check
/// ([lo, hi] strands none). The missing direction rides the displaced rule
/// unchanged - setSplitIndex preserves it - and contributes log 2 to both
/// sides of the prior ratio. [lo, hi] guarantees satisfiability but never
/// occupancy, so a displaced cut can still empty a descendant leaf; the veto
/// resolves that as it does for change.
template <MoveScorableLeafModel L, typename ResidT = double>
double perturbMove(const MoveContext& ctx, const L& leaf, ext_rng* rng,
                   Tree& tree, const ResidT* y, double sigma, bool* stepTaken,
                   int32_t* changedNode = nullptr) {
  *stepTaken = false;

  std::vector<int32_t>& eligible(ctx.scratch.nodeScratch);
  eligible.clear();
  tree.fillNotBottom(0, eligible);
  size_t numEligible = 0;
  for (int32_t i : eligible)
    if (!ctx.data.splitsBySubset(
          static_cast<size_t>(tree.at(i).rule.variableIndex)))
      eligible[numEligible++] = i;
  eligible.resize(numEligible);
  if (eligible.empty()) {
    BARTCORE_CENSUS_NOOP("perturb", tree, invalidNode);
    return -1.0;
  }

  size_t nodeNumber =
    ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, eligible.size());
  int32_t nodeToPerturb = eligible[nodeNumber];
  BARTCORE_CENSUS_SHAPE(tree, nodeToPerturb);

  int32_t variableIndex = tree.at(nodeToPerturb).rule.variableIndex;
  int32_t current = tree.at(nodeToPerturb).rule.splitIndex();
  int32_t lower, upper;
  findGoodOrdinalRules(ctx, tree, nodeToPerturb, variableIndex, &lower, &upper);

  // |W(c)|, the window less the current cut; zero at a degenerate interval
  int32_t forwardLow = std::max(lower, current - perturbWidth);
  int32_t forwardCount = std::min(upper, current + perturbWidth) - forwardLow;
  if (forwardCount <= 0) {
    BARTCORE_CENSUS_NOOP("perturb", tree, nodeToPerturb);
    return -1.0;
  }

  // one draw over W(c), the current cut skipped by shifting the upper half up
  int32_t target = forwardLow + static_cast<int32_t>(
    ext_rng_simulateIntegerUniformInRange(rng, 0, forwardCount));
  if (target >= current) ++target;

  int32_t reverseCount = std::min(upper, target + perturbWidth) -
                         std::max(lower, target - perturbWidth);
  double logProposalCorrection =
    std::log(static_cast<double>(forwardCount)) -
    std::log(static_cast<double>(reverseCount));

  // the node's own split-variable and rule-prior factors cancel exactly, so
  // the pi ratio reduces to the subtree strictly below the perturbed node
  int32_t leftChild = tree.at(nodeToPerturb).leftChild;
  double belowX =
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild) +
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild + 1);
  BranchScore xScore =
    logLikelihoodForBranch(ctx, leaf, tree, nodeToPerturb, y, sigma);

  tree.snapshotSubtree(nodeToPerturb, ctx.scratch.snapshot);

  tree.at(nodeToPerturb).rule.setSplitIndex(target);
  tree.refreshSubtree(ctx.data, nodeToPerturb, y, ctx.weights);

  double belowY =
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild) +
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild + 1);
  BranchScore yScore =
    logLikelihoodForBranch(ctx, leaf, tree, nodeToPerturb, y, sigma);

  // as in changeMove, the veto gates the single exp the acceptance takes
  double xLogL, yLogL;
  resolveVetoRank(xScore, yScore, &xLogL, &yLogL);
  double alpha =
    std::exp((belowY - belowX) + (yLogL - xLogL) + logProposalCorrection);
  alpha = alpha > 1.0 ? 1.0 : alpha;

  if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
    *stepTaken = true;
    if (changedNode != nullptr) *changedNode = nodeToPerturb;
  } else {
    tree.restoreSubtree(ctx.scratch.snapshot);
  }
  BARTCORE_CENSUS_PROPOSAL("perturb", false, *stepTaken, yLogL - xLogL,
                           belowY - belowX, logProposalCorrection);
  BARTCORE_CENSUS_PERTURB(nodeToPerturb, current, target, *stepTaken);
  return alpha;
}

enum class StepType { birth, death, swap, change, perturb };

/// True when the move mixture proposes no structure at all: every structural
/// probability is exactly zero, so the trees stand as they are and only the
/// leaf values, the residual scale and the family's latents keep moving. A
/// sweep reads this once per forest and skips metropolisJumpForTree entirely
/// where it holds - the frozen path draws no uniform for the move choice.
inline bool structureIsFrozen(double birthOrDeathProbability,
                              double swapProbability, double changeProbability,
                              double perturbProbability) {
  return birthOrDeathProbability == 0.0 && swapProbability == 0.0 &&
         changeProbability == 0.0 && perturbProbability == 0.0;
}

/// changedNode, when non-null, receives the index of the node whose subtree an
/// ACCEPTED move repartitioned (the birthed/died node, or the changed, swapped
/// or perturbed subtree root); untouched on rejection or no-op, so gate reads
/// on stepTaken.
///
/// The perturb branch tests at birthOrDeath + swap + perturb and change stays
/// the else, so at a perturb probability of exactly zero the added test IS the
/// swap test in IEEE, fails wherever that one failed, and control reaches
/// changeMove at the same stream position.
template <MoveScorableLeafModel L, typename ResidT = double>
double metropolisJumpForTree(const MoveContext& ctx, const L& leaf, ext_rng* rng,
                             Tree& tree, const ResidT* y, double sigma,
                             bool* stepTaken, StepType* stepType,
                             int32_t* changedNode = nullptr) {
  double alpha;
  double u = ext_rng_simulateContinuousUniform(rng);

  if (u < ctx.birthOrDeathProbability) {
    bool birthed;
    alpha = birthOrDeathMove(ctx, leaf, rng, tree, y, sigma, stepTaken, &birthed,
                             changedNode);
    *stepType = birthed ? StepType::birth : StepType::death;
  } else if (u < ctx.birthOrDeathProbability + ctx.swapProbability) {
    alpha = swapMove(ctx, leaf, rng, tree, y, sigma, stepTaken, changedNode);
    *stepType = StepType::swap;
  } else if (u < ctx.birthOrDeathProbability + ctx.swapProbability +
                   ctx.perturbProbability) {
    alpha = perturbMove(ctx, leaf, rng, tree, y, sigma, stepTaken, changedNode);
    *stepType = StepType::perturb;
  } else {
    alpha = changeMove(ctx, leaf, rng, tree, y, sigma, stepTaken, changedNode);
    *stepType = StepType::change;
  }

  return alpha;
}

}  // namespace bartcore

#endif  // BARTCORE_MOVES_HPP
