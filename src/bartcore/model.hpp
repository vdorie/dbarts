#ifndef BARTCORE_MODEL_HPP
#define BARTCORE_MODEL_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstring>
#include <cfloat>
#include <limits>
#include <memory>
#include <numbers>
#include <type_traits>
#include <vector>

#include <external/random.h>
#include <external/stats.h> // Rf_dnorm4, Rf_pnorm5, Rf_dt for the log-likelihood channel
#include <misc/stats.h>
#include <misc/linearAlgebra.h>

#include "data.hpp"
#include "tree.hpp"

namespace bartcore {

/// Leaf models the conjugate engine can run must be integrable: they expose
/// a closed-form log marginal over their parameters, scored through node
/// context - (tree, working response, weights, node) - which is all the
/// moves ever see. Parameter shape splits the hierarchy below that: scalar
/// models (one parameter per leaf) keep the scalar draw interface,
/// vector models write numParams() doubles per leaf and evaluate fits per
/// observation, function-valued models draw one value per member observation
/// directly into the fit vector (the fits ARE the parameters). Chain code
/// branches on the shape traits at compile time, so a scalar leaf model
/// never instantiates the vector or function paths.
template <typename L>
concept LeafModelCore =
  requires(const L leaf, const Tree& tree, const double* v, double d,
           std::int32_t node) {
    { L::hasVectorParams } -> std::convertible_to<bool>;
    { L::hasFunctionParams } -> std::convertible_to<bool>;
    { leaf.logIntegratedLikelihoodForNode(tree, v, v, d, d, node) }
      -> std::same_as<double>;
  };

template <typename L>
concept ScalarLeafModel = LeafModelCore<L> && !L::hasVectorParams &&
  !L::hasFunctionParams &&
  requires(const L leaf, ext_rng* rng, const Tree& tree, double d,
           std::int32_t node) {
    { leaf.drawFromPosteriorForNode(rng, tree, d, d, node) }
      -> std::same_as<double>;
    { leaf.drawFromPrior(rng, d) } -> std::same_as<double>;
  };

template <typename L>
concept VectorLeafModel = LeafModelCore<L> && L::hasVectorParams &&
  requires(const L leaf, ext_rng* rng, const Tree& tree, const double* v,
           double d, std::int32_t node, double* out, std::size_t i) {
    { leaf.numParams() } -> std::same_as<std::size_t>;
    { leaf.drawFromPosteriorForNode(rng, tree, v, v, d, d, node, out) }
      -> std::same_as<void>;
    { leaf.drawFromPrior(rng, d, out) } -> std::same_as<void>;
    { leaf.fitForObservation(v, i) } -> std::same_as<double>;
    { leaf.fitForTestObservation(v, i) } -> std::same_as<double>;
  };

/// Sufficient statistics from one function-valued leaf draw, feeding the
/// leaf-scale chi-squared update: the standardized sum of squares
/// (f' C^-1 f for a GP draw, the squared value for a constant-fallback draw)
/// in sumSquaredParams, and the draw's parameter count in numParams.
struct FunctionLeafDrawStats {
  double sumSquaredParams = 0.0;
  double numParams = 0.0;
};

/// Function-valued leaves: draws write one value per member observation into
/// the caller's fit vector, so parameters need no per-tree storage. Test
/// rows evaluate through a per-draw cache the chain resets with
/// beginTreeDraw before each tree's parameter sweep.
template <typename L>
concept FunctionLeafModel = LeafModelCore<L> && L::hasFunctionParams &&
  requires(const L leaf, ext_rng* rng, const Tree& tree, const double* v,
           double d, std::int32_t node, double* fits, std::size_t i) {
    { leaf.beginTreeDraw(tree) } -> std::same_as<void>;
    { leaf.drawFromPosteriorForNode(rng, tree, v, v, d, d, node, fits) }
      -> std::same_as<FunctionLeafDrawStats>;
    { leaf.drawFromPriorForNode(rng, tree, d, node, fits) }
      -> std::same_as<FunctionLeafDrawStats>;
    { leaf.fitForTestObservationForNode(tree, node, i) }
      -> std::same_as<double>;
  };

template <typename L>
concept IntegrableLeafModel =
  ScalarLeafModel<L> || VectorLeafModel<L> || FunctionLeafModel<L>;

/// A conjugate scale leaf: one positive variance factor per leaf, scored and
/// drawn through the weighted sum of squared residuals over the node. It
/// shares the marginal shape of a ScalarLeafModel (LeafModelCore, one
/// parameter) but its draw consumes (y, weights) directly - the conjugate
/// statistic is a scale suffstat the node does not cache - rather than the
/// cached location suffstat (sum w, sum wz) the ScalarLeafModel draw reads.
/// The variance forest (heteroscedastic BART, docs/design/heteroscedastic.md)
/// declares it; no location leaf does, and it is not IntegrableLeafModel, so
/// the scale path never instantiates in an existing sampler.
template <typename L>
concept ScaleLeafModel = LeafModelCore<L> && !L::hasVectorParams &&
  !L::hasFunctionParams &&
  requires(const L leaf, ext_rng* rng, const Tree& tree, const double* v,
           double d, std::int32_t node) {
    { leaf.drawFromPosteriorForNode(rng, tree, v, v, d, d, node) }
      -> std::same_as<double>;
    { leaf.drawFromPrior(rng) } -> std::same_as<double>;
  };

/// Leaf models whose birth/death/change/swap structure moves the conjugate
/// machinery can score: any integrable leaf, plus the scale leaf. The scale
/// leaf falls through the SAME per-leaf marginal sum (logLikelihoodForBranch)
/// with a scale statistic in place of a location one; only the leaf parameter
/// DRAW differs, and that is dispatched outside the move.
template <typename L>
concept MoveScorableLeafModel = IntegrableLeafModel<L> || ScaleLeafModel<L>;

/// Optional leaf-model seams the engine dispatches on at compile time; no
/// default leaf declares either, so both compile out and the rng stream is
/// unchanged. A TreeDrawLeafModel draws its whole leaf vector in one coupled
/// sweep (the monotone constrained draw) in place of the default independent
/// per-node draw. A ParamScoringLeafModel scores a structure move against the
/// current leaf-parameter vector M (the tree's persistent block) rather than
/// integrating every leaf out - the generic hook a constrained-conjugate
/// (monotone) or non-conjugate (heteroscedastic) move specializes.
template <typename L>
concept TreeDrawLeafModel =
  requires(const L leaf, ext_rng* rng, const Tree& tree,
           const std::vector<std::int32_t>& bottoms, double d, double* out) {
    { leaf.drawParametersForTree(rng, tree, bottoms, d, d, out) }
      -> std::same_as<void>;
  };

template <typename L>
concept ParamScoringLeafModel =
  requires(const L leaf, const Tree& tree, std::int32_t branchIndex,
           const double* y, const double* weights, double k, double sigma2,
           const double* leafParams) {
    { leaf.logLikelihoodForBranchWithParams(tree, branchIndex, y, weights, k,
                                            sigma2, leafParams) }
      -> std::same_as<double>;
  };

/// Constant Gaussian leaf: mu ~ N(0, (scale / k)^2), Gaussian likelihood.
struct ConstantGaussianLeaf {
  static constexpr bool hasVectorParams = false;
  static constexpr bool hasFunctionParams = false;

  double scale;  // nodeScale / sqrt(numTrees)

  // The (sum w, sum wz) sufficient statistic drives the marginal and the draw
  // directly. The raw response sum of squares a full CGM marginal carries is
  // additive over a node's members, so it cancels in every birth/death MH
  // ratio; dropping it leaves only the mean-explained term below.
  double logIntegratedLikelihood(double k, double residualVariance,
                                 double sumWeights,
                                 double sumWeightedResponse) const {
    // negated so a NaN or a negative total lands here as well; both are
    // unreachable under the bridge's nonnegative-weight validation and a
    // combiner's w m^2, and this guard is all that stands between them and a
    // division by the total below
    if (!(sumWeights > 0.0)) return 0.0;

    double priorPrecision = (k / scale) * (k / scale);
    double posteriorPrecision = sumWeights / residualVariance;
    double mean = sumWeightedResponse / sumWeights;
    double explainedSumOfSquares = sumWeightedResponse * mean;

    double result;
    result  = 0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision));
    result += 0.5 * explainedSumOfSquares / residualVariance;
    result -= 0.5 * ((priorPrecision * mean) * (posteriorPrecision * mean)) /
              (priorPrecision + posteriorPrecision);
    return result;
  }

  double drawFromPosterior(ext_rng* rng, double k, double sumWeights,
                           double sumWeightedResponse,
                           double residualVariance) const {
    double priorPrecision = (k / scale) * (k / scale);
    double posteriorPrecision = sumWeights / residualVariance;
    double posteriorMean = (sumWeightedResponse / residualVariance) /
                           (priorPrecision + posteriorPrecision);
    double posteriorSd = 1.0 / std::sqrt(priorPrecision + posteriorPrecision);
    return posteriorMean + posteriorSd * ext_rng_simulateStandardNormal(rng);
  }

  double drawFromPrior(ext_rng* rng, double k) const {
    return (scale / k) * ext_rng_simulateStandardNormal(rng);
  }

  // The node-context interface the engine drives; delegates to the scalar
  // math above against the node's cached sufficient statistic. y and weights
  // are unused now that the node caches both sums; the residual pointer is
  // templated so the fp32-residual path (const float*) resolves here without a
  // conversion (docs/design/reduced-precision-storage.md sec 3b).
  template <typename ResidT>
  double logIntegratedLikelihoodForNode(const Tree& tree, const ResidT*,
                                        const double*, double k,
                                        double residualVariance,
                                        int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    return logIntegratedLikelihood(k, residualVariance, node.sumWeights,
                                   node.sumWeightedResponse);
  }

  double drawFromPosteriorForNode(ext_rng* rng, const Tree& tree, double k,
                                  double residualVariance,
                                  int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    return drawFromPosterior(rng, k, node.sumWeights,
                             node.sumWeightedResponse, residualVariance);
  }
};

static_assert(ScalarLeafModel<ConstantGaussianLeaf>);

/// Constant variance (scale) leaf for the heteroscedastic variance forest
/// (HBART; Pratola, Chipman, George, McCulloch 2020,
/// docs/design/heteroscedastic.md). Each leaf carries one positive variance
/// factor s^2_k ~ chi^-2(nu', lambda'^2) - the scaled-inverse-chi-squared prior
/// dbarts's sigma prior uses (ChiSquaredScalePrior), here calibrated per
/// variance tree. The conjugate statistic is (n_k, sum_k w_i r_i^2) over the
/// node's index span, where r_i is the caller-supplied scaled residual
/// e_i / sqrt(s^2_{-j}(x_i)) and w_i the case weight; the leaf accumulates its
/// own per call (the LinearGaussianLeaf precedent), so no shared node stat
/// widens and the mean path stays byte-identical. `scale` is on the VARIANCE
/// scale (= ChiSquaredScalePrior::scale = Pratola et al. (2020)'s lambda'^2),
/// so at m' = 1
/// the calibration collapses to the homoscedastic (nu, scale) exactly.
struct ConstantVarianceLeaf {
  static constexpr bool hasVectorParams = false;
  static constexpr bool hasFunctionParams = false;

  double degreesOfFreedom = 3.0;  // nu'
  double scale = 1.0;             // lambda'^2, variance scale

  /// Per-tree calibration from the homoscedastic sigma prior
  /// (nu, scale) over m' variance trees, mean-matching the leaf product's
  /// prior to chi^-2(nu, scale): lambda'^2 = scale^(1/m') and
  /// nu' = 2 / [1 - (1 - 2/nu)^(1/m')]. At m' = 1 this returns (nu, scale).
  static ConstantVarianceLeaf calibrated(double df, double scale,
                                         std::size_t numTrees) {
    double m = static_cast<double>(numTrees);
    double leafScale = std::pow(scale, 1.0 / m);
    double leafDf = 2.0 / (1.0 - std::pow(1.0 - 2.0 / df, 1.0 / m));
    return ConstantVarianceLeaf{leafDf, leafScale};
  }

  /// log integral of prod N(r_i; 0, s^2 / w_i) against the chi^-2(nu, scale)
  /// prior (paper eq. 7), dropping the (2 pi)^{-n/2} (prod w_i)^{1/2} factor
  /// additive over the node's members - conserved under any repartition, the
  /// same reason the constant leaf drops its raw sum of squares. With n
  /// positive-weight rows and ssr = sum w_i r_i^2:
  ///   (nu/2) log(nu scale / 2) - lgamma(nu/2)
  ///     + lgamma((n+nu)/2) - ((n+nu)/2) log((ssr + nu scale)/2).
  /// Empty (n = 0) returns 0, so the empty-leaf veto (moves.hpp) is unchanged.
  double logIntegratedLikelihood(double n, double ssr) const {
    double priorScale = degreesOfFreedom * scale;
    return 0.5 * degreesOfFreedom * std::log(0.5 * priorScale) -
           std::lgamma(0.5 * degreesOfFreedom) +
           std::lgamma(0.5 * (n + degreesOfFreedom)) -
           0.5 * (n + degreesOfFreedom) * std::log(0.5 * (ssr + priorScale));
  }

  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* weights, double, double,
                                        int32_t nodeIndex) const {
    double n, ssr;
    accumulate(tree, y, weights, nodeIndex, n, ssr);
    return logIntegratedLikelihood(n, ssr);
  }

  /// s^2_k | r ~ chi^-2(nu + n, [nu scale + ssr] / (nu + n)) (paper eq. 6),
  /// drawn as (nu scale + ssr) / chisq(nu + n) - ChiSquaredScalePrior's update
  /// per leaf. An empty leaf falls through to the prior draw (no data term).
  double drawFromPosteriorForNode(ext_rng* rng, const Tree& tree,
                                  const double* y, const double* weights,
                                  double, double, int32_t nodeIndex) const {
    double n, ssr;
    accumulate(tree, y, weights, nodeIndex, n, ssr);
    double posteriorScale = degreesOfFreedom * scale + ssr;
    return posteriorScale /
           ext_rng_simulateChiSquared(rng, degreesOfFreedom + n);
  }

  double drawFromPrior(ext_rng* rng) const {
    return degreesOfFreedom * scale /
           ext_rng_simulateChiSquared(rng, degreesOfFreedom);
  }

private:
  /// One pass over the node's index segment: the positive-weight count and the
  /// weighted sum of squared (already scaled) residuals. Zero-weight rows drop
  /// from both, as they drop from ChiSquaredScalePrior's posterior df and SSR.
  static void accumulate(const Tree& tree, const double* y,
                         const double* weights, int32_t nodeIndex, double& n,
                         double& ssr) {
    const Node& node(tree.at(nodeIndex));
    n = 0.0;
    ssr = 0.0;
    for (std::size_t m = node.begin; m < node.end; ++m) {
      std::size_t i = tree.indices[m];
      double w = weights == nullptr ? 1.0 : weights[i];
      if (w <= 0.0) continue;
      n += 1.0;
      ssr += w * y[i] * y[i];
    }
  }
};

static_assert(ScaleLeafModel<ConstantVarianceLeaf>);
static_assert(!ScalarLeafModel<ConstantVarianceLeaf>);
static_assert(!IntegrableLeafModel<ConstantVarianceLeaf>);

// ---- Monotone (mBART) constrained constant leaf ---------------------------
//
// A constant leaf whose forest is monotone in a declared subset of predictors
// (Chipman, George, McCulloch, Shively 2022; docs/design/monotone.md). Each
// tree is constrained so neighboring leaves along a constrained axis are
// ordered; the leaf-parameter draw is a coupled truncated-normal Gibbs sweep
// and the birth/death score is the constrained (truncated) conditional
// marginal (scored under the truncation, not the naive unconstrained marginal
// nor a full integrate-out of the neighbors). It DECLARES the two optional
// leaf seams
// (TreeDrawLeafModel, ParamScoringLeafModel) that the constant leaf does not,
// so it is a separate construction-time instantiation; the unconstrained path
// is byte-identical (design section 8).

inline double gaussianPdf(double z) {
  return std::exp(-0.5 * z * z) / std::sqrt(2.0 * std::numbers::pi);
}
inline double gaussianCdf(double z) {
  return Rf_pnorm5(z, 0.0, 1.0, 1, 0);
}

/// Adaptive Simpson quadrature of a smooth integrand over [a, b]; the callers
/// standardize so the integrand is a standard-normal density times a bounded
/// factor, keeping the recursion shallow.
template <typename F>
double monotoneAdaptiveSimpson(F f, double a, double m, double b, double fa,
                               double fm, double fb, double whole, double tol,
                               int depth) {
  double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
  double flm = f(lm), frm = f(rm);
  double left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
  double right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
  double both = left + right;
  if (depth <= 0 || std::abs(both - whole) <= 15.0 * tol)
    return both + (both - whole) / 15.0;
  return monotoneAdaptiveSimpson(f, a, lm, m, fa, flm, fm, left, 0.5 * tol,
                                 depth - 1) +
         monotoneAdaptiveSimpson(f, m, rm, b, fm, frm, fb, right, 0.5 * tol,
                                 depth - 1);
}
template <typename F>
double monotoneIntegrate(F f, double a, double b, double tol = 1e-12,
                         int maxDepth = 30) {
  if (!(b > a)) return 0.0;
  // an initial uniform split guarantees a narrow peak anywhere in [a, b] is
  // sampled before the adaptive convergence test can stop early on a wide,
  // near-flat interval; each panel is then refined to tolerance.
  const int panels = 16;
  double h = (b - a) / panels;
  double total = 0.0;
  for (int panel = 0; panel < panels; ++panel) {
    double lo = a + panel * h, hi = lo + h;
    double m = 0.5 * (lo + hi);
    double fa = f(lo), fm = f(m), fb = f(hi);
    double whole = (hi - lo) / 6.0 * (fa + 4.0 * fm + fb);
    total += monotoneAdaptiveSimpson(f, lo, m, hi, fa, fm, fb, whole,
                                     tol / panels, maxDepth);
  }
  return total;
}

/// Ordinal code-box of a leaf along one predictor: codes in [lo, hi] reach it.
/// From the ancestor-constrained cut interval (Tree::splitInterval): its left
/// is the low code, its right + 1 the high code.
inline void monotoneLeafBox(const Tree& tree, const ColumnStore& data,
                            std::int32_t leaf, std::int32_t variable,
                            std::int32_t* lo, std::int32_t* hi) {
  std::int32_t left, right;
  tree.splitInterval(data, leaf, variable, &left, &right);
  *lo = left;
  *hi = right + 1;
}

/// The distinct split variables on a leaf's ancestor path (the only axes whose
/// boxes are narrower than the full range, so the only ones an overlap test
/// must check).
inline void monotonePathVariables(const Tree& tree, std::int32_t leaf,
                                  std::vector<std::int32_t>& out) {
  out.clear();
  std::int32_t current = leaf;
  while (tree.at(current).parent != invalidNode) {
    current = tree.at(current).parent;
    std::int32_t v = tree.at(current).rule.variableIndex;
    if (std::find(out.begin(), out.end(), v) == out.end()) out.push_back(v);
  }
}

/// Whether two leaves' boxes overlap on every axis in `axes` except skipAxis.
inline bool monotoneBoxesOverlap(const Tree& tree, const ColumnStore& data,
                                 std::int32_t leafA, std::int32_t leafB,
                                 const std::vector<std::int32_t>& axes,
                                 std::int32_t skipAxis) {
  for (std::int32_t v : axes) {
    if (v == skipAxis) continue;
    std::int32_t loA, hiA, loB, hiB;
    monotoneLeafBox(tree, data, leafA, v, &loA, &hiA);
    monotoneLeafBox(tree, data, leafB, v, &loB, &hiB);
    if (std::max(loA, loB) > std::min(hiA, hiB)) return false;
  }
  return true;
}

/// Scratch for the neighbor walk (per-leaf, reused across the sweep).
struct MonotoneNeighborScratch {
  std::vector<std::int32_t> branch, allBottoms, pathA, pathB, axes;
};

/// Bounds [a, b] on leaf k's value from its neighbors' frozen mu, plus whether
/// k has any neighbor along a constrained axis (c-inflation). j below k along
/// an increasing axis lower-bounds mu_k by mu_j; above upper-bounds it; a
/// decreasing axis (direction -1) flips the roles. A neighbor whose index is
/// in `skip` still marks k constrained but contributes no bound - the touched
/// leaves of a move, integrated out rather than frozen.
inline void monotoneNeighborBounds(const Tree& tree, const ColumnStore& data,
                                   const std::int8_t* directions,
                                   const std::vector<std::int32_t>& bottoms,
                                   std::int32_t k, const double* mu,
                                   const std::int32_t* skip, std::size_t numSkip,
                                   MonotoneNeighborScratch& scratch, double* aOut,
                                   double* bOut, bool* constrained) {
  double a = -HUGE_VAL, b = HUGE_VAL;
  bool hasConstrainedNeighbor = false;
  monotonePathVariables(tree, k, scratch.pathA);
  for (std::int32_t j : bottoms) {
    if (j == k) continue;
    monotonePathVariables(tree, j, scratch.pathB);
    scratch.axes = scratch.pathA;
    for (std::int32_t v : scratch.pathB)
      if (std::find(scratch.axes.begin(), scratch.axes.end(), v) ==
          scratch.axes.end())
        scratch.axes.push_back(v);
    bool jSkipped = false;
    for (std::size_t t = 0; t < numSkip; ++t)
      if (skip[t] == j) jSkipped = true;
    for (std::int32_t i : scratch.axes) {
      if (directions[i] == 0) continue;
      std::int32_t loK, hiK, loJ, hiJ;
      monotoneLeafBox(tree, data, k, i, &loK, &hiK);
      monotoneLeafBox(tree, data, j, i, &loJ, &hiJ);
      bool jBelowK = (hiJ + 1 == loK);
      bool kBelowJ = (hiK + 1 == loJ);
      if (!jBelowK && !kBelowJ) continue;
      if (!monotoneBoxesOverlap(tree, data, k, j, scratch.axes, i)) continue;
      hasConstrainedNeighbor = true;
      if (jSkipped) continue;
      double muj = mu[j];
      if (jBelowK) {
        if (directions[i] > 0) a = std::max(a, muj);
        else                   b = std::min(b, muj);
      } else {
        if (directions[i] > 0) b = std::min(b, muj);
        else                   a = std::max(a, muj);
      }
    }
  }
  *aOut = a;
  *bOut = b;
  *constrained = hasConstrainedNeighbor;
}

/// True when every leaf's value lies within its neighbor bounds - the monotone
/// feasibility invariant. MonotoneConstantGaussianLeaf::drawFromPriorForTree
/// calls it as the acceptance predicate of the constrained prior's rejection
/// sampler, so tol is part of a live draw law and not a test-only slack.
inline bool monotoneTreeIsFeasible(const Tree& tree, const ColumnStore& data,
                                   const std::int8_t* directions,
                                   const double* mu, double tol = 1e-9) {
  MonotoneNeighborScratch s;
  s.allBottoms.clear();
  tree.fillBottom(0, s.allBottoms);
  for (std::int32_t leaf : s.allBottoms) {
    if (tree.at(leaf).numObservations() == 0) continue;
    double a, b;
    bool constrained;
    monotoneNeighborBounds(tree, data, directions, s.allBottoms, leaf, mu,
                           nullptr, 0, s, &a, &b, &constrained);
    if (mu[leaf] < a - tol || mu[leaf] > b + tol) return false;
  }
  return true;
}

struct MonotoneConstantGaussianLeaf {
  static constexpr bool hasVectorParams = false;
  static constexpr bool hasFunctionParams = false;

  double scale = 0.0;
  const ColumnStore* data = nullptr;
  // per-column {-1,0,+1}, owned: the direction spec is copied at construction
  // (the bridge's parsed vector does not outlive the create call) and read
  // every sweep by the geometry
  std::vector<std::int8_t> directions;
  // c^2 = pi / (pi - 1): a constrained leaf uses effective k / c so its post-
  // truncation (skew-normal) marginal variance matches the unconstrained
  // sigma_mu^2 (design section 6).
  double cInflation = 1.0;
  mutable MonotoneNeighborScratch scratch;

  double priorSd(double k, bool constrained) const {
    return (constrained ? cInflation : 1.0) * scale / k;
  }
  void posterior(double sumWeights, double sumWeightedResponse,
                 double residualVariance, double priorStdDev, double* mean,
                 double* stdDev) const {
    double priorPrecision = 1.0 / (priorStdDev * priorStdDev);
    double posteriorPrecision = sumWeights / residualVariance;
    double variance = 1.0 / (priorPrecision + posteriorPrecision);
    *stdDev = std::sqrt(variance);
    *mean = (sumWeightedResponse / residualVariance) * variance;
  }

  // ---- ScalarLeafModel surface (never on the constrained hot path) ---------
  double logIntegratedLikelihood(double k, double residualVariance,
                                 double sumWeights,
                                 double sumWeightedResponse) const {
    return ConstantGaussianLeaf{scale}.logIntegratedLikelihood(
      k, residualVariance, sumWeights, sumWeightedResponse);
  }
  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* w, double k,
                                        double residualVariance,
                                        std::int32_t nodeIndex) const {
    return ConstantGaussianLeaf{scale}.logIntegratedLikelihoodForNode(
      tree, y, w, k, residualVariance, nodeIndex);
  }
  double drawFromPosteriorForNode(ext_rng* rng, const Tree& tree, double k,
                                  double residualVariance,
                                  std::int32_t nodeIndex) const {
    return ConstantGaussianLeaf{scale}.drawFromPosteriorForNode(
      rng, tree, k, residualVariance, nodeIndex);
  }
  double drawFromPrior(ext_rng* rng, double k) const {
    return ConstantGaussianLeaf{scale}.drawFromPrior(rng, k);
  }

  // Draw from N(m, s^2) truncated to (a, b], using the sd-1 primitives; an
  // unbounded side passes through, and a numerically empty interval falls back
  // to the clamped `fallback`.
  static double drawTruncatedNormal(ext_rng* rng, double m, double s, double a,
                                    double b, double fallback) {
    if (!std::isfinite(a) && !std::isfinite(b))
      return m + s * ext_rng_simulateStandardNormal(rng);
    if (!std::isfinite(a))
      return ext_rng_simulateUpperTruncatedNormal(rng, m, s, b);
    if (!std::isfinite(b))
      return ext_rng_simulateLowerTruncatedNormal(rng, m, s, a);
    double z = ext_rng_simulateTruncatedNormalScale1(rng, m / s, a / s, b / s);
    return std::isnan(z) ? std::clamp(fallback, a, b) : s * z;
  }

  // ---- the constrained leaf-parameter draw (TreeDrawLeafModel) -------------
  // One leaf's truncated-normal full conditional, shared by the tree-wide sweep
  // and the birth/death redraws. Empty leaves get mu = 0. `bottoms` is the leaf
  // set monotoneNeighborBounds walks for frozen neighbors; skip/numSkip exclude
  // co-drawn siblings. Fixed k: the truncated draw carries no clean chi-squared
  // statistic (design section 6), so no sumSquaredParams accumulation.
  void drawOneLeaf(ext_rng* rng, const Tree& tree, std::int32_t leaf,
                   const std::vector<std::int32_t>& bottoms,
                   const std::int32_t* skip, std::size_t numSkip, double k,
                   double residualVariance, double* mu) const {
    const Node& node = tree.at(leaf);
    if (node.numObservations() == 0) {
      mu[leaf] = 0.0;
      return;
    }
    double a, b;
    bool constrained;
    monotoneNeighborBounds(tree, *data, directions.data(), bottoms, leaf, mu,
                           skip, numSkip, scratch, &a, &b, &constrained);
    double sd = priorSd(k, constrained);
    double m, s;
    posterior(node.sumWeights, node.sumWeightedResponse, residualVariance, sd,
              &m, &s);
    mu[leaf] = constrained ? drawTruncatedNormal(rng, m, s, a, b, mu[leaf])
                           : m + s * ext_rng_simulateStandardNormal(rng);
  }

  // Sequential single-site truncated-normal Gibbs over the tree's leaves,
  // updating the SURVIVING mu block in place (the move phase keeps it feasible
  // through births/deaths).
  void drawParametersForTree(ext_rng* rng, const Tree& tree,
                             const std::vector<std::int32_t>& bottoms, double k,
                             double residualVariance, double* mu) const {
    for (std::int32_t leaf : bottoms)
      drawOneLeaf(rng, tree, leaf, bottoms, nullptr, 0, k, residualVariance, mu);
  }

  // Exact draw of a tree's leaf vector from the CONSTRAINED prior, by
  // rejection: the leaves' independent (c-inflated) prior marginals, accepted
  // only when the whole vector lands in the monotone cone. A sequential sweep
  // of the truncated full conditionals is NOT the joint truncated law, so
  // rejection is the exact route. Acceptance runs ~1/L! over L leaves chained
  // on a constrained axis - measured 6.5% per try with every axis constrained,
  // 63% with one of five - and prior-drawn trees average 2.5 leaves, so the
  // cap only catches a pathologically deep structure. Empty leaves take mu = 0
  // as everywhere else, flagged by a zero prior sd.
  static constexpr int priorDrawMaxAttempts = 1000000;
  bool drawFromPriorForTree(ext_rng* rng, const Tree& tree,
                            const std::vector<std::int32_t>& bottoms, double k,
                            double* mu) const {
    std::vector<double> priorSds(bottoms.size(), 0.0);
    for (std::size_t i = 0; i < bottoms.size(); ++i) {
      if (tree.at(bottoms[i]).numObservations() == 0) continue;
      double a, b;
      bool constrained;
      monotoneNeighborBounds(tree, *data, directions.data(), bottoms,
                             bottoms[i], mu, nullptr, 0, scratch, &a, &b,
                             &constrained);
      priorSds[i] = priorSd(k, constrained);
    }
    for (int attempt = 0; attempt < priorDrawMaxAttempts; ++attempt) {
      for (std::size_t i = 0; i < bottoms.size(); ++i)
        mu[bottoms[i]] = priorSds[i] == 0.0
          ? 0.0
          : priorSds[i] * ext_rng_simulateStandardNormal(rng);
      if (monotoneTreeIsFeasible(tree, *data, directions.data(), mu))
        return true;
    }
    return false;
  }

  // Redraw of an accepted death's merged leaf: its full conditional
  // truncated to the frozen neighbor bounds. Writing a proper draw (not a
  // deterministic seed) keeps the collapsed move a valid block update so the
  // stationary tree/leaf posterior is exact.
  void redrawAfterDeath(ext_rng* rng, const Tree& tree, std::int32_t leaf,
                        double k, double residualVariance, double* mu) const {
    scratch.allBottoms.clear();
    tree.fillBottom(0, scratch.allBottoms);
    drawOneLeaf(rng, tree, leaf, scratch.allBottoms, nullptr, 0, k,
                residualVariance, mu);
  }

  // Redraw of an accepted birth's two children from their EXACT
  // constrained conditional posterior over the cone {aL<=mu_lower<=bL,
  // aR<=mu_upper<=bR, mu_lower<=mu_upper}: draw mu_upper from its marginal by
  // rejection off the frozen-bound truncated normal (accept in proportion to
  // the mu_lower mass it admits), then mu_lower | mu_upper. Exact, and cheap
  // (the acceptance averages ~1/2). A non-constrained-axis split leaves the two
  // independent, each its own 1-D truncated draw.
  void redrawAfterBirth(ext_rng* rng, const Tree& tree, std::int32_t parent,
                        double k, double residualVariance, double* mu) const {
    scratch.allBottoms.clear();
    tree.fillBottom(0, scratch.allBottoms);
    std::int32_t left = tree.at(parent).leftChild;
    std::int32_t splitVar = tree.at(parent).rule.variableIndex;
    if (directions[splitVar] == 0) {
      redrawLeafFree(rng, tree, left, left + 1, k, residualVariance, mu);
      redrawLeafFree(rng, tree, left + 1, left, k, residualVariance, mu);
      return;
    }
    std::int32_t lower = directions[splitVar] > 0 ? left : left + 1;
    std::int32_t upper = directions[splitVar] > 0 ? left + 1 : left;
    const Node& nodeL = tree.at(lower);
    const Node& nodeR = tree.at(upper);
    if (nodeL.numObservations() == 0 || nodeR.numObservations() == 0) {
      redrawLeafFree(rng, tree, lower, upper, k, residualVariance, mu);
      redrawLeafFree(rng, tree, upper, lower, k, residualVariance, mu);
      return;
    }
    double aL, bL, aR, bR;
    bool cL, cR;
    std::int32_t skipU = upper, skipL = lower;
    monotoneNeighborBounds(tree, *data, directions.data(), scratch.allBottoms,
                           lower, mu, &skipU, 1, scratch, &aL, &bL, &cL);
    monotoneNeighborBounds(tree, *data, directions.data(), scratch.allBottoms,
                           upper, mu, &skipL, 1, scratch, &aR, &bR, &cR);
    double sd = priorSd(k, true);
    double mL, sL, mR, sR;
    posterior(nodeL.sumWeights, nodeL.sumWeightedResponse, residualVariance, sd,
              &mL, &sL);
    posterior(nodeR.sumWeights, nodeR.sumWeightedResponse, residualVariance, sd,
              &mR, &sR);
    double lowR = std::max(aR, aL);
    double hiL = std::isfinite(bR) ? std::min(bL, bR) : bL;  // max reachable min(bL, muR)
    double acceptMax = normalMass(aL, hiL, mL, sL);
    double muUpper = mu[upper], muLower = mu[lower];
    if (acceptMax > 0.0) {
      for (int tries = 0; tries < 100; ++tries) {
        double candidate = drawTruncatedNormal(rng, mR, sR, lowR, bR, muUpper);
        double admitted =
          normalMass(aL, std::isfinite(bL) ? std::min(bL, candidate) : candidate,
                     mL, sL);
        if (ext_rng_simulateContinuousUniform(rng) * acceptMax <= admitted) {
          muUpper = candidate;
          break;
        }
      }
    }
    double upperL = std::isfinite(bL) ? std::min(bL, muUpper) : muUpper;
    muLower = drawTruncatedNormal(rng, mL, sL, aL, upperL, muLower);
    mu[lower] = muLower;
    mu[upper] = muUpper;
  }

  // one child's independent 1-D truncated conditional (unconstrained split axis)
  void redrawLeafFree(ext_rng* rng, const Tree& tree, std::int32_t leaf,
                      std::int32_t sibling, double k, double residualVariance,
                      double* mu) const {
    drawOneLeaf(rng, tree, leaf, scratch.allBottoms, &sibling, 1, k,
                residualVariance, mu);
  }

  // ---- the constrained (truncated) birth/death marginal (ParamScoringLeafModel)
  double logLikelihoodForBranchWithParams(const Tree& tree,
                                          std::int32_t branchIndex,
                                          const double*, const double*,
                                          double k, double residualVariance,
                                          const double* mu) const {
    scratch.branch.clear();
    tree.fillBottom(branchIndex, scratch.branch);
    // The empty-leaf veto is NOT applied here even though this branch owns the
    // whole marginal: it is the caller's branch RANK (moves.hpp), taken over
    // the same leaves for every leaf model, so that two vetoed branches stay
    // comparable (docs/design/empty-leaf-veto.md). What survives below is the
    // FEASIBILITY sentinel, a different -HUGE_VAL: an empty constraint cone
    // has no draw at all. A zero-weight leaf's conditional is its prior
    // truncated to the neighbor bounds, which is finite and is what the
    // constrained forest should sit at where it has nothing to fit.
    scratch.allBottoms.clear();
    tree.fillBottom(0, scratch.allBottoms);

    std::size_t numLeaves = scratch.branch.size();
    if (numLeaves == 2) {
      std::int32_t splitVar = tree.at(branchIndex).rule.variableIndex;
      if (directions[splitVar] != 0) {
        std::int32_t left = tree.at(branchIndex).leftChild;
        std::int32_t lower = directions[splitVar] > 0 ? left : left + 1;
        std::int32_t upper = directions[splitVar] > 0 ? left + 1 : left;
        return twoLeafCoupledLogMarginal(tree, lower, upper, mu, k,
                                         residualVariance);
      }
    }
    // one leaf (birth-old / death-new), two on an unconstrained split, or the
    // >2 of a change/swap (outside the v1 move set): a product of independent
    // 1-D constrained marginals over the frozen outside neighbors. Exact for
    // the birth/death cases; for change/swap it still rejects an infeasible
    // arrangement (any empty bound gives the sentinel).
    double sum = 0.0;
    for (std::int32_t leaf : scratch.branch) {
      double term =
        oneLeafLogMarginal(tree, leaf, mu, k, residualVariance, scratch.branch);
      if (term == -HUGE_VAL) return -HUGE_VAL;
      sum += term;
    }
    return sum;
  }

  double oneLeafLogMarginal(const Tree& tree, std::int32_t leaf,
                            const double* mu, double k, double residualVariance,
                            const std::vector<std::int32_t>& skip) const {
    double a, b;
    bool constrained;
    monotoneNeighborBounds(tree, *data, directions.data(), scratch.allBottoms, leaf, mu,
                           skip.data(), skip.size(), scratch, &a, &b,
                           &constrained);
    if (a > b) return -HUGE_VAL;
    const Node& node = tree.at(leaf);
    double kEff = constrained ? k / cInflation : k;
    double base = ConstantGaussianLeaf{scale}.logIntegratedLikelihood(
      kEff, residualVariance, node.sumWeights, node.sumWeightedResponse);
    if (!constrained) return base;
    double sd = priorSd(k, true);
    double m, s;
    posterior(node.sumWeights, node.sumWeightedResponse, residualVariance, sd,
              &m, &s);
    double postMass = normalMass(a, b, m, s);
    double priorMass = normalMass(a, b, 0.0, sd);
    if (postMass <= 0.0 || priorMass <= 0.0) return -HUGE_VAL;
    return base + std::log(postMass) - std::log(priorMass);
  }

  // The bivariate constrained-axis marginal: the cone is
  // {aL <= mu_lower <= bL, aR <= mu_upper <= bR, mu_lower <= mu_upper}, so the
  // touched-leaf integral and its prior normalizer d_* are each one 1-D
  // quadrature (standardized on the upper leaf), not a grid (design section 4).
  double twoLeafCoupledLogMarginal(const Tree& tree, std::int32_t lower,
                                   std::int32_t upper, const double* mu, double k,
                                   double residualVariance) const {
    double aL, bL, aR, bR;
    bool cL, cR;
    std::int32_t skipUpper = upper, skipLower = lower;
    monotoneNeighborBounds(tree, *data, directions.data(), scratch.allBottoms, lower, mu,
                           &skipUpper, 1, scratch, &aL, &bL, &cL);
    monotoneNeighborBounds(tree, *data, directions.data(), scratch.allBottoms, upper, mu,
                           &skipLower, 1, scratch, &aR, &bR, &cR);
    const Node& nodeL = tree.at(lower);
    const Node& nodeR = tree.at(upper);
    double kEff = k / cInflation;
    double sd = priorSd(k, true);
    double mL, sL, mR, sR;
    posterior(nodeL.sumWeights, nodeL.sumWeightedResponse, residualVariance, sd,
              &mL, &sL);
    posterior(nodeR.sumWeights, nodeR.sumWeightedResponse, residualVariance, sd,
              &mR, &sR);
    double base =
      ConstantGaussianLeaf{scale}.logIntegratedLikelihood(
        kEff, residualVariance, nodeL.sumWeights, nodeL.sumWeightedResponse) +
      ConstantGaussianLeaf{scale}.logIntegratedLikelihood(
        kEff, residualVariance, nodeR.sumWeights, nodeR.sumWeightedResponse);
    double lowR = std::max(aR, aL);
    if (lowR > bR) return -HUGE_VAL;
    double numer = coneProbability(lowR, bR, aL, bL, mL, sL, mR, sR);
    double denom = coneProbability(lowR, bR, aL, bL, 0.0, sd, 0.0, sd);
    if (numer <= 0.0 || denom <= 0.0) return -HUGE_VAL;
    return base + std::log(numer) - std::log(denom);
  }

  // P(a? <= mu_lower <= mu_upper, mu_lower in [aL, bL], mu_upper in [lowR, hiR])
  // for independent mu_lower ~ N(mL, sL^2), mu_upper ~ N(mR, sR^2); standardized
  // on the upper leaf so the integrand is a standard-normal density times a
  // bounded CDF difference.
  double coneProbability(double lowR, double hiR, double aL, double bL,
                         double mL, double sL, double mR, double sR) const {
    double loU = std::max((lowR - mR) / sR, -38.0);
    double hiU = std::isfinite(hiR) ? std::min((hiR - mR) / sR, 38.0) : 38.0;
    double lowerL = std::isfinite(aL) ? gaussianCdf((aL - mL) / sL) : 0.0;
    return monotoneIntegrate(
      [&](double u) {
        double x = mR + sR * u;
        double upper = std::isfinite(bL) ? std::min(bL, x) : x;
        double inner = gaussianCdf((upper - mL) / sL) - lowerL;
        return (inner > 0.0 ? inner : 0.0) * gaussianPdf(u);
      },
      loU, hiU);
  }

  static double normalMass(double a, double b, double mean, double stdDev) {
    double hi = std::isfinite(b) ? gaussianCdf((b - mean) / stdDev) : 1.0;
    double lo = std::isfinite(a) ? gaussianCdf((a - mean) / stdDev) : 0.0;
    return hi - lo;
  }
};

static_assert(ScalarLeafModel<MonotoneConstantGaussianLeaf>);
static_assert(TreeDrawLeafModel<MonotoneConstantGaussianLeaf>);
static_assert(ParamScoringLeafModel<MonotoneConstantGaussianLeaf>);

// ---- Dense linear-algebra primitives shared by the conjugate leaf models ---
//
// LinearGaussianLeaf's ridge normal equations and GPGaussianLeaf's kernel
// solves run the identical dense factorization and triangular solves over
// row-major p x p storage; they live here so both leaves share one copy. The
// amplitude conditional (combiner.hpp) shares the storage convention through
// the square-root-free pair below.

/// In-place lower Cholesky of a symmetric positive-definite matrix; callers
/// guarantee definiteness (a ridge, or a nugget/noise diagonal), so there is
/// no failure path.
inline void choleskyDecompose(double* m, std::size_t p) {
  for (std::size_t j = 0; j < p; ++j) {
    double diagonal = m[j * p + j];
    for (std::size_t a = 0; a < j; ++a)
      diagonal -= m[j * p + a] * m[j * p + a];
    diagonal = std::sqrt(diagonal);
    m[j * p + j] = diagonal;
    for (std::size_t i = j + 1; i < p; ++i) {
      double value = m[i * p + j];
      for (std::size_t a = 0; a < j; ++a)
        value -= m[i * p + a] * m[j * p + a];
      m[i * p + j] = value / diagonal;
    }
  }
}

/// Forward solve L x = b in place, b supplied in x.
inline void solveLowerTriangular(const double* l, std::size_t p, double* x) {
  for (std::size_t i = 0; i < p; ++i) {
    double value = x[i];
    for (std::size_t a = 0; a < i; ++a) value -= l[i * p + a] * x[a];
    x[i] = value / l[i * p + i];
  }
}

/// Back solve L' x = b in place, b supplied in x.
inline void solveLowerTriangularTransposed(const double* l, std::size_t p,
                                           double* x) {
  for (std::size_t i = p; i > 0; --i) {
    double value = x[i - 1];
    for (std::size_t a = i; a < p; ++a) value -= l[a * p + (i - 1)] * x[a];
    x[i - 1] = value / l[(i - 1) * p + (i - 1)];
  }
}

/// In-place square-root-free (L D L') factorization of a symmetric positive-
/// definite matrix over the same row-major p x p storage as the Cholesky
/// above: the pivots d_j overwrite the diagonal and the UNIT lower triangle's
/// strict entries overwrite the strict lower triangle. Callers guarantee
/// definiteness, as there.
///
/// It exists beside choleskyDecompose because the two solve arithmetics differ
/// exactly where the amplitude conditional (combiner.hpp) cannot afford them
/// to. Through L L' a scalar system divides twice by sqrt(d), and
/// (x / sqrt(d)) / sqrt(d) is not x / d - measured, 500345 of a million random
/// (x, d) pairs. Through L D L' the unit triangles contribute nothing at p = 1
/// and nothing off the diagonal at an orthogonal design, so the solve reduces
/// to the ONE division per coordinate the scalar conditional writes, and a
/// q-variate draw over an orthogonal basis reproduces q scalar draws bitwise.
inline void unitLowerDecompose(double* m, std::size_t p) {
  for (std::size_t j = 0; j < p; ++j) {
    double pivot = m[j * p + j];
    for (std::size_t a = 0; a < j; ++a)
      pivot -= m[j * p + a] * m[j * p + a] * m[a * p + a];
    m[j * p + j] = pivot;
    for (std::size_t i = j + 1; i < p; ++i) {
      double value = m[i * p + j];
      for (std::size_t a = 0; a < j; ++a)
        value -= m[i * p + a] * m[j * p + a] * m[a * p + a];
      m[i * p + j] = value / pivot;
    }
  }
}

/// Forward solve L x = b in place against an unitLowerDecompose factorization,
/// b supplied in x. L is UNIT lower triangular, so the diagonal - which holds
/// the pivots, not ones - is read by neither solve, and neither divides.
inline void solveUnitLowerTriangular(const double* m, std::size_t p,
                                     double* x) {
  for (std::size_t i = 0; i < p; ++i) {
    double value = x[i];
    for (std::size_t a = 0; a < i; ++a) value -= m[i * p + a] * x[a];
    x[i] = value;
  }
}

/// Back solve L' x = b in place against the same factorization.
inline void solveUnitLowerTriangularTransposed(const double* m, std::size_t p,
                                               double* x) {
  for (std::size_t i = p; i > 0; --i) {
    double value = x[i - 1];
    for (std::size_t a = i; a < p; ++a) value -= m[a * p + (i - 1)] * x[a];
    x[i - 1] = value;
  }
}

/// Linear Gaussian leaf: each bottom node fits an intercept plus a linear
/// term in q designated ordinal predictor columns, standardized internally
/// to the training mean and sample sd. All q + 1 coefficients are iid
/// N(0, (scale / k)^2), so the marginal over them is a closed-form ridge
/// regression and the conjugate MH moves apply unchanged. The sufficient
/// statistics (U'WU, U'Wz, z'Wz) accumulate per call over the node's index
/// segment - the same O(n_leaf) walk as the constant leaf's variance,
/// widened by (q + 1)^2 - so move rollbacks need no cache invalidation.
/// Missing covariate values enter at the standardized mean (zero); rules on
/// the same column still route the missingness itself. Linear leaves read
/// raw predictor values through ColumnStore::rawColumn: the borrowed matrix
/// normally, gathered copies (with parent-derived standardization constants)
/// on a view built with the designation; a view without them is refused
/// upstream.
struct LinearGaussianLeaf {
  static constexpr bool hasVectorParams = true;
  static constexpr bool hasFunctionParams = false;
  /// Sufficient-statistic scratch lives on the stack, sized for
  /// maxNumCovariates; the factory rejects any designation with
  /// numCovariates > maxNumCovariates.
  static constexpr std::size_t maxNumCovariates = 8;

  double scale = 1.0;  // nodeScale / sqrt(numTrees)

  std::size_t numParams() const { return numCovariates_ + 1; }
  std::size_t numCovariates() const { return numCovariates_; }
  const std::vector<std::size_t>& covariateColumns() const { return columns_; }
  const std::vector<double>& covariateMeans() const { return means_; }
  const std::vector<double>& covariateSds() const { return sds_; }

  /// Gather and standardize the designated columns from the store's raw
  /// training values. Standardization constants come from the observed
  /// values alone; a constant (or all-missing) column keeps sd 1 and
  /// degrades to an extra intercept the ridge prior absorbs.
  void initialize(const ColumnStore& data, const std::size_t* columns,
                  std::size_t numColumns, std::size_t numChains = 1) {
    numCovariates_ = numColumns;
    columns_.assign(columns, columns + numColumns);
    statisticsCacheBudget_ =
      statisticsCacheTotalBudgetBytes / (numChains == 0 ? 1 : numChains);
    reinitialize(data);
  }

  /// Re-derive the standardization constants and covariates for the stored
  /// designation, for whole-data replacement (possibly a new number of
  /// observations): the analogue of setData rebuilding the cut grid. On a
  /// view the raw values come from the gathered copies and the constants
  /// from the parent's full data - the same calibration inheritance as the
  /// copied cut grid, so every fold runs the prior a full-data fit would.
  void reinitialize(const ColumnStore& data) {
    clearStatisticsCache();
    std::size_t numColumns = numCovariates_;
    numObservations_ = data.numObservations;
    means_.assign(numColumns, 0.0);
    sds_.assign(numColumns, 1.0);
    u_.resize(numObservations_ * numColumns);
    for (std::size_t j = 0; j < numColumns; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double mean, sd;
      if (!data.suppliedStandardization(columns_[j], &mean, &sd))
        standardizationMomentsForColumn(column, numObservations_, &mean, &sd);
      means_[j] = mean;
      sds_[j] = sd;
      for (std::size_t i = 0; i < numObservations_; ++i)
        u_[i * numColumns + j] =
          isNA(column[i]) ? 0.0 : (column[i] - mean) / sds_[j];
    }
    rebuildTestCovariates(data);
  }

  /// Regather the training covariates under the existing standardization
  /// constants, for the predictor-mutation surface (including rollback,
  /// which restores the old raw values): the prior's calibration is sticky
  /// across in-place value changes, the way refreshCutsForColumn keeps the
  /// cut count. Whole-data replacement re-initializes instead, refreshing
  /// the constants the way setData rebuilds the cut grid.
  void regatherTrainingCovariates(const ColumnStore& data) {
    clearStatisticsCache();
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      for (std::size_t i = 0; i < numObservations_; ++i)
        u_[i * numCovariates_ + j] =
          isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  /// Drop the crossproduct cache when U'WU's other inputs change: the case
  /// weights (setWeights) or the per-sweep Polya-Gamma refresh of a latent
  /// family. Covariate and whole-data mutations clear it through the two
  /// regather paths above.
  void invalidateStatistics() const { clearStatisticsCache(); }

  /// Regather the test covariates under the training standardization; called
  /// whenever the store's test data changes.
  void rebuildTestCovariates(const ColumnStore& data) {
    numTestObservations_ = data.numTestObservations;
    uTest_.assign(numTestObservations_ * numCovariates_, 0.0);
    if (numTestObservations_ == 0) return;
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawTestColumn(columns_[j]);
      if (column == nullptr) continue;
      double* u = uTest_.data() + j * numTestObservations_;
      for (std::size_t i = 0; i < numTestObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  double fitForObservation(const double* params, std::size_t i) const {
    double result = params[0];
    for (std::size_t j = 0; j < numCovariates_; ++j)
      result += params[j + 1] * u_[i * numCovariates_ + j];
    return result;
  }

  double fitForTestObservation(const double* params, std::size_t i) const {
    double result = params[0];
    for (std::size_t j = 0; j < numCovariates_; ++j)
      result += params[j + 1] * uTest_[i + j * numTestObservations_];
    return result;
  }

  /// log integral of prod N(z_i; b0 + b'u_i, sigma^2 / w_i) against the
  /// N(0, (scale / k)^2 I) prior on (b0, b), dropping the same terms the
  /// constant leaf drops (z'Wz and the pieces depending only on n and w,
  /// all additive over the member set, so the moves' before/after
  /// comparisons share them). With M = U'WU + tau sigma^2 I and
  /// tau = (k / scale)^2:
  ///   0.5 (q+1) log(tau sigma^2) - 0.5 log det M
  ///     + 0.5 b'M^-1 b / sigma^2,  b = U'Wz,
  /// which reduces exactly to the constant leaf's formula at q = 0.
  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* weights, double k,
                                        double residualVariance,
                                        int32_t nodeIndex) const {
    if (tree.at(nodeIndex).numObservations() == 0) return 0.0;

    std::size_t p = numParams();
    double crossproduct[maxStatisticSize], projection[maxNumCovariates + 1];
    accumulateNodeStatistics(tree, y, weights, nodeIndex, crossproduct,
                             projection);

    double ridge = (k / scale) * (k / scale) * residualVariance;
    for (std::size_t a = 0; a < p; ++a) crossproduct[a * p + a] += ridge;
    choleskyDecompose(crossproduct, p);

    double halfLogDet = 0.0;
    for (std::size_t a = 0; a < p; ++a)
      halfLogDet += std::log(crossproduct[a * p + a]);

    // b' M^-1 b through the forward half-solve alone
    solveLowerTriangular(crossproduct, p, projection);
    double quadraticForm = 0.0;
    for (std::size_t a = 0; a < p; ++a)
      quadraticForm += projection[a] * projection[a];

    return 0.5 * static_cast<double>(p) * std::log(ridge) - halfLogDet +
           0.5 * quadraticForm / residualVariance;
  }

  /// (b0, b) | z ~ N(M^-1 b, sigma^2 M^-1): with M = LL' and v = L^-1 b, the
  /// draw is L'^-1 (v + sigma eps) for eps standard normal, drawn in
  /// coordinate order (intercept first). Empty leaves zero the block without
  /// consuming generator draws, like the constant leaf's empty-leaf zero.
  void drawFromPosteriorForNode(ext_rng* rng, const Tree& tree,
                                const double* y, const double* weights,
                                double k, double residualVariance,
                                int32_t nodeIndex, double* out) const {
    std::size_t p = numParams();
    if (tree.at(nodeIndex).numObservations() == 0) {
      for (std::size_t a = 0; a < p; ++a) out[a] = 0.0;
      return;
    }

    double crossproduct[maxStatisticSize];
    accumulateNodeStatistics(tree, y, weights, nodeIndex, crossproduct, out);

    double ridge = (k / scale) * (k / scale) * residualVariance;
    for (std::size_t a = 0; a < p; ++a) crossproduct[a * p + a] += ridge;
    choleskyDecompose(crossproduct, p);

    solveLowerTriangular(crossproduct, p, out);
    double sigma = std::sqrt(residualVariance);
    for (std::size_t a = 0; a < p; ++a)
      out[a] += sigma * ext_rng_simulateStandardNormal(rng);
    solveLowerTriangularTransposed(crossproduct, p, out);
  }

  void drawFromPrior(ext_rng* rng, double k, double* out) const {
    std::size_t p = numParams();
    for (std::size_t a = 0; a < p; ++a)
      out[a] = (scale / k) * ext_rng_simulateStandardNormal(rng);
  }

private:
  static constexpr std::size_t maxStatisticSize =
    (maxNumCovariates + 1) * (maxNumCovariates + 1);

  /// One pass over the node's index segment: crossproduct receives the full
  /// symmetric U'WU (row-major (q+1) x (q+1), leading intercept column) and
  /// projection U'Wz. U'WU is served from the crossproduct cache when the
  /// node's member list matches an entry; the residual-dependent projection
  /// always rescans. A served value is bitwise the fresh scan's, since the
  /// entry was built by the identical fused loop over the same members,
  /// covariates, and weights.
  void accumulateNodeStatistics(const Tree& tree, const double* y,
                                const double* weights, int32_t nodeIndex,
                                double* crossproduct,
                                double* projection) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t p = numParams();
    double row[maxNumCovariates + 1];
    row[0] = 1.0;

    const double* cached = lookupCrossproduct(tree, node, nodeIndex);
    if (cached != nullptr) {
      for (std::size_t a = 0; a < p; ++a) projection[a] = 0.0;
      for (std::size_t m = node.begin; m < node.end; ++m) {
        std::size_t i = tree.indices[m];
        double w = weights == nullptr ? 1.0 : weights[i];
        double z = y[i];
        for (std::size_t j = 0; j < numCovariates_; ++j)
          row[j + 1] = u_[i * numCovariates_ + j];
        for (std::size_t a = 0; a < p; ++a) {
          double scaled = w * row[a];
          projection[a] += scaled * z;
        }
      }
      std::memcpy(crossproduct, cached, p * p * sizeof(double));
      return;
    }

    for (std::size_t a = 0; a < p * p; ++a) crossproduct[a] = 0.0;
    for (std::size_t a = 0; a < p; ++a) projection[a] = 0.0;
    for (std::size_t m = node.begin; m < node.end; ++m) {
      std::size_t i = tree.indices[m];
      double w = weights == nullptr ? 1.0 : weights[i];
      double z = y[i];
      for (std::size_t j = 0; j < numCovariates_; ++j)
        row[j + 1] = u_[i * numCovariates_ + j];
      for (std::size_t a = 0; a < p; ++a) {
        double scaled = w * row[a];
        projection[a] += scaled * z;
        for (std::size_t b = a; b < p; ++b)
          crossproduct[a * p + b] += scaled * row[b];
      }
    }
    for (std::size_t a = 0; a < p; ++a)
      for (std::size_t b = a + 1; b < p; ++b)
        crossproduct[b * p + a] = crossproduct[a * p + b];
    storeCrossproduct(tree, node, nodeIndex, crossproduct);
  }

  /// One leaf's cached U'WU (row-major p x p), tagged with the exact ordered
  /// member list that built it. Coherence: every lookup re-validates by
  /// comparing that list against tree.indices[begin..end], so any structural
  /// move or rejected-move rollback that alters membership fails the compare
  /// and rebuilds - the key is rollback-stable with no per-move hook. U'WU's
  /// other inputs (covariates, weights) are held fixed within a cache
  /// lifetime by clearing wholesale from the two covariate regathers and from
  /// invalidateStatistics (weights); within a sweep the score and draw phases
  /// see identical weights, so a served value is bitwise the fresh scan's.
  struct CachedNodeStatistics {
    std::vector<index_t> members;
    double crossproduct[maxStatisticSize];
  };
  struct TreeStatisticsCache {
    const Tree* tree = nullptr;
    std::vector<CachedNodeStatistics> nodes;  // arena-indexed
  };

  /// Leaves below this rescan; their U'WU is cheap and churns fast.
  static constexpr std::size_t minCachedLeafSize = 32;
  /// Byte ceiling over cached member lists, split across chains at initialize;
  /// when spent, further leaves rescan (still correct, just uncached).
  static constexpr std::size_t statisticsCacheTotalBudgetBytes =
    static_cast<std::size_t>(256) << 20;

  static std::size_t statisticsEntryBytes(std::size_t numObs) {
    return numObs * sizeof(index_t);
  }

  TreeStatisticsCache& statisticsCacheForTree(const Tree& tree) const {
    for (TreeStatisticsCache& cache : statisticsCaches_)
      if (cache.tree == &tree) return cache;
    statisticsCaches_.emplace_back();
    statisticsCaches_.back().tree = &tree;
    return statisticsCaches_.back();
  }

  void clearStatisticsCache() const {
    statisticsCaches_.clear();
    statisticsCacheUsedBytes_ = 0;
  }

  /// The node's cached U'WU when its ordered member list matches, else null.
  const double* lookupCrossproduct(const Tree& tree, const Node& node,
                                   int32_t nodeIndex) const {
    std::size_t numObs = node.numObservations();
    if (numObs < minCachedLeafSize) return nullptr;
    TreeStatisticsCache& cache = statisticsCacheForTree(tree);
    std::size_t index = static_cast<std::size_t>(nodeIndex);
    if (index >= cache.nodes.size()) return nullptr;
    const CachedNodeStatistics& entry = cache.nodes[index];
    if (entry.members.size() == numObs &&
        std::memcmp(entry.members.data(), tree.indices + node.begin,
                    numObs * sizeof(index_t)) == 0)
      return entry.crossproduct;
    return nullptr;
  }

  /// Record the freshly scanned U'WU for the node, subject to the byte budget.
  void storeCrossproduct(const Tree& tree, const Node& node, int32_t nodeIndex,
                         const double* crossproduct) const {
    std::size_t numObs = node.numObservations();
    if (numObs < minCachedLeafSize) return;
    std::size_t p = numParams();
    TreeStatisticsCache& cache = statisticsCacheForTree(tree);
    std::size_t index = static_cast<std::size_t>(nodeIndex);
    if (index >= cache.nodes.size()) cache.nodes.resize(index + 1);
    CachedNodeStatistics& entry = cache.nodes[index];
    std::size_t oldBytes = statisticsEntryBytes(entry.members.size());
    std::size_t newBytes = statisticsEntryBytes(numObs);
    if (statisticsCacheUsedBytes_ - oldBytes + newBytes >
        statisticsCacheBudget_) {
      if (!entry.members.empty()) {
        statisticsCacheUsedBytes_ -= oldBytes;
        entry.members.clear();
      }
      return;
    }
    statisticsCacheUsedBytes_ += newBytes - oldBytes;
    entry.members.assign(tree.indices + node.begin, tree.indices + node.end);
    std::memcpy(entry.crossproduct, crossproduct, p * p * sizeof(double));
  }

  std::size_t numCovariates_ = 0;
  std::size_t numObservations_ = 0;
  std::size_t numTestObservations_ = 0;
  std::vector<std::size_t> columns_;
  std::vector<double> means_, sds_;
  std::vector<double> u_;      // standardized, row-major n x q
  std::vector<double> uTest_;  // standardized, column-major numTest x q
  mutable std::vector<TreeStatisticsCache> statisticsCaches_;
  mutable std::size_t statisticsCacheUsedBytes_ = 0;
  std::size_t statisticsCacheBudget_ = statisticsCacheTotalBudgetBytes;
};

static_assert(VectorLeafModel<LinearGaussianLeaf>);

/// Gaussian-process leaf: each bottom node fits a smooth function of q
/// designated ordinal predictor columns, standardized internally exactly as
/// the linear leaf's. The leaf function over the standardized covariates u
/// is
///   f ~ GP(0, (scale / k)^2 C),
/// with C a squared-exponential correlation kernel (per-column lengthscales
/// theta plus a small fixed nugget for conditioning), so the prior variance
/// ties to k as the constant and linear leaves do and a constant kernel
/// recovers the constant leaf. The marginal over f is the closed form
/// log N(z; 0, sigma^2 W^-1 + (scale / k)^2 C) - one O(n_leaf^3) Cholesky -
/// so the conjugate MH moves apply unchanged. The drawn function values at
/// the leaf's training rows ARE the parameters and land directly in the
/// per-tree fits; test rows evaluate the conditional mean c(x*)' C^-1 f
/// through per-draw cached weights.
///
/// Leaves larger than maxLeafSize score and draw as constant leaves instead
/// (delegating to ConstantGaussianLeaf's exact math). The alternative -
/// vetoing (scoring -infinity) any leaf over the cap - would deadlock: every
/// tree starts as a single root leaf holding all observations, and a birth
/// splitting one over-cap leaf into two over-cap children can never be
/// accepted when both the current single over-cap leaf and the proposed pair
/// of over-cap children score -infinity. The fallback instead makes oversized
/// regions behave exactly
/// like constant-leaf BART - data-informed splits throughout - with the GP
/// refinement switching on once a leaf falls under the cap; the scoring rule
/// is a deterministic function of leaf membership, so the MH comparisons
/// stay coherent.
struct GPGaussianLeaf {
  static constexpr bool hasVectorParams = false;
  static constexpr bool hasFunctionParams = true;
  /// Bounds the replay predictor's stack scratch; the factory rejects any
  /// designation with numCovariates > maxNumCovariates.
  static constexpr std::size_t maxNumCovariates = maxFunctionLeafCovariates;
  /// Conditioning jitter on the correlation diagonal; not a modeling knob.
  static constexpr double nugget = 1.0e-6;

  double scale = 1.0;  // nodeScale / sqrt(numTrees)

  std::size_t numCovariates() const { return numCovariates_; }
  const std::vector<std::size_t>& covariateColumns() const { return columns_; }
  const std::vector<double>& covariateMeans() const { return means_; }
  const std::vector<double>& covariateSds() const { return sds_; }
  const std::vector<double>& lengthscales() const { return lengthscales_; }
  std::size_t maxLeafSize() const { return maxLeafSize_; }

  /// Gather and standardize the designated columns and fix the kernel
  /// lengthscales: the supplied per-column values (standardized scale) when
  /// non-null, otherwise the median pairwise-distance heuristic over the
  /// standardized training values.
  void initialize(const ColumnStore& data, const std::size_t* columns,
                  std::size_t numColumns, const double* lengthscales,
                  std::size_t maxLeafSize, std::size_t numChains = 1) {
    numCovariates_ = numColumns;
    columns_.assign(columns, columns + numColumns);
    maxLeafSize_ = maxLeafSize;
    kernelCacheBudget_ =
      kernelCacheTotalBudgetBytes / (numChains == 0 ? 1 : numChains);
    if (lengthscales != nullptr)
      suppliedLengthscales_.assign(lengthscales, lengthscales + numColumns);
    else
      suppliedLengthscales_.clear();
    reinitialize(data);
  }

  /// Re-derive standardization constants, covariates, and (when not
  /// supplied) lengthscales for the stored designation; the whole-data
  /// analogue of the rebuilt cut grid, sharing the linear leaf's calibration
  /// inheritance on views.
  void reinitialize(const ColumnStore& data) {
    clearKernelCaches();
    std::size_t numColumns = numCovariates_;
    numObservations_ = data.numObservations;
    means_.assign(numColumns, 0.0);
    sds_.assign(numColumns, 1.0);
    u_.resize(numObservations_ * numColumns);
    for (std::size_t j = 0; j < numColumns; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double mean, sd;
      if (!data.suppliedStandardization(columns_[j], &mean, &sd))
        standardizationMomentsForColumn(column, numObservations_, &mean, &sd);
      means_[j] = mean;
      sds_[j] = sd;
      double* u = u_.data() + j * numObservations_;
      for (std::size_t i = 0; i < numObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - mean) / sds_[j];
    }
    lengthscales_.assign(numColumns, 1.0);
    for (std::size_t j = 0; j < numColumns; ++j)
      lengthscales_[j] = suppliedLengthscales_.empty()
        ? medianPairwiseDistance(u_.data() + j * numObservations_,
                                 numObservations_)
        : suppliedLengthscales_[j];
    rebuildTestCovariates(data);
  }

  /// Regather the training covariates under the existing standardization
  /// constants and lengthscales, for the predictor-mutation surface; the
  /// prior's calibration is sticky across in-place value changes.
  void regatherTrainingCovariates(const ColumnStore& data) {
    clearKernelCaches();
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double* u = u_.data() + j * numObservations_;
      for (std::size_t i = 0; i < numObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  /// Regather the test covariates under the training standardization; called
  /// whenever the store's test data changes.
  void rebuildTestCovariates(const ColumnStore& data) {
    numTestObservations_ = data.numTestObservations;
    uTest_.assign(numTestObservations_ * numCovariates_, 0.0);
    if (numTestObservations_ == 0) return;
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawTestColumn(columns_[j]);
      if (column == nullptr) continue;
      double* u = uTest_.data() + j * numTestObservations_;
      for (std::size_t i = 0; i < numTestObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  /// Reset the per-draw prediction cache before one tree's parameter sweep;
  /// draws fill it node by node and the test-fit loop reads it while the
  /// tree's partitions are unchanged.
  void beginTreeDraw(const Tree& tree) const {
    nodeConstant_.assign(tree.nodes.size(), 0.0);
    nodeAlphaOffset_.assign(tree.nodes.size(), -1);
    alphaBuffer_.clear();
    evictStaleKernelEntries(tree);
  }

  /// score = -0.5 log det V + 0.5 log det(sigma^2 W^-1)
  /// - 0.5 (z' V^-1 z - z'Wz / sigma^2) with V = (scale / k)^2 C +
  /// sigma^2 W^-1, dropping the same terms the constant leaf drops (z'Wz is
  /// the model-free part of z'V^-1 z by Woodbury, added back explicitly
  /// here); a constant kernel reduces this to the constant leaf's formula by
  /// Sherman-Morrison, keeping the maxLeafSize fallback below coherent.
  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* weights, double k,
                                        double residualVariance,
                                        int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) return 0.0;
    if (numObs > maxLeafSize_) {
      ConstantGaussianLeaf fallback{scale};
      return fallback.logIntegratedLikelihoodForNode(tree, y, weights, k,
                                                     residualVariance,
                                                     nodeIndex);
    }
    if (weights != nullptr && anyZeroWeight(tree, node, weights))
      return logIntegratedLikelihoodOverPositiveWeights(
        tree, node, y, weights, k, residualVariance, nodeIndex, numObs);

    double s2 = (scale / k) * (scale / k);

    const CachedLeafKernel* entry =
      cachedKernelForNode(tree, node, nodeIndex, numObs);
    if (entry != nullptr) {
      cholV_.resize(numObs * numObs);
      for (std::size_t a = 0; a < numObs * numObs; ++a)
        cholV_[a] = entry->kernel[a] * s2;
    } else {
      gatherLeafCovariates(tree, node, numObs);
      buildKernel(numObs, cholV_);
      for (std::size_t a = 0; a < numObs * numObs; ++a) cholV_[a] *= s2;
    }
    double logDetNoise = 0.0, responseSumOfSquares = 0.0;
    for (std::size_t r = 0; r < numObs; ++r) {
      std::size_t i = tree.indices[node.begin + r];
      double w = weights == nullptr ? 1.0 : weights[i];
      double noise = residualVariance / w;
      logDetNoise += std::log(noise);
      cholV_[r * numObs + r] += noise;
      responseSumOfSquares += w * y[i] * y[i];
    }
    choleskyDecompose(cholV_.data(), numObs);

    double halfLogDetV = 0.0;
    for (std::size_t r = 0; r < numObs; ++r)
      halfLogDetV += std::log(cholV_[r * numObs + r]);

    vectorScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      vectorScratch_[r] = y[tree.indices[node.begin + r]];
    solveLowerTriangular(cholV_.data(), numObs, vectorScratch_.data());
    double quadraticForm = 0.0;
    for (std::size_t r = 0; r < numObs; ++r)
      quadraticForm += vectorScratch_[r] * vectorScratch_[r];

    return -halfLogDetV + 0.5 * logDetNoise -
           0.5 * (quadraticForm - responseSumOfSquares / residualVariance);
  }

  /// Draws f | z into fits at the leaf's member observations by Matheron's
  /// rule - f = f0 + s^2 C V^-1 (z - f0 - e0) with f0 ~ N(0, s^2 C) and
  /// e0 ~ N(0, sigma^2 W^-1) - which has exactly the posterior law and needs
  /// only the two factorizations already in hand. Consumes 2 n_leaf standard
  /// normals (n_leaf + n_positive when zero-weight members are present);
  /// empty leaves consume none and cache a zero constant. Caches the
  /// prediction weights alpha = C^-1 f for the node and returns the chi-k
  /// contribution (f' alpha over n_leaf coordinates).
  FunctionLeafDrawStats drawFromPosteriorForNode(
    ext_rng* rng, const Tree& tree, const double* y, const double* weights,
    double k, double residualVariance, int32_t nodeIndex,
    double* fits) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) {
      setConstantCache(nodeIndex, 0.0);
      return FunctionLeafDrawStats();
    }
    if (numObs > maxLeafSize_) {
      ConstantGaussianLeaf fallback{scale};
      double value = fallback.drawFromPosteriorForNode(rng, tree, k,
                                                       residualVariance,
                                                       nodeIndex);
      for (std::size_t m = node.begin; m < node.end; ++m)
        fits[tree.indices[m]] = value;
      setConstantCache(nodeIndex, value);
      return FunctionLeafDrawStats{value * value, 1.0};
    }
    if (weights != nullptr && anyZeroWeight(tree, node, weights))
      return drawFromPosteriorOverPositiveWeights(
        rng, tree, node, y, weights, k, residualVariance, nodeIndex, fits);

    double s = scale / k, s2 = s * s;

    const double* kernel;
    const double* cholK;
    kernelAndFactorForNode(tree, node, nodeIndex, numObs, &kernel, &cholK);
    cholV_.resize(numObs * numObs);
    for (std::size_t a = 0; a < numObs * numObs; ++a)
      cholV_[a] = s2 * kernel[a];
    for (std::size_t r = 0; r < numObs; ++r) {
      double w =
        weights == nullptr ? 1.0 : weights[tree.indices[node.begin + r]];
      cholV_[r * numObs + r] += residualVariance / w;
    }
    choleskyDecompose(cholV_.data(), numObs);

    // f0 = s L_C eps, drawn first in row order, then e0 row by row
    epsScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      epsScratch_[r] = ext_rng_simulateStandardNormal(rng);
    fScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a <= r; ++a)
        value += cholK[r * numObs + a] * epsScratch_[a];
      fScratch_[r] = s * value;
    }
    vectorScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double w =
        weights == nullptr ? 1.0 : weights[tree.indices[node.begin + r]];
      vectorScratch_[r] = y[tree.indices[node.begin + r]] - fScratch_[r] -
                          std::sqrt(residualVariance / w) *
                            ext_rng_simulateStandardNormal(rng);
    }
    solveLowerTriangular(cholV_.data(), numObs, vectorScratch_.data());
    solveLowerTriangularTransposed(cholV_.data(), numObs, vectorScratch_.data());
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a < numObs; ++a)
        value += kernel[r * numObs + a] * vectorScratch_[a];
      fScratch_[r] += s2 * value;
      fits[tree.indices[node.begin + r]] = fScratch_[r];
    }

    return cacheAlphaForNode(nodeIndex, numObs, cholK);
  }

  /// Prior draw f = s L_C eps into fits, with the same prediction cache and
  /// chi-k accounting as the posterior draw; over-cap leaves draw the
  /// constant leaf's prior. Consumes n_leaf standard normals (one over cap,
  /// none when empty).
  FunctionLeafDrawStats drawFromPriorForNode(ext_rng* rng, const Tree& tree,
                                             double k, int32_t nodeIndex,
                                             double* fits) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) {
      setConstantCache(nodeIndex, 0.0);
      return FunctionLeafDrawStats();
    }
    if (numObs > maxLeafSize_) {
      ConstantGaussianLeaf fallback{scale};
      double value = fallback.drawFromPrior(rng, k);
      for (std::size_t m = node.begin; m < node.end; ++m)
        fits[tree.indices[m]] = value;
      setConstantCache(nodeIndex, value);
      return FunctionLeafDrawStats{value * value, 1.0};
    }

    double s = scale / k;

    const double* kernel;
    const double* cholK;
    kernelAndFactorForNode(tree, node, nodeIndex, numObs, &kernel, &cholK);

    epsScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      epsScratch_[r] = ext_rng_simulateStandardNormal(rng);
    fScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a <= r; ++a)
        value += cholK[r * numObs + a] * epsScratch_[a];
      fScratch_[r] = s * value;
      fits[tree.indices[node.begin + r]] = fScratch_[r];
    }

    return cacheAlphaForNode(nodeIndex, numObs, cholK);
  }

  /// Append one leaf's saved side-channel block (the
  /// computeFunctionBlockOffsets layout) from the live draw cache: the
  /// cached alpha plus the members' plain standardized rows, in member
  /// order; constant-valued nodes (over-cap, empty) append [0, constant].
  /// Valid under the same freshness condition as the test-fit evaluation.
  void appendLeafBlockFromCache(const Tree& tree, int32_t nodeIndex,
                                std::vector<double>& blocks) const {
    std::ptrdiff_t offset =
      nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)];
    if (offset < 0) {
      blocks.push_back(0.0);
      blocks.push_back(nodeConstant_[static_cast<std::size_t>(nodeIndex)]);
      return;
    }
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    blocks.push_back(static_cast<double>(numObs));
    const double* alpha = alphaBuffer_.data() + offset;
    blocks.insert(blocks.end(), alpha, alpha + numObs);
    appendLeafRows(tree, node, numObs, blocks);
  }

  /// The same block recomputed from a tree's persisted fits (the fits ARE
  /// the parameters): over-cap leaves append their constant (fits are
  /// uniform there), gp leaves refactor alpha = C^-1 f against the current
  /// covariates. Serves flatten and live prediction between runs, when no
  /// draw cache is fresh.
  void appendLeafBlock(const Tree& tree, int32_t nodeIndex,
                       const double* fits, std::vector<double>& blocks) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) {
      blocks.push_back(0.0);
      blocks.push_back(0.0);
      return;
    }
    if (numObs > maxLeafSize_) {
      blocks.push_back(0.0);
      blocks.push_back(fits[tree.indices[node.begin]]);
      return;
    }

    gatherLeafCovariates(tree, node, numObs);
    buildKernel(numObs, kernel_);
    cholK_.assign(kernel_.begin(),
                  kernel_.begin() +
                    static_cast<std::ptrdiff_t>(numObs * numObs));
    choleskyDecompose(cholK_.data(), numObs);

    blocks.push_back(static_cast<double>(numObs));
    std::size_t alphaStart = blocks.size();
    for (std::size_t r = 0; r < numObs; ++r)
      blocks.push_back(fits[tree.indices[node.begin + r]]);
    solveLowerTriangular(cholK_.data(), numObs, blocks.data() + alphaStart);
    solveLowerTriangularTransposed(cholK_.data(), numObs,
                             blocks.data() + alphaStart);
    appendLeafRows(tree, node, numObs, blocks);
  }

  /// Conditional mean c(x*)' alpha at one test row from the node's cached
  /// draw; constant-valued nodes (over-cap, empty) return their constant.
  /// Valid only between the node's draw and the next structure change - the
  /// alpha weights pair with the member ordering at draw time.
  double fitForTestObservationForNode(const Tree& tree, int32_t nodeIndex,
                                      std::size_t testIndex) const {
    std::ptrdiff_t offset = nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)];
    if (offset < 0) return nodeConstant_[static_cast<std::size_t>(nodeIndex)];
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    const double* alpha = alphaBuffer_.data() + offset;
    double result = 0.0;
    for (std::size_t r = 0; r < numObs; ++r) {
      std::size_t obs = tree.indices[node.begin + r];
      double distanceSq = 0.0;
      for (std::size_t j = 0; j < numCovariates_; ++j) {
        double difference = (uTest_[testIndex + j * numTestObservations_] -
                             u_[obs + j * numObservations_]) /
                            lengthscales_[j];
        distanceSq += difference * difference;
      }
      result += std::exp(-0.5 * distanceSq) * alpha[r];
    }
    return result;
  }

private:
  /// Median of pairwise absolute differences on the standardized scale,
  /// deterministically subsampled to at most 512 rows; degenerate columns
  /// fall back to 1.
  static double medianPairwiseDistance(const double* u, std::size_t n) {
    std::size_t stride = 1 + (n - 1) / 512;
    std::vector<double> values;
    for (std::size_t i = 0; i < n; i += stride) values.push_back(u[i]);
    if (values.size() < 2) return 1.0;
    std::vector<double> distances;
    distances.reserve(values.size() * (values.size() - 1) / 2);
    for (std::size_t a = 0; a < values.size(); ++a)
      for (std::size_t b = a + 1; b < values.size(); ++b)
        distances.push_back(std::fabs(values[a] - values[b]));
    std::sort(distances.begin(), distances.end());
    std::size_t count = distances.size();
    double median = count % 2 == 1
      ? distances[count / 2]
      : 0.5 * (distances[count / 2 - 1] + distances[count / 2]);
    return median > 0.0 ? median : 1.0;
  }

  /// Leaf rows' standardized covariates pre-divided by the lengthscales,
  /// row-major numObs x q, so the kernel is exp(-0.5 ||row_r - row_c||^2).
  void gatherLeafCovariates(const Tree& tree, const Node& node,
                            std::size_t numObs) const {
    leafU_.resize(numObs * numCovariates_);
    for (std::size_t r = 0; r < numObs; ++r) {
      std::size_t obs = tree.indices[node.begin + r];
      for (std::size_t j = 0; j < numCovariates_; ++j)
        leafU_[r * numCovariates_ + j] =
          u_[obs + j * numObservations_] / lengthscales_[j];
    }
  }

  /// Squared-exponential correlation plus the nugget diagonal over the
  /// gathered leaf rows.
  void buildKernel(std::size_t numObs, std::vector<double>& out) const {
    out.resize(numObs * numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      out[r * numObs + r] = 1.0 + nugget;
      for (std::size_t c = r + 1; c < numObs; ++c) {
        double distanceSq = 0.0;
        for (std::size_t j = 0; j < numCovariates_; ++j) {
          double difference = leafU_[r * numCovariates_ + j] -
                              leafU_[c * numCovariates_ + j];
          distanceSq += difference * difference;
        }
        double value = std::exp(-0.5 * distanceSq);
        out[r * numObs + c] = value;
        out[c * numObs + r] = value;
      }
    }
  }

  void setConstantCache(int32_t nodeIndex, double value) const {
    nodeConstant_[static_cast<std::size_t>(nodeIndex)] = value;
    nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)] = -1;
  }

  /// Append the members' plain standardized covariate rows (row-major
  /// numObs x q, member order; lengthscales NOT baked in - replays divide)
  /// to a saved side-channel block.
  void appendLeafRows(const Tree& tree, const Node& node, std::size_t numObs,
                      std::vector<double>& blocks) const {
    for (std::size_t r = 0; r < numObs; ++r) {
      std::size_t obs = tree.indices[node.begin + r];
      for (std::size_t j = 0; j < numCovariates_; ++j)
        blocks.push_back(u_[obs + j * numObservations_]);
    }
  }

  /// alpha = C^-1 f through the correlation Cholesky, appended to the
  /// per-tree buffer; returns the chi-k contribution f' alpha.
  FunctionLeafDrawStats cacheAlphaForNode(int32_t nodeIndex,
                                          std::size_t numObs,
                                          const double* cholK) const {
    std::size_t offset = alphaBuffer_.size();
    alphaBuffer_.insert(alphaBuffer_.end(), fScratch_.begin(),
                        fScratch_.begin() +
                          static_cast<std::ptrdiff_t>(numObs));
    double* alpha = alphaBuffer_.data() + offset;
    solveLowerTriangular(cholK, numObs, alpha);
    solveLowerTriangularTransposed(cholK, numObs, alpha);
    nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)] =
      static_cast<std::ptrdiff_t>(offset);
    double sumSquared = 0.0;
    for (std::size_t r = 0; r < numObs; ++r)
      sumSquared += fScratch_[r] * alpha[r];
    return FunctionLeafDrawStats{sumSquared,
                                 static_cast<double>(numObs)};
  }

  /// One leaf's cached correlation kernel (nugget included) and its lower
  /// Cholesky factor, tagged with the exact member list that built it. The
  /// kernel depends only on membership order and the fixed lengthscales -
  /// not sigma, k, or the response - so entries survive across sweeps until
  /// an accepted move re-routes observations; every lookup re-validates by
  /// comparing the member list, making a hit bitwise identical to a fresh
  /// build with no invalidation hooks to miss. A vacant slot has an empty
  /// member list.
  struct CachedLeafKernel {
    std::vector<index_t> members;
    std::vector<double> kernel;
    std::vector<double> cholK;  // filled lazily by the first draw
  };
  struct TreeKernelCache {
    const Tree* tree = nullptr;
    std::vector<CachedLeafKernel> nodes;  // arena-indexed
  };

  /// Leaves below this recompute; their kernels are cheap and churn fast.
  static constexpr std::size_t minCachedLeafSize = 32;
  /// Hard byte ceiling over all cached kernels, split evenly across chains
  /// at initialize; when spent, further leaves recompute (still correct,
  /// just uncached).
  static constexpr std::size_t kernelCacheTotalBudgetBytes =
    static_cast<std::size_t>(256) << 20;

  static std::size_t kernelEntryBytes(std::size_t numObs) {
    return numObs * sizeof(index_t) +
           2 * numObs * numObs * sizeof(double);
  }

  TreeKernelCache& kernelCacheForTree(const Tree& tree) const {
    for (TreeKernelCache& cache : kernelCaches_)
      if (cache.tree == &tree) return cache;
    kernelCaches_.emplace_back();
    kernelCaches_.back().tree = &tree;
    return kernelCaches_.back();
  }

  void clearKernelCaches() const {
    kernelCaches_.clear();
    kernelCacheUsedBytes_ = 0;
  }

  /// The validated cache entry for a leaf, building the kernel into a fresh
  /// entry when none matches and the budget allows; null when the leaf is
  /// too small to bother or the budget refuses, in which case the caller
  /// computes into scratch. The returned pointer is valid until the next
  /// cache mutation, which no caller spans.
  CachedLeafKernel* cachedKernelForNode(const Tree& tree, const Node& node,
                                        int32_t nodeIndex,
                                        std::size_t numObs) const {
    if (numObs < minCachedLeafSize) return nullptr;
    TreeKernelCache& cache = kernelCacheForTree(tree);
    std::size_t index = static_cast<std::size_t>(nodeIndex);
    if (index >= cache.nodes.size()) cache.nodes.resize(index + 1);
    CachedLeafKernel& entry = cache.nodes[index];
    if (entry.members.size() == numObs &&
        std::memcmp(entry.members.data(), tree.indices + node.begin,
                    numObs * sizeof(index_t)) == 0)
      return &entry;

    std::size_t oldBytes = kernelEntryBytes(entry.members.size());
    std::size_t newBytes = kernelEntryBytes(numObs);
    if (kernelCacheUsedBytes_ - oldBytes + newBytes > kernelCacheBudget_) {
      if (!entry.members.empty()) {
        kernelCacheUsedBytes_ -= oldBytes;
        entry = CachedLeafKernel();
      }
      return nullptr;
    }
    kernelCacheUsedBytes_ += newBytes - oldBytes;
    entry = CachedLeafKernel();  // release any stale capacity
    entry.members.assign(tree.indices + node.begin, tree.indices + node.end);
    gatherLeafCovariates(tree, node, numObs);
    buildKernel(numObs, entry.kernel);
    return &entry;
  }

  /// The leaf's correlation kernel and lower Cholesky factor for a draw,
  /// served from the cache when possible and computed into the scratch
  /// buffers otherwise; either way the values are bitwise those of a fresh
  /// build.
  void kernelAndFactorForNode(const Tree& tree, const Node& node,
                              int32_t nodeIndex, std::size_t numObs,
                              const double** kernel,
                              const double** cholK) const {
    CachedLeafKernel* entry =
      cachedKernelForNode(tree, node, nodeIndex, numObs);
    if (entry != nullptr) {
      if (entry->cholK.empty()) {
        entry->cholK.assign(entry->kernel.begin(), entry->kernel.end());
        choleskyDecompose(entry->cholK.data(), numObs);
      }
      *kernel = entry->kernel.data();
      *cholK = entry->cholK.data();
      return;
    }
    gatherLeafCovariates(tree, node, numObs);
    buildKernel(numObs, kernel_);
    cholK_.assign(kernel_.begin(),
                  kernel_.begin() +
                    static_cast<std::ptrdiff_t>(numObs * numObs));
    choleskyDecompose(cholK_.data(), numObs);
    *kernel = kernel_.data();
    *cholK = cholK_.data();
  }

  /// True when any member observation carries a zero weight. Those rows
  /// have infinite noise variance, so they contribute no likelihood and
  /// the leaf takes the positive-subset paths below; the constant and
  /// linear leaves get the same behavior for free by multiplying by w.
  bool anyZeroWeight(const Tree& tree, const Node& node,
                     const double* weights) const {
    for (std::size_t m = node.begin; m < node.end; ++m)
      if (weights[tree.indices[m]] == 0.0) return true;
    return false;
  }

  /// The marginal over only the positive-weight members - the exact limit
  /// of infinite noise variance on the rest. Still a deterministic
  /// function of leaf membership and the weights, so MH comparisons stay
  /// coherent; a leaf with no positive-weight members scores 0 like an
  /// empty leaf.
  double logIntegratedLikelihoodOverPositiveWeights(
    const Tree& tree, const Node& node, const double* y,
    const double* weights, double k, double residualVariance,
    int32_t nodeIndex, std::size_t numObs) const {
    positiveScratch_.clear();
    for (std::size_t r = 0; r < numObs; ++r)
      if (weights[tree.indices[node.begin + r]] > 0.0)
        positiveScratch_.push_back(r);
    std::size_t numPos = positiveScratch_.size();
    if (numPos == 0) return 0.0;

    double s2 = (scale / k) * (scale / k);

    const CachedLeafKernel* entry =
      cachedKernelForNode(tree, node, nodeIndex, numObs);
    const double* kernel;
    if (entry != nullptr) {
      kernel = entry->kernel.data();
    } else {
      gatherLeafCovariates(tree, node, numObs);
      buildKernel(numObs, kernel_);
      kernel = kernel_.data();
    }
    cholV_.resize(numPos * numPos);
    for (std::size_t a = 0; a < numPos; ++a)
      for (std::size_t b = 0; b < numPos; ++b)
        cholV_[a * numPos + b] =
          kernel[positiveScratch_[a] * numObs + positiveScratch_[b]] * s2;
    // zero-weight members add nothing to z'Wz, so the positive subset's sum
    // is the full member set's
    double logDetNoise = 0.0, responseSumOfSquares = 0.0;
    for (std::size_t a = 0; a < numPos; ++a) {
      std::size_t i = tree.indices[node.begin + positiveScratch_[a]];
      double w = weights[i];
      double noise = residualVariance / w;
      logDetNoise += std::log(noise);
      cholV_[a * numPos + a] += noise;
      responseSumOfSquares += w * y[i] * y[i];
    }
    choleskyDecompose(cholV_.data(), numPos);

    double halfLogDetV = 0.0;
    for (std::size_t a = 0; a < numPos; ++a)
      halfLogDetV += std::log(cholV_[a * numPos + a]);

    vectorScratch_.resize(numPos);
    for (std::size_t a = 0; a < numPos; ++a)
      vectorScratch_[a] = y[tree.indices[node.begin + positiveScratch_[a]]];
    solveLowerTriangular(cholV_.data(), numPos, vectorScratch_.data());
    double quadraticForm = 0.0;
    for (std::size_t a = 0; a < numPos; ++a)
      quadraticForm += vectorScratch_[a] * vectorScratch_[a];

    return -halfLogDetV + 0.5 * logDetNoise -
           0.5 * (quadraticForm - responseSumOfSquares / residualVariance);
  }

  /// Matheron's rule conditioning only on the positive-weight members: f0
  /// covers every member (numObs normals, row order, as always), e0 and
  /// the correction only the positive rows (numPos more normals), so
  /// zero-weight rows draw from the correct conditional law and the
  /// prediction cache stays defined over the full member list. With no
  /// positive members the draw is the prior draw.
  FunctionLeafDrawStats drawFromPosteriorOverPositiveWeights(
    ext_rng* rng, const Tree& tree, const Node& node, const double* y,
    const double* weights, double k, double residualVariance,
    int32_t nodeIndex, double* fits) const {
    std::size_t numObs = node.numObservations();
    positiveScratch_.clear();
    for (std::size_t r = 0; r < numObs; ++r)
      if (weights[tree.indices[node.begin + r]] > 0.0)
        positiveScratch_.push_back(r);
    std::size_t numPos = positiveScratch_.size();

    double s = scale / k, s2 = s * s;

    const double* kernel;
    const double* cholK;
    kernelAndFactorForNode(tree, node, nodeIndex, numObs, &kernel, &cholK);

    epsScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      epsScratch_[r] = ext_rng_simulateStandardNormal(rng);
    fScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a <= r; ++a)
        value += cholK[r * numObs + a] * epsScratch_[a];
      fScratch_[r] = s * value;
    }

    if (numPos > 0) {
      cholV_.resize(numPos * numPos);
      for (std::size_t a = 0; a < numPos; ++a)
        for (std::size_t b = 0; b < numPos; ++b)
          cholV_[a * numPos + b] =
            s2 * kernel[positiveScratch_[a] * numObs + positiveScratch_[b]];
      for (std::size_t a = 0; a < numPos; ++a) {
        double w = weights[tree.indices[node.begin + positiveScratch_[a]]];
        cholV_[a * numPos + a] += residualVariance / w;
      }
      choleskyDecompose(cholV_.data(), numPos);

      vectorScratch_.resize(numPos);
      for (std::size_t a = 0; a < numPos; ++a) {
        std::size_t r = positiveScratch_[a];
        double w = weights[tree.indices[node.begin + r]];
        vectorScratch_[a] = y[tree.indices[node.begin + r]] - fScratch_[r] -
                            std::sqrt(residualVariance / w) *
                              ext_rng_simulateStandardNormal(rng);
      }
      solveLowerTriangular(cholV_.data(), numPos, vectorScratch_.data());
      solveLowerTriangularTransposed(cholV_.data(), numPos,
                               vectorScratch_.data());
      for (std::size_t r = 0; r < numObs; ++r) {
        double value = 0.0;
        for (std::size_t a = 0; a < numPos; ++a)
          value += kernel[r * numObs + positiveScratch_[a]] *
                   vectorScratch_[a];
        fScratch_[r] += s2 * value;
      }
    }
    for (std::size_t r = 0; r < numObs; ++r)
      fits[tree.indices[node.begin + r]] = fScratch_[r];

    return cacheAlphaForNode(nodeIndex, numObs, cholK);
  }

  /// Drop entries whose node no longer exists, is no longer a bottom node,
  /// or whose membership changed, keeping the budget for live leaves.
  /// Lookups re-validate regardless; this is hygiene, not correctness.
  void evictStaleKernelEntries(const Tree& tree) const {
    for (TreeKernelCache& cache : kernelCaches_) {
      if (cache.tree != &tree) continue;
      for (std::size_t index = 0; index < cache.nodes.size(); ++index) {
        CachedLeafKernel& entry = cache.nodes[index];
        if (entry.members.empty()) continue;
        bool live = index < tree.nodes.size();
        if (live) {
          const Node& node(tree.at(static_cast<int32_t>(index)));
          live = node.isBottom() &&
                 node.numObservations() == entry.members.size() &&
                 std::memcmp(entry.members.data(),
                             tree.indices + node.begin,
                             entry.members.size() * sizeof(index_t)) == 0;
        }
        if (!live) {
          kernelCacheUsedBytes_ -= kernelEntryBytes(entry.members.size());
          entry = CachedLeafKernel();
        }
      }
      break;
    }
  }

  std::size_t numCovariates_ = 0;
  std::size_t numObservations_ = 0;
  std::size_t numTestObservations_ = 0;
  std::size_t maxLeafSize_ = 256;
  std::vector<std::size_t> columns_;
  std::vector<double> means_, sds_;
  std::vector<double> lengthscales_, suppliedLengthscales_;
  std::vector<double> u_;      // standardized, column-major n x q
  std::vector<double> uTest_;  // standardized, column-major numTest x q
  // per-call scratch and the per-tree-draw prediction cache; mutable because
  // scoring and drawing are logically const, and safe because a leaf model
  // instance belongs to a single chain (one instance per chain, never shared
  // across threads)
  mutable std::vector<double> leafU_;    // row-major numObs x q, / theta
  mutable std::vector<double> kernel_, cholK_, cholV_;
  mutable std::vector<double> epsScratch_, fScratch_, vectorScratch_;
  mutable std::vector<std::size_t> positiveScratch_;  // w > 0 member offsets
  mutable std::vector<double> alphaBuffer_;
  mutable std::vector<double> nodeConstant_;         // arena-indexed
  mutable std::vector<std::ptrdiff_t> nodeAlphaOffset_;  // -1 = constant
  // the cross-sweep kernel cache, keyed by tree address (stable for a
  // chain's lifetime) and node arena index
  mutable std::vector<TreeKernelCache> kernelCaches_;
  mutable std::size_t kernelCacheUsedBytes_ = 0;
  std::size_t kernelCacheBudget_ = kernelCacheTotalBudgetBytes;
};

static_assert(FunctionLeafModel<GPGaussianLeaf>);

/// Chipman-George-McCulloch tree structure prior. Split-variable selection
/// is uniform over available variables when splitProbabilities is null;
/// otherwise proportional to the supplied weights restricted to available
/// variables. DART points this at its Gibbs-updated probability vector.
struct CGMTreePrior {
  double base = 0.95;
  double power = 2.0;
  const double* splitProbabilities = nullptr;  // length numPredictors, or null
  // pooled-column scratch (reachable set, pattern draw); mutable because
  // rule scoring and drawing are logically const, and safe because a prior
  // instance belongs to a single chain
  mutable std::vector<std::uint64_t> reachableScratch_;
  mutable std::vector<std::uint64_t> patternScratch_;
  // per-variable availability from Tree::collectAvailableVariables
  mutable std::vector<std::uint8_t> availableScratch_;

  double growthProbability(const Tree& tree, const ColumnStore& data,
                           int32_t nodeIndex) const {
    if (!tree.hasAnyAvailableVariable(data, nodeIndex)) return 0.0;
    return base / std::pow(1.0 + static_cast<double>(tree.depthOf(nodeIndex)), power);
  }

  double splitVariableLogProbability(const Tree& tree, const ColumnStore& data,
                                     int32_t nodeIndex) const {
    availableScratch_.resize(data.numPredictors);
    std::size_t numAvailable =
      tree.collectAvailableVariables(data, nodeIndex, availableScratch_.data());
    if (splitProbabilities == nullptr)
      return -std::log(static_cast<double>(numAvailable));

    double totalProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j)
      if (availableScratch_[j]) totalProbability += splitProbabilities[j];
    return std::log(
      splitProbabilities[tree.at(nodeIndex).rule.variableIndex] / totalProbability);
  }

  /// Categorical rules are uniform over the direction assignments of the
  /// categories reachable at the node that leave neither side empty, both
  /// orientations counted: 2^R - 2 of them. Unreachable categories carry no
  /// probability (their direction bits are pinned to zero by the moves).
  double ruleForVariableLogProbability(const Tree& tree, const ColumnStore& data,
                                       int32_t nodeIndex) const {
    int32_t variableIndex = tree.at(nodeIndex).rule.variableIndex;
    if (data.types[static_cast<size_t>(variableIndex)] ==
        ColumnType::categorical) {
      size_t j = static_cast<size_t>(variableIndex);
      size_t numReachable;
      if (data.columnIsPooled(j)) {
        size_t numWords = maskWordsForCount(data.numCuts[j]);
        reachableScratch_.resize(numWords);
        tree.reachableCategoriesWide(data, nodeIndex, variableIndex,
                                     reachableScratch_.data());
        numReachable = maskPopcount(reachableScratch_.data(), numWords);
      } else {
        numReachable = static_cast<size_t>(std::popcount(
          tree.reachableCategories(data, nodeIndex, variableIndex)));
      }
      // past 54 reachable pow(2, R) - 2 is no longer exact; the closed form
      // is a deterministic function of R, so MH ratios stay consistent
      if (numReachable > 54)
        return -(static_cast<double>(numReachable) * std::log(2.0) +
                 std::log1p(-std::pow(2.0, 1.0 -
                                             static_cast<double>(numReachable))));
      return -std::log(std::pow(2.0, static_cast<double>(numReachable)) - 2.0);
    }
    int32_t left, right;
    tree.splitInterval(data, nodeIndex, variableIndex, &left, &right);
    double logNumRules = std::log(static_cast<double>(right - left + 1));
    // a column with missing values widens every rule by its two-way
    // missing direction (the categorical count absorbs the missing
    // category through the reachable mask instead)
    if (data.hasMissing[static_cast<size_t>(variableIndex)])
      logNumRules += std::log(2.0);
    return -logNumRules;
  }

  /// The prior mass ONE grow-from-root categorical candidate carries, the
  /// commensurable sibling of ruleForVariableLogProbability above; the two
  /// derivations sit together deliberately.
  ///
  /// That builder enumerates the partitions of the P categories PRESENT at a
  /// node, not the 2^R - 2 masks. Grouping the masks by the partition they
  /// induce on the present set leaves 2^(A+1) masks per unordered partition
  /// (A = R - P absent positions, both orientations) over 2^(P-1) - 1
  /// partitions, the residual 2^(A+1) - 2 masks being exactly the empty-child
  /// rules, so the variable's whole enumerable family carries
  ///
  ///     M(R, P) = (2^(P-1) - 1) 2^(A+1) / (2^R - 2)
  ///             = (1 - 2^(1-P)) / (1 - 2^(1-R))
  ///
  /// and a candidate carries that spread evenly over the numEmitted the
  /// enumeration emits. Below the cap numEmitted is 2^(P-1) - 1 and this
  /// reduces (to within an ulp) to the per-partition group mass
  /// (1-P) log 2 - log1p(-exp2(1-R)); above it the sorted-prefix family's P - 1
  /// candidates carry the SAME family total, which is what keeps a variable's
  /// realized split mass continuous in its level count rather than dropping
  /// ~100x at a compile-time constant the user cannot see. The retained
  /// prefixes therefore inherit the whole family's mass, which over-weights a
  /// high-cardinality factor by at most (2^(P-1)-1)/(P-1) in the strong-signal
  /// regime, where a likelihood ratio of exp(O(n_node)) decides anyway; the
  /// exact conditional errs by the same factor the other way in the regime
  /// where the prior IS what decides.
  ///
  /// The log1p/exp2 spelling needs no R > 54 branch: exp2(1-R) underflows to
  /// zero harmlessly where pow(2, R) - 2 stops being exact and then overflows.
  static double categoricalGroupLogProbability(size_t numReachable,
                                               size_t numPresent,
                                               double numEmitted) {
    // meaningless outside this range: it returns a "probability" above one at
    // P < 2, and P > R contradicts every present category being reachable
    assert(numPresent >= 2 && numPresent <= numReachable);
    double R = static_cast<double>(numReachable);
    double P = static_cast<double>(numPresent);
    return std::log1p(-std::exp2(1.0 - P)) - std::log1p(-std::exp2(1.0 - R)) -
           std::log(numEmitted);
  }

  double treeLogProbability(const Tree& tree, const ColumnStore& data,
                            int32_t nodeIndex = 0) const {
    double growth = growthProbability(tree, data, nodeIndex);
    if (tree.at(nodeIndex).isBottom()) return std::log(1.0 - growth);

    double result = std::log(growth);
    result += splitVariableLogProbability(tree, data, nodeIndex);
    result += ruleForVariableLogProbability(tree, data, nodeIndex);
    result += treeLogProbability(tree, data, tree.at(nodeIndex).leftChild);
    result += treeLogProbability(tree, data, tree.at(nodeIndex).leftChild + 1);
    return result;
  }

  int32_t drawSplitVariable(const Tree& tree, const ColumnStore& data,
                            ext_rng* rng, int32_t nodeIndex) const {
    availableScratch_.resize(data.numPredictors);
    std::size_t numGood =
      tree.collectAvailableVariables(data, nodeIndex, availableScratch_.data());

    if (splitProbabilities == nullptr) {
      std::size_t variableNumber =
        ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, numGood);
      std::size_t count = 0;
      for (std::size_t j = 0; j < data.numPredictors; ++j)
        if (availableScratch_[j]) {
          if (count == variableNumber) return static_cast<int32_t>(j);
          ++count;
        }
      return invalidVariable;
    }

    double totalProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j)
      if (availableScratch_[j]) totalProbability += splitProbabilities[j];

    double cutoff = ext_rng_simulateContinuousUniform(rng) * totalProbability;
    double runningProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j) {
      if (availableScratch_[j]) {
        runningProbability += splitProbabilities[j];
        if (runningProbability >= cutoff) return static_cast<int32_t>(j);
      }
    }
    return invalidVariable;  // unreachable with valid probabilities
  }

  /// Build a categorical direction mask from a pattern whose bit k gives
  /// the side of the k-th reachable category (ascending); unreachable
  /// categories stay zero (left), the canonical gauge.
  static std::uint64_t categoryDirectionsForPattern(std::uint64_t reachable,
                                                    std::uint64_t pattern) {
    std::uint64_t directions = 0;
    int bit = 0;
    while (reachable != 0) {
      std::uint32_t category = static_cast<std::uint32_t>(
        std::countr_zero(reachable));
      if ((pattern >> bit) & 1u) directions |= 1ull << category;
      reachable &= reachable - 1;
      ++bit;
    }
    return directions;
  }

  /// Uniform over the direction patterns of numReachable categories that
  /// leave neither side empty, drawn bit by bit with the two all-same
  /// patterns rejected: a single range draw of u * (2^R - 2) has only the
  /// generator's granularity, which for wide masks pins the low pattern
  /// bits to functions of the high ones.
  static std::uint64_t drawCategoryPattern(ext_rng* rng, int numReachable) {
    // 63 categories plus the missing pseudo-category can fill the word
    std::uint64_t allRight =
      numReachable >= 64 ? ~0ull : (1ull << numReachable) - 1ull;
    std::uint64_t pattern;
    do {
      pattern = 0;
      for (int bit = 0; bit < numReachable; ++bit)
        if (ext_rng_simulateBernoulli(rng, 0.5) == 1) pattern |= 1ull << bit;
    } while (pattern == 0 || pattern == allRight);
    return pattern;
  }

  /// The pooled-column extension of drawCategoryPattern: the identical
  /// bit-by-bit scheme, one draw per reachable position ascending, the two
  /// all-same patterns rejected.
  static void drawCategoryPatternWide(ext_rng* rng, size_t numReachable,
                                      std::uint64_t* pattern,
                                      size_t numWords) {
    bool degenerate = true;
    while (degenerate) {
      for (size_t w = 0; w < numWords; ++w) pattern[w] = 0;
      for (size_t bit = 0; bit < numReachable; ++bit)
        if (ext_rng_simulateBernoulli(rng, 0.5) == 1)
          pattern[bit >> 6] |= 1ull << (bit & 63u);
      degenerate = maskIsZero(pattern, numWords) ||
                   maskPopcount(pattern, numWords) == numReachable;
    }
  }

  /// The pooled-column extension of categoryDirectionsForPattern.
  static void categoryDirectionsForPatternWide(const std::uint64_t* reachable,
                                               const std::uint64_t* pattern,
                                               std::uint64_t* directions,
                                               size_t numWords) {
    for (size_t w = 0; w < numWords; ++w) directions[w] = 0;
    size_t bit = 0;
    for (size_t w = 0; w < numWords; ++w) {
      std::uint64_t remaining = reachable[w];
      while (remaining != 0) {
        std::uint32_t category =
          static_cast<std::uint32_t>(std::countr_zero(remaining)) +
          64u * static_cast<std::uint32_t>(w);
        if (((pattern[bit >> 6] >> (bit & 63u)) & 1u) != 0)
          maskSetBit(directions, category);
        remaining &= remaining - 1;
        ++bit;
      }
    }
  }

  /// The tree is non-const because a pooled categorical draw allocates its
  /// mask in the tree's pool; the caller truncates back on rejection.
  Rule drawRuleForVariable(Tree& tree, const ColumnStore& data,
                           ext_rng* rng, int32_t nodeIndex,
                           int32_t variableIndex) const {
    Rule result;
    result.variableIndex = variableIndex;

    if (data.types[static_cast<size_t>(variableIndex)] ==
        ColumnType::categorical) {
      size_t j = static_cast<size_t>(variableIndex);
      if (data.columnIsPooled(j)) {
        size_t numWords = maskWordsForCount(data.numCuts[j]);
        reachableScratch_.resize(numWords);
        tree.reachableCategoriesWide(data, nodeIndex, variableIndex,
                                     reachableScratch_.data());
        size_t numReachable =
          maskPopcount(reachableScratch_.data(), numWords);
        patternScratch_.resize(numWords);
        drawCategoryPatternWide(rng, numReachable, patternScratch_.data(),
                                numWords);
        size_t offset = tree.allocateMask(numWords);
        categoryDirectionsForPatternWide(reachableScratch_.data(),
                                         patternScratch_.data(),
                                         tree.mutableMaskWordsFor(offset),
                                         numWords);
        result.setMaskOffset(offset);
        return result;
      }
      std::uint64_t reachable =
        tree.reachableCategories(data, nodeIndex, variableIndex);
      int numReachable = std::popcount(reachable);
      std::uint64_t pattern = drawCategoryPattern(rng, numReachable);
      result.setCategoryDirections(
        categoryDirectionsForPattern(reachable, pattern));
      return result;
    }

    int32_t left, right;
    tree.splitInterval(data, nodeIndex, variableIndex, &left, &right);
    result.setSplitIndex(static_cast<int32_t>(
      ext_rng_simulateIntegerUniformInRange(rng, left, right + 1)));
    // the missing direction is part of the rule, a symmetric coin drawn
    // only when the column can route a missing value - NA-free data spends
    // no extra generator draws
    if (data.hasMissing[static_cast<size_t>(variableIndex)])
      result.setMissingGoesRight(ext_rng_simulateBernoulli(rng, 0.5) == 1);
    return result;
  }

  Rule drawRuleAndVariable(Tree& tree, const ColumnStore& data,
                           ext_rng* rng, int32_t nodeIndex) const {
    int32_t variableIndex = drawSplitVariable(tree, data, rng, nodeIndex);
    return drawRuleForVariable(tree, data, rng, nodeIndex, variableIndex);
  }
};

/// Stabilized-softmax draw over a grid of log-weights: normalize by the log-
/// sum-exp (subtract the max, exponentiate, sum, divide), then draw a discrete
/// index. `logWeights` is overwritten in place with the normalized
/// probabilities. Shared by the grid full conditionals (DartPrior's alpha,
/// ResidualDfPrior's nu, NBDispersionPrior's r), each of which fills the array
/// with its own log-posterior first.
inline std::size_t drawFromLogWeights(ext_rng* rng, double* logWeights,
                                      std::size_t n) {
  double maxLog = -HUGE_VAL;
  for (std::size_t k = 0; k < n; ++k)
    if (logWeights[k] > maxLog) maxLog = logWeights[k];
  double total = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    logWeights[k] = std::exp(logWeights[k] - maxLog);
    total += logWeights[k];
  }
  for (std::size_t k = 0; k < n; ++k) logWeights[k] /= total;
  return ext_rng_drawFromDiscreteDistribution(rng, logWeights, n);
}

/// DART (Linero 2018): Dirichlet prior over split-variable probabilities.
/// s | counts ~ Dirichlet(alpha/p + m_1, ..., alpha/p + m_p), sampled via
/// normalized gammas. The concentration alpha is optionally sampled on a
/// fixed grid of lambda = alpha / (alpha + rho) with lambda ~ Beta(a, b);
/// all lgamma terms of the grid are precomputed, so the per-iteration alpha
/// update is O(gridSize) multiply-adds.
struct DartPrior {
  double alpha = 1.0;
  bool updateAlpha = true;
  double betaA = 0.5, betaB = 1.0, rho = 0.0;  // rho <= 0 means numPredictors
  // Iterations to hold s uniform before updates begin, so the forest is
  // likelihood-informed when counts first enter the Dirichlet; starting cold
  // can lock s away from the signal variables. The BART package convention
  // is half of burn-in.
  std::size_t updateDelay = 0;
  std::vector<double> probabilities;

  void initialize(std::size_t numPredictors) {
    probabilities.assign(numPredictors, 1.0 / static_cast<double>(numPredictors));
    if (rho <= 0.0) rho = static_cast<double>(numPredictors);
    if (updateAlpha) {
      gridAlpha_.resize(gridSize);
      gridConstant_.resize(gridSize);
      gridWeight_.resize(gridSize);
      double p = static_cast<double>(numPredictors);
      for (std::size_t k = 0; k < gridSize; ++k) {
        double lambda = static_cast<double>(k + 1) / static_cast<double>(gridSize + 1);
        double alphaK = lambda * rho / (1.0 - lambda);
        gridAlpha_[k] = alphaK;
        gridConstant_[k] = std::lgamma(alphaK) - p * std::lgamma(alphaK / p) +
                           (betaA - 1.0) * std::log(lambda) +
                           (betaB - 1.0) * std::log(1.0 - lambda);
      }
    }
  }

  void update(ext_rng* rng, const std::uint32_t* splitCounts) {
    if (numUpdatesSkipped_ < updateDelay) {
      ++numUpdatesSkipped_;
      return;
    }

    std::size_t numPredictors = probabilities.size();
    double p = static_cast<double>(numPredictors);

    double total = 0.0;
    for (std::size_t j = 0; j < numPredictors; ++j) {
      double draw = ext_rng_simulateGamma(
        rng, alpha / p + static_cast<double>(splitCounts[j]), 1.0);
      probabilities[j] = draw > 1.0e-300 ? draw : 1.0e-300;
      total += probabilities[j];
    }
    for (std::size_t j = 0; j < numPredictors; ++j) probabilities[j] /= total;

    if (!updateAlpha) return;

    double sumLogProbabilities = 0.0;
    for (std::size_t j = 0; j < numPredictors; ++j)
      sumLogProbabilities += std::log(probabilities[j]);

    for (std::size_t k = 0; k < gridSize; ++k)
      gridWeight_[k] = gridConstant_[k] + (gridAlpha_[k] / p) * sumLogProbabilities;
    alpha = gridAlpha_[drawFromLogWeights(rng, gridWeight_.data(), gridSize)];
  }

  static constexpr std::size_t gridSize = 1000;

  // State serialization: the delay counter is the only hidden mutable state
  // (alpha and probabilities are public members).
  std::size_t numUpdatesSkipped() const { return numUpdatesSkipped_; }
  void setNumUpdatesSkipped(std::size_t value) { numUpdatesSkipped_ = value; }

private:
  std::vector<double> gridAlpha_, gridConstant_, gridWeight_;
  std::size_t numUpdatesSkipped_ = 0;
};

/// Chi hyperprior on the end-node precision parameter k. Since k ~ chi(nu)
/// gives k^2 ~ Gamma(nu/2, 1/2), the posterior of k^2 is gamma with shape
/// 0.5 (M + nu); a finite prior scale adds 0.5 / scale^2 to the rate.
struct ChiKHyperprior {
  double degreesOfFreedom = 1.5;
  double scale = HUGE_VAL;  // infinite = flat in the rate term

  /// Sentinel ceiling on the sampled k. An improper or weak prior scale leaves
  /// the prior-dominated k Gibbs fixed point with a growth factor above 1, so
  /// k can run away toward a leaf sd that is already statistically zero on the
  /// standardized [-0.5, 0.5] response scale; capping the draw bounds the
  /// excursion. Behavior-neutral outside that runaway regime: healthy draws
  /// sit far below it, so the cap never engages.
  static constexpr double maxDraw = 1.0e6;

  double draw(ext_rng* rng, double sumSquaredParams, double totalNumLeaves,
              double leafScale) const {
    double shape = 0.5 * (totalNumLeaves + degreesOfFreedom);
    // classic form: numTrees * s_sq / nodeScale^2 == s_sq / leafScale^2
    double rate = 0.5 * sumSquaredParams / (leafScale * leafScale);
    if (std::fabs(scale) <= DBL_MAX) rate += 0.5 / (scale * scale);
    double k = std::sqrt(ext_rng_simulateGamma(rng, shape, 1.0 / rate));
    return k > maxDraw ? maxDraw : k;
  }
};

/// Conjugate chi-squared residual variance prior; scale arrives already
/// derived (qchisq(1 - quantile, df) / df, then multiplied by the initial
/// sigma^2 estimate at initialization).
struct ChiSquaredScalePrior {
  double degreesOfFreedom = 3.0;
  double scale = 1.0;

  /// numPositiveWeights is the count of rows with weight > 0 (== numObservations
  /// when unweighted): zero-weight rows drop from the weighted SSR and every
  /// leaf conditional, so they must not inflate the posterior df either.
  double drawSigmaSqFromPosterior(ext_rng* rng, const double* y,
                                  const double* totalFits, const double* weights,
                                  std::size_t numObservations,
                                  std::size_t numPositiveWeights) const {
    double sumOfSquaredResiduals = weights == nullptr
      ? misc_computeSumOfSquaredResiduals(y, numObservations, totalFits)
      : misc_computeWeightedSumOfSquaredResiduals(y, numObservations, weights,
                                                  totalFits);
    double posteriorDegreesOfFreedom =
      degreesOfFreedom + static_cast<double>(numPositiveWeights);
    double posteriorScale = degreesOfFreedom * scale + sumOfSquaredResiduals;
    return posteriorScale /
           ext_rng_simulateChiSquared(rng, posteriorDegreesOfFreedom);
  }
};

/// The leaf model in force, as a value. All four location-scale leaves carry
/// the one `scale` field the calibration surface writes, so the write is total
/// over them; what the reported prior sd MEANS is leaf-specific and only this
/// tag can say which (docs/design/nameable-calibration.md section 3).
enum class LeafModelKind { constant = 0, monotone = 1, linear = 2, gp = 3 };

/// A leaf model type's tag. A property of the type alone, so constexpr; the
/// numbering is the one dbarts.h's dbarts_leaf_model publishes.
template <typename L> constexpr LeafModelKind leafModelKindOf() {
  if constexpr (std::is_same_v<L, MonotoneConstantGaussianLeaf>)
    return LeafModelKind::monotone;
  else if constexpr (std::is_same_v<L, LinearGaussianLeaf>)
    return LeafModelKind::linear;
  else if constexpr (std::is_same_v<L, GPGaussianLeaf>)
    return LeafModelKind::gp;
  else
    return LeafModelKind::constant;
}

/// Response families the sampler can run; gaussian fits the response
/// directly, the binary families fit a latent working response, aft fits
/// log survival times (right-censored observations carry latent log-times),
/// nbinom fits non-negative counts via the Polya-Gamma augmentation (NBResponse,
/// docs/design/negative-binomial.md).
enum class ResponseFamily { gaussian, probit, logistic, aft, ordinal, nbinom };

/// Numerically stable log(1 + exp(x)): the logistic log-likelihood's building
/// block, guarding against overflow for large x.
inline double logOnePlusExp(double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

/// A Polya-Gamma draw PG(b, z) for shape b: the omega-loop's shape-parameterized
/// seam (docs/design/negative-binomial.md section 2). v1 ships the exact integer
/// envelope, so b is integer-valued and PG(b, z) is the exact sum of round(b)
/// independent PG(1, z) variates (the LogisticResponse reduction, PG(n, z) =
/// sum of n PG(1, z)). At b = 1 this is bit-for-bit the shipped PG(1) Devroye
/// stream. A later real-shape mode routes a fractional b to an approximate
/// primitive here without touching the callers.
inline double simulatePolyaGammaShape(ext_rng* rng, double b, double z) {
  long reps = std::lround(b);
  double omega = 0.0;
  for (long c = 0; c < reps; ++c)
    omega += ext_rng_simulatePolyaGamma(rng, z);
  return omega;
}

/// Reshift a Polya-Gamma latent family's working response when only the offset
/// changes: the omega draws and the kappa / omega term stand, so each entry
/// trades the old offset for the new. A null pointer is a zero offset. Shared
/// by LogisticResponse::setOffset and NBResponse::setOffset.
inline void reshiftWorkingForOffset(double* working, const double* oldOffset,
                                    const double* newOffset, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    double unshifted = working[i] + (oldOffset != nullptr ? oldOffset[i] : 0.0);
    working[i] = unshifted - (newOffset != nullptr ? newOffset[i] : 0.0);
  }
}

/// Response models own the working response the backfitting engine sees and
/// any per-iteration latent refresh; concrete classes own their O(n) loops.
/// Response/offset pointers are borrowed: the caller keeps them alive.
class ResponseModel {
public:
  virtual ~ResponseModel() = default;

  /// The engine's working response; contents may change across iterations
  /// for latent-variable models.
  virtual double* workingResponse() = 0;

  /// Per-observation precisions the engine weights tree updates by: the
  /// fixed user weights for gaussian, the current Polya-Gamma draws for
  /// logistic, null (unit weights) for probit. Refreshed alongside the
  /// working response.
  virtual const double* workingWeights() const = 0;

  /// Whether workingWeights() changes across iterations (the latent
  /// Polya-Gamma refresh). Cross-sweep sufficient-statistic caches drop on
  /// each refresh when true; false families weight by fixed user values.
  virtual bool workingWeightsVaryPerSweep() const { return false; }

  /// Called once per iteration after all trees update; sigma is the chain's
  /// current residual sd on the internal scale (1 for the binary families),
  /// which only the grouped decorator's conjugate update consumes.
  virtual void refreshLatents(ext_rng* rng, const double* totalFits,
                              double sigma) = 0;

  /// Residual standard deviation update; fixed-sigma models return sigma.
  virtual double drawSigma(ext_rng* rng, const double* totalFits,
                           double sigma) = 0;

  // Between-sample mutation (the embedded-Gibbs API). sigmaInOut is the
  // chain's current sigma on the internal scale; models that rescale keep it
  // fixed on the original scale across the change. updateScale re-anchors the
  // transform to the new response, as setOffset does; false locks it.
  virtual void setResponse(const double* y, ext_rng* rng,
                           const double* totalFits, bool updateScale,
                           double* sigmaInOut) = 0;
  virtual void setOffset(const double* offset, bool updateScale,
                         double* sigmaInOut) = 0;

  /// Whole-data replacement, possibly changing the number of observations;
  /// latent states are cold-initialized since the fits they condition on are
  /// stale. Models that rescale keep sigmaInOut fixed on the original scale.
  virtual void setData(const double* y, const double* offset,
                       const double* weights, std::size_t numObservations,
                       double* sigmaInOut) = 0;

  /// Replace the case weights (borrowed). For a family whose weights are a
  /// precision this is a bare pointer swap: nothing rescales, and the weighted
  /// residuals enter the next iteration's node statistics and sigma draw. A
  /// family whose augmentation is STATED against the weights - logistic, whose
  /// counts are the Polya-Gamma shape - redraws its latents here instead of
  /// leaving them for the next sweep, which moves its trees first; that is why
  /// the chain generator and the location to draw against arrive too, on the
  /// shape every other mutation setter already has. Families that need neither
  /// take both and drop them. Default: a no-op (the host rejects earlier).
  virtual void setWeights(const double*, ext_rng*, const double*) {}

  /// Whether this family implements the active-row channel. setActiveRows's
  /// own refusal reads it, so the advertised capability and the refusal
  /// cannot disagree.
  virtual bool supportsActiveRows() const { return false; }

  /// Install (or clear, at a null pointer) a per-observation 0/1 active-row
  /// mask: a_i = 0 says row i is not in the data set for this sampler this
  /// sweep. The family composes a into its own workingWeights() - keeping the
  /// user weight pointer, which the two setters hold independently - and skips
  /// the inactive rows' latent draws and family-level sums. The values are
  /// COPIED, so the caller's array is free the moment this returns. The host
  /// has already validated the values and normalized an all-ones mask to a
  /// null pointer; false is the refusal of a family that implements none.
  virtual bool setActiveRows(const double*) { return false; }

  /// Replace the residual-variance prior: re-anchors to the supplied
  /// original-scale sigma estimate exactly as construction does, so a swap
  /// before any run matches creating with the new prior. A no-op for the
  /// fixed-sigma binary families.
  virtual void setSigmaPrior(double /*sigmaEstimate*/, double /*degreesOfFreedom*/,
                             double /*rawScale*/) {}

  virtual const double* latents() const { return nullptr; }

  /// The current training offset (borrowed), or null. Recorded training
  /// fits add it back, keeping them on the original scale.
  virtual const double* offset() const = 0;

  /// State restoration: overwrite the latent state with previously stored
  /// latents() values and rebuild the working response from them. A no-op
  /// for models without latents.
  virtual void restoreLatents(const double*) {}

  /// State serialization of the response transform: gaussian's offset-
  /// adjusted (min, max), whose evolution under setOffset(updateScale) is
  /// otherwise unrecoverable from the data; scale-free families report
  /// (0, 0) and ignore restoration.
  virtual void getScale(double& min, double& max) const { min = max = 0.0; }
  virtual void restoreScale(double /*min*/, double /*max*/) {}

  virtual double initialSigma() const = 0;

  // Sample de-scaling to the original response scale.
  virtual double fitScale() const = 0;
  virtual double fitShift() const = 0;
  virtual double sigmaScale() const = 0;

  /// Test hook: the residual-variance posterior's degrees of freedom,
  /// nu_0 + #{w_i > 0} over THIS model's own precisions. Zero for a family
  /// that draws no sigma.
  virtual double sigmaDegreesOfFreedomForTesting() const { return 0.0; }

  /// Per-observation log-likelihood of the training response under the
  /// current fit, into out (numObservations values). totalFits is the forest
  /// sum on the internal scale and sigma the internal residual sd; each family
  /// applies its own transform and density. The default NaN-fills, so a family
  /// without a defined channel reports "unavailable" rather than a wrong value.
  virtual void computeLogLikelihood(const double* /*totalFits*/,
                                    double /*sigma*/,
                                    std::size_t numObservations,
                                    double* out) const {
    for (std::size_t i = 0; i < numObservations; ++i)
      out[i] = std::numeric_limits<double>::quiet_NaN();
  }

  /// Grouped random effects, when the model carries them (the GroupedResponse
  /// decorator): the per-group intercepts and tau, on the internal scale.
  /// Ungrouped models report none.
  virtual const double* groupEffects() const { return nullptr; }
  virtual std::size_t numGroupEffects() const { return 0; }
  virtual double groupTau() const { return 0.0; }
  /// State restoration of the grouped channel; values are internal-scale
  /// exactly as reported, so restores are exact. A no-op when ungrouped.
  virtual void restoreGroupEffects(const double* /*effects*/,
                                   double /*tau*/) {}

  /// Continuous responses under a Student-t error law (TResponse) carry a
  /// residual degrees of freedom nu the state block serializes alongside the
  /// mixing precisions in latents(); other families carry none, so their
  /// states omit the nu block and a t sampler refuses a state lacking one.
  virtual bool carriesResidualDf() const { return false; }
  virtual double residualDf() const { return 0.0; }
  virtual void restoreResidualDf(double /*nu*/) {}

  /// Ordinal (cumulative-probit) responses (OrdinalResponse) carry a length-
  /// (K-1) cutpoint vector the state block serializes as a by-name "cutpoints"
  /// slot; other families carry none, so their states omit it and an ordinal
  /// sampler refuses a state lacking one. numCutpoints() is the block length,
  /// the residualDf trio's vector analog (a scalar needed no length).
  virtual bool carriesCutpoints() const { return false; }
  virtual const double* cutpoints() const { return nullptr; }
  virtual std::size_t numCutpoints() const { return 0; }
  virtual void restoreCutpoints(const double* /*gamma*/) {}

  /// Count responses under a negative-binomial law (NBResponse) carry a scalar
  /// dispersion r the state block serializes as a by-name "dispersion" slot;
  /// other families carry none, so their states omit it and an NB sampler
  /// refuses a state lacking one (the residualDf scalar analog). The name is
  /// parameterization-neutral and the value real-valued (grid mode stores an
  /// integer-valued double), so a later real-r mode loads and saves with no
  /// state-format change. RESTORE CONTRACT: restoreDispersion MUST run before
  /// restoreLatents, which rebuilds the working response from omega AND r.
  virtual bool carriesDispersion() const { return false; }
  virtual double dispersion() const { return 0.0; }
  virtual void restoreDispersion(double /*r*/) {}
};

class GaussianResponse final : public ResponseModel {
public:
  /// sigmaRawScale is qchisq(1 - quantile, df) / df; offset may be null.
  GaussianResponse(const double* y, const double* offset, const double* weights,
                   std::size_t numObservations, double sigmaEstimate,
                   double sigmaDf, double sigmaRawScale)
    : y_(y), offset_(offset), weights_(weights), userWeights_(weights),
      numObservations_(numObservations),
      numPositiveWeights_(countPositiveWeights(weights, numObservations)) {
    yRescaled_.resize(numObservations);
    rescale();

    initialSigma_ = sigmaEstimate / range_;
    sigmaSqPrior_.degreesOfFreedom = sigmaDf;
    sigmaSqPrior_.scale = initialSigma_ * initialSigma_ * sigmaRawScale;
  }

  double* workingResponse() override { return yRescaled_.data(); }
  const double* workingWeights() const override { return weights_; }
  const double* offset() const override { return offset_; }
  void refreshLatents(ext_rng*, const double*, double) override {}

  double drawSigma(ext_rng* rng, const double* totalFits, double) override {
    return std::sqrt(sigmaSqPrior_.drawSigmaSqFromPosterior(
      rng, yRescaled_.data(), totalFits, weights_, numObservations_,
      numPositiveWeights_));
  }

  double sigmaDegreesOfFreedomForTesting() const override {
    return sigmaSqPrior_.degreesOfFreedom +
           static_cast<double>(numPositiveWeights_);
  }

  void setResponse(const double* y, ext_rng*, const double*, bool updateScale,
                   double* sigmaInOut) override {
    if (updateScale) {
      // sigma and the variance prior stay fixed on the original scale
      double sigmaUnscaled = *sigmaInOut * range_;
      double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

      y_ = y;
      rescale();

      sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
      *sigmaInOut = sigmaUnscaled / range_;
    } else {
      // reuse the existing data scale
      y_ = y;
      std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
      if (offset_ != nullptr)
        misc_subtractVectorsInPlace(offset_, numObservations_, yRescaled_.data());
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
      misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                       1.0 / range_);
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);
    }
  }

  void setData(const double* y, const double* offset, const double* weights,
               std::size_t numObservations, double* sigmaInOut) override {
    // sigma and the variance prior stay fixed on the original scale
    double sigmaUnscaled = *sigmaInOut * range_;
    double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

    y_ = y;
    offset_ = offset;
    numObservations_ = numObservations;
    // the mask is length-n and n may have changed; a data swap drops it
    activeRows_.clear();
    composite_.clear();
    userWeights_ = weights;
    installWeights(weights);
    yRescaled_.resize(numObservations);
    rescale();

    sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
    *sigmaInOut = sigmaUnscaled / range_;
  }

  /// The two setters are absolute and independent: this replaces the borrowed
  /// user weights and recomposes against whatever mask is installed, so the
  /// served precisions are w * a in either call order.
  void setWeights(const double* weights, ext_rng*, const double*) override {
    userWeights_ = weights;
    recomposeActiveRows();
  }

  bool supportsActiveRows() const override { return true; }

  /// a_i multiplies the case weight, so a masked gaussian is exactly
  /// setWeights(w * a) - including the sigma posterior's degrees of freedom,
  /// which installWeights RECOUNTS off the composite rather than leaving at
  /// the unmasked total.
  bool setActiveRows(const double* active) override {
    if (active == nullptr) {
      if (activeRows_.empty()) return true;
      activeRows_.clear();
      installWeights(userWeights_);  // the pre-mask pointer, by identity
      return true;
    }
    activeRows_.assign(active, active + numObservations_);
    composite_.resize(numObservations_);
    recomposeActiveRows();
    return true;
  }

  void setSigmaPrior(double sigmaEstimate, double degreesOfFreedom,
                     double rawScale) override {
    double sigmaInternal = sigmaEstimate / range_;
    sigmaSqPrior_.degreesOfFreedom = degreesOfFreedom;
    sigmaSqPrior_.scale = sigmaInternal * sigmaInternal * rawScale;
  }

  void setOffset(const double* offset, bool updateScale,
                 double* sigmaInOut) override {
    if (updateScale) {
      double sigmaUnscaled = *sigmaInOut * range_;
      double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

      offset_ = offset;
      rescale();

      sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
      *sigmaInOut = sigmaUnscaled / range_;
    } else {
      // reuse the existing data scale
      offset_ = offset;
      std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
      if (offset_ != nullptr)
        misc_subtractVectorsInPlace(offset_, numObservations_, yRescaled_.data());
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
      misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                       1.0 / range_);
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);
    }
  }

  double initialSigma() const override { return initialSigma_; }
  double fitScale() const override { return range_; }
  double fitShift() const override { return range_ * 0.5 + min_; }
  double sigmaScale() const override { return range_; }

  /// dnorm(y_i, mu_i, sigma_i) with mu the original-scale fit (fits carry any
  /// offset) and sigma_i the residual sd scaled by the precision weight:
  /// y | x ~ N(f(x) + offset, sigma^2 / w_i). A row whose composed weight is
  /// zero is not in the model at all, so its entry is NaN - the channel's own
  /// "unavailable" flag - rather than the -Inf an infinite sd would give.
  void computeLogLikelihood(const double* totalFits, double sigma,
                            std::size_t numObservations,
                            double* out) const override {
    double shift = range_ * 0.5 + min_;
    double sigmaOriginal = sigma * range_;
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (weights_ != nullptr && weights_[i] == 0.0) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      double mu = range_ * totalFits[i] + shift +
                  (offset_ != nullptr ? offset_[i] : 0.0);
      double sd = weights_ != nullptr ? sigmaOriginal / std::sqrt(weights_[i])
                                      : sigmaOriginal;
      out[i] = Rf_dnorm4(y_[i], mu, sd, 1);
    }
  }

  void getScale(double& min, double& max) const override {
    min = min_;
    max = max_;
  }

  /// Installs a stored transform, rebuilding the working response under it
  /// and re-anchoring the variance prior on the original scale, exactly as
  /// the mutation methods do; the chain's internal sigma is restored from
  /// the same state afterwards, so no sigmaInOut here.
  void restoreScale(double min, double max) override {
    double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

    min_ = min;
    max_ = max;
    range_ = max_ - min_;
    if (range_ == 0.0) range_ = 1.0;

    std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
    if (offset_ != nullptr)
      misc_subtractVectorsInPlace(offset_, numObservations_,
                                  yRescaled_.data());
    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
    misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                     1.0 / range_);
    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);

    sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
  }

private:
  /// The one write path to the served precisions: every install recounts the
  /// positive ones, so the sigma posterior's df can never lag the weights.
  void installWeights(const double* weights) {
    weights_ = weights;
    numPositiveWeights_ = countPositiveWeights(weights, numObservations_);
  }

  /// c_i = w_i a_i, or the borrowed user pointer BY IDENTITY when no mask is
  /// installed - which is what keeps the fused node-average path reachable on
  /// an unmasked unweighted sampler.
  void recomposeActiveRows() {
    if (activeRows_.empty()) {
      installWeights(userWeights_);
      return;
    }
    for (std::size_t i = 0; i < numObservations_; ++i)
      composite_[i] =
        (userWeights_ != nullptr ? userWeights_[i] : 1.0) * activeRows_[i];
    installWeights(composite_.data());
  }

  /// Zero-weight rows are documented as ignored; the sigma posterior df counts
  /// only positive-weight rows (all n when unweighted).
  static std::size_t countPositiveWeights(const double* weights,
                                          std::size_t numObservations) {
    if (weights == nullptr) return numObservations;
    std::size_t count = 0;
    for (std::size_t i = 0; i < numObservations; ++i)
      if (weights[i] > 0.0) ++count;
    return count;
  }

  void rescale() {
    // an empty response has no min/max to seed from (yRescaled_[0] is OOB); hold
    // the identity scale
    if (numObservations_ == 0) {
      min_ = 0.0;
      max_ = 0.0;
      range_ = 1.0;
      return;
    }
    if (offset_ != nullptr) {
      misc_subtractVectors(offset_, numObservations_, y_, yRescaled_.data());
    } else {
      std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
    }

    min_ = yRescaled_[0];
    max_ = yRescaled_[0];
    for (std::size_t i = 1; i < numObservations_; ++i) {
      if (yRescaled_[i] < min_) min_ = yRescaled_[i];
      if (yRescaled_[i] > max_) max_ = yRescaled_[i];
    }
    range_ = max_ - min_;
    if (range_ == 0.0) range_ = 1.0;

    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
    misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                     1.0 / range_);
    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);
  }

  const double* y_;
  const double* offset_;
  const double* weights_;      // served: userWeights_ by identity, or composite_
  const double* userWeights_;  // borrowed, the pre-mask pointer
  std::size_t numObservations_;
  std::size_t numPositiveWeights_;
  std::vector<double> activeRows_;  // the raw 0/1 mask; empty when none
  std::vector<double> composite_;   // c_i = w_i a_i, served while masked
  std::vector<double> yRescaled_;
  double min_ = 0.0, max_ = 0.0, range_ = 1.0;
  double initialSigma_ = 1.0;
  ChiSquaredScalePrior sigmaSqPrior_;
};

/// Weights are deliberately unsupported: weighted probit has no exact
/// latent-variable form, so latent draws are unweighted by design. The
/// pre-1.0 engine instead scaled the latents by 1 / sqrt(w), which is not a
/// coherent weighted-likelihood model; that scaling was dropped.
class ProbitResponse final : public ResponseModel {
public:
  ProbitResponse(const double* y, const double* offset,
                 std::size_t numObservations)
    : y_(y), offset_(offset), numObservations_(numObservations) {
    latents_.resize(numObservations);
    working_.resize(numObservations);
    // z = 2 y - 1: -1 for failures, 1 for successes.
    misc_setVectorToConstant(latents_.data(), numObservations, -1.0);
    misc_addVectorsInPlaceWithMultiplier(y, numObservations, 2.0,
                                         latents_.data());
    rebuildWorking();
  }

  double* workingResponse() override { return working_.data(); }
  /// The 0/1 mask IS the precision vector here (probit carries no user
  /// weights), and null - the unit-weight fused path - when none is installed.
  const double* workingWeights() const override {
    return activeRows_.empty() ? nullptr : activeRows_.data();
  }
  const double* offset() const override { return offset_; }

  /// An inactive row's truncated-normal draw is SKIPPED, not drawn and
  /// discarded: the primitive is a rejection sampler whose consumption depends
  /// on the bound, so discarding would desynchronize the stream against a
  /// sampler built on the retained rows. Its latents_[i] therefore keeps its
  /// last drawn value - finite, so rebuildWorking and the running residual
  /// never see a NaN - and is stale until the row is active again.
  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double) override {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      if (!activeRows_.empty() && activeRows_[i] == 0.0) continue;
      double mean = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      double sign = 2.0 * y_[i] - 1.0;
      double z =
        sign * ext_rng_simulateLowerTruncatedNormalScale1(rng, sign * mean, 0.0);
      latents_[i] = !std::isnan(z) ? z : sign * DBL_EPSILON;
    }
    rebuildWorking();
  }

  bool supportsActiveRows() const override { return true; }

  bool setActiveRows(const double* active) override {
    if (active == nullptr) activeRows_.clear();
    else activeRows_.assign(active, active + numObservations_);
    return true;
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool, double*) override {
    y_ = y;
    refreshLatents(rng, totalFits, 1.0);
  }

  void setOffset(const double* offset, bool, double*) override {
    offset_ = offset;
    rebuildWorking();
  }

  void setData(const double* y, const double* offset, const double*,
               std::size_t numObservations, double*) override {
    y_ = y;
    offset_ = offset;
    numObservations_ = numObservations;
    activeRows_.clear();  // length-n and n may have changed
    latents_.resize(numObservations);
    working_.resize(numObservations);
    // cold init, z = 2 y - 1
    misc_setVectorToConstant(latents_.data(), numObservations, -1.0);
    misc_addVectorsInPlaceWithMultiplier(y, numObservations, 2.0,
                                         latents_.data());
    rebuildWorking();
  }

  const double* latents() const override { return latents_.data(); }

  void restoreLatents(const double* latents) override {
    std::memcpy(latents_.data(), latents, numObservations_ * sizeof(double));
    rebuildWorking();
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

  /// log dbinom(y_i, 1, Phi(eta_i)) with eta the latent location f(x) + offset,
  /// via the stable log lower/upper normal tail: log Phi(eta) for a success,
  /// log Phi(-eta) for a failure. An inactive row is not in the model, so its
  /// entry is NaN rather than the finite value its fit would still give.
  void computeLogLikelihood(const double* totalFits, double,
                            std::size_t numObservations,
                            double* out) const override {
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (!activeRows_.empty() && activeRows_[i] == 0.0) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      double eta = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      out[i] = Rf_pnorm5(eta, 0.0, 1.0, y_[i] != 0.0 ? 1 : 0, 1);
    }
  }

private:
  void rebuildWorking() {
    std::memcpy(working_.data(), latents_.data(),
                numObservations_ * sizeof(double));
    if (offset_ != nullptr)
      misc_subtractVectorsInPlace(offset_, numObservations_, working_.data());
  }

  const double* y_;
  const double* offset_;
  std::size_t numObservations_;
  std::vector<double> activeRows_;  // the 0/1 mask, served as the weights
  std::vector<double> latents_;
  std::vector<double> working_;
};

/// Ordered categorical response by a cumulative probit (docs/design/ordinal.md).
/// A latent z_i = f(x_i) + o_i + N(0, 1) is thresholded by ordered cutpoints
/// -inf = gamma_0 < gamma_1 < ... < gamma_{K-1} < gamma_K = +inf into category
/// y_i = k iff gamma_{k-1} < z_i <= gamma_k; y_ holds the 1-based category index
/// in {1..K}. Scheme A identification (section 2): sigma fixed at 1 and gamma_1
/// pinned at 0, so K - 2 interior cutpoints are free and K = 2 collapses to
/// probit. Each sweep first updates the free cutpoints by a marginal random-walk
/// Metropolis step against eta = f + offset with the latents integrated out (the
/// Phi-difference likelihood), then redraws z from the doubly-truncated normal on
/// its category interval - cutpoints before latents (section 3). Reuses probit's
/// working = latents - offset, unit working weights, and fixed sigma; the
/// boundary categories route through the same one-sided truncation primitives and
/// NaN fallback probit uses, keeping the K = 2 rng stream bitwise identical.
/// Weights are unsupported for probit's reason (no coherent weighted latent
/// likelihood); the host refuses them.
class OrdinalResponse final : public ResponseModel {
public:
  OrdinalResponse(const double* y, const double* offset,
                  std::size_t numObservations, std::size_t numCategories)
    : y_(y), offset_(offset), numObservations_(numObservations),
      numCategories_(numCategories) {
    latents_.resize(numObservations);
    working_.resize(numObservations);
    gamma_.resize(numCategories - 1);
    proposalScale_.resize(numCategories - 1);
    coldCutpoints();
    computeScales();
    coldLatents();
    rebuildWorking();
  }

  double* workingResponse() override { return working_.data(); }
  /// As probit: the 0/1 mask is the precision vector, null when none is
  /// installed.
  const double* workingWeights() const override {
    return activeRows_.empty() ? nullptr : activeRows_.data();
  }
  const double* offset() const override { return offset_; }

  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double) override {
    updateCutpoints(rng, totalFits);
    drawLatents(rng, totalFits);
    rebuildWorking();
  }

  bool supportsActiveRows() const override { return true; }

  /// Both cutpoint sums restrict to the active rows, so the proposal scale -
  /// a function of the category counts - is recomputed here; the target
  /// (cutpointLogAcceptance) reads the mask directly.
  bool setActiveRows(const double* active) override {
    if (active == nullptr) activeRows_.clear();
    else activeRows_.assign(active, active + numObservations_);
    computeScales();
    return true;
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  /// Embedded-Gibbs y swap: keep the slow-moving cutpoints (a global parameter
  /// the outer sampler wants persisted across a small y perturbation), refresh
  /// the counts the proposal scale derives from, and redraw z under the new
  /// intervals - probit's minimal-disruption re-draw plus a kept-cutpoint clause.
  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool, double*) override {
    y_ = y;
    computeScales();
    drawLatents(rng, totalFits);
    rebuildWorking();
  }

  void setOffset(const double* offset, bool, double*) override {
    offset_ = offset;
    rebuildWorking();
  }

  /// Data swap: everything stale, so cold-init the cutpoints to the default
  /// prior spacing and z to a deterministic point in each category interval,
  /// matching probit/TResponse cold-init on a data swap.
  void setData(const double* y, const double* offset, const double*,
               std::size_t numObservations, double*) override {
    y_ = y;
    offset_ = offset;
    numObservations_ = numObservations;
    activeRows_.clear();  // length-n and n may have changed
    latents_.resize(numObservations);
    working_.resize(numObservations);
    coldCutpoints();
    computeScales();
    coldLatents();
    rebuildWorking();
  }

  const double* latents() const override { return latents_.data(); }

  void restoreLatents(const double* latents) override {
    std::memcpy(latents_.data(), latents, numObservations_ * sizeof(double));
    rebuildWorking();
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

  /// The K - 1 finite cutpoints gamma_1..gamma_{K-1} (gamma_1 pinned at 0), the
  /// by-name "cutpoints" state block: getState reads numCutpoints() values from
  /// cutpoints(), stateIsValid refuses a state of the wrong length, and setState
  /// restoreCutpoints() reinstalls them. The one-based category index rides the
  /// existing latents block. Also read by the component tests.
  bool carriesCutpoints() const override { return true; }
  const double* cutpoints() const override { return gamma_.data(); }
  std::size_t numCutpoints() const override { return gamma_.size(); }
  void restoreCutpoints(const double* gamma) override {
    std::memcpy(gamma_.data(), gamma, gamma_.size() * sizeof(double));
  }

  /// log P(y_i = k | eta, gamma) = log(Phi(gamma_k - eta) - Phi(gamma_{k-1} -
  /// eta)), the cumulative-probit category likelihood generalizing probit's
  /// two-tail form; the +-inf boundary cutpoints give Phi = 1 / 0.
  void computeLogLikelihood(const double* totalFits, double,
                            std::size_t numObservations,
                            double* out) const override {
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (!isActive(i)) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      double eta = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      int k = category(i);
      out[i] = std::log(Phi(cutAt(k) - eta) - Phi(cutAt(k - 1) - eta));
    }
  }

  /// The log Metropolis acceptance (Phi-difference log-likelihood ratio plus
  /// log-prior ratio; symmetric proposal, so no proposal term) for moving free
  /// cutpoint gamma_s to proposal, exposed so a component test can check the
  /// incremental two-category computation against a full-likelihood evaluation.
  double cutpointLogAcceptanceForTesting(const double* totalFits,
                                         std::size_t s, double proposal) const {
    return cutpointLogAcceptance(totalFits, s, proposal);
  }

  /// The three per-sweep kernels in isolation, so a component test can drive
  /// the masked n-row form against the same kernel over the compacted active
  /// rows; computeScales returns the free cutpoints' proposal scales, its only
  /// observable. Nothing else reaches them - refreshLatents runs two at once.
  const double* computeScalesForTesting() {
    computeScales();
    return proposalScale_.data();
  }
  void updateCutpointsForTesting(ext_rng* rng, const double* totalFits) {
    updateCutpoints(rng, totalFits);
  }
  void drawLatentsForTesting(ext_rng* rng, const double* totalFits) {
    drawLatents(rng, totalFits);
  }

private:
  int category(std::size_t i) const { return static_cast<int>(y_[i]); }

  /// Whether row i is in the data set this sweep; true throughout when no
  /// mask is installed.
  bool isActive(std::size_t i) const {
    return activeRows_.empty() || activeRows_[i] != 0.0;
  }

  static double Phi(double x) { return Rf_pnorm5(x, 0.0, 1.0, 1, 0); }

  /// gamma_s with the sentinels gamma_0 = -inf and gamma_K = +inf; the finite
  /// cutpoints live at gamma_[s - 1].
  double cutAt(int s) const {
    if (s <= 0) return -std::numeric_limits<double>::infinity();
    if (s >= static_cast<int>(numCategories_))
      return std::numeric_limits<double>::infinity();
    return gamma_[s - 1];
  }

  /// log of the normal log-gap prior N(log d; mean, sd) times the 1/d Jacobian
  /// of delta = log d, up to an additive constant; the exact-posterior gate
  /// integrates this same pushed-forward density (section 7).
  double logGapTarget(double d) const {
    if (!(d > 0.0)) return -std::numeric_limits<double>::infinity();
    double delta = std::log(d);
    double z = (delta - priorLogGapMean_) / priorLogGapSd_;
    return -0.5 * z * z - delta;
  }

  /// Cutpoints at the default prior spacing exp(mean): gamma_s = (s - 1) *
  /// spacing with gamma_1 = 0 pinned.
  void coldCutpoints() {
    double spacing = std::exp(priorLogGapMean_);
    for (std::size_t s = 1; s < numCategories_; ++s)
      gamma_[s - 1] = static_cast<double>(s - 1) * spacing;
  }

  /// z at a deterministic point in each category interval (interval midpoint, or
  /// one unit past a boundary cutpoint), the ordinal analogue of probit's
  /// z = 2 y - 1; the first sweep's draw replaces it.
  void coldLatents() {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      int k = category(i);
      if (k <= 1) latents_[i] = gamma_[0] - 1.0;
      else if (k >= static_cast<int>(numCategories_))
        latents_[i] = gamma_[numCategories_ - 2] + 1.0;
      else latents_[i] = 0.5 * (gamma_[k - 2] + gamma_[k - 1]);
    }
  }

  /// Fixed count-derived random-walk scale for each free cutpoint gamma_s:
  /// c / sqrt(n_s + n_{s+1}) (section 3), floored so an empty adjacent pair
  /// still moves under the prior. Counts only the active rows, so the proposal
  /// is the retained subsample's.
  void computeScales() {
    std::vector<double> counts(numCategories_ + 1, 0.0);
    for (std::size_t i = 0; i < numObservations_; ++i)
      if (isActive(i)) counts[static_cast<std::size_t>(category(i))] += 1.0;
    for (std::size_t s = 2; s < numCategories_; ++s) {
      double denom = counts[s] + counts[s + 1];
      proposalScale_[s - 1] =
        proposalScaleConstant_ / std::sqrt(denom > 1.0 ? denom : 1.0);
    }
  }

  /// One marginal random-walk Metropolis pass over the free cutpoints gamma_2..
  /// gamma_{K-1} against eta = totalFits + offset, latents integrated out. Plain
  /// symmetric normal proposal; an out-of-order proposal has acceptance 0 and is
  /// rejected before the accept draw; an in-bounds proposal is accepted on the
  /// two-adjacent-cell likelihood ratio times the prior ratio. Consumes no rng
  /// when there are no free cutpoints, so K = 2 stays bitwise probit.
  void updateCutpoints(ext_rng* rng, const double* totalFits) {
    for (std::size_t s = 2; s < numCategories_; ++s) {
      int si = static_cast<int>(s);
      double proposal = gamma_[s - 1] +
        proposalScale_[s - 1] * ext_rng_simulateStandardNormal(rng);
      if (proposal <= cutAt(si - 1) || proposal >= cutAt(si + 1))
        continue;
      double logAccept = cutpointLogAcceptance(totalFits, s, proposal);
      if (std::log(ext_rng_simulateContinuousUniform(rng)) < logAccept)
        gamma_[s - 1] = proposal;
    }
  }

  /// log acceptance for gamma_s -> proposal: only categories s and s + 1 change
  /// (gamma_s is their shared cutpoint), plus the two log-gaps gamma_s touches.
  double cutpointLogAcceptance(const double* totalFits, std::size_t s,
                               double proposal) const {
    int si = static_cast<int>(s);
    double current = gamma_[s - 1];
    double lo = cutAt(si - 1), hi = cutAt(si + 1);
    double logLik = 0.0;
    for (std::size_t i = 0; i < numObservations_; ++i) {
      if (!isActive(i)) continue;  // the target is the subsample's likelihood
      int k = category(i);
      if (k != si && k != si + 1) continue;
      double eta = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      double oldp, newp;
      if (k == si) {
        double base = Phi(lo - eta);
        oldp = Phi(current - eta) - base;
        newp = Phi(proposal - eta) - base;
      } else {
        double base = Phi(hi - eta);
        oldp = base - Phi(current - eta);
        newp = base - Phi(proposal - eta);
      }
      logLik += std::log(newp) - std::log(oldp);
    }
    double logPrior = logGapTarget(proposal - lo) - logGapTarget(current - lo);
    if (s < numCategories_ - 1)  // finite upper gap only below the top cutpoint
      logPrior += logGapTarget(hi - proposal) - logGapTarget(hi - current);
    return logLik + logPrior;
  }

  /// Redraw z_i from N(eta_i, 1) on its category interval (gamma_{y_i-1},
  /// gamma_{y_i}]. Boundary categories keep probit's one-sided rejection
  /// primitives and sign * DBL_EPSILON NaN fallback so K = 2 is bitwise probit;
  /// interior categories use the doubly-truncated primitive with a midpoint
  /// fallback. An inactive row's draw is skipped for probit's reason (the
  /// primitives are rejection samplers), leaving its z stale but finite.
  void drawLatents(ext_rng* rng, const double* totalFits) {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      if (!isActive(i)) continue;
      double mean = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      int k = category(i);
      if (k <= 1) {
        double z =
          ext_rng_simulateUpperTruncatedNormalScale1(rng, mean, gamma_[0]);
        latents_[i] = !std::isnan(z) ? z : -DBL_EPSILON;
      } else if (k >= static_cast<int>(numCategories_)) {
        double z = ext_rng_simulateLowerTruncatedNormalScale1(
          rng, mean, gamma_[numCategories_ - 2]);
        latents_[i] = !std::isnan(z) ? z : DBL_EPSILON;
      } else {
        double lower = gamma_[k - 2], upper = gamma_[k - 1];
        double z =
          ext_rng_simulateTruncatedNormalScale1(rng, mean, lower, upper);
        latents_[i] = !std::isnan(z) ? z : 0.5 * (lower + upper);
      }
    }
  }

  void rebuildWorking() {
    std::memcpy(working_.data(), latents_.data(),
                numObservations_ * sizeof(double));
    if (offset_ != nullptr)
      misc_subtractVectorsInPlace(offset_, numObservations_, working_.data());
  }

  static constexpr double priorLogGapMean_ = 0.0;
  static constexpr double priorLogGapSd_ = 1.5;
  static constexpr double proposalScaleConstant_ = 2.5;

  const double* y_;
  const double* offset_;
  std::size_t numObservations_;
  std::size_t numCategories_;
  std::vector<double> gamma_;          // finite cutpoints gamma_1..gamma_{K-1}
  std::vector<double> proposalScale_;  // per free cutpoint; index 0 unused
  std::vector<double> activeRows_;     // the 0/1 mask, served as the weights
  std::vector<double> latents_;
  std::vector<double> working_;
};

/// Logistic via Polya-Gamma augmentation (Polson, Scott & Windle 2013):
/// given the fit eta_i = f(x_i) + offset_i, omega_i ~ PG(1, eta_i), and
/// kappa_i / omega_i with kappa_i = y_i - 0.5 is conditionally
/// N(eta_i, 1 / omega_i). The backfitting engine therefore sees working
/// response kappa_i / omega_i - offset_i under per-iteration precision
/// weights omega_i, running the same weighted conjugate updates as a
/// weighted gaussian with sigma fixed at 1. Case weights are observation
/// COUNTS, the augmentation's own shape parameter rather than a precision.
/// latents() exposes the omega draws.
class LogisticResponse final : public ResponseModel {
public:
  /// weights are observation counts (positive integers) or null for unit
  /// weights; the host validates integrality. A count-w observation's
  /// Polya-Gamma latent is PG(w, psi), drawn exactly as the sum of w
  /// independent PG(1, psi) variates, so null weights reproduce the
  /// unweighted stream bit for bit.
  LogisticResponse(const double* y, const double* offset,
                   const double* weights, std::size_t numObservations)
    : y_(y), offset_(offset), weights_(weights),
      numObservations_(numObservations) {
    omega_.resize(numObservations);
    working_.resize(numObservations);
    coldStart();
  }

  double* workingResponse() override { return working_.data(); }
  /// The Polya-Gamma draws, or a_i omega_i from a SEPARATE composite while a
  /// mask is installed. The zero is never written into omega_ itself: the
  /// working response divides by it, and 0 * inf in the node kernels is a NaN.
  const double* workingWeights() const override {
    return activeRows_.empty() ? omega_.data() : composite_.data();
  }
  bool workingWeightsVaryPerSweep() const override { return true; }
  const double* offset() const override { return offset_; }

  /// An inactive row's PG draw is SKIPPED, not drawn and discarded: the
  /// primitive is a rejection sampler whose consumption depends on psi, so a
  /// discard would desynchronize the stream against a sampler built on the
  /// retained rows. Its omega_ and working_ keep their last values, which are
  /// finite and positive, and are stale until the row is active again.
  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double) override {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      if (!isActive(i)) continue;
      double offset = offset_ != nullptr ? offset_[i] : 0.0;
      double psi = totalFits[i] + offset;
      long reps = weights_ != nullptr ? std::lround(weights_[i]) : 1L;
      double omega = ext_rng_simulatePolyaGamma(rng, psi);
      for (long c = 1; c < reps; ++c)
        omega += ext_rng_simulatePolyaGamma(rng, psi);
      omega_[i] = omega;
      double weight = weights_ != nullptr ? weights_[i] : 1.0;
      working_[i] = weight * (y_[i] - 0.5) / omega - offset;
    }
    recompose();
  }

  bool supportsActiveRows() const override { return true; }

  /// The mask is NOT redundant with the count weights: a zero count is refused
  /// at creation, no count change is accepted afterwards, and a zero-count row
  /// would still consume one PG variate here.
  bool setActiveRows(const double* active) override {
    if (active == nullptr) {
      activeRows_.clear();
      return true;
    }
    activeRows_.assign(active, active + numObservations_);
    composite_.resize(numObservations_);
    recompose();
    return true;
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool, double*) override {
    y_ = y;
    refreshLatents(rng, totalFits, 1.0);
  }

  /// The counts ARE the Polya-Gamma shape, so a swap restates the model, and
  /// every quantity stated against the old counts is replaced here rather than
  /// at the next sweep, which draws its trees before refreshing. An ACTIVE row
  /// takes the same PG(w, psi) draw refreshLatents makes. An INACTIVE row
  /// draws NOTHING - consuming variates for it would desynchronize the stream
  /// against a sampler built on the retained rows, the property the skip in
  /// refreshLatents protects - and is returned to the deterministic cold start
  /// against its NEW count, so a row that reactivates cannot carry an omega
  /// shaped by counts the sampler no longer holds.
  void setWeights(const double* weights, ext_rng* rng,
                  const double* totalFits) override {
    weights_ = weights;
    for (std::size_t i = 0; i < numObservations_; ++i)
      if (!isActive(i)) coldStartRow(i);
    refreshLatents(rng, totalFits, 1.0);
  }

  void setOffset(const double* offset, bool, double*) override {
    // omega and kappa / omega stand; only the shift into the working
    // response moves
    reshiftWorkingForOffset(working_.data(), offset_, offset, numObservations_);
    offset_ = offset;
  }

  void setData(const double* y, const double* offset, const double* weights,
               std::size_t numObservations, double*) override {
    y_ = y;
    offset_ = offset;
    weights_ = weights;
    numObservations_ = numObservations;
    activeRows_.clear();  // length-n and n may have changed
    composite_.clear();
    omega_.resize(numObservations);
    working_.resize(numObservations);
    coldStart();
  }

  const double* latents() const override { return omega_.data(); }

  void restoreLatents(const double* latents) override {
    std::memcpy(omega_.data(), latents, numObservations_ * sizeof(double));
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double weight = weights_ != nullptr ? weights_[i] : 1.0;
      working_[i] = weight * (y_[i] - 0.5) / omega_[i] -
                    (offset_ != nullptr ? offset_[i] : 0.0);
    }
    recompose();
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

  /// w_i log dbinom(y_i, 1, plogis(eta_i)) with eta the log-odds f(x) + offset
  /// and w_i the integer trial count (1 unweighted): -log(1 + exp(-eta)) for a
  /// success, -log(1 + exp(eta)) for a failure. An inactive row is not in the
  /// model, so its entry is NaN rather than the value its fit would give.
  void computeLogLikelihood(const double* totalFits, double,
                            std::size_t numObservations,
                            double* out) const override {
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (!isActive(i)) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      double eta = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      double ll = y_[i] != 0.0 ? -logOnePlusExp(-eta) : -logOnePlusExp(eta);
      out[i] = (weights_ != nullptr ? weights_[i] : 1.0) * ll;
    }
  }

private:
  /// Whether row i is in the data set this sweep; true throughout when no
  /// mask is installed.
  bool isActive(std::size_t i) const {
    return activeRows_.empty() || activeRows_[i] != 0.0;
  }

  /// c_i = a_i omega_i, the served precisions while a mask is installed; a
  /// no-op otherwise, when omega_ is served by identity.
  void recompose() {
    if (activeRows_.empty()) return;
    for (std::size_t i = 0; i < numObservations_; ++i)
      composite_[i] = activeRows_[i] * omega_[i];
  }

  /// Deterministic cold start, the analogue of probit's z = 2 y - 1: omega
  /// at PG(w, 0)'s mean of w/4, so the working response starts at
  /// 4 (y - 1/2) - offset independent of the weight; real draws replace it
  /// after the first tree sweep.
  void coldStart() {
    for (std::size_t i = 0; i < numObservations_; ++i) coldStartRow(i);
  }

  /// One row of it, also the state a weight swap leaves on a row it does not
  /// draw for.
  void coldStartRow(std::size_t i) {
    omega_[i] = 0.25 * (weights_ != nullptr ? weights_[i] : 1.0);
    working_[i] =
      4.0 * (y_[i] - 0.5) - (offset_ != nullptr ? offset_[i] : 0.0);
  }

  const double* y_;
  const double* offset_;
  const double* weights_;
  std::size_t numObservations_;
  std::vector<double> activeRows_;  // the 0/1 mask; empty when none
  std::vector<double> composite_;   // c_i = a_i omega_i, served while masked
  std::vector<double> omega_;
  std::vector<double> working_;
};

/// The response family for a multinomial (softmax) chain: a thin shell whose
/// single-location seams are vestigial. The K interleaved Polya-Gamma draws,
/// the per-forest working responses, and the level-centering move all live in
/// MultinomialForestCombiner (combiner.hpp), which owns the category labels;
/// this family carries no sigma (fixed, like the binary families) and no
/// latents. refreshLatents and drawSigma are no-ops, and the log-likelihood
/// channel is left NaN-flagged (the default computeLogLikelihood): storeSample
/// scores one forest's fits, which cannot see the K-blend, so
/// logLikelihoodIsDefined() = false is the reporting map's answer (the same
/// choice the Bayesian causal forest (BCF) family makes). fitScale/fitShift
/// are identity so storeSample passes the softmax
/// probabilities the combiner writes through unchanged.
class MultinomialResponse final : public ResponseModel {
public:
  explicit MultinomialResponse(std::size_t numObservations)
    : numObservations_(numObservations) {
    working_.assign(numObservations, 0.0);
  }

  // The combiner owns the working response and precisions; the sweep still
  // reads these, then hands them to formForestResponse, which ignores them.
  double* workingResponse() override { return working_.data(); }
  const double* workingWeights() const override { return nullptr; }
  const double* offset() const override { return nullptr; }

  // Layout guard - combinedFits is K per observation, not one channel: the
  // sweep hands combinedFits' K x n softmax buffer here as a single-location
  // pointer. This no-op ignores it; a future non-no-op refreshLatents must NOT
  // read it as one channel (it is K per observation).
  void refreshLatents(ext_rng*, const double* combinedFits, double) override {
    (void) combinedFits;
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  /// The active-row channel, advertised here and IMPLEMENTED IN THE COMBINER:
  /// this family holds no precisions of its own (workingWeights() is null and
  /// the coupling discards the chain's), so the mask has nothing to compose
  /// into here and lands on the K interleaved Polya-Gamma precisions
  /// MultinomialForestCombiner owns. Chain::setActiveRows forwards it there
  /// after this returns, and a MultinomialResponse is constructed only
  /// alongside that coupling, so the capability the probe reads off the
  /// response is the one the coupling honours. Accepting here rather than
  /// widening the probe keeps SamplerShape::supportsActiveRows derived from the
  /// single predicate the setter refuses on.
  bool supportsActiveRows() const override { return true; }
  bool setActiveRows(const double*) override { return true; }

  // The multinomial response is the borrowed n x K count matrix on the spec,
  // which a flat double* cannot express; the softmax is invariant to a common
  // per-observation shift, so a flat offset points exactly along the null
  // direction and a meaningful one would be n x K too. Both setResponse and
  // setOffset are therefore intentional no-ops that only satisfy the
  // interface. The bridge refuses BOTH mutations on a multi-forest sampler
  // whose coupling does not opt in, and this one does not, so neither no-op is
  // reachable rather than surprising.
  void setResponse(const double*, ext_rng*, const double*, bool,
                   double*) override {}
  void setOffset(const double*, bool, double*) override {}
  void setData(const double*, const double*, const double*,
               std::size_t numObservations, double*) override {
    numObservations_ = numObservations;
    working_.assign(numObservations, 0.0);
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

private:
  std::size_t numObservations_;
  std::vector<double> working_;
};

/// Accelerated failure time, log-normal error (docs/design/survival.md):
///   log T_i = f(x_i) + offset_i + sigma * eps_i,   eps_i ~ N(0, 1).
/// The forest fits log T on the log scale exactly as GaussianResponse, which
/// this contains over a mutable log-time buffer it borrows to that Gaussian
/// as the response. Right-censored observations (status 0) carry a latent
/// log-time redrawn each sweep from the normal truncated below at the observed
/// log censoring time, so it enters the conjugate sigma draw as data;
/// uncensored ones (status 1) are fixed Gaussian data. With no censored
/// observations refreshLatents is a no-op and every hook delegates, so the fit
/// is BIT-IDENTICAL to a Gaussian fit on log T (the reduction gate). Weights
/// are unsupported (a weighted truncated-latent draw is not a coherent
/// likelihood; the host rejects them). State rides the existing `latents`
/// (logT_) and `fit.scale` (the Gaussian scale) blocks unchanged.
class AFTResponse final : public ResponseModel {
public:
  /// logTime holds log of the observed time (event or censoring time);
  /// status[i] != 0 marks an uncensored event, 0 a right-censored observation
  /// whose logTime is its truncation lower bound. offset may be null.
  AFTResponse(const double* logTime, const double* status, const double* offset,
              std::size_t numObservations, double sigmaEstimate, double sigmaDf,
              double sigmaRawScale)
    : numObservations_(numObservations),
      logT_(logTime, logTime + numObservations) {
    for (std::size_t i = 0; i < numObservations; ++i)
      if (status == nullptr || status[i] == 0.0) {
        censoredIndices_.push_back(i);
        censorBound_.push_back(logTime[i]);
      }
    gaussian_ = std::make_unique<GaussianResponse>(
      logT_.data(), offset, nullptr, numObservations, sigmaEstimate, sigmaDf,
      sigmaRawScale);
  }

  double* workingResponse() override { return gaussian_->workingResponse(); }
  const double* workingWeights() const override {
    return gaussian_->workingWeights();
  }
  const double* offset() const override { return gaussian_->offset(); }

  /// Redraw each censored log-time from N(f + offset, sigma^2) truncated below
  /// at its log censoring time, mapping the internal-scale fit back to the log
  /// scale through the public accessors, then rebuild the working response
  /// under the fixed scale. A no-op with no censored observations. An inactive
  /// censored row's redraw is SKIPPED for probit's reason (the truncated-normal
  /// primitive is a rejection sampler), leaving its log-time stale but finite.
  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double sigma) override {
    if (censoredIndices_.empty()) return;
    double scale = gaussian_->fitScale();
    double shift = gaussian_->fitShift();
    double sd = sigma * scale;
    const double* offset = gaussian_->offset();
    for (std::size_t k = 0; k < censoredIndices_.size(); ++k) {
      std::size_t i = censoredIndices_[k];
      if (!isActive(i)) continue;
      double mean = scale * totalFits[i] + shift +
                    (offset != nullptr ? offset[i] : 0.0);
      double draw =
        ext_rng_simulateLowerTruncatedNormal(rng, mean, sd, censorBound_[k]);
      logT_[i] = !std::isnan(draw) ? draw : censorBound_[k];
    }
    rebuildWorking();
  }

  bool supportsActiveRows() const override { return true; }

  /// The mask goes into the contained Gaussian's weights, so the node
  /// statistics and the sigma posterior's degrees of freedom inherit it through
  /// the same recount a masked gaussian gets; the response transform stays the
  /// FULL-data one (rescale spans all n rows), as it does under zero weights.
  bool setActiveRows(const double* active) override {
    if (active == nullptr) activeRows_.clear();
    else activeRows_.assign(active, active + numObservations_);
    return gaussian_->setActiveRows(active);
  }

  double drawSigma(ext_rng* rng, const double* totalFits,
                   double sigma) override {
    return gaussian_->drawSigma(rng, totalFits, sigma);
  }

  double sigmaDegreesOfFreedomForTesting() const override {
    return gaussian_->sigmaDegreesOfFreedomForTesting();
  }

  /// Replaces the observed log-times; the censoring structure is fixed at
  /// creation, so bounds refresh from the new times and the censored latents
  /// redraw against the current fit (the probit pattern).
  void setResponse(const double* logTime, ext_rng* rng, const double* totalFits,
                   bool updateScale, double* sigmaInOut) override {
    std::memcpy(logT_.data(), logTime, numObservations_ * sizeof(double));
    for (std::size_t k = 0; k < censoredIndices_.size(); ++k)
      censorBound_[k] = logT_[censoredIndices_[k]];
    gaussian_->setResponse(logT_.data(), rng, totalFits, updateScale,
                           sigmaInOut);
    if (!censoredIndices_.empty()) refreshLatents(rng, totalFits, *sigmaInOut);
  }

  /// Offset shifts the fit location, not the (offset-independent) log-times, so
  /// the Gaussian re-anchors its scale and rebuilds the working response while
  /// the latents keep their values (redrawn next sweep).
  void setOffset(const double* offset, bool updateScale,
                 double* sigmaInOut) override {
    gaussian_->setOffset(offset, updateScale, sigmaInOut);
  }

  /// Whole-data replacement keeps the censoring structure by index (no new
  /// status is expressible through this signature; the host refuses a public
  /// setData on AFT). Only a length-preserving replacement is coherent; a
  /// shrink drops censored indices that fell off the end.
  void setData(const double* logTime, const double* offset, const double*,
               std::size_t numObservations, double* sigmaInOut) override {
    logT_.assign(logTime, logTime + numObservations);
    numObservations_ = numObservations;
    activeRows_.clear();  // length-n and n may have changed
    while (!censoredIndices_.empty() &&
           censoredIndices_.back() >= numObservations) {
      censoredIndices_.pop_back();
      censorBound_.pop_back();
    }
    for (std::size_t k = 0; k < censoredIndices_.size(); ++k)
      censorBound_[k] = logT_[censoredIndices_[k]];
    gaussian_->setData(logT_.data(), offset, nullptr, numObservations,
                       sigmaInOut);
  }

  void setSigmaPrior(double sigmaEstimate, double degreesOfFreedom,
                     double rawScale) override {
    gaussian_->setSigmaPrior(sigmaEstimate, degreesOfFreedom, rawScale);
  }

  const double* latents() const override { return logT_.data(); }
  void restoreLatents(const double* latents) override {
    std::memcpy(logT_.data(), latents, numObservations_ * sizeof(double));
    rebuildWorking();
  }

  void getScale(double& min, double& max) const override {
    gaussian_->getScale(min, max);
  }
  void restoreScale(double min, double max) override {
    gaussian_->restoreScale(min, max);
  }

  double initialSigma() const override { return gaussian_->initialSigma(); }
  double fitScale() const override { return gaussian_->fitScale(); }
  double fitShift() const override { return gaussian_->fitShift(); }
  double sigmaScale() const override { return gaussian_->sigmaScale(); }

  /// log density of the observed log event time for an event, log survival
  /// past the log censoring bound for a censored observation, both on the log
  /// scale with mu the original-scale fit and sigma the log-scale residual sd.
  /// An inactive row is not in the model, so its entry is NaN.
  void computeLogLikelihood(const double* totalFits, double sigma,
                            std::size_t numObservations,
                            double* out) const override {
    double scale = gaussian_->fitScale();
    double shift = gaussian_->fitShift();
    double sigmaLog = sigma * scale;
    const double* offset = gaussian_->offset();
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (!isActive(i)) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      double mu =
        scale * totalFits[i] + shift + (offset != nullptr ? offset[i] : 0.0);
      out[i] = Rf_dnorm4(logT_[i], mu, sigmaLog, 1);
    }
    for (std::size_t k = 0; k < censoredIndices_.size(); ++k) {
      std::size_t i = censoredIndices_[k];
      if (!isActive(i)) continue;
      double mu =
        scale * totalFits[i] + shift + (offset != nullptr ? offset[i] : 0.0);
      // log P(log T > log C) = log upper normal tail
      out[i] = Rf_pnorm5(censorBound_[k], mu, sigmaLog, 0, 1);
    }
  }

private:
  /// Whether row i is in the data set this sweep; true throughout when no
  /// mask is installed.
  bool isActive(std::size_t i) const {
    return activeRows_.empty() || activeRows_[i] != 0.0;
  }

  /// Rebuild the Gaussian working response from the current log-times under the
  /// fixed scale (updateScale = false), the only side effect the redrawn
  /// latents need.
  void rebuildWorking() {
    double unused = 0.0;
    gaussian_->setResponse(logT_.data(), nullptr, nullptr, false, &unused);
  }

  std::unique_ptr<GaussianResponse> gaussian_;
  std::size_t numObservations_;
  std::vector<double> logT_;             // log survival times; latent when censored
  std::vector<std::size_t> censoredIndices_;
  std::vector<double> censorBound_;      // log censoring time, per censored index
  std::vector<double> activeRows_;       // the 0/1 mask; empty when none
};

/// Sampled residual degrees of freedom nu on a fixed capped grid under a
/// normalized gamma(2, 0.1) prior (docs/design/robust-errors.md section 4).
/// With the mixing precisions lambda_i ~ Gamma(nu/2, nu/2), the log full
/// conditional at grid point nu is
///   n [(nu/2) log(nu/2) - lgamma(nu/2)] + (nu/2)(sumLogLambda - sumLambda)
///     + log p(nu)
/// (the lambda^{-1} density term is constant across the grid, so it drops from
/// the normalized draw). The per-point lgamma and prior constants precompute
/// once, so each draw is O(grid) multiply-adds from the two lambda statistics
/// (the DartPrior idiom). The grid floor of 3 keeps the marginal variance
/// finite.
struct ResidualDfPrior {
  static constexpr std::size_t gridSize = 9;
  static constexpr double grid[gridSize] = {3.0,  4.0,  5.0,  6.0, 8.0,
                                            10.0, 12.0, 15.0, 20.0};
  static constexpr std::size_t medianIndex = 4;  // nu = 8, the cold-start value

  ResidualDfPrior() {
    double priorTotal = 0.0;
    for (std::size_t k = 0; k < gridSize; ++k) {
      double halfNu = 0.5 * grid[k];
      kernel_[k] = halfNu * std::log(halfNu) - std::lgamma(halfNu);
      logPrior_[k] = grid[k] * std::exp(-0.1 * grid[k]);  // gamma(2, 0.1) kernel
      priorTotal += logPrior_[k];
    }
    for (std::size_t k = 0; k < gridSize; ++k)
      logPrior_[k] = std::log(logPrior_[k] / priorTotal);
  }

  /// Draw a grid index from the full conditional given the informative-row
  /// count and the two lambda sufficient statistics.
  std::size_t drawIndex(ext_rng* rng, double numObservations,
                        double sumLogLambda, double sumLambda) {
    for (std::size_t k = 0; k < gridSize; ++k)
      weight_[k] = numObservations * kernel_[k] +
                   0.5 * grid[k] * (sumLogLambda - sumLambda) + logPrior_[k];
    return drawFromLogWeights(rng, weight_.data(), gridSize);
  }

private:
  std::array<double, gridSize> kernel_;
  std::array<double, gridSize> logPrior_;
  std::array<double, gridSize> weight_;
};

/// Student-t continuous errors by the Gaussian scale-mixture augmentation
/// (docs/design/robust-errors.md): sqrt(w_i) r_i / sigma ~ t_nu comes from
///   r_i | lambda_i ~ N(0, sigma^2 / (w_i lambda_i)), lambda_i ~ Gamma(nu/2, nu/2),
/// so a contained GaussianResponse over the same rescaled response sees the
/// composite precision c_i = w_i lambda_i through setWeights,
/// and both its node statistics and its conjugate sigma draw are exact under
/// the mixture with no GaussianResponse edit. refreshLatents redraws lambda
/// each sweep (workingWeightsVaryPerSweep() true, as logistic's omega does) and
/// in grid mode draws nu afterward; fixed mode holds nu. Only
/// computeLogLikelihood overrides, to the marginal t density; latents() exposes
/// lambda for the state block. A zero user weight gives c_i = 0, dropping the
/// row from the sigma df and the nu statistics. sigma is the CONDITIONAL scale
/// (marginal variance sigma^2 nu/(nu - 2)).
class TResponse final : public ResponseModel {
public:
  /// residualDf > 0 fixes nu there; a non-positive or non-finite residualDf
  /// estimates nu on the grid, cold-started at its median. offset and weights
  /// may be null.
  TResponse(const double* y, const double* offset, const double* weights,
            std::size_t numObservations, double sigmaEstimate, double sigmaDf,
            double sigmaRawScale, double residualDf)
    : y_(y), userWeights_(weights), numObservations_(numObservations),
      estimateNu_(!(residualDf > 0.0)),
      nu_(estimateNu_ ? ResidualDfPrior::grid[ResidualDfPrior::medianIndex]
                      : residualDf) {
    lambda_.assign(numObservations, 1.0);
    composite_.resize(numObservations);
    for (std::size_t i = 0; i < numObservations; ++i)
      composite_[i] = userWeights_ != nullptr ? userWeights_[i] : 1.0;
    gaussian_ = std::make_unique<GaussianResponse>(
      y, offset, composite_.data(), numObservations, sigmaEstimate, sigmaDf,
      sigmaRawScale);
  }

  double* workingResponse() override { return gaussian_->workingResponse(); }
  const double* workingWeights() const override {
    return gaussian_->workingWeights();
  }
  bool workingWeightsVaryPerSweep() const override { return true; }
  const double* offset() const override { return gaussian_->offset(); }

  /// Draw each lambda_i from its conjugate Gamma((nu + 1)/2, scale) with
  /// scale = 2 / (nu + w_i r_i^2 / sigma^2) over the internal residual
  /// r_i = z_i - f_i, compose c_i = w_i lambda_i, draw nu from the two
  /// informative-row lambda statistics (grid mode), then hand the composite to
  /// the Gaussian so its weights and sigma draw see the mixture.
  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double sigma) override {
    const double* z = gaussian_->workingResponse();
    double sigmaSq = sigma * sigma;
    double shape = 0.5 * (nu_ + 1.0);
    double sumLogLambda = 0.0, sumLambda = 0.0, numInformative = 0.0;
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double w = userWeights_ != nullptr ? userWeights_[i] : 1.0;
      double a = activeRows_.empty() ? 1.0 : activeRows_[i];
      double r = z[i] - totalFits[i];
      // lambda is drawn for EVERY row: the mask annihilates its value through
      // the composite rather than suppressing the draw, which keeps an
      // inactive row's lambda current for the sweep it reactivates on
      double lambda =
        ext_rng_simulateGamma(rng, shape, 2.0 / (nu_ + w * r * r / sigmaSq));
      lambda_[i] = lambda;
      composite_[i] = w * lambda * a;
      // the nu statistics read the COMPOSED weight, so an inactive row leaves
      // them exactly as a zero user weight does
      if (w * a > 0.0) {
        sumLogLambda += std::log(lambda);
        sumLambda += lambda;
        numInformative += 1.0;
      }
    }
    if (estimateNu_)
      nu_ = ResidualDfPrior::grid[nuPrior_.drawIndex(rng, numInformative,
                                                     sumLogLambda, sumLambda)];
    gaussian_->setWeights(composite_.data(), nullptr, nullptr);
  }

  double drawSigma(ext_rng* rng, const double* totalFits,
                   double sigma) override {
    return gaussian_->drawSigma(rng, totalFits, sigma);
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool updateScale, double* sigmaInOut) override {
    y_ = y;
    gaussian_->setResponse(y, rng, totalFits, updateScale, sigmaInOut);
    coldInit();
  }

  void setOffset(const double* offset, bool updateScale,
                 double* sigmaInOut) override {
    gaussian_->setOffset(offset, updateScale, sigmaInOut);
  }

  void setData(const double* y, const double* offset, const double* weights,
               std::size_t numObservations, double* sigmaInOut) override {
    y_ = y;
    userWeights_ = weights;
    numObservations_ = numObservations;
    activeRows_.clear();  // length-n and n may have changed
    lambda_.assign(numObservations, 1.0);
    composite_.assign(numObservations, 0.0);
    gaussian_->setData(y, offset, weights, numObservations, sigmaInOut);
    coldInit();
  }

  void setWeights(const double* weights, ext_rng*, const double*) override {
    userWeights_ = weights;
    recompose();
  }

  bool supportsActiveRows() const override { return true; }

  /// The mask joins the mixture composite, c_i = w_i lambda_i a_i, so the
  /// contained Gaussian inherits it through the pointer it is already handed -
  /// node statistics, sigma df and all. The lambda draw itself continues at
  /// every row; refreshLatents states why.
  bool setActiveRows(const double* active) override {
    if (active == nullptr) activeRows_.clear();
    else activeRows_.assign(active, active + numObservations_);
    recompose();
    return true;
  }

  void setSigmaPrior(double sigmaEstimate, double degreesOfFreedom,
                     double rawScale) override {
    gaussian_->setSigmaPrior(sigmaEstimate, degreesOfFreedom, rawScale);
  }

  const double* latents() const override { return lambda_.data(); }
  void restoreLatents(const double* latents) override {
    std::memcpy(lambda_.data(), latents, numObservations_ * sizeof(double));
    recompose();
  }

  /// The current residual df and, for the state block, whether it is sampled;
  /// restoreResidualDf reinstalls a stored nu (a grid value in grid mode). The
  /// state block serializes nu whenever carriesResidualDf() is true.
  bool carriesResidualDf() const override { return true; }
  double residualDf() const override { return nu_; }
  bool estimatesResidualDf() const { return estimateNu_; }
  void restoreResidualDf(double nu) override { nu_ = nu; }

  void getScale(double& min, double& max) const override {
    gaussian_->getScale(min, max);
  }
  void restoreScale(double min, double max) override {
    gaussian_->restoreScale(min, max);
  }

  double initialSigma() const override { return gaussian_->initialSigma(); }
  double fitScale() const override { return gaussian_->fitScale(); }
  double fitShift() const override { return gaussian_->fitShift(); }
  double sigmaScale() const override { return gaussian_->sigmaScale(); }

  /// Marginal density of the observed response: a location-scale t_nu with
  /// original-scale location mu = fit (offset carried) and scale
  /// s_i = sigma range_ / sqrt(w_i), so log p = log dt((y - mu)/s_i, nu) - log s_i.
  void computeLogLikelihood(const double* totalFits, double sigma,
                            std::size_t numObservations,
                            double* out) const override {
    double scale = gaussian_->fitScale();
    double shift = gaussian_->fitShift();
    double sigmaOriginal = sigma * scale;
    const double* offset = gaussian_->offset();
    for (std::size_t i = 0; i < numObservations; ++i) {
      // lambda is positive, so the composite vanishes exactly when the user
      // weight or the mask does: not in the model, hence NaN
      if ((userWeights_ != nullptr && userWeights_[i] == 0.0) ||
          (!activeRows_.empty() && activeRows_[i] == 0.0)) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      double mu = scale * totalFits[i] + shift +
                  (offset != nullptr ? offset[i] : 0.0);
      double sd = userWeights_ != nullptr
                    ? sigmaOriginal / std::sqrt(userWeights_[i])
                    : sigmaOriginal;
      out[i] = Rf_dt((y_[i] - mu) / sd, nu_, 1) - std::log(sd);
    }
  }

private:
  /// c_i = w_i lambda_i a_i, then hand the composite to the contained Gaussian
  /// so its weights, node statistics, and sigma draw all see the mixture.
  void recompose() {
    for (std::size_t i = 0; i < numObservations_; ++i)
      composite_[i] =
        (userWeights_ != nullptr ? userWeights_[i] : 1.0) * lambda_[i] *
        (activeRows_.empty() ? 1.0 : activeRows_[i]);
    gaussian_->setWeights(composite_.data(), nullptr, nullptr);
  }

  /// Cold state after a data swap: lambda = 1 (c_i = w_i) and, in grid mode,
  /// nu at the grid median; the next sweep draws both.
  void coldInit() {
    std::fill(lambda_.begin(), lambda_.end(), 1.0);
    if (estimateNu_) nu_ = ResidualDfPrior::grid[ResidualDfPrior::medianIndex];
    recompose();
  }

  const double* y_;
  const double* userWeights_;
  std::size_t numObservations_;
  bool estimateNu_;
  double nu_;
  std::unique_ptr<GaussianResponse> gaussian_;
  std::vector<double> lambda_;      // per-observation mixing precisions
  std::vector<double> activeRows_;  // the 0/1 mask; empty when none
  std::vector<double> composite_;   // c_i = w_i lambda_i a_i, the Gaussian's weights
  ResidualDfPrior nuPrior_;         // grid machinery, used only when estimating
};

/// Sampled negative-binomial dispersion r on a capped positive-integer grid
/// under a normalized gamma(2, 0.1) prior (docs/design/negative-binomial.md
/// sections 2A, 3). This is the r-update seam: a grid full conditional today, a
/// Chinese-restaurant-table (CRT) or real-valued-r strategy later, with
/// integrality assumed only inside it. Under
/// the logit-p parameterization p_i is r-free, so the log full conditional at
/// grid point r_k separates (dropping the r-free normalizer) into
///   L_k + r_k * S + log p(r_k),
/// L_k = sum_c n_c [lgamma(c + r_k) - lgamma(r_k)] over the count histogram n_c,
/// S = sum_i log(1 - p_i). L_k derives ONLY from the fixed counts, so it
/// precomputes once per response (computeKernel); each sweep is one O(n)
/// reduction for S plus O(grid) multiply-adds - the ResidualDfPrior economics,
/// now an exact analogy. The grid is dense where dispersion matters and sparse
/// toward the Poisson-like cap at r = 50.
struct NBDispersionPrior {
  static constexpr std::size_t gridSize = 13;
  static constexpr double grid[gridSize] = {1.0,  2.0,  3.0,  4.0,  5.0,
                                            6.0,  8.0,  10.0, 12.0, 15.0,
                                            20.0, 30.0, 50.0};
  static constexpr std::size_t medianIndex = 6;  // r = 8, the cold-start value

  NBDispersionPrior() {
    double priorTotal = 0.0;
    for (std::size_t k = 0; k < gridSize; ++k) {
      kernel_[k] = 0.0;
      logPrior_[k] = grid[k] * std::exp(-0.1 * grid[k]);  // gamma(2, 0.1) kernel
      priorTotal += logPrior_[k];
    }
    for (std::size_t k = 0; k < gridSize; ++k)
      logPrior_[k] = std::log(logPrior_[k] / priorTotal);
  }

  /// Precompute the count-dependent lgamma kernel L_k from the response: tally
  /// the integer-count histogram n_c, then L_k = sum_c n_c [lgamma(c + r_k) -
  /// lgamma(r_k)]. Called at construction and whenever y changes (setResponse/
  /// setData), the count analogue of a fixed-data precompute. A non-null active
  /// restricts BOTH passes to the rows in the data set, so the kernel is the
  /// retained subsample's; that makes an active-row mask change a rebuild, at
  /// O(n + maxCount * gridSize), the one real per-install cost of the channel.
  void computeKernel(const double* y, std::size_t numObservations,
                     const double* active = nullptr) {
    std::size_t maxCount = 0;
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (active != nullptr && active[i] == 0.0) continue;
      std::size_t c = static_cast<std::size_t>(std::lround(y[i]));
      if (c > maxCount) maxCount = c;
    }
    std::vector<double> histogram(maxCount + 1, 0.0);
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (active != nullptr && active[i] == 0.0) continue;
      histogram[static_cast<std::size_t>(std::lround(y[i]))] += 1.0;
    }
    for (std::size_t k = 0; k < gridSize; ++k) {
      double lgammaR = std::lgamma(grid[k]);
      double L = 0.0;
      for (std::size_t c = 0; c <= maxCount; ++c)
        if (histogram[c] != 0.0)
          L += histogram[c] *
               (std::lgamma(static_cast<double>(c) + grid[k]) - lgammaR);
      kernel_[k] = L;
    }
  }

  /// Draw a grid index from the discrete full conditional given the collapsed
  /// statistic S = sum_i log(1 - p_i).
  std::size_t drawIndex(ext_rng* rng, double sumLog1mP) {
    for (std::size_t k = 0; k < gridSize; ++k)
      weight_[k] = kernel_[k] + grid[k] * sumLog1mP + logPrior_[k];
    return drawFromLogWeights(rng, weight_.data(), gridSize);
  }

  /// The precomputed L_k, exposed so a component test can check the kernel
  /// against a direct lgamma evaluation.
  double kernelValue(std::size_t k) const { return kernel_[k]; }

private:
  std::array<double, gridSize> kernel_;
  std::array<double, gridSize> logPrior_;
  std::array<double, gridSize> weight_;
};

/// Negative-binomial counts by the Polya-Gamma augmentation (Polson-Scott-Windle
/// 2013; Zhou-Li-Dunson-Carin 2012) under the logit-p parameterization
/// (docs/design/negative-binomial.md): the forest fits the log-odds latent
/// psi_i = f(x_i) + offset_i, the count law is y_i ~ NB(r, plogis(psi_i)) with
/// dispersion r, and E[y_i] = r exp(psi_i) so the offset is a log-exposure. The
/// augmentation generalizes LogisticResponse: omega_i ~ PG(y_i + r, psi_i) (a
/// real shape in general, integer while r is integer-valued), kappa_i =
/// (y_i - r)/2, working response
/// z_i = kappa_i/omega_i - offset_i under per-sweep precisions omega_i, sigma
/// fixed at 1. r is a positive integer, fixed (user-supplied) or estimated on
/// the capped grid by the closed-form discrete full conditional
/// (NBDispersionPrior); real-valued r stays deferred to a future extension
/// (docs/design/negative-binomial.md section 2). Counts
/// enter kappa directly, so like the binary families it does not rescale the
/// response. Weights are unsupported (exposure belongs in the offset).
/// latents() exposes the omega draws.
class NBResponse final : public ResponseModel {
public:
  /// dispersion > 0 fixes r there (an integer; the host validates integrality);
  /// a non-positive dispersion estimates r on the grid, cold-started at its
  /// median. offset may be null; weights are unsupported.
  NBResponse(const double* y, const double* offset,
             std::size_t numObservations, double dispersion)
    : y_(y), offset_(offset), numObservations_(numObservations),
      estimateR_(!(dispersion > 0.0)),
      r_(estimateR_ ? NBDispersionPrior::grid[NBDispersionPrior::medianIndex]
                    : dispersion) {
    omega_.resize(numObservations);
    working_.resize(numObservations);
    rPrior_.computeKernel(y, numObservations);
    coldStart();
  }

  double* workingResponse() override { return working_.data(); }
  /// The Polya-Gamma draws, or a_i omega_i from a separate composite while a
  /// mask is installed; as for logistic, the zero never enters omega_ itself.
  const double* workingWeights() const override {
    return activeRows_.empty() ? omega_.data() : composite_.data();
  }
  bool workingWeightsVaryPerSweep() const override { return true; }
  const double* offset() const override { return offset_; }

  /// Per sweep, in the partially-collapsed order (section 5, van Dyk-Park):
  /// (1) update r from its grid full conditional given the fit, COLLAPSED over
  /// omega (it reads only S = sum_i log(1 - p_i), never the omega draws); then
  /// (2) draw omega_i ~ PG(y_i + r_new, psi_i) at the NEW r; then (3) rebuild the
  /// working response with r_new. The reverse order is not invariant - the trees
  /// would consume an omega whose shape carries the stale r. sigma is ignored
  /// (fixed at 1). Fixed-r mode skips step (1) and draws no dispersion variate.
  /// Under a mask both halves of the r update are the subsample's: S sums the
  /// active rows and the count histogram behind L_k was rebuilt over them.
  void refreshLatents(ext_rng* rng, const double* totalFits, double) override {
    if (estimateR_)
      r_ = NBDispersionPrior::grid[rPrior_.drawIndex(
        rng, collapsedStatistic(totalFits))];
    drawOmega(rng, totalFits);
  }

  bool supportsActiveRows() const override { return true; }

  /// Beyond the logistic composition, the dispersion block is the subsample's:
  /// the count-histogram kernel is REBUILT over the active rows here, which is
  /// the channel's one per-install cost.
  bool setActiveRows(const double* active) override {
    if (active == nullptr) {
      activeRows_.clear();
    } else {
      activeRows_.assign(active, active + numObservations_);
      composite_.resize(numObservations_);
      recompose();
    }
    rPrior_.computeKernel(y_, numObservations_, activePointer());
    return true;
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  /// Embedded-Gibbs y swap: keep the slow-moving r (a global the outer sampler
  /// persists across a small y perturbation - the kept-cutpoints/nu clause),
  /// recompute the count-histogram kernel under the new y, and redraw omega and
  /// working against the current fit.
  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool, double*) override {
    y_ = y;
    rPrior_.computeKernel(y, numObservations_, activePointer());
    drawOmega(rng, totalFits);
  }

  void setOffset(const double* offset, bool, double*) override {
    // omega and kappa / omega stand; only the shift into the working response
    // moves (the LogisticResponse setOffset)
    reshiftWorkingForOffset(working_.data(), offset_, offset, numObservations_);
    offset_ = offset;
  }

  /// Data swap: everything stale, so cold-init r to the grid median (or the
  /// held fixed value), rebuild the kernel, and cold-start omega at its
  /// PG(y+r, 0) mean (y_i + r)/4 (the LogisticResponse w/4 generalization); the
  /// first sweep's draw replaces it.
  void setData(const double* y, const double* offset, const double*,
               std::size_t numObservations, double*) override {
    y_ = y;
    offset_ = offset;
    numObservations_ = numObservations;
    activeRows_.clear();  // length-n and n may have changed
    composite_.clear();
    omega_.resize(numObservations);
    working_.resize(numObservations);
    if (estimateR_)
      r_ = NBDispersionPrior::grid[NBDispersionPrior::medianIndex];
    rPrior_.computeKernel(y, numObservations);
    coldStart();
  }

  const double* latents() const override { return omega_.data(); }

  /// Rebuild the working response from the restored omega AND the current r.
  /// RESTORE CONTRACT (section 5): restoreDispersion MUST run before this, since
  /// the rebuild reads r; a restore that installs omega before r rebuilds
  /// against the stale r.
  void restoreLatents(const double* latents) override {
    std::memcpy(omega_.data(), latents, numObservations_ * sizeof(double));
    for (std::size_t i = 0; i < numObservations_; ++i)
      working_[i] = 0.5 * (y_[i] - r_) / omega_[i] -
                    (offset_ != nullptr ? offset_[i] : 0.0);
    recompose();
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

  /// The dispersion r for the by-name "dispersion" state block: getState reads
  /// dispersion() when carriesDispersion(), stateIsValid refuses a non-finite/
  /// non-positive r, and setState restoreDispersion()s it BEFORE restoreLatents
  /// (the restore contract above). estimatesDispersion() is the grid-vs-fixed
  /// flag, the estimatesResidualDf analog.
  bool carriesDispersion() const override { return true; }
  double dispersion() const override { return r_; }
  bool estimatesDispersion() const { return estimateR_; }
  void restoreDispersion(double dispersion) override { r_ = dispersion; }

  /// The dispersion kernel L_k currently installed and the collapsed statistic
  /// S the grid draw reads, so a component test can check a masked rebuild and
  /// a masked reduction against the compacted response's own.
  double dispersionKernelForTesting(std::size_t k) const {
    return rPrior_.kernelValue(k);
  }
  double collapsedStatisticForTesting(const double* totalFits) const {
    return collapsedStatistic(totalFits);
  }

  /// log dnbinom(y_i; r, plogis(eta_i)) with eta the log-odds f(x) + offset:
  ///   lgamma(y+r) - lgamma(r) - lgamma(y+1) + y log p + r log(1 - p),
  /// using the stable log p = -log(1 + e^{-eta}), log(1 - p) = -log(1 + e^{eta}).
  void computeLogLikelihood(const double* totalFits, double,
                            std::size_t numObservations,
                            double* out) const override {
    double lgammaR = std::lgamma(r_);
    for (std::size_t i = 0; i < numObservations; ++i) {
      if (!isActive(i)) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      double eta = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      double y = y_[i];
      out[i] = std::lgamma(y + r_) - lgammaR - std::lgamma(y + 1.0) -
               y * logOnePlusExp(-eta) - r_ * logOnePlusExp(eta);
    }
  }

private:
  bool isActive(std::size_t i) const {
    return activeRows_.empty() || activeRows_[i] != 0.0;
  }
  const double* activePointer() const {
    return activeRows_.empty() ? nullptr : activeRows_.data();
  }

  /// S = sum_i log(1 - p_i) = -sum_i log(1 + e^psi_i), the collapsed statistic
  /// the dispersion grid draw reads, over the ACTIVE rows only.
  double collapsedStatistic(const double* totalFits) const {
    double sumLog1mP = 0.0;
    for (std::size_t i = 0; i < numObservations_; ++i) {
      if (!isActive(i)) continue;
      sumLog1mP -=
        logOnePlusExp(totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0));
    }
    return sumLog1mP;
  }

  /// c_i = a_i omega_i, the served precisions while a mask is installed.
  void recompose() {
    if (activeRows_.empty()) return;
    for (std::size_t i = 0; i < numObservations_; ++i)
      composite_[i] = activeRows_[i] * omega_[i];
  }

  /// Steps (2)-(3) of the sweep at the CURRENT r: draw omega_i ~ PG(y_i + r,
  /// psi_i) and rebuild the working response. Shared by refreshLatents (after
  /// its r step) and setResponse (which keeps r, the ordinal drawLatents
  /// analog). An inactive row's draw is skipped for logistic's reason.
  void drawOmega(ext_rng* rng, const double* totalFits) {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      if (!isActive(i)) continue;
      double offset = offset_ != nullptr ? offset_[i] : 0.0;
      double psi = totalFits[i] + offset;
      double omega = simulatePolyaGammaShape(rng, y_[i] + r_, psi);
      omega_[i] = omega;
      working_[i] = 0.5 * (y_[i] - r_) / omega - offset;
    }
    recompose();
  }

  /// Deterministic cold start: omega at PG(y+r, 0)'s mean (y_i + r)/4, so the
  /// working response starts at 2 (y_i - r)/(y_i + r) - offset independent of the
  /// fit; real draws replace it after the first sweep.
  void coldStart() {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double omega = 0.25 * (y_[i] + r_);
      omega_[i] = omega;
      working_[i] =
        0.5 * (y_[i] - r_) / omega - (offset_ != nullptr ? offset_[i] : 0.0);
    }
  }

  const double* y_;
  const double* offset_;
  std::size_t numObservations_;
  bool estimateR_;
  double r_;
  std::vector<double> activeRows_;  // the 0/1 mask; empty when none
  std::vector<double> composite_;   // c_i = a_i omega_i, served while masked
  std::vector<double> omega_;
  std::vector<double> working_;
  NBDispersionPrior rPrior_;
};

/// Tau priors the in-core grouped sampler supports (rbart.priors); custom
/// R-function priors keep rbart_vi's R loop instead.
enum class TauPriorKind { cauchy, gamma };

/// log p(tau) for the built-in priors, matching rbart.priors evaluated with
/// scale already converted to the engine's internal response scale:
/// dcauchy(tau; 0, scale) or dgamma(tau; shape 2.5, scale), both log.
/// Normalizing constants are kept so component tests compare against R
/// directly; both are scale families, so the internal-scale conversion is
/// exactly a division of the scale parameter.
inline double logTauPrior(TauPriorKind kind, double scale, double tau) {
  if (kind == TauPriorKind::cauchy) {
    double z = tau / scale;
    return -std::log(std::numbers::pi * scale) - std::log1p(z * z);
  }
  return 1.5 * std::log(tau) - tau / scale - 2.5 * std::log(scale) -
         std::lgamma(2.5);
}

/// The R loop's tau log posterior given the group effects:
///   -J log(tau) - 0.5 sum(b^2) / tau^2 + log p(tau)
/// on (0, Inf); anything outside scores -inf.
inline double logTauPosterior(double tau, double numGroups,
                              double sumSquaredEffects, TauPriorKind kind,
                              double priorScale) {
  if (tau <= 0.0 || std::isinf(tau)) return -HUGE_VAL;
  return -numGroups * std::log(tau) -
         0.5 * sumSquaredEffects / (tau * tau) +
         logTauPrior(kind, priorScale, tau);
}

/// One slice-sampling step (Neal 2003: stepping out, then shrinkage) of a
/// log density on (lower, upper); width is the stepping-out unit and the
/// interval endpoints clamp to the boundary, as R's sliceSample does.
/// Stepping out is bounded by Neal's m at 10^4 expansions per side: healthy
/// runs measure under 100 expansions, while an unbounded search chasing a
/// distant mode costs ~mode/width density evaluations (an effective hang
/// once tau has random-walked past ~1e8, reachable through mostly-empty
/// groupings). A capped bracket still contains the current point, so
/// shrinkage samples correctly inside it. Shrinkage terminates in exact
/// arithmetic; its iteration cap turns numeric pathology into keeping the
/// current point.
template <typename LogDensity>
double sliceSampleOnce(ext_rng* rng, const LogDensity& logDensity, double x,
                       double width, double lower, double upper) {
  constexpr int maxStepOut = 10000;
  double logHeight = logDensity(x) - ext_rng_simulateExponential(rng, 1.0);
  double u = ext_rng_simulateContinuousUniform(rng);
  double left = x - u * width;
  double right = left + width;
  int steps = maxStepOut;
  while (steps-- > 0 && left > lower && logDensity(left) > logHeight)
    left -= width;
  if (left < lower) left = lower;
  steps = maxStepOut;
  while (steps-- > 0 && right < upper && logDensity(right) > logHeight)
    right += width;
  if (right > upper) right = upper;

  for (int iteration = 0; iteration < 1000; ++iteration) {
    double proposal =
      left + ext_rng_simulateContinuousUniform(rng) * (right - left);
    if (logDensity(proposal) > logHeight) return proposal;
    if (proposal < x) left = proposal; else right = proposal;
  }
  return x;
}

/// Exact conjugate draw of the grouped scale tau under the half-Cauchy prior
/// tau ~ C+(0, A), via the Makalic-Schmidt (2016) inverse-gamma scale mixture:
///   tau^2 | xi ~ IG(1/2, 1/xi),   xi ~ IG(1/2, 1/A^2)   ==>   tau ~ C+(0, A).
/// With b_j ~ N(0, tau^2), SS = sum b_j^2, J groups, the two full conditionals
/// are exact inverse gammas
///   xi    | tau    ~ IG(1,         1/tau^2 + 1/A^2)
///   tau^2 | b, xi  ~ IG((J + 1)/2, 0.5 SS + 1/xi),
/// drawn as one reciprocal gamma each. This is the drawGlue idiom (chain.hpp):
/// ext_rng_simulateGamma takes a SCALE, so an IG(shape, rate) draw is
/// 1 / Gamma(shape, scale = 1/rate). The update therefore consumes EXACTLY two
/// ext_rng_simulateGamma draws per sweep - a fixed count, unlike the slice
/// sampler's data-dependent step-out/shrinkage; it cannot hang and needs no
/// step-out cap. The auxiliary xi is conditionally independent of everything
/// given tau, so it is redrawn fresh from its conditional against the current
/// tau each sweep and never persisted (no state-format change). tauCurrent is
/// the previous draw the xi conditional conditions on.
inline double drawTauCauchyExactIG(ext_rng* rng, double numGroups,
                                   double sumSquaredEffects, double priorScale,
                                   double tauCurrent) {
  double invASq = 1.0 / (priorScale * priorScale);
  double xi = 1.0 / ext_rng_simulateGamma(
                      rng, 1.0, 1.0 / (1.0 / (tauCurrent * tauCurrent) + invASq));
  double tauSq =
    1.0 / ext_rng_simulateGamma(rng, 0.5 * (numGroups + 1.0),
                                1.0 / (0.5 * sumSquaredEffects + 1.0 / xi));
  return std::sqrt(tauSq);
}

/// One conjugate draw of grouped random intercepts: with working response
/// z_i ~ N(F_i + b_g(i), sigma^2 / w_i) and prior b_j ~ N(0, tau^2), the
/// groups are independent normals
///   prec_j = (sum_{i in j} w_i) / sigma^2 + 1 / tau^2
///   mean_j = (sum_{i in j} w_i (z_i - F_i)) / (sigma^2 prec_j).
/// Null weights are unit weights; a group with no observations draws from
/// its prior through the same formula. Draws land in effects (also used as
/// the residual accumulator) in group order, one standard normal each;
/// weightScratch supplies numGroups doubles of scratch.
inline void drawGroupEffects(ext_rng* rng, const double* z,
                             const double* weights, const double* totalFits,
                             const std::uint32_t* groupIndex,
                             std::size_t numObservations,
                             std::size_t numGroups, double sigma, double tau,
                             double* effects, double* weightScratch) {
  for (std::size_t j = 0; j < numGroups; ++j) {
    effects[j] = 0.0;
    weightScratch[j] = 0.0;
  }
  for (std::size_t i = 0; i < numObservations; ++i) {
    double w = weights == nullptr ? 1.0 : weights[i];
    std::uint32_t j = groupIndex[i];
    weightScratch[j] += w;
    effects[j] += w * (z[i] - totalFits[i]);
  }
  double sigmaSq = sigma * sigma;
  for (std::size_t j = 0; j < numGroups; ++j) {
    double precision = weightScratch[j] / sigmaSq + 1.0 / (tau * tau);
    double mean = effects[j] / (sigmaSq * precision);
    effects[j] =
      mean + ext_rng_simulateStandardNormal(rng) / std::sqrt(precision);
  }
}

/// Grouped random intercepts as a response-model decorator, running
/// rbart_vi's Gibbs blocks in-core over any base family: the trees condition
/// on workingResponse() = base z - b_g(i), so recorded fits stay f(x)-only
/// and the user offset and gaussian scale anchoring belong to the base model
/// untouched. b and tau live on the base model's internal scale; the tau
/// prior's relative scale converts once at construction (both priors are
/// scale families) and reported draws de-scale by sigmaScale(), like sigma.
/// The weighted conjugate update deliberately replaces the R loop's
/// unweighted group means, which is what makes Polya-Gamma logistic weights
/// compose. Group structure is fixed at creation, so the bridge refuses
/// setData; a same-length setResponse or setOffset delegates faithfully and is
/// allowed, except at updateScale = TRUE under a re-anchoring base family,
/// where b, tau and the tau prior scale would silently restate.
class GroupedResponse final : public ResponseModel {
public:
  /// tauPriorScale is the original-scale relative scale (sd(y) continuous,
  /// 0.5 binary); the prior's 2.5 multiplier and the internal-scale
  /// conversion are applied here. Initialization matches the R loop: tau at
  /// a fifth of the relative scale, b drawn from its prior.
  GroupedResponse(std::unique_ptr<ResponseModel> base,
                  std::size_t numObservations,
                  const std::uint32_t* groupIndices, std::size_t numGroups,
                  TauPriorKind tauPriorKind, double tauPriorScale,
                  std::size_t tauSliceSteps, ext_rng* rng)
    : base_(std::move(base)), numObservations_(numObservations),
      groupIndex_(groupIndices, groupIndices + numObservations),
      priorKind_(tauPriorKind), sliceSteps_(tauSliceSteps) {
    double scale = base_->sigmaScale();
    priorScale_ = 2.5 * tauPriorScale / scale;
    tau_ = (tauPriorScale / 5.0) / scale;
    groupEffects_.resize(numGroups);
    for (std::size_t j = 0; j < numGroups; ++j)
      groupEffects_[j] = tau_ * ext_rng_simulateStandardNormal(rng);
    weightScratch_.resize(numGroups);
    workingResponse_.resize(numObservations);
    fitScratch_.resize(numObservations);
    rebuildWorking();
  }

  double* workingResponse() override { return workingResponse_.data(); }
  const double* workingWeights() const override {
    return base_->workingWeights();
  }
  bool workingWeightsVaryPerSweep() const override {
    return base_->workingWeightsVaryPerSweep();
  }
  const double* offset() const override { return base_->offset(); }

  /// rbart_vi's blocks in their order relative to the tree sweep: the base
  /// refreshes its latents against f + b, then b draws conjugately against
  /// the fresh working response and weights, then tau by slice sampling,
  /// keeping the last of sliceSteps_ steps as the R loop does. sigma is the
  /// previous draw, exactly the value the R loop's b update conditions on.
  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double sigma) override {
    shiftFits(totalFits);
    base_->refreshLatents(rng, fitScratch_.data(), sigma);
    drawGroupEffects(rng, base_->workingResponse(), base_->workingWeights(),
                     totalFits, groupIndex_.data(), numObservations_,
                     groupEffects_.size(), sigma, tau_, groupEffects_.data(),
                     weightScratch_.data());

    double sumSquaredEffects = 0.0;
    for (double effect : groupEffects_)
      sumSquaredEffects += effect * effect;
    double numGroups = static_cast<double>(groupEffects_.size());
    if (priorKind_ == TauPriorKind::cauchy) {
      // Exact Makalic-Schmidt two-block conjugate draw: a fixed 2 gamma draws
      // per sweep. The cauchy branch no longer consumes sliceSteps_ or the
      // slice's step-out cap (both stay live for the gamma branch below, which
      // admits no exact draw as parameterized).
      tau_ = drawTauCauchyExactIG(rng, numGroups, sumSquaredEffects,
                                  priorScale_, tau_);
    } else {
      TauPriorKind kind = priorKind_;
      double priorScale = priorScale_;
      auto logDensity = [=](double tau) {
        return logTauPosterior(tau, numGroups, sumSquaredEffects, kind,
                               priorScale);
      };
      for (std::size_t step = 0; step < sliceSteps_; ++step)
        tau_ = sliceSampleOnce(rng, logDensity, tau_, priorScale_, 0.0,
                               HUGE_VAL);
    }

    rebuildWorking();
  }

  /// sigma conditions on f + b through the same shifted scratch.
  double drawSigma(ext_rng* rng, const double* totalFits,
                   double sigma) override {
    shiftFits(totalFits);
    return base_->drawSigma(rng, fitScratch_.data(), sigma);
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool updateScale, double* sigmaInOut) override {
    shiftFits(totalFits);
    base_->setResponse(y, rng, fitScratch_.data(), updateScale, sigmaInOut);
    rebuildWorking();
  }

  /// The base handles any scale re-anchoring; b and tau keep their values
  /// against the new transform, so an updateScale offset change would shift
  /// what they mean on the original scale - which is why the bridge refuses
  /// exactly that pair. rbart_vi's in-core path never calls this.
  void setOffset(const double* offset, bool updateScale,
                 double* sigmaInOut) override {
    base_->setOffset(offset, updateScale, sigmaInOut);
    rebuildWorking();
  }

  /// Group indices are per-observation and fixed at creation, so only a
  /// same-length replacement is coherent; the bridge refuses this whole-data
  /// call outright, at every family and every scale.
  void setData(const double* y, const double* offset, const double* weights,
               std::size_t numObservations, double* sigmaInOut) override {
    base_->setData(y, offset, weights, numObservations, sigmaInOut);
    rebuildWorking();
  }

  /// The weight swap alone, exactly as setResponse delegates the response one:
  /// the base takes it against f + b and refreshes whatever its own
  /// augmentation states against the counts, while b and tau - blocks that
  /// belong to the sweep - do not draw here.
  void setWeights(const double* weights, ext_rng* rng,
                  const double* totalFits) override {
    shiftFits(totalFits);
    base_->setWeights(weights, rng, fitScratch_.data());
    rebuildWorking();
  }

  /// Pure delegation: drawGroupEffects already weights its per-group sums by
  /// workingWeights(), so an inactive row leaves its group's mean and
  /// precision and an all-inactive group falls back to its prior through the
  /// same formula.
  bool supportsActiveRows() const override {
    return base_->supportsActiveRows();
  }
  bool setActiveRows(const double* active) override {
    return base_->setActiveRows(active);
  }

  void setSigmaPrior(double sigmaEstimate, double degreesOfFreedom,
                     double rawScale) override {
    base_->setSigmaPrior(sigmaEstimate, degreesOfFreedom, rawScale);
  }

  const double* latents() const override { return base_->latents(); }

  void restoreLatents(const double* latents) override {
    base_->restoreLatents(latents);
    rebuildWorking();
  }

  // forward the t-residual channel so a grouped Student-t response serializes
  // its nu alongside the lambda latents above
  bool carriesResidualDf() const override { return base_->carriesResidualDf(); }
  double residualDf() const override { return base_->residualDf(); }
  void restoreResidualDf(double nu) override { base_->restoreResidualDf(nu); }

  void getScale(double& min, double& max) const override {
    base_->getScale(min, max);
  }
  void restoreScale(double min, double max) override {
    base_->restoreScale(min, max);
    rebuildWorking();
  }

  double initialSigma() const override { return base_->initialSigma(); }
  double fitScale() const override { return base_->fitScale(); }
  double fitShift() const override { return base_->fitShift(); }
  double sigmaScale() const override { return base_->sigmaScale(); }
  double sigmaDegreesOfFreedomForTesting() const override {
    return base_->sigmaDegreesOfFreedomForTesting();
  }

  /// The per-observation location is f(x_i) + b_g(i); add the group intercept
  /// on the base's internal scale (as shiftFits does) and defer to the base
  /// family's density, so the log-likelihood conditions on the drawn effects.
  void computeLogLikelihood(const double* totalFits, double sigma,
                            std::size_t numObservations,
                            double* out) const override {
    std::vector<double> shifted(numObservations);
    for (std::size_t i = 0; i < numObservations; ++i)
      shifted[i] = totalFits[i] + groupEffects_[groupIndex_[i]];
    base_->computeLogLikelihood(shifted.data(), sigma, numObservations, out);
  }

  const double* groupEffects() const override {
    return groupEffects_.data();
  }
  std::size_t numGroupEffects() const override {
    return groupEffects_.size();
  }
  double groupTau() const override { return tau_; }
  void restoreGroupEffects(const double* effects, double tau) override {
    std::memcpy(groupEffects_.data(), effects,
                groupEffects_.size() * sizeof(double));
    tau_ = tau;
    rebuildWorking();
  }

private:
  void shiftFits(const double* totalFits) {
    for (std::size_t i = 0; i < numObservations_; ++i)
      fitScratch_[i] = totalFits[i] + groupEffects_[groupIndex_[i]];
  }

  void rebuildWorking() {
    const double* z = base_->workingResponse();
    for (std::size_t i = 0; i < numObservations_; ++i)
      workingResponse_[i] = z[i] - groupEffects_[groupIndex_[i]];
  }

  std::unique_ptr<ResponseModel> base_;
  std::size_t numObservations_;
  std::vector<std::uint32_t> groupIndex_;  // per observation, 0..J-1
  std::vector<double> groupEffects_;       // b, internal scale
  double tau_ = 1.0;                       // internal scale
  TauPriorKind priorKind_;
  double priorScale_ = 1.0;                // internal scale, 2.5 applied
  std::size_t sliceSteps_ = 1;
  std::vector<double> workingResponse_;    // base z minus b_g(i)
  std::vector<double> fitScratch_;         // f + b for base updates
  std::vector<double> weightScratch_;      // per-group weight sums
};

}  // namespace bartcore

#endif  // BARTCORE_MODEL_HPP
