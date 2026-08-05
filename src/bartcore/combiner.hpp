#ifndef BARTCORE_COMBINER_HPP
#define BARTCORE_COMBINER_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include <external/random.h>
#include <misc/linearAlgebra.h>

#include "data.hpp"
#include "model.hpp"
#include "moves.hpp"
#include "tree.hpp"

// The multi-forest coupling: the per-forest ensemble Forest<L, ResidT>, the (response,
// precision) pair a combiner forms per forest (ForestResponse), the
// ForestCombiner<L> base a multi-forest Chain delegates to, and its first
// instance BCFForestCombiner<L>. The serializable per-forest and per-chain
// state (ForestStateData, ChainStateData) and BCF's spec/glue structs live here
// too, so a combiner owns its wire format beside its logic.

namespace bartcore {

/// One forest's serializable tree channels: value-encoded flattened live and
/// saved trees, their slope/mask side channels, and the per-forest k. A
/// ChainStateData holds one-or-more (BCF's prognostic and treatment forests);
/// a single-forest state is a length-1 forest vector (docs/design/bcf.md).
struct ForestStateData {
  std::vector<std::vector<FlatNode>> trees;
  std::vector<std::vector<FlatNode>> savedTrees;  // slot-major; empty unless kept
  // vector-parameter leaves only: each tree's slopes, numParams - 1 doubles
  // per leaf in pre-order (the intercept rides in the flattened value);
  // savedTreeParams parallels savedTrees. Empty for scalar leaf models.
  std::vector<std::vector<double>> treeParams;
  std::vector<std::vector<double>> savedTreeParams;
  // wide categorical columns only: each flat tree's mask side channel,
  // parallel to trees/savedTrees. Empty otherwise.
  std::vector<std::vector<std::uint64_t>> treeMasks;
  std::vector<std::vector<std::uint64_t>> savedTreeMasks;
  double k = 2.0;
};

/// Everything a chain's posterior state comprises, in host-exchangeable form:
/// one-or-more forests' trees, plus the chain-shared sigma (original scale),
/// response latents, DART state, the serialized rng, and BCF's glue scalars.
/// Restore rebuilds the rest canonically - partitions from the tree structure
/// and cut points, totalFits by summing the tree fits, the variance prior by
/// re-anchoring through the transform - so a restored chain continues
/// equivalently, not bitwise: the last ulp of the dropped accumulation history
/// is not reproduced.
struct ChainStateData {
  std::vector<ForestStateData> forests;
  double sigma = 1.0;  // original response scale
  // the gaussian response transform at capture; max <= min marks scale-free
  double fitMin = 0.0, fitMax = 0.0;
  std::vector<double> latents;            // empty for gaussian; lambda under t
  // Student-t continuous errors (TResponse) only: the residual df nu at
  // capture. NaN marks absent, so gaussian and every non-t state carry no nu
  // block and a t sampler refuses a state lacking one.
  double residualDf = std::numeric_limits<double>::quiet_NaN();
  // ordinal (cumulative-probit) responses only: the length-(K-1) cutpoint
  // vector at capture. Empty marks absent, so every non-ordinal state carries
  // no cutpoint block and an ordinal sampler refuses a state lacking one.
  std::vector<double> cutpoints;
  // negative-binomial counts (NBResponse) only: the dispersion r at capture
  // (docs/design/negative-binomial.md). NaN marks absent, so every non-NB state
  // carries no dispersion block and an NB sampler refuses a state lacking one.
  // omega rides the latents block above; this is its companion scalar.
  double dispersion = std::numeric_limits<double>::quiet_NaN();
  // grouped samplers only, internal scale so restores are exact
  std::vector<double> groupEffects;
  double groupTau = 0.0;
  std::vector<double> dartProbabilities;  // empty when DART is off
  double dartAlpha = 1.0;
  size_t dartNumUpdatesSkipped = 0;
  std::vector<unsigned char> rngState;
  // BCF combining response's glue (docs/design/bcf.md); false off BCF
  bool hasBCF = false;
  double a = 1.0, aVariance = 1.0, b0 = 0.0, b1 = 1.0;
  // heteroscedastic variance forest (docs/design/heteroscedastic.md): the
  // variance trees' flattened structure, each leaf a POSITIVE scale value on
  // the working scale. Empty for every homoscedastic state; a NEW branch of
  // the flat format (the value semantics and validation differ from a Gaussian
  // leaf - it must be positive). Rebuild recomputes s^2(x) from the product.
  std::vector<std::vector<FlatNode>> varianceTrees;
};

/// One ensemble of the backfitting sampler: its trees, their fits and working
/// residual, the leaf model, the split selector, and the per-forest prior and
/// options. A Chain holds one-or-more (BCF combines a prognostic and a
/// treatment forest); everything a Forest omits - the rng, the response model,
/// sigma, and the shared store - is chain-level and passed in by the methods
/// that drive it.
template <IntegrableLeafModel L, typename ResidT = double>
struct Forest {
  static_assert(std::is_same_v<ResidT, double> || std::is_same_v<ResidT, float>,
                "the residual storage type is fp64 (default) or opt-in fp32");
  // per-forest options: tree count, tree-move probabilities, the k
  // hyperprior, and whether the split selector is DART
  size_t numTrees = 200;
  double birthOrDeathProbability = 0.5;
  double swapProbability = 0.1;
  double changeProbability = 0.4;
  double birthProbability = 0.5;
  bool updateK = false;
  ChiKHyperprior kHyperprior;
  bool useDart = false;

  L leaf;
  CGMTreePrior treePrior;
  DartPrior dart;
  std::vector<double> fixedSplitProbabilities;
  std::vector<std::uint32_t> splitCounts;
  // per-forest split-variable restriction materialized as a 0/1 byte per
  // predictor (1 = splittable); empty leaves every column available and the
  // trees carry a null mask, so the default availability path is unchanged.
  std::vector<std::uint8_t> columnMask;
  // per-forest interaction constraint (max-order + forbidden co-occurrence);
  // inactive by default (build not called), leaving the trees' constraint null
  // and the availability path byte-for-byte unchanged
  // (docs/design/interaction-constraints.md, the per-path max-order +
  // co-occurrence mechanism). HEAP-allocated so
  // the address a tree's interaction_ borrows survives a forests_ vector
  // reallocation/move; a value member relocates with the Forest and dangles
  // every tree's pointer (the BCF multi-forest heap-use-after-free), whereas a
  // heap pointee stays put across the move, exactly as columnMask's data() does.
  std::unique_ptr<InteractionConstraint> interaction;
  // per-tree block-additive constraint (docs/design/interaction-constraints.md,
  // variant A): confine each WHOLE tree to one declared group of predictors so the
  // ensemble is exactly f = sum_G f_G. blockMasks holds numBlocks group-membership
  // rows of numPredictors 0/1 bytes, group-major (row g intersected with any base
  // columnMask - the BCF tau case, where each tree's block row is ANDed with its
  // moderator restriction); blockOfTree[t] is tree t's 0-based group.
  // Each tree's columnMask_ points at its group's row. Both are PLAIN VALUE
  // vectors, NOT unique_ptr: a std::vector's data() heap pointer survives a Forest
  // move (exactly like columnMask above), so the trees' borrowed pointers stay
  // valid - unlike the interaction OBJECT address, which relocates with the Forest.
  // Empty leaves every tree unrestricted, byte-for-byte the default path.
  std::vector<std::uint8_t> blockMasks;
  std::vector<std::size_t> blockOfTree;

  // k is fixed unless updateK; the two accumulators gather the leaf sum of
  // squares and count over a sweep, feeding the k hyperprior draw
  double k = 2.0;
  double kSumSquaredParams = 0.0;
  double kNumLeaves = 0.0;

  std::vector<Tree> trees;
  std::vector<index_t> indexBuffer;
  // dense per-tree fit slab (numObservations x numTrees, tree-major); empty for
  // the constant leaf, which carries muByTree + leafOf instead
  std::vector<double> treeFits;
  std::vector<double> totalFits, totalTestFits;
  // the running per-tree residual; fp64 by default, fp32 under the opt-in
  // storage axis (docs/design/reduced-precision-storage.md sec 3b)
  std::vector<ResidT> treeY;
  std::vector<double> currTestFits;
  std::vector<double> paramByNode;
  // constant leaf only: each live tree's node-indexed leaf values, and the
  // per-observation bottom-node map (tree-major, numObservations x numTrees),
  // so tree t's fit for observation i is muByTree[t][leafOf[t * n + i]]. Empty
  // for vector and function leaves, which keep the dense treeFits slab.
  // leafOf is patched incrementally on accepted moves; a wholesale structure
  // change that leaves fits stale (sampleTreesFromPrior) marks its trees in
  // leafOfStale instead, and the sweep rebuilds before the tree's draw.
  std::vector<std::vector<double>> muByTree;
  std::vector<std::uint32_t> leafOf;
  std::vector<std::uint8_t> leafOfStale;
  // vector-parameter leaves: each live tree's parameter blocks, arena-
  // indexed with stride numParams and kept consistent with the tree's
  // structure between sweeps (fits are no longer constant per leaf, so
  // parameters persist instead of being recovered); savedTreeParams holds
  // the slopes of each saved tree, parallel to savedTrees. Both stay empty
  // for scalar leaf models.
  std::vector<std::vector<double>> paramsByTree;
  MoveScratch scratch;

  // Saved-tree (keepTrees) storage: a circular buffer of capacity slots,
  // each one kept sample's forest in flattened form (slot-major, capacity x
  // numTrees). The slot base is set by the sampler before every run so
  // chains write consistent slots without sharing mutable state.
  std::vector<std::vector<FlatNode>> savedTrees;
  std::vector<std::vector<double>> savedTreeParams;
  // pooled categorical columns (> 63 levels): each saved tree's flattened
  // mask words, parallel to savedTrees; empty otherwise
  std::vector<std::vector<std::uint64_t>> savedTreeMasks;
  size_t savedTreeCapacity = 0;
  size_t savedSlotBase = 0;
};

/// Per-forest calibration for a BCF sampler: the prognostic (mu) and
/// treatment (tau) forests carry bcf's distinct tree counts and priors
/// (docs/design/bcf.md). Node scales are not spec'd: the calibration map
/// derives them from the response sd at construction, and the adaptive
/// magnitude lives in the glue (a for mu, b0/b1 for tau).
struct BCFForestSpec {
  std::size_t numTrees = 200;
  double base = 0.95, power = 2.0;
  double birthOrDeathProbability = 0.5, swapProbability = 0.1,
         changeProbability = 0.4, birthProbability = 0.5;
  // optional split-variable restriction (borrowed 0-based column indices,
  // consumed at construction): the columns this forest may split on. Null or
  // count 0 leaves every column available - the default, byte-for-byte
  // unchanged. BCF's treatment forest carries its moderator subset here; the
  // prognostic forest leaves it empty and reads the full store.
  const std::size_t* columns = nullptr;
  std::size_t numColumns = 0;
  // optional per-forest interaction constraint (docs/design/interaction-
  // constraints.md): a max-order cap (0 = uncapped) and a flat 0-based (a, b)
  // forbidden co-occurrence pair stream (2 per pair; borrowed, consumed at
  // construction). BCF's mu / tau forests carry independent constraints - the
  // calibrated-additivity causal use (an additive-or-low-order tau, a free mu).
  // Both defaults leave the availability path byte-for-byte unchanged.
  std::size_t interactionMaxOrder = 0;
  const std::size_t* interactionForbiddenPairs = nullptr;
  std::size_t interactionNumForbiddenPairs = 0;
  // optional per-forest block-additive constraint (variant A): mu / tau carry
  // independent partitions. blockOfColumn (length numPredictors, borrowed) gives
  // each column's 0-based group (negative = in no block); blockTreeCounts (length
  // numBlocks, borrowed) is the deterministic contiguous per-group tree capacity,
  // summing to numTrees. numBlocks 0 leaves the forest unrestricted; tau's block
  // rows intersect its moderator columnMask at install.
  std::size_t numBlocks = 0;
  const std::int32_t* blockOfColumn = nullptr;
  const std::size_t* blockTreeCounts = nullptr;
};

/// The two model specs plus the treatment vector a BCF chain is built from.
struct BCFSpec {
  BCFForestSpec mu, tau;
  const double* z = nullptr;    // borrowed 0/1 treatment indicator per obs
  double aPriorScale = 2.0;     // half-Cauchy median for the mu scalar a
  double bPriorVariance = 0.5;  // N(0, .) prior variance for b0, b1
  double sdModerate = 1.0;      // treatment effect scale in sd(y) units
  bool updateA = true, updateB = true;  // false fixes the matching glue block
};

/// One category forest's calibration for a multinomial sampler; the K forests
/// are symmetric, so a single spec builds them all (mbart2's convention). Node
/// scale is not spec'd here: the chain constructor sets every forest's leaf
/// scale from nodeScale (the pi*sqrt(3)/sqrt(2) anchor, see
/// docs/design/multinomial.md) and k.
struct MultinomialForestSpec {
  std::size_t numTrees = 200;
  double base = 0.95, power = 2.0;
  double birthOrDeathProbability = 0.5, swapProbability = 0.1,
         changeProbability = 0.4, birthProbability = 0.5;
};

/// The specification a multinomial (softmax) chain is built from: K symmetric
/// category forests over the shared design, the borrowed grouped-count response,
/// and the leaf-scale calibration. The K forests couple through a softmax
/// likelihood with an interleaved one-vs-rest Polya-Gamma augmentation
/// (docs/design/multinomial.md). counts is an n x K nonnegative integer matrix,
/// category-major (column k contiguous, at k*n) to match the combiner's omega_
/// layout; trials is the per-observation trial count n_i = sum_k counts[k*n + i]
/// (>= 1). Single-trial labels enter as a one-hot counts matrix with every
/// trial 1, the exact n_i = 1 reduction the bridge builds.
struct MultinomialSpec {
  std::size_t numCategories = 0;      // K
  const int* counts = nullptr;        // borrowed n x K count matrix, category-major
  const int* trials = nullptr;        // borrowed per-observation trial count n_i
  MultinomialForestSpec forest;       // shared across the K category forests
  // per-forest leaf scale anchor and k: the pairwise log-odds f_ik - f_ij is a
  // difference of two forests, so the per-forest node scale is the logistic
  // pi*sqrt(3) anchor divided by sqrt(2); tau = nodeScale / k is the per-forest
  // total-fit prior sd the level-centering move reads.
  double nodeScale = 3.847649490485592;  // pi*sqrt(3)/sqrt(2)
  double k = 2.0;
};

/// The combining response's glue (docs/design/bcf.md): the prognostic scalar
/// a (half-Cauchy via the scale-mixture auxiliary aVariance) and the
/// treatment scales b0/b1, plus the sweep's per-forest scratch. y = a mu +
/// b_z tau + eps; a Chain holds one only in BCF mode.
struct BCFState {
  const double* z = nullptr;
  double a = 1.0, b0 = 0.0, b1 = 1.0;
  double aVariance = 1.0;
  double aPriorScale = 2.0;
  double bPriorVariance = 0.5;
  bool updateA = true, updateB = true;  // false holds the block at its value
  std::vector<double> combined, forestResponse, forestWeights;
};

/// Forest f's effective response and precision for its own leaf draws: the pair
/// a combiner forms so f's constant-leaf node sums reproduce the residual (y net
/// of the other forests' scaled contributions). Both pointers alias the
/// combiner's per-forest scratch and stay valid only until the next
/// formForestResponse call; a location forest routes into response, a variance
/// forest could later route into weights.
struct ForestResponse {
  const double* response;
  const double* weights;
};

/// The multi-forest coupling a Chain delegates to when it holds more than one
/// forest: per-forest residual formation and the combined per-observation
/// location the response model and sigma draw read. A Chain builds a combiner
/// only in a multi-forest mode (BCF today); a single-forest chain leaves
/// combiner_ null and never pays a virtual call. Templated on the leaf because
/// a combiner's post-combine move reaches Forest<L, ResidT>'s buffers and saved-tree
/// nodes.
///
/// The coupling draw (drawGlue), its post-combine move (afterCombine), the
/// reporting-channel map, and glue (de)serialization are declared here so a
/// subclass owns them without reshaping the base; the base leaves them inert so
/// a combiner that only forms an additive combination need not override them.
template <IntegrableLeafModel L, typename ResidT = double>
struct ForestCombiner {
  virtual ~ForestCombiner() = default;

  /// Forest f's (response, weights) against the residual; forests carries every
  /// forest's current fits, y/w the chain's working response and precisions.
  virtual ForestResponse formForestResponse(std::size_t f,
      const std::vector<Forest<L, ResidT>>& forests, const double* y,
      const double* w) = 0;
  /// The combined per-observation location over all forests; the pointer aliases
  /// combiner scratch, valid only until the next call.
  virtual const double* combinedFits(const std::vector<Forest<L, ResidT>>& forests) = 0;

  /// combinedFits' out-of-sample analog: the combined per-observation TEST
  /// location(s), formed from the forests' totalTestFits, aliasing combiner
  /// scratch valid only until the next call. Base inert (nullptr): it is
  /// reached only for a combiner whose testFitsAreDefined() is true AND whose
  /// numReportedLocations() > 1 (the multinomial softmax blend). BCF leaves
  /// testFitsAreDefined false, and a single-forest chain carries no combiner, so
  /// neither ever calls this - hence the base need not form anything.
  virtual const double* combinedTestFits(const std::vector<Forest<L, ResidT>>&) {
    return nullptr;
  }

  /// Swaps the borrowed treatment vector the coupling reads; inert unless a
  /// subclass carries one. The combiner re-forms its per-forest residuals from
  /// the new vector on the next sweep.
  virtual void setTreatment(const double*) {}

  /// The BCF glue coefficients (a, b0, b1) on the combining response, for the
  /// per-forest reporting path (getBCFGlue); false for a combiner carrying no
  /// such glue, which is how the "no BCF glue" answer reaches the caller.
  virtual bool bcfGlue(double&, double&, double&) const { return false; }

  /// The coupling draw and its likelihood-invariant post-combine move, fired at
  /// the fixed sweep points; inert unless a subclass couples the forests.
  /// afterCombine returns the scale its move applied (1.0 when it makes none) -
  /// the sweep discards it; the component tests read it through the chain.
  virtual void drawGlue(ext_rng*, double, const double*, const double*,
                        const std::vector<Forest<L, ResidT>>&) {}
  virtual double afterCombine(std::vector<Forest<L, ResidT>>&, bool, std::size_t,
                              ext_rng*) { return 1.0; }

  /// A per-forest pre-update hook, fired inside the sweep just before forest f's
  /// tree update, with the partially updated forests (0..f-1 new this sweep,
  /// f..K-1 old). An interleaved coupling draws forest f's latents here against
  /// the CURRENT margins, immediately before formForestResponse(f) reads them;
  /// the base no-op consumes no rng, so every additive combiner (BCF included)
  /// stays bitwise unchanged. Fired in both sweep loops (run and grow-from-root).
  virtual void drawForestGlue(std::size_t, ext_rng*,
                              const std::vector<Forest<L, ResidT>>&) {}

  /// The reporting map: which forest the scalar channels (variable counts, k,
  /// split probabilities) address, and whether the test-fit and log-likelihood
  /// channels are defined. BCF reports forest 0 and leaves those two channels
  /// undefined (no test treatment vector to blend, and no per-observation
  /// location the response model can see to score).
  virtual std::size_t reportedForest() const { return 0; }
  virtual bool testFitsAreDefined() const { return true; }
  virtual bool logLikelihoodIsDefined() const { return true; }

  /// How many per-observation channels combinedFits packs, one per reported
  /// location: 1 for every additive combiner (BCF's a mu + b_z tau is a single
  /// location), K for a model whose combined output is K per-observation values
  /// (softmax probabilities). storeSample loops its training/test writes over
  /// this count, and the run bridge sizes the fits arrays by it; L = 1 leaves
  /// the single-location byte layout every existing path relies on.
  virtual std::size_t numReportedLocations() const { return 1; }

  /// The variable-count reporting set: how many forests the per-sample split-
  /// usage channel records, and which forest slot j addresses. Every additive
  /// combiner reports the single reported forest (count 1, slot 0 =
  /// reportedForest, the exact current channel); a model whose K forests each
  /// carry their own splits (multinomial) reports all K. storeSample loops its
  /// varcount writes over the count, and the run bridge sizes the varcount
  /// array by it; count 1 leaves the single-forest byte layout every existing
  /// path relies on. Distinct from numReportedLocations (softmax output
  /// channels): a future model's fits-location count may diverge from its
  /// forest count, so the varcount axis is keyed on its own forest set.
  virtual std::size_t numVariableCountForests() const { return 1; }
  virtual std::size_t variableCountForest(std::size_t) const {
    return reportedForest();
  }

  /// Whether a whole-response swap (Chain::setResponse at updateScale = false)
  /// is safe on this coupling: true only for a combiner that caches NOTHING
  /// per-forest across sweeps and re-derives every per-forest residual from the
  /// chain's working response on every sweep. Defaults false so a future
  /// multi-forest model stays refused at the bridge until it is audited.
  virtual bool supportsResponseMutation() const { return false; }

  /// Glue (de)serialization into the BCF-shaped state fields; inert unless the
  /// combiner carries glue.
  virtual void serializeGlue(ChainStateData&) const {}
  virtual void restoreGlue(const ChainStateData&) {}
};

/// BCF's combiner (docs/design/bcf.md): a prognostic forest mu (forest 0) and a
/// treatment forest tau (forest 1) combined on a gaussian response as
/// y = a mu + b_z tau + eps. Holds the glue (the scalar a via its half-Cauchy
/// scale mixture aVariance, and the treatment scales b0/b1 over control/treated)
/// and the sweep's per-forest scratch. Constant leaf only, as the whole BCF
/// chain is.
template <IntegrableLeafModel L, typename ResidT = double>
struct BCFForestCombiner : ForestCombiner<L, ResidT> {
  static_assert(!L::hasVectorParams && !L::hasFunctionParams,
                "BCF is a constant-leaf model");

  BCFForestCombiner(const ColumnStore& data, const BCFSpec& spec)
      : data_(data) {
    glue_.z = spec.z;
    glue_.aPriorScale = spec.aPriorScale;
    glue_.bPriorVariance = spec.bPriorVariance;
    glue_.updateA = spec.updateA;
    glue_.updateB = spec.updateB;
  }

  void setTreatment(const double* z) override { glue_.z = z; }

  bool bcfGlue(double& a, double& b0, double& b1) const override {
    a = glue_.a; b0 = glue_.b0; b1 = glue_.b1;
    return true;
  }

  /// Effective response and precision forest f's constant-leaf draws see so that
  /// m_f * f_f explains the residual (y minus the other forests' scaled
  /// contributions): response r_i / m_f, weight w_i m_f^2, which reproduce the
  /// leaf's node sums without touching the leaf math. |m_f| is floored to keep
  /// the division finite in the pathological near-zero-scale case.
  ForestResponse formForestResponse(std::size_t f,
      const std::vector<Forest<L, ResidT>>& forests, const double* y,
      const double* w) override {
    std::size_t n = data_.numObservations;
    glue_.forestResponse.resize(n);
    glue_.forestWeights.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      double resid = y[i];
      for (std::size_t g = 0; g < forests.size(); ++g)
        if (g != f) resid -= forestMultiplier(g, i) * forests[g].totalFits[i];
      double m = forestMultiplier(f, i);
      if (std::fabs(m) < 1.0e-9) m = m < 0.0 ? -1.0e-9 : 1.0e-9;
      glue_.forestResponse[i] = resid / m;
      glue_.forestWeights[i] = (w == nullptr ? 1.0 : w[i]) * m * m;
    }
    return {glue_.forestResponse.data(), glue_.forestWeights.data()};
  }

  const double* combinedFits(const std::vector<Forest<L, ResidT>>& forests) override {
    std::size_t n = data_.numObservations;
    glue_.combined.resize(n);
    const double* mu = forests[0].totalFits.data();
    const double* tau = forests[1].totalFits.data();
    for (std::size_t i = 0; i < n; ++i)
      glue_.combined[i] = glue_.a * mu[i] +
        (glue_.z[i] != 0.0 ? glue_.b1 : glue_.b0) * tau[i];
    return glue_.combined.data();
  }

  /// The glue's Gaussian full conditionals (docs/design/bcf.md): a as the mu
  /// coefficient (prior N(0, aVariance), whose half-Cauchy scale mixture is
  /// refreshed after via an inverse-gamma auxiliary), b0/b1 as the tau
  /// coefficients over control/treated (prior N(0, bPriorVariance)).
  void drawGlue(ext_rng* rng, double sigma, const double* y, const double* w,
                const std::vector<Forest<L, ResidT>>& forests) override {
    std::size_t n = data_.numObservations;
    const double* mu = forests[0].totalFits.data();
    const double* tau = forests[1].totalFits.data();
    double invSigmaSq = 1.0 / (sigma * sigma);

    if (glue_.updateA) {
      double aPrec = 1.0 / glue_.aVariance, aNum = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        double wi = w == nullptr ? 1.0 : w[i];
        double bz = glue_.z[i] != 0.0 ? glue_.b1 : glue_.b0;
        double r = y[i] - bz * tau[i];
        aPrec += wi * mu[i] * mu[i] * invSigmaSq;
        aNum += wi * mu[i] * r * invSigmaSq;
      }
      glue_.a =
        aNum / aPrec + ext_rng_simulateStandardNormal(rng) / std::sqrt(aPrec);

      // t_1 scale mixture: aVariance ~ IG(1/2, scale^2/2) mixes N(0, aVariance)
      // to Cauchy(0, scale), so the conditional's rate carries scale^2, not its
      // inverse
      double rate = 0.5 * glue_.a * glue_.a +
                    0.5 * glue_.aPriorScale * glue_.aPriorScale;
      glue_.aVariance = 1.0 / ext_rng_simulateGamma(rng, 1.0, 1.0 / rate);
    }

    if (glue_.updateB) {
      double bPrec = 1.0 / glue_.bPriorVariance;
      double p0 = bPrec, n0 = 0.0, p1 = bPrec, n1 = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        double wi = w == nullptr ? 1.0 : w[i];
        double r = y[i] - glue_.a * mu[i];
        double prec = wi * tau[i] * tau[i] * invSigmaSq;
        double num = wi * tau[i] * r * invSigmaSq;
        if (glue_.z[i] != 0.0) { p1 += prec; n1 += num; }
        else { p0 += prec; n0 += num; }
      }
      glue_.b0 = n0 / p0 + ext_rng_simulateStandardNormal(rng) / std::sqrt(p0);
      glue_.b1 = n1 / p1 + ext_rng_simulateStandardNormal(rng) / std::sqrt(p1);
    }
  }

  /// Interweaving (ASIS, Yu & Meng 2011) rescale of the prognostic glue ridge.
  /// After the conjugate a draw and the mu leaf draws, jointly rescale the L+1
  /// prognostic-scale coordinates (a, mu_1..mu_L) -> (a/c, c mu_l) along the
  /// likelihood-invariant orbit a mu(x) = (a/c)(c mu(x)), so the move updates
  /// only the amplitude coordinate and preserves the posterior. c = sqrt(v),
  /// v ~ GIG((L-1)/2, M/leafVar, a^2/aVariance) conditioned on the inverse-
  /// gamma auxiliary (exact); L and M are the count and squared sum of the
  /// occupied prognostic leaves. Collapses the slow (a, mu-amplitude) mode
  /// (docs/design/bcf.md). A no-op consuming no rng with a pinned a (updateA
  /// false) or with fewer than two occupied leaves; returns the applied c (1.0
  /// when skipped). record/sampleNum locate the keepTrees saved slot whose mu
  /// leaves, flattened before this move, need the same c so a stored * mu_saved
  /// keeps the identified product.
  double afterCombine(std::vector<Forest<L, ResidT>>& forests, bool record,
                      std::size_t sampleNum, ext_rng* rng) override {
    if (!glue_.updateA) return 1.0;
    Forest<L, ResidT>& forest = forests[0];
    std::size_t n = data_.numObservations;

    // L, M over the occupied prognostic leaves. Recomputed unconditionally:
    // the k-accumulator that would hold these is gated on updateK, which BCF
    // leaves false. A forced-zero empty leaf is not a prior draw, so skip it.
    double M = 0.0;
    std::size_t numLeaves = 0;
    for (std::size_t t = 0; t < forest.numTrees; ++t) {
      Tree& tree = forest.trees[t];
      const std::vector<double>& mu = forest.muByTree[t];
      tree.bottomScratch.clear();
      tree.fillBottom(0, tree.bottomScratch);
      for (int32_t nodeIndex : tree.bottomScratch) {
        const Node& node = tree.at(nodeIndex);
        if (node.numObservations() == 0) continue;
        double value = mu[static_cast<std::size_t>(nodeIndex)];
        M += value * value;
        ++numLeaves;
      }
    }
    if (numLeaves < 2 || !(M > 0.0)) return 1.0;

    // GIG parameters (docs/design/bcf.md): A = M (k/scale)^2, B = a0^2/aVariance
    double a0 = glue_.a;
    double leafPrecision = (forest.k / forest.leaf.scale) *
                           (forest.k / forest.leaf.scale);  // 1 / leafVar
    double gigP = 0.5 * (static_cast<double>(numLeaves) - 1.0);
    double gigA = M * leafPrecision;
    double gigB = a0 * a0 / glue_.aVariance;

    double v = ext_rng_simulateGeneralizedInverseGaussian(rng, gigP, gigA,
                                                          gigB);
    if (!std::isfinite(v) || v <= 0.0) return 1.0;
    double c = std::sqrt(v);
    if (!std::isfinite(c) || c <= 0.0) return 1.0;

    // travel the ridge: a shrinks, the prognostic fits grow by c. Scaling every
    // leaf value scales every gathered fit, so the mu tables carry the rescale.
    glue_.a = a0 / c;
    for (std::size_t t = 0; t < forest.numTrees; ++t)
      misc_scalarMultiplyVectorInPlace(forest.muByTree[t].data(),
                                       forest.muByTree[t].size(), c);
    misc_scalarMultiplyVectorInPlace(forest.totalFits.data(), n, c);
    // aVariance is held: the move conditions on it (ASIS), so refreshing it
    // here re-randomizes the coordinate we just conditioned on and measurably
    // throttles the mixing gain (IACT check, docs/design/bcf.md). The one-sweep
    // lag is benign - the next drawGlue refreshes it | a_new.

    // recorded sweeps carry a live test surface (dead under BCF, but kept
    // self-consistent) and, under keepTrees, this sweep's saved mu slot
    if (record && data_.numTestObservations > 0) {
      misc_scalarMultiplyVectorInPlace(forest.totalTestFits.data(),
                                       data_.numTestObservations, c);
      misc_scalarMultiplyVectorInPlace(forest.currTestFits.data(),
                                       data_.numTestObservations, c);
    }
    if (record && forest.savedTreeCapacity > 0) {
      std::size_t slot =
        (forest.savedSlotBase + sampleNum) % forest.savedTreeCapacity;
      for (std::size_t t = 0; t < forest.numTrees; ++t) {
        std::vector<FlatNode>& flat =
          forest.savedTrees[slot * forest.numTrees + t];
        for (FlatNode& node : flat)
          if (node.variable == invalidVariable) node.value *= c;
      }
    }
    return c;
  }

  /// BCF reports the prognostic forest (forest 0, the base default) but leaves
  /// the test-fit and log-likelihood channels undefined: there is no test
  /// treatment vector to blend a mu + b_z tau off-sample, and the blended
  /// per-observation location is not visible to the response model to score.
  bool testFitsAreDefined() const override { return false; }
  bool logLikelihoodIsDefined() const override { return false; }

  /// BCF admits the scale-pinned response swap, relying on two conditions: the
  /// gaussian response re-maps y through the pinned (min_, range_) and touches
  /// no forest, so both leaf calibrations and the sigma prior stay put; and
  /// this combiner caches nothing per-forest across sweeps - formForestResponse
  /// re-derives every per-forest residual and precision from y each sweep, and
  /// the glue lives here rather than in the response, so it carries over.
  bool supportsResponseMutation() const override { return true; }

  /// The glue scalars into and out of the BCF-shaped wire format. serializeGlue
  /// owns the hasBCF flag (the "carries glue" marker); restoreGlue is a no-op on
  /// a state that carries none, so a mismatched restore leaves the glue at its
  /// constructed values.
  void serializeGlue(ChainStateData& state) const override {
    state.hasBCF = true;
    state.a = glue_.a;
    state.aVariance = glue_.aVariance;
    state.b0 = glue_.b0;
    state.b1 = glue_.b1;
  }
  void restoreGlue(const ChainStateData& state) override {
    if (!state.hasBCF) return;
    glue_.a = state.a;
    glue_.aVariance = state.aVariance;
    glue_.b0 = state.b0;
    glue_.b1 = state.b1;
  }

private:
  /// The scale forest f's constant leaf carries into the combination: a for the
  /// prognostic forest, b_{z_i} for the treatment forest.
  double forestMultiplier(std::size_t f, std::size_t i) const {
    if (f == 0) return glue_.a;
    return glue_.z[i] != 0.0 ? glue_.b1 : glue_.b0;
  }

  const ColumnStore& data_;
  BCFState glue_;
};

/// Softmax over K location-major raw fits (raw[k*n + i] is category k's value
/// for row i) into out in the same layout, log-sum-exp-safe. Shared by the
/// multinomial train/test blends and the K-forest predict replay so the three
/// numerics are one map, not three copies. Safe in place (out == raw): each
/// row's max and normalizer are formed before any of its K entries are
/// overwritten, and distinct rows touch disjoint storage.
inline void softmaxLocationMajor(const double* raw, std::size_t n,
                                 std::size_t K, double* out) {
  for (std::size_t i = 0; i < n; ++i) {
    double maxFit = raw[i];
    for (std::size_t k = 1; k < K; ++k)
      maxFit = std::max(maxFit, raw[k * n + i]);
    double sumExp = 0.0;
    for (std::size_t k = 0; k < K; ++k)
      sumExp += std::exp(raw[k * n + i] - maxFit);
    for (std::size_t k = 0; k < K; ++k)
      out[k * n + i] = std::exp(raw[k * n + i] - maxFit) / sumExp;
  }
}

/// The multinomial (softmax) combiner: K symmetric category forests coupled
/// through a log-linear likelihood, P(y_i = k) = softmax(f_i)_k, with an
/// INTERLEAVED one-vs-rest Polya-Gamma augmentation and a likelihood-invariant
/// level-centering move (docs/design/multinomial.md). Constant leaf only, as
/// the whole chain is.
///
/// Category k's one-vs-rest conditional is a binomial(n_i, .) logistic with
/// linear predictor eta_ik = f_ik - C_ik, where C_ik = log sum_{j != k} exp(f_ij)
/// is the log-sum-exp margin: omega_ik ~ PG(n_i, eta_ik) and forest k sees
/// working response (y_ik - n_i/2)/omega_ik + C_ik under precision omega_ik (the
/// shipped logistic PG machinery, one binomial per category). PG(n_i, .) is the
/// sum of n_i iid PG(1, .) draws, so at n_i = 1 (one-hot single-trial rows) this
/// reduces byte-identically to the label path. The augmentation is INTERLEAVED,
/// not joint: omega_ik is a temporary latent valid only for category k's
/// conditional, so it is drawn against the CURRENT margins in drawForestGlue(k)
/// immediately before forest k's tree update, cycling the categories (Held and
/// Holmes 2006; Polson, Scott and Windle 2013 sec 4). A single post-loop all-K
/// draw would be an invalid Jacobi-style update.
template <IntegrableLeafModel L, typename ResidT = double>
struct MultinomialForestCombiner : ForestCombiner<L, ResidT> {
  static_assert(!L::hasVectorParams && !L::hasFunctionParams,
                "multinomial is a constant-leaf model");

  MultinomialForestCombiner(const ColumnStore& data, const MultinomialSpec& spec)
      : data_(data), numCategories_(spec.numCategories), counts_(spec.counts),
        trials_(spec.trials) {
    std::size_t n = data_.numObservations;
    // n x K omega scratch, cold-started at PG(1, 0)'s mean 1/4 (the logistic
    // seed); drawForestGlue overwrites category k's column each sweep before it
    // is read, so the cold value is only a fallback.
    omega_.assign(n * numCategories_, 0.25);
    margins_.resize(n);
    suffix_.resize(n * numCategories_);
    prefix_.resize(n);
    lastF_ = numCategories_;  // forces the first glue call to rebuild the suffix
    combined_.resize(n * numCategories_);
    forestResponse_.resize(n);
    forestWeights_.resize(n);
  }

  std::size_t numReportedLocations() const override { return numCategories_; }
  /// Each category forest reports its own splits: the per-sample varcount
  /// channel is K forests wide, slot k addressing category k's forest.
  std::size_t numVariableCountForests() const override { return numCategories_; }
  std::size_t variableCountForest(std::size_t j) const override { return j; }
  /// The softmax test blend is well-defined (unlike BCF): the K forests each
  /// accumulate their own totalTestFits in the sweep, and combinedTestFits maps
  /// them through the same softmax combinedFits applies to totalFits. So the
  /// recorded test channel carries the identified K test probabilities.
  bool testFitsAreDefined() const override { return true; }
  bool logLikelihoodIsDefined() const override { return false; }

  /// Interleaved PG draw for category f: form the current margin C_if and draw
  /// omega_if ~ PG(n_i, f_if - C_if) against it, storing both so
  /// formForestResponse(f), called immediately after, reads the SAME margin and
  /// precision. margins_ and omega_'s f-th column are the f-specific handoff.
  /// PG(n_i, psi) is the sum of n_i iid PG(1, psi) draws (the only shipped
  /// sampler); at n_i = 1 the loop is empty, so exactly one PG draw with the
  /// identical psi - the byte-identical single-trial reduction of the label path.
  void drawForestGlue(std::size_t f, ext_rng* rng,
                      const std::vector<Forest<L, ResidT>>& forests) override {
    std::size_t n = data_.numObservations;
    std::size_t K = numCategories_;
    if (f == 0 || f != lastF_ + 1) {
      // Fresh sweep entry (the sampler always starts at f == 0; a direct
      // out-of-order call lands here too). Every j > f fit is OLD at this
      // point, so snapshot the whole suffix by one backward two-way merge from
      // the empty top row: suffix_[g] = LSE over j > g. Seed the prefix (LSE
      // over j < f) from the below-f fits - empty, hence -inf, when f == 0.
      double* top = suffix_.data() + (K - 1) * n;
      for (std::size_t i = 0; i < n; ++i) top[i] = -HUGE_VAL;
      for (std::size_t g = K - 1; g-- > 0; ) {
        const double* next = forests[g + 1].totalFits.data();
        const double* above = suffix_.data() + (g + 1) * n;
        double* here = suffix_.data() + g * n;
        for (std::size_t i = 0; i < n; ++i)
          here[i] = logSumExp2(next[i], above[i]);
      }
      for (std::size_t i = 0; i < n; ++i) prefix_[i] = -HUGE_VAL;
      for (std::size_t g = 0; g < f; ++g) {
        const double* below = forests[g].totalFits.data();
        for (std::size_t i = 0; i < n; ++i)
          prefix_[i] = logSumExp2(prefix_[i], below[i]);
      }
    } else {
      // In-order continuation: category f - 1's tree update landed since its
      // glue call, so its fit is now NEW; fold it into the running prefix (the
      // interleaving's below-f mix). suffix_[f] still holds the OLD above-f mix.
      const double* prev = forests[f - 1].totalFits.data();
      for (std::size_t i = 0; i < n; ++i)
        prefix_[i] = logSumExp2(prefix_[i], prev[i]);
    }
    lastF_ = f;

    const double* fFits = forests[f].totalFits.data();
    const double* suffix = suffix_.data() + f * n;
    double* omega = omega_.data() + f * n;
    for (std::size_t i = 0; i < n; ++i) {
      double margin = logSumExp2(prefix_[i], suffix[i]);
      margins_[i] = margin;
      double psi = fFits[i] - margin;
      double draw = ext_rng_simulatePolyaGamma(rng, psi);
      for (int c = 1; c < trials_[i]; ++c)
        draw += ext_rng_simulatePolyaGamma(rng, psi);
      omega[i] = draw;
    }
  }

  /// Category f's PG working response (y_if - n_i/2)/omega_if + C_if under weight
  /// omega_if, formed from the count matrix the combiner owns (the passed chain y
  /// is ignored). Reads the margin and omega drawForestGlue(f) just stored. At
  /// n_i = 1 the one-hot y_if is in {0, 1} and trials * 0.5 is exactly 0.5, so
  /// this is the byte-identical single-trial reduction.
  ForestResponse formForestResponse(std::size_t f,
      const std::vector<Forest<L, ResidT>>& forests, const double* /*y*/,
      const double* /*w*/) override {
    (void) forests;
    std::size_t n = data_.numObservations;
    const double* omega = omega_.data() + f * n;
    for (std::size_t i = 0; i < n; ++i) {
      double yif = static_cast<double>(counts_[f * n + i]);
      forestResponse_[i] = (yif - trials_[i] * 0.5) / omega[i] + margins_[i];
      forestWeights_[i] = omega[i];
    }
    return {forestResponse_.data(), forestWeights_.data()};
  }

  /// The K softmax probabilities per observation, location-major (channel k at
  /// combined_[k*n + i]); the reported training output. Log-sum-exp-safe.
  const double* combinedFits(const std::vector<Forest<L, ResidT>>& forests) override {
    return blendSoftmax(data_.numObservations,
        [&](std::size_t k) { return forests[k].totalFits.data(); },
        combined_.data());
  }

  /// The K softmax TEST probabilities per test observation, location-major
  /// (channel k at combinedTest_[k*nTest + i]); the reported test output, the
  /// totalTestFits analog of combinedFits. Log-sum-exp-safe. The level-centering
  /// grand shift c (afterCombine) is added uniformly to every totalFits but NOT
  /// to totalTestFits; softmax is invariant to a common shift, so this recovers
  /// the identified test probabilities without it, which is why afterCombine
  /// leaves totalTestFits alone (docs/design/multinomial.md). combinedTest_ is
  /// sized here, not at construction, because numTestObservations may be set
  /// after the combiner is built (setTestPredictors); off the sweep hot path.
  const double* combinedTestFits(
      const std::vector<Forest<L, ResidT>>& forests) override {
    std::size_t nTest = data_.numTestObservations;
    combinedTest_.resize(nTest * numCategories_);
    return blendSoftmax(nTest,
        [&](std::size_t k) { return forests[k].totalTestFits.data(); },
        combinedTest_.data());
  }

  /// The likelihood-invariant LEVEL-CENTERING move (docs/design/multinomial.md):
  /// the softmax is invariant to a common shift of all f_ik, a flat additive
  /// direction the prior pins only weakly, so it mixes as a slow random walk.
  /// This move is a single DATASET-WIDE shift c added to every f_ik at once,
  /// absorbed UNIFORMLY: c/m_k lands on every occupied leaf of every one of
  /// forest k's m_k trees. A dataset-wide shift is the one flat direction a
  /// piecewise-constant forest can represent EXACTLY (add c to every
  /// observation's fit uniformly): it moves only the grand level and leaves
  /// every identified log-odds f_ij - f_ik untouched. A per-observation shift,
  /// by contrast, is not representable by coarse (shared-leaf) trees - the next
  /// backfit projects it, leaking spurious variance into the identified log-odds
  /// and biasing the softmax probabilities (a Jensen bias, confirmed against the
  /// exact gate) - so it is deliberately NOT used.
  ///
  /// The prior lives on LEAF VALUES (per-leaf sd
  /// s_k = leaf.scale_k / (k_k sqrt(m_k))), not on the total fits, so with L_k
  /// and S_k the count and value sum of forest k's OCCUPIED leaves over ALL its
  /// trees, the exact conditional of the uniformly absorbed shift is
  ///   prec = sum_k L_k / (m_k^2 s_k^2),  num = sum_k S_k / (m_k s_k^2),
  ///   c = -num/prec + N(0, 1)/sqrt(prec).
  /// Empty leaves are skipped: they carry no fit (sampleNodeParametersFromPrior
  /// does draw values for them, inert here). Uniform absorption is also the
  /// better mixing device than loading the whole shift onto one tree, whose
  /// constant-leaf conditional then fights it; at the intercept-only
  /// configuration it reduces to an exact independence sampler drawing the level
  /// from its marginal N(0, tau^2/K).
  ///
  /// The shift enters totalFits (which the next sweep's margins read) and the
  /// leaf tables, keeping the residual roll's total = sum-of-tree-fits invariant
  /// TO ROUNDING (m_k * fl(c/m_k) is not exactly c). Shifting muByTree wholesale
  /// is well defined only because a multinomial chain is ConstantGaussianLeaf-
  /// only (a leaf value IS the fit). totalTestFits is deliberately untouched -
  /// the softmax blend is invariant to the common shift. keepTrees is out of
  /// scope this arc, so the saved (flattened) tree leaves are not touched.
  /// Returns 1.0 (no multiplicative scale; the return feeds only BCF's test).
  double afterCombine(std::vector<Forest<L, ResidT>>& forests, bool /*record*/,
                      std::size_t /*sampleNum*/, ext_rng* rng) override {
    std::size_t n = data_.numObservations;
    double prec = 0.0, num = 0.0;
    for (std::size_t k = 0; k < numCategories_; ++k) {
      Forest<L, ResidT>& forest = forests[k];
      double m = static_cast<double>(forest.numTrees);
      if (!(m > 0.0)) continue;
      double s = forest.leaf.scale / (forest.k * std::sqrt(m));
      double invV = s > 0.0 ? 1.0 / (s * s) : 0.0;
      double leafCount = 0.0, leafSum = 0.0;
      // the bottom lists are left in each tree's scratch for the apply pass
      for (std::size_t t = 0; t < numTreesOf(forest); ++t) {
        Tree& tree = forest.trees[t];
        const std::vector<double>& mu = forest.muByTree[t];
        tree.bottomScratch.clear();
        tree.fillBottom(0, tree.bottomScratch);
        for (int32_t nodeIndex : tree.bottomScratch) {
          if (tree.at(nodeIndex).numObservations() == 0) continue;
          leafCount += 1.0;
          leafSum += mu[static_cast<std::size_t>(nodeIndex)];
        }
      }
      prec += leafCount * invV / (m * m);
      num += leafSum * invV / m;
    }
    if (!(prec > 0.0)) return 1.0;
    double c = -num / prec +
               ext_rng_simulateStandardNormal(rng) / std::sqrt(prec);

    for (std::size_t k = 0; k < numCategories_; ++k) {
      Forest<L, ResidT>& forest = forests[k];
      for (std::size_t i = 0; i < n; ++i) forest.totalFits[i] += c;
      if (forest.numTrees == 0) continue;
      double perTree = c / static_cast<double>(forest.numTrees);
      for (std::size_t t = 0; t < numTreesOf(forest); ++t) {
        Tree& tree = forest.trees[t];
        std::vector<double>& mu = forest.muByTree[t];
        for (int32_t nodeIndex : tree.bottomScratch) {
          if (tree.at(nodeIndex).numObservations() == 0) continue;
          mu[static_cast<std::size_t>(nodeIndex)] += perTree;
        }
      }
    }
    return 1.0;
  }

private:
  /// The live tree count the level-centering passes may touch: numTrees clamped
  /// to the tree and leaf-table lengths, so a partially built forest (a
  /// component-test fixture, a forest mid-resize) walks only what exists.
  static std::size_t numTreesOf(const Forest<L, ResidT>& forest) {
    std::size_t count = forest.numTrees;
    if (count > forest.trees.size()) count = forest.trees.size();
    if (count > forest.muByTree.size()) count = forest.muByTree.size();
    return count;
  }

  /// The shared train/test blend: gather the K forests' per-observation fits
  /// (forest k's supplied by fitsOf(k)) into out in location-major order
  /// (channel k at out[k*count + i]), then softmax in place. combinedFits and
  /// combinedTestFits are this same map; only the source vector (totalFits vs
  /// totalTestFits) and length (n vs nTest) differ. The gather is a plain copy
  /// in the same k-outer, i-inner order the two direct loops used, and the
  /// in-place softmax is unchanged, so both callers stay byte-identical.
  template <typename FitsOf>
  const double* blendSoftmax(std::size_t count, FitsOf fitsOf, double* out) {
    for (std::size_t k = 0; k < numCategories_; ++k) {
      const double* fits = fitsOf(k);
      for (std::size_t i = 0; i < count; ++i)
        out[k * count + i] = fits[i];
    }
    softmaxLocationMajor(out, count, numCategories_, out);
    return out;
  }

  /// Stable two-term log-sum-exp. An empty LSE is passed as -HUGE_VAL and
  /// returned through exactly (no exp(-inf)), so the empty-set margins are
  /// bit-exact: C_if for K == 2 is the other category's raw fit.
  static double logSumExp2(double a, double b) {
    if (a == -HUGE_VAL) return b;
    if (b == -HUGE_VAL) return a;
    double hi = a > b ? a : b;
    double lo = a > b ? b : a;
    return hi + std::log1p(std::exp(lo - hi));
  }

  const ColumnStore& data_;
  std::size_t numCategories_;
  const int* counts_;    // borrowed n x K count matrix, category-major (k*n + i)
  const int* trials_;    // borrowed per-observation trial count n_i (>= 1)
  std::vector<double> omega_;          // n x K, category-major (column k at k*n)
  std::vector<double> margins_;        // n; the current forest's C_if handoff
  std::vector<double> suffix_;         // n x K; per-sweep LSE over j > f, old fits
  std::vector<double> prefix_;         // n; running LSE over j < f, new fits
  std::size_t lastF_;                  // prior glue category; detects a fresh sweep
  std::vector<double> combined_;       // n x K softmax probabilities
  std::vector<double> combinedTest_;   // nTest x K softmax test probabilities
  std::vector<double> forestResponse_, forestWeights_;  // n each
};

}  // namespace bartcore

#endif  // BARTCORE_COMBINER_HPP
