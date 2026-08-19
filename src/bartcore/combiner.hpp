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
  /// The leaf prior's other half - the constant leaf is mu ~ N(0, (scale / k)^2)
  /// (model.hpp) - so a continuation that restores k without it would pair a
  /// donor's trees with the destination's calibration. OPTIONAL and append-only:
  /// a live leaf scale is strictly positive, so 0.0 is an unambiguous ABSENT
  /// sentinel and a state written before this block existed restores exactly as
  /// it did then. Both restore paths guard with > 0.0, which also swallows a
  /// non-finite or non-positive value - deliberately the same permissive
  /// posture k has, rather than a new refusal.
  double leafScale = 0.0;
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
  // The amplitude glue a combining response carries (docs/design/bcf.md);
  // false off a coupling that carries none. RAGGED, because a total is not a
  // layout: q = (1, 3) and q = (2, 2) both carry four amplitudes, so the
  // per-forest widths travel with them or a restore silently permutes the
  // blocks. amplitudes is forest-major, block f at the widths' prefix sum;
  // amplitudeVariances is one prior variance per forest, live only where a
  // scale mixture refreshes it.
  bool hasBCF = false;
  std::vector<std::size_t> amplitudeWidths;
  std::vector<double> amplitudes;
  std::vector<double> amplitudeVariances;
  // bcf's (a, aVariance, b0, b1) reading of the same block - the K = 2,
  // q = (1, 2) instance, in the spelling the two-forest fixtures state their
  // glue in and serializeGlue fills alongside the ragged fields. The RAGGED
  // fields are authoritative: restoreGlue reads these four only when
  // amplitudeWidths is empty, which is exactly a state written by hand rather
  // than by a combiner.
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
/// (ForestSpec below) derives them from the family's own latent scale at
/// construction, and the adaptive magnitude lives in the glue (a for mu,
/// b0/b1 for tau).
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

/// One forest's entry in the K-length spec a general combining chain is built
/// from: its tree/structure configuration (BCFForestSpec, unchanged) plus the
/// amplitude channel it enters the combination through.
///
/// The CALIBRATION MAP is stated here rather than in the chain, and is general
/// in K: the node scale is nodeScaleFactor * s / (nodeScaleDivisor * c), with s
/// the FAMILY'S OWN LATENT SCALE (the sample sd of the range-scaled response
/// under gaussian, 1 under probit, pi/sqrt(3) under logistic) and c the median
/// nonzero row norm of this forest's basis. The two ways a magnitude can be
/// carried are the two branches the host picks between - a scale-mixture
/// amplitude carries it in amplitudePriorScale and leaves the node scale at s
/// (bcf's prognostic forest, factor and divisor both 1), a fixed-variance
/// amplitude carries it in the node scale, divided through the half-normal
/// median so the amplitude spread times that scale sits at nodeScaleFactor
/// units of s (bcf's treatment forest, divisor 0.674). The divisor is CARRIED
/// rather than derived from the prior, so a forest declaring a zero scale keeps
/// the node scale its host asked for; and both expressions are written exactly
/// as bcf's were, which is what keeps the two-forest instance bitwise.
///
/// Every scale here is stated PER UNIT OF BASIS ROW NORM, which is what
/// dividing by c means: the induced prior sd on the combined index at row i is
/// sqrt( sum_f nodeScale_f^2 amplitudePriorVariance_f ||B_f(i,.)||^2 ) over the
/// fixed-variance forests, so a basis in other units would otherwise multiply
/// the prior a host asked for. That index is in sd(y) units under gaussian,
/// where a drawn sigma partly absorbs a mis-scaled basis, and in latent sd
/// units under probit and logistic, where sigma is PINNED and nothing does.
///
/// basis is an optional n x numBasisColumns ROW-major amplitude basis, borrowed
/// and COPIED at construction; null synthesizes the dense all-ones column whose
/// contraction is the amplitude itself.
struct ForestSpec {
  BCFForestSpec forest;
  double nodeScaleFactor = 1.0, nodeScaleDivisor = 1.0;
  double amplitudePriorVariance = 1.0;
  double amplitudePriorScale = 0.0;  // half-Cauchy median; 0 = fixed variance
  bool updateAmplitude = true;
  bool ridge = true;
  const double* basis = nullptr;
  std::size_t numBasisColumns = 1;
};

/// The spec a BCF chain is built from. Two readings of one object: the K = 2
/// mu/tau shape below, which is what every fixture and the shipped internal
/// route construct, and the K-LENGTH forests vector, which supersedes it
/// whenever it is non-empty. expandForestSpecs is the thin adapter between
/// them, so there is exactly one shape the chain and the combiner read.
struct BCFSpec {
  BCFForestSpec mu, tau;
  // The response family the K forests are combined under. gaussian, probit and
  // logistic are built; the rest are refused by the factory, which is where the
  // reason for each door lives. It rides the spec rather than the constructor's
  // parameter list so a fixture that names none keeps bcf's own family.
  ResponseFamily family = ResponseFamily::gaussian;
  const double* z = nullptr;    // borrowed 0/1 treatment indicator per obs
  // Half-Cauchy median for the mu scalar a. Deliberately NOT family-aware,
  // unlike the R default it shadows (defaultAmplitudePriorScale, R/model.R,
  // which is 1 under probit and logistic): this initializer is a FIXTURE
  // default that no consumer reaches - createBCFSampler has no flat-C entry
  // point, and the bridge fills spec.forests unconditionally - and the
  // two-forest spelling it belongs to is gaussian bcf's, where 2 is correct.
  // Branching on the family here would put one inside expandForestSpecs, the
  // one adapter the "bcf as the K = 2 instance" contract requires to stay a
  // pure re-spelling.
  double aPriorScale = 2.0;
  double bPriorVariance = 0.5;  // N(0, .) prior variance for b0, b1
  double sdModerate = 1.0;      // treatment effect scale, in units of s above
  bool updateA = true, updateB = true;  // false fixes the matching glue block
  // Whether each forest's amplitude block travels its own likelihood-invariant
  // ASIS ridge after the combination (BCFForestCombiner::afterCombine). mu's is
  // bcf's shipped a-move; tau's is the b-move the general rescale implements.
  // It is OFF here because switching it on consumes a GIG draw per sweep,
  // which re-records bcf-equivalence.
  bool ridgeA = true, ridgeB = false;
  // Whether the shipped K = 2 shape draws its amplitudes through the general
  // q-variate conditional rather than through the two-scalar path it shipped
  // with. The two agree in exact arithmetic and differ only in where the
  // compiler forms fused multiply-adds (BCFForestCombiner::drawGlue), so this
  // is off by default to hold bcf-equivalence bitwise.
  bool generalAmplitudeDraw = false;
  // The K-length reading. EMPTY leaves the mu/tau pair above authoritative and
  // the treatment forest's basis synthesized from z; non-empty supersedes both,
  // and z is then read by nothing.
  std::vector<ForestSpec> forests;
};

/// The thin adapter: bcf's two-forest spelling read as the K-length vector
/// every other layer works in. mu takes the half-Cauchy amplitude (bcf's a,
/// median aPriorScale) over the implicit intercept and leaves its node scale at
/// s; tau takes the fixed-variance amplitude pair (b0, b1) over the treatment
/// indicator basis and carries sdModerate in its node scale. A spec that
/// already carries its own forests vector is returned as it stands.
inline std::vector<ForestSpec> expandForestSpecs(const BCFSpec& spec) {
  if (!spec.forests.empty()) return spec.forests;
  std::vector<ForestSpec> forests(2);
  forests[0].forest = spec.mu;
  forests[0].amplitudePriorScale = spec.aPriorScale;
  forests[0].updateAmplitude = spec.updateA;
  forests[0].ridge = spec.ridgeA;
  forests[1].forest = spec.tau;
  forests[1].nodeScaleFactor = spec.sdModerate;
  forests[1].nodeScaleDivisor = 0.674;  // the half-normal median
  forests[1].amplitudePriorVariance = spec.bPriorVariance;
  forests[1].updateAmplitude = spec.updateB;
  forests[1].ridge = spec.ridgeB;
  forests[1].numBasisColumns = 2;  // synthesized from z by the combiner
  return forests;
}

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
  // borrowed n x K category offset in the same category-major layout, null for
  // none: the latent is f_ik + o_ik, so the offset enters the margins, the
  // working response and the reported softmax, never the leaf values
  const double* offset = nullptr;
  MultinomialForestSpec forest;       // shared across the K category forests
  // per-forest leaf scale anchor and k: the pairwise log-odds f_ik - f_ij is a
  // difference of two forests, so the per-forest node scale is the logistic
  // pi*sqrt(3) anchor divided by sqrt(2); tau = nodeScale / k is the per-forest
  // total-fit prior sd the level-centering move reads.
  double nodeScale = 3.847649490485592;  // pi*sqrt(3)/sqrt(2)
  double k = 2.0;
};

/// One forest's amplitude basis: the n x q matrix B_f whose row contracts with
/// the forest's amplitude vector a_f into the per-observation scalar the forest
/// enters the combination through, m_f(i) = dot(a_f, B_f(i, .)). Stored
/// ROW-major (row i at i * numColumns) because that contraction is the only
/// read: a row is contiguous and the multiplier costs one stream per forest.
///
/// q >= 1 always, and the width is uniform across rows, so there is no implicit
/// all-ones column to special-case - a forest whose multiplier is a plain
/// amplitude carries the ones column densely and reaches it by the same
/// contraction every other forest does. That is what leaves exactly one
/// multiplier path.
struct ForestBasis {
  std::vector<double> values;
  std::size_t numColumns = 0;
};

/// One forest's amplitude prior and the two switches its block carries. The
/// prior is N(0, variance) on every coordinate of the block; a positive
/// halfCauchyScale makes variance a LIVE inverse-gamma auxiliary refreshed
/// after the block's draw, which is bcf's half-Cauchy a (the scale mixture),
/// while a zero leaves it the fixed variance bcf's (b0, b1) carry.
struct ForestAmplitudePrior {
  double variance = 1.0;
  double halfCauchyScale = 0.0;
  bool update = true;  // false holds the whole block at its value
  bool ridge = true;   // whether afterCombine travels this block's ASIS orbit
};

/// The combining response's glue (docs/design/bcf.md): the per-forest amplitude
/// vectors and the bases they contract with, plus the sweep's per-forest
/// scratch. The combined location is a mu + b_z tau - the response mean under
/// gaussian (plus eps), the latent index under probit and logistic; a Chain
/// holds one only in BCF mode.
///
/// The amplitudes are ONE flat vector, forest f's block at amplitudeOffset[f]
/// and as wide as basis[f]: the general shape, of which bcf is the K = 2
/// instance (a) x (b0, b1). The named accessors below are that instance's
/// reading - the same three numbers, in their bcf spelling, for the glue
/// conditionals and the wire format that still speak it. They index THROUGH
/// amplitudeOffset rather than at 0/1/2: forest 1's block starts wherever
/// forest 0's ends, so a prognostic basis wider than one column moves it.
struct BCFState {
  std::vector<ForestBasis> basis;
  /// Per-forest IS-CANONICAL flag: whether basis[f] is still exactly one of the
  /// constructor's synthesized shapes - a dense all-ones column, or a
  /// complementary two-column 0/1 pair. A pure function of the basis VALUES,
  /// recomputed at install and at restore and never serialized, so nothing can
  /// carry a stale answer across a round trip. It selects the draw path:
  /// drawShippedGlue is not a general q = 2 conditional (it never reads
  /// basis[1] at all - it forms two disjoint group accumulators keyed on the
  /// indicator), so a legal continuous two-column basis has to route around it.
  std::vector<std::uint8_t> canonical;
  std::vector<double> amplitudes{1.0, 0.0, 1.0};
  // length K + 1, prefix sums of q_f; seeded with bcf's shipped layout so the
  // accessors below read the constructed amplitudes before any install
  std::vector<std::size_t> amplitudeOffset{0, 1, 3};
  std::vector<ForestAmplitudePrior> prior;
  std::vector<double> combined, forestResponse, forestWeights;
  // per-forest total-fit pointers, refilled at every combinedFits call: the K
  // generalization of the two pointers that blend hoisted out of its loop
  std::vector<const double*> fitsByForest;
  // the q-variate conditional's scratch: the q x q crossproduct (factorized in
  // place), its q-vector of moments (the draw in place), and the row of the
  // weighted design the two accumulate from
  std::vector<double> crossproduct, moments, designRow;

  double* amplitudesOf(std::size_t f) {
    return amplitudes.data() + amplitudeOffset[f];
  }
  const double* amplitudesOf(std::size_t f) const {
    return amplitudes.data() + amplitudeOffset[f];
  }
  std::size_t numAmplitudes(std::size_t f) const {
    return amplitudeOffset[f + 1] - amplitudeOffset[f];
  }

  double& a() { return amplitudes[amplitudeOffset[0]]; }
  double a() const { return amplitudes[amplitudeOffset[0]]; }
  double& b0() { return amplitudes[amplitudeOffset[1]]; }
  double b0() const { return amplitudes[amplitudeOffset[1]]; }
  double& b1() { return amplitudes[amplitudeOffset[1] + 1]; }
  double b1() const { return amplitudes[amplitudeOffset[1] + 1]; }
  double& aVariance() { return prior[0].variance; }
  double aVariance() const { return prior[0].variance; }
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

  /// Forest f's per-observation VETO precisions: the weights half of
  /// formForestResponse, formed alone so a caller holding no glue draw can
  /// still ask which rows forest f's likelihood reaches. The tree prior is
  /// conditioned on exactly that support (docs/design/empty-leaf-veto.md), so
  /// an override owes the vector its own formForestResponse would return, and
  /// may read nothing the response half draws. w is the chain's working
  /// precisions; the pointer aliases combiner scratch, valid until the next
  /// call. Pure rather than inert: a coupling defaulted to the global weights
  /// would draw initial forests from a law its own moves reject.
  virtual const double* formForestVetoWeights(std::size_t f,
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

  /// Installs forest f's own n x numColumns amplitude basis, ROW-major, COPIED;
  /// false, installing nothing, off a coupling that carries no amplitudes or on
  /// a basis this one refuses. It is the SOLE basis-mutation route: synthesis
  /// is construction-only, so there is no second operation whose ordering
  /// against this one would have to be specified. Inert in the base.
  virtual bool setForestBasis(std::size_t, const double*, std::size_t) {
    return false;
  }

  /// Whether this coupling OWNS its response as a count matrix that can be
  /// replaced on a live sampler at fixed n and K. True only for a combiner
  /// whose response is the borrowed counts it reads directly - the chain's y
  /// is not that response, so the response-side conduit
  /// (supportsResponseMutation) cannot express the swap and must not be
  /// widened to it. Defaults false so a future coupling stays refused at the
  /// bridge until it is audited.
  virtual bool supportsCountsMutation() const { return false; }

  /// Swaps the borrowed n x K count matrix and its per-observation trials;
  /// inert unless a subclass carries them. Both are borrowed for the sampler's
  /// lifetime, as the constructed pair is, and both are re-read from scratch on
  /// the next sweep - the coupling caches nothing derived from either.
  virtual void setCounts(const int*, const int*) {}

  /// Swaps the borrowed n x K offset on the coupling's own linear predictor,
  /// clearing it at a null pointer; inert unless a subclass carries one. It is
  /// NOT the response model's offset, which every reported channel adds after
  /// the forests are combined: this one enters BEFORE the combination, which is
  /// the only place a per-category shift means anything under a softmax.
  virtual void setCategoryOffset(const double*) {}

  /// The test-side twin: swaps the borrowed nTest x K offset the combined TEST
  /// location reads, clearing it at a null pointer; inert unless a subclass
  /// carries one. It is a separate object from the train offset - the test rows
  /// are other rows, so neither can be derived from the other - and it enters
  /// the test blend at the same place the train offset enters combinedFits.
  virtual void setCategoryTestOffset(const double*) {}

  /// The amplitude channel on the combining response, for the per-forest
  /// reporting path. totalAmplitudes is the CAPABILITY probe as well as the
  /// length: 0 says this coupling carries no amplitudes at all, which is how
  /// the "no glue" answer reaches the caller, and a coupling that carries any
  /// carries at least one per forest. numForestAmplitudes(f) is forest f's own
  /// width q_f (the vector is ragged), and amplitudes writes all
  /// totalAmplitudes() coordinates forest-major, block f at the widths' prefix
  /// sum.
  virtual std::size_t totalAmplitudes() const { return 0; }
  virtual std::size_t numForestAmplitudes(std::size_t) const { return 0; }
  virtual void amplitudes(double*) const {}

  /// The coupling draw and its likelihood-invariant post-combine move, fired at
  /// the fixed sweep points; inert unless a subclass couples the forests.
  ///
  /// afterCombine's return is a REPORTING channel, not a record of whether it
  /// moved: each override states its own convention, the sweep discards the
  /// value, and only the component tests read it (through
  /// Chain::interweaveGlueRidge). BCF's per-forest rescale returns the scale it
  /// applied to the forest it reports, 1.0 if that one held while another
  /// travelled; the multinomial shift returns 1.0 unconditionally, HAVING
  /// moved, because an additive move has no scale to report. So 1.0 does not
  /// mean the state is unchanged, and no caller may read it that way.
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

  /// Whether the per-draw per-forest reporting channels are defined: every
  /// forest's own function values and the scalars that recombine them into the
  /// per-observation location. True only for a coupling that composes its
  /// forests through such scalars (BCF: a mu + b_z tau, reported as the fits of
  /// both forests plus (a, b0, b1)); false otherwise, so a run over any other
  /// model allocates and computes nothing for either channel. storeSample fills
  /// them from forestTotalFits and bcfGlue, and the run bridge reads this
  /// predicate to decide whether the channels exist at all.
  virtual bool forestReportingIsDefined() const { return false; }

  /// How many per-observation channels combinedFits packs, one per reported
  /// location: 1 for every additive combiner (BCF's a mu + b_z tau is a single
  /// location), K for a model whose combined output is K per-observation values
  /// (softmax probabilities). storeSample loops its training/test writes over
  /// this count, and the run bridge sizes the fits arrays by it; L = 1 leaves
  /// the single-location byte layout every existing path relies on.
  virtual std::size_t numReportedLocations() const { return 1; }

  /// The variable-count reporting set: how many forests the per-sample split-
  /// usage channel CAN record, and which forest slot j addresses. A
  /// single-forest model reports the reported forest alone (count 1, slot 0 =
  /// reportedForest, this base); a model whose forests each carry their own
  /// splits reports all of them, slot j = forest j - multinomial's category
  /// forests and the multi-forest amplitude coupling's alike, so "additive"
  /// does not decide it. This is a CEILING, not the width a run writes:
  /// storeSample loops its varcount writes over the CALLER's declared
  /// Results::numVariableCountForests, which Sampler::run clamps to this count
  /// once and the run bridge sizes the varcount array by, so a caller
  /// declaring 1 (the flat C API, rbart_vi's callback loop) still gets slot 0
  /// and the single-forest byte layout its buffer is sized for. Distinct from
  /// numReportedLocations (softmax output channels): a model's fits-location
  /// count may diverge from its forest count - the amplitude coupling reports
  /// one location and every forest - so the varcount axis is keyed on its own
  /// forest set.
  virtual std::size_t numVariableCountForests() const { return 1; }
  virtual std::size_t variableCountForest(std::size_t) const {
    return reportedForest();
  }

  /// Whether the response-side conduit is safe on this coupling: the bridge
  /// gates the whole response swap, the offset swap (both at
  /// updateScale = false) and the case-weight swap on this one predicate, not
  /// setResponse alone. True only for a combiner that caches NOTHING per-forest
  /// across sweeps and re-derives every per-forest residual and precision from
  /// the chain's working response and weights on every sweep. Defaults false so
  /// a future multi-forest model stays refused at the bridge until it is
  /// audited. Named for the swap it was introduced for; the name does not
  /// narrow it.
  virtual bool supportsResponseMutation() const { return false; }

  /// Whether this coupling admits a caller-supplied per-forest, per-observation
  /// weight; Chain::setForestWeights states the semantics and is the only
  /// consumer. True only for a combiner whose per-forest precisions are a plain
  /// multiplicative factor on forest f's own leaf conditionals, so a row factor
  /// composed into them cannot invalidate anything the coupling carries.
  /// Defaults false so a future multi-forest model stays refused until it is
  /// audited.
  virtual bool supportsForestWeights() const { return false; }

  /// Installs (or clears, at a null pointer) the validated 0/1 active-row mask
  /// on a coupling that owns its OWN per-observation precisions; inert unless a
  /// subclass carries them. An additive coupling needs no override: its
  /// per-forest precisions are formed from the chain's working weights, which
  /// the response model has already composed the mask into.
  ///
  /// GLOBAL by construction - there is no forest index. A row is in the data
  /// set for the sampler or it is not; a coupling whose forests share one
  /// likelihood cannot restrict it forest by forest (MultinomialForestCombiner
  /// states the model reason).
  ///
  /// Chain::setActiveRows owns the single validating and normalizing scan and
  /// calls this after the response model's own install, so the values here are
  /// exactly 0 or 1, the length is the constructed n, and a null pointer means
  /// every row is active. The values are COPIED.
  virtual void setActiveRows(const double*) {}

  /// Glue (de)serialization into the ragged amplitude state fields; inert
  /// unless the combiner carries glue. glueIsValid is the LAYOUT check
  /// stateIsValid routes through: a state carrying the same total over
  /// different per-forest widths - q = (1, 3) against a live q = (2, 2) -
  /// would otherwise be written straight through the live offsets and silently
  /// permute the blocks. Accepting on a state that carries no glue keeps a
  /// mismatched restore a no-op rather than a refusal, as restoreGlue is.
  virtual void serializeGlue(ChainStateData&) const {}
  virtual void restoreGlue(const ChainStateData&) {}
  virtual bool glueIsValid(const ChainStateData&) const { return true; }
};

/// BCF's combiner (docs/design/bcf.md): a prognostic forest mu (forest 0) and a
/// treatment forest tau (forest 1) combined into the location a mu + b_z tau,
/// which is the response mean under gaussian (plus eps) and the latent index
/// under probit and logistic. Holds the glue (the scalar a via its half-Cauchy
/// scale mixture aVariance, and the treatment scales b0/b1 over control/treated)
/// and the sweep's per-forest scratch. Constant leaf only, as the whole BCF
/// chain is.
template <IntegrableLeafModel L, typename ResidT = double>
struct BCFForestCombiner : ForestCombiner<L, ResidT> {
  static_assert(!L::hasVectorParams && !L::hasFunctionParams,
                "BCF is a constant-leaf model");

  /// numForests is the chain's own forest count, which sizes every per-forest
  /// array here; a forest the spec does not reach enters with a plain amplitude
  /// and the default prior. It defaults to bcf's two so a fixture need not
  /// state it.
  ///
  /// Basis SYNTHESIS happens here and nowhere else. Under the K-length spec
  /// each forest takes the basis its entry carries, or the dense all-ones
  /// column when it carries none; under bcf's two-forest spelling forest 1
  /// takes the (1 - z, z) indicator pair z implies. There is deliberately no
  /// mid-life synthesis route: setForestBasis is the only mutator, so the
  /// question of which of two operations wins does not arise.
  BCFForestCombiner(const ColumnStore& data, const BCFSpec& spec,
                    std::size_t numForests = 2)
      : data_(data), numForests_(numForests < 2 ? 2 : numForests),
        generalAmplitudeDraw_(spec.generalAmplitudeDraw) {
    std::size_t n = data_.numObservations;
    std::vector<ForestSpec> forests = expandForestSpecs(spec);
    glue_.basis.resize(numForests_);
    glue_.prior.resize(numForests_);
    for (std::size_t f = 0; f < numForests_; ++f) {
      if (f < forests.size() && forests[f].basis != nullptr) {
        glue_.basis[f].numColumns = forests[f].numBasisColumns;
        glue_.basis[f].values.assign(
          forests[f].basis, forests[f].basis + n * forests[f].numBasisColumns);
      } else {
        glue_.basis[f].numColumns = 1;
        glue_.basis[f].values.assign(n, 1.0);
      }
      if (f < forests.size())
        glue_.prior[f] = {forests[f].amplitudePriorVariance,
                          forests[f].amplitudePriorScale,
                          forests[f].updateAmplitude, forests[f].ridge};
    }
    // bcf's two-forest spelling names its treatment basis by the indicator it
    // contrasts on rather than by the pair; the pair is that indicator's
    // two-level factor basis, whose amplitudes are exactly (b0, b1)
    if (spec.forests.empty() && numForests_ > 1) synthesizeIndicatorBasis(1, spec.z);
    rebuildAmplitudeLayout();
    refreshCanonical();
  }

  /// The SOLE basis-mutation route, at any forest and any width, and therefore
  /// the owner of the guards nothing else is left to apply: the index, a
  /// positive width, and finiteness. A refused install writes nothing.
  ///
  /// Ordering is LAST INSTALL WINS, per forest, and both orderings of a widen
  /// and a swap collapse to that one rule because there is only one operation:
  /// rebuildAmplitudeLayout derives the offsets as a pure prefix sum of the
  /// width vector and carries every block by position, so a widen-then-swap
  /// leaves the widening standing and a swap-then-widen carries the swapped
  /// values to their new offsets. Amplitudes PRESERVE and remap; a
  /// width-preserving install is the bitwise identity on every one of them.
  bool setForestBasis(std::size_t f, const double* values,
                      std::size_t numColumns) override {
    return installForestBasis(f, values, numColumns);
  }

  bool installForestBasis(std::size_t f, const double* values,
                          std::size_t numColumns) {
    std::size_t n = data_.numObservations;
    if (f >= glue_.basis.size() || numColumns == 0 || values == nullptr)
      return false;
    for (std::size_t i = 0; i < n * numColumns; ++i)
      if (!std::isfinite(values[i])) return false;
    glue_.basis[f].numColumns = numColumns;
    glue_.basis[f].values.assign(values, values + n * numColumns);
    rebuildAmplitudeLayout();
    glue_.canonical[f] = basisIsCanonical(f) ? 1 : 0;
    return true;
  }

  /// bcf's three amplitudes (a, b0, b1) off the general ragged channel: the
  /// K = 2, q = (1, 2) instance, kept as a named reading beside the general
  /// one for the conditionals and fixtures that speak it. False on any other
  /// layout, which is how a caller learns it is not looking at bcf.
  bool bcfGlue(double& a, double& b0, double& b1) const {
    if (glue_.amplitudes.size() != 3 || glue_.numAmplitudes(0) != 1)
      return false;
    a = glue_.a(); b0 = glue_.b0(); b1 = glue_.b1();
    return true;
  }

  std::size_t totalAmplitudes() const override {
    return glue_.amplitudes.size();
  }
  std::size_t numForestAmplitudes(std::size_t f) const override {
    return f < glue_.basis.size() ? glue_.numAmplitudes(f) : 0;
  }
  void amplitudes(double* out) const override {
    for (std::size_t j = 0; j < glue_.amplitudes.size(); ++j)
      out[j] = glue_.amplitudes[j];
  }

  /// Multipliers with |m_f| below this are taken as exactly zero:
  /// sqrt(DBL_EPSILON) = 2^-26 = 1.4901161193847656e-08, which R's all.equal,
  /// tinytest::expect_equal and testthat all default to when asked whether two
  /// numbers are the same. Written as a hex literal because 2^-26 is exact.
  static constexpr double zeroMultiplierTolerance = 0x1p-26;

  /// Effective response and precision forest f's constant-leaf draws see so that
  /// m_f * f_f explains the residual (y minus the other forests' scaled
  /// contributions): response r_i / m_f, weight w_i m_f^2, which reproduce the
  /// leaf's node sums without touching the leaf math.
  ///
  /// A multiplier indistinguishable from zero at that tolerance says row i
  /// carries no information about forest f, so the row leaves with exactly zero
  /// weight and exactly zero response. The zero RESPONSE is required rather than
  /// cosmetic: the chain reads this buffer arithmetically when it rolls the
  /// running residual and finalizes total fits, and the node sufficient-statistic
  /// kernels accumulate w * y, where a zero weight against an amplified response
  /// would be 0 * inf. The tolerance doubles as a condition-number cap: the
  /// division amplifies by at most 2^26, whatever the family, the weights and
  /// K. What the amplification does NOT reach is the arithmetic downstream,
  /// which cancels it exactly - the node kernels accumulate sumWeights =
  /// sum_i w_i m_i^2 and sumWeightedResponse = sum_i (w_i m_i^2)(r_i / m_i),
  /// whose exact value is sum_i w_i m_i r_i, and combinedFits and
  /// drawForestAmplitude keep the exact multiplier throughout. (An absolute
  /// precision claim would need a bounded numerator as well, which only
  /// gaussian supplies, by range-anchoring its working response to
  /// [-0.5, 0.5].)
  ///
  /// The snap belongs to this REPARAMETERIZATION, not to the model: combinedFits
  /// and drawGlue keep the exact multiplier, so a snapped row still receives
  /// m_f f_f(x_i) in the combination and still informs the glue draw.
  ///
  /// The residual accumulates FORWARD, subtracting the other forests in
  /// increasing index order from y. That is a CONTRACT once K > 2 - the
  /// direction is bitwise-visible there, as combinedFits' is at K = 2 - and it
  /// is deliberately NOT combinedFits' reverse: this sum has no two-term fused
  /// expression to reproduce, and the amplitude conditional forms the same
  /// residual the same way, so the two per-forest residuals agree by
  /// construction rather than by coincidence.
  ForestResponse formForestResponse(std::size_t f,
      const std::vector<Forest<L, ResidT>>& forests, const double* y,
      const double* w) override {
    std::size_t n = data_.numObservations;
    glue_.forestResponse.resize(n);
    glue_.forestWeights.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      double m = forestMultiplier(f, i);
      if (std::fabs(m) < zeroMultiplierTolerance) {
        glue_.forestResponse[i] = 0.0;
        glue_.forestWeights[i] = 0.0;
        continue;
      }
      double resid = y[i];
      for (std::size_t g = 0; g < forests.size(); ++g)
        if (g != f) resid -= forestMultiplier(g, i) * forests[g].totalFits[i];
      glue_.forestResponse[i] = resid / m;
      glue_.forestWeights[i] = (w == nullptr ? 1.0 : w[i]) * m * m;
    }
    return {glue_.forestResponse.data(), glue_.forestWeights.data()};
  }

  /// The precisions above without the response: w_i m_f(i)^2 under the same
  /// near-zero snap, so a row the reparameterization drops is a row no leaf of
  /// forest f may hold alone. The snap is what makes bcf's creation glue visible
  /// here - b0 = 0 leaves every control row weightless in the treatment forest.
  const double* formForestVetoWeights(std::size_t f, const double* w) override {
    std::size_t n = data_.numObservations;
    glue_.forestWeights.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      double m = forestMultiplier(f, i);
      glue_.forestWeights[i] = std::fabs(m) < zeroMultiplierTolerance
        ? 0.0
        : (w == nullptr ? 1.0 : w[i]) * m * m;
    }
    return glue_.forestWeights.data();
  }

  /// The combined location, sum_f m_f(i) f_f(x_i), accumulated WITHIN the row
  /// and from the LAST forest DOWN.
  ///
  /// The direction is load-bearing and measured, not stylistic. Under fused
  /// multiply-add contraction (clang's default for C++) exactly one product in
  /// a sum escapes its own rounding - the one the closing add absorbs - and the
  /// two-forest blend this replaces, `a mu + b_z tau` as a single expression,
  /// absorbed forest 0's. Seeding the accumulator with the LAST forest leaves
  /// forest 0's product as the only bare multiply in the closing add, so it is
  /// the one contracted and bcf's K = 2 instance stays bitwise. Accumulating
  /// FORWARD absorbs the last forest's product instead and moves ~30% of rows
  /// by one ulp, which walks into the whole trajectory - measured, all 12
  /// bcf-equivalence scenarios red on mu, tau, glue, sigma and train.
  ///
  /// The multiplier here is the EXACT contraction: the reparameterization's
  /// near-zero snap belongs to formForestResponse and must not be shared with
  /// this reader, which still owes a snapped row its m_f f_f.
  const double* combinedFits(const std::vector<Forest<L, ResidT>>& forests) override {
    std::size_t n = data_.numObservations;
    std::size_t numForests = forests.size();
    glue_.combined.resize(n);
    glue_.fitsByForest.resize(numForests);
    for (std::size_t f = 0; f < numForests; ++f)
      glue_.fitsByForest[f] = forests[f].totalFits.data();
    std::size_t last = numForests - 1;
    for (std::size_t i = 0; i < n; ++i) {
      double location = forestMultiplier(last, i) * glue_.fitsByForest[last][i];
      for (std::size_t f = last; f-- > 0;)
        location += forestMultiplier(f, i) * glue_.fitsByForest[f][i];
      glue_.combined[i] = location;
    }
    return glue_.combined.data();
  }

  /// The amplitude conditionals, forest by forest in index order: each block is
  /// drawn jointly from its Gaussian full conditional given the other forests,
  /// and a scale-mixture prior's variance is refreshed immediately after its own
  /// block (docs/design/bcf.md - a as the mu coefficient under the half-Cauchy
  /// mixture, b0/b1 as the tau coefficients over control/treated).
  ///
  /// The SHIPPED K = 2 shape keeps the two-scalar path it landed with, and the
  /// general q-variate conditional draws every other shape. The two are the
  /// SAME conditional in exact arithmetic - measured, all four accumulators
  /// bitwise equal under -ffp-contract=off - and differ only in where the
  /// compiler forms fused multiply-adds: the a block accumulates in one
  /// statement and fuses, while the b block's per-row products are formed
  /// before a branch and accumulated inside it, which fuses unevenly.
  ///
  /// The measured split, against the general loop on the equivalence fixtures,
  /// because this comment is the trigger for deleting the branch below and a
  /// re-measurement has to be able to check it: ALL FOUR PRECISIONS reproduce
  /// bitwise, weighted and unweighted; the divergence is in the two MOMENTS -
  /// unweighted, n1 reproduces and n0 differs; weighted, both differ. No single
  /// accumulation shape reproduces both blocks (21 variants tried), so the
  /// general path CANNOT be bitwise on bcf and the specialized one is kept
  /// until a bcf-equivalence re-record is authorized, at which point
  /// BCFSpec::generalAmplitudeDraw becomes the default and this branch is
  /// deleted.
  ///
  /// Path selection is a per-forest IS-CANONICAL VALUE predicate, not a width
  /// test. The widths alone would admit a continuous two-column basis into the
  /// shipped path, which never reads basis[1] as a DESIGN MATRIX: it tests
  /// column 1 for nonzero as a group key and never multiplies by the stored
  /// values, forming two disjoint group-precision accumulators instead. So on a
  /// 0.25/0.75 pair it would silently draw a different model. A non-canonical
  /// basis at ANY forest therefore forces the general path for the whole draw.
  void drawGlue(ext_rng* rng, double sigma, const double* y, const double* w,
                const std::vector<Forest<L, ResidT>>& forests) override {
    if (forests.size() == 2 && !generalAmplitudeDraw_ && shippedShape())
      drawShippedGlue(rng, sigma, y, w, forests);
    else
      drawAmplitudes(rng, sigma, y, w, forests);
  }

  /// Interweaving (ASIS, Yu & Meng 2011) rescale, per forest, of the amplitude
  /// ridge each forest's block carries. After the amplitude draws and the leaf
  /// draws, forest f's L + q scale coordinates (a_f, mu_1..mu_L) travel
  /// (a_f/c, c mu_l) along the likelihood-invariant orbit
  /// dot(a_f, B_f(i,.)) f_f(x_i) = dot(a_f/c, B_f(i,.)) (c f_f(x_i)), so the
  /// move updates only the amplitude coordinates and preserves the posterior.
  ///
  /// The forests are travelled in index order, one GIG draw each; the blocks
  /// are DISJOINT, so the moves commute and each is an exact Gibbs update given
  /// the rest, and the order is a stream convention rather than a modelling
  /// choice. Returns the scale applied to
  /// the forest this combiner reports - 1.0 if that forest held, which does not
  /// say that no forest moved.
  ///
  /// Instantiated at bcf's prognostic forest (q = 1, the scale mixture) this is
  /// the shipped a-move exactly, and at its treatment forest (q = 2, a fixed
  /// prior variance) it is the b-move - one mechanism, not two, by the same
  /// general exponent rule.
  double afterCombine(std::vector<Forest<L, ResidT>>& forests, bool record,
                      std::size_t sampleNum, ext_rng* rng) override {
    double reported = 1.0;
    for (std::size_t f = 0; f < forests.size(); ++f) {
      const ForestAmplitudePrior& prior = glue_.prior[f];
      if (!prior.update || !prior.ridge) continue;
      double c = rescaleAmplitudeRidge(f, forests[f], record, sampleNum, rng);
      if (f == this->reportedForest()) reported = c;
    }
    return reported;
  }

  /// BCF reports the prognostic forest (forest 0, the base default) but leaves
  /// the test-fit and log-likelihood channels undefined: there is no test
  /// treatment vector to blend a mu + b_z tau off-sample, and the blended
  /// per-observation location is not visible to the response model to score.
  bool testFitsAreDefined() const override { return false; }
  bool logLikelihoodIsDefined() const override { return false; }

  /// Both per-forest channels are defined here: mu and tau are the forests'
  /// own function values, and (a, b0, b1) is exactly the scalar set that
  /// recombines them, so a consumer reads both surfaces off one run.
  bool forestReportingIsDefined() const override { return true; }

  /// Each forest records its own splits, as the softmax coupling's do: the
  /// per-sample varcount channel is numForests_ wide, slot j addressing forest
  /// j, so the prognostic forest keeps slot 0 (the exact bytes a caller
  /// declaring one forest still gets) and the treatment forest follows it.
  /// Chain::numVariableCountForests clamps this to the chain's own forest
  /// count: the constructor rounds numForests_ UP to two, so a combiner built
  /// standalone below that would otherwise report a forest the chain has not
  /// got and storeSample would index past forests_.
  std::size_t numVariableCountForests() const override { return numForests_; }
  std::size_t variableCountForest(std::size_t j) const override { return j; }

  /// BCF admits the scale-pinned response and offset swaps on two conditions,
  /// both of which hold for every family the coupling is built with. FIRST,
  /// the response re-maps y and touches no forest, so every leaf calibration
  /// and the sigma prior stay put: under gaussian that is the re-map through
  /// the pinned (min_, range_), and under a latent family the transform is the
  /// identity and sigma pinned outright, the swap redrawing the latents
  /// against the COMBINED location the chain hands it. SECOND, this combiner
  /// caches nothing per-forest across sweeps - formForestResponse re-derives
  /// every per-forest residual and precision from the chain's y and w each
  /// sweep, and the glue lives here rather than in the response, so it carries
  /// over; read as the WORKING y and w, which the response itself refreshes,
  /// that holds under a latent family unchanged. The same two conditions open
  /// the case-weight swap, which needs no pinned scale of its own: setWeights
  /// does not move the transform, and the map's anchor - the response sd under
  /// gaussian, the link's own error sd otherwise - is unweighted either way.
  bool supportsResponseMutation() const override { return true; }

  /// BCF's per-forest precision is w_i m_f^2, re-derived from y and w every
  /// sweep and read by nothing else, so a caller-supplied row factor composes
  /// into it exactly and carries nothing across sweeps.
  bool supportsForestWeights() const override { return true; }

  /// The amplitude block into and out of the ragged wire format: the per-forest
  /// widths, the flat amplitude vector, and each forest's prior variance.
  /// serializeGlue owns the hasBCF flag (the "carries glue" marker);
  /// restoreGlue is a no-op on a state that carries none, so a mismatched
  /// restore leaves the glue at its constructed values, and glueIsValid has
  /// already refused a state whose LAYOUT differs - a total is not a layout,
  /// and restoreGlue writes through the live offsets.
  ///
  /// The bases themselves are NOT here. They are model configuration, carried
  /// by the host across a re-creation the way the design matrix is, and they
  /// arrive at CONSTRUCTION - so a widening applied after a restore preserves
  /// and remaps the RESTORED amplitudes rather than the constructed ones.
  void serializeGlue(ChainStateData& state) const override {
    state.hasBCF = true;
    std::size_t numForests = glue_.basis.size();
    state.amplitudeWidths.resize(numForests);
    state.amplitudeVariances.resize(numForests);
    for (std::size_t f = 0; f < numForests; ++f) {
      state.amplitudeWidths[f] = glue_.numAmplitudes(f);
      state.amplitudeVariances[f] = glue_.prior[f].variance;
    }
    state.amplitudes = glue_.amplitudes;
    if (numForests == 2 && state.amplitudeWidths[0] == 1 &&
        state.amplitudeWidths[1] == 2) {
      state.a = glue_.a();
      state.aVariance = glue_.aVariance();
      state.b0 = glue_.b0();
      state.b1 = glue_.b1();
    }
  }
  void restoreGlue(const ChainStateData& state) override {
    if (!state.hasBCF) return;
    if (state.amplitudeWidths.empty()) {
      // a hand-written bcf-shaped state: the four named scalars are the whole
      // block, and only on the layout they name
      if (glue_.amplitudes.size() != 3 || glue_.numAmplitudes(0) != 1) return;
      glue_.a() = state.a;
      glue_.aVariance() = state.aVariance;
      glue_.b0() = state.b0;
      glue_.b1() = state.b1;
      return;
    }
    if (!glueIsValid(state)) return;
    glue_.amplitudes = state.amplitudes;
    for (std::size_t f = 0; f < glue_.prior.size(); ++f)
      glue_.prior[f].variance = state.amplitudeVariances[f];
    refreshCanonical();
  }
  bool glueIsValid(const ChainStateData& state) const override {
    if (!state.hasBCF || state.amplitudeWidths.empty()) return true;
    std::size_t numForests = glue_.basis.size();
    if (state.amplitudeWidths.size() != numForests ||
        state.amplitudeVariances.size() != numForests ||
        state.amplitudes.size() != glue_.amplitudes.size())
      return false;
    for (std::size_t f = 0; f < numForests; ++f)
      if (state.amplitudeWidths[f] != glue_.numAmplitudes(f)) return false;
    return true;
  }

private:
  /// The two-scalar conditional the shipped K = 2 shape keeps: a against the
  /// residual net of b_z tau, its scale-mixture variance, then b0 and b1
  /// against the residual net of the NEW a. drawGlue states why this path
  /// survives beside the general one; it is otherwise the general one at
  /// q = 1 and q = 2 over an orthogonal basis, written out.
  ///
  /// The grouping is read from basis[1]'s TREATED column, not from a borrowed
  /// z. There is no borrowed z any more: it had one writer (the retired
  /// synthesis route) and only these two readers, so keeping it while
  /// setForestBasis became the sole mutator would freeze it at its
  /// construction value - a width-preserving swap would install the new pair,
  /// leave the layout unmoved, and then partition by the OLD indicator while
  /// forestMultiplier contracted the NEW basis. The stored column holds
  /// exactly the values z did, so this is bitwise on every shipped route.
  void drawShippedGlue(ext_rng* rng, double sigma, const double* y,
                       const double* w,
                       const std::vector<Forest<L, ResidT>>& forests) {
    std::size_t n = data_.numObservations;
    const double* mu = forests[0].totalFits.data();
    const double* tau = forests[1].totalFits.data();
    const double* treated = glue_.basis[1].values.data();  // row-major, col 1
    double invSigmaSq = 1.0 / (sigma * sigma);

    if (glue_.prior[0].update) {
      double aPrec = 1.0 / glue_.prior[0].variance, aNum = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        double wi = w == nullptr ? 1.0 : w[i];
        double bz = treated[2 * i + 1] != 0.0 ? glue_.b1() : glue_.b0();
        double r = y[i] - bz * tau[i];
        aPrec += wi * mu[i] * mu[i] * invSigmaSq;
        aNum += wi * mu[i] * r * invSigmaSq;
      }
      glue_.a() =
        aNum / aPrec + ext_rng_simulateStandardNormal(rng) / std::sqrt(aPrec);

      // t_1 scale mixture: aVariance ~ IG(1/2, scale^2/2) mixes N(0, aVariance)
      // to Cauchy(0, scale), so the conditional's rate carries scale^2, not its
      // inverse. Gated on a POSITIVE scale, exactly as drawAmplitudes is: a
      // zero scale (reachable - the bridge passes aPriorScale through with no
      // positivity guard) is not a scale mixture at all, and refreshing there
      // would leave the two paths different MODELS rather than the same
      // conditional. Every shipped route sets a positive scale, so the gate is
      // rng-neutral on every baseline.
      if (glue_.prior[0].halfCauchyScale > 0.0) {
        double rate = 0.5 * glue_.a() * glue_.a() +
                      0.5 * glue_.prior[0].halfCauchyScale *
                        glue_.prior[0].halfCauchyScale;
        glue_.prior[0].variance =
          1.0 / ext_rng_simulateGamma(rng, 1.0, 1.0 / rate);
      }
    }

    if (glue_.prior[1].update) {
      double bPrec = 1.0 / glue_.prior[1].variance;
      double p0 = bPrec, n0 = 0.0, p1 = bPrec, n1 = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        double wi = w == nullptr ? 1.0 : w[i];
        double r = y[i] - glue_.a() * mu[i];
        double prec = wi * tau[i] * tau[i] * invSigmaSq;
        double num = wi * tau[i] * r * invSigmaSq;
        if (treated[2 * i + 1] != 0.0) { p1 += prec; n1 += num; }
        else { p0 += prec; n0 += num; }
      }
      glue_.b0() =
        n0 / p0 + ext_rng_simulateStandardNormal(rng) / std::sqrt(p0);
      glue_.b1() =
        n1 / p1 + ext_rng_simulateStandardNormal(rng) / std::sqrt(p1);
    }
  }

  /// The general sweep over the amplitude blocks: forest by forest in INDEX
  /// order, each block's q-variate conditional drawn given the current value of
  /// every other block, so the pass is a Gibbs scan and a block sees the blocks
  /// before it already updated. A scale-mixture prior's variance is refreshed
  /// straight after its own block, from the q-variate inverse-gamma
  /// conditional IG((1 + q)/2, (scale^2 + ||a_f||^2)/2) - at q = 1 the shape is
  /// exactly the shipped path's 1.0.
  void drawAmplitudes(ext_rng* rng, double sigma, const double* y,
                      const double* w,
                      const std::vector<Forest<L, ResidT>>& forests) {
    double invSigmaSq = 1.0 / (sigma * sigma);
    for (std::size_t f = 0; f < forests.size(); ++f) {
      ForestAmplitudePrior& prior = glue_.prior[f];
      if (!prior.update) continue;
      drawForestAmplitude(f, rng, invSigmaSq, y, w, forests);
      if (!(prior.halfCauchyScale > 0.0)) continue;
      std::size_t q = glue_.numAmplitudes(f);
      const double* amplitude = glue_.amplitudesOf(f);
      double squaredNorm = amplitude[0] * amplitude[0];
      for (std::size_t j = 1; j < q; ++j)
        squaredNorm += amplitude[j] * amplitude[j];
      double rate = 0.5 * squaredNorm +
                    0.5 * prior.halfCauchyScale * prior.halfCauchyScale;
      prior.variance = 1.0 / ext_rng_simulateGamma(
        rng, 0.5 * (1.0 + static_cast<double>(q)), 1.0 / rate);
    }
  }

  /// Forest f's amplitude block, drawn jointly from its Gaussian full
  /// conditional. The block's design row is its basis row scaled by the
  /// forest's own fit, x_i = B_f(i,.) f_f(x_i), against the residual y net of
  /// every other forest's contribution, so the conditional's precision is
  /// P = I/priorVar + sum_i w_i x_i x_i' / sigma^2 and its first moment
  /// sum_i w_i x_i r_i / sigma^2. The prior term is what keeps P positive
  /// definite whatever the basis does, which is why the factorization below
  /// needs no failure path.
  ///
  /// ACCUMULATION CONVENTIONS, all three observable once q > 1 or K > 2, and
  /// all three FORWARD here. The residual subtracts the other forests in
  /// increasing index order (formForestResponse's convention, the one a
  /// per-forest residual has always used - and NOT combinedFits' reverse, which
  /// is load-bearing for a different reason: reproducing a two-term fused
  /// expression). The crossproduct and moment sums accumulate over rows in
  /// index order, and the design row contracts forward over columns, exactly as
  /// forestMultiplier does. Each is a bitwise-visible choice; changing any of
  /// them moves draws.
  ///
  /// The factorization is the square-root-free L D L' rather than the Cholesky
  /// the leaf models use, because only its solve reduces to one division per
  /// coordinate: over an ORTHOGONAL basis - bcf's indicator pair, any factor
  /// basis - the unit triangles are exactly identity, so the q-variate draw is
  /// q scalar draws bitwise, in coordinate order, one standard normal each.
  void drawForestAmplitude(std::size_t f, ext_rng* rng, double invSigmaSq,
                           const double* y, const double* w,
                           const std::vector<Forest<L, ResidT>>& forests) {
    std::size_t n = data_.numObservations;
    const ForestBasis& basis = glue_.basis[f];
    std::size_t q = basis.numColumns;
    const double* fits = forests[f].totalFits.data();

    std::vector<double>& crossproduct = glue_.crossproduct;
    std::vector<double>& moments = glue_.moments;
    std::vector<double>& row = glue_.designRow;
    crossproduct.assign(q * q, 0.0);
    moments.assign(q, 0.0);
    row.resize(q);
    double priorPrecision = 1.0 / glue_.prior[f].variance;
    for (std::size_t j = 0; j < q; ++j) crossproduct[j * q + j] = priorPrecision;

    for (std::size_t i = 0; i < n; ++i) {
      double wi = w == nullptr ? 1.0 : w[i];
      double resid = y[i];
      for (std::size_t g = 0; g < forests.size(); ++g)
        if (g != f) resid -= forestMultiplier(g, i) * forests[g].totalFits[i];
      const double* values = basis.values.data() + i * q;
      for (std::size_t j = 0; j < q; ++j) row[j] = values[j] * fits[i];
      for (std::size_t j = 0; j < q; ++j) {
        for (std::size_t k = 0; k <= j; ++k)
          crossproduct[j * q + k] += wi * row[j] * row[k] * invSigmaSq;
        moments[j] += wi * row[j] * resid * invSigmaSq;
      }
    }

    // the draw is L'^-1 (D^-1 L^-1 num + D^-1/2 e): the conditional mean and a
    // normal deviate with covariance P^-1, formed together so the pivot divides
    // the moment exactly once
    unitLowerDecompose(crossproduct.data(), q);
    solveUnitLowerTriangular(crossproduct.data(), q, moments.data());
    for (std::size_t j = 0; j < q; ++j) {
      double pivot = crossproduct[j * q + j];
      moments[j] = moments[j] / pivot +
                   ext_rng_simulateStandardNormal(rng) / std::sqrt(pivot);
    }
    solveUnitLowerTriangularTransposed(crossproduct.data(), q, moments.data());

    double* amplitude = glue_.amplitudesOf(f);
    for (std::size_t j = 0; j < q; ++j) amplitude[j] = moments[j];
  }

  /// Forest f's ASIS ridge: draw c and travel (a_f, leaves) -> (a_f/c, c
  /// leaves). c = sqrt(v), v ~ GIG((L - q)/2, M/leafVar, ||a_f||^2/priorVar)
  /// with L and M the count and squared sum of f's OCCUPIED leaves. The
  /// exponent follows a general rule: rescaling k leaf parameters against d
  /// glue scalars gives p = (k - d)/2, so
  /// q = 1 is the shipped (L - 1)/2 and q = 2 the b-move's (L - 2)/2; the naive
  /// move-map Jacobian's (L - q + 1)/2 is off by one and its prototype rejects
  /// it at KS 1.6e-21 (derivation and prototype at
  /// docs/design/multiplier-combiner.md, "The exponent rule"). B reads the
  /// LIVE prior variance, which for a scale mixture is the auxiliary this
  /// move conditions on (refreshing it here would re-randomize the coordinate
  /// just conditioned on and measurably throttle the mixing gain - IACT 69 ->
  /// 196 on |a|, recorded at the same doc's "The ASIS ridge"); the one-sweep
  /// lag is benign, the next drawGlue refreshing it | a_new.
  ///
  /// A no-op consuming no rng below two occupied leaves or at a zero leaf sum,
  /// returning 1.0. record/sampleNum locate the keepTrees saved slot whose
  /// leaves, flattened before this move, need the same c so a stored
  /// amplitude * leaf keeps the identified product; the recorded test surface
  /// travels for state self-consistency.
  double rescaleAmplitudeRidge(std::size_t f, Forest<L, ResidT>& forest,
                               bool record, std::size_t sampleNum,
                               ext_rng* rng) {
    std::size_t n = data_.numObservations;

    // L, M over the occupied leaves. Recomputed unconditionally: the
    // k-accumulator that would hold these is gated on updateK, which BCF
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

    double* amplitude = glue_.amplitudesOf(f);
    std::size_t q = glue_.numAmplitudes(f);
    double squaredNorm = amplitude[0] * amplitude[0];
    for (std::size_t j = 1; j < q; ++j)
      squaredNorm += amplitude[j] * amplitude[j];
    double leafPrecision = (forest.k / forest.leaf.scale) *
                           (forest.k / forest.leaf.scale);  // 1 / leafVar
    double gigP =
      0.5 * (static_cast<double>(numLeaves) - static_cast<double>(q));
    double gigA = M * leafPrecision;
    double gigB = squaredNorm / glue_.prior[f].variance;

    double v = ext_rng_simulateGeneralizedInverseGaussian(rng, gigP, gigA,
                                                          gigB);
    if (!std::isfinite(v) || v <= 0.0) return 1.0;
    double c = std::sqrt(v);
    if (!std::isfinite(c) || c <= 0.0) return 1.0;

    // travel the ridge: the amplitudes shrink, the forest's fits grow by c.
    // Scaling every leaf value scales every gathered fit, so the leaf tables
    // carry the rescale.
    for (std::size_t j = 0; j < q; ++j) amplitude[j] = amplitude[j] / c;
    for (std::size_t t = 0; t < forest.numTrees; ++t)
      misc_scalarMultiplyVectorInPlace(forest.muByTree[t].data(),
                                       forest.muByTree[t].size(), c);
    misc_scalarMultiplyVectorInPlace(forest.totalFits.data(), n, c);

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

  /// The scale forest f's constant leaf carries into the combination: the
  /// contraction of forest f's amplitude vector with its own basis row,
  /// m_f(i) = dot(a_f, B_f(i, .)). One path, whatever K and whatever q_f - bcf
  /// reaches its two values through it rather than around it: the prognostic
  /// forest's all-ones column contracts to a on every row, and the treatment
  /// forest's (control, treated) indicator pair contracts to b_{z_i}, both
  /// bitwise (a * 1 is a, and the zero-weighted term adds a zero).
  ///
  /// q >= 1 by construction, so the first product seeds the sum and no
  /// zero-initialized accumulator enters the arithmetic. The contraction runs
  /// FORWARD over the columns, which is a CONTRACT and not an implementation
  /// detail: at q > 2 a reassociation moves the multiplier by an ulp and every
  /// reader of it - the blend, the reparameterization, the amplitude
  /// conditional - with it.
  double forestMultiplier(std::size_t f, std::size_t i) const {
    const ForestBasis& basis = glue_.basis[f];
    const double* amplitude =
      glue_.amplitudes.data() + glue_.amplitudeOffset[f];
    const double* row = basis.values.data() + i * basis.numColumns;
    double m = amplitude[0] * row[0];
    for (std::size_t j = 1; j < basis.numColumns; ++j)
      m += amplitude[j] * row[j];
    return m;
  }

  /// SYNTHESIZES forest f's basis from a 0/1 indicator: the two-column
  /// (control, treated) pair, which is that indicator's two-level factor basis
  /// and whose amplitudes are exactly (b0, b1). Called from the CONSTRUCTOR
  /// only - basis synthesis is construction-only, so a widening
  /// cannot be silently reset by a later data swap. A null indicator installs
  /// the all-control basis, matching the constructed b0.
  ///
  /// EVERY per-forest array is sized by the chain's own forest count, not by
  /// bcf's two: combinedFits, the reparameterization and the amplitude
  /// conditional all index basis[f] over the forests they are handed, so a
  /// chain carrying more of them than the spec calibrates must still find a
  /// basis at each one. Those extra forests carry a plain amplitude - the dense
  /// all-ones column, whose contraction is the amplitude itself - and the
  /// default prior.
  void synthesizeIndicatorBasis(std::size_t f, const double* z) {
    std::size_t n = data_.numObservations;
    ForestBasis& basis = glue_.basis[f];
    basis.numColumns = 2;
    basis.values.resize(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
      bool treated = z != nullptr && z[i] != 0.0;
      basis.values[2 * i] = treated ? 0.0 : 1.0;
      basis.values[2 * i + 1] = treated ? 1.0 : 0.0;
    }
  }

  /// Whether basis[f] still holds one of the constructor's synthesized shapes:
  /// a dense all-ones column, or a complementary two-column 0/1 pair (each
  /// entry in {0, 1}, each row summing to 1). Any other width, and any values
  /// off those two shapes, is non-canonical - a MODEL fact, since it is what
  /// selects the amplitude conditional. Recomputed rather than tracked, so
  /// there is no flag to serialize and none to go stale.
  bool basisIsCanonical(std::size_t f) const {
    std::size_t n = data_.numObservations;
    const ForestBasis& basis = glue_.basis[f];
    const double* values = basis.values.data();
    if (basis.numColumns == 1) {
      for (std::size_t i = 0; i < n; ++i)
        if (values[i] != 1.0) return false;
      return true;
    }
    if (basis.numColumns != 2) return false;
    for (std::size_t i = 0; i < n; ++i) {
      double control = values[2 * i], treated = values[2 * i + 1];
      if ((control != 0.0 && control != 1.0) ||
          (treated != 0.0 && treated != 1.0) || control + treated != 1.0)
        return false;
    }
    return true;
  }

  void refreshCanonical() {
    glue_.canonical.resize(glue_.basis.size());
    for (std::size_t f = 0; f < glue_.basis.size(); ++f)
      glue_.canonical[f] = basisIsCanonical(f) ? 1 : 0;
  }

  /// Whether the two bases are still bcf's own - forest 0 the plain intercept,
  /// forest 1 the complementary indicator pair - which is what the two-scalar
  /// draw is written against. The widths are checked alongside the value
  /// predicate because canonical says "one of the two shapes", not which.
  bool shippedShape() const {
    return glue_.basis.size() == 2 && glue_.canonical[0] && glue_.canonical[1] &&
           glue_.basis[0].numColumns == 1 && glue_.basis[1].numColumns == 2;
  }

  /// Re-derives the amplitude layout from the basis widths and carries the
  /// values across it. The offsets are the widths' prefix sums, so nothing can
  /// drift; a layout that did not move leaves the amplitudes untouched (bcf's
  /// mid-life z swap, which must not disturb a live glue draw), and one that
  /// did keeps every block's surviving coordinates at their new offsets and
  /// enters the rest at the neutral amplitude 1.0. The per-forest priors are
  /// resized alongside, so forest f always has both.
  void rebuildAmplitudeLayout() {
    std::size_t numForests = glue_.basis.size();
    glue_.prior.resize(numForests);

    std::vector<std::size_t> offset(numForests + 1);
    offset[0] = 0;
    for (std::size_t f = 0; f < numForests; ++f)
      offset[f + 1] = offset[f] + glue_.basis[f].numColumns;
    if (offset == glue_.amplitudeOffset) return;

    std::vector<double> amplitudes(offset.back(), 1.0);
    std::size_t shared = std::min(numForests, glue_.amplitudeOffset.size() - 1);
    for (std::size_t f = 0; f < shared; ++f) {
      std::size_t width =
        std::min(offset[f + 1] - offset[f], glue_.numAmplitudes(f));
      for (std::size_t j = 0; j < width; ++j)
        amplitudes[offset[f] + j] = glue_.amplitudesOf(f)[j];
    }
    glue_.amplitudeOffset.swap(offset);
    glue_.amplitudes.swap(amplitudes);
  }

  const ColumnStore& data_;
  const std::size_t numForests_;
  const bool generalAmplitudeDraw_;
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
///
/// An optional n x K CATEGORY OFFSET makes the latent f_ik + o_ik. Everything
/// above reads it through rawFits: the margins become C_ik = log sum_{j != k}
/// exp(f_ij + o_ij), the working response gains a - o_if so forest f still
/// estimates its own part, and the reported softmax is taken over the offset
/// fits. It is not the response model's offset, which is added to every
/// reported channel AFTER the forests are combined - past the softmax that
/// would be the wrong side of the nonlinearity, and a flat per-observation
/// shift is the softmax's null direction in any case. The level-centering move
/// is unaffected: its shift is derived from the leaf prior alone, and a common
/// shift of all K is that same null direction. The test rows carry their own
/// nTest x K offset, entering the reported test blend where the train offset
/// enters the reported train blend; it is a separate object, since the test
/// rows are other rows.
template <IntegrableLeafModel L, typename ResidT = double>
struct MultinomialForestCombiner : ForestCombiner<L, ResidT> {
  static_assert(!L::hasVectorParams && !L::hasFunctionParams,
                "multinomial is a constant-leaf model");

  MultinomialForestCombiner(const ColumnStore& data, const MultinomialSpec& spec)
      : data_(data), numCategories_(spec.numCategories), counts_(spec.counts),
        trials_(spec.trials), offset_(spec.offset) {
    std::size_t n = data_.numObservations;
    // the raw slab exists only under an offset; off it every read resolves to
    // the forests' own totalFits and this stays empty
    if (offset_ != nullptr) raw_.resize(n * numCategories_);
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

  /// The response IS the count matrix here, so it is replaceable at fixed n and
  /// K: every quantity derived from it (the margins, omega, the working
  /// response) is per-sweep scratch the next sweep rebuilds, and the leaf
  /// calibration is the data-independent pi*sqrt(3)/sqrt(2) anchor, so nothing
  /// carries an old count forward. The trees do carry over, fitted to the
  /// previous counts, exactly as a single-forest setResponse leaves them.
  bool supportsCountsMutation() const override { return true; }

  /// Installs a replacement count matrix and its per-observation trials, both
  /// borrowed and both sized to the constructed n and K (the host validates the
  /// shape; nothing here can grow a buffer). The next sweep re-reads them: the
  /// trials drive the PG draw count and the counts the working response, and no
  /// scratch is invalidated - drawForestGlue's f == 0 branch rebuilds the whole
  /// suffix/prefix mix, and every omega column is written before it is read.
  /// omega_ is deliberately NOT re-seeded: it holds no count-derived value
  /// between sweeps. A later setState restores trees against WHATEVER counts
  /// the sampler holds then, as a single-forest restore does against the
  /// current y - the counts are data and ride no wire block.
  void setCounts(const int* counts, const int* trials) override {
    counts_ = counts;
    trials_ = trials;
  }

  /// Installs a replacement n x K category offset, borrowed and sized to the
  /// constructed n and K, or clears it at a null pointer. The latent becomes
  /// f_ik + o_ik: the margins C_ik are formed over the OFFSET fits, the working
  /// response subtracts o_if back off (forest f estimates f_if, not the sum),
  /// and the reported softmax is taken over the offset fits, so the offset
  /// shifts the model without ever entering a leaf value. Nothing derived is
  /// cached across the call - the raw slab is rematerialized in full at every
  /// sweep entry and at every reporting point - so a mid-life install or clear
  /// takes effect on the next read with no refresh hook.
  void setCategoryOffset(const double* offset) override {
    offset_ = offset;
    // grows once and is kept: a caller clearing and reinstalling an offset
    // should not pay for the allocation twice, and an unused slab costs nK
    // doubles on a sampler that had one
    if (offset_ != nullptr)
      raw_.resize(data_.numObservations * numCategories_);
  }

  /// Installs a replacement nTest x K test offset, borrowed and sized to the
  /// CURRENT test store, or clears it at a null pointer. It shifts the reported
  /// test probabilities to softmax(f_test + o_test) and nothing else: the test
  /// fits enter no likelihood, so this touches no draw and no working response,
  /// and a sampler that carries only this one is bitwise a sampler with none on
  /// every train channel.
  ///
  /// It is NOT the train offset restricted to the test rows: the test rows are
  /// other rows, and no shape coincidence between them makes one the other. The
  /// host must therefore not resize the test store under an installed offset -
  /// the two would silently stop describing the same rows - which is why the
  /// bridge refuses a test-predictor install while one is here.
  void setCategoryTestOffset(const double* offset) override {
    testOffset_ = offset;
  }

  /// Installs the GLOBAL active-row mask, or clears it at a null pointer. An
  /// inactive row leaves the softmax likelihood entirely: its K interleaved
  /// Polya-Gamma draws are SKIPPED in drawForestGlue and its composed precision
  /// is zero in every category's formForestResponse, so it enters no leaf
  /// sufficient statistic, no branch score and no leaf draw of any forest. It
  /// keeps its leaf occupancy and its reported softmax probabilities, as every
  /// other family's inactive rows keep their fits.
  ///
  /// PER-FOREST masking is refused, permanently and on model grounds rather
  /// than for want of an implementation: category f's margin C_if is a
  /// log-sum-exp over the OTHER K-1 forests, so a row absent from category f's
  /// forest is still in every other category's likelihood, and "row i is out of
  /// category f only" restricts no likelihood at all. Only the global mask is
  /// well-posed here, which is why this takes no forest index.
  ///
  /// The mask is length-n and n is fixed at creation for a multi-forest chain
  /// (whole-data replacement is refused outright), so an installed mask cannot
  /// go stale. A count or offset swap leaves it standing: it names rows, not
  /// responses.
  void setActiveRows(const double* active) override {
    if (active == nullptr) activeRows_.clear();
    else activeRows_.assign(active, active + data_.numObservations);
  }

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
      // Under an offset the whole raw slab is rebuilt first: this is one of the
      // two points that dominate every write to totalFits made outside a sweep
      // (a predictor mutation's tree revalidation), so nothing carries a fit
      // from before it.
      if (offset_ != nullptr) materializeRawFits(forests);
      double* top = suffix_.data() + (K - 1) * n;
      for (std::size_t i = 0; i < n; ++i) top[i] = -HUGE_VAL;
      for (std::size_t g = K - 1; g-- > 0; ) {
        const double* next = rawFits(g + 1, forests);
        const double* above = suffix_.data() + (g + 1) * n;
        double* here = suffix_.data() + g * n;
        for (std::size_t i = 0; i < n; ++i)
          here[i] = logSumExp2(next[i], above[i]);
      }
      for (std::size_t i = 0; i < n; ++i) prefix_[i] = -HUGE_VAL;
      for (std::size_t g = 0; g < f; ++g) {
        const double* below = rawFits(g, forests);
        for (std::size_t i = 0; i < n; ++i)
          prefix_[i] = logSumExp2(prefix_[i], below[i]);
      }
    } else {
      // In-order continuation: category f - 1's tree update landed since its
      // glue call, so its fit is now NEW; fold it into the running prefix (the
      // interleaving's below-f mix). suffix_[f] still holds the OLD above-f mix.
      // Refreshing only column f - 1 is exact WITHIN a sweep: between this call
      // and the previous one the sole totalFits write is that forest's own
      // finalizeTotalFits. It is safe to be lazy here only because both
      // reporting paths rematerialize in full.
      if (offset_ != nullptr) refreshRawColumn(f - 1, forests);
      const double* prev = rawFits(f - 1, forests);
      for (std::size_t i = 0; i < n; ++i)
        prefix_[i] = logSumExp2(prefix_[i], prev[i]);
    }
    lastF_ = f;

    const double* fFits = rawFits(f, forests);
    const double* suffix = suffix_.data() + f * n;
    double* omega = omega_.data() + f * n;
    const double* active = activeRows_.empty() ? nullptr : activeRows_.data();
    for (std::size_t i = 0; i < n; ++i) {
      double margin = logSumExp2(prefix_[i], suffix[i]);
      margins_[i] = margin;
      // An inactive row's K latents are SKIPPED, not drawn and discarded: the
      // Polya-Gamma sampler is a rejection sampler whose consumption depends on
      // psi, so a discard would desynchronize the stream from the compacted
      // model this channel is meant to reproduce. omega keeps its last drawn
      // (or cold-started) value, which is STALE and strictly positive -
      // formForestResponse divides by it, so the zero goes into the composed
      // precision and never in here.
      if (active != nullptr && active[i] == 0.0) continue;
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
  ///
  /// INVARIANT: drawForestGlue(f) immediately precedes this call, for the same
  /// f, on the same fits. margins_ and omega_'s f-th column are that call's
  /// handoff, and nothing else writes them, so a caller that reordered the pair
  /// would form category f's response against another category's margin. Both
  /// sweep loops (Chain::run and growForestFromRoot) fire the pair adjacently.
  /// Not asserted: drawForestGlue sets lastF_ = f unconditionally, so a
  /// lastF_ == f check here is tautological.
  ///
  /// Under a category offset the response carries an extra - o_if: the latent
  /// is f_if + o_if, so the forest is asked for the residual part. This is the
  /// one train-side reader that cannot go through the offset fits, since it
  /// SUBTRACTS the offset rather than adding it. The branch is hoisted out of
  /// the loop, so INVARIANT: the offset-free arm carries neither a per-row test
  /// nor an offset term - its arithmetic is the plain PG reduction, with no
  /// residual dependence on offset_. The subtraction also keeps totalFits
  /// offset-free (finalizeTotalFits sums this response's own tree fits), which
  /// is the invariant the offset fits are built on.
  ForestResponse formForestResponse(std::size_t f,
      const std::vector<Forest<L, ResidT>>& forests, const double* /*y*/,
      const double* /*w*/) override {
    (void) forests;
    std::size_t n = data_.numObservations;
    const double* omega = omega_.data() + f * n;
    if (offset_ == nullptr) {
      for (std::size_t i = 0; i < n; ++i) {
        double yif = static_cast<double>(counts_[f * n + i]);
        forestResponse_[i] = (yif - trials_[i] * 0.5) / omega[i] + margins_[i];
        forestWeights_[i] = omega[i];
      }
    } else {
      const double* o = offset_ + f * n;
      for (std::size_t i = 0; i < n; ++i) {
        double yif = static_cast<double>(counts_[f * n + i]);
        forestResponse_[i] =
          (yif - trials_[i] * 0.5) / omega[i] + margins_[i] - o[i];
        forestWeights_[i] = omega[i];
      }
    }
    // The mask composes into the PRECISION, in a pass of its own so the
    // unmasked loops above stay the statement they were: a_i = 0 zeroes the
    // row's precision in every category, which is what drops it from every
    // leaf sufficient statistic, branch score and leaf draw of forest f. The
    // response is left at its stale-omega value rather than zeroed - nothing
    // reads a zero-precision row's response, and a NaN there would propagate.
    if (!activeRows_.empty())
      for (std::size_t i = 0; i < n; ++i) forestWeights_[i] *= activeRows_[i];
    return {forestResponse_.data(), forestWeights_.data()};
  }

  /// Category f's veto precisions: omega_if under the active-row mask, the same
  /// product the response half composes. The chain's w is ignored here as it is
  /// there - this family carries its precisions in omega, which is strictly
  /// positive whether drawn or cold-started, so the mask is the only zero.
  const double* formForestVetoWeights(std::size_t f,
                                      const double* /*w*/) override {
    std::size_t n = data_.numObservations;
    const double* omega = omega_.data() + f * n;
    for (std::size_t i = 0; i < n; ++i)
      forestWeights_[i] =
        activeRows_.empty() ? omega[i] : omega[i] * activeRows_[i];
    return forestWeights_.data();
  }

  /// The K softmax probabilities per observation, location-major (channel k at
  /// combined_[k*n + i]); the reported training output. Log-sum-exp-safe.
  ///
  /// Under an offset the whole raw slab is rebuilt here, not refreshed
  /// per-column. This is a reporting point, and it is reached from three places
  /// per recorded sweep - after the forest loop, inside storeSample (which runs
  /// after the level-centering move has shifted every totalFits) and in
  /// grow-from-root - so a column left at the previous sweep's fit would be
  /// reported as though it were current. The full pass also absorbs any
  /// totalFits written between sweeps, outside the combiner entirely, by a
  /// predictor mutation's tree revalidation.
  const double* combinedFits(const std::vector<Forest<L, ResidT>>& forests) override {
    if (offset_ != nullptr) materializeRawFits(forests);
    return blendSoftmax(data_.numObservations,
        [&](std::size_t k) { return rawFits(k, forests); },
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
  ///
  /// Under a TEST offset the blend reads the offset test fits instead, through
  /// the same rawTestFits/rawFits split the train side uses. INVARIANT: off a
  /// test offset the pointer the blend gathers IS forest k's own
  /// totalTestFits - rawTest_ is never sized or written, so the reported test
  /// channel is the forests' fits themselves rather than a copy of them.
  /// Rematerializing here rather than at the install is what keeps the two
  /// sides symmetric: totalTestFits is rewritten by every sweep and by a
  /// predictor mutation, and this is the only point that reports it.
  const double* combinedTestFits(
      const std::vector<Forest<L, ResidT>>& forests) override {
    std::size_t nTest = data_.numTestObservations;
    combinedTest_.resize(nTest * numCategories_);
    if (testOffset_ != nullptr) {
      rawTest_.resize(nTest * numCategories_);
      materializeRawTestFits(forests);
    }
    return blendSoftmax(nTest,
        [&](std::size_t k) { return rawTestFits(k, forests); },
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
  /// The prior lives on LEAF VALUES, not on the total fits: the per-leaf sd is
  /// s_k = leaf.scale_k / k_k, and leaf.scale ALREADY carries the 1/sqrt(m_k)
  /// (ConstantGaussianLeaf::scale is nodeScale/sqrt(numTrees)), so dividing by
  /// sqrt(m_k) again here would count that factor twice and draw the shift from
  /// a variance m_k too small. With L_k and S_k the count and value sum of
  /// forest k's OCCUPIED leaves over ALL its trees, the exact conditional of
  /// the uniformly absorbed shift is
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
  /// scope here, so the saved (flattened) tree leaves are not touched.
  /// Returns 1.0 (no multiplicative scale; the return feeds only BCF's test).
  double afterCombine(std::vector<Forest<L, ResidT>>& forests, bool /*record*/,
                      std::size_t /*sampleNum*/, ext_rng* rng) override {
    std::size_t n = data_.numObservations;
    double prec = 0.0, num = 0.0;
    for (std::size_t k = 0; k < numCategories_; ++k) {
      Forest<L, ResidT>& forest = forests[k];
      double m = static_cast<double>(forest.numTrees);
      if (!(m > 0.0)) continue;
      double s = forest.leaf.scale / forest.k;
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
  /// The single definition of "the per-observation fit the softmax sees" for
  /// category k: forest k's own totalFits off an offset, and the offset column
  /// f_k + o_k under one. Every train-side reader goes through it - the suffix
  /// fold, the prefix seed and its continuation, the drawn category's own fits,
  /// and the reported blend - so there is one place where the offset enters the
  /// latent. INVARIANT: off an offset this returns forest k's own totalFits
  /// buffer, so no train-side reader sees a copy and raw_ stays empty; the
  /// offset costs an offset-free run neither an allocation nor an addition.
  const double* rawFits(std::size_t k,
                        const std::vector<Forest<L, ResidT>>& forests) const {
    if (offset_ == nullptr) return forests[k].totalFits.data();
    return raw_.data() + k * data_.numObservations;
  }

  /// Rewrites category k's offset fits from the forest's current totalFits.
  /// Never called off an offset.
  void refreshRawColumn(std::size_t k,
                        const std::vector<Forest<L, ResidT>>& forests) {
    std::size_t n = data_.numObservations;
    const double* fits = forests[k].totalFits.data();
    const double* o = offset_ + k * n;
    double* out = raw_.data() + k * n;
    for (std::size_t i = 0; i < n; ++i) out[i] = fits[i] + o[i];
  }

  /// Rewrites all K columns. The two callers (a fresh sweep entry and the
  /// reported blend) are what make the in-sweep per-column refresh safe.
  void materializeRawFits(const std::vector<Forest<L, ResidT>>& forests) {
    for (std::size_t k = 0; k < numCategories_; ++k)
      refreshRawColumn(k, forests);
  }

  /// rawFits' test twin: forest k's own totalTestFits off a test offset, and
  /// the offset column f_test_k + o_test_k under one. The single definition of
  /// "the test fit the reported softmax sees", on rawFits' invariant: off a
  /// test offset it returns the forest's own buffer, never a copy of it.
  const double* rawTestFits(
      std::size_t k, const std::vector<Forest<L, ResidT>>& forests) const {
    if (testOffset_ == nullptr) return forests[k].totalTestFits.data();
    return rawTest_.data() + k * data_.numTestObservations;
  }

  /// Rewrites all K test columns from the forests' current totalTestFits. The
  /// one caller is the reported test blend, so nothing here is ever stale;
  /// never called off a test offset.
  void materializeRawTestFits(const std::vector<Forest<L, ResidT>>& forests) {
    std::size_t nTest = data_.numTestObservations;
    for (std::size_t k = 0; k < numCategories_; ++k) {
      const double* fits = forests[k].totalTestFits.data();
      const double* o = testOffset_ + k * nTest;
      double* out = rawTest_.data() + k * nTest;
      for (std::size_t i = 0; i < nTest; ++i) out[i] = fits[i] + o[i];
    }
  }

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
  const double* offset_; // borrowed n x K category offset, same layout; or null
  const double* testOffset_ = nullptr;  // borrowed nTest x K test one, or null
  std::vector<double> raw_;  // n x K fits + offset; empty and unread off one
  std::vector<double> rawTest_;  // nTest x K test fits + test offset; likewise
  std::vector<double> omega_;          // n x K, category-major (column k at k*n)
  std::vector<double> margins_;        // n; the current forest's C_if handoff
  std::vector<double> suffix_;         // n x K; per-sweep LSE over j > f, old fits
  std::vector<double> prefix_;         // n; running LSE over j < f, new fits
  std::size_t lastF_;                  // prior glue category; detects a fresh sweep
  std::vector<double> combined_;       // n x K softmax probabilities
  std::vector<double> combinedTest_;   // nTest x K softmax test probabilities
  std::vector<double> forestResponse_, forestWeights_;  // n each
  std::vector<double> activeRows_;     // the global 0/1 mask; empty when none
};

}  // namespace bartcore

#endif  // BARTCORE_COMBINER_HPP
