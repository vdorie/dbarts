#ifndef R_INTERFACE_BARTCORE_COMMON_HPP
#define R_INTERFACE_BARTCORE_COMMON_HPP

// bartcore bridge internals shared between the R interface
// (R_interface_bartcore.cpp) and the flat C API (C_interface.cpp);
// definitions live in R_interface_bartcore.cpp

#include <cstddef> // size_t
#include <cstdint> // int32_t
#include <memory>  // unique_ptr
#include <vector>

#include <external/Rinternals.h> // SEXP
#include <external/random.h>     // ext_rng

#include "bartcore/bartcore.hpp"

namespace bartcore_bridge {

struct BartcoreHolder {
  std::unique_ptr<bartcore::SamplerBase> sampler;
  std::vector<ext_rng*> rngs; // one per chain
  bool keepTrainingFits;

  // a sampler created over a data handle owns its row-sliced vectors; the
  // usual creation path leaves these empty and borrows from R instead.
  // Default-initialized so a creation site lists only the three fields above
  // and every owned buffer defaults to empty (no partial-aggregate hazard).
  std::vector<double> ownedResponse{}, ownedWeights{}, ownedOffset{},
                      ownedTestOffset{};

  // one owned per-forest observation weight vector per forest, borrowed by the
  // chains. The NESTING is load-bearing: a flat numForests * n buffer would
  // dangle every installed pointer the moment it were resized, whereas growing
  // the outer vector MOVES the inner vectors without relocating their heap
  // storage, so an already-installed pointer survives.
  std::vector<std::vector<double>> ownedForestWeights{};

  // multinomial holds an owned copy of the category-major n x K count matrix
  // and the per-observation trials the combiner borrows for the sampler's
  // lifetime (the single-trial label entry builds a one-hot counts matrix)
  std::vector<int> ownedCounts{}, ownedTrials{};

  // and, when one is installed, an owned copy of the n x K category offset the
  // combiner borrows in the same layout. Empty means none: the bridge reads
  // emptiness as the sampler's "carries a category offset" answer, since the
  // engine's own probe is a capability rather than a state.
  std::vector<double> ownedCategoryOffset{};

  // the test-side twin, nTest x K over the CURRENT test rows, likewise borrowed
  // and likewise empty for none. A separate buffer because it describes
  // different rows: nothing derives one from the other, which is why replacing
  // the test rows under one is refused rather than reinterpreted.
  std::vector<double> ownedCategoryTestOffset{};

  ~BartcoreHolder() {
    for (std::size_t c = rngs.size(); c > 0; --c)
      if (rngs[c - 1] != NULL) ext_rng_destroy(rngs[c - 1]);
  }
};

/// Parses the R specification objects, builds the per-chain rngs (drawing
/// through R's stream when unseeded), and creates the sampler; prints the
/// initial summary under a verbose control. Raises R errors on invalid
/// specifications after cleaning up.
BartcoreHolder* createHolder(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                             const char* familyName);

/// The complete sampler state as a serializable R object of class
/// "bartcoreState"; unprotected on return.
SEXP storeState(bartcore::SamplerBase& sampler);

/// Restores a state object into a sampler of matching shape; raises R
/// errors on malformed or inconsistent states. currentPredictors is the
/// call-time predictor matrix a cross-grid restore re-quantizes from (data@x,
/// or the retained creation spec's @x); null for CSC/mixed stores and for a
/// same-spec continuation, which re-quantizes nothing.
void setState(bartcore::SamplerBase& sampler, SEXP stateExpr,
              const double* currentPredictors);

/// A data.frame of tree structure over 0-based index arrays; unprotected on
/// return. Reads saved trees unless useLiveTrees (sample indices are then
/// ignored); newdata, when non-null, is a borrowed predictor view - dense or
/// CSC-backed - whose newdataNumRows rows replay for the n column. Only the
/// dense flavor is materialized, and its caller already had it.
/// trainingReplay supplies the training predictors the saved
/// trees replay when no newdata is given (the engine keeps no matrix): the
/// caller passes the current data@x, or null to leave saved-tree counts
/// unpopulated. caller labels range errors.
SEXP getTrees(bartcore::SamplerBase& sampler, const std::size_t* chainIndices,
              std::size_t numChainIndices, const std::size_t* sampleIndices,
              std::size_t numSampleIndices, const std::size_t* treeIndices,
              std::size_t numTreeIndices, bool useLiveTrees,
              const bartcore::PredictorSource* newdata,
              std::size_t newdataNumRows, const double* trainingReplay,
              std::size_t trainingReplayNumRows, std::size_t forestIndex,
              const char* caller);

/// Errors unless every replacement value for a categorical column is an
/// existing category code; ordinal columns pass through.
void validateColumnValues(const bartcore::ColumnStore& store,
                          std::size_t column, const double* values,
                          std::size_t numValues);

/// True on a multi-forest sampler (BCF, numForests >= 2), the condition
/// refuseMultiForestMutation raises on. Stated apart from the raiser because
/// the flat C API answers this refusal through a return value rather than a
/// longjmp; one predicate, two channels.
bool isMultiForest(const bartcore::SamplerBase& sampler);

/// Errors on a multi-forest sampler (BCF, numForests >= 2): a whole-data or
/// whole-model mutation (setData, setModel) rebuilds or reprices only
/// forest 0, leaving the rest against stale data or an uncalibrated prior.
/// setResponse and setWeights are governed by
/// refuseMultiForestResponseMutation below, not this guard. caller labels
/// the error.
void refuseMultiForestMutation(const bartcore::SamplerBase& sampler,
                               const char* caller);

/// The response-side conduits refuseMultiForestResponseMutation covers. They
/// carry one rule and differ only in what a refusal names and in whether the
/// conduit has a scale to pin (weights do not).
enum class ResponseConduit { response, offset, weights };

/// True where a coupling fixes every response-side conduit at creation, which
/// is refuseMultiForestResponseMutation's first arm alone: the sampler caches
/// something per forest across sweeps rather than re-deriving it. Takes the
/// shape rather than the sampler so a caller reading it for other reasons pays
/// for one shape(). The second arm - a scale update against a pinned transform
/// - is NOT covered here: it depends on the conduit and on updateScale, and it
/// stays a raise on both surfaces.
bool responseConduitIsFixed(const bartcore::SamplerShape& shape);

/// Errors on a multi-forest sampler (numForests >= 2) that either fixes
/// conduit at creation - its coupling caches something per forest across
/// sweeps rather than re-deriving it (!supportsResponseMutation) - or was asked
/// to move the response transform its per-forest leaf calibrations are stated
/// against (any updateScale but FALSE, which the scaleless weight conduit
/// ignores). The one multi-forest mutation family that is opt-in rather than
/// refused; both the R bridge and the flat C API guard with this, so the two
/// surfaces cannot state different rules. A coupling whose response is its own
/// count matrix (SamplerShape::supportsCountsMutation) is still refused here -
/// no conduit of this shape can express that response - but the refusal names
/// the channel that can. caller labels the error.
void refuseMultiForestResponseMutation(const bartcore::SamplerBase& sampler,
                                       const char* caller,
                                       ResponseConduit conduit,
                                       int updateScale);

/// Errors when a heteroscedastic sampler's response-side conduit was asked to
/// re-anchor the response transform (updateScale = TRUE): the variance forest's
/// scale leaf is calibrated once, against the transform fixed at creation, and
/// nothing re-states it, so a re-anchored transform leaves s^2(x) measured on
/// the old scale and the fit runs away with getSigmas() reading unchanged. The
/// fifth sigma door beyond the transactional ones; updateScale = FALSE pins the
/// transform and stays allowed. Both the R bridge and the flat C API guard with
/// this. caller labels the error.
void refuseVarianceForestScaleUpdate(const bartcore::SamplerBase& sampler,
                                     const char* caller,
                                     ResponseConduit conduit, int updateScale);

/// Errors when a grouped (random-intercept) sampler's response-side conduit was
/// asked to re-anchor the response transform (updateScale = TRUE) under a base
/// family that has one to re-anchor: the group effects b and their scale tau are
/// held on the base model's INTERNAL scale, and the decorator converts neither,
/// so a re-anchored transform silently restates both on the original scale while
/// sigma and the residual prior ARE converted. gaussian (Student-t included,
/// which reports gaussian) and aft carry such a transform; probit and logistic
/// have none, so a grouped binary sampler takes updateScale = TRUE as the no-op
/// it already is. updateScale = FALSE pins the transform and is the supported
/// grouped response and offset swap; only setData, which may move the group
/// structure out from under the fixed per-observation indices, stays refused
/// outright. Both the R bridge and the flat C API guard with this. caller labels
/// the error.
void refuseGroupedScaleUpdate(const bartcore::SamplerBase& sampler,
                              const char* caller, ResponseConduit conduit,
                              int updateScale);

/// Errors on a weight vector its family cannot carry: probit outright (a
/// weighted probit has no tractable latent-variable form), logistic on any
/// element that is not a positive integer, since its weights are observation
/// counts and PG(w, psi) is the sum of w PG(1, psi) draws, and gaussian on any
/// element that is not finite and non-negative, since a case weight is a
/// precision multiplier there and a negative one subtracts information from a
/// leaf's sufficient statistics. numObservations sizes the per-family scan; a
/// null weights pointer is a no-op. Called at creation and on every conduit
/// that installs weights afterwards, so one text states the rule wherever it
/// is enforced, and no caller repeats it.
void enforceBinaryWeightPolicy(bartcore::ResponseFamily family,
                               const double* weights,
                               std::size_t numObservations);

/// Errors on a post-creation case-weight change under a family that carries no
/// weights at all: probit, ordinal, aft and nbinom. gaussian and logistic are
/// accepted - a logistic swap is a model change with a defined meaning, since
/// the counts are the Polya-Gamma shape and the engine redraws the latents
/// against them - and are held to the creation policy by
/// enforceBinaryWeightPolicy above. The message names the actual family rather
/// than "a binary response", which is the only thing an aft, ordinal or nbinom
/// caller can act on. The R bridge raises with this; the flat C API answers the
/// same condition through familyCarriesNoWeights below, so the two surfaces
/// state one rule in two channels.
void refuseBinaryWeightChange(const bartcore::SamplerBase& sampler);

/// True on a family that carries no case weights at all - probit, ordinal, aft
/// and nbinom - which is exactly refuseBinaryWeightChange's condition. The
/// value rule for the families that DO carry weights (logistic counts) is
/// enforceBinaryWeightPolicy's and stays a raise.
bool familyCarriesNoWeights(const bartcore::SamplerBase& sampler);

/// Errors on a response value outside the family's support, the post-creation
/// half of the rule the R surface (R/spec.R) enforces when the sampler is
/// built: 0/1 for probit and logistic, an integer category index in [1, K] for
/// ordinal, a finite non-negative integer count no larger than 1e6 for nbinom.
/// gaussian and aft (log survival times) constrain nothing, so they pass
/// through. Memory safety, not only modelling: NBDispersionPrior::computeKernel
/// sizes its count histogram from lround(max y), which a negative element
/// underflows into a ~1.8e19 allocation and an unbounded positive one grows
/// linearly out of memory. numCategories is K for ordinal and ignored
/// otherwise; a null y is a no-op. Called at creation (so the flat C API, which
/// has no R surface ahead of it, states the same rule) and on every conduit
/// that swaps y. caller labels the error.
void validateResponseSupport(bartcore::ResponseFamily family,
                             std::size_t numCategories, const double* y,
                             std::size_t numObservations, const char* caller);

/// Errors on a multi-forest sampler (numForests >= 2) whose test fits are
/// undefined, its amplitudes having no off-sample basis to multiply: the
/// combined location is ill-defined off the training rows, so the engine would
/// fall back to the bare first forest and silently misreport. Gated on
/// testFitsAreDefined rather than the forest count, so a multi-forest model
/// whose test blend IS defined (the multinomial softmax) passes through.
/// caller labels the error.
void refuseUndefinedTestFits(const bartcore::SamplerBase& sampler,
                             const char* caller);

/// True where that blend is undefined, the condition the raiser above states.
/// The flat C API's test-side entries answer it through their return value
/// instead of unwinding, so a host that did not build the sampler can probe it.
bool testFitsAreUndefined(const bartcore::SamplerBase& sampler);

/// Errors on a saved-tree read against a store holding no RECORDED draws -
/// keepTrees on, nothing run since the last reset (creation, setTreeStorage,
/// a warm-start install). Its slots hold zero-leaf trees, which replay as a
/// legitimate constant surface, so the read is refused rather than answered.
/// A sampler without tree storage passes through: it replays its live trees.
/// caller labels the error.
void refuseEmptyTreeStore(const bartcore::SamplerBase& sampler,
                          const char* caller);

/// Errors on a sampler whose residual sd is structurally pinned - a family
/// that fixes it by definition (probit, logistic, multinomial, ordinal,
/// nbinom) or a heteroscedastic variance forest - where a change would
/// silently rescale every leaf posterior precision. caller labels the error.
void refusePinnedSigmaChange(const bartcore::SamplerBase& sampler,
                             const char* caller);

/// True where that residual sd is structurally pinned, the disjunction the
/// raiser above states. A heteroscedastic sampler answers gaussian, so this is
/// the only thing that predicts the refusal without attempting the write.
bool sigmaIsPinned(const bartcore::SamplerBase& sampler);

/// Errors when an active-row mask carries an element that is neither 0 nor 1,
/// a fractional value being a weighted likelihood the latent families have no
/// coherent form for. Recoverable - a different mask would have worked - so it
/// is an error on both surfaces rather than the engine's capability false,
/// which every surface reads as "this family implements no mask". A null
/// pointer (the clear) passes through.
void refuseNonBinaryMask(const double* active, std::size_t numObservations);

/// Errors when a CSC-backed column declares a reference level against a
/// store-ORDINAL column: a reference is the code the column's IMPLICIT rows
/// take, which only a categorical column has, so the engine would read those
/// rows as the quantized zero while the caller's own densification read them
/// as the declared level. columnSources is the view's source map (a negative
/// entry names CSC column ~v) and storeTypes is indexed by SOURCE column, so a
/// subset mutation passes the types of the columns it names; a null
/// referenceMeta (nothing declared) passes through.
void refuseCscReferenceAgainstStore(const bartcore::ColumnType* storeTypes,
                                    const std::int32_t* columnSources,
                                    std::size_t numColumns,
                                    const int* referenceMeta,
                                    std::size_t numSparseColumns);

/// Errors when a designated leaf covariate column arrives CSC-backed in a test
/// view: a leaf covariate reads contiguous raw values, which CSC storage does
/// not serve. The store entrances answer this with setTestData's false return;
/// a read-only replay builds no store, so it checks the view itself.
void refuseSparseLeafCovariate(const bartcore::SamplerShape& shape,
                               const bartcore::PredictorSource& source);

/// Errors on a test view's categorical code outside the STORE's fixed
/// category counts - the training-side bound the view's author cannot see,
/// since its own declared K is the caller's, not the sampler's. Covers a
/// dense-backed column's slice, a CSC-backed column's stored codes, and the
/// reference code its implicit rows read. Run after the view is assembled and
/// before any store change: creation-time validation is long past by then, and
/// setTestData would otherwise quantize an unbounded code into the test store.
void validateTestContainerAgainstStore(const bartcore::ColumnStore& store,
                                       const bartcore::PredictorSource& view);

/// What one per-observation augmentation step reads, from either surface. fit
/// is the location WITHOUT the offset - the engine's own totalFits convention -
/// so the linear predictor psi = fit + offset is formed inside and a null
/// offset is zero; a working response reads its latent through the same member.
/// weights are the logistic counts (null is unit), cutpoints the ordinal's
/// K - 1, and the three scalars are read only by the family whose law names
/// them.
struct AugmentationInputs {
  std::size_t numObservations;
  const double* fit;
  const double* y;
  const double* weights;
  const double* offset;
  const double* cutpoints;
  std::size_t numCutpoints;
  double sigma;
  double dispersion;
  double df;
};

/// The family an augmentation helper's name selects: the law it draws AND the
/// arm validateResponseSupport states the response's support with. "student"
/// has no member of its own - a Student-t residual distribution is a gaussian
/// response under a scale mixture - so it maps to gaussian, whose support arm
/// constrains nothing and whose law IS the mixture. Errors on any other name.
bartcore::ResponseFamily augmentationFamily(const char* name);

/// One draw per observation into result, each the law its engine site draws,
/// from a per-call generator seeded off R's own stream: main R thread only, and
/// the bracket is taken here rather than by the caller. Both surfaces reach the
/// stream through this one function, so a flat and an R call under one seed
/// agree bitwise. caller labels the error a generator that cannot be created
/// raises; every argument the law reads is the caller's to have validated.
void drawAugmentation(bartcore::ResponseFamily family,
                      const AugmentationInputs& in, double* result,
                      const char* caller);

/// The quantity a host regresses on given the drawn latent, read off the same
/// engine sites: a Polya-Gamma family divides its kappa by the drawn precision,
/// Student-t regresses on y (its latent weights the row instead), and a
/// location family on the latent - each less the offset. Draws nothing.
void computeWorkingResponse(bartcore::ResponseFamily family,
                            const AugmentationInputs& in, const double* latent,
                            double* result);

} // namespace bartcore_bridge

#endif // R_INTERFACE_BARTCORE_COMMON_HPP
