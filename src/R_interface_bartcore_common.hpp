#ifndef R_INTERFACE_BARTCORE_COMMON_HPP
#define R_INTERFACE_BARTCORE_COMMON_HPP

// bartcore bridge internals shared between the R interface
// (R_interface_bartcore.cpp) and the flat C API (C_interface.cpp);
// definitions live in R_interface_bartcore.cpp

#include <cstddef> // size_t
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

  // BCF holds an owned copy of the 0/1 treatment the chains borrow, so
  // setTreatment installs a replacement without a protection slot
  std::vector<double> ownedTreatment{};

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

/// A BCF two-forest sampler (docs/design/bcf.md): the model spec supplies the
/// prognostic forest, bcfParams (length 8: tau tree count, base, power;
/// aPriorScale; sdModerate; bPriorVariance; updateA and updateB flags) the
/// treatment forest and glue, z the 0/1 treatment. moderators holds the
/// treatment forest's optional 1-based column restriction, and the four
/// trailing expressions the per-forest interactions() and blocks() lists (mu
/// is the prognostic forest, tau the treatment forest); each may be null.
/// Gaussian only; raises R errors otherwise. The public creation route reaches
/// the same build through createHolder, which reads the same pieces off a
/// control attribute.
BartcoreHolder* createBCFHolder(SEXP controlExpr, SEXP modelExpr,
                                SEXP dataExpr, SEXP zExpr, SEXP bcfParamsExpr,
                                SEXP moderatorsExpr, SEXP muInteractionsExpr,
                                SEXP tauInteractionsExpr, SEXP muBlocksExpr,
                                SEXP tauBlocksExpr);

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

/// Warm start: seed the sampler's live forests from a donor "bartcoreState"
/// over the same predictors. samplesExpr, when non-null, maps each chain to a
/// 1-based donor-sample index; NULL spreads chains across the donor pool.
/// Raises R errors on malformed states or an incompatible donor.
void installForests(bartcore::SamplerBase& sampler, SEXP donorStateExpr,
                    SEXP samplesExpr);

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

/// Errors on a multi-forest sampler (BCF, numForests >= 2): a whole-data or
/// whole-model mutation (setData, setResponse, setWeights, setModel)
/// rebuilds or reprices only forest 0, leaving the rest against stale data
/// or an uncalibrated prior. caller labels the error.
void refuseMultiForestMutation(const bartcore::SamplerBase& sampler,
                               const char* caller);

/// The response-side conduits refuseMultiForestResponseMutation covers. They
/// carry one rule and differ only in what a refusal names and in whether the
/// conduit has a scale to pin (weights do not).
enum class ResponseConduit { response, offset, weights };

/// Errors on a multi-forest sampler (numForests >= 2) that either fixes
/// conduit at creation - its coupling caches something per forest across
/// sweeps rather than re-deriving it (!supportsResponseMutation) - or was asked
/// to move the response transform its per-forest leaf calibrations are stated
/// against (any updateScale but FALSE, which the scaleless weight conduit
/// ignores). The one multi-forest mutation family that is opt-in rather than
/// refused; both the R bridge and the flat C API guard with this, so the two
/// surfaces cannot state different rules. caller labels the error.
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

/// Errors on a response value outside the family's support, the post-creation
/// half of the rule the R surface (R/spec.R) enforces when the sampler is
/// built: 0/1 for probit and logistic, an integer category index in [1, K] for
/// ordinal, a finite non-negative integer count for nbinom. gaussian and aft
/// (log survival times) constrain nothing, so they pass through. Memory safety,
/// not only modelling: NBDispersionPrior::computeKernel sizes its count
/// histogram from lround(max y), which a negative element underflows into a
/// ~1.8e19 allocation. numCategories is K for ordinal and ignored otherwise; a
/// null y is a no-op. Called at creation (so the flat C API, which has no R
/// surface ahead of it, states the same rule) and on every conduit that swaps
/// y. caller labels the error.
void validateResponseSupport(bartcore::ResponseFamily family,
                             std::size_t numCategories, const double* y,
                             std::size_t numObservations, const char* caller);

/// Errors on a BCF sampler (numForests >= 2 and its test fits undefined),
/// whose test surface has no meaning without a test treatment vector: a blend
/// a * mu + b_z * tau is ill-defined, so the engine would fall back to the bare
/// prognostic forest and silently misreport. Gated on testFitsAreDefined rather
/// than the forest count, so a multi-forest model whose test blend IS defined
/// (the multinomial softmax) passes through. caller labels the error.
void refuseBCFTestSurface(const bartcore::SamplerBase& sampler,
                          const char* caller);

/// Errors on a sampler whose residual sd is structurally pinned - a family
/// that fixes it by definition (probit, logistic, multinomial, ordinal,
/// nbinom) or a heteroscedastic variance forest - where a change would
/// silently rescale every leaf posterior precision. caller labels the error.
void refusePinnedSigmaChange(const bartcore::SamplerBase& sampler,
                             const char* caller);

/// Errors when the sampler holds a variance forest: it lives outside
/// forests_, which the transactional revalidate/rebuild helpers loop, so an
/// accepted transactional change would leave s^2(x) routing observations by
/// the predictors it was built with. The two per-observation sessions, which
/// have no force variant, call this directly; the transactional entries reach
/// it through refuseMultiForestTransactionalUpdate. caller labels the error.
void refuseVarianceForestPredictorMutation(
    const bartcore::SamplerBase& sampler, const char* caller);

/// Errors on the whole-matrix and subset transactional predictor paths
/// (setPredictor and updatePredictor without forceUpdate) for a sampler the
/// revalidation cannot cover: a variance forest, which sits outside forests_.
/// Multi-forest samplers are covered and pass. caller labels the error.
void refuseMultiForestTransactionalUpdate(const bartcore::SamplerBase& sampler,
                                          const char* caller,
                                          bool forcedUpdate = false);

} // namespace bartcore_bridge

#endif // R_INTERFACE_BARTCORE_COMMON_HPP
