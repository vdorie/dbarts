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
/// treatment forest and glue, z the 0/1 treatment. Gaussian only; raises R
/// errors otherwise.
BartcoreHolder* createBCFHolder(SEXP controlExpr, SEXP modelExpr,
                                SEXP dataExpr, SEXP zExpr, SEXP bcfParamsExpr);

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

/// Errors on a sampler whose residual sd is structurally pinned - a family
/// that fixes it by definition (probit, logistic, multinomial, ordinal,
/// nbinom) or a heteroscedastic variance forest - where a change would
/// silently rescale every leaf posterior precision. caller labels the error.
void refusePinnedSigmaChange(const bartcore::SamplerBase& sampler,
                             const char* caller);

/// Errors when the sampler holds a variance forest: it lives outside
/// forests_, which the transactional revalidate/rebuild helpers loop, so an
/// accepted transactional change would leave s^2(x) routing observations by
/// the predictors it was built with. caller labels the error.
void refuseVarianceForestPredictorMutation(
    const bartcore::SamplerBase& sampler, const char* caller);

/// Errors on the transactional predictor paths (setPredictor and
/// updatePredictor without forceUpdate, and the per-observation sessions,
/// which have no force variant) for a sampler a single-forest revalidation
/// cannot safely cover: a variance forest (always, through
/// refuseVarianceForestPredictorMutation), or numForests >= 2 when
/// !forcedUpdate. caller labels the error.
void refuseMultiForestTransactionalUpdate(const bartcore::SamplerBase& sampler,
                                          const char* caller,
                                          bool forcedUpdate = false);

} // namespace bartcore_bridge

#endif // R_INTERFACE_BARTCORE_COMMON_HPP
