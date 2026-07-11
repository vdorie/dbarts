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
  // usual creation path leaves these empty and borrows from R instead
  std::vector<double> ownedResponse, ownedWeights, ownedOffset,
                      ownedTestOffset;

  // BCF holds an owned copy of the 0/1 treatment the chains borrow, so
  // setTreatment installs a replacement without a protection slot
  std::vector<double> ownedTreatment;

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
/// ignored); newdata, when non-null, replays its newdataNumRows rows for
/// the n column. trainingReplay supplies the training predictors the saved
/// trees replay when no newdata is given (the engine keeps no matrix): the
/// caller passes the current data@x, or null to leave saved-tree counts
/// unpopulated. caller labels range errors.
SEXP getTrees(bartcore::SamplerBase& sampler, const std::size_t* chainIndices,
              std::size_t numChainIndices, const std::size_t* sampleIndices,
              std::size_t numSampleIndices, const std::size_t* treeIndices,
              std::size_t numTreeIndices, bool useLiveTrees,
              const double* newdata, std::size_t newdataNumRows,
              const double* trainingReplay, std::size_t trainingReplayNumRows,
              const char* caller);

/// Errors unless every replacement value for a categorical column is an
/// existing category code; ordinal columns pass through.
void validateColumnValues(const bartcore::ColumnStore& store,
                          std::size_t column, const double* values,
                          std::size_t numValues);

} // namespace bartcore_bridge

#endif // R_INTERFACE_BARTCORE_COMMON_HPP
