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
/// errors on malformed or inconsistent states.
void setState(bartcore::SamplerBase& sampler, SEXP stateExpr);

/// A data.frame of tree structure over 0-based index arrays; unprotected on
/// return. Reads saved trees unless useLiveTrees (sample indices are then
/// ignored); newdata, when non-null, replays its newdataNumRows rows for
/// the n column. caller labels range errors.
SEXP getTrees(bartcore::SamplerBase& sampler, const std::size_t* chainIndices,
              std::size_t numChainIndices, const std::size_t* sampleIndices,
              std::size_t numSampleIndices, const std::size_t* treeIndices,
              std::size_t numTreeIndices, bool useLiveTrees,
              const double* newdata, std::size_t newdataNumRows,
              const char* caller);

/// Errors unless every replacement value for a categorical column is an
/// existing category code; ordinal columns pass through.
void validateColumnValues(const bartcore::ColumnStore& store,
                          std::size_t column, const double* values,
                          std::size_t numValues);

} // namespace bartcore_bridge

#endif // R_INTERFACE_BARTCORE_COMMON_HPP
