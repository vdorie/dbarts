#include "R_interface_rbart.hpp"

#include <cstddef> // size_t

#include <external/Rinternals.h> // SEXP

#include <external/stats.h> // ext_cumulativeProbabilityOfNormal

#include <rc/util.h>

using std::size_t;

extern "C" {

/// Posterior mean of a random-intercept fit: for each observation, the BART
/// draws plus the intercept its own group drew, averaged over draws.
/// \p groupByExpr carries one 1-based index into ranef's group dimension per
/// observation, and both loops below index ranef through it with no fallback,
/// so it must be exactly as long as there are observations and hold no NA and
/// nothing outside [1, numGroups].
SEXP rbart_getFitted(SEXP yhatExpr, SEXP ranefExpr, SEXP groupByExpr, SEXP responseIsBinaryExpr) {
  SEXP ranefDimsExpr = PROTECT(rc_getDims(ranefExpr));
  SEXP yhatDimsExpr = PROTECT(rc_getDims(yhatExpr));

  const int* ranefDims = INTEGER(ranefDimsExpr);
  const int* yhatDims = INTEGER(yhatDimsExpr);

  const double* yhat = REAL(yhatExpr);
  const double* ranef = REAL(ranefExpr);

  const int* groupBy = INTEGER(groupByExpr);

  size_t n;
  size_t numTotalSamples;
  size_t numGroups;

  bool responseIsBinary = INTEGER(responseIsBinaryExpr)[0] != 0;

  if (rc_getLength(yhatDimsExpr) == 2) {
    // chains were combined or only one exists;
    // ranef: (n.chains * n.samples) x n.groups
    // yhat:  (n.chains * n.samples) x n.obs

    n = static_cast<size_t>(yhatDims[1]);
    numTotalSamples = static_cast<size_t>(ranefDims[0]);
    numGroups = static_cast<size_t>(ranefDims[1]);
  } else {
    n = static_cast<size_t>(yhatDims[2]);
    numTotalSamples = static_cast<size_t>(ranefDims[0] * ranefDims[1]);
    numGroups = static_cast<size_t>(ranefDims[2]);
  }

  // Enforced in one pass rather than per element inside the loops below;
  // NA_INTEGER is INT_MIN and so fails the lower bound with no test of its
  // own. The R callers refuse an unmatched group label first - this backstops
  // a direct .Call.
  if (rc_getLength(groupByExpr) != n)
    Rf_error("group index length must match the number of observations");
  for (size_t i = 0; i < n; ++i)
    if (groupBy[i] < 1 || static_cast<size_t>(groupBy[i]) > numGroups)
      Rf_error(
        "group index must be a non-NA index into the 'ranef' group dimension");

  SEXP resultExpr = PROTECT(rc_newReal(n));
  double* result = REAL(resultExpr);
  
  if (responseIsBinary) {
    for (size_t i = 0; i < n; ++i) {
      result[i] = 0.0;
      for (size_t j = 0; j < numTotalSamples; ++j) {
        result[i] += ext_cumulativeProbabilityOfNormal(
          yhat[j + i * numTotalSamples] + ranef[j + (groupBy[i] - 1) * numTotalSamples],
          0.0,
          1.0);
      }
      result[i] /= static_cast<double>(numTotalSamples);
    }
  } else {
    for (size_t i = 0; i < n; ++i) {
      result[i] = 0.0;
      for (size_t j = 0; j < numTotalSamples; ++j) {
        result[i] += yhat[j + i * numTotalSamples] + ranef[j + (groupBy[i] - 1) * numTotalSamples];
      }
      result[i] /= static_cast<double>(numTotalSamples);
    }
  }

  UNPROTECT(3);

  return resultExpr;
}
} // extern "C"

