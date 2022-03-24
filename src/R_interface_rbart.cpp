#include "R_interface_rbart.hpp"

#include <cstddef> // size_t

#include <external/Rinternals.h> // SEXP

#include <external/stats.h> // ext_cumulativeProbabilityOfNormal

#include <rc/util.h>

using std::size_t;

extern "C" {

SEXP rbart_getFitted(SEXP yhatExpr, SEXP ranefExpr, SEXP groupByExpr, SEXP responseIsBinaryExpr) {
  SEXP ranefDimsExpr = rc_getDims(ranefExpr);
  SEXP yhatDimsExpr =  rc_getDims(yhatExpr);

  const int* ranefDims = INTEGER(ranefDimsExpr);
  const int* yhatDims = INTEGER(yhatDimsExpr);

  const double* yhat = REAL(yhatExpr);
  const double* ranef = REAL(ranefExpr);

  const int* groupBy = INTEGER(groupByExpr);

  size_t n;
  // size_t q;
  size_t numTotalSamples;

  bool responseIsBinary = INTEGER(responseIsBinaryExpr)[0] != 0;

  if (rc_getLength(yhatDimsExpr) == 2) {
    // chains were combined or only one exists;
    // ranef: (n.chains * n.samples) x n.groups
    // yhat:  (n.chains * n.samples) x n.obs

    // q = static_cast<size_t>(ranefDims[1]);
    n = static_cast<size_t>(yhatDims[1]);
    numTotalSamples = static_cast<size_t>(ranefDims[0]);
  } else {
    // q = static_cast<size_t>(ranefDims[2]);
    n = static_cast<size_t>(yhatDims[2]);
    numTotalSamples = static_cast<size_t>(ranefDims[0] * ranefDims[1]);
  }
  
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

  UNPROTECT(1);

  return resultExpr;
}
} // extern "C"

