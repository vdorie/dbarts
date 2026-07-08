// .Call shim exposing the bartcore sampler to R for the statistical
// equivalence gate (benchmarks/R/equivalence.R). Not part of the package
// build; compiled standalone via R CMD SHLIB (see benchmarks/R/bartcore-shim.R).

#include <cstddef>
#include <cstdint>
#include <vector>

// R's remapped api (length, error, ...) collides with the standard library
// headers bartcore pulls in
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

#include <misc/simd.h>
#include <external/random.h>

#include <bartcore/bartcore.hpp>

using std::size_t;

namespace {

const double* optionalReal(SEXP x) {
  return Rf_isNull(x) ? nullptr : REAL(x);
}

}  // namespace

extern "C" SEXP bartcore_fit(SEXP xExpr, SEXP yExpr, SEXP xTestExpr,
                             SEXP weightsExpr, SEXP offsetExpr, SEXP binaryExpr,
                             SEXP sigestExpr, SEXP sigdfExpr, SEXP sigScaleExpr,
                             SEXP ntreeExpr, SEXP kExpr, SEXP nodeScaleExpr,
                             SEXP numcutExpr, SEXP usequantsExpr,
                             SEXP ndpostExpr, SEXP nskipExpr,
                             SEXP splitprobsExpr, SEXP dartExpr) {
  static bool simdInitialized = false;
  if (!simdInitialized) {
    misc_simd_init();
    simdInitialized = true;
  }

  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  size_t n = static_cast<size_t>(INTEGER(dims)[0]);
  size_t p = static_cast<size_t>(INTEGER(dims)[1]);

  size_t numTest = 0;
  const double* xTest = nullptr;
  if (!Rf_isNull(xTestExpr)) {
    SEXP testDims = Rf_getAttrib(xTestExpr, R_DimSymbol);
    numTest = static_cast<size_t>(INTEGER(testDims)[0]);
    xTest = REAL(xTestExpr);
  }

  bool binary = Rf_asLogical(binaryExpr) == TRUE;
  size_t numTrees = static_cast<size_t>(Rf_asInteger(ntreeExpr));
  size_t ndpost = static_cast<size_t>(Rf_asInteger(ndpostExpr));
  size_t nskip = static_cast<size_t>(Rf_asInteger(nskipExpr));

  bartcore::SamplerOptions options;
  options.numTrees = numTrees;
  options.k = Rf_asReal(kExpr);
  options.nodeScale = Rf_asReal(nodeScaleExpr);
  options.maxNumCuts = static_cast<std::uint32_t>(Rf_asInteger(numcutExpr));
  options.useQuantiles = Rf_asLogical(usequantsExpr) == TRUE;
  options.splitProbabilities = optionalReal(splitprobsExpr);
  if (!Rf_isNull(dartExpr)) {
    // dartExpr: numeric c(alpha, updateAlpha, a, b, rho)
    const double* dartParams = REAL(dartExpr);
    options.useDart = true;
    options.dart.alpha = dartParams[0];
    options.dart.updateAlpha = dartParams[1] != 0.0;
    options.dart.betaA = dartParams[2];
    options.dart.betaB = dartParams[3];
    options.dart.rho = dartParams[4];
    options.dart.updateDelay = nskip / 2;  // BART-package startdart convention
  }

  GetRNGstate();
  ext_rng* rng = ext_rng_createDefault(true);
  if (rng == nullptr) {
    PutRNGstate();
    Rf_error("could not create rng");
  }

  bartcore::ConstantLeafSampler sampler(
    REAL(xExpr), REAL(yExpr), n, p, optionalReal(weightsExpr),
    optionalReal(offsetExpr),
    binary ? bartcore::ResponseFamily::probit
           : bartcore::ResponseFamily::gaussian,
    Rf_asReal(sigestExpr), Rf_asReal(sigdfExpr), Rf_asReal(sigScaleExpr),
    options, &rng);
  if (numTest > 0) sampler.setTestPredictors(xTest, numTest);

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP sigmaExpr = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(ndpost)));
  SEXP trainExpr = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n),
                                          static_cast<int>(ndpost)));
  SEXP testExpr = numTest > 0
    ? PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numTest),
                             static_cast<int>(ndpost)))
    : PROTECT(R_NilValue);
  SEXP varcountExpr = PROTECT(Rf_allocMatrix(INTSXP, static_cast<int>(p),
                                             static_cast<int>(ndpost)));

  std::vector<std::uint32_t> variableCounts(p * ndpost);

  bartcore::Results results;
  results.sigma = REAL(sigmaExpr);
  results.trainingFits = REAL(trainExpr);
  results.testFits = numTest > 0 ? REAL(testExpr) : nullptr;
  results.variableCounts = variableCounts.data();

  sampler.run(nskip, ndpost, results);

  ext_rng_destroy(rng);
  PutRNGstate();

  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < p * ndpost; ++i)
    varcountOut[i] = static_cast<int>(variableCounts[i]);

  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, trainExpr);
  SET_VECTOR_ELT(resultExpr, 2, testExpr);
  SET_VECTOR_ELT(resultExpr, 3, varcountExpr);

  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("yhat.train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("yhat.test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(6);
  return resultExpr;
}
