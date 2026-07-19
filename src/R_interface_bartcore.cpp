#include "config.hpp"
#include "R_interface_bartcore.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <external/Rinternals.h>
#include <R_ext/Random.h> // GetRNGstate, PutRNGstate
#include <R_ext/Utils.h>  // R_CheckUserInterrupt

#include <external/random.h>
#include <external/stats.h> // ext_quantileOfChiSquared and its inverse

#include <rc/bounds.h>
#include <rc/util.h>

#include "bartcore/bartcore.hpp"

#include "R_interface_bartcore_common.hpp"

using std::size_t;
using std::uint32_t;
using bartcore_bridge::BartcoreHolder;
using bartcore_bridge::validateColumnValues;

namespace {

// The external pointer's protection slot pins the vectors the sampler
// borrows, one fixed slot per borrowable so replacements do not accumulate.
// Predictors are not among them: the engine owns its quantized codes and
// borrows the raw only for the duration of a build or re-quantize call, so the
// R data object (PROT_DATA at creation, plus the live sampler$data the R
// methods hold) is the sole predictor GC anchor - no PROT_PREDICTORS slot.
enum {
  PROT_DATA = 0,
  PROT_RESPONSE,
  PROT_OFFSET,
  PROT_TEST_OFFSET,
  PROT_WEIGHTS,
  PROT_COUNT
};

void retain(SEXP ptrExpr, int slot, SEXP value) {
  SET_VECTOR_ELT(R_ExternalPtrProtected(ptrExpr), slot, value);
}

void holderFinalizer(SEXP ptrExpr) {
  BartcoreHolder* holder =
    static_cast<BartcoreHolder*>(R_ExternalPtrAddr(ptrExpr));
  if (holder == NULL) return;
  delete holder;
  R_ClearExternalPtr(ptrExpr);
}

BartcoreHolder& holderFromExpression(SEXP ptrExpr) {
  BartcoreHolder* holder =
    static_cast<BartcoreHolder*>(R_ExternalPtrAddr(ptrExpr));
  if (holder == NULL)
    Rf_error("bartcore function called on NULL external pointer");
  return *holder;
}

// validates a column-major matrix of predictors against the store: matching
// column count and representable categorical codes; returns the row count
size_t validatePredictorMatrix(const bartcore::SamplerBase& sampler,
                               SEXP xExpr, const char* caller) {
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  if (!Rf_isReal(xExpr) || Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[1]) != sampler.numPredictors())
    Rf_error("%s requires a numeric matrix with matching columns", caller);
  size_t numRows = static_cast<size_t>(INTEGER(dims)[0]);
  for (size_t j = 0; j < sampler.numPredictors(); ++j)
    validateColumnValues(sampler.data(), j, REAL(xExpr) + j * numRows,
                         numRows);
  return numRows;
}

SEXP getListElement(SEXP listExpr, const char* name) {
  SEXP namesExpr = Rf_getAttrib(listExpr, R_NamesSymbol);
  if (Rf_isNull(namesExpr)) return R_NilValue;
  for (R_xlen_t i = 0; i < Rf_xlength(listExpr); ++i)
    if (std::strcmp(CHAR(STRING_ELT(namesExpr, i)), name) == 0)
      return VECTOR_ELT(listExpr, i);
  return R_NilValue;
}

// The subset of the R specification objects (dbartsControl, dbartsData,
// dbartsModel) the engine consumes. Pointers borrow from the expressions;
// error paths may leak the parse vectors - Rf_error longjmps past
// destructors, so the leak is deliberate and bounded by the R error.

struct ParsedControl {
  bool responseIsBinary = false;
  // ordinal (cumulative-probit) response shape (docs/design/ordinal.md): the
  // ordered category count K, the third response shape beside the binary and
  // continuous ones. K >= 2 selects the ordinal family; 0 (absent) is a
  // non-ordinal response, so every existing family parses unchanged. The R
  // surface (C3) attaches the count; parseControl reads it here.
  size_t numOrdinalCategories = 0;
  // negative-binomial count response shape (docs/design/negative-binomial.md
  // section 4): the fourth response shape beside the binary, ordinal, and
  // continuous ones. countResponse marks it (a count is none of the others);
  // dispersion is the r spec, the residualDf sign convention - a positive value
  // fixes r there (an integer), a non-positive value estimates r on the grid.
  // The R surface (C3) attaches both through one control attribute; parseControl
  // reads them here. Absent (the default) leaves a non-count response, so every
  // existing family parses byte-for-byte unchanged.
  bool countResponse = false;
  double dispersion = NA_REAL;
  bool verbose = false;
  bool keepTrainingFits = true;
  bool useQuantiles = false;
  bool keepTrees = false;
  size_t defaultNumSamples = 0;
  size_t numTrees = 75;
  size_t numChains = 1;
  size_t numThreads = 1;
  uint32_t treeThinningRate = 1;
  uint32_t printEvery = 100;
  uint32_t printCutoffs = 0;
  bool haveRngSeed = false;
  std::uint_least32_t rngSeed = 0;
};

// Parsed outputs of an R-side mixed test container (a dbartsMixedMatrix as
// x.test). denseAssembly and the index vectors are owned here (the test store
// copies them at build, so they need only outlive the ensuing build call);
// the CSC pointers borrow the container's sparse part, valid while it stays
// protected. Held by value in ParsedData and in the mutation entry points'
// unwind-protected scopes, so an error jump destroys it and frees its buffers.
struct ParsedTestContainer {
  std::vector<double> denseAssembly;
  const double* mixedDenseValues = NULL;
  std::vector<std::int32_t> columnSources;
  const int* cscColumnPointers = NULL;
  const int* cscRowIndices = NULL;
  const double* cscValues = NULL;
  std::vector<bartcore::xint_t> cscReferenceCodes;
  size_t numTestObservations = 0;
};

struct ParsedData {
  const double* y = NULL;
  const double* x = NULL;
  size_t numObservations = 0;
  size_t numPredictors = 0;
  std::vector<bartcore::ColumnType> columnTypes;
  bool anyCategorical = false;
  // set when x is a dgCMatrix: x stays null and the CSC slots are borrowed
  // (docs/design/sparse-columns.md)
  bool xIsSparse = false;
  const int* cscColumnPointers = NULL;
  const int* cscRowIndices = NULL;
  const double* cscValues = NULL;
  // set when x is the R-side mixed container: x stays null, the CSC slots
  // borrow the container's sparse part, mixedDenseValues its dense part,
  // and columnSources maps each column in engine convention (a nonnegative
  // dense column or the complement of a CSC column)
  bool xIsMixed = false;
  const double* mixedDenseValues = NULL;
  std::vector<std::int32_t> columnSources;
  // per sparse column of a mixed container: the reference level's 0-based
  // code and the level count K, borrowed from the container (null unless it
  // carries the metadata). Resolved per predictor into cscReferenceCodes /
  // cscCategoryCounts once varTypes marks the CSC-backed categorical columns.
  const int* cscReferenceMeta = NULL;
  const int* cscCategoryCountMeta = NULL;
  size_t numSparseColumns = 0;
  std::vector<std::uint32_t> cscCategoryCounts;
  std::vector<bartcore::xint_t> cscReferenceCodes;
  // the dense columnar container's transiently assembled block (x points
  // into it); owned here so it lives exactly as long as the parse result
  std::vector<double> denseAssembly;
  const double* x_test = NULL;
  size_t numTestObservations = 0;
  // set when x.test is the R-side mixed container: x_test stays null and the
  // parsed container carries the transient sources, which the test store copies
  // (no holder ownership, unlike the training block the store borrows)
  bool testIsMixed = false;
  ParsedTestContainer testContainer;
  const double* weights = NULL;
  const double* offset = NULL;
  const double* testOffset = NULL;
  double sigmaEstimate = 1.0;
  std::vector<uint32_t> maxNumCuts;
};

struct ParsedModel {
  double birthOrDeathProbability = 0.5;
  double swapProbability = 0.1;
  double changeProbability = 0.4;
  double birthProbability = 0.5;
  double nodeScale = 0.5;
  double power = 2.0;
  double base = 0.95;
  const double* splitProbabilities = NULL;
  bool updateK = false;
  double k = 2.0;
  double kDf = 1.5;
  double kScale = HUGE_VAL;
  bool sigmaIsFixed = false;
  double fixedSigmaSq = 1.0; // the residual variance fixed() holds
  double sigmaDf = 3.0;
  double sigmaQuantile = 0.9;
  // the chi-squared prior's scale before anchoring to the sigma estimate
  double sigmaRawScale = 1.0;
  // Student-t residual degrees of freedom (docs/design/robust-errors.md): NaN
  // is the Gaussian default, a positive value fixes nu, 0 estimates it on the
  // grid. The R surface (C3) attaches it; absent it stays NaN.
  double residualDf = NA_REAL;
  // per-predictor monotone directions in {-1, 0, +1}, narrowed from the R
  // integer spec; empty when no constraint is declared
  std::vector<std::int8_t> monotoneDirections;
  // a linear or gp node prior's designated covariate columns (0-based);
  // empty for the constant leaf
  std::vector<size_t> leafCovariateColumns;
  // gp node priors only: selects the function-valued leaf model over the
  // linear one; lengthscales empty for the median-distance heuristic
  bool gpLeaves = false;
  std::vector<double> gpLengthscales;
  size_t gpMaxLeafSize = 256;
};

// Slot temporaries feed rc_ validators and further attribute reads, both of
// which can allocate; reprotecting each read keeps the temporary live over
// its use window without per-site protection stack churn.
#define REPROTECT_SLOT(target, parent, name, index) \
  REPROTECT((target) = Rf_getAttrib((parent), Rf_install(name)), (index))

// Rf_error longjmps past C++ destructors, so an entry point that holds owning
// C++ containers across an Rf_error-capable call (the rc_ validators error
// internally) would leak them. Running that owning scope as a heap-held
// closure under R_UnwindProtect destroys the closure - and the containers it
// captures by value - on the error jump as well as the normal return, so
// deliberate error paths free their buffers instead of leaking. This is the
// only place the longjmp constraint is stated; the wrapped scopes below just
// capture their owning locals into the closure.
template <typename Body>
SEXP unwindProtect(Body body) {
  // allocate the continuation before the closure: R_MakeUnwindCont can longjmp
  // on OOM, and a closure allocated first would leak past that jump
  SEXP continuation = PROTECT(R_MakeUnwindCont());
  Body* held = new Body(std::move(body));
  SEXP result = R_UnwindProtect(
    [](void* p) -> SEXP { return (*static_cast<Body*>(p))(); }, held,
    [](void* p, Rboolean) { delete static_cast<Body*>(p); }, held,
    continuation);
  UNPROTECT(1);
  return result;
}

void parseControl(ParsedControl& control, SEXP controlExpr) {
  SEXP slotExpr;
  PROTECT_INDEX slotIndex;
  PROTECT_WITH_INDEX(R_NilValue, &slotIndex);

  REPROTECT_SLOT(slotExpr, controlExpr, "binary", slotIndex);
  control.responseIsBinary =
    rc_getBool(slotExpr, "binary response signifier", RC_LENGTH | RC_GEQ,
               rc_asRLength(1), RC_END);

  REPROTECT_SLOT(slotExpr, controlExpr, "verbose", slotIndex);
  control.verbose = rc_getBool(slotExpr, "verbose", RC_LENGTH | RC_GEQ,
                               rc_asRLength(1), RC_END);

  REPROTECT_SLOT(slotExpr, controlExpr, "keepTrainingFits", slotIndex);
  control.keepTrainingFits =
    rc_getBool(slotExpr, "keep training fits", RC_LENGTH | RC_EQ,
               rc_asRLength(1), RC_END);

  REPROTECT_SLOT(slotExpr, controlExpr, "useQuantiles", slotIndex);
  control.useQuantiles = rc_getBool(slotExpr, "use quantiles",
                                    RC_LENGTH | RC_EQ, rc_asRLength(1),
                                    RC_END);

  REPROTECT_SLOT(slotExpr, controlExpr, "keepTrees", slotIndex);
  control.keepTrees = rc_getBool(slotExpr, "keep trees", RC_LENGTH | RC_EQ,
                                 rc_asRLength(1), RC_END);

  REPROTECT_SLOT(slotExpr, controlExpr, "n.samples", slotIndex);
  control.defaultNumSamples = static_cast<size_t>(
    rc_getInt(slotExpr, "number of samples", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END));

  REPROTECT_SLOT(slotExpr, controlExpr, "n.trees", slotIndex);
  control.numTrees = static_cast<size_t>(
    rc_getInt(slotExpr, "number of trees", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END));

  REPROTECT_SLOT(slotExpr, controlExpr, "n.chains", slotIndex);
  control.numChains = static_cast<size_t>(
    rc_getInt(slotExpr, "number of chains", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END));

  REPROTECT_SLOT(slotExpr, controlExpr, "n.threads", slotIndex);
  control.numThreads = static_cast<size_t>(
    rc_getInt(slotExpr, "number of threads", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END));

  REPROTECT_SLOT(slotExpr, controlExpr, "n.thin", slotIndex);
  control.treeThinningRate = static_cast<uint32_t>(
    rc_getInt(slotExpr, "tree thinning rate", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END));

  REPROTECT_SLOT(slotExpr, controlExpr, "printEvery", slotIndex);
  int i_temp = rc_getInt(slotExpr, "print every", RC_LENGTH | RC_EQ,
                         rc_asRLength(1), RC_VALUE | RC_GEQ, 1,
                         RC_NA | RC_YES, RC_END);
  if (i_temp != NA_INTEGER) control.printEvery = static_cast<uint32_t>(i_temp);

  REPROTECT_SLOT(slotExpr, controlExpr, "printCutoffs", slotIndex);
  i_temp = rc_getInt(slotExpr, "print cutoffs", RC_LENGTH | RC_EQ,
                     rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES,
                     RC_END);
  if (i_temp == NA_INTEGER) i_temp = 0;
  control.printCutoffs = static_cast<uint32_t>(i_temp);

  // rngKind and rngNormalKind are classic-only; the R side refuses them
  // before creation, so they are not read here
  REPROTECT_SLOT(slotExpr, controlExpr, "rngSeed", slotIndex);
  if (rc_getLength(slotExpr) != 1)
    Rf_error("slot 'rngSeed' must be of length 1");
  control.haveRngSeed = INTEGER(slotExpr)[0] != NA_INTEGER;
  if (control.haveRngSeed)
    control.rngSeed = static_cast<std::uint_least32_t>(INTEGER(slotExpr)[0]);

  // optional ordinal shape (docs/design/ordinal.md): an integer K >= 2 the R
  // surface (C3) attaches for an ordered-factor response, read raw and guarded
  // like resid.df. Absent (the default) leaves a non-ordinal response, so every
  // existing family parses byte-for-byte unchanged.
  SEXP ordinalExpr =
    Rf_getAttrib(controlExpr, Rf_install("bartcore.n.categories"));
  if (Rf_isInteger(ordinalExpr) && Rf_xlength(ordinalExpr) == 1 &&
      INTEGER(ordinalExpr)[0] >= 2)
    control.numOrdinalCategories =
      static_cast<size_t>(INTEGER(ordinalExpr)[0]);

  // optional count shape (docs/design/negative-binomial.md section 4): a length-1
  // real the R surface (C3) attaches for a count response, guarded like
  // bartcore.n.categories. Its presence marks the count shape; its value is the
  // dispersion spec (positive fixes r, non-positive estimates on the grid).
  // Absent (the default) leaves a non-count response, so every existing family
  // parses byte-for-byte unchanged.
  SEXP dispersionExpr =
    Rf_getAttrib(controlExpr, Rf_install("bartcore.dispersion"));
  if (Rf_isReal(dispersionExpr) && Rf_xlength(dispersionExpr) == 1) {
    control.countResponse = true;
    control.dispersion = REAL(dispersionExpr)[0];
  }

  UNPROTECT(1);
}

struct CscSlots {
  const int* pointers = NULL;
  const int* rows = NULL;
  const double* values = NULL;
  size_t numColumns = 0;
};

// Borrowed slots of a dgCMatrix; the structure is validated here because
// the engine trusts it (rows unique, ascending, and in range per column).
CscSlots parseCscMatrix(SEXP matrixExpr, size_t numObservations) {
  CscSlots result;
  SEXP dimExpr = PROTECT(Rf_getAttrib(matrixExpr, Rf_install("Dim")));
  if (!Rf_isInteger(dimExpr) || rc_getLength(dimExpr) != 2)
    Rf_error("malformed sparse predictor matrix");
  if (static_cast<size_t>(INTEGER(dimExpr)[0]) != numObservations)
    Rf_error("number of rows of 'x' must equal length of 'y'");
  result.numColumns = static_cast<size_t>(INTEGER(dimExpr)[1]);

  SEXP pointersExpr = PROTECT(Rf_getAttrib(matrixExpr, Rf_install("p")));
  SEXP rowsExpr = PROTECT(Rf_getAttrib(matrixExpr, Rf_install("i")));
  SEXP valuesExpr = PROTECT(Rf_getAttrib(matrixExpr, Rf_install("x")));
  if (!Rf_isInteger(pointersExpr) || !Rf_isInteger(rowsExpr) ||
      !Rf_isReal(valuesExpr) ||
      static_cast<size_t>(rc_getLength(pointersExpr)) !=
        result.numColumns + 1)
    Rf_error("malformed sparse predictor matrix");
  const int* pointers = INTEGER(pointersExpr);
  const int* rows = INTEGER(rowsExpr);
  size_t numNonzero = static_cast<size_t>(rc_getLength(rowsExpr));
  if (pointers[0] != 0 || pointers[result.numColumns] < 0 ||
      static_cast<size_t>(pointers[result.numColumns]) != numNonzero ||
      static_cast<size_t>(rc_getLength(valuesExpr)) != numNonzero)
    Rf_error("malformed sparse predictor matrix");
  for (size_t j = 0; j < result.numColumns; ++j) {
    if (pointers[j + 1] < pointers[j])
      Rf_error("malformed sparse predictor matrix");
    for (int k = pointers[j]; k < pointers[j + 1]; ++k) {
      if (rows[k] < 0 ||
          static_cast<size_t>(rows[k]) >= numObservations ||
          (k > pointers[j] && rows[k] <= rows[k - 1]))
        Rf_error("malformed sparse predictor matrix");
    }
  }
  result.pointers = pointers;
  result.rows = rows;
  result.values = REAL(valuesExpr);
  UNPROTECT(4);
  return result;
}

// Code one dense container column into \p target as the zero-based doubles the
// retained coded matrix held: a factor's 1-based codes become
// as.integer(column) - 1 (NA preserved), a real column is copied verbatim. Any
// other type raises \p malformedMessage (the side's predictor-vs-test wording).
void codeDenseColumn(double* target, SEXP columnExpr, size_t numRows,
                     const char* malformedMessage) {
  if (Rf_isFactor(columnExpr)) {
    const int* codes = INTEGER(columnExpr);
    for (size_t i = 0; i < numRows; ++i)
      target[i] = codes[i] == NA_INTEGER
        ? NA_REAL : static_cast<double>(codes[i] - 1);
  } else if (Rf_isReal(columnExpr)) {
    std::memcpy(target, REAL(columnExpr), numRows * sizeof(double));
  } else {
    Rf_error("%s", malformedMessage);
  }
}

// Map a mixed container's 1-based column map into engine column sources: a
// positive k names dense column k (stored as k - 1), a negative -k names CSC
// column k (stored as the engine's ~(k - 1), the negative value itself). An
// out-of-range entry raises \p malformedMessage (the side's wording).
void mapColumnSources(std::vector<std::int32_t>& out, const int* map,
                      size_t numPredictors, size_t numDenseColumns,
                      size_t numCscColumns, const char* malformedMessage) {
  out.resize(numPredictors);
  for (size_t j = 0; j < numPredictors; ++j) {
    if (map[j] > 0 && static_cast<size_t>(map[j]) <= numDenseColumns)
      out[j] = map[j] - 1;
    else if (map[j] < 0 && static_cast<size_t>(-map[j]) <= numCscColumns)
      out[j] = map[j];
    else
      Rf_error("%s", malformedMessage);
  }
}

// Assemble a mixed container's dense block and column-source map, the shared
// core of both mixed call sites: code each dense column into \p denseAssembly
// via codeDenseColumn, publish its base through \p denseValues (null if empty),
// and fill \p columnSources via mapColumnSources. \p rowLengthMessage and
// \p malformedMessage carry the side's wording. The CSC slots and reference
// resolution stay at each call site (sparse requirement / resolve order differ).
void parseMixedContainerBlock(SEXP denseExpr, const int* map,
                              size_t numPredictors, size_t numRows,
                              size_t numCscColumns,
                              std::vector<double>& denseAssembly,
                              const double*& denseValues,
                              std::vector<std::int32_t>& columnSources,
                              const char* rowLengthMessage,
                              const char* malformedMessage) {
  size_t numDenseColumns = Rf_isNull(denseExpr)
    ? 0 : static_cast<size_t>(rc_getLength(denseExpr));
  denseAssembly.resize(numDenseColumns * numRows);
  for (size_t k = 0; k < numDenseColumns; ++k) {
    SEXP columnExpr = VECTOR_ELT(denseExpr, static_cast<R_xlen_t>(k));
    if (static_cast<size_t>(rc_getLength(columnExpr)) != numRows)
      Rf_error("%s", rowLengthMessage);
    codeDenseColumn(denseAssembly.data() + k * numRows, columnExpr, numRows,
                    malformedMessage);
  }
  denseValues = numDenseColumns > 0 ? denseAssembly.data() : NULL;
  mapColumnSources(columnSources, map, numPredictors, numDenseColumns,
                   numCscColumns, malformedMessage);
}

// Resolve each CSC-backed categorical column's reference code (and, for the
// training side via non-null \p categoryCountsOut, its level count K) from the
// per-sparse-column metadata in level order. \p boundMessage covers the
// redundant source-range check (unreachable: mapColumnSources already bounds
// it) and \p referenceMessage the [0, K) validation; both carry the side's
// wording.
void resolveCscCategoricalReferences(
    const bartcore::ColumnType* columnTypes,
    const std::int32_t* columnSources, size_t numPredictors,
    const int* referenceMeta, const int* categoryCountMeta,
    size_t numSparseColumns,
    std::vector<bartcore::xint_t>& referenceCodesOut,
    std::vector<std::uint32_t>* categoryCountsOut, const char* boundMessage,
    const char* referenceMessage) {
  referenceCodesOut.assign(numPredictors, 0);
  if (categoryCountsOut != nullptr)
    categoryCountsOut->assign(numPredictors, 0);
  for (size_t j = 0; j < numPredictors; ++j) {
    if (columnTypes[j] != bartcore::ColumnType::categorical ||
        columnSources[j] >= 0)
      continue;
    size_t source = static_cast<size_t>(~columnSources[j]);
    if (source >= numSparseColumns)
      Rf_error("%s", boundMessage);
    int count = categoryCountMeta[source];
    int reference = referenceMeta[source];
    if (count <= 0 || reference == NA_INTEGER || reference < 0 ||
        reference >= count)
      Rf_error("%s", referenceMessage);
    if (categoryCountsOut != nullptr)
      (*categoryCountsOut)[j] = static_cast<std::uint32_t>(count);
    referenceCodesOut[j] = static_cast<bartcore::xint_t>(reference);
  }
}

// Parse an R-side mixed test container (a dbartsMixedMatrix as x.test) against
// the training cut grid the engine already holds: assemble the transient dense
// block (factors carrying zero-based codes), gather the CSC slices, and resolve
// each CSC-backed categorical column's reference code (the code its implicit
// rows take). numPredictors fixes the column count the container's map must
// match; columnTypes[j] marks whether predictor j is categorical. The store
// copies everything, so the outputs need only outlive the ensuing build call.
// Shared by the creation parse and the setTestPredictor/AndOffset mutations.
void parseTestContainer(ParsedTestContainer& out, SEXP containerExpr,
                        size_t numPredictors,
                        const bartcore::ColumnType* columnTypes) {
  SEXP denseExpr = PROTECT(getListElement(containerExpr, "dense"));
  SEXP sparseExpr = PROTECT(getListElement(containerExpr, "sparse"));
  SEXP mapExpr = PROTECT(getListElement(containerExpr, "map"));
  if (!Rf_isInteger(mapExpr) ||
      static_cast<size_t>(rc_getLength(mapExpr)) != numPredictors)
    Rf_error("number of columns in 'x.test' must equal that of 'x'");
  SEXP numObsExpr = getListElement(containerExpr, "numObservations");
  if (!Rf_isInteger(numObsExpr) || rc_getLength(numObsExpr) != 1 ||
      INTEGER(numObsExpr)[0] < 0)
    Rf_error("malformed mixed test container");
  size_t numTest = static_cast<size_t>(INTEGER(numObsExpr)[0]);
  out.numTestObservations = numTest;
  const int* map = INTEGER(mapExpr);

  bool hasSparse = Rf_inherits(sparseExpr, "dgCMatrix");
  if (!Rf_isNull(sparseExpr) && !hasSparse)
    Rf_error("malformed mixed test container");
  if (!Rf_isNull(denseExpr) && TYPEOF(denseExpr) != VECSXP)
    Rf_error("malformed mixed test container");
  CscSlots csc;
  if (hasSparse) csc = parseCscMatrix(sparseExpr, numTest);
  size_t numCscColumns = hasSparse ? csc.numColumns : 0;

  parseMixedContainerBlock(denseExpr, map, numPredictors, numTest,
                           numCscColumns, out.denseAssembly,
                           out.mixedDenseValues, out.columnSources,
                           "number of rows of 'x.test' columns must match",
                           "malformed mixed test container");
  if (hasSparse) {
    out.cscColumnPointers = csc.pointers;
    out.cscRowIndices = csc.rows;
    out.cscValues = csc.values;
  }

  // resolve the reference code per CSC-backed categorical test column (the code
  // its implicit rows take); the container's per-sparse-column metadata is
  // already in level order
  bool anyTestCscCategorical = false;
  for (size_t j = 0; j < numPredictors; ++j)
    if (columnTypes[j] == bartcore::ColumnType::categorical &&
        out.columnSources[j] < 0) {
      anyTestCscCategorical = true;
      break;
    }
  if (anyTestCscCategorical) {
    SEXP referenceExpr = getListElement(containerExpr, "sparseReference");
    SEXP categoryCountExpr =
      getListElement(containerExpr, "sparseCategoryCount");
    if (!Rf_isInteger(referenceExpr) || !Rf_isInteger(categoryCountExpr) ||
        static_cast<size_t>(rc_getLength(referenceExpr)) != numCscColumns ||
        static_cast<size_t>(rc_getLength(categoryCountExpr)) != numCscColumns)
      Rf_error("sparse categorical test predictor columns require reference "
               "metadata");
    resolveCscCategoricalReferences(
      columnTypes, out.columnSources.data(), numPredictors,
      INTEGER(referenceExpr), INTEGER(categoryCountExpr), numCscColumns,
      out.cscReferenceCodes, nullptr, "malformed mixed test container",
      "sparse categorical test predictor columns require reference "
      "metadata");
  }
  UNPROTECT(3);
}

// Route a parsed test container to the engine's typed test store (against the
// training cut grid, owning its raw). Returns false, store untouched, when a
// designated leaf covariate column would be CSC-backed - sparse storage serves
// no dense raw test covariate; the caller raises the leaf-covariate error.
bool installTestContainer(bartcore::SamplerBase& sampler,
                          const ParsedTestContainer& parsed) {
  return sampler.setTestData(
    parsed.mixedDenseValues, parsed.cscColumnPointers, parsed.cscRowIndices,
    parsed.cscValues, parsed.columnSources.data(),
    parsed.cscReferenceCodes.empty() ? NULL : parsed.cscReferenceCodes.data(),
    parsed.numTestObservations);
}

void parseData(ParsedData& data, SEXP dataExpr) {
  SEXP slotExpr;
  PROTECT_INDEX slotIndex;
  PROTECT_WITH_INDEX(R_NilValue, &slotIndex);

  REPROTECT_SLOT(slotExpr, dataExpr, "y", slotIndex);
  if (!Rf_isReal(slotExpr)) Rf_error("y must be of type real");
  if (rc_getLength(slotExpr) <= 0)
    Rf_error("length of y must be greater than 0");
  data.y = REAL(slotExpr);
  data.numObservations = rc_getLength(slotExpr);

  REPROTECT_SLOT(slotExpr, dataExpr, "x", slotIndex);
  if (Rf_inherits(slotExpr, "dgCMatrix")) {
    CscSlots csc = parseCscMatrix(slotExpr, data.numObservations);
    data.numPredictors = csc.numColumns;
    data.xIsSparse = true;
    data.cscColumnPointers = csc.pointers;
    data.cscRowIndices = csc.rows;
    data.cscValues = csc.values;
    data.x = NULL;
  } else if (Rf_inherits(slotExpr, "dbartsMixedMatrix")) {
    SEXP denseExpr = PROTECT(getListElement(slotExpr, "dense"));
    SEXP sparseExpr = PROTECT(getListElement(slotExpr, "sparse"));
    SEXP mapExpr = PROTECT(getListElement(slotExpr, "map"));
    if (!Rf_isInteger(mapExpr) || rc_getLength(mapExpr) == 0)
      Rf_error("malformed mixed predictor container");
    if (TYPEOF(denseExpr) == VECSXP && Rf_isNull(sparseExpr)) {
      // the dense columnar flavor (R/mixedMatrix.R): per-column vectors,
      // factors carrying their integer codes, no sparse part. Assemble the
      // transient contiguous block - the exact doubles the retained cbind
      // held - and take the plain dense path from here; build quantizes
      // into owned codes and retains nothing.
      data.numPredictors = rc_getLength(mapExpr);
      const int* map = INTEGER(mapExpr);
      size_t numDenseColumns = static_cast<size_t>(rc_getLength(denseExpr));
      data.denseAssembly.resize(data.numPredictors * data.numObservations);
      for (size_t j = 0; j < data.numPredictors; ++j) {
        if (map[j] < 1 || static_cast<size_t>(map[j]) > numDenseColumns)
          Rf_error("malformed mixed predictor container");
        SEXP columnExpr =
          VECTOR_ELT(denseExpr, static_cast<R_xlen_t>(map[j] - 1));
        if (static_cast<size_t>(rc_getLength(columnExpr)) !=
            data.numObservations)
          Rf_error("number of rows of 'x' must equal length of 'y'");
        codeDenseColumn(
          data.denseAssembly.data() + j * data.numObservations, columnExpr,
          data.numObservations, "malformed mixed predictor container");
      }
      data.x = data.denseAssembly.data();
    } else {
      // the mixed flavor: a per-column dense list (factors carrying their
      // integer codes, or NULL for no dense columns), a dgCMatrix, and a
      // 1-based map - positive k names dense column k, negative -k sparse
      // column k, the engine's ~(k - 1). Assemble the transient block - the
      // exact doubles the retained cbind held - which the holder/handle owns
      // for the store's lifetime; the store borrows dense slices of it.
      if (!Rf_inherits(sparseExpr, "dgCMatrix"))
        Rf_error("malformed mixed predictor container");
      if (!Rf_isNull(denseExpr) && TYPEOF(denseExpr) != VECSXP)
        Rf_error("malformed mixed predictor container");
      CscSlots csc = parseCscMatrix(sparseExpr, data.numObservations);

      data.numPredictors = rc_getLength(mapExpr);
      const int* map = INTEGER(mapExpr);
      parseMixedContainerBlock(
        denseExpr, map, data.numPredictors, data.numObservations,
        csc.numColumns, data.denseAssembly, data.mixedDenseValues,
        data.columnSources,
        "number of rows of 'x' must equal length of 'y'",
        "malformed mixed predictor container");
      data.xIsMixed = true;
      data.cscColumnPointers = csc.pointers;
      data.cscRowIndices = csc.rows;
      data.cscValues = csc.values;
      data.x = NULL;

      // borrow the per-sparse-column reference metadata (a CSC-backed
      // categorical column needs it); it rides the container, so stays valid
      // while dataExpr is protected
      SEXP referenceExpr = getListElement(slotExpr, "sparseReference");
      SEXP categoryCountExpr = getListElement(slotExpr, "sparseCategoryCount");
      if (Rf_isInteger(referenceExpr) && Rf_isInteger(categoryCountExpr) &&
          static_cast<size_t>(rc_getLength(referenceExpr)) ==
            csc.numColumns &&
          static_cast<size_t>(rc_getLength(categoryCountExpr)) ==
            csc.numColumns) {
        data.cscReferenceMeta = INTEGER(referenceExpr);
        data.cscCategoryCountMeta = INTEGER(categoryCountExpr);
        data.numSparseColumns = csc.numColumns;
      }
    }
    UNPROTECT(3);
  } else {
    if (!Rf_isReal(slotExpr)) Rf_error("x must be of type real");
    rc_assertDimConstraints(slotExpr, "dimensions of x", RC_LENGTH | RC_EQ,
                            rc_asRLength(2), RC_VALUE | RC_EQ,
                            static_cast<int>(data.numObservations), RC_END);
    int* dims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
    data.x = REAL(slotExpr);
    data.numPredictors = static_cast<size_t>(dims[1]);
  }

  REPROTECT_SLOT(slotExpr, dataExpr, "varTypes", slotIndex);
  rc_assertIntConstraints(slotExpr, "variable types", RC_LENGTH | RC_EQ,
                          rc_asRLength(data.numPredictors), RC_END);
  int* i_variableTypes = INTEGER(slotExpr);
  data.columnTypes.assign(data.numPredictors, bartcore::ColumnType::ordinal);
  for (size_t j = 0; j < data.numPredictors; ++j) {
    if (i_variableTypes[j] == 0) continue; // 0 encodes ordinal
    data.columnTypes[j] = bartcore::ColumnType::categorical;
    data.anyCategorical = true;
  }
  if (data.xIsSparse && data.anyCategorical)
    Rf_error("sparse predictor matrices must be entirely ordinal");
  if (data.xIsMixed && data.anyCategorical) {
    // a CSC-backed categorical column reaches the engine only with its level
    // count K and reference code, resolved here per predictor from the
    // container's per-sparse-column metadata; without it, refuse cleanly
    bool anyCscCategorical = false;
    for (size_t j = 0; j < data.numPredictors; ++j)
      if (data.columnTypes[j] == bartcore::ColumnType::categorical &&
          data.columnSources[j] < 0) {
        anyCscCategorical = true;
        break;
      }
    if (anyCscCategorical) {
      if (data.cscReferenceMeta == NULL || data.cscCategoryCountMeta == NULL)
        Rf_error("sparse categorical predictor columns require reference "
                 "metadata");
      resolveCscCategoricalReferences(
        data.columnTypes.data(), data.columnSources.data(),
        data.numPredictors, data.cscReferenceMeta, data.cscCategoryCountMeta,
        data.numSparseColumns, data.cscReferenceCodes,
        &data.cscCategoryCounts, "malformed mixed predictor container",
        "sparse categorical predictor columns require reference "
        "metadata");
    }
  }

  REPROTECT_SLOT(slotExpr, dataExpr, "x.test", slotIndex);
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.x_test = NULL;
    data.numTestObservations = 0;
  } else if (Rf_inherits(slotExpr, "dbartsMixedMatrix")) {
    // the mixed test container parses against the training cut grid; the store
    // copies its raw, so the parsed sources ride in data (freed on an error
    // jump with it) rather than being pinned by the holder
    parseTestContainer(data.testContainer, slotExpr, data.numPredictors,
                       data.columnTypes.data());
    data.testIsMixed = true;
    data.numTestObservations = data.testContainer.numTestObservations;
  } else {
    if (!Rf_isReal(slotExpr)) Rf_error("x.test must be of type real");
    rc_assertDimConstraints(slotExpr, "dimensions of x.test",
                            RC_LENGTH | RC_EQ, rc_asRLength(2), RC_NA,
                            RC_VALUE | RC_EQ,
                            static_cast<int>(data.numPredictors), RC_END);
    int* testDims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
    data.x_test = REAL(slotExpr);
    data.numTestObservations = static_cast<size_t>(testDims[0]);
  }

  REPROTECT_SLOT(slotExpr, dataExpr, "weights", slotIndex);
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.weights = NULL;
  } else {
    rc_assertDoubleConstraints(slotExpr, "weights", RC_LENGTH | RC_EQ,
                               rc_asRLength(data.numObservations), RC_END);
    data.weights = REAL(slotExpr);
  }

  REPROTECT_SLOT(slotExpr, dataExpr, "offset", slotIndex);
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.offset = NULL;
  } else {
    rc_assertDoubleConstraints(slotExpr, "offset", RC_LENGTH | RC_EQ,
                               rc_asRLength(data.numObservations), RC_END);
    data.offset = REAL(slotExpr);
  }

  REPROTECT_SLOT(slotExpr, dataExpr, "offset.test", slotIndex);
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.testOffset = NULL;
  } else {
    rc_assertDoubleConstraints(slotExpr, "test offset", RC_LENGTH | RC_EQ,
                               rc_asRLength(data.numTestObservations),
                               RC_END);
    data.testOffset = REAL(slotExpr);
  }

  REPROTECT_SLOT(slotExpr, dataExpr, "sigma", slotIndex);
  data.sigmaEstimate = rc_getDouble(
    slotExpr, "sigma estimate", RC_LENGTH | RC_EQ, rc_asRLength(1),
    RC_NA | RC_YES, RC_VALUE | RC_GT, 0.0, RC_END);

  REPROTECT_SLOT(slotExpr, dataExpr, "n.cuts", slotIndex);
  rc_assertIntConstraints(slotExpr, "maximum number of cuts",
                          RC_LENGTH | RC_EQ,
                          rc_asRLength(data.numPredictors), RC_END);
  int* i_maxNumCuts = INTEGER(slotExpr);
  data.maxNumCuts.resize(data.numPredictors);
  for (size_t j = 0; j < data.numPredictors; ++j)
    data.maxNumCuts[j] = static_cast<uint32_t>(i_maxNumCuts[j]);

  UNPROTECT(1);
}

void parseModel(ParsedModel& model, SEXP modelExpr, size_t numPredictors) {
  SEXP slotExpr;
  PROTECT_INDEX slotIndex;
  PROTECT_WITH_INDEX(R_NilValue, &slotIndex);

  REPROTECT_SLOT(slotExpr, modelExpr, "p.birth_death", slotIndex);
  // a monotone forest is birth/death-only (p.birth_death == 1); permit it
  model.birthOrDeathProbability = rc_getDouble(
    slotExpr, "probability of birth/death rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LEQ, 1.0, RC_END);

  REPROTECT_SLOT(slotExpr, modelExpr, "p.swap", slotIndex);
  model.swapProbability = rc_getDouble(
    slotExpr, "probability of swap rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  REPROTECT_SLOT(slotExpr, modelExpr, "p.change", slotIndex);
  model.changeProbability = rc_getDouble(
    slotExpr, "probability of change rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  if (std::fabs(model.birthOrDeathProbability + model.swapProbability +
                model.changeProbability - 1.0) >= 1.0e-10)
    Rf_error("rule proposal probabilities must sum to 1.0");

  REPROTECT_SLOT(slotExpr, modelExpr, "p.birth", slotIndex);
  model.birthProbability = rc_getDouble(
    slotExpr, "probability of birth in birth/death rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  REPROTECT_SLOT(slotExpr, modelExpr, "node.scale", slotIndex);
  model.nodeScale = rc_getDouble(
    slotExpr, "scale of node prior", RC_LENGTH | RC_EQ, rc_asRLength(1),
    RC_VALUE | RC_GT, 0.0, RC_END);

  // monotone (mBART) directions ride the model as a per-predictor integer
  // attribute the R surface resolves; absent or all-zero keeps the
  // unconstrained constant leaf (docs/design/monotone.md section 8)
  REPROTECT_SLOT(slotExpr, modelExpr, "monotone", slotIndex);
  if (!Rf_isNull(slotExpr) && rc_getLength(slotExpr) > 0) {
    if (static_cast<size_t>(rc_getLength(slotExpr)) != numPredictors)
      Rf_error("length of monotone directions must equal number of predictors");
    if (!Rf_isInteger(slotExpr))
      Rf_error("monotone directions must be resolved integer signs");
    const int* directions = INTEGER(slotExpr);
    model.monotoneDirections.resize(numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) {
      if (directions[j] < -1 || directions[j] > 1)
        Rf_error("monotone directions must be -1, 0, or 1");
      model.monotoneDirections[j] = static_cast<std::int8_t>(directions[j]);
    }
  }

  // linear and gp node priors designate leaf covariate columns, resolved
  // R-side to 1-based model matrix indices; every other node prior is the
  // constant leaf and carries nothing beyond node.scale/node.hyperprior.
  // gp priors add per-column lengthscales (NULL for the median-distance
  // heuristic; the R side validates and recycles) and the leaf-size cap.
  SEXP nodePriorExpr =
    PROTECT(Rf_getAttrib(modelExpr, Rf_install("node.prior")));
  bool isLinearPrior = !Rf_isNull(nodePriorExpr) &&
                       Rf_inherits(nodePriorExpr, "dbartsLinearPrior");
  bool isGPPrior = !Rf_isNull(nodePriorExpr) &&
                   Rf_inherits(nodePriorExpr, "dbartsGPPrior");
  if (isLinearPrior || isGPPrior) {
    SEXP columnsExpr =
      PROTECT(Rf_getAttrib(nodePriorExpr, Rf_install("columns")));
    if (!Rf_isInteger(columnsExpr) || Rf_xlength(columnsExpr) < 1)
      Rf_error("node prior columns must be resolved integer indices");
    R_xlen_t numColumns = Rf_xlength(columnsExpr);
    model.leafCovariateColumns.resize(static_cast<size_t>(numColumns));
    for (R_xlen_t j = 0; j < numColumns; ++j) {
      int column = INTEGER(columnsExpr)[j];
      if (column < 1 || static_cast<size_t>(column) > numPredictors)
        Rf_error("node prior column out of range");
      model.leafCovariateColumns[static_cast<size_t>(j)] =
        static_cast<size_t>(column - 1);
    }
    UNPROTECT(1);
  }
  if (isGPPrior) {
    model.gpLeaves = true;
    // a NULL slot arrives as S4's pseudo-NULL symbol, not R_NilValue, so
    // test positively for the resolved numeric vector
    SEXP lengthscaleExpr =
      PROTECT(Rf_getAttrib(nodePriorExpr, Rf_install("lengthscale")));
    if (Rf_isReal(lengthscaleExpr)) {
      if (static_cast<size_t>(Rf_xlength(lengthscaleExpr)) !=
          model.leafCovariateColumns.size())
        Rf_error("gp node prior lengthscales must be resolved per column");
      const double* lengthscales = REAL(lengthscaleExpr);
      for (size_t j = 0; j < model.leafCovariateColumns.size(); ++j)
        if (!(lengthscales[j] > 0.0))
          Rf_error("gp node prior lengthscales must be positive");
      model.gpLengthscales.assign(
        lengthscales, lengthscales + model.leafCovariateColumns.size());
    }
    SEXP maxLeafSizeExpr =
      PROTECT(Rf_getAttrib(nodePriorExpr, Rf_install("max.leaf.size")));
    int maxLeafSize = rc_getInt(
      maxLeafSizeExpr, "gp node prior maximum leaf size", RC_LENGTH | RC_EQ,
      rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END);
    model.gpMaxLeafSize = static_cast<size_t>(maxLeafSize);
    UNPROTECT(2);
  }

  SEXP priorExpr;
  PROTECT_INDEX priorIndex;
  PROTECT_WITH_INDEX(R_NilValue, &priorIndex);

  REPROTECT_SLOT(priorExpr, modelExpr, "tree.prior", priorIndex);
  REPROTECT_SLOT(slotExpr, priorExpr, "power", slotIndex);
  model.power = rc_getDouble(slotExpr, "tree prior power", RC_LENGTH | RC_EQ,
                             rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);

  REPROTECT_SLOT(slotExpr, priorExpr, "base", slotIndex);
  model.base = rc_getDouble(slotExpr, "tree prior base", RC_LENGTH | RC_EQ,
                            rc_asRLength(1), RC_VALUE | RC_GT, 0.0,
                            RC_VALUE | RC_LT, 1.0, RC_END);

  REPROTECT_SLOT(slotExpr, priorExpr, "splitProbabilities", slotIndex);
  if (rc_getLength(slotExpr) == 0) {
    model.splitProbabilities = NULL;
  } else {
    model.splitProbabilities = REAL(slotExpr);
    if (rc_getLength(slotExpr) != numPredictors)
      Rf_error("length of split probabilities must equal number of "
               "predictors");
    double totalProbability = 0.0;
    for (size_t j = 0; j < numPredictors; ++j) {
      if (model.splitProbabilities[j] < 0.0)
        Rf_error("split probabilities must be non-negative");
      totalProbability += model.splitProbabilities[j];
    }
    if (std::fabs(totalProbability - 1.0) >= 1.0e-10)
      Rf_error("split probabilities must sum to 1.0");
  }

  REPROTECT_SLOT(priorExpr, modelExpr, "node.hyperprior", priorIndex);
  const char* classStr = CHAR(STRING_ELT(rc_getClass(priorExpr), 0));
  if (std::strcmp(classStr, "dbartsChiHyperprior") == 0) {
    model.updateK = true;
    model.kDf = rc_getDouble(
      Rf_getAttrib(priorExpr, Rf_install("degreesOfFreedom")),
      "degreesOfFreedom", RC_LENGTH | RC_EQ, rc_asRLength(1),
      RC_VALUE | RC_GT, 0.0, RC_END);
    model.kScale = rc_getDouble(
      Rf_getAttrib(priorExpr, Rf_install("scale")), "scale",
      RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
  } else if (std::strcmp(classStr, "dbartsFixedHyperprior") == 0) {
    model.k = rc_getDouble(Rf_getAttrib(priorExpr, Rf_install("k")), "k",
                           RC_LENGTH | RC_EQ, rc_asRLength(1),
                           RC_VALUE | RC_GT, 0.0, RC_END);
  } else {
    Rf_error("unsupported k prior type '%s'", classStr);
  }

  REPROTECT_SLOT(priorExpr, modelExpr, "resid.prior", priorIndex);
  classStr = CHAR(STRING_ELT(rc_getClass(priorExpr), 0));
  if (std::strcmp(classStr, "dbartsChiSqPrior") == 0) {
    model.sigmaDf = rc_getDouble(
      Rf_getAttrib(priorExpr, Rf_install("df")),
      "sigma prior degrees of freedom", RC_LENGTH | RC_EQ, rc_asRLength(1),
      RC_VALUE | RC_GT, 0.0, RC_END);
    model.sigmaQuantile = rc_getDouble(
      Rf_getAttrib(priorExpr, Rf_install("quantile")), "sigma prior quantile",
      RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0,
      RC_VALUE | RC_LT, 1.0, RC_END);
    model.sigmaRawScale =
      ext_quantileOfChiSquared(1.0 - model.sigmaQuantile, model.sigmaDf) /
      model.sigmaDf;
  } else if (std::strcmp(classStr, "dbartsFixedPrior") == 0) {
    model.sigmaIsFixed = true;
    model.fixedSigmaSq = rc_getDouble(
      Rf_getAttrib(priorExpr, Rf_install("value")),
      "residual variance prior fixed value", RC_LENGTH | RC_EQ,
      rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
  } else {
    Rf_error("unsupported residual variance prior type '%s'", classStr);
  }

  // Student-t residual df (docs/design/robust-errors.md): an optional numeric
  // model attribute the R surface (C3) attaches; absent it stays NaN, the
  // Gaussian law. Read raw here - the family cross-check and sign policy live
  // in parseSamplerSpecification, once the family is known.
  SEXP residDfExpr = Rf_getAttrib(modelExpr, Rf_install("resid.df"));
  if (Rf_isReal(residDfExpr) && Rf_xlength(residDfExpr) == 1)
    model.residualDf = REAL(residDfExpr)[0];

  UNPROTECT(3);
}

// The classic engine's BARTFit::printInitialSummary, line for line, printed
// at creation under verbose; the classic version's quantile and scale lines
// reduce at creation to expressions over the raw prior scale and the
// response range, used here directly.
void printInitialSummary(const ParsedControl& control,
                         const ParsedModel& model, const ParsedData& data,
                         const bartcore::SamplerBase& sampler) {
  if (control.responseIsBinary)
    ext_printf("\nRunning BART with binary y\n\n");
  else
    ext_printf("\nRunning BART with numeric y\n\n");

  ext_printf("number of trees: %lu\n",
             static_cast<unsigned long>(control.numTrees));
  ext_printf("number of chains: %lu, default number of threads %lu\n",
             static_cast<unsigned long>(control.numChains),
             static_cast<unsigned long>(control.numThreads));
  ext_printf("tree thinning rate: %u\n", control.treeThinningRate);

  ext_printf("Prior:\n");
  if (!model.updateK) {
    ext_printf("\tk prior fixed to %f\n", model.k);
  } else {
    ext_printf("\tprior on k: chi with %f degrees of freedom and %f scale\n",
               model.kDf, model.kScale);
  }
  if (!control.responseIsBinary) {
    if (model.sigmaIsFixed) {
      // the classic engine's FixedPrior::print
      ext_printf("\tresidual variance prior fixed to %f\n",
                 model.fixedSigmaSq);
    } else {
      ext_printf("\tdegrees of freedom in sigma prior: %f\n", model.sigmaDf);
      ext_printf("\tquantile in sigma prior: %f\n", model.sigmaQuantile);
      double sigmaInternal = data.sigmaEstimate / sampler.fitScale();
      ext_printf("\tscale in sigma prior: %f\n",
                 sigmaInternal * sigmaInternal * model.sigmaRawScale);
    }
  }

  ext_printf("\tpower and base for tree prior: %f %f\n", model.power,
             model.base);
  if (model.splitProbabilities != NULL) {
    ext_printf("\ttree split probabilities: %f",
               model.splitProbabilities[0]);
    size_t printLength = 5 < data.numPredictors ? 5 : data.numPredictors;
    for (size_t i = 1; i < printLength; ++i)
      ext_printf(", %f", model.splitProbabilities[i]);
    ext_printf("\n");
  }
  ext_printf("\tuse quantiles for rule cut points: %s\n",
             control.useQuantiles ? "true" : "false");
  ext_printf(
    "\tproposal probabilities: birth/death %.2f, swap %.2f, change %.2f; "
    "birth %.2f\n",
    model.birthOrDeathProbability, model.swapProbability,
    model.changeProbability, model.birthProbability);

  ext_printf("data:\n");
  ext_printf("\tnumber of training observations: %lu\n",
             static_cast<unsigned long>(data.numObservations));
  ext_printf("\tnumber of test observations: %lu\n",
             static_cast<unsigned long>(data.numTestObservations));
  ext_printf("\tnumber of explanatory variables: %lu\n",
             static_cast<unsigned long>(data.numPredictors));
  if (!control.responseIsBinary)
    ext_printf("\tinit sigma: %f, curr sigma: %f\n", data.sigmaEstimate,
               sampler.sigma(0));
  if (data.weights != NULL) ext_printf("\tusing observation weights\n");
  ext_printf("\n");

  const bartcore::ColumnStore& store(sampler.data());
  ext_printf("Cutoff rules c in x<=c vs x>c\n");
  ext_printf("Number of cutoffs: (var: number of possible c):\n");
  for (size_t j = 0; j < store.numPredictors; ++j) {
    ext_printf("(%lu: %u) ", static_cast<unsigned long>(j + 1),
               store.numCuts[j]);
    if ((j + 1) % 5 == 0) ext_printf("\n");
  }
  ext_printf("\n");
  if (control.printCutoffs > 0) {
    ext_printf("cutoffs:\n");
    for (size_t j = 0; j < store.numPredictors; ++j) {
      if (store.types[j] == bartcore::ColumnType::categorical) continue;
      ext_printf("x(%lu) cutoffs: ", static_cast<unsigned long>(j + 1));

      size_t k;
      for (k = 0;
           k < store.numCuts[j] - 1 && k < control.printCutoffs - 1;
           ++k) {
        ext_printf("%f", store.cutPoints[j][k]);
        if ((k + 1) % 5 == 0) ext_printf("\n\t");
      }
      if (k > 2 && k == control.printCutoffs && k < store.numCuts[j] - 1)
        ext_printf("...");

      ext_printf("%f", store.cutPoints[j][store.numCuts[j] - 1]);
      ext_printf("\n");
    }
  }

  if (data.offset != NULL ||
      (data.numTestObservations > 0 && data.testOffset != NULL)) {
    ext_printf("offsets:\n");

    if (data.offset != NULL) {
      ext_printf("\treg : %.2f", data.offset[0]);
      for (size_t i = 1;
           i < (5 < data.numObservations ? 5 : data.numObservations); ++i)
        ext_printf(" %.2f", data.offset[i]);
      ext_printf("\n");
    }
    if (data.numTestObservations > 0 && data.testOffset != NULL) {
      ext_printf("\ttest: %.2f", data.testOffset[0]);
      for (size_t i = 1;
           i < (5 < data.numTestObservations ? 5 : data.numTestObservations);
           ++i)
        ext_printf(" %.2f", data.testOffset[i]);
      ext_printf("\n");
    }
  }
}

// family selects the response model, keyed on the response shape. For a binary
// response: "" or "probit" give the classic probit latents, "logistic" the
// Polya-Gamma sampler. For a K-level ordinal response (control's category
// count): "" or "ordinal" give the cumulative probit. For a continuous
// response: "" or "gaussian" fits ordinary BART, "aft" the log-normal survival
// model (censored latents via a control attribute). Each shape refuses the
// others' family names by name.
bartcore::ResponseFamily resolveFamily(const ParsedControl& control,
                                       const char* familyName) {
  if (control.responseIsBinary) {
    if (std::strcmp(familyName, "logistic") == 0)
      return bartcore::ResponseFamily::logistic;
    if (familyName[0] == '\0' || std::strcmp(familyName, "probit") == 0)
      return bartcore::ResponseFamily::probit;
    Rf_error("unrecognized response family for a binary response");
  }
  // ordinal is the K-level categorical shape (docs/design/ordinal.md): the
  // cumulative probit is the only family defined on it
  if (control.numOrdinalCategories >= 2) {
    if (familyName[0] == '\0' || std::strcmp(familyName, "ordinal") == 0)
      return bartcore::ResponseFamily::ordinal;
    Rf_error("an ordinal (K-level) response supports only family \"ordinal\"");
  }
  // nbinom is the count shape (docs/design/negative-binomial.md section 4): the
  // negative binomial is the only family defined on it
  if (control.countResponse) {
    if (familyName[0] == '\0' || std::strcmp(familyName, "nbinom") == 0)
      return bartcore::ResponseFamily::nbinom;
    Rf_error("a count response supports only family \"nbinom\"");
  }
  // continuous shape: "ordinal" needs the ordered categorical response above
  if (std::strcmp(familyName, "ordinal") == 0)
    Rf_error("family \"ordinal\" requires an ordered categorical response");
  // "nbinom" needs the count response (a non-negative integer y the R surface
  // marks) above
  if (std::strcmp(familyName, "nbinom") == 0)
    Rf_error("family \"nbinom\" requires a non-negative integer (count) response");
  // aft fits a continuous response (log survival times) with truncated latents
  // for right-censored observations; the status arrives on a control attribute
  if (std::strcmp(familyName, "aft") == 0)
    return bartcore::ResponseFamily::aft;
  if (familyName[0] != '\0' && std::strcmp(familyName, "gaussian") != 0)
    Rf_error("response families other than gaussian and aft require a binary "
             "response");
  return bartcore::ResponseFamily::gaussian;
}

// Binary weight policy at sampler creation: a probit has no tractable
// weighted latent-variable form and is refused; logistic treats weights as
// observation counts (its PG(w, psi) latent is the sum of w PG(1, psi)
// draws), so they must be positive integers; gaussian accepts any positive
// weight and is validated elsewhere. The R layer mirrors this, so these
// errors backstop direct-API consumers.
void enforceBinaryWeightPolicy(bartcore::ResponseFamily family,
                               const double* weights,
                               size_t numObservations) {
  if (weights == NULL) return;
  if (family == bartcore::ResponseFamily::probit)
    Rf_error("probit models do not support weights: a weighted probit has no "
             "tractable latent-variable form; use family = \"logistic\" for "
             "weighted binary regression, or model the latents directly");
  if (family == bartcore::ResponseFamily::logistic)
    for (size_t i = 0; i < numObservations; ++i)
      if (!(weights[i] > 0.0) || weights[i] != std::floor(weights[i]))
        Rf_error("logistic weights are observation counts and must be "
                 "positive integers; drop zero-count rows, and use a gaussian "
                 "model for continuous weights");
}

// The post-creation half of the policy: weight changes on binary-response
// samplers are refused outright, since probit never supports weights and
// logistic weights enter the latent construction at creation.
void refuseBinaryWeightChange(const bartcore::SamplerBase& sampler) {
  if (sampler.family() != bartcore::ResponseFamily::gaussian)
    Rf_error("weights on a binary response cannot be set after creation: "
             "probit does not support weights, and logistic weights "
             "(observation counts) are fixed when the sampler is created");
}

// Raw training values of column j, when a dense source serves them: the
// x matrix, or the mixed container's dense slice via the source map.
const double* rawTrainingColumn(const ParsedData& data, size_t j) {
  if (data.xIsMixed)
    return data.columnSources[j] >= 0
      ? data.mixedDenseValues +
        static_cast<size_t>(data.columnSources[j]) * data.numObservations
      : NULL;
  return data.x != NULL ? data.x + j * data.numObservations : NULL;
}

// Raw test values of column j, when a dense source serves them: the x.test
// matrix, or the mixed test container's dense slice. Null for a CSC-backed
// test column (coded R-side against the training levels, nothing to scan).
const double* rawParsedTestColumn(const ParsedData& data, size_t j) {
  if (data.testIsMixed)
    return data.testContainer.columnSources[j] >= 0
      ? data.testContainer.mixedDenseValues +
        static_cast<size_t>(data.testContainer.columnSources[j]) *
          data.numTestObservations
      : NULL;
  return data.x_test != NULL ? data.x_test + j * data.numTestObservations
                             : NULL;
}

void validateCategoricalPredictors(const ParsedData& data) {
  if (data.xIsSparse) return;  // parseData enforced all-ordinal
  for (size_t j = 0; j < data.numPredictors; ++j) {
    if (data.columnTypes[j] != bartcore::ColumnType::categorical) continue;
    // a CSC-backed categorical column has no dense slice to scan; its codes
    // came from the R surface's level table and parseData carried its K and
    // reference code
    if (data.xIsMixed && data.columnSources[j] < 0) continue;
    const double* column = rawTrainingColumn(data, j);
    for (size_t i = 0; i < data.numObservations; ++i) {
      double value = column[i];
      if (bartcore::isNA(value)) continue;  // the reserved missing category
      if (value < 0.0 ||
          value >= static_cast<double>(bartcore::maxCategories) ||
          value != std::floor(value))
        Rf_error("categorical predictors must hold integer category codes "
                 "in [0, 65535)");
    }
  }
  if (data.anyCategorical && data.numTestObservations > 0) {
    // test codes must also be representable; category counts come from the
    // training columns
    for (size_t j = 0; j < data.numPredictors; ++j) {
      if (data.columnTypes[j] != bartcore::ColumnType::categorical) continue;
      // a CSC-backed categorical column has no dense slice to bound the codes
      // against; both its training and test codes came from the R level table
      if (data.xIsMixed && data.columnSources[j] < 0) continue;
      const double* testColumn = rawParsedTestColumn(data, j);
      if (testColumn == NULL) continue;  // CSC-backed test column: R-side coded
      const double* column = rawTrainingColumn(data, j);
      double maxValue = 0.0;
      for (size_t i = 0; i < data.numObservations; ++i) {
        double value = column[i];
        if (!bartcore::isNA(value) && value > maxValue) maxValue = value;
      }
      for (size_t i = 0; i < data.numTestObservations; ++i) {
        double value = testColumn[i];
        if (bartcore::isNA(value)) continue;
        if (value < 0.0 || value > maxValue || value != std::floor(value))
          Rf_error("categorical test predictors must hold existing category "
                   "codes");
      }
    }
  }
}

bartcore::SamplerOptions optionsFromParsed(const ParsedControl& control,
                                           const ParsedModel& model,
                                           const ParsedData& data,
                                           SEXP modelExpr, bool sigmaIsFixed) {
  bartcore::SamplerOptions options;
  options.numTrees = control.numTrees;
  options.sigmaIsFixed = sigmaIsFixed;
  // the gaussian construction path reads this: finite selects Student-t errors
  // (docs/design/robust-errors.md), NaN keeps the Gaussian law
  options.residualDf = model.residualDf;
  // ordinal construction reads this: K >= 2 selects OrdinalResponse with a K-1
  // cutpoint vector (docs/design/ordinal.md), 0 leaves a non-ordinal response
  options.numCategories = control.numOrdinalCategories;
  // nbinom construction reads this: a positive value fixes the dispersion r, a
  // non-positive value estimates it on the grid (docs/design/negative-binomial.md);
  // NaN (a non-count response) is ignored
  options.dispersion = control.dispersion;
  options.k = model.k;
  options.nodeScale = model.nodeScale;
  options.base = model.base;
  options.power = model.power;
  options.birthOrDeathProbability = model.birthOrDeathProbability;
  options.swapProbability = model.swapProbability;
  options.changeProbability = model.changeProbability;
  options.birthProbability = model.birthProbability;
  options.maxNumCutsPerVariable = data.maxNumCuts.data(); // copied at build
  options.useQuantiles = control.useQuantiles;
  options.columnTypes =
    data.anyCategorical ? data.columnTypes.data() : NULL;
  if (data.xIsSparse || data.xIsMixed) {
    options.cscColumnPointers = data.cscColumnPointers;
    options.cscRowIndices = data.cscRowIndices;
    options.cscValues = data.cscValues;
  }
  if (data.xIsMixed) {
    options.mixedDenseValues = data.mixedDenseValues;
    options.columnSources = data.columnSources.data();  // consumed at build
    if (!data.cscCategoryCounts.empty()) {
      options.cscCategoryCounts = data.cscCategoryCounts.data();
      options.cscReferenceCodes = data.cscReferenceCodes.data();
    }
  }
  options.splitProbabilities = model.splitProbabilities; // copied by ctor
  options.monotoneDirections = model.monotoneDirections.empty()
    ? NULL : model.monotoneDirections.data();  // consumed at construction
  options.leafCovariateColumns = model.leafCovariateColumns.empty()
    ? NULL : model.leafCovariateColumns.data();  // consumed at construction
  options.numLeafCovariates = model.leafCovariateColumns.size();
  options.gpLeaves = model.gpLeaves;
  options.gpLengthscales = model.gpLengthscales.empty()
    ? NULL : model.gpLengthscales.data();  // consumed at construction
  options.gpMaxLeafSize = model.gpMaxLeafSize;

  // the generic slot parse above reads the CGM structure a DART prior
  // contains; the Dirichlet configuration comes off the R object directly
  SEXP treePriorExpr =
    PROTECT(Rf_getAttrib(modelExpr, Rf_install("tree.prior")));
  if (Rf_inherits(treePriorExpr, "dbartsDartPrior")) {
    options.useDart = true;
    options.dart.betaA = rc_getDouble(
      Rf_getAttrib(treePriorExpr, Rf_install("a")), "dart a",
      RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
    options.dart.betaB = rc_getDouble(
      Rf_getAttrib(treePriorExpr, Rf_install("b")), "dart b",
      RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
    double rho = REAL(Rf_getAttrib(treePriorExpr, Rf_install("rho")))[0];
    options.dart.rho = ISNAN(rho) ? 0.0 : rho;  // <= 0 means numPredictors
    options.dart.alpha = rc_getDouble(
      Rf_getAttrib(treePriorExpr, Rf_install("alpha")), "dart alpha",
      RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
    options.dart.updateAlpha = rc_getBool(
      Rf_getAttrib(treePriorExpr, Rf_install("update.alpha")), "dart update.alpha",
      RC_LENGTH | RC_EQ, rc_asRLength(1), RC_NA | RC_NO, RC_END);
    options.dart.updateDelay = static_cast<size_t>(rc_getDouble(
      Rf_getAttrib(treePriorExpr, Rf_install("update.delay")), "dart update.delay",
      RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0.0, RC_END));
  }
  UNPROTECT(1);
  options.updateK = model.updateK;
  options.kHyperprior.degreesOfFreedom = model.kDf;
  options.kHyperprior.scale = model.kScale;
  options.numChains = control.numChains;
  options.numThreads = control.numThreads;
  options.numThin = control.treeThinningRate;
  options.keepTrees = control.keepTrees;
  options.numSamplesToStore = control.defaultNumSamples;
  options.verbose = control.verbose;
  options.printEvery = control.printEvery;

  return options;
}

// Every chain gets its own Mersenne twister, so worker threads never touch
// the R API and results do not depend on the thread count. A control rngSeed
// makes results reproducible without R's stream: it seeds a dedicated
// generator that hands each chain its seed, so a single-chain run with seed
// S reproduces chain 0 of any multi-chain run with seed S. Without one,
// chain seeds are drawn from R's stream, so a set.seed() beforehand
// suffices; sampling itself never advances R's stream.
std::vector<ext_rng*> createChainRngs(const ParsedControl& control,
                                      size_t numChains) {
  bool haveSeed = control.haveRngSeed;
  std::vector<ext_rng*> rngs(numChains, static_cast<ext_rng*>(NULL));
  bool rngFailed = false;
  if (haveSeed) {
    ext_rng* seedGenerator =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    rngFailed = seedGenerator == NULL ||
      ext_rng_setSeed(seedGenerator, control.rngSeed) != 0;
    for (size_t c = 0; c < numChains && !rngFailed; ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      rngFailed = rngs[c] == NULL ||
        ext_rng_setSeed(rngs[c], static_cast<std::uint_least32_t>(
          ext_rng_simulateUnsignedIntegerUniformInRange(
            seedGenerator, 0, static_cast<std::uint_least32_t>(-1)))) != 0;
    }
    if (seedGenerator != NULL) ext_rng_destroy(seedGenerator);
  } else {
    GetRNGstate();
    for (size_t c = 0; c < numChains && !rngFailed; ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      rngFailed = rngs[c] == NULL ||
        ext_rng_setSeed(rngs[c], static_cast<std::uint_least32_t>(
                                   unif_rand() * 4294967295.0)) != 0;
    }
    PutRNGstate();
  }
  if (rngFailed) {
    for (size_t c = rngs.size(); c > 0; --c)
      if (rngs[c - 1] != NULL) ext_rng_destroy(rngs[c - 1]);
    Rf_error("could not allocate rng");
  }
  return rngs;
}

// rbart_vi's in-core Gibbs path passes its grouping through an internal
// attribute on the control object: 1-based per-observation group indices,
// the group count, the built-in tau prior's name and original-scale
// relative scale, and the slice-step count. Public surfaces never set the
// attribute, and only full creation reads it (setControl ignores it; the
// group structure is fixed at creation). groupIndices outlives the options
// borrow: the chains copy it during construction.
void applyGroupAttribute(SEXP controlExpr, size_t numObservations,
                         bartcore::SamplerOptions& options,
                         std::vector<std::uint32_t>& groupIndices) {
  SEXP groupsExpr = Rf_getAttrib(controlExpr, Rf_install("bartcore.groups"));
  if (Rf_isNull(groupsExpr)) return;

  SEXP indicesExpr = getListElement(groupsExpr, "indices");
  SEXP numGroupsExpr = getListElement(groupsExpr, "n.groups");
  SEXP priorExpr = getListElement(groupsExpr, "prior");
  SEXP scaleExpr = getListElement(groupsExpr, "rel.scale");
  SEXP stepsExpr = getListElement(groupsExpr, "n.steps");
  if (!Rf_isInteger(indicesExpr) ||
      static_cast<size_t>(Rf_xlength(indicesExpr)) != numObservations ||
      !Rf_isInteger(numGroupsExpr) || Rf_xlength(numGroupsExpr) != 1 ||
      !Rf_isString(priorExpr) || Rf_xlength(priorExpr) != 1 ||
      !Rf_isReal(scaleExpr) || Rf_xlength(scaleExpr) != 1 ||
      !Rf_isInteger(stepsExpr) || Rf_xlength(stepsExpr) != 1)
    Rf_error("malformed grouped random effects specification");

  int numGroups = INTEGER(numGroupsExpr)[0];
  if (numGroups < 1)
    Rf_error("grouped random effects require at least one group");
  groupIndices.resize(numObservations);
  for (size_t i = 0; i < numObservations; ++i) {
    int index = INTEGER(indicesExpr)[i];
    if (index < 1 || index > numGroups)
      Rf_error("group indices must be in [1, number of groups]");
    groupIndices[i] = static_cast<std::uint32_t>(index - 1);
  }

  const char* priorName = CHAR(STRING_ELT(priorExpr, 0));
  if (std::strcmp(priorName, "cauchy") == 0) {
    options.tauPriorKind = bartcore::TauPriorKind::cauchy;
  } else if (std::strcmp(priorName, "gamma") == 0) {
    options.tauPriorKind = bartcore::TauPriorKind::gamma;
  } else {
    Rf_error("unrecognized tau prior for grouped random effects");
  }

  double relScale = REAL(scaleExpr)[0];
  if (!(relScale > 0.0)) Rf_error("tau prior scale must be positive");
  int numSteps = INTEGER(stepsExpr)[0];
  if (numSteps < 1) Rf_error("tau slice steps must be at least 1");

  options.groupIndices = groupIndices.data();
  options.numGroups = static_cast<size_t>(numGroups);
  options.tauPriorScale = relScale;
  options.tauSliceSteps = static_cast<size_t>(numSteps);
}

// AFT survival status arrives on an internal control attribute alongside the
// log-time response (the y creation argument): a per-observation numeric
// vector, 1 for an uncensored event, 0 for a right-censored observation. Only
// full creation reads it, and only for family aft. status outlives the options
// borrow: the response copies it during construction. The public surface sets
// the attribute from a Surv or (time, status) ingest (docs/design/survival.md).
void applySurvivalAttribute(SEXP controlExpr, size_t numObservations,
                            bartcore::ResponseFamily family,
                            bartcore::SamplerOptions& options,
                            std::vector<double>& status) {
  SEXP statusExpr = Rf_getAttrib(controlExpr, Rf_install("bartcore.survival"));
  if (family != bartcore::ResponseFamily::aft) {
    if (!Rf_isNull(statusExpr))
      Rf_error("survival status supplied for a non-aft response family");
    return;
  }
  if (Rf_isNull(statusExpr))
    Rf_error("aft (survival) models require a status vector");
  if (!Rf_isReal(statusExpr) ||
      static_cast<size_t>(Rf_xlength(statusExpr)) != numObservations)
    Rf_error("survival status must be a numeric vector of length equal to the "
             "number of observations");
  const double* raw = REAL(statusExpr);
  status.resize(numObservations);
  for (size_t i = 0; i < numObservations; ++i) {
    if (raw[i] != 0.0 && raw[i] != 1.0)
      Rf_error("survival status must be 0 (censored) or 1 (event)");
    status[i] = raw[i];
  }
  options.survivalStatus = status.data();
}

// A sampler created over a data handle holds no raw predictor values, so
// the raw-x mutation surface has nothing to work from; a CSC-built sampler
// holds only borrowed slices and refuses the same surface by design
// (docs/design/sparse-columns.md).
void refusePredictorMutation(const bartcore::SamplerBase& sampler,
                             const char* caller) {
  const bartcore::ColumnStore& data = sampler.data();
  if (data.acceptsNewRawPredictors()) return;
  if (data.builtFromCsc)
    Rf_error("%s: sparse predictors fix the design at creation; make a new "
             "sampler instead", caller);
  Rf_error("%s requires a sampler that owns its predictors; data-handle "
           "views hold none", caller);
}

// Unlike views, CSC-built samplers re-quantize from their retained slices,
// so cut installation and state restore stay available on them.
void refuseRequantizeWithoutSource(const bartcore::SamplerBase& sampler,
                                   const char* caller) {
  if (!sampler.data().hasRequantizeSource())
    Rf_error("%s requires a sampler that owns its predictors; data-handle "
             "views hold none", caller);
}

// The single-forest test-fit and prediction surface has no meaning under BCF:
// with two forests and no test treatment vector, a blend a * mu + b_z * tau is
// ill-defined, so the engine would fall back to the bare prognostic forest and
// silently misreport. Reject test data and out-of-sample prediction on a BCF
// sampler; consumers recombine per forest via getForestFits + the BCF glue
// (docs/design/bcf.md). Gated on testFitsAreDefined() rather than the forest
// count so a multi-forest model whose test blend IS defined (multinomial
// softmax over the K forests' totalTestFits) is allowed through; only BCF, the
// one multi-forest model that leaves the channel undefined, is refused.
void refuseBCFTestSurface(const bartcore::SamplerBase& sampler,
                          const char* caller) {
  if (sampler.numForests() >= 2 && !sampler.testFitsAreDefined())
    Rf_error("%s: a BCF sampler carries no test treatment vector, so its test "
             "fits are undefined; predict per forest with getForestFits and "
             "the BCF glue instead", caller);
}

// A whole-data mutation (setData, setResponse, setWeights) rebuilds only
// forest 0: applyNewData and the response/latent refresh touch forests_[0], so
// on a multi-forest sampler (BCF, and any future multi-forest model) the other
// forests would keep fits against the old data. Refuse it; a multi-forest
// sampler fixes its data at creation, as grouped/sparse/aft samplers do.
// setTreatment, the one supported multi-forest data swap, routes through the
// combiner and stays allowed.
void refuseMultiForestMutation(const bartcore::SamplerBase& sampler,
                               const char* caller) {
  if (sampler.numForests() >= 2)
    Rf_error("%s: a multi-forest sampler fixes its data at creation; make a "
             "new sampler instead", caller);
}

// The transactional predictor paths (setPredictor and updatePredictor without
// forceUpdate, and the per-observation sessions, which have no force variant)
// validate and rebuild through revalidateAllChains, which revalidates only the
// primary forest - an accepted change would leave a multi-forest sampler's
// other forests routed against stale codes. The FORCE paths refresh every
// forest (forceRefreshTrees) and stay available: a forced whole-matrix
// setPredictor is the supported multi-forest predictor swap (the bartCause
// propensity pattern).
void refuseMultiForestTransactionalUpdate(const bartcore::SamplerBase& sampler,
                                          const char* caller) {
  if (sampler.numForests() >= 2)
    Rf_error("%s: a transactional predictor update validates only the primary "
             "forest of a multi-forest sampler; use setPredictor with "
             "forceUpdate = TRUE or make a new sampler instead", caller);
}

// A built column store (cuts + codes) shared by row-subset view samplers
// (public-surface.md section 5; internal). The external pointer's
// protection slot pins the data expression whose x the store borrows.
struct DataHandle {
  bartcore::ColumnStore store;
  // a mixed store borrows dense slices of a transiently assembled block; the
  // handle owns it so views can gather from the parent after creation returns
  std::vector<double> ownedMixedDense;
};

void dataHandleFinalizer(SEXP ptrExpr) {
  DataHandle* handle = static_cast<DataHandle*>(R_ExternalPtrAddr(ptrExpr));
  if (handle == NULL) return;
  delete handle;
  R_ClearExternalPtr(ptrExpr);
}

DataHandle& dataHandleFromExpression(SEXP ptrExpr) {
  DataHandle* handle = static_cast<DataHandle*>(R_ExternalPtrAddr(ptrExpr));
  if (handle == NULL)
    Rf_error("data handle function called on NULL external pointer");
  return *handle;
}

// Roots a freshly allocated result column in its protected container and
// hands it back, so run-result assembly needs no per-column PROTECT.
SEXP installResult(SEXP resultExpr, int slot, SEXP value) {
  SET_VECTOR_ELT(resultExpr, slot, value);
  return value;
}

// Shared parse/validate prologue of the two sampler-creation paths: fills
// the parsed views from the specification objects, resolves the response
// family, applies the fixed-sigma variance semantics, and enforces the
// binary weight policy at creation.
bartcore::ResponseFamily parseSamplerSpecification(
    SEXP controlExpr, SEXP modelExpr, SEXP dataExpr, const char* familyName,
    ParsedControl& control, ParsedModel& model, ParsedData& data,
    bool& sigmaIsFixed) {
  parseControl(control, controlExpr);
  parseData(data, dataExpr);
  parseModel(model, modelExpr, data.numPredictors);

  bartcore::ResponseFamily family = resolveFamily(control, familyName);
  sigmaIsFixed = !control.responseIsBinary && model.sigmaIsFixed;
  // documented semantics: fixed(value) holds the residual variance at
  // value, so sigma enters as sqrt(value) and is never drawn
  if (sigmaIsFixed) data.sigmaEstimate = std::sqrt(model.fixedSigmaSq);
  enforceBinaryWeightPolicy(family, data.weights, data.numObservations);
  // a weighted truncated-latent draw is not a coherent likelihood; AFT v1
  // rejects weights (docs/design/survival.md)
  if (family == bartcore::ResponseFamily::aft && data.weights != NULL)
    Rf_error("aft (survival) models do not support case weights");
  // ordinal refuses weights for probit's reason: a weighted truncated-normal
  // latent likelihood is not a coherent model (docs/design/ordinal.md)
  if (family == bartcore::ResponseFamily::ordinal && data.weights != NULL)
    Rf_error("ordinal models do not support weights: a weighted truncated-"
             "normal latent likelihood is not a coherent model");
  // nbinom refuses weights in v1 (docs/design/negative-binomial.md section 4):
  // the usual count "weight" is exposure, which belongs in the offset as a
  // log-exposure term, not in observation replication
  if (family == bartcore::ResponseFamily::nbinom && data.weights != NULL)
    Rf_error("nbinom (count) models do not support weights: exposure belongs in "
             "the offset as a log-exposure term");
  // Student-t continuous errors (docs/design/robust-errors.md): a finite
  // resid.df selects the scale-mixture error law - a positive value fixes the
  // degrees of freedom, 0 estimates them on the grid. Only the gaussian family
  // carries it; NaN (absent) keeps the Gaussian law.
  if (std::isfinite(model.residualDf)) {
    if (family != bartcore::ResponseFamily::gaussian)
      Rf_error("Student-t residuals (resid.df) require a continuous gaussian "
               "response");
    if (model.residualDf < 0.0)
      Rf_error("resid.df must be positive to fix the degrees of freedom, or 0 "
               "to estimate them");
  }
  return family;
}

} // namespace

namespace bartcore_bridge {

void validateColumnValues(const bartcore::ColumnStore& store, size_t column,
                          const double* values, size_t numValues) {
  if (store.types[column] != bartcore::ColumnType::categorical) return;
  for (size_t i = 0; i < numValues; ++i) {
    if (!store.categoricalValueIsValid(column, values[i]))
      Rf_error("categorical predictor values must be existing category codes");
  }
}

// family resolution: see resolveFamily() above.
BartcoreHolder* createHolder(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                             const char* familyName) {
  BartcoreHolder* holder = nullptr;
  unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                 model = ParsedModel{},
                 groupIndices = std::vector<std::uint32_t>{},
                 survivalStatus = std::vector<double>{},
                 rngs = std::vector<ext_rng*>{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    bartcore::ResponseFamily family = parseSamplerSpecification(
      controlExpr, modelExpr, dataExpr, familyName, control, model, data,
      sigmaIsFixed);
    validateCategoricalPredictors(data);

    bartcore::SamplerOptions options =
      optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);

    // grouped random intercepts (rbart_vi's in-core path) arrive on an
    // internal control attribute; the chains copy the indices at construction
    applyGroupAttribute(controlExpr, data.numObservations, options,
                        groupIndices);
    // grouped ordinal is a recorded but unbuilt door (docs/design/ordinal.md
    // section 8): the cutpoint block and the group block are not yet shown to
    // interleave, so refuse the composition here, the host backstop the R
    // surface (rbart_vi, C3) mirrors
    if (family == bartcore::ResponseFamily::ordinal && options.numGroups > 0)
      Rf_error("grouped random effects are not supported for ordinal responses");
    // grouped nbinom is a recorded but unbuilt door (docs/design/negative-binomial.md
    // section 7): the dispersion block and the group block are not yet shown to
    // interleave, so refuse the composition here, the backstop rbart_vi (C3) mirrors
    if (family == bartcore::ResponseFamily::nbinom && options.numGroups > 0)
      Rf_error("grouped random effects are not supported for count (nbinom) responses");
    // AFT survival status arrives the same way; the response copies it
    applySurvivalAttribute(controlExpr, data.numObservations, family, options,
                           survivalStatus);

    rngs = createChainRngs(control, options.numChains);

    // dispatches on the leaf model: a linear node prior's designated columns
    // select the linear-leaf instantiation, everything else the constant leaf
    std::unique_ptr<bartcore::SamplerBase> sampler = bartcore::createSampler(
      data.x, data.y, data.numObservations, data.numPredictors, data.weights,
      data.offset, family, data.sigmaEstimate, model.sigmaDf,
      model.sigmaRawScale, options, rngs.data());
    if (sampler == NULL) {
      // R-side resolution validates first, so only an invariant breach lands
      for (ext_rng* rng : rngs) if (rng != NULL) ext_rng_destroy(rng);
      Rf_error("invalid leaf covariate designation");
    }

    if (data.numTestObservations > 0) {
      if (data.testIsMixed) {
        // the container's raw is copied into the test store; a leaf covariate
        // that would land on a CSC-backed test column is refused (sparse
        // storage serves no dense raw test covariate)
        if (!installTestContainer(*sampler, data.testContainer)) {
          sampler.reset();  // borrows the rngs, so tear it down before them
          for (ext_rng* rng : rngs) if (rng != NULL) ext_rng_destroy(rng);
          Rf_error("a leaf covariate column cannot be a sparse test column; "
                   "supply it as a dense test column");
        }
      } else {
        sampler->setTestPredictors(data.x_test, data.numTestObservations);
      }
      sampler->setTestOffset(data.testOffset);
    }

    if (control.verbose) printInitialSummary(control, model, data, *sampler);

    holder = new BartcoreHolder{std::move(sampler), std::move(rngs),
                                control.keepTrainingFits, {}, {}, {}, {}, {}, {}, {}};
    // the mixed store borrows dense slices of the transiently assembled block;
    // the holder owns it so they outlive this call (a vector move preserves
    // the buffer address the store cached)
    if (data.xIsMixed)
      holder->ownedMixedDense = std::move(data.denseAssembly);
    return R_NilValue;
  });
  return holder;
}

BartcoreHolder* createBCFHolder(SEXP controlExpr, SEXP modelExpr,
                                SEXP dataExpr, SEXP zExpr,
                                SEXP bcfParamsExpr, SEXP moderatorsExpr) {
  if (!Rf_isReal(bcfParamsExpr) || Rf_xlength(bcfParamsExpr) != 8)
    Rf_error("bcf parameters must be a length-8 numeric vector");

  BartcoreHolder* holder = nullptr;
  unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                 model = ParsedModel{}, rngs = std::vector<ext_rng*>{},
                 z = std::vector<double>{},
                 moderators = std::vector<size_t>{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    bartcore::ResponseFamily family = parseSamplerSpecification(
      controlExpr, modelExpr, dataExpr, "", control, model, data, sigmaIsFixed);
    validateCategoricalPredictors(data);
    if (family != bartcore::ResponseFamily::gaussian)
      Rf_error("BCF requires a continuous (gaussian) response");
    if (data.x == NULL)
      Rf_error("BCF requires dense predictors");
    if (static_cast<size_t>(Rf_xlength(zExpr)) != data.numObservations)
      Rf_error("treatment length must match the number of observations");

    z.resize(data.numObservations);
    for (size_t i = 0; i < data.numObservations; ++i)
      z[i] = REAL(zExpr)[i] != 0.0 ? 1.0 : 0.0;

    bartcore::SamplerOptions options =
      optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);
    rngs = createChainRngs(control, options.numChains);

    const double* p = REAL(bcfParamsExpr);
    bartcore::BCFSpec spec;
    // node scales come from the calibration map, not the host model
    spec.mu.numTrees = options.numTrees;
    spec.mu.base = model.base;
    spec.mu.power = model.power;
    spec.tau.numTrees = static_cast<size_t>(p[0]);
    spec.tau.base = p[1];
    spec.tau.power = p[2];
    spec.aPriorScale = p[3];
    spec.sdModerate = p[4];
    spec.bPriorVariance = p[5];
    spec.updateA = p[6] != 0.0;
    spec.updateB = p[7] != 0.0;
    spec.z = z.data();

    // the treatment forest's optional moderator restriction: 1-based column
    // indices resolved R-side, or NULL for no restriction (tau reads the full
    // store). Consumed at construction, so the buffer need only outlive the
    // sampler build.
    if (!Rf_isNull(moderatorsExpr)) {
      if (!Rf_isInteger(moderatorsExpr))
        Rf_error("bcf moderators must be resolved integer column indices");
      R_xlen_t numModerators = Rf_xlength(moderatorsExpr);
      moderators.resize(static_cast<size_t>(numModerators));
      for (R_xlen_t j = 0; j < numModerators; ++j) {
        int column = INTEGER(moderatorsExpr)[j];
        if (column < 1 || static_cast<size_t>(column) > data.numPredictors)
          Rf_error("bcf moderator column out of range");
        moderators[static_cast<size_t>(j)] = static_cast<size_t>(column - 1);
      }
      spec.tau.columns = moderators.data();
      spec.tau.numColumns = moderators.size();
    }

    std::unique_ptr<bartcore::SamplerBase> sampler = bartcore::createBCFSampler(
      data.x, data.y, data.numObservations, data.numPredictors, data.weights,
      data.offset, data.sigmaEstimate, model.sigmaDf, model.sigmaRawScale,
      options, spec, rngs.data());

    holder = new BartcoreHolder{std::move(sampler), std::move(rngs),
                                control.keepTrainingFits, {}, {}, {}, {}, {}, {}, {}};
    // moving z keeps its buffer, so the chains' borrowed z stays valid
    holder->ownedTreatment = std::move(z);
    return R_NilValue;
  });
  return holder;
}

// The parse and validation both multinomial entries share: parse the sampler
// spec, refuse the response combinations the single-trial softmax cannot carry
// (case weights, a test offset, a mixed test store), and require dense
// predictors and K >= 2. Runs before any rng is created, so a refusal needs no
// cleanup. The predictors ride the data object; the data's response is ignored
// (the counts are the response, borrowed separately).
static void parseMultinomialData(SEXP controlExpr, SEXP modelExpr,
                                 SEXP dataExpr, ParsedControl& control,
                                 ParsedModel& model, ParsedData& data,
                                 bool& sigmaIsFixed, size_t numCategories) {
  parseSamplerSpecification(controlExpr, modelExpr, dataExpr, "", control, model,
                            data, sigmaIsFixed);
  validateCategoricalPredictors(data);
  if (data.x == NULL)
    Rf_error("multinomial requires dense predictors");
  if (data.weights != NULL)
    Rf_error("multinomial (softmax) models do not support case weights");
  if (data.testOffset != NULL)
    Rf_error("multinomial (softmax) models do not support a test offset");
  if (data.numTestObservations > 0 && data.testIsMixed)
    Rf_error("multinomial (softmax) models require a dense test matrix");
  if (numCategories < 2)
    Rf_error("multinomial requires at least two categories");
}

// Builds the K-forest multinomial sampler shared by both entries: sizes the
// options, creates the per-chain rngs, builds the sampler from the count-native
// spec, and sets the test predictors. counts is the borrowed category-major
// n x K matrix and trials the per-observation trial counts; both the label
// entry (one-hot, unit trials) and the count entry route through here so their
// draw streams are the one code path. Leaf-scale and k follow the multinomial
// calibration, not the host node prior (a gaussian default).
static std::unique_ptr<bartcore::SamplerBase> buildMultinomialSampler(
    const ParsedControl& control, const ParsedModel& model,
    const ParsedData& data, SEXP modelExpr, bool sigmaIsFixed,
    size_t numCategories, const int* counts, const int* trials,
    std::vector<ext_rng*>& rngs) {
  bartcore::SamplerOptions options =
    optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);
  rngs = createChainRngs(control, options.numChains);

  bartcore::MultinomialSpec spec;
  spec.numCategories = numCategories;
  spec.counts = counts;
  spec.trials = trials;
  // the K=2 pairwise-log-odds anchor (pi*sqrt(3)/sqrt(2)) and the host k; the
  // host node scale (a gaussian default) is deliberately not read here
  spec.k = model.k;
  spec.forest.numTrees = options.numTrees;
  spec.forest.base = model.base;
  spec.forest.power = model.power;
  spec.forest.birthOrDeathProbability = model.birthOrDeathProbability;
  spec.forest.swapProbability = model.swapProbability;
  spec.forest.changeProbability = model.changeProbability;
  spec.forest.birthProbability = model.birthProbability;

  std::unique_ptr<bartcore::SamplerBase> sampler =
    bartcore::createMultinomialSampler(data.x, data.numObservations,
                                       data.numPredictors, options, spec,
                                       rngs.data());

  // test-at-creation: the K forests each accumulate their own totalTestFits in
  // the sweep, and storeSample blends them into the K softmax test
  // probabilities (chain testFitsAreDefined() is true). buildTest copies the
  // dense test values, so nothing is pinned; parseMultinomialData ruled out a
  // test offset and a mixed test store.
  if (data.numTestObservations > 0)
    sampler->setTestPredictors(data.x_test, data.numTestObservations);
  return sampler;
}

// A K-forest multinomial (softmax) sampler; internal, constant-leaf only
// (docs/design/multinomial.md). The single-trial label entry: the category codes
// (0..K-1, one per observation) become a one-hot n x K count matrix with every
// trial 1, the exact n_i = 1 reduction of the count-native combiner.
BartcoreHolder* createMultinomialHolder(SEXP controlExpr, SEXP modelExpr,
                                        SEXP dataExpr, SEXP labelsExpr,
                                        SEXP numCategoriesExpr) {
  BartcoreHolder* holder = nullptr;
  unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                 model = ParsedModel{}, rngs = std::vector<ext_rng*>{},
                 counts = std::vector<int>{},
                 trials = std::vector<int>{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    size_t numCategories = static_cast<size_t>(Rf_asInteger(numCategoriesExpr));
    parseMultinomialData(controlExpr, modelExpr, dataExpr, control, model, data,
                         sigmaIsFixed, numCategories);

    if (!Rf_isInteger(labelsExpr) ||
        static_cast<size_t>(Rf_xlength(labelsExpr)) != data.numObservations)
      Rf_error("multinomial labels must be an integer code per observation");
    size_t n = data.numObservations;
    // one-hot category-major counts (column k contiguous at k*n) with unit
    // trials: at n_i = 1 the combiner's PG summing loop and (y - n_i/2) working
    // response reduce byte-identically to the label path
    counts.assign(n * numCategories, 0);
    trials.assign(n, 1);
    const int* src = INTEGER(labelsExpr);
    for (size_t i = 0; i < n; ++i) {
      if (src[i] < 0 || static_cast<size_t>(src[i]) >= numCategories)
        Rf_error("multinomial label out of range 0..K-1");
      counts[static_cast<size_t>(src[i]) * n + i] = 1;
    }

    std::unique_ptr<bartcore::SamplerBase> sampler = buildMultinomialSampler(
      control, model, data, modelExpr, sigmaIsFixed, numCategories,
      counts.data(), trials.data(), rngs);

    holder = new BartcoreHolder{std::move(sampler), std::move(rngs),
                                control.keepTrainingFits, {}, {}, {}, {}, {}, {}, {}};
    // moving keeps the buffers, so the combiner's borrowed counts/trials stay valid
    holder->ownedCounts = std::move(counts);
    holder->ownedTrials = std::move(trials);
    return R_NilValue;
  });
  return holder;
}

// A K-forest multinomial (softmax) sampler over a GROUPED-COUNT response: Y is
// an n x K nonnegative integer matrix, category-major (R column-major = the
// combiner's counts_ layout, so the buffer copies directly), and the trials
// n_i = sum_k Y_ik must be >= 1 (an empty row carries no information; a PG(0, .)
// point mass at 0 would break the working response). Same engine as the label
// entry, count-native (docs/design/multinomial.md).
BartcoreHolder* createMultinomialCountsHolder(SEXP controlExpr, SEXP modelExpr,
                                              SEXP dataExpr, SEXP countsExpr,
                                              SEXP numCategoriesExpr) {
  BartcoreHolder* holder = nullptr;
  unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                 model = ParsedModel{}, rngs = std::vector<ext_rng*>{},
                 counts = std::vector<int>{},
                 trials = std::vector<int>{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    size_t numCategories = static_cast<size_t>(Rf_asInteger(numCategoriesExpr));
    parseMultinomialData(controlExpr, modelExpr, dataExpr, control, model, data,
                         sigmaIsFixed, numCategories);

    size_t n = data.numObservations;
    if (!Rf_isInteger(countsExpr) ||
        static_cast<size_t>(Rf_xlength(countsExpr)) != n * numCategories)
      Rf_error("multinomial counts must be an n x K integer matrix");
    const int* src = INTEGER(countsExpr);
    counts.assign(src, src + n * numCategories);
    trials.assign(n, 0);
    for (size_t k = 0; k < numCategories; ++k)
      for (size_t i = 0; i < n; ++i) {
        int y = counts[k * n + i];
        if (y < 0) Rf_error("multinomial counts must be nonnegative");
        trials[i] += y;
      }
    for (size_t i = 0; i < n; ++i)
      if (trials[i] < 1)
        Rf_error("every multinomial count row must have at least one trial");

    std::unique_ptr<bartcore::SamplerBase> sampler = buildMultinomialSampler(
      control, model, data, modelExpr, sigmaIsFixed, numCategories,
      counts.data(), trials.data(), rngs);

    holder = new BartcoreHolder{std::move(sampler), std::move(rngs),
                                control.keepTrainingFits, {}, {}, {}, {}, {}, {}, {}};
    // moving keeps the buffers, so the combiner's borrowed counts/trials stay valid
    holder->ownedCounts = std::move(counts);
    holder->ownedTrials = std::move(trials);
    return R_NilValue;
  });
  return holder;
}

} // namespace bartcore_bridge

extern "C" {

// The external pointer's protection slot pins everything the sampler
// borrows: the data expression at creation, and any replacement vectors the
// setters install later. family is as createHolder's.
SEXP bartcore_create(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                     SEXP familyExpr) {
  const char* familyName =
    Rf_isNull(familyExpr) ? "" : CHAR(STRING_ELT(familyExpr, 0));
  // Register the finalizer on a null-address pointer, then build and install:
  // R_SetExternalPtrAddr does not allocate, so an OOM longjmp cannot strand the
  // holder between its construction and its finalizer.
  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  BartcoreHolder* holder = bartcore_bridge::createHolder(
    controlExpr, modelExpr, dataExpr, familyName);
  R_SetExternalPtrAddr(result, holder);

  UNPROTECT(2);
  return result;
}

// Builds the two-layer store (cuts + codes) once for sharing across
// row-subset samplers: control contributes useQuantiles, data contributes
// x, the column types, and n.cuts. Internal, with no serialization;
// see public-surface.md section 5.
SEXP bartcore_createDataHandle(SEXP controlExpr, SEXP dataExpr,
                               SEXP leafCovariateColumnsExpr) {
  return unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                        gatherColumns =
                          std::vector<size_t>{}]() mutable -> SEXP {
    parseControl(control, controlExpr);
    parseData(data, dataExpr);
    validateCategoricalPredictors(data);

    // 1-based indices of the columns whose raw values a view's leaf model will
    // read, or NULL for none. The handle serves views whose leaf specs are
    // unknown at creation, so the caller declares the union its views need; an
    // undeclared column is left ungathered and its view designation refused
    // downstream (rawColumn null). Dense builds only: a mixed build's
    // dense-backed columns already serve raw, and a sparse build serves none.
    if (!Rf_isNull(leafCovariateColumnsExpr)) {
      if (!Rf_isInteger(leafCovariateColumnsExpr))
        Rf_error("leaf covariate columns must be an integer vector or NULL");
      R_xlen_t numGather = Rf_xlength(leafCovariateColumnsExpr);
      gatherColumns.resize(static_cast<size_t>(numGather));
      for (R_xlen_t j = 0; j < numGather; ++j) {
        int column = INTEGER(leafCovariateColumnsExpr)[j];
        if (column < 1 || static_cast<size_t>(column) > data.numPredictors)
          Rf_error("leaf covariate column out of range");
        gatherColumns[static_cast<size_t>(j)] = static_cast<size_t>(column - 1);
      }
    }

    DataHandle* handle = new DataHandle;
    if (data.xIsMixed) {
      handle->store.buildMixed(data.mixedDenseValues, data.cscColumnPointers,
                               data.cscRowIndices, data.cscValues,
                               data.columnSources.data(), data.numObservations,
                               data.numPredictors, data.maxNumCuts.data(), 0,
                               control.useQuantiles,
                               data.anyCategorical ? data.columnTypes.data()
                                                   : NULL,
                               data.cscCategoryCounts.empty()
                                 ? NULL : data.cscCategoryCounts.data(),
                               data.cscReferenceCodes.empty()
                                 ? NULL : data.cscReferenceCodes.data());
      // the store borrows dense slices of the transiently assembled block;
      // the handle owns it so views can gather from the parent afterwards
      handle->ownedMixedDense = std::move(data.denseAssembly);
    } else if (data.xIsSparse) {
      handle->store.buildFromCsc(data.cscColumnPointers, data.cscRowIndices,
                                 data.cscValues, data.numObservations,
                                 data.numPredictors, data.maxNumCuts.data(), 0,
                                 control.useQuantiles);
    } else {
      // the handle owns raw only for the declared leaf-covariate columns; a
      // view then gathers any of those it designates from parent.rawColumn
      handle->store.build(data.x, data.numObservations, data.numPredictors,
                          data.maxNumCuts.data(), control.useQuantiles,
                          data.anyCategorical ? data.columnTypes.data() : NULL,
                          gatherColumns.empty() ? NULL : gatherColumns.data(),
                          gatherColumns.size());
    }

    SEXP result = PROTECT(R_MakeExternalPtr(handle, R_NilValue, dataExpr));
    R_RegisterCFinalizerEx(result, dataHandleFinalizer,
                           static_cast<Rboolean>(FALSE));
    UNPROTECT(1);
    return result;
  });
}

// A sampler over a row-subset view of a data handle: the view copies the
// handle's cut grid and gathers the subset's codes, so folds bin
// identically to the full data. dataExpr is the full data object the
// handle was built from; the response vectors are sliced here by the
// 1-based trainRows and owned by the holder, and a test offset is gathered
// from the regular offset at testRows (xbart's fold semantics). The result
// refuses the raw-x mutation surface (setPredictor and friends, setData,
// setCutPoints, setState).
SEXP bartcore_createFromHandle(SEXP controlExpr, SEXP modelExpr,
                               SEXP dataExpr, SEXP handleExpr,
                               SEXP trainRowsExpr, SEXP testRowsExpr,
                               SEXP familyExpr, SEXP columnsExpr) {
  const bartcore::ColumnStore& parent =
    dataHandleFromExpression(handleExpr).store;
  const char* familyName =
    Rf_isNull(familyExpr) ? "" : CHAR(STRING_ELT(familyExpr, 0));

  return unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                        model = ParsedModel{},
                        trainRows = std::vector<size_t>{},
                        testRows = std::vector<size_t>{},
                        response = std::vector<double>{},
                        weights = std::vector<double>{},
                        offset = std::vector<double>{},
                        testOffset = std::vector<double>{},
                        store = bartcore::ColumnStore{},
                        columns = std::vector<size_t>{},
                        viewLeafCovariates = std::vector<size_t>{},
                        rngs = std::vector<ext_rng*>{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    bartcore::ResponseFamily family = parseSamplerSpecification(
      controlExpr, modelExpr, dataExpr, familyName, control, model, data,
      sigmaIsFixed);

    if (data.numObservations != parent.numObservations ||
        data.numPredictors != parent.numPredictors)
      Rf_error("data does not match the shape the handle was built from");

    if (!Rf_isInteger(trainRowsExpr) || Rf_xlength(trainRowsExpr) == 0)
      Rf_error("trainRows must be a non-empty integer vector");
    if (!Rf_isNull(testRowsExpr) && !Rf_isInteger(testRowsExpr))
      Rf_error("testRows must be an integer vector or NULL");
    size_t numTrainRows = static_cast<size_t>(Rf_xlength(trainRowsExpr));
    size_t numTestRows = Rf_isNull(testRowsExpr)
      ? 0 : static_cast<size_t>(Rf_xlength(testRowsExpr));

    trainRows.resize(numTrainRows);
    testRows.resize(numTestRows);
    const int* i_trainRows = INTEGER(trainRowsExpr);
    for (size_t i = 0; i < numTrainRows; ++i) {
      if (i_trainRows[i] < 1 ||
          static_cast<size_t>(i_trainRows[i]) > parent.numObservations)
        Rf_error("train row out of range");
      trainRows[i] = static_cast<size_t>(i_trainRows[i] - 1);
    }
    for (size_t i = 0; i < numTestRows; ++i) {
      int row = INTEGER(testRowsExpr)[i];
      if (row < 1 || static_cast<size_t>(row) > parent.numObservations)
        Rf_error("test row out of range");
      testRows[i] = static_cast<size_t>(row - 1);
    }

    response.resize(numTrainRows);
    for (size_t i = 0; i < numTrainRows; ++i)
      response[i] = data.y[trainRows[i]];
    if (data.weights != NULL) {
      weights.resize(numTrainRows);
      for (size_t i = 0; i < numTrainRows; ++i)
        weights[i] = data.weights[trainRows[i]];
    }
    if (data.offset != NULL) {
      offset.resize(numTrainRows);
      for (size_t i = 0; i < numTrainRows; ++i)
        offset[i] = data.offset[trainRows[i]];
      if (numTestRows > 0) {
        testOffset.resize(numTestRows);
        for (size_t i = 0; i < numTestRows; ++i)
          testOffset[i] = data.offset[testRows[i]];
      }
    }

    // an optional column subset: 1-based indices into the handle's columns,
    // or NULL for the full set. The view spans only these, binning each
    // identically to the parent. Consumed at construction, so the buffer need
    // only outlive the sampler build.
    if (!Rf_isNull(columnsExpr)) {
      if (!Rf_isInteger(columnsExpr) || Rf_xlength(columnsExpr) == 0)
        Rf_error("view columns must be a non-empty integer vector or NULL");
      R_xlen_t numColumns = Rf_xlength(columnsExpr);
      columns.resize(static_cast<size_t>(numColumns));
      for (R_xlen_t j = 0; j < numColumns; ++j) {
        int column = INTEGER(columnsExpr)[j];
        if (column < 1 || static_cast<size_t>(column) > parent.numPredictors)
          Rf_error("view column out of range");
        columns[static_cast<size_t>(j)] = static_cast<size_t>(column - 1);
      }
    }

    // leaf covariates address the view's own columns; translate the
    // parent-space designation onto the subset (identity for a full-span view)
    const size_t* gatherColumns = model.leafCovariateColumns.empty()
      ? NULL : model.leafCovariateColumns.data();
    size_t numGatherColumns = model.leafCovariateColumns.size();
    if (!columns.empty() && numGatherColumns > 0) {
      viewLeafCovariates.reserve(numGatherColumns);
      for (size_t parentColumn : model.leafCovariateColumns) {
        size_t viewColumn = columns.size();
        for (size_t j = 0; j < columns.size(); ++j)
          if (columns[j] == parentColumn) { viewColumn = j; break; }
        if (viewColumn == columns.size())
          Rf_error("leaf covariate column absent from the view's columns");
        viewLeafCovariates.push_back(viewColumn);
      }
      gatherColumns = viewLeafCovariates.data();
      numGatherColumns = viewLeafCovariates.size();
    }

    // a linear node prior's designated columns have the view gather their raw
    // values, with standardization constants from the handle's full data - the
    // same calibration inheritance as the copied cut grid
    store.buildFromParent(parent, trainRows.data(), numTrainRows,
                          testRows.data(), numTestRows, gatherColumns,
                          numGatherColumns,
                          columns.empty() ? NULL : columns.data(),
                          columns.size());

    bartcore::SamplerOptions options =
      optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);
    if (!viewLeafCovariates.empty()) {
      options.leafCovariateColumns = viewLeafCovariates.data();
      options.numLeafCovariates = viewLeafCovariates.size();
    }
    rngs = createChainRngs(control, options.numChains);

    std::unique_ptr<bartcore::SamplerBase> sampler =
      bartcore::createSamplerOverStore(
        std::move(store), response.data(),
        data.weights != NULL ? weights.data() : NULL,
        data.offset != NULL ? offset.data() : NULL, family, data.sigmaEstimate,
        model.sigmaDf, model.sigmaRawScale, options, rngs.data());
    if (sampler == NULL) {
      // R-side resolution validates column legality first, so what lands here
      // is a covariate the handle did not gather raw for at creation
      for (ext_rng* rng : rngs) if (rng != NULL) ext_rng_destroy(rng);
      Rf_error("invalid leaf covariate designation; a view's leaf covariate "
               "columns must be among those the data handle was told to "
               "gather raw values for at creation");
    }

    BartcoreHolder* holder = new BartcoreHolder{
      std::move(sampler), std::move(rngs), control.keepTrainingFits,
      {}, {}, {}, {}, {}, {}, {}};
    // moving the vectors keeps their buffers, so the chains' borrowed
    // pointers stay valid for the holder's lifetime
    holder->ownedResponse = std::move(response);
    holder->ownedWeights = std::move(weights);
    holder->ownedOffset = std::move(offset);
    holder->ownedTestOffset = std::move(testOffset);
    if (!holder->ownedTestOffset.empty())
      holder->sampler->setTestOffset(holder->ownedTestOffset.data());

    SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
    SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
    SEXP result = PROTECT(R_MakeExternalPtr(holder, R_NilValue, protExpr));
    R_RegisterCFinalizerEx(result, holderFinalizer,
                           static_cast<Rboolean>(FALSE));

    UNPROTECT(2);
    return result;
  });
}

// A BCF two-forest sampler; internal, gaussian only (docs/design/bcf.md).
// The model spec is the prognostic forest, bcfParams the treatment forest and
// glue, z the 0/1 treatment.
SEXP bartcore_createBCF(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                        SEXP zExpr, SEXP bcfParamsExpr, SEXP moderatorsExpr) {
  // null-address first, then install; see bartcore_create
  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  BartcoreHolder* holder = bartcore_bridge::createBCFHolder(
    controlExpr, modelExpr, dataExpr, zExpr, bcfParamsExpr, moderatorsExpr);
  R_SetExternalPtrAddr(result, holder);

  UNPROTECT(2);
  return result;
}

// A K-forest multinomial (softmax) sampler; internal, constant-leaf only
// (docs/design/multinomial.md). labels are the 0..K-1 category codes.
SEXP bartcore_createMultinomial(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                                SEXP labelsExpr, SEXP numCategoriesExpr) {
  // null-address first, then install; see bartcore_create
  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  BartcoreHolder* holder = bartcore_bridge::createMultinomialHolder(
    controlExpr, modelExpr, dataExpr, labelsExpr, numCategoriesExpr);
  R_SetExternalPtrAddr(result, holder);

  UNPROTECT(2);
  return result;
}

// The grouped-count entry: counts is an n x K nonnegative integer matrix, the
// trials n_i = sum_k counts_ik derived C-side (docs/design/multinomial.md).
SEXP bartcore_createMultinomialCounts(SEXP controlExpr, SEXP modelExpr,
                                      SEXP dataExpr, SEXP countsExpr,
                                      SEXP numCategoriesExpr) {
  // null-address first, then install; see bartcore_create
  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  BartcoreHolder* holder = bartcore_bridge::createMultinomialCountsHolder(
    controlExpr, modelExpr, dataExpr, countsExpr, numCategoriesExpr);
  R_SetExternalPtrAddr(result, holder);

  UNPROTECT(2);
  return result;
}

SEXP bartcore_setTreatment(SEXP ptrExpr, SEXP zExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t n = holder.sampler->numObservations();
  if (holder.sampler->numForests() < 2)
    Rf_error("bartcore_setTreatment requires a BCF sampler");
  if (static_cast<size_t>(Rf_xlength(zExpr)) != n)
    Rf_error("treatment length must match the number of observations");
  holder.ownedTreatment.resize(n);
  for (size_t i = 0; i < n; ++i)
    holder.ownedTreatment[i] = REAL(zExpr)[i] != 0.0 ? 1.0 : 0.0;
  holder.sampler->setTreatment(holder.ownedTreatment.data());
  return R_NilValue;
}

// The glue on the combining response, one column {a, b0, b1} per chain.
SEXP bartcore_getBCFGlue(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numChains = holder.sampler->numChains();
  SEXP result =
    PROTECT(Rf_allocMatrix(REALSXP, 3, static_cast<int>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    if (!holder.sampler->bcfGlue(c, REAL(result) + 3 * c)) {
      UNPROTECT(1);
      Rf_error("bartcore_getBCFGlue requires a BCF sampler");
    }
  UNPROTECT(1);
  return result;
}

// One forest's internal-scale function values, numObservations x numChains
// (forest 0 prognostic, 1 treatment).
SEXP bartcore_getForestFits(SEXP ptrExpr, SEXP forestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t forestIndex = static_cast<size_t>(Rf_asInteger(forestExpr));
  if (forestIndex >= holder.sampler->numForests())
    Rf_error("forest index out of range");
  size_t n = holder.sampler->numObservations();
  size_t numChains = holder.sampler->numChains();
  SEXP result = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n),
                                       static_cast<int>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    holder.sampler->forestTotalFits(c, forestIndex, REAL(result) + c * n);
  UNPROTECT(1);
  return result;
}

SEXP bartcore_getForestVariableCounts(SEXP ptrExpr, SEXP forestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t forestIndex = static_cast<size_t>(Rf_asInteger(forestExpr));
  if (forestIndex >= holder.sampler->numForests())
    Rf_error("forest index out of range");
  size_t numPredictors = holder.sampler->numPredictors();
  size_t numChains = holder.sampler->numChains();
  // counts alias R integers: variable-use counts never approach 2^31
  SEXP result = PROTECT(Rf_allocMatrix(INTSXP, static_cast<int>(numPredictors),
                                       static_cast<int>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    holder.sampler->forestVariableCounts(
      c, forestIndex,
      reinterpret_cast<uint32_t*>(INTEGER(result)) + c * numPredictors);
  UNPROTECT(1);
  return result;
}

// R_CheckUserInterrupt longjmps when an interrupt is pending; running it
// through R_ToplevelExec catches that jump so the sampler can join its worker
// threads before the interrupt becomes an error (a bare longjmp would strand
// them). Must be called only on the main R thread. R_ToplevelExec returns
// FALSE when the wrapped call jumped, i.e. when an interrupt was pending.
static void bartcore_checkInterrupt(void*) { R_CheckUserInterrupt(); }
static bool bartcore_userInterrupted() {
  return R_ToplevelExec(bartcore_checkInterrupt, nullptr) == FALSE;
}

SEXP bartcore_run(SEXP ptrExpr, SEXP numBurnInExpr, SEXP numSamplesExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  size_t numBurnIn = static_cast<size_t>(Rf_asInteger(numBurnInExpr));
  size_t numSamples = static_cast<size_t>(Rf_asInteger(numSamplesExpr));

  size_t numObservations = sampler.numObservations();
  size_t numPredictors = sampler.numPredictors();
  size_t numTestObservations = sampler.numTestObservations();
  size_t numChains = sampler.numChains();
  int numSamplesInt = static_cast<int>(numSamples);
  int numChainsInt = static_cast<int>(numChains);
  // a multi-location combiner (multinomial: K softmax channels) inserts a
  // location dimension between the observations and the samples; L = 1 keeps
  // the exact n x numSamples (x numChains) shape and byte layout every other
  // model relies on.
  size_t numLocations = sampler.numReportedLocations();
  auto allocFitsArray = [&](size_t leadingDim) -> SEXP {
    int dims[4];
    int numDims = 0;
    dims[numDims++] = static_cast<int>(leadingDim);
    dims[numDims++] = static_cast<int>(numLocations);
    dims[numDims++] = numSamplesInt;
    if (numChains > 1) dims[numDims++] = numChainsInt;
    R_xlen_t total = 1;
    for (int d = 0; d < numDims; ++d) total *= dims[d];
    SEXP arr = PROTECT(Rf_allocVector(REALSXP, total));
    SEXP dimExpr = Rf_allocVector(INTSXP, numDims);
    for (int d = 0; d < numDims; ++d) INTEGER(dimExpr)[d] = dims[d];
    Rf_setAttrib(arr, R_DimSymbol, dimExpr);
    UNPROTECT(1);
    return arr;
  };
  // the varcount channel widens on its own forest axis (multinomial: K category
  // forests), inserting a forest dimension between the predictors and the
  // samples exactly as the fits seam inserts locations; count 1 keeps the exact
  // numPredictors x numSamples (x numChains) INTSXP shape and byte layout every
  // other model relies on.
  size_t numVCForests = sampler.numVariableCountForests();
  auto allocVarcountArray = [&]() -> SEXP {
    int dims[4];
    int numDims = 0;
    dims[numDims++] = static_cast<int>(numPredictors);
    dims[numDims++] = static_cast<int>(numVCForests);
    dims[numDims++] = numSamplesInt;
    if (numChains > 1) dims[numDims++] = numChainsInt;
    R_xlen_t total = 1;
    for (int d = 0; d < numDims; ++d) total *= dims[d];
    SEXP arr = PROTECT(Rf_allocVector(INTSXP, total));
    SEXP dimExpr = Rf_allocVector(INTSXP, numDims);
    for (int d = 0; d < numDims; ++d) INTEGER(dimExpr)[d] = dims[d];
    Rf_setAttrib(arr, R_DimSymbol, dimExpr);
    UNPROTECT(1);
    return arr;
  };

  if (numSamples == 0) {
    bartcore::Results empty;
    GetRNGstate();
    bool cancelled = sampler.run(numBurnIn, 0, empty, bartcore_userInterrupted);
    PutRNGstate();
    if (cancelled) Rf_error("sampler run interrupted");
    return R_NilValue;
  }

  // ordinal reports its K-1 cutpoints in an extra channel appended after ranef;
  // every other family carries none, so the list keeps its 8 slots and byte-for-
  // byte layout. numCutpoints == 0 off ordinal.
  size_t numCutpoints = sampler.numCutpoints();
  bool hasCutpoints = numCutpoints > 0;
  int numResultSlots = hasCutpoints ? 9 : 8;

  // several chains add a trailing chain dimension, as the classic engine's
  // results do. Every column roots in the protected container the moment it
  // is allocated (installResult), so there is no hand-counted PROTECT stack
  // to keep in sync with the slot list.
  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, numResultSlots));
  SEXP sigmaExpr = installResult(
    resultExpr, 0,
    numChains == 1
      ? Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples))
      : Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  SEXP trainExpr = installResult(
    resultExpr, 1,
    !holder.keepTrainingFits
      ? R_NilValue
      : numLocations > 1
        ? allocFitsArray(numObservations)
        : numChains == 1
          ? Rf_allocMatrix(REALSXP, static_cast<int>(numObservations),
                           numSamplesInt)
          : Rf_alloc3DArray(REALSXP, static_cast<int>(numObservations),
                            numSamplesInt, numChainsInt));
  SEXP testExpr = installResult(
    resultExpr, 2,
    numTestObservations == 0
      ? R_NilValue
      : numLocations > 1
        ? allocFitsArray(numTestObservations)
        : numChains == 1
          ? Rf_allocMatrix(REALSXP, static_cast<int>(numTestObservations),
                           numSamplesInt)
          : Rf_alloc3DArray(REALSXP, static_cast<int>(numTestObservations),
                            numSamplesInt, numChainsInt));
  SEXP varcountExpr = installResult(
    resultExpr, 3,
    numVCForests > 1
      ? allocVarcountArray()
      : numChains == 1
        ? Rf_allocMatrix(INTSXP, static_cast<int>(numPredictors), numSamplesInt)
        : Rf_alloc3DArray(INTSXP, static_cast<int>(numPredictors),
                          numSamplesInt, numChainsInt));
  SEXP kExpr = installResult(
    resultExpr, 4,
    !sampler.kIsSampled()
      ? R_NilValue
      : numChains == 1
        ? Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples))
        : Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  SEXP varprobsExpr = installResult(
    resultExpr, 5,
    !sampler.usesDart()
      ? R_NilValue
      : numChains == 1
        ? Rf_allocMatrix(REALSXP, static_cast<int>(numPredictors),
                         numSamplesInt)
        : Rf_alloc3DArray(REALSXP, static_cast<int>(numPredictors),
                          numSamplesInt, numChainsInt));
  size_t numGroups = sampler.numGroups();
  SEXP tauExpr = installResult(
    resultExpr, 6,
    numGroups == 0
      ? R_NilValue
      : numChains == 1
        ? Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples))
        : Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  SEXP ranefExpr = installResult(
    resultExpr, 7,
    numGroups == 0
      ? R_NilValue
      : numChains == 1
        ? Rf_allocMatrix(REALSXP, static_cast<int>(numGroups), numSamplesInt)
        : Rf_alloc3DArray(REALSXP, static_cast<int>(numGroups),
                          numSamplesInt, numChainsInt));
  SEXP cutpointsExpr = !hasCutpoints
    ? R_NilValue
    : installResult(
        resultExpr, 8,
        numChains == 1
          ? Rf_allocMatrix(REALSXP, static_cast<int>(numCutpoints),
                           numSamplesInt)
          : Rf_alloc3DArray(REALSXP, static_cast<int>(numCutpoints),
                            numSamplesInt, numChainsInt));

  std::vector<std::uint32_t> variableCounts(numPredictors * numVCForests *
                                            numSamples * numChains);

  bartcore::Results results;
  results.sigma = REAL(sigmaExpr);
  results.trainingFits = holder.keepTrainingFits ? REAL(trainExpr) : NULL;
  results.testFits = numTestObservations > 0 ? REAL(testExpr) : NULL;
  results.variableCounts = variableCounts.data();
  results.k = sampler.kIsSampled() ? REAL(kExpr) : NULL;
  results.splitProbabilities = sampler.usesDart() ? REAL(varprobsExpr) : NULL;
  results.tau = numGroups > 0 ? REAL(tauExpr) : NULL;
  results.groupEffects = numGroups > 0 ? REAL(ranefExpr) : NULL;
  // one per-observation fits channel per reported location; 1 for every
  // additive model, K for multinomial. The location stride drives the
  // chain-major slabbing (multiple chains).
  results.numReportedLocations = numLocations;
  // one varcount slab per reported forest; 1 for every additive model, K for
  // multinomial's per-category counts. Drives the chain-major varcount stride.
  results.numVariableCountForests = numVCForests;
  // K-1 cutpoints per sample for ordinal, none otherwise; drives the chain-major
  // cutpoint stride
  results.cutpoints = hasCutpoints ? REAL(cutpointsExpr) : NULL;
  results.numCutpoints = numCutpoints;

  GetRNGstate();
  bool cancelled = sampler.run(numBurnIn, numSamples, results,
                               bartcore_userInterrupted);
  PutRNGstate();
  if (cancelled) {
    std::vector<std::uint32_t>().swap(variableCounts);  // free before longjmp
    Rf_error("sampler run interrupted");
  }

  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < numPredictors * numVCForests * numSamples * numChains;
       ++i)
    varcountOut[i] = static_cast<int>(variableCounts[i]);
  // free the copied-out counts before the names block, which can OOM-longjmp
  std::vector<std::uint32_t>().swap(variableCounts);

  // named as the classic engine's run results are, so the engines are
  // drop-in replacements for each other; varprobs, tau, and ranef are
  // bartcore extensions. The names vector roots through the container's
  // attribute before the mkChar allocations fill it.
  SEXP namesExpr = Rf_allocVector(STRSXP, numResultSlots);
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));
  SET_STRING_ELT(namesExpr, 5, Rf_mkChar("varprobs"));
  SET_STRING_ELT(namesExpr, 6, Rf_mkChar("tau"));
  SET_STRING_ELT(namesExpr, 7, Rf_mkChar("ranef"));
  if (hasCutpoints) SET_STRING_ELT(namesExpr, 8, Rf_mkChar("cutpoints"));

  UNPROTECT(1);
  return resultExpr;
}

// rbart_vi's custom-prior Gibbs sampler: one run with a per-sweep R closure in
// place of a run(0, 1) per kept sample. results is a named list of caller-owned
// per-sweep buffers the engine fills (sigma, train, test, k, varprobs reals; an
// integer varcount whose storage aliases the engine's uint32 slots, valid
// because variable-use counts never approach 2^31); a null slot skips that
// channel. Single chain only, so the closure runs inline.
//
// No GetRNGstate/PutRNGstate: the chain's generator never touches R's stream,
// while the closure draws from it (the random intercepts and the tau slice
// sampler), so R must own .Random.seed throughout. The closure is evaluated
// under R_tryEval so an error cannot longjmp across Chain::run's C++ frames; it
// becomes a cooperative stop, re-raised in R from the closure's own handler.
SEXP bartcore_runWithCallback(SEXP ptrExpr, SEXP numBurnInExpr,
                              SEXP numSamplesExpr, SEXP resultsExpr,
                              SEXP callbackExpr, SEXP rhoExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  if (sampler.numChains() != 1)
    Rf_error("dbarts_bartcore_runWithCallback requires a single chain");
  if (!Rf_isFunction(callbackExpr)) Rf_error("callback must be a function");

  size_t numBurnIn = static_cast<size_t>(Rf_asInteger(numBurnInExpr));
  size_t numSamples = static_cast<size_t>(Rf_asInteger(numSamplesExpr));

  SEXP sigmaExpr = getListElement(resultsExpr, "sigma");
  SEXP trainExpr = getListElement(resultsExpr, "train");
  SEXP testExpr = getListElement(resultsExpr, "test");
  SEXP varcountExpr = getListElement(resultsExpr, "varcount");
  SEXP kExpr = getListElement(resultsExpr, "k");
  SEXP varprobsExpr = getListElement(resultsExpr, "varprobs");

  bartcore::Results results;
  results.sigma = Rf_isNull(sigmaExpr) ? NULL : REAL(sigmaExpr);
  results.trainingFits = Rf_isNull(trainExpr) ? NULL : REAL(trainExpr);
  results.testFits = Rf_isNull(testExpr) ? NULL : REAL(testExpr);
  results.variableCounts = Rf_isNull(varcountExpr)
    ? NULL : reinterpret_cast<uint32_t*>(INTEGER(varcountExpr));
  results.k = Rf_isNull(kExpr) ? NULL : REAL(kExpr);
  results.splitProbabilities = Rf_isNull(varprobsExpr) ? NULL : REAL(varprobsExpr);
  // single chain here, so the location stride only shapes the fits buffers the
  // caller allocated; 1 for every model today (n x numSamples)
  results.numReportedLocations = sampler.numReportedLocations();
  // rbart_vi's caller-owned varcount buffer is single-forest (rbart is never
  // multinomial); 1 keeps the aliased numPredictors-per-sweep layout
  results.numVariableCountForests = sampler.numVariableCountForests();

  bool callbackErrored = false;  // an error escaped the closure (R_tryEval)
  bool closureStopped = false;   // the closure returned TRUE (self-caught stop)
  bartcore::SweepCallback onSweep =
    [&](size_t, size_t sweepIndex, bool) -> bool {
      SEXP call = PROTECT(Rf_lang2(
        callbackExpr, Rf_ScalarInteger(static_cast<int>(sweepIndex))));
      int errorOccurred = 0;
      SEXP res = R_tryEval(call, rhoExpr, &errorOccurred);
      bool stop = errorOccurred == 0 && res != R_NilValue &&
                  Rf_asLogical(res) == TRUE;
      UNPROTECT(1);
      if (errorOccurred) { callbackErrored = true; return true; }
      closureStopped = stop;
      return stop;
    };

  bool cancelled = sampler.run(numBurnIn, numSamples, results,
                               bartcore_userInterrupted, onSweep);
  if (callbackErrored)
    Rf_error("error evaluating the rbart_vi sweep callback");
  if (cancelled && !closureStopped) Rf_error("sampler run interrupted");
  return R_NilValue;
}

SEXP bartcore_sampleTreesFromPrior(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  GetRNGstate();
  holder.sampler->sampleTreesFromPrior();
  PutRNGstate();
  return R_NilValue;
}

SEXP bartcore_sampleNodeParametersFromPrior(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  GetRNGstate();
  holder.sampler->sampleNodeParametersFromPrior();
  PutRNGstate();
  return R_NilValue;
}

SEXP bartcore_growFromRoot(SEXP ptrExpr, SEXP numSweepsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  // constant leaf only in v1 (the scan's closed-form marginal); the R5 method
  // refuses first, this is the engine-side backstop
  if (holder.sampler->numLeafCovariates() != 0 ||
      holder.sampler->usesFunctionLeaves())
    Rf_error("grow-from-root warm start is only available for the "
             "constant-leaf model");
  int numSweeps = Rf_asInteger(numSweepsExpr);
  if (numSweeps == NA_INTEGER || numSweeps <= 0)
    Rf_error("n.sweeps must be a positive integer");
  GetRNGstate();
  holder.sampler->growFromRoot(static_cast<size_t>(numSweeps));
  PutRNGstate();
  return R_NilValue;
}

SEXP bartcore_setOffset(SEXP ptrExpr, SEXP offsetExpr, SEXP updateScaleExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (!Rf_isNull(offsetExpr) &&
      (!Rf_isReal(offsetExpr) ||
       static_cast<size_t>(Rf_xlength(offsetExpr)) !=
         holder.sampler->numObservations()))
    Rf_error("length of replacement offset is not equal to number of observations");
  const double* offset = Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr);
  holder.sampler->setOffset(offset, Rf_asLogical(updateScaleExpr) == TRUE);
  retain(ptrExpr, PROT_OFFSET, offsetExpr);
  return R_NilValue;
}

SEXP bartcore_setResponse(SEXP ptrExpr, SEXP yExpr, SEXP updateScaleExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refuseMultiForestMutation(*holder.sampler, "bartcore_setResponse");
  if (holder.sampler->numGroups() > 0)
    Rf_error("grouped random effects fix the response at creation; make a "
             "new sampler instead");
  if (!Rf_isReal(yExpr) ||
      static_cast<size_t>(Rf_xlength(yExpr)) !=
        holder.sampler->numObservations())
    Rf_error("y must be of length equal to %lu",
             static_cast<unsigned long>(holder.sampler->numObservations()));
  GetRNGstate(); // probit latent redraw
  holder.sampler->setResponse(REAL(yExpr),
                              Rf_asLogical(updateScaleExpr) == TRUE);
  PutRNGstate();
  retain(ptrExpr, PROT_RESPONSE, yExpr);
  return R_NilValue;
}

SEXP bartcore_setSigma(SEXP ptrExpr, SEXP sigmaExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  holder.sampler->setSigma(Rf_asReal(sigmaExpr));
  return R_NilValue;
}

SEXP bartcore_setData(SEXP ptrExpr, SEXP dataExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  refusePredictorMutation(sampler, "bartcore_setData");
  refuseMultiForestMutation(sampler, "bartcore_setData");
  if (sampler.numGroups() > 0)
    Rf_error("grouped random effects fix the data at creation; make a new "
             "sampler instead");
  if (sampler.family() == bartcore::ResponseFamily::aft)
    Rf_error("aft (survival) models fix the censoring structure at creation; "
             "make a new sampler instead");

  if (!Rf_inherits(dataExpr, "dbartsData"))
    Rf_error("'data' argument to bartcore_setData not of class 'dbartsData'");

  return unwindProtect([&, data = ParsedData{}]() mutable -> SEXP {
    parseData(data, dataExpr);
    if (data.xIsSparse || data.xIsMixed)
      Rf_error("bartcore setData requires a dense predictor matrix; sparse "
               "predictors fix the design at creation");
    if (data.testIsMixed)
      Rf_error("bartcore setData requires a dense test matrix; a sparse test "
               "set fixes the design at creation");

    if (data.numPredictors != sampler.numPredictors())
      Rf_error("bartcore setData requires the same predictors");
    if (data.weights != NULL) refuseBinaryWeightChange(sampler);
    for (size_t j = 0; j < data.numPredictors; ++j) {
      bool wasCategorical = sampler.data().types[j] ==
                            bartcore::ColumnType::categorical;
      bool isCategorical = data.columnTypes[j] ==
                           bartcore::ColumnType::categorical;
      if (isCategorical != wasCategorical)
        Rf_error("bartcore setData requires the same predictor types");
      if (!wasCategorical) continue;
      // category counts are fixed at creation; new values must be existing
      // codes, in the training and test data both
      for (size_t i = 0; i < data.numObservations; ++i)
        if (!sampler.data().categoricalValueIsValid(
              j, data.x[i + j * data.numObservations]))
          Rf_error("categorical predictor values must be existing category "
                   "codes");
      for (size_t i = 0; i < data.numTestObservations; ++i)
        if (!sampler.data().categoricalValueIsValid(
              j, data.x_test[i + j * data.numTestObservations]))
          Rf_error("categorical predictor values must be existing category "
                   "codes");
    }

    sampler.setData(data.x, data.y, data.numObservations, data.weights,
                    data.offset, data.x_test, data.numTestObservations,
                    data.testOffset);

    // the new spec is the creation contract and the call-time raw source; the
    // engine re-quantized it and retains no predictor pointer
    retain(ptrExpr, PROT_DATA, dataExpr);
    retain(ptrExpr, PROT_RESPONSE, R_NilValue);
    retain(ptrExpr, PROT_OFFSET, R_NilValue);

    return R_NilValue;
  });
}

SEXP bartcore_setTestPredictor(SEXP ptrExpr, SEXP xTestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (Rf_isNull(xTestExpr)) {
    // removal: back to the no-test-data state, offset included
    holder.sampler->setTestPredictors(NULL, 0);
    holder.sampler->setTestOffset(NULL);
    retain(ptrExpr, PROT_TEST_OFFSET, R_NilValue);
    return R_NilValue;
  }
  refuseBCFTestSurface(*holder.sampler, "bartcore_setTestPredictor");
  if (Rf_inherits(xTestExpr, "dbartsMixedMatrix"))
    // whole-object replacement by a mixed/sparse container: parse against the
    // training cut grid, then rebuild the typed test store. Both the parse and
    // the leaf-covariate refusal precede any store change, so a rejected
    // container leaves the sampler in its prior state. The parsed sources own
    // buffers, so the unwind-protected scope frees them on the error jump.
    return unwindProtect([&, parsed = ParsedTestContainer{}]() mutable -> SEXP {
      parseTestContainer(parsed, xTestExpr, holder.sampler->numPredictors(),
                         holder.sampler->data().types.data());
      if (holder.sampler->data().testOffset != NULL &&
          parsed.numTestObservations != holder.sampler->numTestObservations())
        Rf_error("test offset length would no longer match; set the predictors "
                 "and offset together");
      if (!installTestContainer(*holder.sampler, parsed))
        Rf_error("a leaf covariate column cannot be a sparse test column; "
                 "supply it as a dense test column");
      return R_NilValue;
    });
  SEXP dims = Rf_getAttrib(xTestExpr, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[1]) != holder.sampler->numPredictors())
    Rf_error("bartcore_setTestPredictor requires a matrix with matching columns");
  size_t numTestObservations = static_cast<size_t>(INTEGER(dims)[0]);
  if (holder.sampler->data().testOffset != NULL &&
      numTestObservations != holder.sampler->numTestObservations())
    Rf_error("test offset length would no longer match; set the predictors "
             "and offset together");
  for (size_t j = 0; j < holder.sampler->numPredictors(); ++j)
    validateColumnValues(holder.sampler->data(), j,
                         REAL(xTestExpr) + j * numTestObservations,
                         numTestObservations);
  // buildTest copies the test values into owned storage, so nothing is pinned
  holder.sampler->setTestPredictors(REAL(xTestExpr), numTestObservations);
  return R_NilValue;
}

// offsetExpr may be null (clear); length must match the current test data.
SEXP bartcore_setTestOffset(SEXP ptrExpr, SEXP offsetExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (Rf_isNull(offsetExpr)) {
    holder.sampler->setTestOffset(NULL);
    retain(ptrExpr, PROT_TEST_OFFSET, R_NilValue);
    return R_NilValue;
  }
  refuseBCFTestSurface(*holder.sampler, "bartcore_setTestOffset");
  if (holder.sampler->numTestObservations() == 0)
    Rf_error("cannot set a test offset without test predictors");
  if (!Rf_isReal(offsetExpr) ||
      static_cast<size_t>(Rf_xlength(offsetExpr)) !=
        holder.sampler->numTestObservations())
    Rf_error("length of test offset must equal number of test observations");
  holder.sampler->setTestOffset(REAL(offsetExpr));
  retain(ptrExpr, PROT_TEST_OFFSET, offsetExpr);
  return R_NilValue;
}

// The combined form the row count may change through; offsetExpr null
// clears the offset alongside the new predictors. A null test matrix
// removes the test data (its offset must be null as well).
SEXP bartcore_setTestPredictorAndOffset(SEXP ptrExpr, SEXP xTestExpr,
                                        SEXP offsetExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (Rf_isNull(xTestExpr)) {
    if (!Rf_isNull(offsetExpr))
      Rf_error("when test matrix is NULL, test offset must be as well");
    holder.sampler->setTestPredictors(NULL, 0);
    holder.sampler->setTestOffset(NULL);
    retain(ptrExpr, PROT_TEST_OFFSET, R_NilValue);
    return R_NilValue;
  }
  refuseBCFTestSurface(*holder.sampler, "bartcore_setTestPredictorAndOffset");
  if (Rf_inherits(xTestExpr, "dbartsMixedMatrix"))
    // whole-object container replacement plus its offset; parse and the
    // leaf-covariate refusal precede the store change, so a rejected container
    // leaves the sampler untouched
    return unwindProtect([&, parsed = ParsedTestContainer{}]() mutable -> SEXP {
      parseTestContainer(parsed, xTestExpr, holder.sampler->numPredictors(),
                         holder.sampler->data().types.data());
      if (!Rf_isNull(offsetExpr) &&
          (!Rf_isReal(offsetExpr) ||
           static_cast<size_t>(Rf_xlength(offsetExpr)) !=
             parsed.numTestObservations))
        Rf_error("length of test offset must equal number of rows in test "
                 "matrix");
      if (!installTestContainer(*holder.sampler, parsed))
        Rf_error("a leaf covariate column cannot be a sparse test column; "
                 "supply it as a dense test column");
      holder.sampler->setTestOffset(Rf_isNull(offsetExpr) ? NULL
                                                          : REAL(offsetExpr));
      retain(ptrExpr, PROT_TEST_OFFSET, offsetExpr);
      return R_NilValue;
    });
  SEXP dims = Rf_getAttrib(xTestExpr, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[1]) != holder.sampler->numPredictors())
    Rf_error("bartcore_setTestPredictorAndOffset requires a matrix with "
             "matching columns");
  size_t numTestObservations = static_cast<size_t>(INTEGER(dims)[0]);
  if (!Rf_isNull(offsetExpr) &&
      (!Rf_isReal(offsetExpr) ||
       static_cast<size_t>(Rf_xlength(offsetExpr)) != numTestObservations))
    Rf_error("length of test offset must equal number of rows in test matrix");
  for (size_t j = 0; j < holder.sampler->numPredictors(); ++j)
    validateColumnValues(holder.sampler->data(), j,
                         REAL(xTestExpr) + j * numTestObservations,
                         numTestObservations);

  // buildTest owns the test values; only the borrowed test offset is pinned
  holder.sampler->setTestPredictors(REAL(xTestExpr), numTestObservations);
  holder.sampler->setTestOffset(Rf_isNull(offsetExpr) ? NULL
                                                      : REAL(offsetExpr));
  retain(ptrExpr, PROT_TEST_OFFSET, offsetExpr);
  return R_NilValue;
}

// Case weights, like the classic engine's setWeights: a pointer swap with
// nothing rescaled; refused for binary responses, whose reference-engine
// weighting was incorrect and was stripped rather than ported.
SEXP bartcore_setWeights(SEXP ptrExpr, SEXP weightsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refuseMultiForestMutation(*holder.sampler, "bartcore_setWeights");
  refuseBinaryWeightChange(*holder.sampler);
  if (!Rf_isReal(weightsExpr) ||
      static_cast<size_t>(Rf_xlength(weightsExpr)) !=
        holder.sampler->numObservations())
    Rf_error("length of weights must equal number of observations");
  const double* weights = REAL(weightsExpr);
  for (size_t i = 0; i < holder.sampler->numObservations(); ++i)
    if (!(weights[i] >= 0.0))
      Rf_error("weights must be non-negative");
  holder.sampler->setWeights(weights);
  retain(ptrExpr, PROT_WEIGHTS, weightsExpr);
  return R_NilValue;
}

// Between-run reconfiguration, the classic engine's setControl. The R side
// refuses changes to the engine, rng, and cut settings; chain and tree
// counts shape live storage, so they are re-checked here.
SEXP bartcore_setControl(SEXP ptrExpr, SEXP controlExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  if (!Rf_inherits(controlExpr, "dbartsControl"))
    Rf_error("'control' argument to bartcore_setControl not of class "
             "'dbartsControl'");

  ParsedControl control;
  parseControl(control, controlExpr);

  if (control.numChains != sampler.numChains())
    Rf_error("the bartcore engine cannot change the number of chains of an "
             "existing sampler");
  if (control.numTrees != sampler.numTrees())
    Rf_error("the bartcore engine cannot change the number of trees of an "
             "existing sampler");

  holder.keepTrainingFits = control.keepTrainingFits;
  sampler.setNumThreads(control.numThreads);
  sampler.setNumThin(control.treeThinningRate);
  sampler.setVerbose(control.verbose, control.printEvery);
  sampler.setTreeStorage(control.keepTrees, control.defaultNumSamples);

  return R_NilValue;
}

/// Prior replacement, the classic engine's setModel; installing a model
/// before any run matches creating with it.
SEXP bartcore_setModel(SEXP ptrExpr, SEXP modelExpr, SEXP controlExpr,
                       SEXP dataExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  if (!Rf_inherits(modelExpr, "dbartsModel"))
    Rf_error("'model' argument to bartcore_setModel not of class "
             "'dbartsModel'");

  (void) controlExpr; // arity fixed by the call table; nothing read from it

  return unwindProtect([&, model = ParsedModel{}]() mutable -> SEXP {
    parseModel(model, modelExpr, sampler.numPredictors());

    // the leaf model is a template instantiation: the designation and its
    // kind are fixed at creation, so a replacement prior must carry the same
    bool designationMatches =
      model.leafCovariateColumns.size() == sampler.numLeafCovariates() &&
      model.gpLeaves == sampler.usesFunctionLeaves();
    for (size_t j = 0; designationMatches &&
         j < model.leafCovariateColumns.size(); ++j)
      designationMatches =
        model.leafCovariateColumns[j] == sampler.leafCovariateColumns()[j];
    if (!designationMatches)
      Rf_error("the leaf covariate designation is fixed when a sampler is "
               "created; make a new sampler instead");

    // gaussian and aft both draw sigma from a variance prior; the binary
    // families fix it, so setModel leaves their sigma alone
    bool drawsSigma = sampler.family() == bartcore::ResponseFamily::gaussian ||
                      sampler.family() == bartcore::ResponseFamily::aft;

    bartcore::ModelParameters parameters;
    parameters.base = model.base;
    parameters.power = model.power;
    parameters.splitProbabilities = model.splitProbabilities;
    parameters.birthOrDeathProbability = model.birthOrDeathProbability;
    parameters.swapProbability = model.swapProbability;
    parameters.changeProbability = model.changeProbability;
    parameters.birthProbability = model.birthProbability;
    parameters.nodeScale = model.nodeScale;
    parameters.updateK = model.updateK;
    if (parameters.updateK) {
      parameters.kHyperprior.degreesOfFreedom = model.kDf;
      parameters.kHyperprior.scale = model.kScale;
    } else {
      parameters.k = model.k;
    }
    if (drawsSigma) {
      if (model.sigmaIsFixed) {
        // documented semantics: fixed(value) holds the residual variance
        parameters.sigmaIsFixed = true;
        parameters.sigmaEstimate = std::sqrt(model.fixedSigmaSq);
      } else {
        parameters.sigmaEstimate = rc_getDouble(
          Rf_getAttrib(dataExpr, Rf_install("sigma")), "sigma estimate",
          RC_LENGTH | RC_EQ, rc_asRLength(1), RC_NA | RC_YES,
          RC_VALUE | RC_GT, 0.0, RC_END);
        parameters.sigmaDf = model.sigmaDf;
        parameters.sigmaRawScale = model.sigmaRawScale;
      }
    }

    // split probabilities are copied per chain before the model goes away
    sampler.setModel(parameters);

    return R_NilValue;
  });
}

SEXP bartcore_isValidPointer(SEXP ptrExpr) {
  return Rf_ScalarLogical(R_ExternalPtrAddr(ptrExpr) != NULL ? TRUE : FALSE);
}

SEXP bartcore_getSigmas(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numChains = holder.sampler->numChains();
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    REAL(result)[c] = holder.sampler->sigma(c);
  UNPROTECT(1);
  return result;
}

SEXP bartcore_getSumsOfSquaredResiduals(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numChains = holder.sampler->numChains();
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    REAL(result)[c] = holder.sampler->sumOfSquaredResiduals(c);
  UNPROTECT(1);
  return result;
}

// Predictor mutation. The sampler borrows a full replacement matrix
// (R/bartcore.R retains it on success); column and per-observation updates
// write in place into the matrix the sampler currently borrows, aliasing the
// R-side data like the classic engine does.

SEXP bartcore_setPredictor(SEXP ptrExpr, SEXP xExpr, SEXP forceUpdateExpr,
                           SEXP updateCutPointsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refusePredictorMutation(*holder.sampler, "bartcore_setPredictor");
  if (Rf_asLogical(forceUpdateExpr) != TRUE)
    refuseMultiForestTransactionalUpdate(*holder.sampler,
                                         "bartcore_setPredictor");
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[0]) != holder.sampler->numObservations() ||
      static_cast<size_t>(INTEGER(dims)[1]) != holder.sampler->numPredictors())
    Rf_error("bartcore_setPredictor requires a matrix with matching dimensions");

  for (size_t j = 0; j < holder.sampler->numPredictors(); ++j)
    validateColumnValues(holder.sampler->data(), j,
                         REAL(xExpr) + j * holder.sampler->numObservations(),
                         holder.sampler->numObservations());

  bartcore::PredictorUpdateResult result = holder.sampler->setPredictor(
    REAL(xExpr), Rf_asLogical(forceUpdateExpr) == TRUE,
    Rf_asLogical(updateCutPointsExpr) == TRUE);
  if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
    Rf_error("number of induced cut points in new predictor less than "
             "previous: old splits would be invalid");
  // the engine re-quantized xExpr into owned codes and retains no pointer; the
  // R method reassigns sampler$data@x, which keeps the current matrix alive
  return Rf_ScalarLogical(
    result == bartcore::PredictorUpdateResult::accepted ? TRUE : FALSE);
}

SEXP bartcore_updatePredictor(SEXP ptrExpr, SEXP xExpr, SEXP columnsExpr,
                              SEXP forceUpdateExpr, SEXP updateCutPointsExpr) {
  return unwindProtect([&, columns = std::vector<size_t>{}]() mutable -> SEXP {
    BartcoreHolder& holder(holderFromExpression(ptrExpr));
    refusePredictorMutation(*holder.sampler, "bartcore_updatePredictor");
    if (Rf_asLogical(forceUpdateExpr) != TRUE)
      refuseMultiForestTransactionalUpdate(*holder.sampler,
                                           "bartcore_updatePredictor");
    size_t numObservations = holder.sampler->numObservations();
    size_t numPredictors = holder.sampler->numPredictors();

    size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
    if (numColumns == 0 ||
        static_cast<size_t>(Rf_xlength(xExpr)) != numObservations * numColumns)
      Rf_error("bartcore_updatePredictor requires numObservations values per column");

    columns.resize(numColumns);
    for (size_t k = 0; k < numColumns; ++k) {
      int column = INTEGER(columnsExpr)[k];
      if (column < 1 || static_cast<size_t>(column) > numPredictors)
        Rf_error("bartcore_updatePredictor column out of range");
      columns[k] = static_cast<size_t>(column - 1);
    }

    for (size_t k = 0; k < numColumns; ++k)
      validateColumnValues(holder.sampler->data(), columns[k],
                           REAL(xExpr) + k * numObservations, numObservations);

    bartcore::PredictorUpdateResult result = holder.sampler->updatePredictor(
      REAL(xExpr), columns.data(), numColumns,
      Rf_asLogical(forceUpdateExpr) == TRUE,
      Rf_asLogical(updateCutPointsExpr) == TRUE);
    if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
      Rf_error("number of induced cut points in new predictor less than "
               "previous: old splits would be invalid");
    return Rf_ScalarLogical(
      result == bartcore::PredictorUpdateResult::accepted ? TRUE : FALSE);
  });
}

SEXP bartcore_setCutPoints(SEXP ptrExpr, SEXP cutPointsExpr,
                           SEXP columnsExpr, SEXP currentPredictorsExpr) {
  // the by-value scratch is freed by the wrapper on the Rf_error jump
  return unwindProtect([&, cutPoints = std::vector<const double*>{},
                        numCutPoints = std::vector<std::uint32_t>{},
                        columns = std::vector<size_t>{}]() mutable -> SEXP {
    BartcoreHolder& holder(holderFromExpression(ptrExpr));
    refuseRequantizeWithoutSource(*holder.sampler, "bartcore_setCutPoints");
    size_t numPredictors = holder.sampler->numPredictors();
    // dense columns re-quantize from the supplied data@x; CSC/mixed columns
    // read their retained slices, so a non-matrix source is passed as null
    const double* currentPredictors =
      Rf_isReal(currentPredictorsExpr) ? REAL(currentPredictorsExpr) : NULL;

    size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
    if (numColumns == 0 ||
        static_cast<size_t>(Rf_xlength(cutPointsExpr)) != numColumns)
      Rf_error("bartcore_setCutPoints requires one cut point vector per column");

    cutPoints.resize(numColumns);
    numCutPoints.resize(numColumns);
    columns.resize(numColumns);
    for (size_t k = 0; k < numColumns; ++k) {
      int column = INTEGER(columnsExpr)[k];
      if (column < 1 || static_cast<size_t>(column) > numPredictors)
        Rf_error("bartcore_setCutPoints column out of range");
      columns[k] = static_cast<size_t>(column - 1);
      if (holder.sampler->data().types[columns[k]] ==
          bartcore::ColumnType::categorical)
        Rf_error("cannot set cut points for a categorical predictor");

      SEXP cutsExpr = VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(k));
      R_xlen_t numCuts = Rf_xlength(cutsExpr);
      if (numCuts > 65535)  // codes must fit xint_t, including numCuts itself
        Rf_error("bartcore_setCutPoints cut point vector too long");
      const double* cuts = REAL(cutsExpr);
      for (R_xlen_t i = 1; i < numCuts; ++i)
        if (cuts[i] <= cuts[i - 1])
          Rf_error("bartcore_setCutPoints requires strictly increasing cut "
                   "points");
      cutPoints[k] = cuts;
      numCutPoints[k] = static_cast<std::uint32_t>(numCuts);
    }

    holder.sampler->setCutPoints(cutPoints.data(), numCutPoints.data(),
                                 columns.data(), numColumns,
                                 currentPredictors);
    return R_NilValue;
  });
}

SEXP bartcore_updatePredictorPerObservation(SEXP ptrExpr, SEXP xExpr,
                                            SEXP columnExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refusePredictorMutation(
    *holder.sampler, "bartcore_updatePredictorPerObservation");
  refuseMultiForestTransactionalUpdate(
    *holder.sampler, "bartcore_updatePredictorPerObservation");
  size_t numObservations = holder.sampler->numObservations();

  if (static_cast<size_t>(Rf_xlength(xExpr)) != numObservations)
    Rf_error("bartcore_updatePredictorPerObservation requires one value per "
             "observation");
  int column = Rf_asInteger(columnExpr);
  if (column < 1 ||
      static_cast<size_t>(column) > holder.sampler->numPredictors())
    Rf_error("bartcore_updatePredictorPerObservation column out of range");
  validateColumnValues(holder.sampler->data(),
                       static_cast<size_t>(column - 1), REAL(xExpr),
                       numObservations);

  std::unique_ptr<bool[]> installed(new bool[numObservations]);

  GetRNGstate();  // scan-order permutation
  bool treesAreValid = holder.sampler->updatePredictorPerObservation(
    REAL(xExpr), static_cast<size_t>(column - 1), installed.get());
  PutRNGstate();

  // The sequential guard admits no empty leaves, so an invalid rebuild is an
  // internal invariant violation; fail loudly rather than leave the sampler
  // in an invalid state.
  if (!treesAreValid) {
    installed.reset();  // free before longjmp
    Rf_error("bartcore updatePredictorPerObservation produced a tree with an "
             "empty leaf");
  }

  SEXP result = PROTECT(
    Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(numObservations)));
  for (size_t i = 0; i < numObservations; ++i)
    LOGICAL(result)[i] = installed[i] ? TRUE : FALSE;
  UNPROTECT(1);
  return result;
}

SEXP bartcore_updatePredictorPerObservationJointly(SEXP ptrsExpr, SEXP xExpr,
                                                   SEXP columnsExpr) {
  return unwindProtect(
    [&, samplers = std::vector<bartcore::SamplerBase*>{},
     columns = std::vector<size_t>{}]() mutable -> SEXP {
    size_t numSamplers = static_cast<size_t>(Rf_xlength(ptrsExpr));
    if (numSamplers == 0 ||
        static_cast<size_t>(Rf_xlength(columnsExpr)) != numSamplers)
      Rf_error("bartcore_updatePredictorPerObservationJointly requires one "
               "column per sampler");

    samplers.resize(numSamplers);
    columns.resize(numSamplers);
    for (size_t k = 0; k < numSamplers; ++k) {
      BartcoreHolder& holder(
        holderFromExpression(VECTOR_ELT(ptrsExpr, static_cast<R_xlen_t>(k))));
      refusePredictorMutation(
        *holder.sampler, "bartcore_updatePredictorPerObservationJointly");
      refuseMultiForestTransactionalUpdate(
        *holder.sampler, "bartcore_updatePredictorPerObservationJointly");
      samplers[k] = holder.sampler.get();
      int column = INTEGER(columnsExpr)[k];
      if (column < 1 ||
          static_cast<size_t>(column) > samplers[k]->numPredictors())
        Rf_error("bartcore_updatePredictorPerObservationJointly column out of "
                 "range");
      columns[k] = static_cast<size_t>(column - 1);
    }

    size_t numObservations = samplers[0]->numObservations();
    for (size_t k = 1; k < numSamplers; ++k)
      if (samplers[k]->numObservations() != numObservations)
        Rf_error("bartcore_updatePredictorPerObservationJointly requires "
                 "index-aligned samplers");
    if (static_cast<size_t>(Rf_xlength(xExpr)) != numObservations)
      Rf_error("bartcore_updatePredictorPerObservationJointly requires one "
               "value per observation");
    for (size_t k = 0; k < numSamplers; ++k)
      validateColumnValues(samplers[k]->data(), columns[k], REAL(xExpr),
                           numObservations);

    std::unique_ptr<bool[]> installed(new bool[numObservations]);

    GetRNGstate();  // scan-order permutation
    bool treesAreValid = bartcore::updatePredictorPerObservationJointly(
      samplers.data(), numSamplers, REAL(xExpr), columns.data(),
      installed.get());
    PutRNGstate();

    if (!treesAreValid) {
      installed.reset();  // free before longjmp
      Rf_error("bartcore updatePredictorPerObservationJointly produced a tree "
               "with an empty leaf");
    }

    SEXP result = PROTECT(
      Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(numObservations)));
    for (size_t i = 0; i < numObservations; ++i)
      LOGICAL(result)[i] = installed[i] ? TRUE : FALSE;
    UNPROTECT(1);
    return result;
  });
}

// State serialization. The returned object is engine-specific and opaque:
// one list per chain (flattened trees, the saved-tree buffer, accumulated
// fits, per-tree observation orderings, internal-scale sigma, k, the
// response transform, latents, dart state, and the serialized rng), with
// the store's cut points and the saved-tree write position as attributes.
// A sampler restored from it over the same data continues bitwise
// identically to one that was never stored; the stored transform makes
// that hold even when setOffset(updateScale) moved the scale after
// creation.

// flattened trees as parallel R vectors: concatenated 1-based-with-(-1)
// variables, tagged 8-byte payloads, per-tree node counts, and flags (the
// missing direction plus the payload tag). The payload is a raw word rather
// than a REALSXP because an inline categorical mask is a uint64 whose bits a
// numeric might normalize; leaf and ordinal payloads are the doubles they were
static void storeFlatTrees(SEXP chainExpr, int variablesSlot, int valuesSlot,
                           int sizesSlot, int flagsSlot,
                           const std::vector<std::vector<bartcore::FlatNode>>& trees) {
  R_xlen_t totalNumNodes = 0;
  for (const std::vector<bartcore::FlatNode>& tree : trees)
    totalNumNodes += static_cast<R_xlen_t>(tree.size());

  SET_VECTOR_ELT(chainExpr, variablesSlot,
                 Rf_allocVector(INTSXP, totalNumNodes));
  SET_VECTOR_ELT(chainExpr, valuesSlot,
                 Rf_allocVector(RAWSXP, totalNumNodes * 8));
  SET_VECTOR_ELT(chainExpr, sizesSlot,
                 Rf_allocVector(INTSXP, static_cast<R_xlen_t>(trees.size())));
  SET_VECTOR_ELT(chainExpr, flagsSlot, Rf_allocVector(RAWSXP, totalNumNodes));

  int* variables = INTEGER(VECTOR_ELT(chainExpr, variablesSlot));
  Rbyte* values = RAW(VECTOR_ELT(chainExpr, valuesSlot));
  int* sizes = INTEGER(VECTOR_ELT(chainExpr, sizesSlot));
  Rbyte* flags = RAW(VECTOR_ELT(chainExpr, flagsSlot));
  R_xlen_t offset = 0;
  for (size_t t = 0; t < trees.size(); ++t) {
    sizes[t] = static_cast<int>(trees[t].size());
    for (const bartcore::FlatNode& node : trees[t]) {
      variables[offset] = node.variable >= 0 ? node.variable + 1
                                             : node.variable;
      std::memcpy(values + offset * 8, &node.value, 8);
      flags[offset] = node.flags;
      ++offset;
    }
  }
}

// linear-leaf slope arrays as one concatenated R vector; each tree's block
// holds numSlopes doubles per leaf in pre-order, splittable by the tree
// sizes ((m + 1) / 2 leaves for m records)
static void storeTreeParams(SEXP chainExpr, int paramsSlot,
                            const std::vector<std::vector<double>>& params) {
  R_xlen_t totalNumParams = 0;
  for (const std::vector<double>& treeParams : params)
    totalNumParams += static_cast<R_xlen_t>(treeParams.size());

  SET_VECTOR_ELT(chainExpr, paramsSlot,
                 Rf_allocVector(REALSXP, totalNumParams));
  double* out = REAL(VECTOR_ELT(chainExpr, paramsSlot));
  for (const std::vector<double>& treeParams : params) {
    std::memcpy(out, treeParams.data(), treeParams.size() * sizeof(double));
    out += treeParams.size();
  }
}

// pooled categorical mask channels as one raw vector, each word written
// explicitly little-endian (raw bytes serialize portably where NaN payloads
// in a numeric might not), concatenated across trees; splittable by walking
// each tree's rules against the store's category counts
static void storeTreeMasks(
    SEXP chainExpr, int masksSlot,
    const std::vector<std::vector<std::uint64_t>>& masks) {
  R_xlen_t totalNumWords = 0;
  for (const std::vector<std::uint64_t>& treeMasks : masks)
    totalNumWords += static_cast<R_xlen_t>(treeMasks.size());

  SET_VECTOR_ELT(chainExpr, masksSlot,
                 Rf_allocVector(RAWSXP, totalNumWords * 8));
  Rbyte* out = RAW(VECTOR_ELT(chainExpr, masksSlot));
  for (const std::vector<std::uint64_t>& treeMasks : masks)
    for (std::uint64_t word : treeMasks)
      for (int b = 0; b < 8; ++b)
        *out++ = static_cast<Rbyte>((word >> (8 * b)) & 0xFFu);
}

// the number of side-channel words a flattened tree's rules occupy
static size_t maskWordsForFlatTree(
    const std::vector<bartcore::FlatNode>& tree,
    const bartcore::ColumnStore& store) {
  size_t numWords = 0;
  for (const bartcore::FlatNode& node : tree) {
    if (node.variable < 0) continue;
    size_t j = static_cast<size_t>(node.variable);
    if (j < store.numPredictors && store.columnIsPooled(j))
      numWords += bartcore::maskWordsForCount(store.numCuts[j]);
  }
  return numWords;
}

static bool readTreeMasks(
    SEXP masksExpr, const std::vector<std::vector<bartcore::FlatNode>>& trees,
    const bartcore::ColumnStore& store,
    std::vector<std::vector<std::uint64_t>>& masks,
    const char** errorMessage) {
  R_xlen_t totalNumWords = 0;
  for (const std::vector<bartcore::FlatNode>& tree : trees)
    totalNumWords += static_cast<R_xlen_t>(maskWordsForFlatTree(tree, store));

  if (Rf_isNull(masksExpr) || TYPEOF(masksExpr) != RAWSXP ||
      Rf_xlength(masksExpr) != totalNumWords * 8) {
    *errorMessage = "malformed category masks in bartcore state";
    return false;
  }

  const Rbyte* in = RAW(masksExpr);
  masks.resize(trees.size());
  for (size_t t = 0; t < trees.size(); ++t) {
    size_t numWords = maskWordsForFlatTree(trees[t], store);
    masks[t].resize(numWords);
    for (size_t w = 0; w < numWords; ++w) {
      std::uint64_t word = 0;
      for (int b = 0; b < 8; ++b)
        word |= static_cast<std::uint64_t>(*in++) << (8 * b);
      masks[t][w] = word;
    }
  }
  return true;
}

static bool readTreeParams(
    SEXP paramsExpr,
    const std::vector<std::vector<bartcore::FlatNode>>& trees,
    size_t numSlopes, std::vector<std::vector<double>>& params,
    const char** errorMessage) {
  R_xlen_t totalNumParams = 0;
  for (const std::vector<bartcore::FlatNode>& tree : trees)
    totalNumParams +=
      static_cast<R_xlen_t>(((tree.size() + 1) / 2) * numSlopes);

  if (Rf_isNull(paramsExpr) || !Rf_isReal(paramsExpr) ||
      Rf_xlength(paramsExpr) != totalNumParams) {
    *errorMessage = "malformed leaf parameters in bartcore state";
    return false;
  }

  const double* values = REAL(paramsExpr);
  params.resize(trees.size());
  for (size_t t = 0; t < trees.size(); ++t) {
    size_t numParams = ((trees[t].size() + 1) / 2) * numSlopes;
    params[t].assign(values, values + numParams);
    values += numParams;
  }
  return true;
}

// function-valued live-tree parameters: one fits slab of numObservations
// doubles per tree, in observation order
static bool readFunctionTreeParams(
    SEXP paramsExpr, size_t numTrees, size_t numObservations,
    std::vector<std::vector<double>>& params, const char** errorMessage) {
  if (Rf_isNull(paramsExpr) || !Rf_isReal(paramsExpr) ||
      static_cast<size_t>(Rf_xlength(paramsExpr)) !=
        numTrees * numObservations) {
    *errorMessage = "malformed leaf parameters in bartcore state";
    return false;
  }
  const double* values = REAL(paramsExpr);
  params.resize(numTrees);
  for (size_t t = 0; t < numTrees; ++t) {
    params[t].assign(values, values + numObservations);
    values += numObservations;
  }
  return true;
}

// function-valued saved side channels: variable-length per-leaf blocks (see
// computeFunctionBlockOffsets), split by walking each tree's leaf count;
// the engine's stateIsValid re-validates against the rebuilt trees
static bool readFunctionSavedParams(
    SEXP paramsExpr, const std::vector<std::vector<bartcore::FlatNode>>& trees,
    size_t numCovariates, std::vector<std::vector<double>>& params,
    const char** errorMessage) {
  if (Rf_isNull(paramsExpr) || !Rf_isReal(paramsExpr)) {
    *errorMessage = "malformed leaf parameters in bartcore state";
    return false;
  }
  const double* values = REAL(paramsExpr);
  size_t length = static_cast<size_t>(Rf_xlength(paramsExpr));
  size_t cursor = 0;
  params.resize(trees.size());
  for (size_t t = 0; t < trees.size(); ++t) {
    size_t numLeaves = (trees[t].size() + 1) / 2;
    size_t start = cursor;
    for (size_t b = 0; b < numLeaves; ++b) {
      if (cursor >= length) {
        *errorMessage = "malformed leaf parameters in bartcore state";
        return false;
      }
      double count = values[cursor];
      if (!(count >= 0.0) || count != std::floor(count) || count > 1.0e8) {
        *errorMessage = "malformed leaf parameters in bartcore state";
        return false;
      }
      size_t width = count == 0.0
        ? 2
        : 1 + static_cast<size_t>(count) * (1 + numCovariates);
      if (cursor + width > length) {
        *errorMessage = "malformed leaf parameters in bartcore state";
        return false;
      }
      cursor += width;
    }
    params[t].assign(values + start, values + cursor);
  }
  if (cursor != length) {
    *errorMessage = "malformed leaf parameters in bartcore state";
    return false;
  }
  return true;
}

// the inverse; errorMessage is set instead of erroring so callers can clean
// up C++ state first. A pooled node's numMaskWords is re-derived from the
// store here so raw-value replay of a restored saved tree has its clamp
static bool readFlatTrees(SEXP variablesExpr, SEXP valuesExpr, SEXP sizesExpr,
                          SEXP flagsExpr, const bartcore::ColumnStore& store,
                          std::vector<std::vector<bartcore::FlatNode>>& trees,
                          const char** errorMessage) {
  if (!Rf_isInteger(variablesExpr) || TYPEOF(valuesExpr) != RAWSXP ||
      !Rf_isInteger(sizesExpr) ||
      Rf_xlength(valuesExpr) != Rf_xlength(variablesExpr) * 8 ||
      (!Rf_isNull(flagsExpr) &&
       (TYPEOF(flagsExpr) != RAWSXP ||
        Rf_xlength(flagsExpr) != Rf_xlength(variablesExpr)))) {
    *errorMessage = "malformed trees in bartcore state";
    return false;
  }
  const Rbyte* flags = Rf_isNull(flagsExpr) ? NULL : RAW(flagsExpr);
  R_xlen_t numTrees = Rf_xlength(sizesExpr);
  const int* variables = INTEGER(variablesExpr);
  const Rbyte* values = RAW(valuesExpr);
  const int* sizes = INTEGER(sizesExpr);

  R_xlen_t totalNumNodes = 0;
  for (R_xlen_t t = 0; t < numTrees; ++t) {
    if (sizes[t] < 1) {
      *errorMessage = "malformed trees in bartcore state";
      return false;
    }
    totalNumNodes += sizes[t];
  }
  if (totalNumNodes != Rf_xlength(variablesExpr)) {
    *errorMessage = "malformed trees in bartcore state";
    return false;
  }

  trees.resize(static_cast<size_t>(numTrees));
  R_xlen_t offset = 0;
  for (R_xlen_t t = 0; t < numTrees; ++t) {
    trees[static_cast<size_t>(t)].resize(static_cast<size_t>(sizes[t]));
    for (int k = 0; k < sizes[t]; ++k) {
      bartcore::FlatNode& node(trees[static_cast<size_t>(t)]
                                    [static_cast<size_t>(k)]);
      int variable = variables[offset];
      node.variable = variable >= 1 ? variable - 1 : bartcore::invalidVariable;
      std::memcpy(&node.value, values + offset * 8, 8);
      node.flags = flags != NULL ? static_cast<std::uint8_t>(flags[offset]) : 0;
      if (node.variable >= 0 &&
          static_cast<size_t>(node.variable) < store.numPredictors &&
          store.columnIsPooled(static_cast<size_t>(node.variable)))
        node.numMaskWords = static_cast<std::uint32_t>(
          bartcore::maskWordsForCount(
            store.numCuts[static_cast<size_t>(node.variable)]));
      ++offset;
    }
  }
  return true;
}

SEXP bartcore_storeState(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  return bartcore_bridge::storeState(*holder.sampler);
}

SEXP bartcore_setState(SEXP ptrExpr, SEXP stateExpr,
                       SEXP currentPredictorsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  // restoring cut points re-quantizes from raw values, which views lack
  refuseRequantizeWithoutSource(*holder.sampler, "bartcore_setState");
  // a cross-grid restore re-quantizes dense columns from the supplied data@x;
  // a same-spec continuation skips per column, so a null source is harmless
  const double* currentPredictors =
    Rf_isReal(currentPredictorsExpr) ? REAL(currentPredictorsExpr) : NULL;
  bartcore_bridge::setState(*holder.sampler, stateExpr, currentPredictors);
  return R_NilValue;
}

SEXP bartcore_installForests(SEXP ptrExpr, SEXP donorStateExpr,
                             SEXP samplesExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refuseRequantizeWithoutSource(*holder.sampler, "bartcore_installForests");
  bartcore_bridge::installForests(*holder.sampler, donorStateExpr, samplesExpr);
  return R_NilValue;
}

// Fits for new data on the original response scale (binary responses give
// the latent scale, as the classic engine does). With keepTrees the saved
// trees produce numTestObservations x numSamples (x numChains) fits; without,
// the live trees produce a single set per chain. A multi-location combiner
// (multinomial: K softmax channels) inserts the K dimension between the rows
// and the samples, matching the run's test channel. offset, when non-null, is
// added to every sample's fits (refused for a multi-location surface, whose
// output is probabilities rather than an additive latent scale).
SEXP bartcore_predict(SEXP ptrExpr, SEXP xTestExpr, SEXP offsetExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  // predict() sums only the prognostic forest, so a BCF prediction would drop
  // the treatment forest and the glue; refuse it for the same reason recorded
  // BCF test fits are undefined.
  refuseBCFTestSurface(sampler, "bartcore_predict");

  size_t numTestObservations =
    validatePredictorMatrix(sampler, xTestExpr, "bartcore_predict");
  if (numTestObservations == 0) Rf_error("bartcore_predict requires rows");

  const double* offset = NULL;
  if (!Rf_isNull(offsetExpr)) {
    if (!Rf_isReal(offsetExpr) ||
        static_cast<size_t>(Rf_xlength(offsetExpr)) != numTestObservations)
      Rf_error("bartcore_predict offset must have one value per row");
    offset = REAL(offsetExpr);
  }

  size_t capacity = sampler.savedTreeCapacity();
  size_t numChains = sampler.numChains();
  size_t numSamples = capacity > 0 ? capacity : 1;
  size_t numLocations = sampler.numReportedLocations();

  // a multi-location (multinomial softmax) surface reports probabilities, not
  // an additive latent scale, so an offset has no meaning there
  if (offset != NULL && numLocations > 1)
    Rf_error("bartcore_predict: an offset is undefined for a multi-location "
             "(multinomial softmax) predict surface");

  int numTestInt = static_cast<int>(numTestObservations);
  SEXP resultExpr;
  if (numLocations > 1) {
    // numTestObservations x numLocations x numSamples (x numChains): the K
    // location dimension inserts between the rows and the slots, so the shape
    // and byte layout match the run's test channel exactly (allocFitsArray)
    int dims[4];
    int numDims = 0;
    dims[numDims++] = numTestInt;
    dims[numDims++] = static_cast<int>(numLocations);
    dims[numDims++] = static_cast<int>(numSamples);
    if (numChains > 1) dims[numDims++] = static_cast<int>(numChains);
    R_xlen_t total = 1;
    for (int d = 0; d < numDims; ++d) total *= dims[d];
    resultExpr = PROTECT(Rf_allocVector(REALSXP, total));
    SEXP dimExpr = Rf_allocVector(INTSXP, numDims);
    for (int d = 0; d < numDims; ++d) INTEGER(dimExpr)[d] = dims[d];
    Rf_setAttrib(resultExpr, R_DimSymbol, dimExpr);
  } else if (capacity > 0) {
    resultExpr = numChains == 1
      ? PROTECT(Rf_allocMatrix(REALSXP, numTestInt,
                               static_cast<int>(capacity)))
      : PROTECT(Rf_alloc3DArray(REALSXP, numTestInt,
                                static_cast<int>(capacity),
                                static_cast<int>(numChains)));
  } else {
    resultExpr = numChains == 1
      ? PROTECT(Rf_allocVector(REALSXP,
                               static_cast<R_xlen_t>(numTestObservations)))
      : PROTECT(Rf_allocMatrix(REALSXP, numTestInt,
                               static_cast<int>(numChains)));
  }

  sampler.predict(REAL(xTestExpr), numTestObservations, REAL(resultExpr));

  if (offset != NULL) {
    for (size_t slab = 0; slab < numSamples * numChains; ++slab)
      misc_addVectorsInPlace(offset, numTestObservations,
                             REAL(resultExpr) + slab * numTestObservations);
  }

  UNPROTECT(1);
  return resultExpr;
}

// index conversion around bartcore_bridge::getTrees, which describes the
// data.frame produced
SEXP bartcore_getTrees(SEXP ptrExpr, SEXP chainNumsExpr, SEXP sampleNumsExpr,
                       SEXP treeNumsExpr, SEXP currentExpr, SEXP newdataExpr,
                       SEXP trainingDataExpr, SEXP forestExpr) {
  return unwindProtect([&, chainIndices = std::vector<size_t>{},
                        sampleIndices = std::vector<size_t>{},
                        treeIndices = std::vector<size_t>{}]() mutable -> SEXP {
    BartcoreHolder& holder(holderFromExpression(ptrExpr));
    bartcore::SamplerBase& sampler(*holder.sampler);

    // forest addressing follows getForestFits: 0-based, unconverted
    size_t forestIndex = static_cast<size_t>(Rf_asInteger(forestExpr));
    if (forestIndex >= sampler.numForests())
      Rf_error("bartcore_getTrees forest index out of range");

    bool useLiveTrees = Rf_asLogical(currentExpr) == TRUE;
    bool useSaved = sampler.savedTreeCapacity() > 0 && !useLiveTrees;

    chainIndices.resize(static_cast<size_t>(Rf_xlength(chainNumsExpr)));
    for (size_t i = 0; i < chainIndices.size(); ++i) {
      int chainNum = INTEGER(chainNumsExpr)[i];
      if (chainNum < 1) Rf_error("bartcore_getTrees chain number out of range");
      chainIndices[i] = static_cast<size_t>(chainNum - 1);
    }
    if (useSaved) {
      sampleIndices.resize(static_cast<size_t>(Rf_xlength(sampleNumsExpr)));
      for (size_t i = 0; i < sampleIndices.size(); ++i) {
        int sampleNum = INTEGER(sampleNumsExpr)[i];
        if (sampleNum < 1)
          Rf_error("bartcore_getTrees sample number out of range");
        sampleIndices[i] = static_cast<size_t>(sampleNum - 1);
      }
    }
    treeIndices.resize(static_cast<size_t>(Rf_xlength(treeNumsExpr)));
    for (size_t i = 0; i < treeIndices.size(); ++i) {
      int treeNum = INTEGER(treeNumsExpr)[i];
      if (treeNum < 1) Rf_error("bartcore_getTrees tree number out of range");
      treeIndices[i] = static_cast<size_t>(treeNum - 1);
    }

    const double* newdata = NULL;
    size_t newdataNumRows = 0;
    if (!Rf_isNull(newdataExpr)) {
      newdataNumRows =
        validatePredictorMatrix(sampler, newdataExpr, "bartcore_getTrees");
      newdata = REAL(newdataExpr);
    }

    // saved-tree replay reads the current training predictors the R method
    // supplies (data@x); the engine keeps no matrix. A dense matrix serves; a
    // sparse/absent source leaves saved-tree counts unpopulated (a view sampler)
    const double* trainingReplay = NULL;
    size_t trainingReplayNumRows = 0;
    if (useSaved && newdata == NULL && Rf_isReal(trainingDataExpr)) {
      SEXP dims = Rf_getAttrib(trainingDataExpr, R_DimSymbol);
      if (!Rf_isNull(dims) && Rf_xlength(dims) == 2 &&
          static_cast<size_t>(INTEGER(dims)[1]) == sampler.numPredictors()) {
        trainingReplay = REAL(trainingDataExpr);
        trainingReplayNumRows = static_cast<size_t>(INTEGER(dims)[0]);
      }
    }

    return bartcore_bridge::getTrees(
      sampler, chainIndices.data(), chainIndices.size(), sampleIndices.data(),
      sampleIndices.size(), treeIndices.data(), treeIndices.size(),
      useLiveTrees, newdata, newdataNumRows, trainingReplay,
      trainingReplayNumRows, forestIndex, "bartcore_getTrees");
  });
}

SEXP bartcore_printTrees(SEXP ptrExpr, SEXP chainNumsExpr, SEXP sampleNumsExpr,
                         SEXP treeNumsExpr) {
  return unwindProtect([&, chainIndices = std::vector<size_t>{},
                        sampleIndices = std::vector<size_t>{},
                        treeIndices = std::vector<size_t>{}]() mutable -> SEXP {
    BartcoreHolder& holder(holderFromExpression(ptrExpr));
    bartcore::SamplerBase& sampler(*holder.sampler);

    size_t capacity = sampler.savedTreeCapacity();

    if (Rf_isNull(chainNumsExpr)) {
      for (size_t i = 0; i < sampler.numChains(); ++i) chainIndices.push_back(i);
    } else {
      for (R_xlen_t i = 0; i < Rf_xlength(chainNumsExpr); ++i) {
        int chainNum = INTEGER(chainNumsExpr)[i];
        if (chainNum < 1 || static_cast<size_t>(chainNum) > sampler.numChains())
          Rf_error("bartcore_printTrees chain number out of range");
        chainIndices.push_back(static_cast<size_t>(chainNum - 1));
      }
    }
    if (capacity > 0) {
      if (Rf_isNull(sampleNumsExpr)) {
        for (size_t i = 0; i < capacity; ++i) sampleIndices.push_back(i);
      } else {
        for (R_xlen_t i = 0; i < Rf_xlength(sampleNumsExpr); ++i) {
          int sampleNum = INTEGER(sampleNumsExpr)[i];
          if (sampleNum < 1 || static_cast<size_t>(sampleNum) > capacity)
            Rf_error("bartcore_printTrees sample number out of range");
          sampleIndices.push_back(static_cast<size_t>(sampleNum - 1));
        }
      }
    }
    if (Rf_isNull(treeNumsExpr)) {
      for (size_t i = 0; i < sampler.numTrees(); ++i) treeIndices.push_back(i);
    } else {
      for (R_xlen_t i = 0; i < Rf_xlength(treeNumsExpr); ++i) {
        int treeNum = INTEGER(treeNumsExpr)[i];
        if (treeNum < 1 || static_cast<size_t>(treeNum) > sampler.numTrees())
          Rf_error("bartcore_printTrees tree number out of range");
        treeIndices.push_back(static_cast<size_t>(treeNum - 1));
      }
    }

    sampler.printTrees(chainIndices.data(), chainIndices.size(),
                       sampleIndices.data(), sampleIndices.size(),
                       treeIndices.data(), treeIndices.size());

    return R_NilValue;
  });
}

// resultExpr, when non-null, is a preallocated numeric filled in place (the
// classic engine's storeLatents contract, which rbart_vi relies on).
SEXP bartcore_getLatents(SEXP ptrExpr, SEXP resultExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (holder.sampler->latents(0) == NULL) return R_NilValue;

  size_t numObservations = holder.sampler->numObservations();
  size_t numChains = holder.sampler->numChains();

  SEXP result;
  if (!Rf_isNull(resultExpr)) {
    if (!Rf_isReal(resultExpr) ||
        static_cast<size_t>(Rf_xlength(resultExpr)) !=
          numObservations * numChains)
      Rf_error("preallocated latents must be a numeric of length equal to "
               "the number of observations times the number of chains");
    result = PROTECT(resultExpr);
  } else if (numChains == 1) {
    result =
      PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numObservations)));
  } else {
    result = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numObservations),
                                    static_cast<int>(numChains)));
  }
  for (size_t c = 0; c < numChains; ++c)
    std::memcpy(REAL(result) + c * numObservations,
                holder.sampler->latents(c), numObservations * sizeof(double));
  UNPROTECT(1);
  return result;
}

} // extern "C"

namespace bartcore_bridge {

// Provenance stamp for the on-disk layout storeState/setState exchange.
// Version 2 tags each flat node and stores its payload as a raw word - an
// inline categorical mask no longer bit-casts through a double - so the side
// channel holds only pooled masks (past 63 categories). It also drops the
// accumulation-history slots (total.fits, indices) and the variance prior's
// internal scale (the third fit.scale element): restore rebuilds them from the
// trees, and sigma rides on the original response scale. Version 3 gives each
// chain's tree channels a forest dimension (docs/design/bcf.md): the
// tree/saved/param/mask/k slots move into a per-chain "forests" list, and BCF
// adds a "bcf" glue slot (a, aVariance, b0, b1). A single-forest state is a
// length-1 forest list.
//
// Registry rule for evolving the format (docs/design/public-surface.md 2):
// block names are APPEND-ONLY and a shipped name's on-disk encoding is FROZEN.
// A new capability adds a NEW optional block name; setState reads blocks by
// name (getListElement), defaults an absent OPTIONAL block, and refuses -
// naming the block - only when a REQUIRED (or config-conditionally-required)
// block is missing. So an ADDITIVE block addition does NOT bump the version
// (an old reader ignores the unknown name; a new reader defaults it), and MUST
// NOT bump minReadableStateFormatVersion. Only a non-additive change to an
// existing block's encoding - one that cannot be expressed as a new name -
// bumps both.
static const int stateFormatVersion = 1;

// The oldest ENCODING this reader still understands: additive block additions
// leave it here; only a non-additive encoding change raises it. Currently 1:
// the 1.0-0 encoding is the FIRST shipped format, so the pre-release
// development increments (the forests-list/BCF restructure included) are
// collapsed into it - no released reader or writer ever saw them. Pre-1.0
// states are not a compat target and cannot structurally reach the by-name
// reader (they lack the "forests" block); a state with no version attribute
// reads as 0 and is refused at the floor.
static const int minReadableStateFormatVersion = 1;

SEXP storeState(bartcore::SamplerBase& sampler) {
  bartcore::SamplerStateData state;
  sampler.getState(state);

  size_t numChains = state.chains.size();
  size_t numObservations = sampler.numObservations();

  // per-forest tree channels (a length-1 list off BCF); the k rides here too
  enum {
    FSLOT_TREE_VARS = 0, FSLOT_TREE_VALUES, FSLOT_TREE_SIZES, FSLOT_TREE_FLAGS,
    FSLOT_TREE_PARAMS, FSLOT_TREE_MASKS,
    FSLOT_SAVED_VARS, FSLOT_SAVED_VALUES, FSLOT_SAVED_SIZES, FSLOT_SAVED_FLAGS,
    FSLOT_SAVED_PARAMS, FSLOT_SAVED_MASKS,
    FSLOT_K, FSLOT_COUNT
  };
  static const char* forestSlotNames[FSLOT_COUNT] = {
    "tree.vars", "tree.values", "tree.sizes", "tree.flags", "tree.params",
    "tree.masks",
    "saved.vars", "saved.values", "saved.sizes", "saved.flags",
    "saved.params", "saved.masks",
    "k"
  };

  enum {
    SLOT_FORESTS = 0, SLOT_SIGMA, SLOT_FIT_SCALE, SLOT_LATENTS,
    SLOT_RANEF, SLOT_TAU,
    SLOT_DART_PROBABILITIES, SLOT_DART_ALPHA, SLOT_DART_UPDATES_SKIPPED,
    SLOT_RNG_STATE, SLOT_BCF, SLOT_RESID_DF, SLOT_CUTPOINTS, SLOT_DISPERSION,
    SLOT_COUNT
  };
  static const char* slotNames[SLOT_COUNT] = {
    "forests", "sigma", "fit.scale",
    "latents", "ranef", "tau",
    "dart.probabilities", "dart.alpha", "dart.updates.skipped",
    "rng.state", "bcf", "resid.df", "cutpoints", "dispersion"
  };

  SEXP resultExpr =
    PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(numChains)));
  SEXP slotNamesExpr = PROTECT(Rf_allocVector(STRSXP, SLOT_COUNT));
  for (int slot = 0; slot < SLOT_COUNT; ++slot)
    SET_STRING_ELT(slotNamesExpr, slot, Rf_mkChar(slotNames[slot]));
  SEXP forestSlotNamesExpr = PROTECT(Rf_allocVector(STRSXP, FSLOT_COUNT));
  for (int slot = 0; slot < FSLOT_COUNT; ++slot)
    SET_STRING_ELT(forestSlotNamesExpr, slot, Rf_mkChar(forestSlotNames[slot]));

  for (size_t c = 0; c < numChains; ++c) {
    const bartcore::ChainStateData& chainState(state.chains[c]);
    SEXP chainExpr = PROTECT(Rf_allocVector(VECSXP, SLOT_COUNT));
    Rf_setAttrib(chainExpr, R_NamesSymbol, slotNamesExpr);

    size_t numForests = chainState.forests.size();
    SEXP forestsExpr =
      PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(numForests)));
    for (size_t f = 0; f < numForests; ++f) {
      const bartcore::ForestStateData& fs(chainState.forests[f]);
      SEXP forestExpr = PROTECT(Rf_allocVector(VECSXP, FSLOT_COUNT));
      Rf_setAttrib(forestExpr, R_NamesSymbol, forestSlotNamesExpr);

      storeFlatTrees(forestExpr, FSLOT_TREE_VARS, FSLOT_TREE_VALUES,
                     FSLOT_TREE_SIZES, FSLOT_TREE_FLAGS, fs.trees);
      if (!fs.treeParams.empty())
        storeTreeParams(forestExpr, FSLOT_TREE_PARAMS, fs.treeParams);
      if (!fs.treeMasks.empty())
        storeTreeMasks(forestExpr, FSLOT_TREE_MASKS, fs.treeMasks);
      if (!fs.savedTrees.empty()) {
        storeFlatTrees(forestExpr, FSLOT_SAVED_VARS, FSLOT_SAVED_VALUES,
                       FSLOT_SAVED_SIZES, FSLOT_SAVED_FLAGS, fs.savedTrees);
        if (!fs.savedTreeParams.empty())
          storeTreeParams(forestExpr, FSLOT_SAVED_PARAMS, fs.savedTreeParams);
        if (!fs.savedTreeMasks.empty())
          storeTreeMasks(forestExpr, FSLOT_SAVED_MASKS, fs.savedTreeMasks);
      }
      SET_VECTOR_ELT(forestExpr, FSLOT_K, Rf_ScalarReal(fs.k));
      SET_VECTOR_ELT(forestsExpr, static_cast<R_xlen_t>(f), forestExpr);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(chainExpr, SLOT_FORESTS, forestsExpr);
    UNPROTECT(1);

    SET_VECTOR_ELT(chainExpr, SLOT_SIGMA, Rf_ScalarReal(chainState.sigma));

    SET_VECTOR_ELT(chainExpr, SLOT_FIT_SCALE, Rf_allocVector(REALSXP, 2));
    REAL(VECTOR_ELT(chainExpr, SLOT_FIT_SCALE))[0] = chainState.fitMin;
    REAL(VECTOR_ELT(chainExpr, SLOT_FIT_SCALE))[1] = chainState.fitMax;

    if (!chainState.latents.empty()) {
      SET_VECTOR_ELT(chainExpr, SLOT_LATENTS,
                     Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                               numObservations)));
      std::memcpy(REAL(VECTOR_ELT(chainExpr, SLOT_LATENTS)),
                  chainState.latents.data(),
                  numObservations * sizeof(double));
    }

    if (!chainState.groupEffects.empty()) {
      SET_VECTOR_ELT(chainExpr, SLOT_RANEF,
                     Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                      chainState.groupEffects.size())));
      std::memcpy(REAL(VECTOR_ELT(chainExpr, SLOT_RANEF)),
                  chainState.groupEffects.data(),
                  chainState.groupEffects.size() * sizeof(double));
      SET_VECTOR_ELT(chainExpr, SLOT_TAU,
                     Rf_ScalarReal(chainState.groupTau));
    }

    if (!chainState.dartProbabilities.empty()) {
      SET_VECTOR_ELT(chainExpr, SLOT_DART_PROBABILITIES,
                     Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                      chainState.dartProbabilities.size())));
      std::memcpy(REAL(VECTOR_ELT(chainExpr, SLOT_DART_PROBABILITIES)),
                  chainState.dartProbabilities.data(),
                  chainState.dartProbabilities.size() * sizeof(double));
      SET_VECTOR_ELT(chainExpr, SLOT_DART_ALPHA,
                     Rf_ScalarReal(chainState.dartAlpha));
      SET_VECTOR_ELT(chainExpr, SLOT_DART_UPDATES_SKIPPED,
                     Rf_ScalarInteger(static_cast<int>(
                       chainState.dartNumUpdatesSkipped)));
    }

    SET_VECTOR_ELT(chainExpr, SLOT_RNG_STATE,
                   Rf_allocVector(INTSXP, static_cast<R_xlen_t>(
                                    chainState.rngState.size() / sizeof(int))));
    std::memcpy(INTEGER(VECTOR_ELT(chainExpr, SLOT_RNG_STATE)),
                chainState.rngState.data(), chainState.rngState.size());

    if (chainState.hasBCF) {
      SET_VECTOR_ELT(chainExpr, SLOT_BCF, Rf_allocVector(REALSXP, 4));
      double* bcf = REAL(VECTOR_ELT(chainExpr, SLOT_BCF));
      bcf[0] = chainState.a;
      bcf[1] = chainState.aVariance;
      bcf[2] = chainState.b0;
      bcf[3] = chainState.b1;
    }

    // the t-only residual df; lambda already rode the latents slot above. A
    // gaussian (or any non-t) chain carries NaN and writes no block.
    if (std::isfinite(chainState.residualDf))
      SET_VECTOR_ELT(chainExpr, SLOT_RESID_DF,
                     Rf_ScalarReal(chainState.residualDf));

    // the ordinal-only length-(K-1) cutpoint vector; z already rode the latents
    // slot above. A non-ordinal chain leaves it empty and writes no block, so
    // old and other-family states omit the slot.
    if (!chainState.cutpoints.empty()) {
      SET_VECTOR_ELT(chainExpr, SLOT_CUTPOINTS,
                     Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                      chainState.cutpoints.size())));
      std::memcpy(REAL(VECTOR_ELT(chainExpr, SLOT_CUTPOINTS)),
                  chainState.cutpoints.data(),
                  chainState.cutpoints.size() * sizeof(double));
    }

    // the nbinom-only dispersion r; omega already rode the latents slot above. A
    // non-count chain carries NaN and writes no block, so old and other-family
    // states omit the slot.
    if (std::isfinite(chainState.dispersion))
      SET_VECTOR_ELT(chainExpr, SLOT_DISPERSION,
                     Rf_ScalarReal(chainState.dispersion));

    SET_VECTOR_ELT(resultExpr, static_cast<R_xlen_t>(c), chainExpr);
    UNPROTECT(1);
  }

  SEXP cutPointsExpr = PROTECT(Rf_allocVector(
    VECSXP, static_cast<R_xlen_t>(state.cutPoints.size())));
  for (size_t j = 0; j < state.cutPoints.size(); ++j) {
    if (state.cutPoints[j].empty()) continue;  // categorical column
    SET_VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(j),
                   Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                             state.cutPoints[j].size())));
    std::memcpy(REAL(VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(j))),
                state.cutPoints[j].data(),
                state.cutPoints[j].size() * sizeof(double));
  }
  Rf_setAttrib(resultExpr, Rf_install("cutPoints"), cutPointsExpr);
  Rf_setAttrib(resultExpr, Rf_install("currentSampleNum"),
               Rf_ScalarInteger(static_cast<int>(state.currentSampleNum)));
  Rf_setAttrib(resultExpr, Rf_install("formatVersion"),
               Rf_ScalarInteger(stateFormatVersion));
  Rf_setAttrib(resultExpr, Rf_install("packageVersion"),
               Rf_mkString(PACKAGE_VERSION));
  SEXP classExpr = PROTECT(Rf_mkString("bartcoreState"));
  Rf_setAttrib(resultExpr, R_ClassSymbol, classExpr);

  UNPROTECT(5);
  return resultExpr;
}

void setState(bartcore::SamplerBase& sampler, SEXP stateExpr,
              const double* currentPredictors) {
  if (!Rf_inherits(stateExpr, "bartcoreState"))
    Rf_error("'state' must be a bartcore state object");

  SEXP formatVersionExpr =
    PROTECT(Rf_getAttrib(stateExpr, Rf_install("formatVersion")));
  int formatVersion = Rf_isInteger(formatVersionExpr) &&
      Rf_xlength(formatVersionExpr) == 1 ? INTEGER(formatVersionExpr)[0] : 0;
  UNPROTECT(1);
  // An encoding floor, not an equality gate: a newer release's additive blocks
  // are read by name and defaulted when absent, so a state written by any
  // release at or past the floor loads; only a genuinely older, incompatible
  // encoding is refused (registry rule at minReadableStateFormatVersion).
  if (formatVersion < minReadableStateFormatVersion) {
    SEXP packageVersionExpr =
      Rf_getAttrib(stateExpr, Rf_install("packageVersion"));
    const char* packageVersion = Rf_isString(packageVersionExpr) &&
        Rf_xlength(packageVersionExpr) == 1 ?
      CHAR(STRING_ELT(packageVersionExpr, 0)) : "unknown";
    Rf_error("state encoding version %d (written by dbarts %s) predates the "
             "oldest this dbarts (%d) can read; re-fit with this version",
             formatVersion, packageVersion, minReadableStateFormatVersion);
  }

  if (static_cast<size_t>(Rf_xlength(stateExpr)) != sampler.numChains())
    Rf_error("'state' length must equal number of chains");

  // Rf_error longjmps past destructors, so parse with an error accumulator,
  // free the C++ state, and error at the end.
  const char* errorMessage = NULL;

  // Block-name refusals: a REQUIRED (or config-conditionally-required) block
  // is named on absence vs malformation. The buffer outlives every break and
  // the final Rf_error below, which longjmps past this frame.
  char blockError[96];
  auto missingBlock = [&blockError](const char* name) -> const char* {
    std::snprintf(blockError, sizeof blockError,
                  "bartcore state is missing required block '%s'", name);
    return blockError;
  };
  auto malformedBlock = [&blockError](const char* name) -> const char* {
    std::snprintf(blockError, sizeof blockError,
                  "bartcore state block '%s' is malformed", name);
    return blockError;
  };

  bartcore::SamplerStateData state;
  state.chains.resize(sampler.numChains());

  SEXP cutPointsExpr = Rf_getAttrib(stateExpr, Rf_install("cutPoints"));
  if (Rf_isNull(cutPointsExpr) ||
      static_cast<size_t>(Rf_xlength(cutPointsExpr)) !=
        sampler.numPredictors()) {
    errorMessage = "malformed cut points in bartcore state";
  } else {
    state.cutPoints.resize(sampler.numPredictors());
    for (size_t j = 0; j < sampler.numPredictors(); ++j) {
      SEXP cutsExpr = VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(j));
      if (Rf_isNull(cutsExpr)) continue;
      if (!Rf_isReal(cutsExpr)) {
        errorMessage = "malformed cut points in bartcore state";
        break;
      }
      state.cutPoints[j].assign(
        REAL(cutsExpr), REAL(cutsExpr) + Rf_xlength(cutsExpr));
    }
  }

  SEXP sampleNumExpr =
    PROTECT(Rf_getAttrib(stateExpr, Rf_install("currentSampleNum")));
  if (errorMessage == NULL) {
    if (!Rf_isInteger(sampleNumExpr) || Rf_xlength(sampleNumExpr) != 1 ||
        INTEGER(sampleNumExpr)[0] < 0)
      errorMessage = "malformed sample number in bartcore state";
    else
      state.currentSampleNum = static_cast<size_t>(INTEGER(sampleNumExpr)[0]);
  }
  UNPROTECT(1);

  for (size_t c = 0; c < sampler.numChains() && errorMessage == NULL; ++c) {
    SEXP chainExpr = VECTOR_ELT(stateExpr, static_cast<R_xlen_t>(c));
    bartcore::ChainStateData& chainState(state.chains[c]);

    SEXP forestsExpr = getListElement(chainExpr, "forests");
    if (Rf_isNull(forestsExpr)) {
      errorMessage = missingBlock("forests");
      break;
    }
    if (TYPEOF(forestsExpr) != VECSXP) {
      errorMessage = malformedBlock("forests");
      break;
    }
    size_t numForests = static_cast<size_t>(Rf_xlength(forestsExpr));
    chainState.forests.resize(numForests);
    // every forest of a chain shares the leaf shape, so these dispatch flags
    // are chain-level
    size_t numLeafCovariates = sampler.numLeafCovariates();
    for (size_t f = 0; f < numForests && errorMessage == NULL; ++f) {
      SEXP forestExpr = VECTOR_ELT(forestsExpr, static_cast<R_xlen_t>(f));
      bartcore::ForestStateData& fs(chainState.forests[f]);

      if (Rf_isNull(getListElement(forestExpr, "tree.vars"))) {
        errorMessage = missingBlock("tree.vars");
        break;
      }
      if (!readFlatTrees(getListElement(forestExpr, "tree.vars"),
                         getListElement(forestExpr, "tree.values"),
                         getListElement(forestExpr, "tree.sizes"),
                         getListElement(forestExpr, "tree.flags"),
                         sampler.data(), fs.trees, &errorMessage))
        break;
      SEXP savedSizesExpr = getListElement(forestExpr, "saved.sizes");
      if (!Rf_isNull(savedSizesExpr) &&
          !readFlatTrees(getListElement(forestExpr, "saved.vars"),
                         getListElement(forestExpr, "saved.values"),
                         savedSizesExpr,
                         getListElement(forestExpr, "saved.flags"),
                         sampler.data(), fs.savedTrees, &errorMessage))
        break;

      // linear-leaf states must carry their slope arrays; function-valued
      // states carry fits slabs and variable-length saved blocks instead
      if (sampler.usesFunctionLeaves()) {
        if (Rf_isNull(getListElement(forestExpr, "tree.params"))) {
          errorMessage = missingBlock("tree.params");
          break;
        }
        if (!readFunctionTreeParams(getListElement(forestExpr, "tree.params"),
                                    fs.trees.size(), sampler.numObservations(),
                                    fs.treeParams, &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readFunctionSavedParams(getListElement(forestExpr, "saved.params"),
                                     fs.savedTrees, numLeafCovariates,
                                     fs.savedTreeParams, &errorMessage))
          break;
      } else if (numLeafCovariates > 0) {
        if (Rf_isNull(getListElement(forestExpr, "tree.params"))) {
          errorMessage = missingBlock("tree.params");
          break;
        }
        if (!readTreeParams(getListElement(forestExpr, "tree.params"),
                            fs.trees, numLeafCovariates, fs.treeParams,
                            &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeParams(getListElement(forestExpr, "saved.params"),
                            fs.savedTrees, numLeafCovariates,
                            fs.savedTreeParams, &errorMessage))
          break;
      }

      // pooled-categorical states must carry their mask channels
      if (sampler.data().hasPooledCategorical) {
        if (Rf_isNull(getListElement(forestExpr, "tree.masks"))) {
          errorMessage = missingBlock("tree.masks");
          break;
        }
        if (!readTreeMasks(getListElement(forestExpr, "tree.masks"), fs.trees,
                           sampler.data(), fs.treeMasks, &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeMasks(getListElement(forestExpr, "saved.masks"),
                           fs.savedTrees, sampler.data(), fs.savedTreeMasks,
                           &errorMessage))
          break;
      }

      SEXP kExpr = getListElement(forestExpr, "k");
      if (Rf_isNull(kExpr)) {
        errorMessage = missingBlock("k");
        break;
      }
      if (!Rf_isReal(kExpr) || Rf_xlength(kExpr) != 1) {
        errorMessage = malformedBlock("k");
        break;
      }
      fs.k = REAL(kExpr)[0];
    }
    if (errorMessage != NULL) break;

    SEXP sigmaExpr = getListElement(chainExpr, "sigma");
    if (Rf_isNull(sigmaExpr)) {
      errorMessage = missingBlock("sigma");
      break;
    }
    if (!Rf_isReal(sigmaExpr) || Rf_xlength(sigmaExpr) != 1) {
      errorMessage = malformedBlock("sigma");
      break;
    }
    chainState.sigma = REAL(sigmaExpr)[0];

    SEXP fitScaleExpr = getListElement(chainExpr, "fit.scale");
    if (Rf_isNull(fitScaleExpr)) {
      errorMessage = missingBlock("fit.scale");
      break;
    }
    if (!Rf_isReal(fitScaleExpr) || Rf_xlength(fitScaleExpr) != 2) {
      errorMessage = malformedBlock("fit.scale");
      break;
    }
    chainState.fitMin = REAL(fitScaleExpr)[0];
    chainState.fitMax = REAL(fitScaleExpr)[1];

    SEXP latentsExpr = getListElement(chainExpr, "latents");
    if (!Rf_isNull(latentsExpr)) {
      if (!Rf_isReal(latentsExpr)) {
        errorMessage = "malformed latents in bartcore state";
        break;
      }
      chainState.latents.assign(
        REAL(latentsExpr), REAL(latentsExpr) + Rf_xlength(latentsExpr));
    }

    SEXP ranefExpr = getListElement(chainExpr, "ranef");
    if (!Rf_isNull(ranefExpr)) {
      SEXP tauExpr = getListElement(chainExpr, "tau");
      if (!Rf_isReal(ranefExpr) || !Rf_isReal(tauExpr) ||
          Rf_xlength(tauExpr) != 1) {
        errorMessage = "malformed grouped effects in bartcore state";
        break;
      }
      chainState.groupEffects.assign(
        REAL(ranefExpr), REAL(ranefExpr) + Rf_xlength(ranefExpr));
      chainState.groupTau = REAL(tauExpr)[0];
    }

    SEXP dartProbabilitiesExpr =
      getListElement(chainExpr, "dart.probabilities");
    if (!Rf_isNull(dartProbabilitiesExpr)) {
      SEXP dartAlphaExpr = getListElement(chainExpr, "dart.alpha");
      SEXP dartSkippedExpr =
        getListElement(chainExpr, "dart.updates.skipped");
      if (!Rf_isReal(dartProbabilitiesExpr) || !Rf_isReal(dartAlphaExpr) ||
          Rf_xlength(dartAlphaExpr) != 1 || !Rf_isInteger(dartSkippedExpr) ||
          Rf_xlength(dartSkippedExpr) != 1 || INTEGER(dartSkippedExpr)[0] < 0) {
        errorMessage = "malformed dart state in bartcore state";
        break;
      }
      chainState.dartProbabilities.assign(
        REAL(dartProbabilitiesExpr),
        REAL(dartProbabilitiesExpr) + Rf_xlength(dartProbabilitiesExpr));
      chainState.dartAlpha = REAL(dartAlphaExpr)[0];
      chainState.dartNumUpdatesSkipped =
        static_cast<size_t>(INTEGER(dartSkippedExpr)[0]);
    }

    SEXP rngStateExpr = getListElement(chainExpr, "rng.state");
    if (!Rf_isNull(rngStateExpr)) {
      if (!Rf_isInteger(rngStateExpr)) {
        errorMessage = "malformed rng state in bartcore state";
        break;
      }
      chainState.rngState.resize(
        static_cast<size_t>(Rf_xlength(rngStateExpr)) * sizeof(int));
      std::memcpy(chainState.rngState.data(), INTEGER(rngStateExpr),
                  chainState.rngState.size());
    }

    SEXP bcfExpr = getListElement(chainExpr, "bcf");
    if (!Rf_isNull(bcfExpr)) {
      if (!Rf_isReal(bcfExpr) || Rf_xlength(bcfExpr) != 4) {
        errorMessage = "malformed bcf glue in bartcore state";
        break;
      }
      chainState.hasBCF = true;
      chainState.a = REAL(bcfExpr)[0];
      chainState.aVariance = REAL(bcfExpr)[1];
      chainState.b0 = REAL(bcfExpr)[2];
      chainState.b1 = REAL(bcfExpr)[3];
    }

    // additive t-only block: absent (an old or gaussian state) leaves the NaN
    // default, which stateIsValid refuses only for a t sampler
    SEXP residDfExpr = getListElement(chainExpr, "resid.df");
    if (!Rf_isNull(residDfExpr)) {
      if (!Rf_isReal(residDfExpr) || Rf_xlength(residDfExpr) != 1) {
        errorMessage = malformedBlock("resid.df");
        break;
      }
      chainState.residualDf = REAL(residDfExpr)[0];
    }

    // additive ordinal-only block: absent (an old or non-ordinal state) leaves
    // the vector empty, which stateIsValid refuses only for an ordinal sampler
    SEXP cutpointsExpr = getListElement(chainExpr, "cutpoints");
    if (!Rf_isNull(cutpointsExpr)) {
      if (!Rf_isReal(cutpointsExpr)) {
        errorMessage = malformedBlock("cutpoints");
        break;
      }
      chainState.cutpoints.assign(
        REAL(cutpointsExpr), REAL(cutpointsExpr) + Rf_xlength(cutpointsExpr));
    }

    // additive nbinom-only block: absent (an old or non-count state) leaves the
    // NaN default, which stateIsValid refuses only for an NB sampler
    SEXP dispersionExpr = getListElement(chainExpr, "dispersion");
    if (!Rf_isNull(dispersionExpr)) {
      if (!Rf_isReal(dispersionExpr) || Rf_xlength(dispersionExpr) != 1) {
        errorMessage = malformedBlock("dispersion");
        break;
      }
      chainState.dispersion = REAL(dispersionExpr)[0];
    }
  }

  bool restored =
    errorMessage == NULL && sampler.setState(state, currentPredictors);
  {
    bartcore::SamplerStateData empty;
    std::swap(state, empty);  // free before a potential longjmp
  }
  if (errorMessage != NULL) Rf_error("%s", errorMessage);
  if (!restored)
    Rf_error("state is not consistent with this sampler");
}

// Parses a "bartcoreState" donor into a SamplerStateData for a warm start,
// validating flat trees against the destination sampler's data. Only the
// channels a warm start consumes are read (trees, leaf params, masks, k,
// sigma, the fit scale, DART, and BCF glue); latents, group effects, and the
// rng are left for the destination to redraw. Function-leaf donors seed from
// their live trees, so their saved channel is skipped. The donor's own chain
// count is honored (a short donor may seed many chains). Returns an error
// string, or NULL, rather than longjmping so the caller can free state first.
static const char* readWarmStartState(SEXP stateExpr,
                                      bartcore::SamplerBase& sampler,
                                      bartcore::SamplerStateData& state) {
  size_t numChains = static_cast<size_t>(Rf_xlength(stateExpr));
  if (numChains == 0) return "warm-start donor holds no chains";

  SEXP cutPointsExpr = Rf_getAttrib(stateExpr, Rf_install("cutPoints"));
  if (Rf_isNull(cutPointsExpr) ||
      static_cast<size_t>(Rf_xlength(cutPointsExpr)) != sampler.numPredictors())
    return "malformed cut points in warm-start donor";
  state.cutPoints.resize(sampler.numPredictors());
  for (size_t j = 0; j < sampler.numPredictors(); ++j) {
    SEXP cutsExpr = VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(j));
    if (Rf_isNull(cutsExpr)) continue;
    if (!Rf_isReal(cutsExpr)) return "malformed cut points in warm-start donor";
    state.cutPoints[j].assign(REAL(cutsExpr),
                              REAL(cutsExpr) + Rf_xlength(cutsExpr));
  }

  const char* errorMessage = NULL;
  state.chains.resize(numChains);
  size_t numLeafCovariates = sampler.numLeafCovariates();
  bool functionLeaves = sampler.usesFunctionLeaves();
  for (size_t c = 0; c < numChains && errorMessage == NULL; ++c) {
    SEXP chainExpr = VECTOR_ELT(stateExpr, static_cast<R_xlen_t>(c));
    bartcore::ChainStateData& chainState(state.chains[c]);

    SEXP forestsExpr = getListElement(chainExpr, "forests");
    if (Rf_isNull(forestsExpr) || TYPEOF(forestsExpr) != VECSXP) {
      errorMessage = "malformed forests in warm-start donor";
      break;
    }
    size_t numForests = static_cast<size_t>(Rf_xlength(forestsExpr));
    chainState.forests.resize(numForests);
    for (size_t f = 0; f < numForests && errorMessage == NULL; ++f) {
      SEXP forestExpr = VECTOR_ELT(forestsExpr, static_cast<R_xlen_t>(f));
      bartcore::ForestStateData& fs(chainState.forests[f]);

      if (!readFlatTrees(getListElement(forestExpr, "tree.vars"),
                         getListElement(forestExpr, "tree.values"),
                         getListElement(forestExpr, "tree.sizes"),
                         getListElement(forestExpr, "tree.flags"),
                         sampler.data(), fs.trees, &errorMessage))
        break;
      SEXP savedSizesExpr = getListElement(forestExpr, "saved.sizes");
      if (!functionLeaves && !Rf_isNull(savedSizesExpr) &&
          !readFlatTrees(getListElement(forestExpr, "saved.vars"),
                         getListElement(forestExpr, "saved.values"),
                         savedSizesExpr, getListElement(forestExpr,
                         "saved.flags"), sampler.data(), fs.savedTrees,
                         &errorMessage))
        break;

      if (functionLeaves) {
        if (!readFunctionTreeParams(getListElement(forestExpr, "tree.params"),
                                    fs.trees.size(), sampler.numObservations(),
                                    fs.treeParams, &errorMessage))
          break;
      } else if (numLeafCovariates > 0) {
        if (!readTreeParams(getListElement(forestExpr, "tree.params"),
                            fs.trees, numLeafCovariates, fs.treeParams,
                            &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeParams(getListElement(forestExpr, "saved.params"),
                            fs.savedTrees, numLeafCovariates,
                            fs.savedTreeParams, &errorMessage))
          break;
      }

      if (sampler.data().hasPooledCategorical) {
        if (!readTreeMasks(getListElement(forestExpr, "tree.masks"), fs.trees,
                           sampler.data(), fs.treeMasks, &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeMasks(getListElement(forestExpr, "saved.masks"),
                           fs.savedTrees, sampler.data(), fs.savedTreeMasks,
                           &errorMessage))
          break;
      }

      SEXP kExpr = getListElement(forestExpr, "k");
      if (!Rf_isReal(kExpr) || Rf_xlength(kExpr) != 1) {
        errorMessage = "malformed parameters in warm-start donor";
        break;
      }
      fs.k = REAL(kExpr)[0];
    }
    if (errorMessage != NULL) break;

    SEXP sigmaExpr = getListElement(chainExpr, "sigma");
    SEXP fitScaleExpr = getListElement(chainExpr, "fit.scale");
    if (!Rf_isReal(sigmaExpr) || Rf_xlength(sigmaExpr) != 1 ||
        !Rf_isReal(fitScaleExpr) || Rf_xlength(fitScaleExpr) != 2) {
      errorMessage = "malformed parameters in warm-start donor";
      break;
    }
    chainState.sigma = REAL(sigmaExpr)[0];
    chainState.fitMin = REAL(fitScaleExpr)[0];
    chainState.fitMax = REAL(fitScaleExpr)[1];

    SEXP dartProbabilitiesExpr =
      getListElement(chainExpr, "dart.probabilities");
    if (!Rf_isNull(dartProbabilitiesExpr)) {
      SEXP dartAlphaExpr = getListElement(chainExpr, "dart.alpha");
      SEXP dartSkippedExpr = getListElement(chainExpr, "dart.updates.skipped");
      if (!Rf_isReal(dartProbabilitiesExpr) || !Rf_isReal(dartAlphaExpr) ||
          Rf_xlength(dartAlphaExpr) != 1 || !Rf_isInteger(dartSkippedExpr) ||
          Rf_xlength(dartSkippedExpr) != 1 || INTEGER(dartSkippedExpr)[0] < 0) {
        errorMessage = "malformed dart state in warm-start donor";
        break;
      }
      chainState.dartProbabilities.assign(
        REAL(dartProbabilitiesExpr),
        REAL(dartProbabilitiesExpr) + Rf_xlength(dartProbabilitiesExpr));
      chainState.dartAlpha = REAL(dartAlphaExpr)[0];
      chainState.dartNumUpdatesSkipped =
        static_cast<size_t>(INTEGER(dartSkippedExpr)[0]);
    }

    SEXP bcfExpr = getListElement(chainExpr, "bcf");
    if (!Rf_isNull(bcfExpr)) {
      if (!Rf_isReal(bcfExpr) || Rf_xlength(bcfExpr) != 4) {
        errorMessage = "malformed bcf glue in warm-start donor";
        break;
      }
      chainState.hasBCF = true;
      chainState.a = REAL(bcfExpr)[0];
      chainState.aVariance = REAL(bcfExpr)[1];
      chainState.b0 = REAL(bcfExpr)[2];
      chainState.b1 = REAL(bcfExpr)[3];
    }
  }
  return errorMessage;
}

void installForests(bartcore::SamplerBase& sampler, SEXP donorStateExpr,
                    SEXP samplesExpr) {
  if (!Rf_inherits(donorStateExpr, "bartcoreState"))
    Rf_error("'warm.start' must supply a bartcore state object");

  SEXP formatVersionExpr =
    PROTECT(Rf_getAttrib(donorStateExpr, Rf_install("formatVersion")));
  int formatVersion = Rf_isInteger(formatVersionExpr) &&
      Rf_xlength(formatVersionExpr) == 1 ? INTEGER(formatVersionExpr)[0] : 0;
  UNPROTECT(1);
  if (formatVersion < minReadableStateFormatVersion)
    Rf_error("warm-start donor encoding version %d predates the oldest this "
             "dbarts (%d) can read; re-fit the donor",
             formatVersion, minReadableStateFormatVersion);

  bartcore::SamplerStateData donor;
  const char* errorMessage = readWarmStartState(donorStateExpr, sampler, donor);

  // The mapping vectors live in a scope that closes before any Rf_error, so
  // the longjmp cannot leak them.
  bartcore::WarmStartResult result = bartcore::WarmStartResult::ok;
  {
    // Donor pool of (chain, slot): every saved slot when the donor kept trees,
    // otherwise each donor chain's live forest (slot -1).
    std::vector<std::pair<size_t, int>> pool;
    if (errorMessage == NULL) {
      for (size_t dc = 0; dc < donor.chains.size(); ++dc) {
        const bartcore::ForestStateData& f0 = donor.chains[dc].forests[0];
        if (!f0.savedTrees.empty() && !f0.trees.empty()) {
          size_t capacity = f0.savedTrees.size() / f0.trees.size();
          for (size_t s = 0; s < capacity; ++s)
            pool.emplace_back(dc, static_cast<int>(s));
        } else {
          pool.emplace_back(dc, -1);
        }
      }
      if (pool.empty()) errorMessage = "warm-start donor holds no samples";
    }

    size_t numChains = sampler.numChains();
    std::vector<std::pair<size_t, int>> sampleMap;
    if (errorMessage == NULL) {
      sampleMap.resize(numChains);
      if (Rf_isNull(samplesExpr)) {
        // spread the chains across the pool, so many chains from one donor
        // draw overdispersed starts rather than the same forest
        for (size_t c = 0; c < numChains; ++c)
          sampleMap[c] = pool[(c * pool.size()) / numChains];
      } else if (!Rf_isInteger(samplesExpr) ||
                 static_cast<size_t>(Rf_xlength(samplesExpr)) != numChains) {
        errorMessage = "'samples' must be an integer vector, one per chain";
      } else {
        for (size_t c = 0; c < numChains && errorMessage == NULL; ++c) {
          int idx = INTEGER(samplesExpr)[c];
          if (idx < 1 || static_cast<size_t>(idx) > pool.size())
            errorMessage = "'samples' entries must index the donor pool";
          else
            sampleMap[c] = pool[static_cast<size_t>(idx) - 1];
        }
      }
    }

    if (errorMessage == NULL)
      result = sampler.installForests(donor, sampleMap);
  }
  {
    bartcore::SamplerStateData empty;
    std::swap(donor, empty);  // free before a potential longjmp
  }
  if (errorMessage != NULL) Rf_error("%s", errorMessage);
  switch (result) {
    case bartcore::WarmStartResult::ok:
      break;
    case bartcore::WarmStartResult::gridMismatch:
      Rf_error("warm-start donor was fit on a different cut grid (predictors, "
               "n.cuts, or useQuantiles differ); cross-grid warm starts are "
               "not supported");
    case bartcore::WarmStartResult::dartMismatch:
      Rf_error("DART state transfers only between two DART fits; the donor and "
               "destination disagree on dart");
    case bartcore::WarmStartResult::shapeMismatch:
      Rf_error("warm-start donor is not shape-compatible with this sampler "
               "(number of trees, forests, or predictors differ)");
  }
}

// Column-major gather of the requested trees, plus the flags that decide
// which data.frame columns exist.
struct GatheredTrees {
  bool includeChain, includeSample, anyCategorical, anyMissing;
  size_t numSlopes;
  std::vector<int> chain, sample, tree, count, variable, missing;
  std::vector<double> value;
  // every categorical rule reports its mask decoded, one L/R per observed
  // category (value is NA); the R wrapper pads to the declared level count
  std::vector<std::string> directions;
  // one column per leaf covariate, NA on internal nodes
  std::vector<std::vector<double>> slopes;
};

// Walks the requested (chain, sample, tree) triples, flattening live trees
// or reading saved records, and appends one row per node. replayData, when
// non-null, re-routes its rows through each tree for the n column.
void gatherTrees(bartcore::SamplerBase& sampler, const size_t* chainIndices,
                 size_t numChainIndices, const size_t* sampleIndices,
                 size_t numSampleIndices, const size_t* treeIndices,
                 size_t numTreeIndices, bool useSaved,
                 const double* replayData, size_t replayNumRows,
                 size_t forestIndex, GatheredTrees& out) {
  const bartcore::ColumnStore& store(sampler.data());

  std::vector<bartcore::FlatNode> liveNodes;
  std::vector<double> liveSlopes;
  std::vector<std::uint64_t> liveMasks;
  std::vector<std::uint32_t> counts;
  std::vector<size_t> replayIndices(replayNumRows);
  std::string directionsScratch;
  bool functionLeaves = sampler.usesFunctionLeaves();
  // the mask side channel exists only for pooled columns (past 63 levels)
  bool anyPooled = false;
  for (size_t j = 0; j < store.numPredictors; ++j)
    if (store.columnIsPooled(j)) { anyPooled = true; break; }

  for (size_t i = 0; i < numChainIndices; ++i) {
    size_t chainNum = chainIndices[i];
    for (size_t j = 0; j < numSampleIndices; ++j) {
      size_t sampleNum = useSaved ? sampleIndices[j] : 0;
      for (size_t k = 0; k < numTreeIndices; ++k) {
        size_t treeNum = treeIndices[k];

        const std::vector<bartcore::FlatNode>* nodes;
        const std::vector<double>* slopes = NULL;
        const std::vector<std::uint64_t>* masks = NULL;
        if (useSaved) {
          nodes = &sampler.savedTree(chainNum, sampleNum, treeNum, forestIndex);
          if (out.numSlopes > 0)
            slopes = &sampler.savedTreeSlopes(chainNum, sampleNum, treeNum,
                                              forestIndex);
          if (anyPooled)
            masks = &sampler.savedTreeMasks(chainNum, sampleNum, treeNum,
                                            forestIndex);
        } else {
          sampler.flattenTree(chainNum, treeNum, liveNodes, counts,
                              out.numSlopes > 0 ? &liveSlopes : NULL,
                              anyPooled ? &liveMasks : NULL, forestIndex);
          nodes = &liveNodes;
          if (out.numSlopes > 0) slopes = &liveSlopes;
          if (anyPooled) masks = &liveMasks;
        }
        if (replayData != NULL) {
          counts.resize(nodes->size());
          for (size_t l = 0; l < replayNumRows; ++l) replayIndices[l] = l;
          bartcore::countFlatObservationsBelow(
            nodes->data(), replayData, replayNumRows, replayIndices.data(), 0,
            replayNumRows, counts.data(), masks != NULL ? masks->data() : NULL);
        } else if (useSaved) {
          // saved trees carry no counts; a null replay source (a sparse or
          // absent creation spec) leaves the n column zeroed rather than
          // reading whichever tree last filled the shared scratch
          counts.assign(nodes->size(), 0);
        }

        size_t leafNum = 0;
        for (size_t l = 0; l < nodes->size(); ++l) {
          out.chain.push_back(static_cast<int>(chainNum + 1));
          out.sample.push_back(static_cast<int>(sampleNum + 1));
          out.tree.push_back(static_cast<int>(treeNum + 1));
          out.count.push_back(static_cast<int>(counts[l]));
          const bartcore::FlatNode& node((*nodes)[l]);
          out.variable.push_back(
            node.variable >= 0 ? node.variable + 1 : node.variable);
          bool isCategoricalRule =
            node.variable >= 0 &&
            store.types[static_cast<size_t>(node.variable)] ==
              bartcore::ColumnType::categorical;
          // a categorical rule's payload is a mask, not a data value; the
          // directions column carries its decode and its value is NA
          out.value.push_back(
            isCategoricalRule ||
                (functionLeaves && node.variable == bartcore::invalidVariable)
              ? NA_REAL : node.value);
          if (out.anyCategorical) {
            if (isCategoricalRule) {
              size_t variable = static_cast<size_t>(node.variable);
              directionsScratch.clear();
              if (store.columnIsPooled(variable)) {
                const std::uint64_t* words = masks->data() + node.maskOffset;
                for (std::uint32_t c = 0; c < store.numCuts[variable]; ++c)
                  directionsScratch.push_back(
                    bartcore::maskTestBit(words, c) ? 'R' : 'L');
              } else {
                for (std::uint32_t c = 0; c < store.numCuts[variable]; ++c)
                  directionsScratch.push_back(
                    ((node.mask >> c) & 1u) != 0 ? 'R' : 'L');
              }
              out.directions.push_back(directionsScratch);
            } else {
              out.directions.push_back(std::string());
            }
          }
          if (out.anyMissing)
            out.missing.push_back(
              node.variable >= 0 &&
                  store.hasMissing[static_cast<size_t>(node.variable)]
                ? static_cast<int>(node.flags & bartcore::flatMissingGoesRight)
                : NA_INTEGER);
          if (out.numSlopes > 0) {
            bool isLeaf = node.variable == bartcore::invalidVariable;
            for (size_t s = 0; s < out.numSlopes; ++s)
              out.slopes[s].push_back(
                isLeaf ? (*slopes)[leafNum * out.numSlopes + s] : NA_REAL);
            if (isLeaf) ++leafNum;
          }
        }
      }
    }
  }
}

// Appends an integer or numeric column to the data frame under
// construction, copying from the gathered vector.
template <typename T>
void emitTreeColumn(SEXP resultExpr, SEXP namesExpr, R_xlen_t columnNum,
                    const char* name, const std::vector<T>& column) {
  R_xlen_t length = static_cast<R_xlen_t>(column.size());
  if constexpr (std::is_same_v<T, int>) {
    SET_VECTOR_ELT(resultExpr, columnNum, Rf_allocVector(INTSXP, length));
    std::memcpy(INTEGER(VECTOR_ELT(resultExpr, columnNum)), column.data(),
                column.size() * sizeof(int));
  } else {
    SET_VECTOR_ELT(resultExpr, columnNum, Rf_allocVector(REALSXP, length));
    std::memcpy(REAL(VECTOR_ELT(resultExpr, columnNum)), column.data(),
                column.size() * sizeof(double));
  }
  SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar(name));
}

// Builds the classic-format data.frame from a gather: ([chain,] [sample,]
// tree, n, var, value[, directions][, missing][, beta.*]).
SEXP emitTreeDataFrame(const GatheredTrees& gathered) {
  R_xlen_t totalNumNodes = static_cast<R_xlen_t>(gathered.value.size());
  R_xlen_t numColumns = 4 + (gathered.includeChain ? 1 : 0) +
                        (gathered.includeSample ? 1 : 0) +
                        (gathered.anyCategorical ? 1 : 0) +
                        (gathered.anyMissing ? 1 : 0) +
                        static_cast<R_xlen_t>(gathered.numSlopes);

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, numColumns));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, numColumns));

  R_xlen_t columnNum = 0;
  if (gathered.includeChain)
    emitTreeColumn(resultExpr, namesExpr, columnNum++, "chain",
                   gathered.chain);
  if (gathered.includeSample)
    emitTreeColumn(resultExpr, namesExpr, columnNum++, "sample",
                   gathered.sample);
  emitTreeColumn(resultExpr, namesExpr, columnNum++, "tree", gathered.tree);
  emitTreeColumn(resultExpr, namesExpr, columnNum++, "n", gathered.count);
  emitTreeColumn(resultExpr, namesExpr, columnNum++, "var", gathered.variable);
  emitTreeColumn(resultExpr, namesExpr, columnNum++, "value", gathered.value);
  if (gathered.anyCategorical) {
    SET_VECTOR_ELT(resultExpr, columnNum,
                   Rf_allocVector(STRSXP, totalNumNodes));
    SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar("directions"));
    SEXP directionsExpr = VECTOR_ELT(resultExpr, columnNum);
    for (R_xlen_t l = 0; l < totalNumNodes; ++l)
      SET_STRING_ELT(directionsExpr, l,
                     gathered.directions[static_cast<size_t>(l)].empty()
                       ? NA_STRING
                       : Rf_mkChar(gathered.directions[static_cast<size_t>(l)]
                                     .c_str()));
    ++columnNum;
  }
  if (gathered.anyMissing)
    emitTreeColumn(resultExpr, namesExpr, columnNum++, "missing",
                   gathered.missing);
  // generically named here; the R wrapper renames to beta.<column name>
  for (size_t s = 0; s < gathered.numSlopes; ++s) {
    char slopeName[32];
    snprintf(slopeName, sizeof(slopeName), "beta.%lu",
             static_cast<unsigned long>(s + 1));
    emitTreeColumn(resultExpr, namesExpr, columnNum++, slopeName,
                   gathered.slopes[s]);
  }

  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  SEXP rowNamesExpr = PROTECT(Rf_allocVector(STRSXP, totalNumNodes));
  char buffer[32];
  for (R_xlen_t l = 0; l < totalNumNodes; ++l) {
    snprintf(buffer, sizeof(buffer), "%ld", static_cast<long>(l + 1));
    SET_STRING_ELT(rowNamesExpr, l, Rf_mkChar(buffer));
  }
  Rf_setAttrib(resultExpr, R_RowNamesSymbol, rowNamesExpr);

  SEXP classExpr = PROTECT(Rf_mkString("data.frame"));
  Rf_setAttrib(resultExpr, R_ClassSymbol, classExpr);

  UNPROTECT(4);
  return resultExpr;
}

// A data.frame of tree structure in the classic engine's format: pre-order
// rows of ([chain,] [sample,] tree, n, var, value), var 1-based with -1
// marking leaves, an ordinal rule's value its cut point, leaf values on the
// engine's internal response scale. A categorical rule carries no data value:
// its value is NA and a 'directions' character column - present whenever the
// store has any categorical column - carries one L/R per observed category
// (NA on ordinal rules and leaves; the R wrapper pads to the declared level
// count). When any predictor contains missing values a 'missing' integer
// column reports each rule's missing direction (0 left, 1 right; NA on leaves
// and on columns without missing values). Saved trees replay the training
// predictors for n unless newdata is supplied; live trees report their own
// counts. Validates, gathers, and emits, in that order.
SEXP getTrees(bartcore::SamplerBase& sampler, const size_t* chainIndices,
              size_t numChainIndices, const size_t* sampleIndices,
              size_t numSampleIndices, const size_t* treeIndices,
              size_t numTreeIndices, bool useLiveTrees, const double* newdata,
              size_t newdataNumRows, const double* trainingReplay,
              size_t trainingReplayNumRows, size_t forestIndex,
              const char* caller) {
  const bartcore::ColumnStore& store(sampler.data());

  bool useSaved = sampler.savedTreeCapacity() > 0 && !useLiveTrees;
  if (!useSaved) numSampleIndices = 1;

  for (size_t i = 0; i < numChainIndices; ++i) {
    if (chainIndices[i] >= sampler.numChains())
      Rf_error("%s chain number out of range", caller);
  }
  if (useSaved) {
    for (size_t i = 0; i < numSampleIndices; ++i) {
      if (sampleIndices[i] >= sampler.savedTreeCapacity())
        Rf_error("%s sample number out of range", caller);
    }
  }
  for (size_t i = 0; i < numTreeIndices; ++i) {
    if (treeIndices[i] >= sampler.numTreesInForest(forestIndex))
      Rf_error("%s tree number out of range", caller);
  }

  // saved trees carry no counts of their own and replay the training rows the
  // caller supplies (the engine keeps no matrix); newdata replays its rows
  // through live and saved trees alike
  const double* replayData = NULL;
  size_t replayNumRows = 0;
  if (newdata != NULL) {
    replayData = newdata;
    replayNumRows = newdataNumRows;
  } else if (useSaved) {
    replayData = trainingReplay;
    replayNumRows = trainingReplayNumRows;
  }

  bool anyMissing = false, anyCategorical = false;
  for (size_t j = 0; j < store.numPredictors; ++j) {
    if (store.hasMissing[j]) anyMissing = true;
    if (store.types[j] == bartcore::ColumnType::categorical)
      anyCategorical = true;
  }

  // function-valued (gp) leaves report no per-leaf coefficients - the
  // function rides prediction only - and their leaves' values print NA
  // (a whole function per row does not fit a data frame)
  size_t numSlopes =
    sampler.usesFunctionLeaves() ? 0 : sampler.numLeafCovariates();

  GatheredTrees gathered{sampler.numChains() > 1, useSaved,
                         anyCategorical, anyMissing, numSlopes,
                         {}, {}, {}, {}, {}, {}, {}, {}, {}};
  gathered.slopes.resize(numSlopes);
  gatherTrees(sampler, chainIndices, numChainIndices, sampleIndices,
              numSampleIndices, treeIndices, numTreeIndices, useSaved,
              replayData, replayNumRows, forestIndex, gathered);

  return emitTreeDataFrame(gathered);
}

} // namespace bartcore_bridge
