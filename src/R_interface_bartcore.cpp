#include "config.hpp"
#include "R_interface_bartcore.hpp"

#include <climits> // INT_MAX
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <memory>
#include <numbers>
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
using bartcore_bridge::AugmentationInputs;
using bartcore_bridge::augmentationFamily;
using bartcore_bridge::BartcoreHolder;
using bartcore_bridge::computeWorkingResponse;
using bartcore_bridge::drawAugmentation;
using bartcore_bridge::enforceBinaryWeightPolicy;
using bartcore_bridge::refuseUndefinedTestFits;
using bartcore_bridge::refuseBinaryWeightChange;
using bartcore_bridge::refuseCscReferenceAgainstStore;
using bartcore_bridge::refuseEmptyTreeStore;
using bartcore_bridge::refuseGroupedScaleUpdate;
using bartcore_bridge::refuseMultiForestMutation;
using bartcore_bridge::refuseMultiForestResponseMutation;
using bartcore_bridge::refuseNonBinaryMask;
using bartcore_bridge::refusePinnedSigmaChange;
using bartcore_bridge::refuseSparseLeafCovariate;
using bartcore_bridge::refuseVarianceForestScaleUpdate;
using bartcore_bridge::ResponseConduit;
using bartcore_bridge::validateColumnValues;
using bartcore_bridge::validateResponseSupport;
using bartcore_bridge::validateTestContainerAgainstStore;

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

// Refusal shared by every storage = "single" gate: fp32 residual storage is
// v1-scoped to the gaussian constant-leaf path.
static const char storageSingleUnsupportedMessage[] =
  "storage = \"single\" is currently supported only for gaussian models "
  "with constant leaves";

// Tolerance for a user-supplied probability vector that must sum to 1.0:
// sqrt(DBL_EPSILON) = 2^-26 = 1.4901161193847656e-08, which R's all.equal,
// tinytest::expect_equal and testthat all default to when asked whether two
// numbers are the same; shares combiner.hpp's zeroMultiplierTolerance
// derivation. Written as a hex literal because 2^-26 is exact.
static constexpr double sumToOneTolerance = 0x1p-26;

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

// The bartcore_create* entry points share this holder-installation boilerplate:
// register the finalizer on a null-address external pointer first
// (R_SetExternalPtrAddr does not allocate, so no OOM longjmp can strand the
// holder between its construction and its finalizer), build the holder, then
// install its address. \p buildHolder returns the holder; its own parse errors
// longjmp past this before any address is installed. The protection vector pins
// dataExpr (PROT_DATA); the setters fill the remaining slots later.
template <typename BuildHolder>
SEXP createExternalHolder(SEXP dataExpr, BuildHolder buildHolder) {
  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  BartcoreHolder* holder = buildHolder();
  R_SetExternalPtrAddr(result, holder);

  UNPROTECT(2);
  return result;
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
  size_t numPredictors = sampler.shape().numPredictors;
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  if (!Rf_isReal(xExpr) || Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[1]) != numPredictors)
    Rf_error("%s: requires a numeric matrix with matching columns", caller);
  size_t numRows = static_cast<size_t>(INTEGER(dims)[0]);
  for (size_t j = 0; j < numPredictors; ++j)
    validateColumnValues(sampler.data(), j, REAL(xExpr) + j * numRows,
                         numRows);
  return numRows;
}

// The subset of the R specification objects (dbartsControl, dbartsData,
// dbartsModel) the engine consumes. Pointers borrow from the expressions;
// error paths may leak the parse vectors - Rf_error longjmps past
// destructors, so the leak is deliberate and bounded by the R error.

struct ParsedControl {
  bool responseIsBinary = false;
  // ordinal (cumulative-probit) response shape: the ordered category count
  // K, the third response shape beside the binary and continuous ones. K >=
  // 2 selects the ordinal family; 0 (absent) is a non-ordinal response, so
  // every existing family parses unchanged. The R surface attaches the count
  // (the bartcore.n.categories control attribute); parseControl reads it
  // here.
  size_t numOrdinalCategories = 0;
  // negative-binomial count response shape: the fourth response shape
  // beside the binary, ordinal, and continuous ones. countResponse marks it
  // (a count is none of the others);
  // dispersion is the r spec, the residualDf sign convention - a positive value
  // fixes r there (an integer), a non-positive value estimates r on the grid.
  // The dbarts()/bart2 surface attaches both through one control attribute
  // (bartcore.dispersion); parseControl reads them here. Absent (the default) leaves a non-count response, so every
  // existing family parses byte-for-byte unchanged.
  bool countResponse = false;
  double dispersion = NA_REAL;
  bool verbose = false;
  bool keepTrainingFits = true;
  bool useQuantiles = false;
  bool keepTrees = false;
  // opt-in fp32 running residual (control@storage == "single"); the
  // createSampler gate refuses it for anything but a gaussian constant-leaf
  // model
  bool fp32Residual = false;
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
  std::vector<std::int32_t> columnSources;
  std::vector<bartcore::xint_t> cscReferenceCodes;
  // the borrowed view the test build takes, over the buffers above and the
  // container's sparse part; published once each is final
  bartcore::PredictorSource view;
};

// A sparse mutation argument's parsed storage: the borrowed view, the buffers
// it points into, the store types of the columns it names, and the dense block
// the engine's re-quantize reads. Held by value in the entrance's
// unwindProtect frame, so the validation jump frees it.
struct ParsedMutationSource {
  std::vector<double> denseAssembly;
  std::vector<std::int32_t> columnSources;
  std::vector<bartcore::xint_t> referenceCodes;
  std::vector<bartcore::ColumnType> storeTypes;
  std::vector<double> block;
  // the argument's per-sparse-column reference metadata, borrowed from the
  // container; NULL for a bare dgCMatrix, whose implicit rows are the zero its
  // own storage means
  const int* referenceMeta = NULL;
  const int* categoryCountMeta = NULL;
  size_t numSparseColumns = 0;
  bartcore::PredictorSource view;
};

struct ParsedData {
  const double* y = NULL;
  size_t numObservations = 0;
  size_t numPredictors = 0;
  std::vector<bartcore::ColumnType> columnTypes;
  bool anyCategorical = false;
  // The borrowed view every training entrance reads: a plain matrix and the
  // dense columnar container leave it unmapped over denseAssembly or the R
  // matrix, while a dgCMatrix and the mixed container map each column onto the
  // CSC triple (negative) or the dense block (nonnegative). Published field by
  // field below, as each owning buffer becomes final.
  bartcore::PredictorSource predictors;
  // the storage FLAVOR the parse saw, kept for the refusal texts and the
  // CSC code-validation memory-safety gates that must not key on the map
  bool xIsSparse = false;
  bool xIsMixed = false;
  std::vector<std::int32_t> columnSources;
  // per sparse column of a mixed container: the reference level's 0-based
  // code and the level count K, borrowed from the container (null unless it
  // carries a sparse block, whose metadata the parse validates). Resolved per
  // predictor into cscReferenceCodes / categoryCounts once varTypes marks the
  // CSC-backed categorical columns.
  const int* cscReferenceMeta = NULL;
  const int* cscCategoryCountMeta = NULL;
  size_t numSparseColumns = 0;
  // per predictor, the level count its host declares (0 = infer from the
  // codes): a CSC-backed categorical column's from the container metadata
  // above, a dense-backed one's from the attr(x, "factor.levels") table the
  // model-matrix builders attach. Empty when nothing declares one.
  std::vector<std::uint32_t> categoryCounts;
  std::vector<bartcore::xint_t> cscReferenceCodes;
  // a columnar container's transiently assembled dense block (the view's
  // denseValues either flavor); owned here so it lives exactly as long as the
  // parse result, which is all either build needs - an unmapped build
  // quantizes into owned codes and a mapped one copies the block
  std::vector<double> denseAssembly;
  // the test view, unmapped over the x.test matrix or filled from the parsed
  // container below; the test store copies its raw either way
  bartcore::PredictorSource testPredictors;
  size_t numTestObservations = 0;
  bool testIsMixed = false;
  ParsedTestContainer testContainer;
  const double* weights = NULL;
  const double* offset = NULL;
  const double* testOffset = NULL;
  // the per-forest amplitude bases a multi-forest fit's amplitudes
  // multiply, EMPTY for every single-forest model. Entry f is forest f's
  // n x basisColumns[f] matrix as R holds it - COLUMN-major -
  // borrowed from the data object and transposed into the engine's row-major
  // contract at the build; a null entry leaves forest f the dense all-ones
  // column the engine synthesizes.
  std::vector<const double*> bases;
  std::vector<size_t> basisColumns;
  double sigmaEstimate = 1.0;
  std::vector<uint32_t> maxNumCuts;
};

struct ParsedModel {
  double birthOrDeathProbability = 0.5;
  double swapProbability = 0.1;
  double changeProbability = 0.4;
  double birthProbability = 0.5;
  double nodeScale = 0.5;
  // the model's prior.scale slot: the NAMED leaf calibration in response
  // units, NA when unnamed. Finite, the single-forest engine converts it
  // against the response transform and it overrides nodeScale; the two-forest
  // and multinomial creation paths refuse it by name below.
  double priorScale = NA_REAL;
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
  // Student-t residual degrees of freedom: NaN is the Gaussian default, a
  // positive value fixes nu, 0 estimates it on the grid. The dbarts()/bart2
  // surface attaches it (the model's resid.df attribute); absent it stays
  // NaN.
  double residualDf = NA_REAL;
  // per-predictor monotone directions in {-1, 0, +1}, narrowed from the R
  // integer spec; empty when no constraint is declared
  std::vector<std::int8_t> monotoneDirections;
  // per-forest interaction constraint: interactionMaxOrder caps the
  // distinct split variables on any path (0 = uncapped);
  // interactionForbiddenPairs is a flat 0-based (a, b) stream (two indices
  // per forbidden co-occurrence). Both defaults leave the path free.
  size_t interactionMaxOrder = 0;
  std::vector<size_t> interactionForbiddenPairs;
  // per-forest block-additive constraint (per-tree column grouping): blockOfColumn is the 0-based
  // group index per predictor (negative = in no block), blockTreeCounts the
  // per-group tree capacity summing to numTrees. Empty leaves every tree
  // unrestricted, byte-for-byte the default path.
  std::vector<std::int32_t> blockOfColumn;
  std::vector<size_t> blockTreeCounts;
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

// Attach \p value under the attribute named \p name. Interning a name the
// process has not seen allocates, and C++ leaves the evaluation order of a
// call's arguments unspecified, so a freshly allocated value passed as a
// sibling argument of Rf_install could be collected before it is attached;
// rooting it here orders the two.
void setAttribByName(SEXP target, const char* name, SEXP value) {
  PROTECT(value);
  Rf_setAttrib(target, Rf_install(name), value);
  UNPROTECT(1);
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

  // storage precision: "single" opts the running residual into fp32; anything
  // else (default "double", or an old control lacking the slot) keeps the
  // fp64 engine
  REPROTECT_SLOT(slotExpr, controlExpr, "storage", slotIndex);
  control.fp32Residual =
    Rf_isString(slotExpr) && rc_getLength(slotExpr) >= 1 &&
    STRING_ELT(slotExpr, 0) != NA_STRING &&
    std::strcmp(CHAR(STRING_ELT(slotExpr, 0)), "single") == 0;

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
  int printValue = rc_getInt(slotExpr, "print every", RC_LENGTH | RC_EQ,
                             rc_asRLength(1), RC_VALUE | RC_GEQ, 1,
                             RC_NA | RC_YES, RC_END);
  if (printValue != NA_INTEGER)
    control.printEvery = static_cast<uint32_t>(printValue);

  REPROTECT_SLOT(slotExpr, controlExpr, "printCutoffs", slotIndex);
  printValue = rc_getInt(slotExpr, "print cutoffs", RC_LENGTH | RC_EQ,
                         rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES,
                         RC_END);
  if (printValue == NA_INTEGER) printValue = 0;
  control.printCutoffs = static_cast<uint32_t>(printValue);

  // rngKind and rngNormalKind are not consumed by this engine; the R side
  // refuses them before creation, so they are not read here
  REPROTECT_SLOT(slotExpr, controlExpr, "seed", slotIndex);
  if (rc_getLength(slotExpr) != 1)
    Rf_error("slot 'seed' must be of length 1");
  control.haveRngSeed = INTEGER(slotExpr)[0] != NA_INTEGER;
  if (control.haveRngSeed)
    control.rngSeed = static_cast<std::uint_least32_t>(INTEGER(slotExpr)[0]);

  // optional ordinal shape: an integer K >= 2 the R surface attaches for an
  // ordered-factor response, read raw and guarded like resid.df. Absent (the
  // default) leaves a non-ordinal response, so every existing family parses
  // byte-for-byte unchanged.
  // An attribute of the rooted controlExpr cannot be collected, so the PROTECT
  // is redundant to that rooting; it is what the PROTECT-balance analyzer reads,
  // which sees only that Rf_isInteger may allocate (it tests for a factor).
  SEXP ordinalExpr =
    PROTECT(Rf_getAttrib(controlExpr, Rf_install("bartcore.n.categories")));
  if (Rf_isInteger(ordinalExpr) && Rf_xlength(ordinalExpr) == 1 &&
      INTEGER(ordinalExpr)[0] >= 2)
    control.numOrdinalCategories =
      static_cast<size_t>(INTEGER(ordinalExpr)[0]);
  UNPROTECT(1);

  // optional count shape: a length-1 real the R surface attaches for a
  // count response, guarded like bartcore.n.categories. Its presence marks
  // the count shape; its value is the dispersion spec (positive fixes r,
  // non-positive estimates on the grid). Absent (the default) leaves a
  // non-count response, so every existing family parses byte-for-byte
  // unchanged.
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
// \p rowCountMessage carries the caller's wording for a disagreeing row count:
// the default is creation's, which must not surface from a mutation call.
CscSlots parseCscMatrix(SEXP matrixExpr, size_t numObservations,
                        const char* rowCountMessage =
                          "number of rows of 'x' must equal length of 'y'") {
  CscSlots result;
  SEXP dimExpr = PROTECT(Rf_getAttrib(matrixExpr, Rf_install("Dim")));
  if (!Rf_isInteger(dimExpr) || rc_getLength(dimExpr) != 2)
    Rf_error("malformed sparse predictor matrix");
  if (static_cast<size_t>(INTEGER(dimExpr)[0]) != numObservations)
    Rf_error("%s", rowCountMessage);
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
    // NA_INTEGER is INT_MIN, whose negation overflows, so the CSC arm bounds
    // the entry before it flips the sign rather than after
    if (map[j] > 0 && static_cast<size_t>(map[j]) <= numDenseColumns)
      out[j] = map[j] - 1;
    else if (map[j] < 0 && map[j] >= -static_cast<std::int64_t>(numCscColumns))
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

// Validate a mixed container's per-sparse-column reference metadata and borrow
// it: both fields must be integer vectors of one entry per CSC column, in the
// sparse block's column order - how every consumer indexes them. Shared by all
// three container funnels (creation and either mutation), each passing its own
// \p malformedMessage, so a container carrying a sparse block is refused where
// it arrives rather than taken at one entrance and refused at the next. The
// borrowed pointers ride the container, valid as long as it is protected.
void requireCscReferenceMeta(SEXP containerExpr, size_t numCscColumns,
                             const int*& referenceMeta,
                             const int*& categoryCountMeta,
                             const char* malformedMessage) {
  SEXP referenceExpr = rc_getListElement(containerExpr, "sparseReference");
  SEXP categoryCountExpr =
    rc_getListElement(containerExpr, "sparseCategoryCount");
  if (!Rf_isInteger(referenceExpr) || !Rf_isInteger(categoryCountExpr) ||
      static_cast<size_t>(rc_getLength(referenceExpr)) != numCscColumns ||
      static_cast<size_t>(rc_getLength(categoryCountExpr)) != numCscColumns)
    Rf_error("%s", malformedMessage);
  referenceMeta = INTEGER(referenceExpr);
  categoryCountMeta = INTEGER(categoryCountExpr);
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
  SEXP denseExpr = PROTECT(rc_getListElement(containerExpr, "dense"));
  SEXP sparseExpr = PROTECT(rc_getListElement(containerExpr, "sparse"));
  SEXP mapExpr = PROTECT(rc_getListElement(containerExpr, "map"));
  if (!Rf_isInteger(mapExpr) ||
      static_cast<size_t>(rc_getLength(mapExpr)) != numPredictors)
    Rf_error("number of columns in 'x.test' must equal that of 'x'");
  SEXP numObsExpr = rc_getListElement(containerExpr, "numObservations");
  if (!Rf_isInteger(numObsExpr) || rc_getLength(numObsExpr) != 1 ||
      INTEGER(numObsExpr)[0] < 0)
    Rf_error("malformed mixed test container");
  size_t numTest = static_cast<size_t>(INTEGER(numObsExpr)[0]);
  out.view.numRows = numTest;
  out.view.numColumns = numPredictors;
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
                           out.view.denseValues, out.columnSources,
                           "number of rows of 'x.test' columns must match",
                           "malformed mixed test container");
  out.view.columnSources = out.columnSources.data();
  if (hasSparse) {
    out.view.cscColumnPointers = csc.pointers;
    out.view.cscRowIndices = csc.rows;
    out.view.cscValues = csc.values;
  }

  // Reference metadata is present whenever the container carries a sparse
  // block, independent of whether any CSC column happens to be
  // store-categorical: a non-NA reference against a store-ORDINAL column is
  // malformed at every funnel (rule 2), so the refusal below must not be
  // gated on a categorical column existing among this container's columns.
  const int* referenceMeta = nullptr;
  const int* categoryCountMeta = nullptr;
  if (hasSparse)
    requireCscReferenceMeta(containerExpr, numCscColumns, referenceMeta,
                            categoryCountMeta,
                            "sparse categorical test predictor columns require "
                            "reference metadata");
  // a non-NA reference against a store-ORDINAL column is refused
  // here, for every parseTestContainer caller (creation, setTestPredictor,
  // setTestPredictorAndOffset) alike, since columnTypes is the STORE's types
  // at mutation and the being-built types at creation
  refuseCscReferenceAgainstStore(columnTypes, out.columnSources.data(),
                                 numPredictors, referenceMeta, numCscColumns);

  // resolve the reference code per CSC-backed categorical test column (the code
  // its implicit rows take); the container's per-sparse-column metadata is
  // already in level order
  bool anyTestCscCategorical = false;
  for (size_t j = 0; j < numPredictors; ++j)
    if (columnTypes[j] == bartcore::ColumnType::categorical &&
        out.view.sourceOf(j) < 0) {
      anyTestCscCategorical = true;
      break;
    }
  if (anyTestCscCategorical) {
    resolveCscCategoricalReferences(
      columnTypes, out.columnSources.data(), numPredictors, referenceMeta,
      categoryCountMeta, numCscColumns, out.cscReferenceCodes, nullptr,
      "malformed mixed test container",
      "sparse categorical test predictor columns require reference "
      "metadata");
    out.view.referenceCodes = out.cscReferenceCodes.data();
  }
  UNPROTECT(3);
}

// Whether a test-side argument carries its values sparsely, so it takes the
// resident container parse rather than the dense matrix one.
bool testSourceIsSparse(SEXP xTestExpr) {
  return Rf_inherits(xTestExpr, "dbartsMixedMatrix") ||
         Rf_inherits(xTestExpr, "dgCMatrix");
}

// A sparse test set at a read-only entrance (predict, tree replay): a mixed
// container parses as at every other test funnel, and a bare dgCMatrix is the
// all-CSC spelling of the same container - every column CSC-backed, implicit
// rows the zero its own storage means, no reference metadata to declare.
void parseTestSource(ParsedTestContainer& out, SEXP xTestExpr,
                     size_t numPredictors,
                     const bartcore::ColumnType* storeTypes) {
  if (!Rf_inherits(xTestExpr, "dgCMatrix")) {
    parseTestContainer(out, xTestExpr, numPredictors, storeTypes);
    return;
  }
  SEXP dimExpr = PROTECT(Rf_getAttrib(xTestExpr, Rf_install("Dim")));
  if (!Rf_isInteger(dimExpr) || rc_getLength(dimExpr) != 2)
    Rf_error("malformed sparse predictor matrix");
  size_t numTest = static_cast<size_t>(INTEGER(dimExpr)[0]);
  UNPROTECT(1);
  CscSlots csc = parseCscMatrix(xTestExpr, numTest);
  if (csc.numColumns != numPredictors)
    Rf_error("number of columns in 'x.test' must equal that of 'x'");
  out.columnSources.resize(numPredictors);
  for (size_t j = 0; j < numPredictors; ++j)
    out.columnSources[j] = ~static_cast<std::int32_t>(j);
  out.view.numRows = numTest;
  out.view.numColumns = numPredictors;
  out.view.columnSources = out.columnSources.data();
  out.view.cscColumnPointers = csc.pointers;
  out.view.cscRowIndices = csc.rows;
  out.view.cscValues = csc.values;
}

// Mutation-side wording for the shared parse helpers; the creation-flavored
// texts they default to describe a design being built, not one being replaced.
const char* const mutationContainerMessage =
  "malformed mixed predictor container";
const char* const mutationReferenceMessage =
  "sparse categorical predictor columns require reference metadata";

// Parse a sparse mutation argument - a bare dgCMatrix or a dbartsMixedMatrix -
// into a borrowed view of \p numRows rows over the \p numColumns predictor
// columns it replaces. Returns false with \p out untouched for anything else,
// leaving the entrance to refuse it; \p shapeMessage is the entrance's own
// wording for a disagreeing shape.
bool parseMutationSource(ParsedMutationSource& out, SEXP xExpr, size_t numRows,
                         size_t numColumns, const char* shapeMessage) {
  if (Rf_inherits(xExpr, "dgCMatrix")) {
    CscSlots csc = parseCscMatrix(xExpr, numRows, shapeMessage);
    if (csc.numColumns != numColumns) Rf_error("%s", shapeMessage);
    out.columnSources.resize(numColumns);
    for (size_t j = 0; j < numColumns; ++j)
      out.columnSources[j] = ~static_cast<std::int32_t>(j);
    out.view.cscColumnPointers = csc.pointers;
    out.view.cscRowIndices = csc.rows;
    out.view.cscValues = csc.values;
  } else if (Rf_inherits(xExpr, "dbartsMixedMatrix")) {
    SEXP denseExpr = PROTECT(rc_getListElement(xExpr, "dense"));
    SEXP sparseExpr = PROTECT(rc_getListElement(xExpr, "sparse"));
    SEXP mapExpr = PROTECT(rc_getListElement(xExpr, "map"));
    if (!Rf_isInteger(mapExpr) ||
        static_cast<size_t>(rc_getLength(mapExpr)) != numColumns)
      Rf_error("%s", shapeMessage);
    bool hasSparse = Rf_inherits(sparseExpr, "dgCMatrix");
    if ((!Rf_isNull(sparseExpr) && !hasSparse) ||
        (!Rf_isNull(denseExpr) && TYPEOF(denseExpr) != VECSXP))
      Rf_error("%s", mutationContainerMessage);
    CscSlots csc;
    if (hasSparse) csc = parseCscMatrix(sparseExpr, numRows, shapeMessage);
    parseMixedContainerBlock(denseExpr, INTEGER(mapExpr), numColumns, numRows,
                             csc.numColumns, out.denseAssembly,
                             out.view.denseValues, out.columnSources,
                             shapeMessage, mutationContainerMessage);
    if (hasSparse) {
      out.view.cscColumnPointers = csc.pointers;
      out.view.cscRowIndices = csc.rows;
      out.view.cscValues = csc.values;
      requireCscReferenceMeta(xExpr, csc.numColumns, out.referenceMeta,
                              out.categoryCountMeta, mutationContainerMessage);
      out.numSparseColumns = csc.numColumns;
    }
    UNPROTECT(3);
  } else {
    return false;
  }
  out.view.numRows = numRows;
  out.view.numColumns = numColumns;
  out.view.columnSources = out.columnSources.data();
  return true;
}

// Materialize a parsed mutation argument into the dense block the engine's
// re-quantize takes, under the STORE's implicit rule: \p storeTypes gives the
// store's type of each column the argument names. A declared reference against
// a non-categorical store column is refused before anything is read.
const double* materializeMutationSource(
    ParsedMutationSource& parsed, const bartcore::ColumnType* storeTypes) {
  size_t numColumns = parsed.view.numColumns;
  refuseCscReferenceAgainstStore(storeTypes, parsed.columnSources.data(),
                                 numColumns, parsed.referenceMeta,
                                 parsed.numSparseColumns);
  if (parsed.referenceMeta != NULL) {
    resolveCscCategoricalReferences(
      storeTypes, parsed.columnSources.data(), numColumns,
      parsed.referenceMeta, parsed.categoryCountMeta, parsed.numSparseColumns,
      parsed.referenceCodes, nullptr, mutationContainerMessage,
      mutationReferenceMessage);
    parsed.view.referenceCodes = parsed.referenceCodes.data();
  }
  parsed.block.resize(parsed.view.numRows * numColumns);
  bartcore::materializePredictorSource(parsed.view, storeTypes, 0,
                                       parsed.view.numRows,
                                       parsed.block.data());
  return parsed.block.data();
}

// Route a parsed test container to the engine's typed test store (against the
// training cut grid, owning its raw). Returns false, store untouched, when a
// designated leaf covariate column would be CSC-backed - sparse storage serves
// no dense raw test covariate; the caller raises the leaf-covariate error.
bool installTestContainer(bartcore::SamplerBase& sampler,
                          const ParsedTestContainer& parsed) {
  return sampler.setTestData(parsed.view);
}

// Take each dense-backed categorical column's declared level count from the
// per-column level tables attr(x, "factor.levels") carries (one list element
// per predictor, NULL where the column is not a factor). The gate is
// varTypes: an ordered factor carries a level table too, but enters as an
// ordinal column, whose grid is cut points rather than categories. Absent -
// no attribute, a wrong-length one, or a NULL entry - leaves the count 0, the
// spelling resolveCscCategoricalReferences already uses for "infer from the
// observed codes". CSC-backed columns are skipped: their container declares
// its own K, which the caller has already resolved.
void readDeclaredCategoryCounts(ParsedData& data, SEXP factorLevelsExpr) {
  if (data.categoryCounts.empty())
    data.categoryCounts.assign(data.numPredictors, 0);
  if (TYPEOF(factorLevelsExpr) != VECSXP ||
      static_cast<size_t>(rc_getLength(factorLevelsExpr)) != data.numPredictors)
    return;
  for (size_t j = 0; j < data.numPredictors; ++j) {
    if (data.columnTypes[j] != bartcore::ColumnType::categorical) continue;
    if (data.predictors.sourceOf(j) < 0) continue;
    R_xlen_t numLevels =
      Rf_xlength(VECTOR_ELT(factorLevelsExpr, static_cast<R_xlen_t>(j)));
    // a table wider than the code type can hold declares nothing usable; the
    // inferred count then applies and the training bound stays maxCategories
    if (numLevels > 0 &&
        static_cast<size_t>(numLevels) <= bartcore::maxCategories)
      data.categoryCounts[j] = static_cast<std::uint32_t>(numLevels);
  }
}

// Reads a list of per-forest amplitude bases into the parsed data: entry f is
// forest f's n x q_f numeric matrix as R holds it, COLUMN-major and BORROWED,
// and a NULL entry is a forest that declares none (the dense all-ones column
// the engine synthesizes). A null or empty list leaves the channel empty,
// which is what marks a fit single-forest. Shared by the data-object slot and
// the internal creation entry, so the two read the same shape.
void readForestBases(SEXP basesExpr, ParsedData& data) {
  data.bases.clear();
  data.basisColumns.clear();
  if (Rf_isNull(basesExpr) || Rf_xlength(basesExpr) == 0) return;
  if (!Rf_isNewList(basesExpr))
    Rf_error("'bases' must be a list of per-forest basis matrices");
  R_xlen_t numForests = Rf_xlength(basesExpr);
  data.bases.resize(static_cast<size_t>(numForests), NULL);
  data.basisColumns.resize(static_cast<size_t>(numForests), 0);
  for (R_xlen_t f = 0; f < numForests; ++f) {
    SEXP basisExpr = VECTOR_ELT(basesExpr, f);
    if (Rf_isNull(basisExpr)) continue;
    if (!Rf_isReal(basisExpr))
      Rf_error("a 'bases' element must be a numeric matrix");
    R_xlen_t length = Rf_xlength(basisExpr);
    if (length == 0 || data.numObservations == 0 ||
        length % static_cast<R_xlen_t>(data.numObservations) != 0)
      Rf_error("a 'bases' element must have as many rows as 'y' has elements");
    data.bases[static_cast<size_t>(f)] = REAL(basisExpr);
    data.basisColumns[static_cast<size_t>(f)] =
      static_cast<size_t>(length) / data.numObservations;
  }
}

void parseData(ParsedData& data, SEXP dataExpr) {
  SEXP slotExpr;
  PROTECT_INDEX slotIndex;
  PROTECT_WITH_INDEX(R_NilValue, &slotIndex);

  REPROTECT_SLOT(slotExpr, dataExpr, "y", slotIndex);
  if (!Rf_isReal(slotExpr)) Rf_error("'y' must be numeric");
  if (rc_getLength(slotExpr) <= 0)
    Rf_error("length of y must be greater than 0");
  data.y = REAL(slotExpr);
  data.numObservations = rc_getLength(slotExpr);
  // The hot gather index is stored as bartcore::index_t (uint32); a subscript
  // n - 1 must fit. Unreachable in practice (predictor codes would be
  // terabytes first), but must not silently truncate.
  if (data.numObservations > static_cast<size_t>(UINT32_MAX))
    Rf_error("number of observations exceeds the %u index limit", UINT32_MAX);

  REPROTECT_SLOT(slotExpr, dataExpr, "x", slotIndex);
  // the per-column level tables the model-matrix builders attach (attr
  // "factor.levels", R/utility.R): one list element per predictor, NULL where
  // the column is not a factor, and present on a plain matrix, the dense
  // columnar container, and the mixed container alike. Read once here, while
  // slotExpr is still x, and resolved into declared category counts below,
  // where varTypes says which columns are categorical. Protected outright
  // rather than leaned on dataExpr's own protection, since the slot is
  // re-pointed at every field parsed after it.
  SEXP factorLevelsExpr =
    PROTECT(Rf_getAttrib(slotExpr, Rf_install("factor.levels")));
  if (Rf_inherits(slotExpr, "dgCMatrix")) {
    CscSlots csc = parseCscMatrix(slotExpr, data.numObservations);
    data.numPredictors = csc.numColumns;
    data.xIsSparse = true;
    // the all-CSC map: column j is CSC column j, the engine's ~j. Spelling the
    // identity out lets the bare matrix take the same mapped path a container
    // does, so nothing keys a rule on "no map means dense"
    data.columnSources.resize(data.numPredictors);
    for (size_t j = 0; j < data.numPredictors; ++j)
      data.columnSources[j] = ~static_cast<std::int32_t>(j);
    data.predictors.columnSources = data.columnSources.data();
    data.predictors.cscColumnPointers = csc.pointers;
    data.predictors.cscRowIndices = csc.rows;
    data.predictors.cscValues = csc.values;
  } else if (Rf_inherits(slotExpr, "dbartsMixedMatrix")) {
    SEXP denseExpr = PROTECT(rc_getListElement(slotExpr, "dense"));
    SEXP sparseExpr = PROTECT(rc_getListElement(slotExpr, "sparse"));
    SEXP mapExpr = PROTECT(rc_getListElement(slotExpr, "map"));
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
      data.predictors.denseValues = data.denseAssembly.data();
    } else {
      // the mixed flavor: a per-column dense list (factors carrying their
      // integer codes, or NULL for no dense columns), a dgCMatrix, and a
      // 1-based map - positive k names dense column k, negative -k sparse
      // column k, the engine's ~(k - 1). Assemble the transient block - the
      // exact doubles the retained cbind held - which the mapped build copies
      // into the store, so it need only outlive the ensuing build call.
      if (!Rf_inherits(sparseExpr, "dgCMatrix"))
        Rf_error("malformed mixed predictor container");
      if (!Rf_isNull(denseExpr) && TYPEOF(denseExpr) != VECSXP)
        Rf_error("malformed mixed predictor container");
      CscSlots csc = parseCscMatrix(sparseExpr, data.numObservations);

      data.numPredictors = rc_getLength(mapExpr);
      const int* map = INTEGER(mapExpr);
      parseMixedContainerBlock(
        denseExpr, map, data.numPredictors, data.numObservations,
        csc.numColumns, data.denseAssembly, data.predictors.denseValues,
        data.columnSources,
        "number of rows of 'x' must equal length of 'y'",
        "malformed mixed predictor container");
      data.xIsMixed = true;
      data.predictors.columnSources = data.columnSources.data();
      data.predictors.cscColumnPointers = csc.pointers;
      data.predictors.cscRowIndices = csc.rows;
      data.predictors.cscValues = csc.values;

      // the per-sparse-column reference metadata a CSC-backed categorical
      // column needs, borrowed from the container (protected via dataExpr)
      requireCscReferenceMeta(slotExpr, csc.numColumns, data.cscReferenceMeta,
                              data.cscCategoryCountMeta,
                              "malformed mixed predictor container");
      data.numSparseColumns = csc.numColumns;
    }
    UNPROTECT(3);
  } else {
    if (!Rf_isReal(slotExpr)) Rf_error("'x' must be numeric");
    rc_assertDimConstraints(slotExpr, "dimensions of x", RC_LENGTH | RC_EQ,
                            rc_asRLength(2), RC_VALUE | RC_EQ,
                            static_cast<int>(data.numObservations), RC_END);
    int* dims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
    data.predictors.denseValues = REAL(slotExpr);
    data.numPredictors = static_cast<size_t>(dims[1]);
  }
  data.predictors.numRows = data.numObservations;
  data.predictors.numColumns = data.numPredictors;

  REPROTECT_SLOT(slotExpr, dataExpr, "varTypes", slotIndex);
  rc_assertIntConstraints(slotExpr, "variable types", RC_LENGTH | RC_EQ,
                          rc_asRLength(data.numPredictors), RC_END);
  int* variableTypes = INTEGER(slotExpr);
  data.columnTypes.assign(data.numPredictors, bartcore::ColumnType::ordinal);
  for (size_t j = 0; j < data.numPredictors; ++j) {
    if (variableTypes[j] == 0) continue; // 0 encodes ordinal
    data.columnTypes[j] = bartcore::ColumnType::categorical;
    data.anyCategorical = true;
  }
  data.predictors.columnTypes =
    data.anyCategorical ? data.columnTypes.data() : NULL;
  // keyed on the parse-time flavor, never on the map: a bare dgCMatrix carries
  // an all-CSC map exactly as a container's sparse columns do, and only the
  // former lacks the reference metadata a categorical column needs
  if (data.xIsSparse && data.anyCategorical)
    Rf_error("sparse predictor matrices must be entirely ordinal");
  if (data.xIsMixed && data.anyCategorical) {
    // a CSC-backed categorical column reaches the engine only with its level
    // count K and reference code, resolved here per predictor from the
    // container's per-sparse-column metadata - which the parse above required
    // of every container carrying a sparse block, so it is present here
    bool anyCscCategorical = false;
    for (size_t j = 0; j < data.numPredictors; ++j)
      if (data.columnTypes[j] == bartcore::ColumnType::categorical &&
          data.columnSources[j] < 0) {
        anyCscCategorical = true;
        break;
      }
    if (anyCscCategorical) {
      resolveCscCategoricalReferences(
        data.columnTypes.data(), data.columnSources.data(),
        data.numPredictors, data.cscReferenceMeta, data.cscCategoryCountMeta,
        data.numSparseColumns, data.cscReferenceCodes,
        &data.categoryCounts, "malformed mixed predictor container",
        "sparse categorical predictor columns require reference "
        "metadata");
      data.predictors.referenceCodes = data.cscReferenceCodes.data();
    }
  }
  // the dense-backed categorical columns' declared counts, on top of whatever
  // the CSC resolution above settled (which stays authoritative for the
  // columns it owns - a container declares its own K)
  if (data.anyCategorical)
    readDeclaredCategoryCounts(data, factorLevelsExpr);
  data.predictors.categoryCounts =
    data.categoryCounts.empty() ? NULL : data.categoryCounts.data();

  REPROTECT_SLOT(slotExpr, dataExpr, "x.test", slotIndex);
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.numTestObservations = 0;
  } else if (Rf_inherits(slotExpr, "dbartsMixedMatrix")) {
    // the mixed test container parses against the training cut grid; the store
    // copies its raw, so the parsed sources ride in data (freed on an error
    // jump with it) rather than being pinned by the holder
    parseTestContainer(data.testContainer, slotExpr, data.numPredictors,
                       data.columnTypes.data());
    data.testIsMixed = true;
    data.testPredictors = data.testContainer.view;
    data.numTestObservations = data.testPredictors.numRows;
  } else {
    if (!Rf_isReal(slotExpr)) Rf_error("'x.test' must be numeric");
    rc_assertDimConstraints(slotExpr, "dimensions of x.test",
                            RC_LENGTH | RC_EQ, rc_asRLength(2), RC_NA,
                            RC_VALUE | RC_EQ,
                            static_cast<int>(data.numPredictors), RC_END);
    int* testDims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
    data.numTestObservations = static_cast<size_t>(testDims[0]);
    data.testPredictors = bartcore::densePredictorSource(
      REAL(slotExpr), data.numTestObservations, data.numPredictors);
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

  // the per-forest amplitude bases ride the data object beside the weights they
  // mirror; absent (the usual case) leaves the list empty and every
  // single-forest path unchanged
  REPROTECT_SLOT(slotExpr, dataExpr, "bases", slotIndex);
  readForestBases(rc_isS4Null(slotExpr) ? R_NilValue : slotExpr, data);

  REPROTECT_SLOT(slotExpr, dataExpr, "sigma", slotIndex);
  data.sigmaEstimate = rc_getDouble(
    slotExpr, "sigma estimate", RC_LENGTH | RC_EQ, rc_asRLength(1),
    RC_NA | RC_YES, RC_VALUE | RC_GT, 0.0, RC_END);

  REPROTECT_SLOT(slotExpr, dataExpr, "n.cuts", slotIndex);
  rc_assertIntConstraints(slotExpr, "maximum number of cuts",
                          RC_LENGTH | RC_EQ,
                          rc_asRLength(data.numPredictors), RC_END);
  int* maxNumCuts = INTEGER(slotExpr);
  data.maxNumCuts.resize(data.numPredictors);
  for (size_t j = 0; j < data.numPredictors; ++j)
    data.maxNumCuts[j] = static_cast<uint32_t>(maxNumCuts[j]);

  UNPROTECT(2);  // the slot index and the factor level tables
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
                model.changeProbability - 1.0) >= sumToOneTolerance)
    Rf_error("rule proposal probabilities must sum to 1.0");

  REPROTECT_SLOT(slotExpr, modelExpr, "p.birth", slotIndex);
  model.birthProbability = rc_getDouble(
    slotExpr, "probability of birth in birth/death rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  REPROTECT_SLOT(slotExpr, modelExpr, "node.scale", slotIndex);
  model.nodeScale = rc_getDouble(
    slotExpr, "scale of node prior", RC_LENGTH | RC_EQ, rc_asRLength(1),
    RC_VALUE | RC_GT, 0.0, RC_END);

  // the named calibration, response units. NA is the unnamed default, which
  // the engine reads off the same non-finite test that lets a NaN through;
  // anything else must be usable as a divisor, so infinity is refused here
  // alongside the R validity method's own check.
  REPROTECT_SLOT(slotExpr, modelExpr, "prior.scale", slotIndex);
  model.priorScale = rc_getDouble(
    slotExpr, "named prior scale", RC_LENGTH | RC_EQ, rc_asRLength(1),
    RC_NA | RC_YES, RC_VALUE | RC_GT, 0.0, RC_END);
  if (!ISNAN(model.priorScale) && !std::isfinite(model.priorScale))
    Rf_error("named prior scale must be NA or a positive finite number");

  // monotone (mBART) directions ride the model as a per-predictor integer
  // attribute the R surface resolves; absent or all-zero keeps the
  // unconstrained constant leaf
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

  // interaction constraint: two model attributes the R surface resolves - a
  // scalar max order and a 2 x k integer matrix of 0-based forbidden pairs
  // (column-major, so INTEGER() is the flat (a, b) pair stream the engine
  // reads). Absent leaves the availability path byte-for-byte unchanged.
  REPROTECT_SLOT(slotExpr, modelExpr, "interaction.max.order", slotIndex);
  if (!Rf_isNull(slotExpr) && rc_getLength(slotExpr) > 0)
    model.interactionMaxOrder = static_cast<size_t>(rc_getInt(
      slotExpr, "interaction max order", RC_LENGTH | RC_EQ, rc_asRLength(1),
      RC_VALUE | RC_GEQ, 0, RC_END));

  REPROTECT_SLOT(slotExpr, modelExpr, "interaction.forbidden", slotIndex);
  if (!Rf_isNull(slotExpr) && rc_getLength(slotExpr) > 0) {
    if (!Rf_isInteger(slotExpr))
      Rf_error("interaction forbidden pairs must be resolved integer indices");
    R_xlen_t numEntries = rc_getLength(slotExpr);
    if (numEntries % 2 != 0)
      Rf_error("interaction forbidden pairs must come in (a, b) pairs");
    const int* forbidden = INTEGER(slotExpr);
    model.interactionForbiddenPairs.resize(static_cast<size_t>(numEntries));
    for (R_xlen_t j = 0; j < numEntries; ++j) {
      if (forbidden[j] < 0 ||
          static_cast<size_t>(forbidden[j]) >= numPredictors)
        Rf_error("interaction forbidden pair column out of range");
      model.interactionForbiddenPairs[static_cast<size_t>(j)] =
        static_cast<size_t>(forbidden[j]);
    }
  }

  // block-additive constraint (per-tree column grouping): two resolved model attributes - a
  // per-column 0-based group index (length numPredictors; negative = no block)
  // and a per-group tree capacity. Absent leaves every tree unrestricted,
  // byte-for-byte unchanged. The R surface validates the total disjoint partition
  // and the capacity sum; here we range-check defensively.
  REPROTECT_SLOT(slotExpr, modelExpr, "block.of.column", slotIndex);
  if (!Rf_isNull(slotExpr) && rc_getLength(slotExpr) > 0) {
    if (!Rf_isInteger(slotExpr))
      Rf_error("block column groups must be resolved integer indices");
    if (static_cast<size_t>(rc_getLength(slotExpr)) != numPredictors)
      Rf_error("block column groups must have one entry per predictor");
    const int* groups = INTEGER(slotExpr);
    model.blockOfColumn.assign(groups, groups + numPredictors);
  }

  REPROTECT_SLOT(slotExpr, modelExpr, "block.tree.counts", slotIndex);
  if (!Rf_isNull(slotExpr) && rc_getLength(slotExpr) > 0) {
    if (!Rf_isInteger(slotExpr))
      Rf_error("block tree counts must be resolved integers");
    R_xlen_t numBlocks = rc_getLength(slotExpr);
    const int* counts = INTEGER(slotExpr);
    model.blockTreeCounts.resize(static_cast<size_t>(numBlocks));
    for (R_xlen_t g = 0; g < numBlocks; ++g) {
      if (counts[g] <= 0) Rf_error("block tree counts must be positive");
      model.blockTreeCounts[static_cast<size_t>(g)] =
        static_cast<size_t>(counts[g]);
    }
    for (int32_t group : model.blockOfColumn)
      if (group >= static_cast<int32_t>(numBlocks))
        Rf_error("block column group index out of range");
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
    if (std::fabs(totalProbability - 1.0) >= sumToOneTolerance)
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

  // Student-t residual df: an optional numeric model attribute (resid.df)
  // the R surface attaches; absent it stays NaN, the Gaussian law. Read raw
  // here - the family cross-check and sign policy live in
  // parseSamplerSpecification, once the family is known.
  SEXP residDfExpr = Rf_getAttrib(modelExpr, Rf_install("resid.df"));
  if (Rf_isReal(residDfExpr) && Rf_xlength(residDfExpr) == 1)
    model.residualDf = REAL(residDfExpr)[0];

  UNPROTECT(3);
}

// The creation-time verbose summary printed under control.verbose; the
// quantile and scale lines reduce at creation to expressions over the raw
// prior scale and the response range, used here directly.
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
      // the fixed residual variance prior's summary line
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
      // numCuts is unsigned and every index below is taken relative to the
      // last cut, so a column carrying none prints nothing rather than wraps
      if (store.numCuts[j] == 0) {
        ext_printf("\n");
        continue;
      }

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
// response: "" or "probit" give the standard probit latents, "logistic" the
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
  // ordinal is the K-level categorical shape: the cumulative probit is the
  // only family defined on it
  if (control.numOrdinalCategories >= 2) {
    if (familyName[0] == '\0' || std::strcmp(familyName, "ordinal") == 0)
      return bartcore::ResponseFamily::ordinal;
    Rf_error("an ordinal (K-level) response supports only family \"ordinal\"");
  }
  // nbinom is the count shape: the negative binomial is the only family
  // defined on it
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

// Column j's dense slice of a parsed view, when a dense source serves it: the
// plain matrix's own column, or the block column the source map names. Null
// for a CSC-backed column (coded by its container, nothing contiguous to scan)
// and for a view with no values at all.
const double* rawViewColumn(const bartcore::PredictorSource& view, size_t j) {
  std::int32_t source = view.sourceOf(j);
  if (source < 0 || view.denseValues == NULL) return NULL;
  return view.denseValues + static_cast<size_t>(source) * view.numRows;
}

// The two refusal texts the categorical entrances share - the training side
// reporting representability, the test side membership in the training levels.
// The dense entrances have always spelled them this way and the CSC views
// reuse them, so a refusal reads the same wherever it fires.
const char* const categoricalTrainingMessage =
  "categorical predictors must hold integer category codes in [0, 65535)";
const char* const categoricalTestMessage =
  "categorical test predictors must hold existing category codes";

// The stored codes of a CSC-backed column of a parsed container, found through
// its engine-convention source (the complement of a sparse column's position).
struct ParsedCscCodes {
  const double* values = NULL;
  size_t numValues = 0;
};

ParsedCscCodes parsedCscCodes(const int* columnPointers, const double* values,
                              std::int32_t source) {
  size_t k = static_cast<size_t>(~source);
  size_t begin = static_cast<size_t>(columnPointers[k]);
  return { values + begin,
           static_cast<size_t>(columnPointers[k + 1]) - begin };
}

// Refuse any code outside [0, bound): a code must be integral, and only the
// reserved missing value is exempt. Cold path shared by every categorical
// entrance, so it takes the caller's refusal text rather than deciding one.
void refuseInvalidCategoryCodes(const double* values, size_t numValues,
                                double bound, const char* message) {
  for (size_t i = 0; i < numValues; ++i) {
    double value = values[i];
    if (bartcore::isNA(value)) continue;  // the reserved missing category
    if (value < 0.0 || value >= bound || value != std::floor(value))
      Rf_error("%s", message);
  }
}

// Column j's host-declared level count, or 0 where nothing declares one.
double declaredCategoryCount(const ParsedData& data, size_t j) {
  return data.categoryCounts.empty()
    ? 0.0 : static_cast<double>(data.categoryCounts[j]);
}

// The category count column j's test codes must fall under, keyed on the view
// the TRAINING side presents rather than on its storage kind: the K its host
// declared - the container's for a CSC-backed column, the factor.levels table's
// for a dense one - but never below the max + 1 its own codes reach. Mirrors
// buildCutsForColumn (data.hpp) exactly, so the bound is the numCuts the store
// is about to hold.
double trainingCategoryBound(const ParsedData& data, size_t j) {
  double declared = declaredCategoryCount(data, j);
  if (data.predictors.sourceOf(j) < 0) return declared;
  const double* column = rawViewColumn(data.predictors, j);
  double maxValue = -1.0;
  for (size_t i = 0; i < data.numObservations; ++i) {
    double value = column[i];
    if (!bartcore::isNA(value) && value > maxValue) maxValue = value;
  }
  double inferred = maxValue < 0.0 ? 0.0 : maxValue + 1.0;
  return declared > inferred ? declared : inferred;
}

// Bound every categorical code the creation parse is about to ingest. Runs
// before the store exists, so the training-side count is reconstructed rather
// than read off numCuts: a CSC-backed training column's stored codes must lie
// in its declared K, and a dense one's in its own declared K where its host
// supplied a level table, else need only be representable (the count is then
// inferred from them). Each test view is then bounded by that count, whatever
// backs it - the x.test matrix, a container's dense slice, a container's CSC
// slice, or the reference code its implicit rows read. An unbounded code would
// mis-bin, shift past a tree's category mask, or over-read a pooled bitmap.
void validateCategoricalPredictors(const ParsedData& data) {
  for (size_t j = 0; j < data.numPredictors; ++j) {
    if (data.columnTypes[j] != bartcore::ColumnType::categorical) continue;
    if (data.predictors.sourceOf(j) < 0) {
      // a CSC-backed column stores only its non-reference codes; those must lie
      // in the K its container declared - the count the store will take - as
      // parseData already required of its reference code
      ParsedCscCodes stored = parsedCscCodes(data.predictors.cscColumnPointers,
                                             data.predictors.cscValues,
                                             data.predictors.sourceOf(j));
      refuseInvalidCategoryCodes(stored.values, stored.numValues,
                                 declaredCategoryCount(data, j),
                                 categoricalTrainingMessage);
      continue;
    }
    // a declared level table bounds a dense training column outright; without
    // one the column's own codes fix the count, so they need only be
    // representable
    double declared = declaredCategoryCount(data, j);
    refuseInvalidCategoryCodes(
      rawViewColumn(data.predictors, j), data.numObservations,
      declared > 0.0 ? declared : static_cast<double>(bartcore::maxCategories),
      categoricalTrainingMessage);
  }
  if (!data.anyCategorical || data.numTestObservations == 0) return;
  for (size_t j = 0; j < data.numPredictors; ++j) {
    if (data.columnTypes[j] != bartcore::ColumnType::categorical) continue;
    double bound = trainingCategoryBound(data, j);
    if (data.testPredictors.sourceOf(j) < 0) {
      ParsedCscCodes stored = parsedCscCodes(
        data.testPredictors.cscColumnPointers, data.testPredictors.cscValues,
        data.testPredictors.sourceOf(j));
      refuseInvalidCategoryCodes(stored.values, stored.numValues, bound,
                                 categoricalTestMessage);
      double reference =
        static_cast<double>(data.testPredictors.referenceCodeOf(j));
      refuseInvalidCategoryCodes(&reference, 1, bound, categoricalTestMessage);
      continue;
    }
    const double* testColumn = rawViewColumn(data.testPredictors, j);
    if (testColumn == NULL) continue;  // no test view of this column
    refuseInvalidCategoryCodes(testColumn, data.numTestObservations, bound,
                               categoricalTestMessage);
  }
}

bartcore::SamplerOptions optionsFromParsed(const ParsedControl& control,
                                           const ParsedModel& model,
                                           const ParsedData& data,
                                           SEXP modelExpr, bool sigmaIsFixed) {
  bartcore::SamplerOptions options;
  options.numTrees = control.numTrees;
  options.sigmaIsFixed = sigmaIsFixed;
  // the gaussian construction path reads this: finite selects Student-t
  // errors, NaN keeps the Gaussian law
  options.residualDf = model.residualDf;
  // ordinal construction reads this: K >= 2 selects OrdinalResponse with a
  // K-1 threshold vector, 0 leaves a non-ordinal response
  options.numCategories = control.numOrdinalCategories;
  // nbinom construction reads this: a positive value fixes the dispersion r,
  // a non-positive value estimates it on the grid; NaN (a non-count
  // response) is ignored
  options.dispersion = control.dispersion;
  options.k = model.k;
  options.nodeScale = model.nodeScale;
  // NA_REAL is a NaN, so the engine's isfinite test reads "unnamed" from it
  options.priorScale = model.priorScale;
  options.base = model.base;
  options.power = model.power;
  options.birthOrDeathProbability = model.birthOrDeathProbability;
  options.swapProbability = model.swapProbability;
  options.changeProbability = model.changeProbability;
  options.birthProbability = model.birthProbability;
  options.maxNumCutsPerVariable = data.maxNumCuts.data(); // copied at build
  options.useQuantiles = control.useQuantiles;
  // one borrowed view carries the storage and the typing channel (types,
  // declared level counts, CSC reference codes); consumed at build
  options.predictors = data.predictors;
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

  // per-forest interaction constraint (single-forest path): the resolved max
  // order and flat 0-based forbidden-pair stream, consumed at construction
  options.interactionMaxOrder = model.interactionMaxOrder;
  options.interactionForbiddenPairs = model.interactionForbiddenPairs.empty()
    ? NULL : model.interactionForbiddenPairs.data();
  options.interactionNumForbiddenPairs =
    model.interactionForbiddenPairs.size() / 2;

  // per-forest block-additive constraint (per-tree column grouping, single-forest path): the
  // per-column group and per-group tree capacity, consumed at construction. The
  // capacity must sum to the tree count (the R surface guarantees it; a defensive
  // backstop here since a mismatch would misassign trees to groups).
  if (!model.blockTreeCounts.empty()) {
    size_t sum = 0;
    for (size_t c : model.blockTreeCounts) sum += c;
    if (sum != options.numTrees)
      Rf_error("block tree counts must sum to the number of trees");
    options.numBlocks = model.blockTreeCounts.size();
    options.blockOfColumn = model.blockOfColumn.data();
    options.blockTreeCounts = model.blockTreeCounts.data();
  }

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
  // the gaussian constant-leaf gate is enforced at each sampler-construction
  // site (this parse cannot see the family/leaf); the factory mints the fp32
  // instantiation only when this flag survives that gate
  options.fp32Residual = control.fp32Residual;

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

  SEXP indicesExpr = rc_getListElement(groupsExpr, "indices");
  SEXP numGroupsExpr = rc_getListElement(groupsExpr, "n.groups");
  SEXP priorExpr = rc_getListElement(groupsExpr, "prior");
  SEXP scaleExpr = rc_getListElement(groupsExpr, "rel.scale");
  SEXP stepsExpr = rc_getListElement(groupsExpr, "n.steps");
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
// the attribute from a Surv or (time, status) ingest.
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

// The heteroscedastic variance forest arrives on an internal control
// attribute `bartcore.variance`: a list carrying the
// tree count, the tree-structure prior (base/power), and optional 1-based
// column indices (the `variance = ~ subset` selector; null spans all mean
// predictors). Absent leaves the fit homoscedastic. The factory refuses the
// combination for non-gaussian or non-constant-leaf models; the R surface
// mirrors that refusal. columns outlives the options borrow (captured by value
// in createHolder's lambda).
void applyVarianceAttributes(SEXP controlExpr, size_t numPredictors,
                             bartcore::SamplerOptions& options,
                             std::vector<std::size_t>& columns) {
  SEXP varExpr = Rf_getAttrib(controlExpr, Rf_install("bartcore.variance"));
  if (Rf_isNull(varExpr)) return;

  SEXP nTreesExpr = rc_getListElement(varExpr, "n.trees");
  SEXP baseExpr = rc_getListElement(varExpr, "base");
  SEXP powerExpr = rc_getListElement(varExpr, "power");
  SEXP columnsExpr = rc_getListElement(varExpr, "columns");
  if (!Rf_isInteger(nTreesExpr) || Rf_xlength(nTreesExpr) != 1 ||
      !Rf_isReal(baseExpr) || Rf_xlength(baseExpr) != 1 ||
      !Rf_isReal(powerExpr) || Rf_xlength(powerExpr) != 1)
    Rf_error("malformed variance forest specification");

  int numTrees = INTEGER(nTreesExpr)[0];
  if (numTrees < 1) Rf_error("variance n.trees must be a positive integer");
  double base = REAL(baseExpr)[0], power = REAL(powerExpr)[0];
  if (!(base > 0.0 && base < 1.0)) Rf_error("variance base must be in (0, 1)");
  if (!(power > 0.0)) Rf_error("variance power must be positive");
  options.numVarianceTrees = static_cast<size_t>(numTrees);
  options.varianceBase = base;
  options.variancePower = power;

  if (!Rf_isNull(columnsExpr)) {
    if (!Rf_isInteger(columnsExpr))
      Rf_error("variance columns must be resolved integer indices");
    R_xlen_t numColumns = Rf_xlength(columnsExpr);
    columns.resize(static_cast<size_t>(numColumns));
    for (R_xlen_t j = 0; j < numColumns; ++j) {
      int column = INTEGER(columnsExpr)[j];
      if (column < 1 || static_cast<size_t>(column) > numPredictors)
        Rf_error("variance column index out of range");
      columns[static_cast<size_t>(j)] = static_cast<size_t>(column - 1);
    }
    options.varianceForestColumns = columns.data();
    options.numVarianceForestColumns = columns.size();
  }
}

// Read a resolved interactions() list (the resolveInteractions output: a
// max.order scalar and a 2 x k integer matrix of 0-based forbidden pairs) into
// one forest's spec. A NULL list leaves the forest unconstrained. The pair
// stream is appended to `storage` (which must outlive the sampler build and be
// distinct per forest, so the borrowed pointer stays stable) and borrowed by
// the spec.
void applyForestInteractions(SEXP listExpr, size_t numPredictors,
                             bartcore::ForestStructureSpec& spec,
                             std::vector<size_t>& storage) {
  if (Rf_isNull(listExpr)) return;
  SEXP maxOrderExpr = rc_getListElement(listExpr, "max.order");
  if (!Rf_isNull(maxOrderExpr) && Rf_xlength(maxOrderExpr) == 1) {
    int order = Rf_asInteger(maxOrderExpr);
    if (order > 0) spec.interactionMaxOrder = static_cast<size_t>(order);
  }
  SEXP forbiddenExpr = rc_getListElement(listExpr, "forbidden");
  if (!Rf_isNull(forbiddenExpr) && Rf_xlength(forbiddenExpr) > 0) {
    if (!Rf_isInteger(forbiddenExpr))
      Rf_error("forest interaction forbidden pairs must be resolved integers");
    R_xlen_t numEntries = Rf_xlength(forbiddenExpr);
    if (numEntries % 2 != 0)
      Rf_error("forest interaction forbidden pairs must come in (a, b) pairs");
    const int* forbidden = INTEGER(forbiddenExpr);
    storage.resize(static_cast<size_t>(numEntries));
    for (R_xlen_t j = 0; j < numEntries; ++j) {
      if (forbidden[j] < 0 || static_cast<size_t>(forbidden[j]) >= numPredictors)
        Rf_error("forest interaction forbidden pair column out of range");
      storage[static_cast<size_t>(j)] = static_cast<size_t>(forbidden[j]);
    }
    spec.interactionForbiddenPairs = storage.data();
    spec.interactionNumForbiddenPairs = storage.size() / 2;
  }
}

// Read a resolved blocks() list (the resolveBlocks output: a per-column 0-based
// group index and a per-group tree capacity) into one forest's spec. A NULL
// list leaves the forest unrestricted. The group and capacity buffers are stored
// in `groupStore` / `countStore` (distinct per forest, outliving the sampler
// build) and borrowed by the spec. The capacity must sum to the forest's tree
// count (the R surface guarantees it; a defensive backstop here).
void applyForestBlocks(SEXP listExpr, size_t numPredictors, size_t numTrees,
                       bartcore::ForestStructureSpec& spec,
                       std::vector<std::int32_t>& groupStore,
                       std::vector<size_t>& countStore) {
  if (Rf_isNull(listExpr)) return;
  SEXP groupExpr = rc_getListElement(listExpr, "block.of.column");
  SEXP countExpr = rc_getListElement(listExpr, "block.tree.counts");
  if (Rf_isNull(groupExpr) || Rf_isNull(countExpr)) return;
  if (!Rf_isInteger(groupExpr) || !Rf_isInteger(countExpr))
    Rf_error("forest block spec must hold resolved integers");
  if (static_cast<size_t>(Rf_xlength(groupExpr)) != numPredictors)
    Rf_error("forest block column groups must have one entry per predictor");
  R_xlen_t numBlocks = Rf_xlength(countExpr);
  const int* counts = INTEGER(countExpr);
  size_t sum = 0;
  countStore.resize(static_cast<size_t>(numBlocks));
  for (R_xlen_t g = 0; g < numBlocks; ++g) {
    if (counts[g] <= 0) Rf_error("forest block tree counts must be positive");
    countStore[static_cast<size_t>(g)] = static_cast<size_t>(counts[g]);
    sum += static_cast<size_t>(counts[g]);
  }
  if (sum != numTrees)
    Rf_error("forest block tree counts must sum to the forest's tree count");
  const int* groups = INTEGER(groupExpr);
  groupStore.assign(groups, groups + numPredictors);
  for (int32_t group : groupStore)
    if (group >= static_cast<int32_t>(numBlocks))
      Rf_error("forest block column group index out of range");
  spec.numBlocks = static_cast<size_t>(numBlocks);
  spec.blockOfColumn = groupStore.data();
  spec.blockTreeCounts = countStore.data();
}

/// Backing storage for the pointers an AmplitudeSpec borrows: each forest's
/// column subset, its forbidden-interaction pair stream, and its block
/// partition. One instance per sampler build, outliving it. RAGGED - one entry
/// per forest of the spec, so nothing is keyed on bcf's two.
struct AmplitudeSpecStorage {
  std::vector<std::vector<size_t>> columns, pairs, blockCounts;
  std::vector<std::vector<std::int32_t>> blockGroups;
  // the per-forest bases, transposed from R's column-major into the engine's
  // row-major contract; borrowed by the spec until the chains copy them
  std::vector<std::vector<double>> bases;
};

// One element of a per-forest list, or R_NilValue when the list is null or
// too short. The lists the control attribute carries are all K-length and
// parallel, but one that stops short leaves the rest at their defaults rather
// than erroring twice for the same omission.
SEXP forestListElement(SEXP listExpr, size_t f) {
  if (Rf_isNull(listExpr) || !Rf_isNewList(listExpr)) return R_NilValue;
  if (static_cast<R_xlen_t>(f) >= Rf_xlength(listExpr)) return R_NilValue;
  return VECTOR_ELT(listExpr, static_cast<R_xlen_t>(f));
}

// Fill the K-length forest spec from the resolved R objects. The FIRST forest
// takes its tree count and structure prior from the host model, since the
// fit's own control@n.trees and tree prior are its declaration; every forest
// takes the rest from its own length-8 params vector - tree count, base,
// power, the node-scale factor and divisor the calibration map reads, the
// amplitude prior's variance and half-Cauchy scale, and the amplitude update
// flag. The per-forest vars/interactions()/blocks() lists are optional and
// parallel. Shared by the internal entry, which passes these as arguments, and
// the public creation path, which reads them off a control attribute - so the
// two build the same sampler.
//
// A forest travels its amplitude's likelihood-invariant ASIS ridge exactly
// when its prior is a SCALE MIXTURE. That reproduces bcf - the a-move on, the
// b-move off - and it is the honest general rule: the fixed-variance ridge is
// a held-off door, whose cost is a GIG draw per sweep and whose acceptance
// gate has not been run.
void applyAmplitudeSpec(SEXP paramsExpr, SEXP varsExpr, SEXP interactionsExpr,
                        SEXP blocksExpr, const ParsedModel& model,
                        size_t numTrees, size_t numPredictors,
                        bartcore::AmplitudeSpec& spec,
                        AmplitudeSpecStorage& storage) {
  if (Rf_isNull(paramsExpr) || !Rf_isNewList(paramsExpr) ||
      Rf_xlength(paramsExpr) < 2)
    Rf_error("forest parameters must be a list of at least two per-forest "
             "parameter vectors");
  size_t numForests = static_cast<size_t>(Rf_xlength(paramsExpr));
  spec.forests.assign(numForests, bartcore::ForestSpec{});
  storage.columns.resize(numForests);
  storage.pairs.resize(numForests);
  storage.blockCounts.resize(numForests);
  storage.blockGroups.resize(numForests);

  for (size_t f = 0; f < numForests; ++f) {
    SEXP forestParamsExpr = VECTOR_ELT(paramsExpr, static_cast<R_xlen_t>(f));
    if (!Rf_isReal(forestParamsExpr) || Rf_xlength(forestParamsExpr) != 8)
      Rf_error("forest parameters must be a length-8 numeric vector per forest");
    const double* params = REAL(forestParamsExpr);
    bartcore::ForestSpec& forest = spec.forests[f];
    // node scales come from the calibration map, not the host model
    forest.forest.numTrees =
      f == 0 ? numTrees : static_cast<size_t>(params[0]);
    forest.forest.base = f == 0 ? model.base : params[1];
    forest.forest.power = f == 0 ? model.power : params[2];
    forest.nodeScaleFactor = params[3];
    forest.nodeScaleDivisor = params[4];
    forest.amplitudePriorVariance = params[5];
    forest.amplitudePriorScale = params[6];
    forest.updateAmplitude = params[7] != 0.0;
    forest.ridge = forest.amplitudePriorScale > 0.0;

    // the forest's optional column restriction: 1-based indices resolved
    // R-side, or NULL for an unrestricted forest reading the full store.
    // Consumed at construction, so the buffer need only outlive the build.
    SEXP forestVarsExpr = forestListElement(varsExpr, f);
    if (!Rf_isNull(forestVarsExpr)) {
      if (!Rf_isInteger(forestVarsExpr))
        Rf_error("forest variables must be resolved integer column indices");
      R_xlen_t numColumns = Rf_xlength(forestVarsExpr);
      storage.columns[f].resize(static_cast<size_t>(numColumns));
      for (R_xlen_t j = 0; j < numColumns; ++j) {
        int column = INTEGER(forestVarsExpr)[j];
        if (column < 1 || static_cast<size_t>(column) > numPredictors)
          Rf_error("forest variable column out of range");
        storage.columns[f][static_cast<size_t>(j)] =
          static_cast<size_t>(column - 1);
      }
      forest.forest.columns = storage.columns[f].data();
      forest.forest.numColumns = storage.columns[f].size();
    }

    // independent per-forest interaction constraints (the calibrated-additivity
    // causal use) and block-additive partitions: each forest resolves its own
    // prior R-side over the columns it may split on, and the engine intersects
    // a restricted forest's block rows with its column mask at install
    applyForestInteractions(forestListElement(interactionsExpr, f),
                            numPredictors, forest.forest, storage.pairs[f]);
    applyForestBlocks(forestListElement(blocksExpr, f), numPredictors,
                      forest.forest.numTrees, forest.forest,
                      storage.blockGroups[f], storage.blockCounts[f]);
  }
}

// Transposes the data object's per-forest bases from R's COLUMN-major layout
// into the engine's ROW-major contract and borrows them onto the spec, which
// the chains copy at construction. A forest the list does not reach, or whose
// entry is null, keeps the dense all-ones column the engine synthesizes.
void applyForestBases(const ParsedData& data, bartcore::AmplitudeSpec& spec,
                      AmplitudeSpecStorage& storage) {
  size_t n = data.numObservations;
  storage.bases.resize(spec.forests.size());
  for (size_t f = 0; f < spec.forests.size() && f < data.bases.size(); ++f) {
    const double* values = data.bases[f];
    if (values == NULL) continue;
    size_t q = data.basisColumns[f];
    std::vector<double>& rowMajor = storage.bases[f];
    rowMajor.resize(n * q);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < q; ++j) {
        double value = values[j * n + i];
        if (!std::isfinite(value))
          Rf_error("a forest basis value is not finite");
        rowMajor[i * q + j] = value;
      }
    spec.forests[f].basis = rowMajor.data();
    spec.forests[f].numBasisColumns = q;
  }
}

// The multi-forest twin of applyVarianceAttributes: the forests' configuration
// arrives on the internal control attribute `bartcore.forests`, a list
// carrying the K-length params list and the parallel per-forest
// vars/interactions()/blocks() lists. Absent leaves the fit single-forest and
// returns false, so every existing creation path reads one attribute and moves
// on. The bases themselves ride the data object; createHolder cross-checks the
// two halves.
bool applyForestAttributes(SEXP controlExpr, const ParsedModel& model,
                           size_t numTrees, size_t numPredictors,
                           bartcore::AmplitudeSpec& spec,
                           AmplitudeSpecStorage& storage) {
  SEXP forestsExpr = Rf_getAttrib(controlExpr, Rf_install("bartcore.forests"));
  if (Rf_isNull(forestsExpr)) return false;
  applyAmplitudeSpec(rc_getListElement(forestsExpr, "params"),
                     rc_getListElement(forestsExpr, "vars"),
                     rc_getListElement(forestsExpr, "interactions"),
                     rc_getListElement(forestsExpr, "blocks"), model, numTrees,
                     numPredictors, spec, storage);
  return true;
}

// The families a K-forest coupling admits, and why each door is shut: null for
// gaussian, probit and logistic, the reason otherwise. Admission turns on the
// calibration map having a latent scale to state its node scales against, and
// on the family's own parameter block being shown to interleave with the
// amplitude block. The two creation routes ask this separately rather than
// sharing a gate, since only one runs the whole offender cascade below.
const char* refusedAmplitudeFamilyReason(bartcore::ResponseFamily family) {
  switch (family) {
  case bartcore::ResponseFamily::gaussian:
  case bartcore::ResponseFamily::probit:
  case bartcore::ResponseFamily::logistic:
    return NULL;
  case bartcore::ResponseFamily::aft:
    return "an AFT (survival) response: it draws sigma, which the coupling "
           "pins, and its censoring status reaches no K-forest creation path";
  case bartcore::ResponseFamily::ordinal:
    return "an ordinal response: its threshold block is not shown to interleave "
           "with the amplitude block";
  case bartcore::ResponseFamily::nbinom:
    return "a count (nbinom) response: its dispersion block is not shown to "
           "interleave with the amplitude block";
  }
  return "this response family";
}

// The node scale the R surface writes for a family that names none, mirroring
// defaultNodeScale() in R/model.R. The K-forest gate reads it so "a non-default
// node scale" means the same thing under every family: the literal 0.5 it used
// to compare against is GAUSSIAN'S, inlined while the coupling was gaussian.
double defaultNodeScale(bartcore::ResponseFamily family) {
  switch (family) {
  case bartcore::ResponseFamily::probit:
  case bartcore::ResponseFamily::ordinal:
    return 3.0;
  case bartcore::ResponseFamily::logistic:
  case bartcore::ResponseFamily::nbinom:
    return std::numbers::pi * std::sqrt(3.0);
  // gaussian and aft state the scale in response units, below. Every
  // enumerator is listed and there is no default arm: a family added without
  // one must fail the build rather than silently take gaussian's 0.5, and its
  // arm must be added to R's defaultNodeScale in the same change.
  case bartcore::ResponseFamily::gaussian:
  case bartcore::ResponseFamily::aft:
    break;
  }
  return 0.5;
}

// Every option the amplitude chain constructor does not read, refused rather
// than dropped in silence: the calibration map fixes every forest's leaf scale
// and k, buildSpecifiedForest takes no DART or split probabilities, the
// constant leaf is the single instantiation, the grouped decorator and the
// variance forest are built only by the single-forest constructor, the cut cap
// and the test surface are left undefined, and the gaussian response law is
// not the Student-t mixture. The R surface refuses the same list ahead of this
// backstop, which is what a direct dbarts.h consumer meets.
void refuseUnsupportedAmplitudeComposition(
    bartcore::ResponseFamily family, const ParsedModel& model,
    const ParsedData& data, const bartcore::SamplerOptions& options) {
  if (const char* refused = refusedAmplitudeFamilyReason(family))
    Rf_error("a treatment forest does not support %s", refused);
  if (!data.predictors.isDenseBlock())
    Rf_error("a treatment forest requires dense predictors");
  const char* offender = NULL;
  if (options.useDart) offender = "a DART tree prior";
  else if (options.splitProbabilities != NULL) offender = "split probabilities";
  else if (!model.monotoneDirections.empty()) offender = "monotone constraints";
  else if (options.numLeafCovariates != 0 || options.gpLeaves)
    offender = "a linear or Gaussian-process node prior";
  else if (model.updateK) offender = "a k hyperprior";
  else if (model.k != 2.0) offender = "a non-default k";
  else if (model.nodeScale != defaultNodeScale(family))
    offender = "a non-default node scale";
  // the node-scale gate above does not fire on a model that names its
  // calibration in response units instead, and the calibration map would drop
  // it in silence, so it is its own offender
  else if (std::isfinite(model.priorScale)) offender = "a named 'prior.scale'";
  else if (model.birthOrDeathProbability != 0.5 ||
           model.swapProbability != 0.1 || model.changeProbability != 0.4 ||
           model.birthProbability != 0.5)
    offender = "non-default proposal probabilities";
  else if (std::isfinite(model.residualDf)) offender = "Student-t residuals";
  else if (options.numGroups > 0) offender = "grouped random effects";
  else if (options.numVarianceTrees > 0) offender = "a variance forest";
  else if (options.fp32Residual) offender = "single-precision storage";
  else if (data.numTestObservations > 0) offender = "test predictors";
  if (offender == NULL) {
    // a uniform cut cap is the only one the amplitude chain keeps; it nulls the
    // per-variable vector outright
    for (size_t j = 1; j < data.numPredictors; ++j)
      if (data.maxNumCuts[j] != data.maxNumCuts[0]) {
        offender = "per-column cut counts";
        break;
      }
  }
  if (offender != NULL)
    Rf_error("a treatment forest does not support %s; drop it or fit a "
             "single-forest model", offender);
}

// A sampler created over a data handle holds no raw predictor values, so
// the raw-x mutation surface has nothing to work from; a CSC-built sampler
// holds only borrowed slices and refuses the same surface by design.
void refusePredictorMutation(const bartcore::SamplerBase& sampler,
                             const char* caller) {
  const bartcore::ColumnStore& data = sampler.data();
  if (data.acceptsNewRawPredictors()) return;
  if (data.builtFromCsc)
    Rf_error("%s: sparse predictors fix the design at creation; make a new "
             "sampler instead", caller);
  Rf_error("%s: requires a sampler that owns its predictors; data-handle "
           "views hold none", caller);
}

// The predictor-MUTATION surface (setPredictor, updatePredictor, and the
// per-observation sessions) re-quantizes columns from a new dense column, so
// unlike whole-data setData it runs on CSC/mixed stores too (they retain
// their slices; a mutated column repoints at engine-owned nonzeros). Cut
// installation and state restore re-quantize from those same retained
// slices, so they share this guard rather than carrying one of their own.
// Only a data-handle view, which keeps
// no raw source at all, is refused here; a per-observation caller additionally
// refuses a CSC-backed target column by name below.
void refuseMutationOnView(const bartcore::SamplerBase& sampler,
                          const char* caller) {
  if (!sampler.data().hasRequantizeSource())
    Rf_error("%s: requires a sampler that owns its predictors; data-handle "
             "views hold none", caller);
}

// storeSample adds the test offset to every reported channel AFTER the forests
// are blended, so on the one multi-forest model whose test blend IS defined
// (the multinomial softmax) it would shift the reported probabilities off the
// simplex rather than shift any latent. parseMultinomialData refuses a test
// offset at creation; refuse the post-creation install to match. Ordered after
// refuseUndefinedTestFits at every call site, so an amplitude coupling keeps
// its own message.
//
// The flat vector stays refused on the softmax coupling FOREVER, and truthfully
// - after the blend it leaves the simplex, and before it a common
// per-observation shift is the softmax's own null direction. What the coupling
// does carry is the per-category matrix, so the message names it rather than
// leaving a caller to conclude the test rows admit no offset at all. The
// generic wording stays verbatim for every other multi-forest coupling.
void refuseMultiForestTestOffset(const bartcore::SamplerBase& sampler,
                                 const char* caller) {
  bartcore::SamplerShape shape = sampler.shape();
  if (shape.numForests < 2) return;
  if (shape.supportsCountsMutation)
    Rf_error("%s: a flat test offset is added to every reported channel after "
             "the categories are blended, which would move the reported "
             "probabilities off the simplex, and before the blend a common "
             "per-observation shift is inert; write an nTest x K matrix "
             "through the category test offset channel", caller);
  Rf_error("%s: a multi-forest sampler adds a test offset after its forests "
           "are blended, which would move the reported values off the "
           "softmax scale; it carries no test offset", caller);
}

// Validates a category offset argument in place and returns its values,
// borrowed from R for the call's duration; a null expression yields null. The
// shared check behind every category-offset entrance, train, test and predict:
// anything non-null must be a real matrix of exactly n x K. The DIMENSIONS are
// checked, not the length, for the reason the counts entrance checks them - a
// transposed matrix carries exactly n*K entries and would shift every cell by
// another cell's offset. Every entry must be FINITE: an infinity propagates
// through the log-sum-exp margin into a NaN for every category at that row, and
// an NA is not a shift.
const double* validateCategoryOffset(SEXP offsetExpr, size_t n, size_t K,
                                     const char* caller) {
  if (Rf_isNull(offsetExpr)) return NULL;
  SEXP dimsExpr = Rf_getAttrib(offsetExpr, R_DimSymbol);
  if (!Rf_isReal(offsetExpr) || Rf_xlength(dimsExpr) != 2 ||
      static_cast<size_t>(INTEGER(dimsExpr)[0]) != n ||
      static_cast<size_t>(INTEGER(dimsExpr)[1]) != K)
    Rf_error("%s: requires a real matrix of %lu observations x %lu categories",
             caller, static_cast<unsigned long>(n),
             static_cast<unsigned long>(K));
  const double* src = REAL(offsetExpr);
  for (size_t j = 0; j < n * K; ++j)
    if (!R_finite(src[j]))
      Rf_error("%s: requires every category offset entry to be finite",
               caller);
  return src;
}

// The same into scratch, for the entrances that INSTALL one: the copy is the
// caller's to swap in - the combiner borrows what it is handed - so nothing is
// retained from R. A null expression clears (out empty).
void parseCategoryOffset(SEXP offsetExpr, size_t n, size_t K,
                         std::vector<double>& out, const char* caller) {
  out.clear();
  const double* src = validateCategoryOffset(offsetExpr, n, K, caller);
  if (src == NULL) return;
  out.assign(src, src + n * K);
}

// The test twin, whose row count is the CURRENT test store's: without test rows
// there is nothing for a per-test-row offset to describe, and accepting one
// would leave it silently unread. The engine sizes nothing off it, so this is
// the only place the shape is pinned to the test store.
void parseCategoryTestOffset(SEXP offsetExpr, size_t nTest, size_t K,
                             std::vector<double>& out, const char* caller) {
  out.clear();
  if (Rf_isNull(offsetExpr)) return;
  if (nTest == 0)
    Rf_error("%s: requires test data: it carries one row per test row",
             caller);
  parseCategoryOffset(offsetExpr, nTest, K, out, caller);
}

// Every test-side channel now carries its own per-category offset - the
// reported test blend through the resident nTest x K one, and predict through
// the per-call matrix it takes - so a sampler holding a train category offset
// has a test surface that can express one. What is NOT expressible is deriving
// one from the other: the test rows are other rows, so an unset test offset
// means zero there, and predict refuses rather than guess (both stated at their
// entrances below). Installing test predictors under a resident test offset is
// refused for the same reason - the offset describes rows that are being
// replaced.
void refuseStaleCategoryTestOffset(bool hasCategoryTestOffset,
                                   const char* caller) {
  if (!hasCategoryTestOffset) return;
  Rf_error("%s: this sampler carries an nTest x K category test offset, one "
           "row per CURRENT test row; clear it with a null category test "
           "offset before replacing those rows", caller);
}

// A data handle is a built column store (cuts + codes) shared by row-subset
// view samplers; internal, held by external
// pointer. The pointer's protection slot pins the data expression whose x the
// store borrows.
void dataHandleFinalizer(SEXP ptrExpr) {
  bartcore::ColumnStore* store =
    static_cast<bartcore::ColumnStore*>(R_ExternalPtrAddr(ptrExpr));
  if (store == NULL) return;
  delete store;
  R_ClearExternalPtr(ptrExpr);
}

bartcore::ColumnStore& dataHandleFromExpression(SEXP ptrExpr) {
  bartcore::ColumnStore* store =
    static_cast<bartcore::ColumnStore*>(R_ExternalPtrAddr(ptrExpr));
  if (store == NULL)
    Rf_error("data handle function called on NULL external pointer");
  return *store;
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
  bartcore_bridge::enforceBinaryWeightPolicy(family, data.weights,
                                            data.numObservations);
  // the response-support rule, stated once for creation and every mutation
  // conduit: the R surface checks it first, so this is a no-op there and the
  // real gate for the flat C API, which has no R layer ahead of it
  bartcore_bridge::validateResponseSupport(family,
                                           control.numOrdinalCategories, data.y,
                                           data.numObservations,
                                           "sampler creation");
  // a weighted truncated-latent draw is not a coherent likelihood; AFT v1
  // rejects weights
  if (family == bartcore::ResponseFamily::aft && data.weights != NULL)
    Rf_error("aft (survival) models do not support case weights");
  // ordinal refuses weights for probit's reason: a weighted truncated-normal
  // latent likelihood is not a coherent model
  if (family == bartcore::ResponseFamily::ordinal && data.weights != NULL)
    Rf_error("ordinal models do not support weights: a weighted truncated-"
             "normal latent likelihood is not a coherent model");
  // nbinom refuses weights in v1: the usual count "weight" is exposure,
  // which belongs in the offset as a log-exposure term, not in observation
  // replication
  if (family == bartcore::ResponseFamily::nbinom && data.weights != NULL)
    Rf_error("nbinom (count) models do not support weights: exposure belongs in "
             "the offset as a log-exposure term");
  // Student-t continuous errors: a finite resid.df selects the
  // scale-mixture error law - a positive value fixes the degrees of
  // freedom, 0 estimates them on the grid. Only the gaussian family carries
  // it; NaN (absent) keeps the Gaussian law.
  if (std::isfinite(model.residualDf)) {
    if (family != bartcore::ResponseFamily::gaussian)
      Rf_error("a continuous gaussian response is required for Student-t "
               "residuals ('resid.df')");
    if (model.residualDf < 0.0)
      Rf_error("'resid.df' must be positive to fix the degrees of freedom, "
               "or 0 to estimate them");
  }
  return family;
}

// The starting residual sd is a REQUIRED input of every creation path that
// hands one to the engine: it seeds sigma and calibrates the residual-variance
// prior's scale (initialSigma_ and sigmaSqPrior_.scale), so an unresolved
// NA_real_ makes both NaN and every draw after them NaN - silently, and past
// repair, since setSigma writes the value but not the poisoned prior scale.
// dbartsData's @sigma defaults to NA and only the R resolution layer fills it
// (estimateStartingSigma), so this is the gate for a caller that arrives with a
// raw specification - the flat C API, which has no layer ahead of it. Keyed on
// the family alone: gaussian and aft read sigmaEstimate at chain construction,
// buildVarianceForest included - it calibrates the variance forest's own
// scale leaf from the same value before pinning the scalar sigma at 1, so a
// heteroscedastic sampler needs a resolved @sigma too - while the
// fixed-unit-scale families (probit, logistic, ordinal, nbinom) never read
// the value and keep the NA their R path leaves them.
void requireResolvedSigmaEstimate(bartcore::ResponseFamily family,
                                  double sigmaEstimate) {
  if (family != bartcore::ResponseFamily::gaussian &&
      family != bartcore::ResponseFamily::aft)
    return;
  if (std::isfinite(sigmaEstimate) && sigmaEstimate > 0.0) return;
  Rf_error("sampler creation: a gaussian or aft sampler calibrates its "
           "residual variance prior from a starting estimate of sigma, which "
           "the data specification leaves unresolved (data@sigma is NA); "
           "supply a positive value, as dbarts() and dbartsSpec() do");
}

} // namespace

namespace bartcore_bridge {

// A whole-data or whole-model mutation (setData, setWeights, setModel)
// rebuilds or reprices only forest 0: applyNewData, the response/latent
// refresh, and Chain::setModel's prior installation all touch forests_[0]
// alone, so on a multi-forest sampler (BCF, and any future multi-forest model)
// the other forests would keep fits - or, for setModel, a leaf scale - against
// the old data or an uncalibrated prior. For setData the guard is also memory
// safety, not just staleness: every per-forest amplitude basis is COPIED at the
// width and observation count it was installed at, and combinedFits, the
// reparameterization and the amplitude conditional all index it per row, so a
// per-forest applyNewData that grew n would over-read every basis and a
// fixed-n one would silently re-pair the old bases with new rows - a lifted
// refusal must take the bases in the same call.
// Refuse it; a multi-forest sampler
// fixes its data and prior at creation, as grouped/sparse/aft samplers do.
// setForestBasis, the one supported multi-forest data swap, routes through the
// combiner and stays allowed; setResponse is opt-in and scale-pinned rather
// than refused, and carries its own condition in
// refuseMultiForestResponseMutation below.
// External linkage: the flat C API (C_interface.cpp) answers the same
// condition on its own setTestOffset entry, through the predicate below.
bool isMultiForest(const bartcore::SamplerBase& sampler) {
  return sampler.shape().numForests >= 2;
}

void refuseMultiForestMutation(const bartcore::SamplerBase& sampler,
                               const char* caller) {
  if (!isMultiForest(sampler)) return;
  Rf_error("%s: a multi-forest sampler fixes its data at creation; make a "
           "new sampler instead", caller);
}

// The response, offset and weight conduits are one rule under three names. The
// coupling must re-derive every per-forest residual AND precision from y and w
// each sweep, caching nothing across sweeps (supportsResponseMutation), and the
// response transform must stay where the per-forest leaf calibrations were
// stated against it: updateScale = TRUE would re-anchor it while both forests
// keep the old calibration, silently decalibrating them. NA is not FALSE here -
// only an explicit FALSE takes the permitted branch. The offset conduit is the
// response-side swap under a different pointer (for a gaussian response
// setOffset(yBuild - yNew, FALSE) re-maps through the pinned transform exactly
// as setResponse(yNew, FALSE) does), so it carries both conditions; a coupling
// that does not opt in has no offset semantics at all, since the K-forest
// softmax is invariant to a common per-observation shift. Weights carry no
// scale - setWeights never moves the transform, and BCF's scaledResponseSd is
// unweighted - so only the opt-in half applies to them.
// External linkage: the flat C API reuses this guard's SCALE arm on its own
// setResponse and setOffset entries, and answers the conduit arm - a fixed
// property of the sampler, so no argument would have worked - through the
// predicate below rather than by unwinding.
bool responseConduitIsFixed(const bartcore::SamplerShape& shape) {
  return shape.numForests >= 2 && !shape.supportsResponseMutation;
}

void refuseMultiForestResponseMutation(const bartcore::SamplerBase& sampler,
                                       const char* caller,
                                       ResponseConduit conduit,
                                       int updateScale) {
  bartcore::SamplerShape shape = sampler.shape();
  if (responseConduitIsFixed(shape)) {
    // A coupling that owns its response as a count matrix does not fix it at
    // creation - it swaps through its own entry, and an integer case weight is
    // exactly the row-wise count replication that matrix already expresses - so
    // naming the channel is the whole repair. The offset conduit names the
    // other half: a flat offset really is inert under a softmax, but an n x K
    // one is not, and it has its own entry too. Both this bridge and the flat C
    // API call this helper, so the hint lands on both surfaces, which is the
    // point: a guard that exists so two surfaces cannot state different rules
    // must not be made to state two.
    if (shape.supportsCountsMutation) {
      if (conduit == ResponseConduit::offset)
        Rf_error("%s: this sampler's offset is its n x K category matrix, "
                 "which a flat offset cannot express; write it through the "
                 "category offset channel", caller);
      // Weights get their own arm rather than falling through to the response
      // one: they are not refused because a flat vector cannot carry the
      // matrix, but on MODEL grounds, and a caller told to swap the counts
      // instead would not learn that an integer weight already IS that swap.
      if (conduit == ResponseConduit::weights)
        Rf_error("%s: this sampler's case weights are already expressed by its "
                 "n x K count matrix - an integer weight is row-wise count "
                 "replication - and a non-integer one has no exact "
                 "augmentation sampler", caller);
      Rf_error("%s: this sampler's response is its n x K count matrix; replace "
               "it through the counts channel", caller);
    }
    const char* fixed = "fixes its response at creation";
    if (conduit == ResponseConduit::offset) fixed = "carries no offset";
    else if (conduit == ResponseConduit::weights)
      fixed = "fixes its case weights at creation";
    Rf_error("%s: this multi-forest sampler %s; make a new sampler instead",
             caller, fixed);
  }
  if (shape.numForests >= 2 && conduit != ResponseConduit::weights &&
      updateScale != FALSE)
    Rf_error("%s: a multi-forest sampler supports %s swap only with "
             "updateScale = FALSE, which pins the response transform its "
             "per-forest leaf calibrations are stated against", caller,
             conduit == ResponseConduit::response ? "a response" : "an offset");
}

// The heteroscedastic analogue of the multi-forest scale pin above, and the
// fifth sigma door: a variance forest's scale leaf is calibrated once, at
// creation, against the response transform in force then, and no path re-states
// it. updateScale = TRUE re-anchors that transform under the calibration, so
// every s^2(x) the forest reports is measured on a scale the model no longer
// uses and the fit runs away while getSigmas() - which reads the pinned sigma,
// not the forest - shows nothing. There is no algebra that would rescue it
// short of recalibrating the scale leaf, so refuse; updateScale = FALSE pins
// the transform and is the supported heteroscedastic response swap. Weights
// carry no transform and never reach here.
// External linkage: the flat C API reuses this guard on its own setResponse and
// setOffset entries.
void refuseVarianceForestScaleUpdate(const bartcore::SamplerBase& sampler,
                                     const char* caller,
                                     ResponseConduit conduit, int updateScale) {
  if (conduit == ResponseConduit::weights || updateScale == FALSE) return;
  if (!sampler.shape().hasVarianceForest) return;
  Rf_error("%s: a heteroscedastic sampler's variance forest is calibrated "
           "against the response transform fixed at creation, so %s swap is "
           "supported only with updateScale = FALSE, which pins it", caller,
           conduit == ResponseConduit::response ? "a response" : "an offset");
}

// The grouped analogue of the two scale pins above, and the one place the
// random-intercept decorator is not scale-transparent: GroupedResponse holds b,
// tau and the tau prior scale on the BASE MODEL'S INTERNAL scale (its class
// comment, and the constructor's single division by sigmaScale), and its
// setResponse/setOffset delegate to the base and rebuild the working response
// without touching any of the three. At updateScale = TRUE a re-anchoring base
// recomputes its range and converts exactly sigma and the residual prior scale,
// so b and tau silently come to mean something else on the original scale while
// nothing reports the move - the same defect class as the two guards above.
// Keyed on the family for refusePinnedSigmaChange's reason and with its exact
// two-way shape: gaussian and aft are the families whose transform is derived
// from the data, ResponseFamily reports gaussian for a Student-t sampler (so it
// is covered here without a member of its own), and probit and logistic have a
// transform fixed by the link that updateScale does not touch at all, leaving a
// grouped binary sampler's TRUE the documented no-op it already is. Anything
// but FALSE is refused, the sibling guard's condition-keying: the two surfaces
// convert to the engine's bool differently (the R bridge on == TRUE, the flat
// API on != 0), so only the value both read as "pin it" is let through.
// Weights carry no transform and never reach here. External linkage: the flat
// C API reuses this guard on its own setResponse and setOffset entries.
void refuseGroupedScaleUpdate(const bartcore::SamplerBase& sampler,
                              const char* caller, ResponseConduit conduit,
                              int updateScale) {
  if (conduit == ResponseConduit::weights || updateScale == FALSE) return;
  bartcore::SamplerShape shape = sampler.shape();
  if (shape.numGroups == 0) return;
  if (shape.family != bartcore::ResponseFamily::gaussian &&
      shape.family != bartcore::ResponseFamily::aft) return;
  Rf_error("%s: a grouped sampler holds its random intercepts b and their "
           "scale tau against the response transform fixed at creation and "
           "converts neither, so %s swap is supported only with "
           "updateScale = FALSE, which pins it", caller,
           conduit == ResponseConduit::response ? "a response" : "an offset");
}

// The weight policy, stated once for creation and every mutation conduit: a
// probit has no tractable weighted latent-variable form and is refused;
// logistic treats weights as observation counts (its PG(w, psi) latent is the
// sum of w PG(1, psi) draws), so they must be positive integers; gaussian
// takes any finite non-negative weight. The R layer mirrors this, so these
// errors backstop direct-API consumers, and the mutation entries reuse it
// rather than stating a second text.
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
  // a gaussian weight enters the leaf sufficient statistics as a precision, so
  // a negative one subtracts information and NaN/Inf poisons the sum - both
  // fit silently rather than erroring. !(w >= 0.0) catches NaN. The O(n) scan
  // is at a setter, and the flat entrance has no R layer ahead of it to have
  // scanned already
  if (family == bartcore::ResponseFamily::gaussian)
    for (size_t i = 0; i < numObservations; ++i)
      if (!(weights[i] >= 0.0) || !std::isfinite(weights[i]))
        Rf_error("gaussian case weights must be finite and non-negative: a "
                 "weight is a precision multiplier, sd_i = sigma / sqrt(w_i)");
}

// The post-creation half of that policy: gaussian and logistic accept a weight
// change - the logistic one is a model change with a defined meaning, since
// the counts are its Polya-Gamma shape and the engine redraws the latents
// against them - while probit, ordinal, aft and nbinom carry no weights to
// change under any surface. Naming the actual family matters: the one message
// this used to carry told an aft, ordinal or nbinom caller about "a binary
// response" they had not asked for.
// External linkage: the flat C API answers the same condition on
// dbarts_sampler_setWeights, through the predicate below, rather than
// dropping a probit/ordinal/aft/nbinom weight change silently as it once did.
bool familyCarriesNoWeights(const bartcore::SamplerBase& sampler) {
  bartcore::ResponseFamily family = sampler.shape().family;
  return family != bartcore::ResponseFamily::gaussian &&
         family != bartcore::ResponseFamily::logistic;
}

void refuseBinaryWeightChange(const bartcore::SamplerBase& sampler) {
  if (!familyCarriesNoWeights(sampler)) return;
  bartcore::ResponseFamily family = sampler.shape().family;
  const char* name = "probit";
  if (family == bartcore::ResponseFamily::aft) name = "aft (survival)";
  else if (family == bartcore::ResponseFamily::ordinal) name = "ordinal";
  else if (family == bartcore::ResponseFamily::nbinom) name = "nbinom (count)";
  Rf_error("%s models do not support case weights, so none can be set after "
           "creation; fit a gaussian model for a weighted likelihood", name);
}

// The largest count any surface accepts for nbinom. The bound is an ALLOCATION
// bound: NBDispersionPrior::computeKernel sizes its count histogram as
// maxCount + 1 doubles, 8 bytes per unit of the largest count, so y = 1e9 asks
// for 8 GB where no R error can be raised, while this bound pins the request
// at 8 * (1e6 + 1) = 8 MB and the kernel's rebuild at 13e6 multiply-adds -
// trivially safe on every host, and orders of magnitude past any count a
// negative-binomial regression is a sensible model for. The exact Polya-Gamma
// augmentation's O(y + r) draw cost stays a recorded family cost above and
// below the bound; it is not what this refuses.
constexpr double maximumCount = 1.0e6;

// The post-creation half of the response-support policy the R surface applies
// at creation (R/spec.R): mutation must accept exactly what creation does, or a
// swap walks the sampler off its family's support. Two harms, both confirmed:
// probit/ordinal latents drawn against an out-of-support y are silently garbage
// (a non-0/1 y drives probit latents into the hundreds), and an nbinom count
// that is negative underflows NBDispersionPrior::computeKernel's
// static_cast<size_t>(lround(y)) into a ~1.8e19 histogram allocation - an
// uncatchable crash, not an error. gaussian and aft impose nothing (aft's y is
// a log survival time, any real), so they pass through; multinomial counts are
// not reachable by this conduit. Magnitude is bounded for the same allocation
// reason the sign is (see maximumCount), at creation and at every mutation
// alike. A multinomial sampler reports the logistic family, but the
// multi-forest response guard refuses every conduit that reaches this ahead of
// it, so its counts are never read against the binary rule.
// External linkage: the creation prologue and the flat C API both call this, so
// creation and every mutation conduit state one rule.
void validateResponseSupport(bartcore::ResponseFamily family,
                             size_t numCategories, const double* y,
                             size_t numObservations, const char* caller) {
  if (y == NULL) return;
  switch (family) {
  case bartcore::ResponseFamily::probit:
  case bartcore::ResponseFamily::logistic:
    for (size_t i = 0; i < numObservations; ++i)
      if (y[i] != 0.0 && y[i] != 1.0)
        Rf_error("%s: a binary (probit/logistic) response must be coded 0 or "
                 "1; fit family = \"gaussian\" for a continuous response",
                 caller);
    break;
  case bartcore::ResponseFamily::ordinal:
    for (size_t i = 0; i < numObservations; ++i)
      if (!std::isfinite(y[i]) || y[i] != std::floor(y[i]) || y[i] < 1.0 ||
          y[i] > static_cast<double>(numCategories))
        Rf_error("%s: an ordinal response must be an integer category index "
                 "in [1, %lu], the coding fixed when the sampler was created; "
                 "make a new sampler to change the category set", caller,
                 static_cast<unsigned long>(numCategories));
    break;
  case bartcore::ResponseFamily::nbinom:
    for (size_t i = 0; i < numObservations; ++i) {
      if (!std::isfinite(y[i]) || y[i] < 0.0 || y[i] != std::floor(y[i]))
        Rf_error("%s: family \"nbinom\" requires a non-negative integer "
                 "(count) response", caller);
      if (y[i] > maximumCount)
        Rf_error("%s: family \"nbinom\" requires counts no larger than %.0f; "
                 "the dispersion grid's count histogram is sized from the "
                 "largest count, so a larger one allocates without bound",
                 caller, maximumCount);
    }
    break;
  // gaussian and aft constrain nothing (any real y). Every enumerator is
  // listed and there is no default arm: a family added without one would get
  // NO validation, which is exactly the silent-garbage and huge-allocation
  // pair above, so it must fail the build instead.
  case bartcore::ResponseFamily::gaussian:
  case bartcore::ResponseFamily::aft:
    break;
  }
}

// The single-forest test-fit and prediction surface has no meaning on a
// coupling whose amplitudes have no off-sample basis to multiply: the blend
// sum_f dot(a_f, B_f(i,.)) f_f(x_i) is ill-defined off the training rows, so
// the engine would fall back to the bare first forest and silently misreport.
// Reject test data and out-of-sample prediction there; consumers recombine per
// forest via getForestFits + the amplitude glue. Gated on
// testFitsAreDefined rather than the forest count so a multi-forest model
// whose test blend IS defined (multinomial softmax over the K forests'
// totalTestFits) is allowed through; only the couplings that leave the channel
// undefined are refused. External linkage: the flat C API answers the same
// condition on its own predict and test-predictor entries, through the
// predicate below.
bool testFitsAreUndefined(const bartcore::SamplerBase& sampler) {
  bartcore::SamplerShape shape = sampler.shape();
  return shape.numForests >= 2 && !shape.testFitsAreDefined;
}

void refuseUndefinedTestFits(const bartcore::SamplerBase& sampler,
                             const char* caller) {
  if (!testFitsAreUndefined(sampler)) return;
  Rf_error("%s: this sampler's forest amplitudes have no off-sample basis, "
           "so its combined test fits are undefined; replay the forests "
           "separately and recombine them with the amplitude glue instead",
           caller);
}

void refuseEmptyTreeStore(const bartcore::SamplerBase& sampler,
                          const char* caller) {
  bartcore::SamplerShape shape = sampler.shape();
  if (shape.savedTreeCapacity > 0 && shape.numSavedDraws == 0)
    Rf_error("%s: the saved-tree store holds no recorded draws; run the "
             "sampler with keepTrees on before reading from it", caller);
}

// A sampler whose residual sd is not a free parameter has none to set:
// Chain::setSigma installs the value unconditionally and the constant-leaf
// draws divide by sigma_ * sigma_, but sigmaIsFixed_ gates off the redraw that
// would correct it, so the value silently rescales every leaf posterior
// precision for the sampler's life. Two structurally pinned cases: a family
// whose sigmaScale() is 1 by the model definition (probit and logistic - the
// mark multinomial gives itself - and the latent-augmented ordinal/nbinom),
// and a heteroscedastic variance forest, which IS the residual variance
// (buildVarianceForest leaves family_ gaussian, so the family test alone would
// miss it). Keyed on the family, NOT on sigmaIsFixed_: a gaussian sampler with
// resid.prior = fixed() pins sigma too, and driving it per sweep is the
// supported outer-Gibbs conditioning idiom.
// External linkage: the flat C API answers the same condition on
// dbarts_sampler_setSigma, through the predicate below.
bool sigmaIsPinned(const bartcore::SamplerBase& sampler) {
  bartcore::SamplerShape shape = sampler.shape();
  return shape.hasVarianceForest ||
         (shape.family != bartcore::ResponseFamily::gaussian &&
          shape.family != bartcore::ResponseFamily::aft);
}

void refusePinnedSigmaChange(const bartcore::SamplerBase& sampler,
                             const char* caller) {
  if (!sigmaIsPinned(sampler)) return;
  if (sampler.shape().hasVarianceForest)
    Rf_error("%s: this sampler's variance forest owns the residual scale; "
             "there is no single sigma to set", caller);
  // the only other disjunct, so no second test
  Rf_error("%s: this response family fixes the residual standard deviation "
           "by definition; only gaussian and aft samplers carry a sigma to "
           "set", caller);
}

// The value half of the active-row contract, apart from the engine's
// capability bool: a mask element that is neither 0 nor 1 is recoverable - a
// different mask would have worked - so it raises on both surfaces, leaving
// the engine's false to mean "this family implements no mask" alone. The scan
// is the engine's own (bartcore/chain.hpp), restated here so a surface can
// refuse before the install rather than after it.
void refuseNonBinaryMask(const double* active, size_t numObservations) {
  if (active == NULL) return;
  for (size_t i = 0; i < numObservations; ++i)
    if (active[i] != 0.0 && active[i] != 1.0)
      Rf_error("active rows must be exactly 0 or 1: a fractional value is a "
               "weighted likelihood, which the latent families have no "
               "coherent form for");
}

void validateColumnValues(const bartcore::ColumnStore& store, size_t column,
                          const double* values, size_t numValues) {
  if (store.types[column] != bartcore::ColumnType::categorical) return;
  for (size_t i = 0; i < numValues; ++i) {
    if (!store.categoricalValueIsValid(column, values[i]))
      Rf_error("categorical predictor values must be existing category codes");
  }
}

// A sparse column's declared reference level is the code its IMPLICIT rows
// take, which only a categorical column has: whatever the container declares,
// the engine reads an ordinal column's implicit rows as the quantized zero, so
// a reference against one is malformed rather than an alternative implicit
// value - and the container's own densification would read that column
// differently from the engine. \p storeTypes is indexed by SOURCE column, so a
// subset mutation passes the types of the columns it names. External linkage:
// promoted to the bridge block so the flat C API's source-shaped entries state
// this rule once rather than restating it.
void refuseCscReferenceAgainstStore(const bartcore::ColumnType* storeTypes,
                                    const std::int32_t* columnSources,
                                    size_t numColumns, const int* referenceMeta,
                                    size_t numSparseColumns) {
  if (referenceMeta == NULL) return;
  for (size_t j = 0; j < numColumns; ++j) {
    if (columnSources[j] >= 0 ||
        storeTypes[j] == bartcore::ColumnType::categorical)
      continue;
    size_t source = static_cast<size_t>(~columnSources[j]);
    if (source < numSparseColumns && referenceMeta[source] != NA_INTEGER)
      Rf_error("a sparse predictor column may declare a reference level only "
               "for a categorical predictor");
  }
}

// A designated leaf covariate reads contiguous raw values, which CSC storage
// does not serve. The test-store entrances answer this with setTestData's false
// return; a read-only replay builds no store, so it checks the view itself and
// raises the same text. External linkage: promoted to the bridge block so the
// flat C API's source-shaped entries state this rule once rather than
// restating it.
void refuseSparseLeafCovariate(const bartcore::SamplerShape& shape,
                               const bartcore::PredictorSource& source) {
  for (size_t k = 0; k < shape.numLeafCovariates; ++k)
    if (source.sourceOf(shape.leafCovariateColumns[k]) < 0)
      Rf_error("a leaf covariate column cannot be a sparse test column; "
               "supply it as a dense test column");
}

// Bound a parsed test view's categorical codes against the STORE's fixed
// category counts - the training-side bound the view's author cannot see,
// since its own declared K is the caller's, not the sampler's. Covers a
// dense-backed column's slice, a CSC-backed column's stored codes, and the
// reference code its implicit rows read. The container mutation entrances run
// this after parseTestContainer and before any store change: creation-time
// validation is long past there, and setTestData would otherwise quantize an
// unbounded code straight into the test store. Takes the bare view rather
// than the R-side parse it came from, which no header can name. External
// linkage: promoted to the bridge block so the flat C API's source-shaped
// entries state this rule once rather than restating it.
void validateTestContainerAgainstStore(
    const bartcore::ColumnStore& store,
    const bartcore::PredictorSource& view) {
  for (size_t j = 0; j < store.numPredictors; ++j) {
    if (store.types[j] != bartcore::ColumnType::categorical) continue;
    double bound = static_cast<double>(store.numCuts[j]);
    if (view.sourceOf(j) >= 0) {
      refuseInvalidCategoryCodes(rawViewColumn(view, j), view.numRows, bound,
                                 categoricalTestMessage);
      continue;
    }
    ParsedCscCodes stored = parsedCscCodes(view.cscColumnPointers,
                                           view.cscValues, view.sourceOf(j));
    refuseInvalidCategoryCodes(stored.values, stored.numValues, bound,
                               categoricalTestMessage);
    double reference = static_cast<double>(view.referenceCodeOf(j));
    refuseInvalidCategoryCodes(&reference, 1, bound, categoricalTestMessage);
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
                 varianceColumns = std::vector<std::size_t>{},
                 amplitudeSpec = bartcore::AmplitudeSpec{},
                 amplitudeStorage = AmplitudeSpecStorage{},
                 rngs = std::vector<ext_rng*>{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    bartcore::ResponseFamily family = parseSamplerSpecification(
      controlExpr, modelExpr, dataExpr, familyName, control, model, data,
      sigmaIsFixed);
    validateCategoricalPredictors(data);
    requireResolvedSigmaEstimate(family, data.sigmaEstimate);

    bartcore::SamplerOptions options =
      optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);

    // grouped random intercepts (rbart_vi's in-core path) arrive on an
    // internal control attribute; the chains copy the indices at construction
    applyGroupAttribute(controlExpr, data.numObservations, options,
                        groupIndices);
    // grouped ordinal is a recorded but unbuilt door: the threshold block and
    // the group block are not yet shown to interleave, so refuse the
    // composition here, the host backstop the R surface (rbart_vi) mirrors
    if (family == bartcore::ResponseFamily::ordinal && options.numGroups > 0)
      Rf_error("grouped random effects are not supported for ordinal responses");
    // grouped nbinom is a recorded but unbuilt door: the dispersion block
    // and the group block are not yet shown to interleave, so refuse the
    // composition here, the backstop rbart_vi mirrors
    if (family == bartcore::ResponseFamily::nbinom && options.numGroups > 0)
      Rf_error("grouped random effects are not supported for count (nbinom) responses");
    // AFT survival status arrives the same way; the response copies it
    applySurvivalAttribute(controlExpr, data.numObservations, family, options,
                           survivalStatus);
    // the heteroscedastic variance forest arrives on a control attribute; the
    // factory refuses it for non-gaussian or non-constant-leaf models
    applyVarianceAttributes(controlExpr, data.numPredictors, options,
                            varianceColumns);
    // grouped random effects and a variance forest is an unadjudicated
    // composition: the group block draws b at the scalar sigma a variance
    // forest pins at 1, so the effects condition on a residual variance the
    // model does not have. Backstop for the entrances that skip the R
    // surface's own refusal.
    if (options.numGroups > 0 && options.numVarianceTrees > 0)
      Rf_error("grouped random effects are not supported with a "
               "heteroscedastic variance forest");

    // opt-in fp32 residual (storage = "single") is v1-scoped to the gaussian
    // PLAIN constant-leaf path; refuse
    // every other model rather than silently ignore the request. The monotone
    // constrained leaf is a distinct (double-only) instantiation the factory
    // dispatches ahead of the fp32 branch, so refuse it here too (it is not yet
    // R-exposed, but this keeps the flag from being silently dropped).
    if (options.fp32Residual) {
      bool monotoneActive = false;
      for (std::int8_t d : model.monotoneDirections)
        if (d != 0) { monotoneActive = true; break; }
      if (family != bartcore::ResponseFamily::gaussian ||
          options.numLeafCovariates != 0 || options.gpLeaves || monotoneActive)
        Rf_error("%s", storageSingleUnsupportedMessage);
      if (options.numVarianceTrees > 0)
        Rf_error("storage = \"single\" is not supported with a "
                 "heteroscedastic variance forest");
    }

    // The multi-forest model is selected by the two halves of a spec: the
    // per-forest bases on the data object and the forests'
    // configuration on the control attribute. Cross-check them in BOTH
    // directions, so a half stripped in transit (a setControl that drops the
    // attribute, a data object rebuilt without the slot) is a loud refusal
    // naming the missing piece rather than a silent single-forest fit.
    // the spec's forests vector owns its storage, so it rides the closure
    // beside the buffers it borrows: the refusals below jump past every local
    bool carriesAmplitudes =
      applyForestAttributes(controlExpr, model, options.numTrees,
                            data.numPredictors, amplitudeSpec,
                            amplitudeStorage);
    if (carriesAmplitudes && data.bases.empty())
      Rf_error("a basis forest was configured but the data carry no forest "
               "bases");
    if (!carriesAmplitudes && !data.bases.empty())
      Rf_error("the data carry forest bases but no basis forest was "
               "configured; build the sampler through dbarts() or dbartsSpec()");
    if (carriesAmplitudes) {
      refuseUnsupportedAmplitudeComposition(family, model, data, options);
      amplitudeSpec.family = family;
      if (data.bases.size() != amplitudeSpec.forests.size())
        Rf_error("the data carry %lu forest bases but %lu forests were "
                 "configured",
                 static_cast<unsigned long>(data.bases.size()),
                 static_cast<unsigned long>(amplitudeSpec.forests.size()));
      // the chains COPY every basis at construction, so this row-major
      // transposition need only outlive the build
      applyForestBases(data, amplitudeSpec, amplitudeStorage);
    }

    rngs = createChainRngs(control, options.numChains);

    // dispatches on the leaf model: a linear node prior's designated columns
    // select the linear-leaf instantiation, everything else the constant leaf
    std::unique_ptr<bartcore::SamplerBase> sampler =
      carriesAmplitudes
        ? bartcore::createAmplitudeSampler(
            data.predictors.denseValues, data.y, data.numObservations,
            data.numPredictors, data.weights, data.offset, data.sigmaEstimate,
            model.sigmaDf, model.sigmaRawScale, options, amplitudeSpec,
            rngs.data())
        : bartcore::createSampler(
            data.predictors.denseValues, data.y, data.numObservations,
            data.numPredictors, data.weights, data.offset, family,
            data.sigmaEstimate, model.sigmaDf, model.sigmaRawScale, options,
            rngs.data());
    if (sampler == NULL) {
      // the R surface refuses these first, so only an entrance that skips it
      // (the flat C API) reaches this
      for (ext_rng* rng : rngs) if (rng != NULL) ext_rng_destroy(rng);
      Rf_error("%s",
               carriesAmplitudes
                 ? "invalid forest specification"
                 : "invalid sampler specification: either the leaf covariate "
                   "designation is not one the leaf models can take, or a "
                   "variance forest is combined with a non-gaussian family, "
                   "Student-t residuals, leaf covariates, or a monotone "
                   "constraint");
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
        sampler->setTestPredictors(data.testPredictors.denseValues,
                               data.numTestObservations);
      }
      sampler->setTestOffset(data.testOffset);
    }

    if (control.verbose) printInitialSummary(control, model, data, *sampler);

    holder = new BartcoreHolder{std::move(sampler), std::move(rngs),
                                control.keepTrainingFits};
    // one empty per-forest weight slot per forest keeps the engine's
    // pass-through until a caller installs one; nothing else is owned here,
    // since every basis was copied into the chains at construction
    if (carriesAmplitudes)
      holder->ownedForestWeights.resize(holder->sampler->shape().numForests);
    return R_NilValue;
  });
  return holder;
}

/// A K-forest amplitude sampler reached from the internal bcf constructor:
/// bases supplies each forest's amplitude basis and
/// bcfParams the K-length list of length-8 parameter vectors (tree count, base,
/// power; the node-scale factor and divisor; the amplitude prior's variance and
/// half-Cauchy scale; the amplitude update flag), the first forest taking its
/// tree count and structure prior from the model spec instead. vars holds each
/// forest's optional 1-based column restriction, and the two trailing
/// expressions the per-forest interactions() and blocks() lists; each may be
/// null. Family admission is refusedAmplitudeFamilyReason's. The public
/// creation route reaches the same build through createHolder, which reads the
/// same pieces off a control attribute, and is the only route the flat C API
/// has here.
BartcoreHolder* createBCFHolder(SEXP controlExpr, SEXP modelExpr,
                                SEXP dataExpr, SEXP basesExpr,
                                SEXP bcfParamsExpr, SEXP varsExpr,
                                SEXP interactionsExpr, SEXP blocksExpr) {
  BartcoreHolder* holder = nullptr;
  unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                 model = ParsedModel{}, rngs = std::vector<ext_rng*>{},
                 spec = bartcore::AmplitudeSpec{},
                 storage = AmplitudeSpecStorage{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    // the family is read off the model this route was handed rather than
    // threaded beside it, which is what makes the two agree structurally:
    // the R-level route already resolves dbartsModel@family away from
    // "auto" and maps that lone remaining case to the bridge's own ""
    // dispatch, so every
    // existing caller derives exactly what it derived before
    SEXP familyExpr = Rf_getAttrib(modelExpr, Rf_install("family"));
    const char* familyName =
      Rf_isString(familyExpr) && rc_getLength(familyExpr) >= 1 &&
          STRING_ELT(familyExpr, 0) != NA_STRING &&
          std::strcmp(CHAR(STRING_ELT(familyExpr, 0)), "auto") != 0
        ? CHAR(STRING_ELT(familyExpr, 0))
        : "";
    bartcore::ResponseFamily family =
      parseSamplerSpecification(controlExpr, modelExpr, dataExpr, familyName,
                                control, model, data, sigmaIsFixed);
    validateCategoricalPredictors(data);
    if (const char* refused = refusedAmplitudeFamilyReason(family))
      Rf_error("a treatment forest does not support %s", refused);
    requireResolvedSigmaEstimate(family, data.sigmaEstimate);
    if (!data.predictors.isDenseBlock())
      Rf_error("a treatment forest requires dense predictors");

    bartcore::SamplerOptions options =
      optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);
    if (options.fp32Residual)
      Rf_error("%s", storageSingleUnsupportedMessage);

    spec.family = family;
    applyAmplitudeSpec(bcfParamsExpr, varsExpr, interactionsExpr, blocksExpr,
                       model, options.numTrees, data.numPredictors, spec,
                       storage);
    // the bases arrive as an argument here rather than on the data object, so
    // they are read into the same ParsedData channel the public route fills
    // and transposed by the same helper
    readForestBases(basesExpr, data);
    if (data.bases.size() != spec.forests.size())
      Rf_error("%lu forest bases were given but %lu forests were configured",
               static_cast<unsigned long>(data.bases.size()),
               static_cast<unsigned long>(spec.forests.size()));
    applyForestBases(data, spec, storage);

    rngs = createChainRngs(control, options.numChains);

    std::unique_ptr<bartcore::SamplerBase> sampler =
      bartcore::createAmplitudeSampler(
        data.predictors.denseValues, data.y, data.numObservations,
        data.numPredictors, data.weights, data.offset, data.sigmaEstimate,
        model.sigmaDf, model.sigmaRawScale, options, spec, rngs.data());
    // the factory returns null on a composition it cannot build; storing that
    // unchecked would hand back a live external pointer wrapping a null
    // sampler, which every entry dereferences
    if (sampler == NULL) {
      for (ext_rng* rng : rngs) if (rng != NULL) ext_rng_destroy(rng);
      Rf_error("invalid basis forest specification");
    }

    holder = new BartcoreHolder{std::move(sampler), std::move(rngs),
                                control.keepTrainingFits};
    // one empty per-forest weight slot per forest; nothing is installed until
    // a caller asks, so the engine keeps its pass-through
    holder->ownedForestWeights.resize(holder->sampler->shape().numForests);
    return R_NilValue;
  });
  return holder;
}

// The parse and validation both multinomial entries share: parse the sampler
// spec, refuse the response combinations the single-trial softmax cannot carry
// (case weights, an offset, a test offset, a mixed test store), and require
// dense predictors and K >= 2. Runs before any rng is created, so a refusal
// needs no cleanup. The predictors ride the data object; the data's response
// is ignored (the counts are the response, borrowed separately).
static void parseMultinomialData(SEXP controlExpr, SEXP modelExpr,
                                 SEXP dataExpr, ParsedControl& control,
                                 ParsedModel& model, ParsedData& data,
                                 bool& sigmaIsFixed, size_t numCategories) {
  parseSamplerSpecification(controlExpr, modelExpr, dataExpr, "", control, model,
                            data, sigmaIsFixed);
  validateCategoricalPredictors(data);
  if (!data.predictors.isDenseBlock())
    Rf_error("multinomial requires dense predictors");
  // an integer weight is exactly the row-wise count replication the counts
  // response already expresses, and a non-integer one would need a
  // real-shape PG(w_i n_i, .) draw with no exact sampler
  if (data.weights != NULL)
    Rf_error("multinomial (softmax) models do not support case weights: an "
             "integer weight is already expressible as row-wise count "
             "replication in the response, and a non-integer one has no "
             "exact augmentation sampler");
  // the softmax is invariant to a common per-observation shift, so the HOST
  // data object's flat offset points exactly along the null direction and is
  // identically inert; the meaningful one is the n x K category offset the
  // creation entries take separately
  if (data.offset != NULL)
    Rf_error("multinomial (softmax) models do not support a flat offset; the "
             "category offset is an n x K matrix");
  // same reasoning as the train-side flat offset above: the meaningful test
  // offset is the nTest x K category test offset, an internal-channel-only
  // capability (this creation entry's own categoryTestOffsetExpr argument, or
  // bartcore_setCategoryTestOffset after) - never this HOST data object's flat
  // one, which is the same null direction the train side is
  if (data.testOffset != NULL)
    Rf_error("multinomial (softmax) models do not support a flat test "
             "offset; the category test offset is an nTest x K matrix "
             "(bartcore_setCategoryTestOffset)");
  if (data.numTestObservations > 0 && data.testIsMixed)
    Rf_error("multinomial (softmax) models require a dense test matrix");
  if (numCategories < 2)
    Rf_error("multinomial requires at least two categories");
}

// Builds the K-forest multinomial sampler shared by both entries: sizes the
// options, creates the per-chain rngs, builds the sampler from the count-native
// spec, and sets the test predictors and their offset. counts is the borrowed
// category-major n x K matrix and trials the per-observation trial counts; both
// the label entry (one-hot, unit trials) and the count entry route through here
// so their draw streams are the one code path. offset is the borrowed n x K
// category offset in the same layout, or null, and testOffset its nTest x K
// test twin. Leaf-scale and k follow the multinomial calibration, not the host
// node prior (a gaussian default).
static std::unique_ptr<bartcore::SamplerBase> buildMultinomialSampler(
    const ParsedControl& control, const ParsedModel& model,
    const ParsedData& data, SEXP modelExpr, bool sigmaIsFixed,
    size_t numCategories, const int* counts, const int* trials,
    const double* offset, const double* testOffset,
    std::vector<ext_rng*>& rngs) {
  bartcore::SamplerOptions options =
    optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);
  if (options.fp32Residual)
    Rf_error("%s", storageSingleUnsupportedMessage);
  // the K category forests take their leaf scale from the softmax calibration
  // map below, never from the host node prior, so a named calibration has
  // nowhere to land; refuse it rather than drop it. This is the first
  // node-scale-class refusal on this path - the host's own node.scale is
  // deliberately not read, and carries a gaussian default no user chose.
  if (std::isfinite(model.priorScale))
    Rf_error("a multinomial forest does not support a named 'prior.scale'; "
             "its leaf scale comes from the softmax calibration map");
  rngs = createChainRngs(control, options.numChains);

  bartcore::MultinomialSpec spec;
  spec.numCategories = numCategories;
  spec.counts = counts;
  spec.trials = trials;
  // borrowed for the sampler's lifetime, as the counts are; the caller owns the
  // buffer and the combiner reads it every sweep
  spec.offset = offset;
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
    bartcore::createMultinomialSampler(data.predictors.denseValues,
                                       data.numObservations,
                                       data.numPredictors, options, spec,
                                       rngs.data());
  // the factory returns null on a composition it cannot build - a variance
  // forest, whose precision channel the softmax's own augmentation owns.
  // Storing that unchecked would hand back a live external pointer wrapping a
  // null sampler, which every entry dereferences. The R surface refuses it
  // first; this is the backstop for the callers that reach here without one.
  if (sampler == NULL) {
    for (ext_rng* rng : rngs)
      if (rng != NULL) ext_rng_destroy(rng);
    Rf_error("a multinomial (softmax) model does not support a "
             "heteroscedastic variance forest");
  }

  // test-at-creation: the K forests each accumulate their own totalTestFits in
  // the sweep, and storeSample blends them into the K softmax test
  // probabilities (chain testFitsAreDefined() is true). buildTest copies the
  // dense test values, so nothing is pinned; parseMultinomialData ruled out the
  // host object's FLAT test offset (the softmax's null direction) and a mixed
  // test store. The category test offset installs after the rows it describes.
  if (data.numTestObservations > 0) {
    sampler->setTestPredictors(data.testPredictors.denseValues,
                               data.numTestObservations);
    if (testOffset != NULL) sampler->setCategoryTestOffset(testOffset);
  }
  return sampler;
}

// A K-forest multinomial (softmax) sampler over a GROUPED-COUNT response: Y is
// an n x K nonnegative integer matrix, category-major (R column-major = the
// combiner's counts_ layout, so the buffer copies directly), and the trials
// n_i = sum_k Y_ik must be >= 1 (an empty row carries no information; a PG(0, .)
// point mass at 0 would break the working response). Same engine as the label
// entry, count-native. categoryOffsetExpr is the optional n x K category
// offset (null for none) and categoryTestOffsetExpr its
// optional nTest x K test twin, which requires test rows to describe.
//
// numCategories is taken as a value rather than an expression because it is
// the count matrix's own column count, read by the one caller
// (createMultinomialDataHolder) off countsExpr's dim attribute before this
// is reached - so it is already known to match and needs no re-check here.
BartcoreHolder* createMultinomialCountsHolder(SEXP controlExpr, SEXP modelExpr,
                                              SEXP dataExpr, SEXP countsExpr,
                                              size_t numCategories,
                                              SEXP categoryOffsetExpr,
                                              SEXP categoryTestOffsetExpr) {
  BartcoreHolder* holder = nullptr;
  unwindProtect([&, control = ParsedControl{}, data = ParsedData{},
                 model = ParsedModel{}, rngs = std::vector<ext_rng*>{},
                 counts = std::vector<int>{}, trials = std::vector<int>{},
                 offset = std::vector<double>{},
                 testOffset = std::vector<double>{}]() mutable -> SEXP {
    bool sigmaIsFixed;
    parseMultinomialData(controlExpr, modelExpr, dataExpr, control, model, data,
                         sigmaIsFixed, numCategories);

    size_t n = data.numObservations;
    // The DIMENSIONS are checked, not just the length, exactly as
    // bartcore_setCounts checks them: a transposed matrix has exactly n*K
    // entries and would place every count in the wrong cell. Each refusal
    // reports the matrix's ACTUAL shape, which is what tells a caller that
    // theirs is transposed. dims is an attribute of the rooted countsExpr, so
    // the PROTECT is redundant to that rooting and is what the PROTECT-balance
    // analyzer reads; the refusals below longjmp past it
    SEXP dimsExpr = PROTECT(Rf_getAttrib(countsExpr, R_DimSymbol));
    if (!Rf_isInteger(countsExpr) || Rf_xlength(dimsExpr) != 2)
      Rf_error("multinomial counts must be an n x K integer matrix");
    int numRowsGiven = INTEGER(dimsExpr)[0];
    int numColumnsGiven = INTEGER(dimsExpr)[1];
    UNPROTECT(1);
    if (static_cast<size_t>(numRowsGiven) != n)
      Rf_error("multinomial counts must have one row per observation: %lu are "
               "needed and the matrix is %d x %d",
               static_cast<unsigned long>(n), numRowsGiven, numColumnsGiven);
    const int* src = INTEGER(countsExpr);
    counts.assign(src, src + n * numCategories);
    trials.assign(n, 0);
    for (size_t k = 0; k < numCategories; ++k)
      for (size_t i = 0; i < n; ++i) {
        int y = counts[k * n + i];
        if (y < 0) Rf_error("multinomial counts must be non-negative");
        // the trials are int, as the combiner's PG loop counter is; a row sum
        // that overflows would wrap into a negative or absurd draw count
        if (trials[i] > INT_MAX - y)
          Rf_error("multinomial count row sums must fit in an integer");
        trials[i] += y;
      }
    for (size_t i = 0; i < n; ++i)
      if (trials[i] < 1)
        Rf_error("every multinomial count row must have at least one trial");
    parseCategoryOffset(categoryOffsetExpr, n, numCategories, offset,
                        "multinomial category offset");
    parseCategoryTestOffset(categoryTestOffsetExpr, data.numTestObservations,
                            numCategories, testOffset,
                            "multinomial category test offset");

    std::unique_ptr<bartcore::SamplerBase> sampler = buildMultinomialSampler(
      control, model, data, modelExpr, sigmaIsFixed, numCategories,
      counts.data(), trials.data(),
      offset.empty() ? NULL : offset.data(),
      testOffset.empty() ? NULL : testOffset.data(), rngs);

    holder = new BartcoreHolder{std::move(sampler), std::move(rngs),
                                control.keepTrainingFits};
    // moving keeps the buffers, so the combiner's borrowed counts/trials and
    // both category offsets stay valid
    holder->ownedCounts = std::move(counts);
    holder->ownedTrials = std::move(trials);
    holder->ownedCategoryOffset = std::move(offset);
    holder->ownedCategoryTestOffset = std::move(testOffset);
    return R_NilValue;
  });
  return holder;
}

/// Reads a data-object slot whose class union admits NULL, returning
/// R_NilValue for an unset one. An S4 slot holding NULL reads back as the null
/// SYMBOL rather than as R_NilValue, which every "is this absent" test must
/// therefore ask rc_isS4Null about first.
SEXP readNullableSlot(SEXP objectExpr, const char* name) {
  SEXP slotExpr = Rf_getAttrib(objectExpr, Rf_install(name));
  return rc_isS4Null(slotExpr) ? R_NilValue : slotExpr;
}

/// The PUBLIC spec route's multinomial creation: the count response and
/// both category offsets ride the data object's own
/// slots rather than arriving beside it, so a sampler re-created from its own
/// (control, model, data) triple - which is what getPointer does after a save
/// and load, and what setState does on a dead pointer - finds them exactly
/// where creation found them, with no R-side mirroring to keep in step. K is
/// the count matrix's own column count; no control attribute carries it. An
/// object serialized before the slots existed reads them as null, which is
/// what routes it to the ordinary single-forest path instead of here.
BartcoreHolder* createMultinomialDataHolder(SEXP controlExpr, SEXP modelExpr,
                                            SEXP dataExpr) {
  // an unset slot of a NULL class union reads back as the S4 null SYMBOL, not
  // as R_NilValue, so every read here is normalized the way parseData
  // normalizes 'offset' and 'weights'
  SEXP countsExpr = readNullableSlot(dataExpr, "counts");
  if (Rf_isNull(countsExpr))
    Rf_error("family \"multinomial\" requires a data object carrying an "
             "n x K count matrix");
  SEXP dimsExpr = Rf_getAttrib(countsExpr, R_DimSymbol);
  if (!Rf_isInteger(countsExpr) || Rf_xlength(dimsExpr) != 2)
    Rf_error("multinomial counts must be an n x K integer matrix");
  return createMultinomialCountsHolder(
    controlExpr, modelExpr, dataExpr, countsExpr,
    static_cast<size_t>(INTEGER(dimsExpr)[1]),
    readNullableSlot(dataExpr, "offset.category"),
    readNullableSlot(dataExpr, "offset.category.test"));
}

/// Warm start: seed the sampler's live forests from a donor "bartcoreState"
/// over the same predictors. samplesExpr, when non-null, maps each chain to a
/// 1-based donor-sample index; NULL spreads chains across the donor pool.
/// Raises R errors on malformed states or an incompatible donor. Declared here
/// because its one caller precedes the state readers it is defined among.
void installForests(bartcore::SamplerBase& sampler, SEXP donorStateExpr,
                    SEXP samplesExpr);

} // namespace bartcore_bridge

extern "C" {

// The external pointer's protection slot pins everything the sampler
// borrows: the data expression at creation, and any replacement vectors the
// setters install later. family is as createHolder's.
SEXP bartcore_create(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                     SEXP familyExpr) {
  const char* familyName =
    Rf_isNull(familyExpr) ? "" : CHAR(STRING_ELT(familyExpr, 0));
  // The multinomial (softmax) family is the one whose response is not the data
  // object's y - it is the n x K count matrix beside it - so it takes its own
  // factory rather than a ResponseModel resolveFamily could name. Dispatched
  // here, at the single creation entry every R-level route reaches, so a
  // sampler re-created from its own triple after a save and load comes back as
  // the same K-forest engine.
  if (std::strcmp(familyName, "multinomial") == 0)
    return createExternalHolder(dataExpr, [&]() {
      return bartcore_bridge::createMultinomialDataHolder(controlExpr,
                                                          modelExpr, dataExpr);
    });
  return createExternalHolder(dataExpr, [&]() {
    return bartcore_bridge::createHolder(controlExpr, modelExpr, dataExpr,
                                         familyName);
  });
}

// Builds the two-layer store (cuts + codes) once for sharing across
// row-subset samplers: control contributes useQuantiles, data contributes
// x, the column types, and n.cuts. Internal, with no serialization.
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

    bartcore::ColumnStore* handle = new bartcore::ColumnStore;
    // the handle owns raw only for the declared leaf-covariate columns of a
    // DENSE build; a mapped build's dense-backed columns already serve raw
    // from the store's own block, and its CSC-backed ones serve none
    if (data.predictors.isMapped()) gatherColumns.clear();
    handle->build(data.predictors, data.maxNumCuts.data(), 0,
                  control.useQuantiles,
                  gatherColumns.empty() ? NULL : gatherColumns.data(),
                  gatherColumns.size());

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
  const bartcore::ColumnStore& parent = dataHandleFromExpression(handleExpr);
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
    requireResolvedSigmaEstimate(family, data.sigmaEstimate);

    if (data.numObservations != parent.numObservations ||
        data.numPredictors != parent.numPredictors)
      Rf_error("data does not match the shape the handle was built from");

    if (!Rf_isInteger(trainRowsExpr) || Rf_xlength(trainRowsExpr) == 0)
      Rf_error("'trainRows' must be a non-empty integer vector");
    if (!Rf_isNull(testRowsExpr) && !Rf_isInteger(testRowsExpr))
      Rf_error("'testRows' must be an integer vector or NULL");
    size_t numTrainRows = static_cast<size_t>(Rf_xlength(trainRowsExpr));
    size_t numTestRows = Rf_isNull(testRowsExpr)
      ? 0 : static_cast<size_t>(Rf_xlength(testRowsExpr));

    trainRows.resize(numTrainRows);
    testRows.resize(numTestRows);
    for (size_t i = 0; i < numTrainRows; ++i) {
      int row = INTEGER(trainRowsExpr)[i];
      if (row < 1 || static_cast<size_t>(row) > parent.numObservations)
        Rf_error("train row out of range");
      trainRows[i] = static_cast<size_t>(row - 1);
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
    // the store-view factory keeps the fp64 residual (fp32 is minted only by
    // createSampler's gaussian constant-leaf branch); refuse rather than
    // silently ignore a "single" request
    if (options.fp32Residual)
      Rf_error("storage = \"single\" is not supported for this sampler");
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
      std::move(sampler), std::move(rngs), control.keepTrainingFits};
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

// A K-forest combining sampler; internal. The model spec is forest 0 - and
// carries the family, which is read off its own slot rather
// than passed beside it - bases the per-forest amplitude bases (a null entry
// leaving that forest the implicit intercept), and bcfParams the K-length
// per-forest parameter list.
SEXP bartcore_createBCF(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                        SEXP basesExpr, SEXP bcfParamsExpr, SEXP varsExpr,
                        SEXP interactionsExpr, SEXP blocksExpr) {
  return createExternalHolder(dataExpr, [&]() {
    return bartcore_bridge::createBCFHolder(controlExpr, modelExpr, dataExpr,
                                            basesExpr, bcfParamsExpr, varsExpr,
                                            interactionsExpr, blocksExpr);
  });
}

// Replaces a multinomial (softmax) sampler's response: the n x K
// category-major count matrix and the trials n_i = sum_k y_ik it implies. The
// trees carry over, fitted to the previous counts, so the swap is the response
// mutation every other family reaches through setResponse - which cannot serve
// here, since the combiner names the chain's y out and reads these counts
// directly.
//
// n and K are fixed: every combiner buffer and every forest allocation is
// sized by them, and K is the forest count, which no live sampler can change.
// A length mismatch is refused naming both.
//
// COST: the sweep draws n_i Polya-Gamma variates per observation per category
// (combiner.hpp drawForestGlue sums n_i PG(1, psi)), so a swap to large row
// sums multiplies sweep cost by mean(n_i) - a data property, not a defect, but
// one a caller replacing single-trial labels with grouped counts should know.
//
// Validation is total before anything is installed: the combiner BORROWS
// ownedCounts.data(), so an in-place write would BE the mutation and an error
// midway would leave half-new counts. The scratch is built and checked whole,
// then swapped in, then the engine is pointed at the fresh buffers. The
// scratch rides the unwindProtect closure, so an error frees it.
SEXP bartcore_setCounts(SEXP ptrExpr, SEXP countsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  // The capability probe comes FIRST, and it is not a forest count: BCF also
  // carries several forests and owns no counts, while a future coupling that
  // did would have to opt in here.
  if (!shape.supportsCountsMutation)
    Rf_error("bartcore_setCounts: requires a multinomial (softmax) sampler");
  size_t n = shape.numObservations;
  // K for a counts-owning coupling; the reported-location count IS the category
  // count, so there is no second field to keep in step with it
  size_t K = shape.numReportedLocations;

  return unwindProtect([&, counts = std::vector<int>{},
                        trials = std::vector<int>{}]() mutable -> SEXP {
    // the strict TYPE check creation applies, kept rather than a value test:
    // a double matrix of whole numbers is still the wrong buffer to borrow.
    // NA_INTEGER is INT_MIN, so the nonnegativity test below catches NA too.
    // The DIMENSIONS are checked, not just the length: a transposed matrix has
    // exactly n*K entries and would install every count in the wrong cell
    // dims is an attribute of the rooted countsExpr, so the PROTECT is
    // redundant to that rooting and is what the PROTECT-balance analyzer reads;
    // the refusal below longjmps, which unwinds the stack past it
    SEXP dimsExpr = PROTECT(Rf_getAttrib(countsExpr, R_DimSymbol));
    if (!Rf_isInteger(countsExpr) || Rf_xlength(dimsExpr) != 2 ||
        static_cast<size_t>(INTEGER(dimsExpr)[0]) != n ||
        static_cast<size_t>(INTEGER(dimsExpr)[1]) != K)
      Rf_error("bartcore_setCounts: requires an integer matrix of %lu "
               "observations x %lu categories",
               static_cast<unsigned long>(n), static_cast<unsigned long>(K));
    UNPROTECT(1);
    const int* src = INTEGER(countsExpr);
    counts.assign(src, src + n * K);
    trials.assign(n, 0);
    for (size_t k = 0; k < K; ++k)
      for (size_t i = 0; i < n; ++i) {
        int y = counts[k * n + i];
        if (y < 0) Rf_error("multinomial counts must be non-negative");
        // the trials are int, as the combiner's PG loop counter is; a row sum
        // that overflows would wrap into a negative or absurd draw count
        if (trials[i] > INT_MAX - y)
          Rf_error("multinomial count row sums must fit in an integer");
        trials[i] += y;
      }
    for (size_t i = 0; i < n; ++i)
      if (trials[i] < 1)
        Rf_error("every multinomial count row must have at least one trial");

    holder.ownedCounts.swap(counts);
    holder.ownedTrials.swap(trials);
    holder.sampler->setCounts(holder.ownedCounts.data(),
                              holder.ownedTrials.data());
    return R_NilValue;
  });
}

// Installs (or clears, at NULL) a multinomial sampler's n x K category offset:
// the latent becomes f_ik + o_ik, so the offset enters the margins, the working
// responses and the reported softmax, and never a leaf value. It is the
// response-side counterpart of bartcore_setCounts, and NOT bartcore_setOffset,
// whose vector the response model adds to every reported channel after the
// forests are blended - past the softmax, where a shift means nothing a
// modeller wants (a common per-observation shift is the softmax's own null
// direction, which is why the flat entry stays refused here).
//
// n and K are fixed, as they are for the counts; a length that happens to match
// from a transposed matrix is refused on dimensions. Validation is total before
// anything is installed, for the reason the counts entrance states: the
// combiner BORROWS ownedCategoryOffset.data(), so an in-place write would BE
// the mutation. Clearing returns the sampler to the exact offset-free path,
// pointer for pointer.
//
// This one shifts the TRAIN latent only. The test rows are other rows, so they
// carry their own offset (bartcore_setCategoryTestOffset), and predict takes
// its own per call; neither is derived from this one.
SEXP bartcore_setCategoryOffset(SEXP ptrExpr, SEXP offsetExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  // the capability probe comes FIRST and is not a forest count, exactly as the
  // counts entrance's is: the coupling that owns an n x K count response is the
  // one with an n x K linear predictor to shift
  if (!shape.supportsCountsMutation)
    Rf_error("bartcore_setCategoryOffset: requires a multinomial (softmax) "
             "sampler");
  size_t n = shape.numObservations;
  size_t K = shape.numReportedLocations;

  return unwindProtect([&, offset = std::vector<double>{}]() mutable -> SEXP {
    parseCategoryOffset(offsetExpr, n, K, offset, "bartcore_setCategoryOffset");
    holder.ownedCategoryOffset.swap(offset);
    const double* installed = holder.ownedCategoryOffset.empty()
      ? NULL : holder.ownedCategoryOffset.data();
    holder.sampler->setCategoryOffset(installed);
    return R_NilValue;
  });
}

// Installs (or clears, at NULL) a multinomial sampler's nTest x K category test
// offset: the reported test channel becomes softmax(f_test + o_test), formed
// where the train blend forms softmax(f + o). Nothing else moves - the test
// fits enter no likelihood - so a run under this offset alone is bitwise a run
// without one on every train channel.
//
// The shape is pinned to the CURRENT test rows, which is also why installing
// test predictors is refused while one is here: the offset describes the rows
// being replaced, and no row-count coincidence makes it describe the new ones.
// Clearing first is the explicit way through. Validation is total before
// anything is installed, as at every borrowed-buffer entrance.
SEXP bartcore_setCategoryTestOffset(SEXP ptrExpr, SEXP offsetExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  // the capability probe comes FIRST and is not a forest count: what makes a
  // per-category test shift meaningful is the softmax test blend, which arrived
  // with the count response
  if (!shape.supportsCountsMutation)
    Rf_error("bartcore_setCategoryTestOffset: requires a multinomial (softmax) "
             "sampler");
  size_t nTest = shape.numTestObservations;
  size_t K = shape.numReportedLocations;

  return unwindProtect([&, offset = std::vector<double>{}]() mutable -> SEXP {
    parseCategoryTestOffset(offsetExpr, nTest, K, offset,
                            "bartcore_setCategoryTestOffset");
    holder.ownedCategoryTestOffset.swap(offset);
    const double* installed = holder.ownedCategoryTestOffset.empty()
      ? NULL : holder.ownedCategoryTestOffset.data();
    holder.sampler->setCategoryTestOffset(installed);
    return R_NilValue;
  });
}

// Forest addressing on the entry points that name one: 0-based and
// unconverted, so the bound is the sampler's own forest count.
static size_t forestIndexFrom(SEXP forestExpr,
                              const bartcore::SamplerShape& shape) {
  size_t forestIndex = static_cast<size_t>(Rf_asInteger(forestExpr));
  if (forestIndex >= shape.numForests) Rf_error("forest index out of range");
  return forestIndex;
}

SEXP bartcore_setForestBasis(SEXP ptrExpr, SEXP forestExpr, SEXP basisExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t n = shape.numObservations;
  // A capability probe, not a forest count: a basis is defined only as what the
  // amplitudes multiply, and a K-forest multinomial (K >= 2) defeats a
  // numForests test. totalAmplitudes is 0 off a combiner and 0 for the
  // multinomial combiner, so the single-forest case this already covered keeps
  // its refusal.
  if (holder.sampler->totalAmplitudes() == 0)
    Rf_error("a forest basis requires a sampler whose forests carry "
             "amplitudes");
  size_t forestIndex = forestIndexFrom(forestExpr, shape);
  if (!Rf_isReal(basisExpr) || Rf_xlength(basisExpr) == 0 || n == 0 ||
      Rf_xlength(basisExpr) % static_cast<R_xlen_t>(n) != 0)
    Rf_error("basis length must be a multiple of the number of observations");
  size_t numColumns = static_cast<size_t>(Rf_xlength(basisExpr)) / n;
  return unwindProtect([&, rowMajor = std::vector<double>{}]() mutable -> SEXP {
    // R holds a matrix COLUMN-major and the engine reads a basis ROW-major (the
    // contraction with the forest's amplitude vector is its only read), so the
    // transposition happens here, once, on the way in
    rowMajor.resize(n * numColumns);
    const double* values = REAL(basisExpr);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < numColumns; ++j) {
        double value = values[j * n + i];
        if (!R_FINITE(value)) Rf_error("a forest basis value is not finite");
        rowMajor[i * numColumns + j] = value;
      }
    holder.sampler->setForestBasis(forestIndex, rowMajor.data(), numColumns);
    return R_NilValue;
  });
}

// A caller-supplied per-forest, per-observation weight: a multiplicative
// precision factor on that forest's own leaf conditionals, composing with the
// case weights (bartcore/chain.hpp states the semantics and its two edges).
SEXP bartcore_setForestWeights(SEXP ptrExpr, SEXP forestExpr,
                               SEXP weightsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  // The capability probe comes FIRST, and it is not a forest count: a K-forest
  // multinomial carries several forests and admits no such weight.
  if (!shape.supportsForestWeights) {
    // This is the only per-forest, per-observation channel a caller can reach,
    // so it is where an attempt to restrict a softmax sampler's rows to one
    // category arrives - and that is refused on model grounds rather than for
    // want of an implementation, permanently. A multinomial mask is GLOBAL
    // ($setActiveRows); the coupling's own setActiveRows carries the argument.
    if (shape.supportsCountsMutation)
      Rf_error("a multinomial mask applies to every category: the softmax "
               "margin reads all K forests, so a row cannot leave one "
               "category's likelihood alone");
    Rf_error("bartcore_setForestWeights: requires a sampler that carries "
             "forest amplitudes");
  }
  size_t forestIndex = forestIndexFrom(forestExpr, shape);
  size_t n = shape.numObservations;
  if (!Rf_isReal(weightsExpr) ||
      static_cast<size_t>(Rf_xlength(weightsExpr)) != n)
    Rf_error("forest weight length must match the number of observations");
  // defense in depth: dbartsSampler$setForestWeights already
  // enforces this; load-bearing only for a direct .Call
  // bypassing it
  const double* weights = REAL(weightsExpr);
  for (size_t i = 0; i < n; ++i)
    if (!R_FINITE(weights[i]) || weights[i] < 0.0)
      Rf_error("forest weights must be finite and non-negative");
  // Copy rather than retain: PROT_COUNT is a fixed enum, so a per-forest
  // multiplicity cannot take a protection slot. Growing the outer vector never
  // relocates an inner vector's storage, so an installed pointer survives.
  holder.ownedForestWeights.resize(shape.numForests);
  std::vector<double>& owned = holder.ownedForestWeights[forestIndex];
  owned.assign(weights, weights + n);
  holder.sampler->setForestWeights(forestIndex, owned.data());
  return R_NilValue;
}

// A per-observation 0/1 mask saying which rows are in the data set this sweep:
// an inactive row leaves every sufficient statistic, every family-level
// parameter update and its own latent draw, while keeping its leaf occupancy
// and its fitted value (bartcore/chain.hpp states the semantics). The
// capability probe comes FIRST and never switches on the family - Student-t
// reports as gaussian - then the length; the exact-{0,1} scan and the all-ones
// normalization are the engine's, so every surface inherits them. NULL clears.
// The values are copied into the sampler, so nothing is retained here.
SEXP bartcore_setActiveRows(SEXP ptrExpr, SEXP activeExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  // Every shipped family implements the channel, multinomial included (there
  // the mask is GLOBAL: it lands on the softmax coupling's K interleaved
  // precisions, since the response holds none of its own). The probe stays as
  // the contract's first step - a family that does not override the base
  // refusal is refused here rather than silently ignored, which is the
  // r-c-division defect this channel's shape exists to avoid.
  if (!shape.supportsActiveRows)
    Rf_error("active-row masking is not implemented for this response family");
  if (activeExpr == R_NilValue) {
    holder.sampler->setActiveRows(NULL);
    return R_NilValue;
  }
  size_t n = shape.numObservations;
  if (!Rf_isReal(activeExpr) ||
      static_cast<size_t>(Rf_xlength(activeExpr)) != n)
    Rf_error("active row length must match the number of observations");
  // the value refusal, stated once for both surfaces; with the capability
  // probed above, no other reason to refuse remains
  refuseNonBinaryMask(REAL(activeExpr), n);
  holder.sampler->setActiveRows(REAL(activeExpr));
  return R_NilValue;
}

// The sampler's forest count, the R twin of dbarts_sampler_numForests. The R5
// readers stack their per-forest reads at forest = NULL and need the bound; a
// capability probe cannot supply it, since a plain single-forest sampler
// answers no to every one of them.
SEXP bartcore_numForests(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  return Rf_ScalarInteger(
    static_cast<int>(holder.sampler->shape().numForests));
}

// One forest's amplitudes, its own q_f x numChains matrix, or - at a NULL
// forest - the whole vector stacked forest-major, sum_f q_f x numChains, which
// is the shape the run's own glue channel carries. The vector is RAGGED, forest
// f being as wide as its basis, so a caller wanting one block names it.
SEXP bartcore_getForestAmplitudes(SEXP ptrExpr, SEXP forestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t total = holder.sampler->totalAmplitudes();
  if (total == 0)
    Rf_error("forest amplitudes require a sampler whose forests carry them");
  bool wholeVector = Rf_isNull(forestExpr);
  size_t forestIndex = wholeVector ? 0 : forestIndexFrom(forestExpr, shape);
  size_t offset = 0;
  if (!wholeVector)
    for (size_t f = 0; f < forestIndex; ++f)
      offset += holder.sampler->numForestAmplitudes(f);
  size_t numAmplitudes =
    wholeVector ? total : holder.sampler->numForestAmplitudes(forestIndex);
  size_t numChains = shape.numChains;
  SEXP result = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numAmplitudes),
                                       static_cast<int>(numChains)));
  std::vector<double> amplitudes(total);
  for (size_t c = 0; c < numChains; ++c) {
    holder.sampler->amplitudes(c, amplitudes.data());
    for (size_t j = 0; j < numAmplitudes; ++j)
      REAL(result)[c * numAmplitudes + j] = amplitudes[offset + j];
  }
  UNPROTECT(1);
  return result;
}

// One forest's internal-scale function values, numObservations x numChains
// (forest 0 prognostic, 1 treatment).
SEXP bartcore_getForestFits(SEXP ptrExpr, SEXP forestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t forestIndex = forestIndexFrom(forestExpr, shape);
  size_t n = shape.numObservations;
  size_t numChains = shape.numChains;
  SEXP result = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n),
                                       static_cast<int>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    holder.sampler->forestTotalFits(c, forestIndex, REAL(result) + c * n);
  UNPROTECT(1);
  return result;
}

// The combined per-observation location on the RESPONSE scale and WITHOUT the
// offset, numObservations x numChains - what the recorded training channel
// carries with the offset folded in, which is the read a host driving an outer
// block one sweep at a time needs offset-free. The shape refusal is the
// ENGINE'S, mapped here to a named error: a second test on this side would
// make Chain::fitsWithoutOffset's own branch unreachable from R and from
// tests/cpp alike, which is untested dead code.
SEXP bartcore_getFitsWithoutOffset(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t n = shape.numObservations;
  size_t numChains = shape.numChains;
  SEXP result = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n),
                                       static_cast<int>(numChains)));
  for (size_t c = 0; c < numChains; ++c) {
    if (!holder.sampler->fitsWithoutOffset(c, REAL(result) + c * n)) {
      UNPROTECT(1);
      Rf_error("getFitsWithoutOffset: a multi-location sampler (multinomial) "
               "reports one softmax probability channel per category and not a "
               "single additive location, so it has no offset-free fit; use "
               "predict(x) for the current state's probabilities, which "
               "reports the saved samples instead when keepTrees is TRUE");
    }
  }
  UNPROTECT(1);
  return result;
}

static const char* leafModelName(bartcore::LeafModelKind kind) {
  switch (kind) {
  case bartcore::LeafModelKind::monotone: return "monotone";
  case bartcore::LeafModelKind::linear: return "linear";
  case bartcore::LeafModelKind::gp: return "gp";
  case bartcore::LeafModelKind::constant: break;
  }
  return "constant";
}

// One forest's leaf-prior calibration in RESPONSE units, one ROW per chain -
// the chains carry their own transforms and their own drawn k, so a flattened
// summary would hide a divergence this surface exists to show. The leaf-model
// tag rides as an attribute because it is a property of the sampler, not of a
// chain, and it qualifies what the reported prior sd means: an equality for
// the constant leaf, a stated bound for the other three.
SEXP bartcore_getCalibration(SEXP ptrExpr, SEXP forestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t forestIndex = forestIndexFrom(forestExpr, shape);
  // the seven reported since the reader shipped, then the five calibration-map
  // quantities: NaN on a forest with no map entry, and the two amplitude
  // columns exclusive per forest on one that has
  const char* const columnNames[] = {
    "prior.scale", "prior.sd", "prior.mean", "k",
    "k.has.hyperprior", "response.scale", "response.shift",
    "amplitude.prior.variance", "amplitude.prior.scale",
    "node.scale.factor", "node.scale.divisor", "basis.row.norm"
  };
  size_t numColumns = sizeof columnNames / sizeof columnNames[0];
  size_t numChains = shape.numChains;
  SEXP resultExpr = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numChains),
                                           static_cast<int>(numColumns)));
  double* result = REAL(resultExpr);
  for (size_t c = 0; c < numChains; ++c) {
    bartcore::ForestCalibration calibration =
      holder.sampler->forestCalibration(c, forestIndex);
    result[c] = calibration.priorScale;
    result[c + numChains] = calibration.priorSd;
    result[c + 2 * numChains] = calibration.priorMean;
    result[c + 3 * numChains] = calibration.k;
    result[c + 4 * numChains] = calibration.kHasHyperprior ? 1.0 : 0.0;
    result[c + 5 * numChains] = calibration.responseScale;
    result[c + 6 * numChains] = calibration.responseShift;
    result[c + 7 * numChains] = calibration.amplitudePriorVariance;
    result[c + 8 * numChains] = calibration.amplitudePriorScale;
    result[c + 9 * numChains] = calibration.nodeScaleFactor;
    result[c + 10 * numChains] = calibration.nodeScaleDivisor;
    result[c + 11 * numChains] = calibration.basisRowNorm;
  }
  SEXP dimNamesExpr = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimNamesExpr, 0, R_NilValue);
  SEXP columnNamesExpr =
    PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(numColumns)));
  for (size_t j = 0; j < numColumns; ++j)
    SET_STRING_ELT(columnNamesExpr, static_cast<R_xlen_t>(j),
                   Rf_mkChar(columnNames[j]));
  SET_VECTOR_ELT(dimNamesExpr, 1, columnNamesExpr);
  Rf_setAttrib(resultExpr, R_DimNamesSymbol, dimNamesExpr);
  setAttribByName(resultExpr, "leaf.model",
                  Rf_mkString(leafModelName(shape.leafModel)));
  UNPROTECT(3);
  return resultExpr;
}

// Restates one forest's leaf prior on every chain, response units. Two error
// channels, as the flat entry has: a CAPABILITY answer (no such forest, or a
// combiner owns the calibration) and a MALFORMED VALUE. Neither the response
// transform, k, sigma nor the tree prior moves, and a write reproducing what
// is in force is skipped bitwise inside the engine, so a read-then-write
// cannot perturb a draw.
SEXP bartcore_setCalibration(SEXP ptrExpr, SEXP forestExpr,
                             SEXP priorScaleExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t forestIndex = forestIndexFrom(forestExpr, shape);
  double priorScale = Rf_asReal(priorScaleExpr);
  if (!std::isfinite(priorScale) || priorScale <= 0.0)
    Rf_error("'prior.scale' must be a positive finite number");
  // The map is named so the refusal points at what fixed the scale rather than
  // at the refusal itself. Only the K = 2 coupling is a two-forest map; above
  // it the same map spans K forests, and naming it two-forest would be false.
  if (!holder.sampler->setForestPriorScale(forestIndex, priorScale))
    Rf_error("this forest's leaf scale comes from the %s, which owns both "
             "halves of its calibration; make a new sampler instead",
             shape.supportsCountsMutation ? "softmax calibration map"
             : shape.numForests == 2      ? "two-forest calibration map"
                                          : "multi-forest calibration map");
  return R_NilValue;
}

SEXP bartcore_getForestVariableCounts(SEXP ptrExpr, SEXP forestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t forestIndex = forestIndexFrom(forestExpr, shape);
  size_t numPredictors = shape.numPredictors;
  size_t numChains = shape.numChains;
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

  bartcore::SamplerShape shape = sampler.shape();
  size_t numObservations = shape.numObservations;
  size_t numPredictors = shape.numPredictors;
  size_t numTestObservations = shape.numTestObservations;
  size_t numChains = shape.numChains;
  int numSamplesInt = static_cast<int>(numSamples);
  int numChainsInt = static_cast<int>(numChains);
  // One recorded channel: its leading dimensions, then the samples, then a
  // trailing chain dimension when there is more than one chain. A per-draw
  // scalar on a single chain lands as a bare vector - one dimension carries no
  // dim attribute - which is the shape and byte layout every consumer relies
  // on.
  auto allocChannel = [&](SEXPTYPE type,
                          std::initializer_list<size_t> leading) -> SEXP {
    int dims[4];
    int numDims = 0;
    for (size_t extent : leading) dims[numDims++] = static_cast<int>(extent);
    dims[numDims++] = numSamplesInt;
    if (numChains > 1) dims[numDims++] = numChainsInt;
    R_xlen_t total = 1;
    for (int d = 0; d < numDims; ++d) total *= dims[d];
    SEXP arr = PROTECT(Rf_allocVector(type, total));
    if (numDims > 1) {
      SEXP dimExpr = Rf_allocVector(INTSXP, numDims);
      for (int d = 0; d < numDims; ++d) INTEGER(dimExpr)[d] = dims[d];
      Rf_setAttrib(arr, R_DimSymbol, dimExpr);
    }
    UNPROTECT(1);
    return arr;
  };
  // one scalar per draw
  auto allocScalarChannel = [&]() { return allocChannel(REALSXP, {}); };
  // A channel that widens on an inner per-sample axis carries that axis
  // between its leading dimension and the samples. An inner extent of 1 drops
  // the axis, keeping the exact leadingDim x numSamples (x numChains) shape
  // every unwidened model relies on.
  auto allocWideChannel = [&](SEXPTYPE type, size_t leadingDim,
                              size_t innerDim) -> SEXP {
    return innerDim > 1 ? allocChannel(type, {leadingDim, innerDim})
                        : allocChannel(type, {leadingDim});
  };
  // a multi-location combiner (multinomial: K softmax channels) inserts a
  // location dimension between the observations and the samples
  size_t numLocations = shape.numReportedLocations;
  // the varcount channel widens on its own forest axis (K for multinomial's
  // category forests and for a multi-forest amplitude model's, prognostic
  // first), inserting a forest dimension between the predictors and the
  // samples exactly as the fits seam inserts locations
  size_t numVCForests = shape.numVariableCountForests;

  if (numSamples == 0) {
    bartcore::Results empty;
    GetRNGstate();
    bool cancelled = sampler.run(numBurnIn, 0, empty, bartcore_userInterrupted);
    PutRNGstate();
    if (cancelled) Rf_error("sampler run interrupted");
    return R_NilValue;
  }

  // ordinal reports its K-1 thresholds in an extra channel appended after ranef;
  // every other family carries none, so the list keeps its 8 slots and byte-for-
  // byte layout. numOrdinalThresholds == 0 off ordinal.
  size_t numOrdinalThresholds = shape.numOrdinalThresholds;
  bool hasOrdinalThresholds = numOrdinalThresholds > 0;
  // heteroscedastic samplers append s.train (+ s.test when test rows exist) as
  // a separately-typed variance channel; gaussian-only, so mutually exclusive
  // with the ordinal threshold slot
  bool hasVariance = shape.hasVarianceForest;
  // a coupling that composes its forests through scalar glue (BCF) appends the
  // per-forest fits and that glue, so one run reports both surfaces and every
  // draw's recombination; every other model appends neither slot and the engine
  // computes nothing for them
  bool hasForestReporting = shape.forestReportingIsDefined;
  size_t numForests = shape.numForests;
  // an nbinom sampler appends its per-draw dispersion r right after the ordinal
  // threshold slot, so every later conditional slot shifts by it; no family
  // carries both, but the arithmetic composes regardless of that
  bool hasDispersion = shape.carriesDispersion;
  // a Student-t error law appends its per-draw df nu next, on the same
  // arithmetic; no response carries both, but the count composes regardless
  bool hasResidualDf = shape.carriesResidualDf;
  int numResultSlots = 8 + (hasOrdinalThresholds ? 1 : 0) +
                       (hasDispersion ? 1 : 0) + (hasResidualDf ? 1 : 0) +
                       (hasVariance ? 2 : 0) + (hasForestReporting ? 2 : 0);

  // several chains add a trailing chain dimension. Every column roots in the
  // protected container the moment it is allocated (installChannel), so there
  // is no hand-counted PROTECT stack to keep in sync with the slot list.
  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, numResultSlots));
  // The names vector roots through the container's attribute before the mkChar
  // allocations fill it. Slot order IS install order: a channel claims the next
  // slot and names it in the same step, so the conditional channels need no
  // slot arithmetic and the names cannot drift from the values.
  SEXP namesExpr = Rf_allocVector(STRSXP, numResultSlots);
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);
  int slot = 0;
  auto installChannel = [&](const char* name, SEXP value) -> SEXP {
    // the value roots first: it is unprotected until it is in the container,
    // and mkChar allocates. That ordering is what keeps it alive; the PROTECT
    // is redundant to it and is what the PROTECT-balance analyzer reads, which
    // does not model rooting by container assignment.
    PROTECT(value);
    SET_VECTOR_ELT(resultExpr, slot, value);
    SET_STRING_ELT(namesExpr, slot++, Rf_mkChar(name));
    UNPROTECT(1);
    return value;
  };

  SEXP sigmaExpr = installChannel("sigma", allocScalarChannel());
  SEXP trainExpr = installChannel(
    "train", !holder.keepTrainingFits
               ? R_NilValue
               : allocWideChannel(REALSXP, numObservations, numLocations));
  SEXP testExpr = installChannel(
    "test", numTestObservations == 0
              ? R_NilValue
              : allocWideChannel(REALSXP, numTestObservations, numLocations));
  SEXP varcountExpr = installChannel(
    "varcount", allocWideChannel(INTSXP, numPredictors, numVCForests));
  SEXP kExpr =
    installChannel("k", !shape.kIsSampled ? R_NilValue : allocScalarChannel());
  SEXP varprobsExpr = installChannel(
    "varprobs",
    !shape.usesDart ? R_NilValue : allocChannel(REALSXP, {numPredictors}));
  size_t numGroups = shape.numGroups;
  SEXP tauExpr =
    installChannel("tau", numGroups == 0 ? R_NilValue : allocScalarChannel());
  SEXP ranefExpr = installChannel(
    "ranef", numGroups == 0 ? R_NilValue : allocChannel(REALSXP, {numGroups}));
  SEXP ordinalThresholdsExpr = !hasOrdinalThresholds
    ? R_NilValue
    : installChannel("thresholds",
                     allocChannel(REALSXP, {numOrdinalThresholds}));
  // one scalar per draw, so the dispersion channel takes sigma's own shape
  SEXP dispersionExpr = !hasDispersion
    ? R_NilValue
    : installChannel("dispersion", allocScalarChannel());
  // likewise one scalar per draw, so the df channel takes sigma's shape too
  SEXP residualDfExpr = !hasResidualDf
    ? R_NilValue
    : installChannel("resid.df", allocScalarChannel());
  // the test half keeps its slot and its name with no test rows to fill it
  SEXP varianceTrainExpr = R_NilValue;
  SEXP varianceTestExpr = R_NilValue;
  if (hasVariance) {
    varianceTrainExpr =
      installChannel("variance", allocChannel(REALSXP, {numObservations}));
    varianceTestExpr = installChannel(
      "varianceTest", numTestObservations == 0
                        ? R_NilValue
                        : allocChannel(REALSXP, {numTestObservations}));
  }
  SEXP forestFitsExpr = R_NilValue;
  SEXP glueExpr = R_NilValue;
  if (hasForestReporting) {
    // the forest axis is always present here (the channel exists only for a
    // multi-forest coupling), so this is n x numForests x numSamples
    // (x numChains)
    forestFitsExpr = installChannel(
      "forestFits", allocChannel(REALSXP, {numObservations, numForests}));
    // the glue axis is the RAGGED amplitude vector, sum_f q_f long,
    // forest-major within a draw; bcf's q = (1, 2) makes it the 3 rows
    // (a, b0, b1) it shipped with
    glueExpr =
      installChannel("glue", allocChannel(REALSXP, {shape.numAmplitudes}));
  }

  std::vector<std::uint32_t> variableCounts(numPredictors * numVCForests *
                                            numSamples * numChains);

  bartcore::Results results;
  results.sigma = REAL(sigmaExpr);
  results.trainingFits = holder.keepTrainingFits ? REAL(trainExpr) : NULL;
  results.testFits = numTestObservations > 0 ? REAL(testExpr) : NULL;
  results.variableCounts = variableCounts.data();
  results.k = shape.kIsSampled ? REAL(kExpr) : NULL;
  results.splitProbabilities = shape.usesDart ? REAL(varprobsExpr) : NULL;
  results.tau = numGroups > 0 ? REAL(tauExpr) : NULL;
  results.groupEffects = numGroups > 0 ? REAL(ranefExpr) : NULL;
  // one per-observation fits channel per reported location; 1 for every
  // additive model, K for multinomial. The location stride drives the
  // chain-major slabbing (multiple chains).
  results.numReportedLocations = numLocations;
  // one varcount slab per reported forest; 1 for a single-forest model, K for
  // multinomial's per-category counts and for a multi-forest amplitude model's
  // per-forest counts. Drives the chain-major varcount stride. The R run bridge
  // is the surface that OPTS IN to the widened channel; a caller leaving the
  // field at its default 1 (the flat C API) keeps the single prognostic slab.
  results.numVariableCountForests = numVCForests;
  // K-1 thresholds per sample for ordinal, none otherwise; drives the
  // chain-major threshold stride
  results.ordinalThresholds =
    hasOrdinalThresholds ? REAL(ordinalThresholdsExpr) : NULL;
  results.numOrdinalThresholds = numOrdinalThresholds;
  // the dispersion r each draw is conditioned on; null off nbinom, which is the
  // guard storeSample's write shares
  results.dispersion = hasDispersion ? REAL(dispersionExpr) : NULL;
  // the residual df nu each draw is conditioned on; null off a Student-t error
  // law, the guard storeSample's write shares
  results.residualDf = hasResidualDf ? REAL(residualDfExpr) : NULL;
  results.varianceFits = hasVariance ? REAL(varianceTrainExpr) : NULL;
  results.varianceTestFits =
    (hasVariance && numTestObservations > 0) ? REAL(varianceTestExpr) : NULL;
  // each forest's own internal-scale fits and the (a, b0, b1) that recombines
  // them, per draw; both null unless the coupling defines the channels
  results.forestFits = hasForestReporting ? REAL(forestFitsExpr) : NULL;
  results.glue = hasForestReporting ? REAL(glueExpr) : NULL;

  GetRNGstate();
  bool cancelled = sampler.run(numBurnIn, numSamples, results,
                               bartcore_userInterrupted);
  PutRNGstate();
  if (cancelled) {
    std::vector<std::uint32_t>().swap(variableCounts);  // free before longjmp
    Rf_error("sampler run interrupted");
  }

  // nothing past the copy-out allocates, so the counts need no early free
  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < numPredictors * numVCForests * numSamples * numChains;
       ++i)
    varcountOut[i] = static_cast<int>(variableCounts[i]);

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
  bartcore::SamplerShape shape = sampler.shape();
  if (shape.numChains != 1)
    Rf_error("bartcore_runWithCallback: requires a single chain");
  if (!Rf_isFunction(callbackExpr)) Rf_error("callback must be a function");

  size_t numBurnIn = static_cast<size_t>(Rf_asInteger(numBurnInExpr));
  size_t numSamples = static_cast<size_t>(Rf_asInteger(numSamplesExpr));

  SEXP sigmaExpr = rc_getListElement(resultsExpr, "sigma");
  SEXP trainExpr = rc_getListElement(resultsExpr, "train");
  SEXP testExpr = rc_getListElement(resultsExpr, "test");
  SEXP varcountExpr = rc_getListElement(resultsExpr, "varcount");
  SEXP kExpr = rc_getListElement(resultsExpr, "k");
  SEXP varprobsExpr = rc_getListElement(resultsExpr, "varprobs");

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
  results.numReportedLocations = shape.numReportedLocations;
  // rbart_vi's caller-owned varcount buffer is single-slab (R/rbart.R sizes it
  // numPredictors x n.samples), so the count is PINNED to 1 here rather than
  // read off the shape: the layout is then true at this site whatever a future
  // slice widens upstream, and the two guards that keep a multi-forest sampler
  // off this path - setOffset's BCF refusal on the R-loop path, and the grouped
  // refusal in R/spec.R on the in-core one - stop being load-bearing for
  // MEMORY SAFETY. A multi-forest sampler that ever reached here would report
  // its prognostic forest, exactly as the flat C API does.
  results.numVariableCountForests = 1;

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
  bartcore::SamplerShape shape = holder.sampler->shape();
  // constant leaf only in v1 (the scan's closed-form marginal); the R-level
  // method refuses first, this is the engine-side backstop
  if (shape.numLeafCovariates != 0 || shape.usesFunctionLeaves)
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
  bartcore::SamplerShape shape = holder.sampler->shape();
  int updateScale = Rf_asLogical(updateScaleExpr);
  refuseMultiForestResponseMutation(*holder.sampler, "bartcore_setOffset",
                                    ResponseConduit::offset, updateScale);
  refuseVarianceForestScaleUpdate(*holder.sampler, "bartcore_setOffset",
                                  ResponseConduit::offset, updateScale);
  refuseGroupedScaleUpdate(*holder.sampler, "bartcore_setOffset",
                           ResponseConduit::offset, updateScale);
  if (!Rf_isNull(offsetExpr) &&
      (!Rf_isReal(offsetExpr) ||
       static_cast<size_t>(Rf_xlength(offsetExpr)) != shape.numObservations))
    Rf_error("length of replacement offset is not equal to number of observations");
  const double* offset = Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr);
  holder.sampler->setOffset(offset, updateScale == TRUE);
  retain(ptrExpr, PROT_OFFSET, offsetExpr);
  return R_NilValue;
}

SEXP bartcore_setResponse(SEXP ptrExpr, SEXP yExpr, SEXP updateScaleExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  int updateScale = Rf_asLogical(updateScaleExpr);
  refuseMultiForestResponseMutation(*holder.sampler, "bartcore_setResponse",
                                    ResponseConduit::response, updateScale);
  refuseVarianceForestScaleUpdate(*holder.sampler, "bartcore_setResponse",
                                  ResponseConduit::response, updateScale);
  refuseGroupedScaleUpdate(*holder.sampler, "bartcore_setResponse",
                           ResponseConduit::response, updateScale);
  if (!Rf_isReal(yExpr) ||
      static_cast<size_t>(Rf_xlength(yExpr)) != shape.numObservations)
    Rf_error("y must be of length equal to %lu",
             static_cast<unsigned long>(shape.numObservations));
  // support before install: the engine's latent refresh consumes y immediately
  validateResponseSupport(shape.family, shape.numOrdinalThresholds + 1,
                          REAL(yExpr), shape.numObservations,
                          "bartcore_setResponse");
  GetRNGstate(); // probit latent redraw
  holder.sampler->setResponse(REAL(yExpr), updateScale == TRUE);
  PutRNGstate();
  retain(ptrExpr, PROT_RESPONSE, yExpr);
  return R_NilValue;
}

SEXP bartcore_setSigma(SEXP ptrExpr, SEXP sigmaExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refusePinnedSigmaChange(*holder.sampler, "bartcore_setSigma");
  holder.sampler->setSigma(Rf_asReal(sigmaExpr));
  return R_NilValue;
}

SEXP bartcore_setData(SEXP ptrExpr, SEXP dataExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  bartcore::SamplerShape shape = sampler.shape();
  refusePredictorMutation(sampler, "bartcore_setData");
  refuseMultiForestMutation(sampler, "bartcore_setData");
  if (shape.numGroups > 0)
    Rf_error("grouped random effects fix the data at creation; make a new "
             "sampler instead");
  if (shape.family == bartcore::ResponseFamily::aft)
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

    if (data.numPredictors != shape.numPredictors)
      Rf_error("bartcore setData requires the same predictors");
    // the whole-data conduit carries the weight policy too, and had NO
    // backstop of its own: it feeds LogisticResponse::setData's cold start
    // directly, where a zero or negative count becomes a phantom row carrying
    // a full PG(1, psi) precision. Replacement data given WITHOUT weights is
    // unweighted, exactly as at creation and as it is for gaussian - for
    // logistic that reads as single-trial rows, unit counts.
    if (data.weights != NULL) {
      refuseBinaryWeightChange(sampler);
      enforceBinaryWeightPolicy(shape.family, data.weights,
                                data.numObservations);
    }
    // the whole-data conduit swaps y too, so it carries the same support rule
    validateResponseSupport(shape.family, shape.numOrdinalThresholds + 1,
                            data.y, data.numObservations,
                            "bartcore setData");
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
              j, data.predictors.denseValues[i + j * data.numObservations]))
          Rf_error("categorical predictor values must be existing category "
                   "codes");
      for (size_t i = 0; i < data.numTestObservations; ++i)
        if (!sampler.data().categoricalValueIsValid(
              j,
              data.testPredictors.denseValues[i +
                                              j * data.numTestObservations]))
          Rf_error("categorical predictor values must be existing category "
                   "codes");
    }

    sampler.setData(data.predictors.denseValues, data.y, data.numObservations,
                    data.weights, data.offset,
                    data.testPredictors.denseValues, data.numTestObservations,
                    data.testOffset);
    // a family whose augmentation is STATED against the counts takes them
    // through the weight conduit as well, so the latents are drawn against the
    // counts now in force rather than left at the deterministic cold start a
    // data swap installs - a null weight vector being unit counts. Without
    // this a sampler created with counts and handed weightless data would
    // carry omega = 1/4 into the next sweep's tree moves.
    if (shape.family == bartcore::ResponseFamily::logistic)
      sampler.setWeights(data.weights);

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
  bartcore::SamplerShape shape = holder.sampler->shape();
  if (Rf_isNull(xTestExpr)) {
    // removal: back to the no-test-data state, both offsets included
    refuseStaleCategoryTestOffset(!holder.ownedCategoryTestOffset.empty(),
                                  "bartcore_setTestPredictor");
    holder.sampler->setTestPredictors(NULL, 0);
    holder.sampler->setTestOffset(NULL);
    retain(ptrExpr, PROT_TEST_OFFSET, R_NilValue);
    return R_NilValue;
  }
  refuseUndefinedTestFits(*holder.sampler, "bartcore_setTestPredictor");
  refuseStaleCategoryTestOffset(!holder.ownedCategoryTestOffset.empty(),
                                "bartcore_setTestPredictor");
  if (Rf_inherits(xTestExpr, "dbartsMixedMatrix"))
    // whole-object replacement by a mixed/sparse container: parse against the
    // training cut grid, then rebuild the typed test store. Both the parse and
    // the leaf-covariate refusal precede any store change, so a rejected
    // container leaves the sampler in its prior state. The parsed sources own
    // buffers, so the unwind-protected scope frees them on the error jump.
    return unwindProtect([&, parsed = ParsedTestContainer{}]() mutable -> SEXP {
      parseTestContainer(parsed, xTestExpr, shape.numPredictors,
                         holder.sampler->data().types.data());
      validateTestContainerAgainstStore(holder.sampler->data(), parsed.view);
      if (holder.sampler->data().testOffset != NULL &&
          parsed.view.numRows != shape.numTestObservations)
        Rf_error("test offset length would no longer match; set the predictors "
                 "and offset together");
      if (!installTestContainer(*holder.sampler, parsed))
        Rf_error("a leaf covariate column cannot be a sparse test column; "
                 "supply it as a dense test column");
      return R_NilValue;
    });
  SEXP dims = Rf_getAttrib(xTestExpr, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[1]) != shape.numPredictors)
    Rf_error("bartcore_setTestPredictor: requires a matrix with matching columns");
  size_t numTestObservations = static_cast<size_t>(INTEGER(dims)[0]);
  if (holder.sampler->data().testOffset != NULL &&
      numTestObservations != shape.numTestObservations)
    Rf_error("test offset length would no longer match; set the predictors "
             "and offset together");
  for (size_t j = 0; j < shape.numPredictors; ++j)
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
  refuseUndefinedTestFits(*holder.sampler, "bartcore_setTestOffset");
  refuseMultiForestTestOffset(*holder.sampler, "bartcore_setTestOffset");
  size_t numTestObservations = holder.sampler->shape().numTestObservations;
  if (numTestObservations == 0)
    Rf_error("cannot set a test offset without test predictors");
  if (!Rf_isReal(offsetExpr) ||
      static_cast<size_t>(Rf_xlength(offsetExpr)) != numTestObservations)
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
  size_t numPredictors = holder.sampler->shape().numPredictors;
  if (Rf_isNull(xTestExpr)) {
    if (!Rf_isNull(offsetExpr))
      Rf_error("when test matrix is NULL, test offset must be as well");
    refuseStaleCategoryTestOffset(!holder.ownedCategoryTestOffset.empty(),
                                  "bartcore_setTestPredictorAndOffset");
    holder.sampler->setTestPredictors(NULL, 0);
    holder.sampler->setTestOffset(NULL);
    retain(ptrExpr, PROT_TEST_OFFSET, R_NilValue);
    return R_NilValue;
  }
  refuseUndefinedTestFits(*holder.sampler,
                          "bartcore_setTestPredictorAndOffset");
  refuseStaleCategoryTestOffset(!holder.ownedCategoryTestOffset.empty(),
                                "bartcore_setTestPredictorAndOffset");
  if (!Rf_isNull(offsetExpr))
    refuseMultiForestTestOffset(*holder.sampler,
                                "bartcore_setTestPredictorAndOffset");
  if (Rf_inherits(xTestExpr, "dbartsMixedMatrix"))
    // whole-object container replacement plus its offset; parse and the
    // leaf-covariate refusal precede the store change, so a rejected container
    // leaves the sampler untouched
    return unwindProtect([&, parsed = ParsedTestContainer{}]() mutable -> SEXP {
      parseTestContainer(parsed, xTestExpr, numPredictors,
                         holder.sampler->data().types.data());
      validateTestContainerAgainstStore(holder.sampler->data(), parsed.view);
      if (!Rf_isNull(offsetExpr) &&
          (!Rf_isReal(offsetExpr) ||
           static_cast<size_t>(Rf_xlength(offsetExpr)) != parsed.view.numRows))
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
      static_cast<size_t>(INTEGER(dims)[1]) != numPredictors)
    Rf_error("bartcore_setTestPredictorAndOffset: requires a matrix with "
             "matching columns");
  size_t numTestObservations = static_cast<size_t>(INTEGER(dims)[0]);
  if (!Rf_isNull(offsetExpr) &&
      (!Rf_isReal(offsetExpr) ||
       static_cast<size_t>(Rf_xlength(offsetExpr)) != numTestObservations))
    Rf_error("length of test offset must equal number of rows in test matrix");
  for (size_t j = 0; j < numPredictors; ++j)
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

// Case weights: a pointer swap with nothing rescaled for gaussian, and for
// logistic a model change - the counts are the Polya-Gamma shape, so the
// engine redraws the latents against them off the chain generator, never R's
// stream, which is why no GetRNGstate bracket appears here.
// refuseBinaryWeightChange refuses the families whose likelihood has no weight
// slot at all.
SEXP bartcore_setWeights(SEXP ptrExpr, SEXP weightsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t numObservations = shape.numObservations;
  refuseMultiForestResponseMutation(*holder.sampler, "bartcore_setWeights",
                                    ResponseConduit::weights, FALSE);
  refuseBinaryWeightChange(*holder.sampler);
  if (!Rf_isReal(weightsExpr) ||
      static_cast<size_t>(Rf_xlength(weightsExpr)) != numObservations)
    Rf_error("length of weights must equal number of observations");
  const double* weights = REAL(weightsExpr);
  // the same values creation accepts, or the swap walks the sampler off its
  // family's weight policy: a logistic count that is zero, negative or
  // fractional is silently rounded by the PG draw's lround and leaves a row
  // carrying a full PG(1, psi) precision it has no observation for
  enforceBinaryWeightPolicy(shape.family, weights, numObservations);
  holder.sampler->setWeights(weights);
  retain(ptrExpr, PROT_WEIGHTS, weightsExpr);
  return R_NilValue;
}

// Between-run reconfiguration. The R side refuses changes to the engine, rng,
// and cut settings; chain and tree counts shape live storage, so they are
// re-checked here.
SEXP bartcore_setControl(SEXP ptrExpr, SEXP controlExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  bartcore::SamplerShape shape = sampler.shape();

  if (!Rf_inherits(controlExpr, "dbartsControl"))
    Rf_error("'control' argument to bartcore_setControl not of class "
             "'dbartsControl'");

  ParsedControl control;
  parseControl(control, controlExpr);

  if (control.numChains != shape.numChains)
    Rf_error("the bartcore engine cannot change the number of chains of an "
             "existing sampler");
  if (control.numTrees != shape.numTrees)
    Rf_error("the bartcore engine cannot change the number of trees of an "
             "existing sampler");

  holder.keepTrainingFits = control.keepTrainingFits;
  sampler.setNumThreads(control.numThreads);
  sampler.setNumThin(control.treeThinningRate);
  sampler.setVerbose(control.verbose, control.printEvery);
  sampler.setTreeStorage(control.keepTrees, control.defaultNumSamples);

  return R_NilValue;
}

/// Prior replacement; installing a model before any run matches creating
/// with it.
SEXP bartcore_setModel(SEXP ptrExpr, SEXP modelExpr, SEXP dataExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  bartcore::SamplerShape shape = sampler.shape();
  refuseMultiForestMutation(sampler, "bartcore_setModel");

  if (!Rf_inherits(modelExpr, "dbartsModel"))
    Rf_error("'model' argument to bartcore_setModel not of class "
             "'dbartsModel'");

  return unwindProtect([&, model = ParsedModel{}]() mutable -> SEXP {
    parseModel(model, modelExpr, shape.numPredictors);

    // the leaf model is a template instantiation: the designation and its
    // kind are fixed at creation, so a replacement prior must carry the same
    bool designationMatches =
      model.leafCovariateColumns.size() == shape.numLeafCovariates &&
      model.gpLeaves == shape.usesFunctionLeaves;
    for (size_t j = 0; designationMatches &&
         j < model.leafCovariateColumns.size(); ++j)
      designationMatches =
        model.leafCovariateColumns[j] == shape.leafCovariateColumns[j];
    if (!designationMatches)
      Rf_error("the leaf covariate designation is fixed when a sampler is "
               "created; make a new sampler instead");

    // gaussian and aft both draw sigma from a variance prior; the binary
    // families fix it, so setModel leaves their sigma alone
    bool drawsSigma = shape.family == bartcore::ResponseFamily::gaussian ||
                      shape.family == bartcore::ResponseFamily::aft;

    bartcore::ModelParameters parameters;
    parameters.base = model.base;
    parameters.power = model.power;
    parameters.splitProbabilities = model.splitProbabilities;
    parameters.birthOrDeathProbability = model.birthOrDeathProbability;
    parameters.swapProbability = model.swapProbability;
    parameters.changeProbability = model.changeProbability;
    parameters.birthProbability = model.birthProbability;
    parameters.nodeScale = model.nodeScale;
    // carried so the install re-derives the named calibration against the
    // transform in force; without it $setModel(sampler$model) - a no-op round
    // trip - would revert to the family-keyed node scale
    parameters.priorScale = model.priorScale;
    parameters.updateK = model.updateK;
    if (parameters.updateK) {
      parameters.kHyperprior.degreesOfFreedom = model.kDf;
      parameters.kHyperprior.scale = model.kScale;
    } else {
      parameters.k = model.k;
    }
    if (drawsSigma) {
      // the prior half travels in both arms: a homoscedastic chain reads it
      // only when sigma is drawn, but a heteroscedastic one recalibrates its
      // scale leaf from the whole triple whatever fixed() says about sigma
      parameters.sigmaDf = model.sigmaDf;
      parameters.sigmaRawScale = model.sigmaRawScale;
      if (model.sigmaIsFixed) {
        // documented semantics: fixed(value) holds the residual variance
        parameters.sigmaIsFixed = true;
        parameters.sigmaEstimate = std::sqrt(model.fixedSigmaSq);
      } else {
        parameters.sigmaEstimate = rc_getDouble(
          Rf_getAttrib(dataExpr, Rf_install("sigma")), "sigma estimate",
          RC_LENGTH | RC_EQ, rc_asRLength(1), RC_NA | RC_YES,
          RC_VALUE | RC_GT, 0.0, RC_END);
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
  size_t numChains = holder.sampler->shape().numChains;
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    REAL(result)[c] = holder.sampler->sigma(c);
  UNPROTECT(1);
  return result;
}

// The dispersion r in force on each chain - the mid-sweep read of the scalar
// the recorded dispersion channel stores once per kept draw, without
// serializing state. NULL, rather than a filler value, off a family that
// carries no dispersion: the caller then tests for the channel instead of
// comparing against a number that would mean nothing.
SEXP bartcore_getDispersion(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerShape shape = holder.sampler->shape();
  if (!shape.carriesDispersion) return R_NilValue;
  size_t numChains = shape.numChains;
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    REAL(result)[c] = holder.sampler->dispersion(c);
  UNPROTECT(1);
  return result;
}

SEXP bartcore_getSumsOfSquaredResiduals(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numChains = holder.sampler->shape().numChains;
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    REAL(result)[c] = holder.sampler->sumOfSquaredResiduals(c);
  UNPROTECT(1);
  return result;
}

// Predictor mutation. The sampler borrows a full replacement matrix
// (R/bartcore.R retains it on success); column and per-observation updates
// write in place into the matrix the sampler currently borrows, so they
// alias the R-side data.

SEXP bartcore_setPredictor(SEXP ptrExpr, SEXP xExpr, SEXP forceUpdateExpr,
                           SEXP updateCutPointsExpr) {
  // the materialized block is owned across validateColumnValues, whose refusal
  // longjmps past its destructor
  return unwindProtect([&, parsed = ParsedMutationSource{}]() mutable -> SEXP {
    const char* shapeMessage =
      "bartcore_setPredictor requires a matrix with matching dimensions";
    BartcoreHolder& holder(holderFromExpression(ptrExpr));
    bartcore::SamplerShape shape = holder.sampler->shape();
    // no shape guard left: the two-phase transaction covers every forest and
    // the variance forest, so an unforced call vetoes or rolls back rather
    // than accepting a change one ensemble would misroute
    refuseMutationOnView(*holder.sampler, "bartcore_setPredictor");
    const double* values;
    if (Rf_isReal(xExpr)) {
      SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
      if (Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
          static_cast<size_t>(INTEGER(dims)[0]) != shape.numObservations ||
          static_cast<size_t>(INTEGER(dims)[1]) != shape.numPredictors)
        Rf_error("%s", shapeMessage);
      values = REAL(xExpr);
    } else {
      if (!parseMutationSource(parsed, xExpr, shape.numObservations,
                               shape.numPredictors, shapeMessage))
        Rf_error("%s", shapeMessage);
      values = materializeMutationSource(parsed,
                                         holder.sampler->data().types.data());
    }

    for (size_t j = 0; j < shape.numPredictors; ++j)
      validateColumnValues(holder.sampler->data(), j,
                           values + j * shape.numObservations,
                           shape.numObservations);

    bartcore::PredictorUpdateResult result = holder.sampler->setPredictor(
      values, Rf_asLogical(forceUpdateExpr) == TRUE,
      Rf_asLogical(updateCutPointsExpr) == TRUE);
    if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
      Rf_error("number of induced cut points in new predictor less than "
               "previous: old splits would be invalid");
    // the engine re-quantized the values into owned codes and retains no
    // pointer; the R method reassigns sampler$data@x, which keeps the current
    // source alive
    return Rf_ScalarLogical(
      result == bartcore::PredictorUpdateResult::accepted ? TRUE : FALSE);
  });
}

SEXP bartcore_updatePredictor(SEXP ptrExpr, SEXP xExpr, SEXP columnsExpr,
                              SEXP forceUpdateExpr, SEXP updateCutPointsExpr) {
  return unwindProtect([&, columns = std::vector<size_t>{},
                        parsed = ParsedMutationSource{}]() mutable -> SEXP {
    const char* shapeMessage =
      "bartcore_updatePredictor requires numObservations values per column";
    BartcoreHolder& holder(holderFromExpression(ptrExpr));
    // unguarded, as bartcore_setPredictor above
    refuseMutationOnView(*holder.sampler, "bartcore_updatePredictor");
    bartcore::SamplerShape shape = holder.sampler->shape();
    size_t numObservations = shape.numObservations;
    size_t numPredictors = shape.numPredictors;

    size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
    if (numColumns == 0 ||
        (Rf_isReal(xExpr) && static_cast<size_t>(Rf_xlength(xExpr)) !=
                               numObservations * numColumns))
      Rf_error("%s", shapeMessage);

    columns.resize(numColumns);
    for (size_t k = 0; k < numColumns; ++k) {
      int column = INTEGER(columnsExpr)[k];
      if (column < 1 || static_cast<size_t>(column) > numPredictors)
        Rf_error("bartcore_updatePredictor: column out of range");
      columns[k] = static_cast<size_t>(column - 1);
    }

    const double* values;
    if (Rf_isReal(xExpr)) {
      values = REAL(xExpr);
    } else {
      if (!parseMutationSource(parsed, xExpr, numObservations, numColumns,
                               shapeMessage))
        Rf_error("%s", shapeMessage);
      // the store's type of each column the argument names, in argument order
      parsed.storeTypes.resize(numColumns);
      for (size_t k = 0; k < numColumns; ++k)
        parsed.storeTypes[k] = holder.sampler->data().types[columns[k]];
      values = materializeMutationSource(parsed, parsed.storeTypes.data());
    }

    for (size_t k = 0; k < numColumns; ++k)
      validateColumnValues(holder.sampler->data(), columns[k],
                           values + k * numObservations, numObservations);

    bartcore::PredictorUpdateResult result = holder.sampler->updatePredictor(
      values, columns.data(), numColumns,
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
    refuseMutationOnView(*holder.sampler, "bartcore_setCutPoints");
    size_t numPredictors = holder.sampler->shape().numPredictors;
    // dense columns re-quantize from the supplied data@x; CSC/mixed columns
    // read their retained slices, so a non-matrix source is passed as null
    const double* currentPredictors =
      Rf_isReal(currentPredictorsExpr) ? REAL(currentPredictorsExpr) : NULL;

    size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
    if (numColumns == 0 ||
        static_cast<size_t>(Rf_xlength(cutPointsExpr)) != numColumns)
      Rf_error("bartcore_setCutPoints: requires one cut point vector per column");

    cutPoints.resize(numColumns);
    numCutPoints.resize(numColumns);
    columns.resize(numColumns);
    for (size_t k = 0; k < numColumns; ++k) {
      int column = INTEGER(columnsExpr)[k];
      if (column < 1 || static_cast<size_t>(column) > numPredictors)
        Rf_error("bartcore_setCutPoints: column out of range");
      columns[k] = static_cast<size_t>(column - 1);
      if (holder.sampler->data().types[columns[k]] ==
          bartcore::ColumnType::categorical)
        Rf_error("cannot set cut points for a categorical predictor");

      SEXP cutsExpr = VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(k));
      R_xlen_t numCuts = Rf_xlength(cutsExpr);
      if (numCuts > 65535)  // codes must fit xint_t, including numCuts itself
        Rf_error("bartcore_setCutPoints: cut point vector too long");
      // an ordinal column with no cut point is not a state the store can hold:
      // its own validator refuses one, so the sampler could not restore itself
      if (numCuts < 1)
        Rf_error("bartcore_setCutPoints: requires at least one cut point per "
                 "column");
      const double* cuts = REAL(cutsExpr);
      for (R_xlen_t i = 1; i < numCuts; ++i)
        if (cuts[i] <= cuts[i - 1])
          Rf_error("bartcore_setCutPoints: requires strictly increasing cut "
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
  bartcore::SamplerShape shape = holder.sampler->shape();
  // the session's cell guard caches every forest AND the variance forest,
  // pruned to the trees the column can move, so every sampler shape takes this
  // entry: a row installs only where no leaf of any tree would empty
  refuseMutationOnView(
    *holder.sampler, "bartcore_updatePredictorPerObservation");
  size_t numObservations = shape.numObservations;

  if (static_cast<size_t>(Rf_xlength(xExpr)) != numObservations)
    Rf_error("bartcore_updatePredictorPerObservation: requires one value per "
             "observation");
  int column = Rf_asInteger(columnExpr);
  if (column < 1 || static_cast<size_t>(column) > shape.numPredictors)
    Rf_error("bartcore_updatePredictorPerObservation: column out of range");
  // per-observation replacement writes one cell at a time, which a sparse
  // column's rank storage cannot take without an O(nnz) shift per cell; the
  // sparse mutation path is whole-column, so a CSC-backed target is refused by
  // name (replace it wholesale with updatePredictor). Dense-backed columns of
  // a mixed store - the IRT latent-in-a-sparse-design case - stay open.
  if (holder.sampler->data().columnIsCscBacked(static_cast<size_t>(column - 1)))
    Rf_error("bartcore_updatePredictorPerObservation: per-observation updates "
             "require a dense-backed column; a sparse column fixes its nonzero "
             "pattern per cell - replace the whole column with updatePredictor");
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
      Rf_error("bartcore_updatePredictorPerObservationJointly: requires one "
               "column per sampler");

    samplers.resize(numSamplers);
    columns.resize(numSamplers);
    for (size_t k = 0; k < numSamplers; ++k) {
      BartcoreHolder& holder(
        holderFromExpression(VECTOR_ELT(ptrsExpr, static_cast<R_xlen_t>(k))));
      // unguarded, as the single-sampler entry above
      refuseMutationOnView(
        *holder.sampler, "bartcore_updatePredictorPerObservationJointly");
      samplers[k] = holder.sampler.get();
      int column = INTEGER(columnsExpr)[k];
      if (column < 1 ||
          static_cast<size_t>(column) > samplers[k]->shape().numPredictors)
        Rf_error("bartcore_updatePredictorPerObservationJointly: column out of "
                 "range");
      // per-observation cell writes need a dense-backed target (see the
      // single-sampler entry point)
      if (samplers[k]->data().columnIsCscBacked(static_cast<size_t>(column - 1)))
        Rf_error("bartcore_updatePredictorPerObservationJointly: per-"
                 "observation updates require a dense-backed column; replace a "
                 "sparse column wholesale with updatePredictor");
      columns[k] = static_cast<size_t>(column - 1);
    }

    size_t numObservations = samplers[0]->shape().numObservations;
    for (size_t k = 1; k < numSamplers; ++k)
      if (samplers[k]->shape().numObservations != numObservations)
        Rf_error("bartcore_updatePredictorPerObservationJointly: requires "
                 "index-aligned samplers");
    if (static_cast<size_t>(Rf_xlength(xExpr)) != numObservations)
      Rf_error("bartcore_updatePredictorPerObservationJointly: requires one "
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

// State serialization. The returned object is engine-specific and opaque: one
// list per chain (flattened live and saved trees with each forest's k and leaf
// scale, original-scale sigma, the response transform at capture, latents and
// their per-family companions, group, glue and dart state, and the serialized
// rng), with the store's cut points, the saved-tree write position, the draws
// it is read against, and the format version as attributes. Restore reinstalls
// the captured transform - a scale setOffset(updateScale) moved after creation
// survives - and rebuilds the rest canonically, partitions and fits from the
// trees, so a restored sampler continues EQUIVALENTLY, not bitwise.

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
      // 1.0e8 is a defensive sanity cap so a corrupt count cannot drive the
      // width computation below into an absurd allocation
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
  refuseMutationOnView(*holder.sampler, "bartcore_setState");
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
  refuseMutationOnView(*holder.sampler, "bartcore_installForests");
  bartcore_bridge::installForests(*holder.sampler, donorStateExpr, samplesExpr);
  return R_NilValue;
}

// The shared tail of both predict entrances: allocate the surface's result
// shape, replay the borrowed view into it, and add the offset. The view's rows
// are the result's rows either flavor, so the dense and resident-sparse paths
// differ only in how the view was built. numThreads is the replay's
// per-call worker count; the offset add and the variance clone below are
// serial passes over the same output, so the bridge's own serial share is one
// pass against the replay's numTrees.
static SEXP predictFromSource(bartcore::SamplerBase& sampler,
                              const bartcore::SamplerShape& shape,
                              const bartcore::PredictorSource& source,
                              SEXP offsetExpr, size_t numThreads) {
  size_t numTestObservations = source.numRows;
  if (numTestObservations == 0) Rf_error("bartcore_predict: requires rows");

  refuseEmptyTreeStore(sampler, "bartcore_predict");

  size_t capacity = shape.savedTreeCapacity;
  size_t numChains = shape.numChains;
  // the draw axis is what the store has RECORDED, oldest first, not its
  // capacity: a run short of capacity reports the draws it made
  size_t numSamples = capacity > 0 ? shape.numSavedDraws : 1;
  size_t numLocations = shape.numReportedLocations;

  // The offset takes the shape of the surface. A multi-location (multinomial
  // softmax) surface reports probabilities, not an additive latent scale, so
  // the offset is a PER-CATEGORY nNew x K matrix entering the raw fits before
  // the blend - the only place a shift means anything under a softmax. A flat
  // vector stays refused there, and truthfully: after the blend it would move
  // the values off the simplex, and before it a common per-observation shift is
  // the softmax's own null direction. Every other surface takes the flat form,
  // added to the replayed fits below.
  const double* offset = NULL;
  const double* categoryOffset = NULL;
  if (!Rf_isNull(offsetExpr)) {
    if (numLocations > 1) {
      if (Rf_isNull(Rf_getAttrib(offsetExpr, R_DimSymbol)))
        Rf_error("bartcore_predict: a flat offset is undefined for a "
                 "multi-location (multinomial softmax) predict surface, whose "
                 "offset must be a per-category matrix, one row per predicted "
                 "row");
      categoryOffset = validateCategoryOffset(offsetExpr, numTestObservations,
                                              numLocations,
                                              "bartcore_predict offset");
    } else {
      if (!Rf_isReal(offsetExpr) ||
          static_cast<size_t>(Rf_xlength(offsetExpr)) != numTestObservations)
        Rf_error("bartcore_predict: offset must have one value per row");
      offset = REAL(offsetExpr);
    }
  }

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
                               static_cast<int>(numSamples)))
      : PROTECT(Rf_alloc3DArray(REALSXP, numTestInt,
                                static_cast<int>(numSamples),
                                static_cast<int>(numChains)));
  } else {
    resultExpr = numChains == 1
      ? PROTECT(Rf_allocVector(REALSXP,
                               static_cast<R_xlen_t>(numTestObservations)))
      : PROTECT(Rf_allocMatrix(REALSXP, numTestInt,
                               static_cast<int>(numChains)));
  }

  sampler.predict(source, numTestObservations, categoryOffset, numThreads,
                  REAL(resultExpr));

  if (offset != NULL) {
    for (size_t slab = 0; slab < numSamples * numChains; ++slab)
      misc_addVectorsInPlace(offset, numTestObservations,
                             REAL(resultExpr) + slab * numTestObservations);
  }

  // heteroscedastic: the variance surface s^2(x) rides alongside f(x) as a
  // named list (mean, variance), same shape (gaussian, single location). Off a
  // variance forest the bare mean array returns, backward-compatible. Predict
  // needs saved trees, so a null-capacity variance forest has nothing to replay.
  if (shape.hasVarianceForest && capacity > 0) {
    SEXP varianceExpr = PROTECT(Rf_duplicate(resultExpr));  // clone the shape
    sampler.predictVariance(source, numTestObservations, numThreads,
                            REAL(varianceExpr));
    SEXP listExpr = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(listExpr, 0, resultExpr);
    SET_VECTOR_ELT(listExpr, 1, varianceExpr);
    SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(namesExpr, 0, Rf_mkChar("mean"));
    SET_STRING_ELT(namesExpr, 1, Rf_mkChar("variance"));
    Rf_setAttrib(listExpr, R_NamesSymbol, namesExpr);
    UNPROTECT(4);  // resultExpr, varianceExpr, listExpr, namesExpr
    return listExpr;
  }

  UNPROTECT(1);
  return resultExpr;
}

// Fits for new data on the original response scale (a binary response gives
// the latent scale). With keepTrees the saved trees produce
// numTestObservations x numSamples (x numChains) fits; without, the live
// trees produce a single set per chain. A multi-location combiner
// (multinomial: K softmax channels) inserts the K dimension between the rows
// and the samples, matching the run's test channel. offset, when non-null, is
// added to every sample's fits - or, on a multi-location surface, is the
// nNew x K matrix each replay adds to its raw fits before the softmax, since
// probabilities admit no additive shift. A dense matrix, a dgCMatrix, and a
// dbartsMixedMatrix all predict; only the dense one is materialized, and it was
// already.
SEXP bartcore_predict(SEXP ptrExpr, SEXP xTestExpr, SEXP offsetExpr,
                      SEXP nThreadsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  bartcore::SamplerShape shape = sampler.shape();

  // a per-call worker count, not a stored one; the engine partitions the
  // replay by (chain, draw) and reports the same numbers at every value
  size_t numThreads = static_cast<size_t>(
    rc_getInt(nThreadsExpr, "number of threads", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_NA | RC_NO, RC_END));

  // predict() sums only the first forest, so an amplitude coupling's
  // prediction would drop every other forest and the glue; refuse it for the
  // same reason its recorded test fits are undefined.
  refuseUndefinedTestFits(sampler, "bartcore_predict");
  // The rows are the CALLER's, so the offset must be too: predict never reads
  // the sampler's resident category offsets, whose rows are other rows, and a
  // row-count coincidence between them is not consent. A sampler that models a
  // per-category shift therefore cannot answer here without one being named -
  // and an all-zero matrix is how a caller asks for the offset-free surface.
  // Either resident offset triggers the refusal: a train-only offset predicts
  // the offset-free surface just as wrongly as a test-only one would.
  if ((!holder.ownedCategoryOffset.empty() ||
       !holder.ownedCategoryTestOffset.empty()) && Rf_isNull(offsetExpr))
    Rf_error("bartcore_predict: this sampler carries an n x K category offset, "
             "and the predicted rows are not its rows, so their offset cannot "
             "be inferred; pass one per predicted row (an all-zero matrix for "
             "the offset-free surface)");

  // a sparse test set replays resident, off the view the container parses to;
  // the parse owns buffers, so the unwind-protected scope frees them on the
  // error jump
  if (testSourceIsSparse(xTestExpr))
    return unwindProtect([&, parsed = ParsedTestContainer{}]() mutable -> SEXP {
      parseTestSource(parsed, xTestExpr, shape.numPredictors,
                      sampler.data().types.data());
      validateTestContainerAgainstStore(sampler.data(), parsed.view);
      refuseSparseLeafCovariate(shape, parsed.view);
      return predictFromSource(sampler, shape, parsed.view, offsetExpr,
                               numThreads);
    });

  size_t numTestObservations =
    validatePredictorMatrix(sampler, xTestExpr, "bartcore_predict");
  return predictFromSource(
    sampler, shape,
    bartcore::densePredictorSource(REAL(xTestExpr), numTestObservations,
                                   shape.numPredictors),
    offsetExpr, numThreads);
}

// The allocation/replay tail of bartcore_predictPerForest. It cannot reuse
// predictFromSource: the forest margin changes the result's dim vector, and
// there is no offset to add. An amplitude coupling's location is
// shift + sum_f (B_f a_f) * (scale * f_f), so the response transform AND the
// offset both belong to the combination, never to one forest's own total.
//
// The result is numTestObservations x numForests x numSamples (x numChains) -
// the layout the run's own per-forest channel carries, with numSamples the
// saved capacity, or 1 off keepTrees, where the live trees replay instead.
static SEXP predictPerForestFromSource(bartcore::SamplerBase& sampler,
                                       const bartcore::SamplerShape& shape,
                                       const bartcore::PredictorSource& source,
                                       size_t numThreads) {
  size_t numTestObservations = source.numRows;
  if (numTestObservations == 0)
    Rf_error("bartcore_predictPerForest: requires rows");

  refuseEmptyTreeStore(sampler, "bartcore_predictPerForest");

  size_t capacity = shape.savedTreeCapacity;
  size_t numSamples = capacity > 0 ? shape.numSavedDraws : 1;

  int dims[4];
  int numDims = 0;
  dims[numDims++] = static_cast<int>(numTestObservations);
  dims[numDims++] = static_cast<int>(shape.numForests);
  dims[numDims++] = static_cast<int>(numSamples);
  if (shape.numChains > 1) dims[numDims++] = static_cast<int>(shape.numChains);
  R_xlen_t total = 1;
  for (int d = 0; d < numDims; ++d) total *= dims[d];
  SEXP resultExpr = PROTECT(Rf_allocVector(REALSXP, total));
  SEXP dimExpr = Rf_allocVector(INTSXP, numDims);
  for (int d = 0; d < numDims; ++d) INTEGER(dimExpr)[d] = dims[d];
  Rf_setAttrib(resultExpr, R_DimSymbol, dimExpr);

  sampler.predictPerForest(source, numTestObservations, numThreads,
                           REAL(resultExpr));
  UNPROTECT(1);
  return resultExpr;
}

// Each forest's own RAW fits for new data, on the forests' INTERNAL scale and
// without any offset - the off-sample twin of bartcore_getForestFits. Gated on
// forestReportingIsDefined, the same predicate the run's forestFits channel
// exists under, so an ordinary single-forest sampler and a multinomial one are
// both refused: the latter's raw per-category fits ARE defined, but its
// reported quantity is a softmax probability and its off-sample surface is
// bartcore_predict, which reports that. Sparse test sources replay resident,
// exactly as they do there.
SEXP bartcore_predictPerForest(SEXP ptrExpr, SEXP xTestExpr, SEXP offsetExpr,
                               SEXP nThreadsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  bartcore::SamplerShape shape = sampler.shape();

  size_t numThreads = static_cast<size_t>(
    rc_getInt(nThreadsExpr, "number of threads", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_NA | RC_NO, RC_END));

  if (!shape.forestReportingIsDefined)
    Rf_error("bartcore_predictPerForest: this sampler reports no per-forest "
             "fits; only a coupling that composes its forests through scalar "
             "amplitude glue carries them");
  // A per-forest raw fit takes no offset: an offset shifts the COMBINATION,
  // exactly as the response transform's shift does, and neither is any one
  // forest's. Refused rather than ignored, so a caller who means to shift the
  // recombination is told where the shift belongs.
  if (!Rf_isNull(offsetExpr))
    Rf_error("bartcore_predictPerForest: a per-forest fit takes no offset, "
             "which shifts the recombination rather than any one forest's own "
             "total; add it there instead");

  if (testSourceIsSparse(xTestExpr))
    return unwindProtect([&, parsed = ParsedTestContainer{}]() mutable -> SEXP {
      parseTestSource(parsed, xTestExpr, shape.numPredictors,
                      sampler.data().types.data());
      validateTestContainerAgainstStore(sampler.data(), parsed.view);
      refuseSparseLeafCovariate(shape, parsed.view);
      return predictPerForestFromSource(sampler, shape, parsed.view,
                                        numThreads);
    });

  size_t numTestObservations =
    validatePredictorMatrix(sampler, xTestExpr, "bartcore_predictPerForest");
  return predictPerForestFromSource(
    sampler, shape,
    bartcore::densePredictorSource(REAL(xTestExpr), numTestObservations,
                                   shape.numPredictors),
    numThreads);
}

// The shape of the last out-of-sample replay's fan-out: list(resolved,
// n.workers, worker), the last a per-slab worker index in the same (chain,
// draw) order the replay's output is laid out in. It exists so a test can prove that a thread
// argument reached the engine and that the partition covers every slab exactly
// once, without measuring time; no R surface reads it, and the numbers it
// reports never change an answer.
SEXP bartcore_lastPredictPartition(void) {
  const bartcore::PredictPartitionChannel& channel = bartcore::predictPartition;
  size_t numSlabs = channel.workerForSlab.size();
  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(resultExpr, 0,
                 Rf_ScalarInteger(static_cast<int>(channel.resolvedThreads)));
  SET_VECTOR_ELT(resultExpr, 1,
                 Rf_ScalarInteger(static_cast<int>(channel.numWorkers)));
  SEXP workerExpr =
    PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(numSlabs)));
  for (size_t i = 0; i < numSlabs; ++i)
    INTEGER(workerExpr)[i] = static_cast<int>(channel.workerForSlab[i]);
  SET_VECTOR_ELT(resultExpr, 2, workerExpr);
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("resolved"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("n.workers"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("worker"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);
  UNPROTECT(3);
  return resultExpr;
}

// Test-only companion: replaces the derived traversal cutoff, so a fixture
// small enough to run in a test still fans out and an identity-across-thread-
// counts assertion is not vacuous. 0 restores the derived constant. Returns
// the value it replaced.
SEXP bartcore_setPredictParallelCutoff(SEXP cutoffExpr) {
  size_t previous = bartcore::predictPartition.cutoffOverride;
  bartcore::predictPartition.cutoffOverride = static_cast<size_t>(
    rc_getInt(cutoffExpr, "predict parallel cutoff", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_NO, RC_END));
  return Rf_ScalarInteger(static_cast<int>(previous));
}

// index conversion around bartcore_bridge::getTrees, which describes the
// data.frame produced
SEXP bartcore_getTrees(SEXP ptrExpr, SEXP chainNumsExpr, SEXP sampleNumsExpr,
                       SEXP treeNumsExpr, SEXP currentExpr, SEXP newdataExpr,
                       SEXP trainingDataExpr, SEXP forestExpr) {
  return unwindProtect([&, chainIndices = std::vector<size_t>{},
                        sampleIndices = std::vector<size_t>{},
                        treeIndices = std::vector<size_t>{},
                        parsed = ParsedTestContainer{}]() mutable -> SEXP {
    BartcoreHolder& holder(holderFromExpression(ptrExpr));
    bartcore::SamplerBase& sampler(*holder.sampler);
    bartcore::SamplerShape shape = sampler.shape();

    // forest addressing follows getForestFits: 0-based, unconverted
    size_t forestIndex = static_cast<size_t>(Rf_asInteger(forestExpr));
    if (forestIndex >= shape.numForests)
      Rf_error("bartcore_getTrees: forest index out of range");

    bool useLiveTrees = Rf_asLogical(currentExpr) == TRUE;
    bool useSaved = shape.savedTreeCapacity > 0 && !useLiveTrees;

    chainIndices.resize(static_cast<size_t>(Rf_xlength(chainNumsExpr)));
    for (size_t i = 0; i < chainIndices.size(); ++i) {
      int chainNum = INTEGER(chainNumsExpr)[i];
      if (chainNum < 1) Rf_error("bartcore_getTrees: chain number out of range");
      chainIndices[i] = static_cast<size_t>(chainNum - 1);
    }
    if (useSaved) {
      // a null sample list is every RECORDED draw, oldest first, as for
      // printTrees: the caller cannot count them, the engine can
      if (Rf_isNull(sampleNumsExpr)) {
        for (size_t i = 0; i < shape.numSavedDraws; ++i)
          sampleIndices.push_back(i);
      } else {
        sampleIndices.resize(static_cast<size_t>(Rf_xlength(sampleNumsExpr)));
        for (size_t i = 0; i < sampleIndices.size(); ++i) {
          int sampleNum = INTEGER(sampleNumsExpr)[i];
          if (sampleNum < 1)
            Rf_error("bartcore_getTrees: sample number out of range");
          sampleIndices[i] = static_cast<size_t>(sampleNum - 1);
        }
      }
    }
    treeIndices.resize(static_cast<size_t>(Rf_xlength(treeNumsExpr)));
    for (size_t i = 0; i < treeIndices.size(); ++i) {
      int treeNum = INTEGER(treeNumsExpr)[i];
      if (treeNum < 1) Rf_error("bartcore_getTrees: tree number out of range");
      treeIndices[i] = static_cast<size_t>(treeNum - 1);
    }

    // newdata replays through the trees exactly as predict routes it: dense
    // from its own block, sparse from the container's rank bitmaps
    const bartcore::PredictorSource* newdata = NULL;
    size_t newdataNumRows = 0;
    bartcore::PredictorSource newdataView;
    if (testSourceIsSparse(newdataExpr)) {
      parseTestSource(parsed, newdataExpr, shape.numPredictors,
                      sampler.data().types.data());
      validateTestContainerAgainstStore(sampler.data(), parsed.view);
      refuseSparseLeafCovariate(shape, parsed.view);
      newdata = &parsed.view;
      newdataNumRows = parsed.view.numRows;
    } else if (!Rf_isNull(newdataExpr)) {
      newdataNumRows =
        validatePredictorMatrix(sampler, newdataExpr, "bartcore_getTrees");
      newdataView = bartcore::densePredictorSource(
        REAL(newdataExpr), newdataNumRows, shape.numPredictors);
      newdata = &newdataView;
    }

    // saved-tree replay reads the current training predictors the R method
    // supplies (data@x); the engine keeps no matrix. A dense matrix serves; a
    // sparse/absent source leaves saved-tree counts unpopulated (a view sampler)
    const double* trainingReplay = NULL;
    size_t trainingReplayNumRows = 0;
    if (useSaved && newdata == NULL && Rf_isReal(trainingDataExpr)) {
      SEXP dims = Rf_getAttrib(trainingDataExpr, R_DimSymbol);
      if (!Rf_isNull(dims) && Rf_xlength(dims) == 2 &&
          static_cast<size_t>(INTEGER(dims)[1]) == shape.numPredictors) {
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
    bartcore::SamplerShape shape = sampler.shape();

    size_t capacity = shape.savedTreeCapacity;
    if (capacity > 0) refuseEmptyTreeStore(sampler, "bartcore_printTrees");

    if (Rf_isNull(chainNumsExpr)) {
      for (size_t i = 0; i < shape.numChains; ++i) chainIndices.push_back(i);
    } else {
      for (R_xlen_t i = 0; i < Rf_xlength(chainNumsExpr); ++i) {
        int chainNum = INTEGER(chainNumsExpr)[i];
        if (chainNum < 1 || static_cast<size_t>(chainNum) > shape.numChains)
          Rf_error("bartcore_printTrees: chain number out of range");
        chainIndices.push_back(static_cast<size_t>(chainNum - 1));
      }
    }
    if (capacity > 0) {
      // sample numbers are DRAWS, oldest first, exactly as getTrees reads them
      size_t numDraws = shape.numSavedDraws;
      if (Rf_isNull(sampleNumsExpr)) {
        for (size_t i = 0; i < numDraws; ++i) sampleIndices.push_back(i);
      } else {
        for (R_xlen_t i = 0; i < Rf_xlength(sampleNumsExpr); ++i) {
          int sampleNum = INTEGER(sampleNumsExpr)[i];
          if (sampleNum < 1 || static_cast<size_t>(sampleNum) > numDraws)
            Rf_error("bartcore_printTrees: sample number out of range");
          sampleIndices.push_back(static_cast<size_t>(sampleNum - 1));
        }
      }
    }
    if (Rf_isNull(treeNumsExpr)) {
      for (size_t i = 0; i < shape.numTrees; ++i) treeIndices.push_back(i);
    } else {
      for (R_xlen_t i = 0; i < Rf_xlength(treeNumsExpr); ++i) {
        int treeNum = INTEGER(treeNumsExpr)[i];
        if (treeNum < 1 || static_cast<size_t>(treeNum) > shape.numTrees)
          Rf_error("bartcore_printTrees: tree number out of range");
        treeIndices.push_back(static_cast<size_t>(treeNum - 1));
      }
    }

    // forest 0: shape.numTrees, which the tree range check above reads, is
    // forest 0's count on a multi-forest sampler
    sampler.printTrees(chainIndices.data(), chainIndices.size(),
                       sampleIndices.data(), sampleIndices.size(),
                       treeIndices.data(), treeIndices.size(), 0, false);

    return R_NilValue;
  });
}

// resultExpr, when non-null, is a preallocated numeric filled in place rather
// than a fresh allocation, which is what rbart_vi's per-sweep loop relies on.
SEXP bartcore_getLatents(SEXP ptrExpr, SEXP resultExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (holder.sampler->latents(0) == NULL) return R_NilValue;

  bartcore::SamplerShape shape = holder.sampler->shape();
  size_t numObservations = shape.numObservations;
  size_t numChains = shape.numChains;

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

// ---- the exported augmentation helpers, dbartsDrawLatents and
// dbartsWorkingResponse: the per-sweep augmentation draws, run against R's own
// stream on a caller's fit rather than inside a sampler. The laws themselves
// live in bartcore_bridge, which the flat C API's wrapped forms call too.

/// The inputs both helpers read off R; the flat forms fill the same struct
/// from their own arguments.
static AugmentationInputs augmentationInputs(SEXP fitExpr, SEXP yExpr,
                                             SEXP weightsExpr, SEXP offsetExpr,
                                             SEXP ordinalThresholdsExpr,
                                             SEXP sigmaExpr,
                                             SEXP dispersionExpr, SEXP dfExpr) {
  AugmentationInputs in;
  in.numObservations = static_cast<size_t>(Rf_xlength(fitExpr));
  in.fit = REAL(fitExpr);
  in.y = REAL(yExpr);
  in.weights = Rf_isNull(weightsExpr) ? NULL : REAL(weightsExpr);
  in.offset = Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr);
  in.ordinalThresholds =
    Rf_isNull(ordinalThresholdsExpr) ? NULL : REAL(ordinalThresholdsExpr);
  in.numOrdinalThresholds = Rf_isNull(ordinalThresholdsExpr)
    ? 0 : static_cast<size_t>(Rf_xlength(ordinalThresholdsExpr));
  in.sigma = Rf_asReal(sigmaExpr);       // a null scalar reads as NA, which
  in.dispersion = Rf_asReal(dispersionExpr);  // only an arm that ignores it
  in.df = Rf_asReal(dfExpr);                  // ever sees
  return in;
}

// R/augmentation.R has validated every length, every family's applicable
// arguments and the scalars' ranges; the response support is stated HERE, by
// the same function every conduit that swaps a y calls.
SEXP bartcore_drawLatents(SEXP familyExpr, SEXP fitExpr, SEXP yExpr,
                          SEXP weightsExpr, SEXP offsetExpr, SEXP sigmaExpr,
                          SEXP dispersionExpr, SEXP ordinalThresholdsExpr,
                          SEXP dfExpr) {
  bartcore::ResponseFamily family =
    augmentationFamily(CHAR(STRING_ELT(familyExpr, 0)));
  AugmentationInputs in =
    augmentationInputs(fitExpr, yExpr, weightsExpr, offsetExpr,
                       ordinalThresholdsExpr, sigmaExpr, dispersionExpr,
                       dfExpr);
  validateResponseSupport(family, in.numOrdinalThresholds + 1, in.y,
                          in.numObservations, "dbartsDrawLatents");
  // everything that can longjmp runs BEFORE the generator exists, the result
  // vector included, so the draw loop cannot strand it
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(in.numObservations)));
  using RF = bartcore::ResponseFamily;
  bool precision = family == RF::logistic || family == RF::nbinom ||
                   family == RF::gaussian;
  drawAugmentation(family, in, REAL(result), "dbartsDrawLatents");
  setAttribByName(result, "quantity",
                  Rf_mkString(precision ? "precision" : "location"));
  UNPROTECT(1);
  return result;
}

SEXP bartcore_workingResponse(SEXP familyExpr, SEXP latentExpr, SEXP yExpr,
                              SEXP weightsExpr, SEXP offsetExpr,
                              SEXP dispersionExpr) {
  bartcore::ResponseFamily family =
    augmentationFamily(CHAR(STRING_ELT(familyExpr, 0)));
  AugmentationInputs in =
    augmentationInputs(latentExpr, yExpr, weightsExpr, offsetExpr, R_NilValue,
                       R_NilValue, dispersionExpr, R_NilValue);
  // the ordinal working response is the latent less the offset, so y never
  // enters it and there is no category count here to state its support against
  if (family != bartcore::ResponseFamily::ordinal)
    validateResponseSupport(family, 0, in.y, in.numObservations,
                            "dbartsWorkingResponse");
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(in.numObservations)));
  computeWorkingResponse(family, in, REAL(latentExpr), REAL(result));
  UNPROTECT(1);
  return result;
}

} // extern "C"

namespace bartcore_bridge {

// ---- the augmentation cores, called by both surfaces: the R helpers
// dbartsDrawLatents/dbartsWorkingResponse above and the flat
// dbarts_drawLatents/dbarts_workingResponse (C_interface.cpp). The laws are
// RESTATED here rather than shared with the response models, which run against
// a different generator.

bartcore::ResponseFamily augmentationFamily(const char* name) {
  using RF = bartcore::ResponseFamily;
  if (std::strcmp(name, "probit") == 0) return RF::probit;
  if (std::strcmp(name, "logistic") == 0) return RF::logistic;
  if (std::strcmp(name, "ordinal") == 0) return RF::ordinal;
  if (std::strcmp(name, "aft") == 0) return RF::aft;
  if (std::strcmp(name, "nbinom") == 0) return RF::nbinom;
  if (std::strcmp(name, "student") != 0)
    Rf_error("unrecognized augmentation family \"%s\"", name);
  return RF::gaussian;
}

/// A per-call generator drawing R's stream, held BY an unwindProtect closure:
/// its cleanup runs this deleter on Rf_error's longjmp as well as on the
/// normal return, where a bare local's destructor would not. Never cached -
/// native mode holds no state of its own.
struct NativeRngDeleter {
  void operator()(ext_rng* rng) const { ext_rng_destroy(rng); }
};
typedef std::unique_ptr<ext_rng, NativeRngDeleter> NativeRng;

/// One draw per observation, each the law its engine site draws
/// (ProbitResponse::refreshLatents, OrdinalResponse::drawLatents,
/// AFTResponse::refreshLatents, LogisticResponse::refreshLatents,
/// NBResponse::drawOmega, TResponse::refreshLatents), NaN fallbacks included.
static void drawAugmentationLaws(ext_rng* rng, bartcore::ResponseFamily family,
                                 const AugmentationInputs& in, double* result) {
  using RF = bartcore::ResponseFamily;
  for (size_t i = 0; i < in.numObservations; ++i) {
    double psi = in.fit[i] + (in.offset != NULL ? in.offset[i] : 0.0);
    switch (family) {
    case RF::probit: {
      double sign = 2.0 * in.y[i] - 1.0;
      double z = sign *
        ext_rng_simulateLowerTruncatedNormalScale1(rng, sign * psi, 0.0);
      result[i] = !std::isnan(z) ? z : sign * DBL_EPSILON;
      break;
    }
    case RF::ordinal: {
      // category k lies in (ordinalThresholds[k - 2],
      // ordinalThresholds[k - 1]], the boundary categories one-sided;
      // numOrdinalThresholds is K - 1
      int k = static_cast<int>(std::lround(in.y[i]));
      if (k <= 1) {
        double z = ext_rng_simulateUpperTruncatedNormalScale1(
          rng, psi, in.ordinalThresholds[0]);
        result[i] = !std::isnan(z) ? z : -DBL_EPSILON;
      } else if (static_cast<size_t>(k) > in.numOrdinalThresholds) {
        double z = ext_rng_simulateLowerTruncatedNormalScale1(
          rng, psi, in.ordinalThresholds[in.numOrdinalThresholds - 1]);
        result[i] = !std::isnan(z) ? z : DBL_EPSILON;
      } else {
        double lower = in.ordinalThresholds[k - 2];
        double upper = in.ordinalThresholds[k - 1];
        double z = ext_rng_simulateTruncatedNormalScale1(rng, psi, lower, upper);
        result[i] = !std::isnan(z) ? z : 0.5 * (lower + upper);
      }
      break;
    }
    case RF::aft: { // the imputed log survival time of a row censored at y
      double draw =
        ext_rng_simulateLowerTruncatedNormal(rng, psi, in.sigma, in.y[i]);
      result[i] = !std::isnan(draw) ? draw : in.y[i];
      break;
    }
    case RF::logistic: {
      long reps = in.weights != NULL ? std::lround(in.weights[i]) : 1L;
      double omega = ext_rng_simulatePolyaGamma(rng, psi);
      for (long c = 1; c < reps; ++c)
        omega += ext_rng_simulatePolyaGamma(rng, psi);
      result[i] = omega;
      break;
    }
    case RF::nbinom:
      result[i] =
        bartcore::simulatePolyaGammaShape(rng, in.y[i] + in.dispersion, psi);
      break;
    case RF::gaussian: { // the Student-t scale mixer
      double residual = in.y[i] - psi;
      result[i] = ext_rng_simulateGamma(
        rng, 0.5 * (in.df + 1.0),
        2.0 / (in.df + residual * residual / (in.sigma * in.sigma)));
      break;
    }
    }
  }
}

void drawAugmentation(bartcore::ResponseFamily family,
                      const AugmentationInputs& in, double* result,
                      const char* caller) {
  // the caller has already allocated result, so nothing that can longjmp runs
  // while the generator is alive
  GetRNGstate();
  unwindProtect([&, rng = NativeRng(ext_rng_createDefault(true))]()
                  mutable -> SEXP {
    if (rng == nullptr) {
      PutRNGstate();
      Rf_error("%s: could not create a random number generator", caller);
    }
    drawAugmentationLaws(rng.get(), family, in, result);
    PutRNGstate();
    return R_NilValue;
  });
}

void computeWorkingResponse(bartcore::ResponseFamily family,
                            const AugmentationInputs& in, const double* latent,
                            double* result) {
  using RF = bartcore::ResponseFamily;
  for (size_t i = 0; i < in.numObservations; ++i) {
    // the latent families' working response IS the drawn latent; the switch
    // names every enumerator and carries no default arm, so a family added
    // without one must fail the build rather than take that reading silently
    double value = latent[i];
    switch (family) {
    case RF::logistic:
      value = (in.weights != NULL ? in.weights[i] : 1.0) * (in.y[i] - 0.5) /
        latent[i];
      break;
    case RF::nbinom:
      value = 0.5 * (in.y[i] - in.dispersion) / latent[i];
      break;
    case RF::gaussian: // the Student-t mixer's working response is y itself
      value = in.y[i];
      break;
    case RF::probit:
    case RF::ordinal:
    case RF::aft:
      break;
    }
    result[i] = value - (in.offset != NULL ? in.offset[i] : 0.0);
  }
}

// Provenance stamp for the on-disk layout storeState/setState exchange.
//
// Registry rule for evolving the format: block names are APPEND-ONLY and a
// shipped name's on-disk encoding is FROZEN. A new capability adds a NEW
// optional block name; setState reads blocks by name (rc_getListElement),
// defaults an absent OPTIONAL block, and refuses - naming the block - only
// when a REQUIRED (or config-conditionally-required) block is missing. So an
// ADDITIVE block addition does NOT bump the version (an old reader ignores
// the unknown name; a new reader defaults it), and MUST NOT bump
// minReadableStateFormatVersion. Only a non-additive change to an existing
// block's encoding - one that cannot be expressed as a new name - bumps
// both. RENAMING a block is such a change, and a silent one: the reader
// would find the old name absent and default an OPTIONAL block rather than
// error, leaving the amplitudes at their construction values. Version 2 was
// the glue block's rename from "bcf", and the floor moved with it so a state
// carrying the old name is refused by version before any block is read.
// Version 3 adds the saved-tree store's recorded-draw count beside its write
// cursor. It is REQUIRED and not defaultable: a version-2 state carries no
// count, and every value one could infer is a misread - capacity promotes a
// partly filled store to full and replays slots nothing wrote, 0 discards a
// full one - so the floor moved with it too.
//
// The rule is stated over block names but governs the TOP-LEVEL ATTRIBUTES the
// same way and for the same reason: they are read by name too, and a reader
// ignores one it does not know. "weights.digest" is such an addition - a state
// written before it carries none, and setState then behaves as it did before
// the attribute existed. Making it REQUIRED behind a floor bump would buy no
// compatibility and orphan in-flight states for nothing.
static const int stateFormatVersion = 3;

// The oldest ENCODING this reader still understands: additive block additions
// leave it here; only a non-additive encoding change raises it. Currently 3:
// the 1.0-0 encoding is the FIRST shipped format, so the pre-release
// development increments (the forests-list restructure included) are collapsed
// into version 1, which the glue rename and then the recorded-draw count
// supersede - no released reader or writer ever saw any of the three. Pre-1.0 states are not a compat target; a state
// with no version attribute reads as 0 and is refused at the floor.
static const int minReadableStateFormatVersion = 3;

// The weights a state was stored under, as the 64-bit digest the engine
// computes over their bytes, little-endian into 8 raw bytes. It travels at
// TOP level rather than per chain because weights are chain-invariant. Its one
// consumer is setState, which compares it against the DESTINATION's own live
// digest: equal means the stored latents were shaped by the weights now in
// force and install unchanged, different means they were not.
static const std::size_t weightsDigestBytes = 8;

SEXP encodeWeightsDigest(std::uint64_t digest) {
  SEXP result = Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(
                                         weightsDigestBytes));
  Rbyte* bytes = RAW(result);
  for (std::size_t i = 0; i < weightsDigestBytes; ++i)
    bytes[i] = static_cast<Rbyte>((digest >> (8 * i)) & 0xffULL);
  return result;
}

std::uint64_t decodeWeightsDigest(const Rbyte* bytes) {
  std::uint64_t digest = 0;
  for (std::size_t i = 0; i < weightsDigestBytes; ++i)
    digest |= static_cast<std::uint64_t>(bytes[i]) << (8 * i);
  return digest;
}

SEXP storeState(bartcore::SamplerBase& sampler) {
  bartcore::SamplerStateData state;
  sampler.getState(state);

  size_t numChains = state.chains.size();
  size_t numObservations = sampler.shape().numObservations;

  // per-forest tree channels (a length-1 list off BCF); the k rides here too
  enum {
    FSLOT_TREE_VARS = 0, FSLOT_TREE_VALUES, FSLOT_TREE_SIZES, FSLOT_TREE_FLAGS,
    FSLOT_TREE_PARAMS, FSLOT_TREE_MASKS,
    FSLOT_SAVED_VARS, FSLOT_SAVED_VALUES, FSLOT_SAVED_SIZES, FSLOT_SAVED_FLAGS,
    FSLOT_SAVED_PARAMS, FSLOT_SAVED_MASKS,
    FSLOT_K, FSLOT_LEAF_SCALE, FSLOT_COUNT
  };
  // append-only here too: leaf.scale is the leaf prior's scale factor, k's
  // other half (mu ~ N(0, (scale / k)^2)), added AFTER k and read as OPTIONAL,
  // so a state written before it existed decodes as absent.
  static const char* forestSlotNames[FSLOT_COUNT] = {
    "tree.vars", "tree.values", "tree.sizes", "tree.flags", "tree.params",
    "tree.masks",
    "saved.vars", "saved.values", "saved.sizes", "saved.flags",
    "saved.params", "saved.masks",
    "k", "leaf.scale"
  };

  // append-only slot registry: new blocks go before SLOT_COUNT and do NOT bump
  // the format version (an old state simply lacks the name and decodes as
  // empty). The variance.* block is the heteroscedastic variance forest's
  // flattened trees.
  enum {
    SLOT_FORESTS = 0, SLOT_SIGMA, SLOT_FIT_SCALE, SLOT_LATENTS,
    SLOT_RANEF, SLOT_TAU,
    SLOT_DART_PROBABILITIES, SLOT_DART_ALPHA, SLOT_DART_UPDATES_SKIPPED,
    SLOT_RNG_STATE, SLOT_GLUE, SLOT_RESID_DF, SLOT_THRESHOLDS, SLOT_DISPERSION,
    SLOT_VARIANCE_VARS, SLOT_VARIANCE_VALUES, SLOT_VARIANCE_SIZES,
    SLOT_VARIANCE_FLAGS, SLOT_VARIANCE_MASKS,
    SLOT_VARIANCE_SAVED_VARS, SLOT_VARIANCE_SAVED_VALUES,
    SLOT_VARIANCE_SAVED_SIZES, SLOT_VARIANCE_SAVED_FLAGS,
    SLOT_VARIANCE_SAVED_MASKS,
    SLOT_COUNT
  };
  static const char* slotNames[SLOT_COUNT] = {
    "forests", "sigma", "fit.scale",
    "latents", "ranef", "tau",
    "dart.probabilities", "dart.alpha", "dart.updates.skipped",
    "rng.state", "glue", "resid.df", "thresholds", "dispersion",
    "variance.vars", "variance.values", "variance.sizes", "variance.flags",
    "variance.masks",
    "variance.saved.vars", "variance.saved.values", "variance.saved.sizes",
    "variance.saved.flags", "variance.saved.masks"
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
      SET_VECTOR_ELT(forestExpr, FSLOT_LEAF_SCALE, Rf_ScalarReal(fs.leafScale));
      SET_VECTOR_ELT(forestsExpr, static_cast<R_xlen_t>(f), forestExpr);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(chainExpr, SLOT_FORESTS, forestsExpr);
    UNPROTECT(1);

    // heteroscedastic variance forest: its flattened LIVE trees ride four
    // chain-level slots and its SAVED (keepTrees) buffer four more, each with a
    // pooled-categorical mask channel - a scale leaf carries no leaf params,
    // but nothing confines a variance tree to narrow columns. Empty off a
    // variance forest, so the slots stay R_NilValue.
    if (!chainState.varianceTrees.empty()) {
      storeFlatTrees(chainExpr, SLOT_VARIANCE_VARS, SLOT_VARIANCE_VALUES,
                     SLOT_VARIANCE_SIZES, SLOT_VARIANCE_FLAGS,
                     chainState.varianceTrees);
      if (!chainState.varianceTreeMasks.empty())
        storeTreeMasks(chainExpr, SLOT_VARIANCE_MASKS,
                       chainState.varianceTreeMasks);
    }
    if (!chainState.savedVarianceTrees.empty()) {
      storeFlatTrees(chainExpr, SLOT_VARIANCE_SAVED_VARS,
                     SLOT_VARIANCE_SAVED_VALUES, SLOT_VARIANCE_SAVED_SIZES,
                     SLOT_VARIANCE_SAVED_FLAGS,
                     chainState.savedVarianceTrees);
      if (!chainState.savedVarianceTreeMasks.empty())
        storeTreeMasks(chainExpr, SLOT_VARIANCE_SAVED_MASKS,
                       chainState.savedVarianceTreeMasks);
    }

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

    // the amplitude glue, SELF-DESCRIBING because the vector is ragged:
    // [K, q_1..q_K, sum q_f amplitudes, K prior variances]. A total is not a
    // layout - q = (1, 3) and q = (2, 2) carry four amplitudes each - so the
    // widths travel with the values, and both readers reconstruct the layout
    // from them rather than from a length.
    if (chainState.hasAmplitudes) {
      size_t numForests = chainState.amplitudeWidths.size();
      size_t numAmplitudes = chainState.amplitudes.size();
      SET_VECTOR_ELT(chainExpr, SLOT_GLUE,
                     Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                       1 + 2 * numForests + numAmplitudes)));
      double* glue = REAL(VECTOR_ELT(chainExpr, SLOT_GLUE));
      glue[0] = static_cast<double>(numForests);
      for (size_t f = 0; f < numForests; ++f)
        glue[1 + f] = static_cast<double>(chainState.amplitudeWidths[f]);
      for (size_t j = 0; j < numAmplitudes; ++j)
        glue[1 + numForests + j] = chainState.amplitudes[j];
      for (size_t f = 0; f < numForests; ++f)
        glue[1 + numForests + numAmplitudes + f] =
          chainState.amplitudeVariances[f];
    }

    // the t-only residual df; lambda already rode the latents slot above. A
    // gaussian (or any non-t) chain carries NaN and writes no block.
    if (std::isfinite(chainState.residualDf))
      SET_VECTOR_ELT(chainExpr, SLOT_RESID_DF,
                     Rf_ScalarReal(chainState.residualDf));

    // the ordinal-only length-(K-1) threshold vector; z already rode the latents
    // slot above. A non-ordinal chain leaves it empty and writes no block, so
    // old and other-family states omit the slot.
    if (!chainState.ordinalThresholds.empty()) {
      SET_VECTOR_ELT(chainExpr, SLOT_THRESHOLDS,
                     Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                      chainState.ordinalThresholds.size())));
      std::memcpy(REAL(VECTOR_ELT(chainExpr, SLOT_THRESHOLDS)),
                  chainState.ordinalThresholds.data(),
                  chainState.ordinalThresholds.size() * sizeof(double));
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
  setAttribByName(resultExpr, "cutPoints", cutPointsExpr);
  setAttribByName(resultExpr, "currentSampleNum",
                  Rf_ScalarInteger(static_cast<int>(state.currentSampleNum)));
  // the cursor alone does not say how much of the ring is real; both travel
  setAttribByName(resultExpr, "recordedDraws",
                  Rf_ScalarInteger(static_cast<int>(state.recordedDraws)));
  setAttribByName(resultExpr, "formatVersion",
                  Rf_ScalarInteger(stateFormatVersion));
  // the weights themselves do not ride the state; their digest does, so a
  // restore can tell whether the latents it carries belong to the weights the
  // destination holds
  setAttribByName(resultExpr, "weights.digest",
                  encodeWeightsDigest(sampler.weightsDigest()));
  setAttribByName(resultExpr, "packageVersion", Rf_mkString(PACKAGE_VERSION));
  SEXP classExpr = PROTECT(Rf_mkString("bartcoreState"));
  Rf_setAttrib(resultExpr, R_ClassSymbol, classExpr);

  UNPROTECT(5);
  return resultExpr;
}

// Reads the self-describing amplitude glue block - [K, q_1..q_K, sum q_f
// amplitudes, K prior variances] - into a chain state; false on anything that
// does not parse as that layout, which each caller reports in its own words.
// The LAYOUT is validated here and the layout's agreement with the live
// sampler at Chain::stateIsValid, because a total is not a layout: q = (1, 3)
// and q = (2, 2) are both four amplitudes, and restoreGlue writes through the
// live offsets.
bool readAmplitudeGlue(SEXP glueExpr, bartcore::ChainStateData& chainState) {
  if (!Rf_isReal(glueExpr) || Rf_xlength(glueExpr) < 1) return false;
  const double* glue = REAL(glueExpr);
  double numForestsValue = glue[0];
  if (!(numForestsValue >= 1.0) ||
      numForestsValue != std::floor(numForestsValue))
    return false;
  size_t numForests = static_cast<size_t>(numForestsValue);
  if (static_cast<size_t>(Rf_xlength(glueExpr)) < 1 + 2 * numForests)
    return false;
  size_t numAmplitudes = 0;
  chainState.amplitudeWidths.resize(numForests);
  for (size_t f = 0; f < numForests; ++f) {
    double width = glue[1 + f];
    if (!(width >= 1.0) || width != std::floor(width)) return false;
    chainState.amplitudeWidths[f] = static_cast<size_t>(width);
    numAmplitudes += chainState.amplitudeWidths[f];
  }
  if (static_cast<size_t>(Rf_xlength(glueExpr)) !=
      1 + 2 * numForests + numAmplitudes)
    return false;
  chainState.hasAmplitudes = true;
  chainState.amplitudes.assign(glue + 1 + numForests,
                               glue + 1 + numForests + numAmplitudes);
  chainState.amplitudeVariances.assign(
    glue + 1 + numForests + numAmplitudes,
    glue + 1 + 2 * numForests + numAmplitudes);
  return true;
}

/// The single refusal both tree-install entries report when an incoming tree
/// splits on a column the recipient forest's mask forbids. setState and
/// installForests run the one predicate, so they speak with the one voice: a
/// state either entry refuses is refused by the other, in the same words.
static const char* const columnMaskMismatchMessage =
  "warm-start donor holds a tree that splits on a variable outside "
  "this forest's allowed column set; the donor's fit is "
  "incompatible with the column restriction (a forest's own "
  "column subset or a restricted variance forest) in force here";

void setState(bartcore::SamplerBase& sampler, SEXP stateExpr,
              const double* currentPredictors) {
  bartcore::SamplerShape shape = sampler.shape();
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

  if (static_cast<size_t>(Rf_xlength(stateExpr)) != shape.numChains)
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
  state.chains.resize(shape.numChains);

  // Whether the stored latents were shaped by other weights than the ones in
  // force here. ABSENT (a state written before the attribute existed) is not a
  // mismatch: it reconciles nothing and behaves exactly as this reader did
  // before the attribute.
  bool weightsDiffer = false;
  SEXP weightsDigestExpr =
    Rf_getAttrib(stateExpr, Rf_install("weights.digest"));
  if (!Rf_isNull(weightsDigestExpr)) {
    if (TYPEOF(weightsDigestExpr) != RAWSXP ||
        static_cast<size_t>(Rf_xlength(weightsDigestExpr)) !=
          weightsDigestBytes)
      errorMessage = "malformed weights digest in bartcore state";
    else
      weightsDiffer = decodeWeightsDigest(RAW(weightsDigestExpr)) !=
        sampler.weightsDigest();
  }

  SEXP cutPointsExpr = Rf_getAttrib(stateExpr, Rf_install("cutPoints"));
  if (errorMessage == NULL &&
      (Rf_isNull(cutPointsExpr) ||
       static_cast<size_t>(Rf_xlength(cutPointsExpr)) !=
         shape.numPredictors)) {
    errorMessage = "malformed cut points in bartcore state";
  } else if (errorMessage == NULL) {
    state.cutPoints.resize(shape.numPredictors);
    for (size_t j = 0; j < shape.numPredictors; ++j) {
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

  SEXP recordedDrawsExpr =
    PROTECT(Rf_getAttrib(stateExpr, Rf_install("recordedDraws")));
  if (errorMessage == NULL) {
    if (!Rf_isInteger(recordedDrawsExpr) ||
        Rf_xlength(recordedDrawsExpr) != 1 ||
        INTEGER(recordedDrawsExpr)[0] < 0)
      errorMessage = "malformed recorded draw count in bartcore state";
    else
      state.recordedDraws =
        static_cast<size_t>(INTEGER(recordedDrawsExpr)[0]);
  }
  UNPROTECT(1);

  for (size_t c = 0; c < shape.numChains && errorMessage == NULL; ++c) {
    SEXP chainExpr = VECTOR_ELT(stateExpr, static_cast<R_xlen_t>(c));
    bartcore::ChainStateData& chainState(state.chains[c]);

    SEXP forestsExpr = rc_getListElement(chainExpr, "forests");
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
    size_t numLeafCovariates = shape.numLeafCovariates;
    for (size_t f = 0; f < numForests && errorMessage == NULL; ++f) {
      SEXP forestExpr = VECTOR_ELT(forestsExpr, static_cast<R_xlen_t>(f));
      bartcore::ForestStateData& fs(chainState.forests[f]);

      if (Rf_isNull(rc_getListElement(forestExpr, "tree.vars"))) {
        errorMessage = missingBlock("tree.vars");
        break;
      }
      if (!readFlatTrees(rc_getListElement(forestExpr, "tree.vars"),
                         rc_getListElement(forestExpr, "tree.values"),
                         rc_getListElement(forestExpr, "tree.sizes"),
                         rc_getListElement(forestExpr, "tree.flags"),
                         sampler.data(), fs.trees, &errorMessage))
        break;
      SEXP savedSizesExpr = rc_getListElement(forestExpr, "saved.sizes");
      if (!Rf_isNull(savedSizesExpr) &&
          !readFlatTrees(rc_getListElement(forestExpr, "saved.vars"),
                         rc_getListElement(forestExpr, "saved.values"),
                         savedSizesExpr,
                         rc_getListElement(forestExpr, "saved.flags"),
                         sampler.data(), fs.savedTrees, &errorMessage))
        break;

      // linear-leaf states must carry their slope arrays; function-valued
      // states carry fits slabs and variable-length saved blocks instead
      if (shape.usesFunctionLeaves) {
        if (Rf_isNull(rc_getListElement(forestExpr, "tree.params"))) {
          errorMessage = missingBlock("tree.params");
          break;
        }
        if (!readFunctionTreeParams(
              rc_getListElement(forestExpr, "tree.params"), fs.trees.size(),
              shape.numObservations, fs.treeParams, &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readFunctionSavedParams(
              rc_getListElement(forestExpr, "saved.params"), fs.savedTrees,
              numLeafCovariates, fs.savedTreeParams, &errorMessage))
          break;
      } else if (numLeafCovariates > 0) {
        if (Rf_isNull(rc_getListElement(forestExpr, "tree.params"))) {
          errorMessage = missingBlock("tree.params");
          break;
        }
        if (!readTreeParams(rc_getListElement(forestExpr, "tree.params"),
                            fs.trees, numLeafCovariates, fs.treeParams,
                            &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeParams(rc_getListElement(forestExpr, "saved.params"),
                            fs.savedTrees, numLeafCovariates,
                            fs.savedTreeParams, &errorMessage))
          break;
      }

      // pooled-categorical states must carry their mask channels
      if (sampler.data().hasPooledCategorical) {
        if (Rf_isNull(rc_getListElement(forestExpr, "tree.masks"))) {
          errorMessage = missingBlock("tree.masks");
          break;
        }
        if (!readTreeMasks(
              rc_getListElement(forestExpr, "tree.masks"), fs.trees,
              sampler.data(), fs.treeMasks, &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeMasks(rc_getListElement(forestExpr, "saved.masks"),
                           fs.savedTrees, sampler.data(), fs.savedTreeMasks,
                           &errorMessage))
          break;
      }

      SEXP kExpr = rc_getListElement(forestExpr, "k");
      if (Rf_isNull(kExpr)) {
        errorMessage = missingBlock("k");
        break;
      }
      if (!Rf_isReal(kExpr) || Rf_xlength(kExpr) != 1) {
        errorMessage = malformedBlock("k");
        break;
      }
      fs.k = REAL(kExpr)[0];

      // OPTIONAL, append-only: a state written before the block existed lacks
      // the name and keeps the 0.0 ABSENT sentinel, restoring exactly as it
      // did then. Presence need not be uniform across a chain's forests. The
      // type/length check names the block; a present but non-finite or
      // non-positive VALUE is NOT refused - it falls through the restore
      // paths' > 0.0 guard as absent, matching k's permissive posture rather
      // than inventing a stricter one.
      SEXP leafScaleExpr = rc_getListElement(forestExpr, "leaf.scale");
      if (!Rf_isNull(leafScaleExpr)) {
        if (!Rf_isReal(leafScaleExpr) || Rf_xlength(leafScaleExpr) != 1) {
          errorMessage = malformedBlock("leaf.scale");
          break;
        }
        fs.leafScale = REAL(leafScaleExpr)[0];
      }
    }
    if (errorMessage != NULL) break;

    SEXP sigmaExpr = rc_getListElement(chainExpr, "sigma");
    if (Rf_isNull(sigmaExpr)) {
      errorMessage = missingBlock("sigma");
      break;
    }
    if (!Rf_isReal(sigmaExpr) || Rf_xlength(sigmaExpr) != 1) {
      errorMessage = malformedBlock("sigma");
      break;
    }
    chainState.sigma = REAL(sigmaExpr)[0];

    // heteroscedastic variance forest: optional blocks, absent (empty) off a
    // variance state. stateIsValid refuses a variance sampler lacking them and
    // a homoscedastic sampler carrying them (the additive-block contract), and
    // pairs the saved buffer's size against the live capacity. readFlatTrees
    // and readTreeMasks carry no name channel, so their generic refusals are
    // re-stated here against the block that failed.
    SEXP varianceVarsExpr = rc_getListElement(chainExpr, "variance.vars");
    if (!Rf_isNull(varianceVarsExpr)) {
      if (!readFlatTrees(varianceVarsExpr,
                         rc_getListElement(chainExpr, "variance.values"),
                         rc_getListElement(chainExpr, "variance.sizes"),
                         rc_getListElement(chainExpr, "variance.flags"),
                         sampler.data(), chainState.varianceTrees,
                         &errorMessage)) {
        errorMessage = malformedBlock("variance.vars");
        break;
      }
      if (sampler.data().hasPooledCategorical &&
          !readTreeMasks(rc_getListElement(chainExpr, "variance.masks"),
                         chainState.varianceTrees, sampler.data(),
                         chainState.varianceTreeMasks, &errorMessage)) {
        errorMessage = malformedBlock("variance.masks");
        break;
      }
    }
    SEXP varianceSavedVarsExpr =
      rc_getListElement(chainExpr, "variance.saved.vars");
    if (!Rf_isNull(varianceSavedVarsExpr)) {
      if (!readFlatTrees(varianceSavedVarsExpr,
                         rc_getListElement(chainExpr, "variance.saved.values"),
                         rc_getListElement(chainExpr, "variance.saved.sizes"),
                         rc_getListElement(chainExpr, "variance.saved.flags"),
                         sampler.data(), chainState.savedVarianceTrees,
                         &errorMessage)) {
        errorMessage = malformedBlock("variance.saved.vars");
        break;
      }
      if (sampler.data().hasPooledCategorical &&
          !readTreeMasks(rc_getListElement(chainExpr, "variance.saved.masks"),
                         chainState.savedVarianceTrees, sampler.data(),
                         chainState.savedVarianceTreeMasks, &errorMessage)) {
        errorMessage = malformedBlock("variance.saved.masks");
        break;
      }
    }

    SEXP fitScaleExpr = rc_getListElement(chainExpr, "fit.scale");
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

    SEXP latentsExpr = rc_getListElement(chainExpr, "latents");
    if (!Rf_isNull(latentsExpr)) {
      if (!Rf_isReal(latentsExpr)) {
        errorMessage = "malformed latents in bartcore state";
        break;
      }
      chainState.latents.assign(
        REAL(latentsExpr), REAL(latentsExpr) + Rf_xlength(latentsExpr));
    }

    SEXP ranefExpr = rc_getListElement(chainExpr, "ranef");
    if (!Rf_isNull(ranefExpr)) {
      SEXP tauExpr = rc_getListElement(chainExpr, "tau");
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
      rc_getListElement(chainExpr, "dart.probabilities");
    if (!Rf_isNull(dartProbabilitiesExpr)) {
      SEXP dartAlphaExpr = rc_getListElement(chainExpr, "dart.alpha");
      SEXP dartSkippedExpr =
        rc_getListElement(chainExpr, "dart.updates.skipped");
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

    SEXP rngStateExpr = rc_getListElement(chainExpr, "rng.state");
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

    SEXP glueExpr = rc_getListElement(chainExpr, "glue");
    if (!Rf_isNull(glueExpr)) {
      if (!readAmplitudeGlue(glueExpr, chainState)) {
        errorMessage = "malformed amplitude glue in bartcore state";
        break;
      }
    }

    // additive t-only block: absent (an old or gaussian state) leaves the NaN
    // default, which stateIsValid refuses only for a t sampler
    SEXP residDfExpr = rc_getListElement(chainExpr, "resid.df");
    if (!Rf_isNull(residDfExpr)) {
      if (!Rf_isReal(residDfExpr) || Rf_xlength(residDfExpr) != 1) {
        errorMessage = malformedBlock("resid.df");
        break;
      }
      chainState.residualDf = REAL(residDfExpr)[0];
    }

    // additive ordinal-only block: absent (an old or non-ordinal state) leaves
    // the vector empty, which stateIsValid refuses only for an ordinal sampler
    SEXP ordinalThresholdsExpr = rc_getListElement(chainExpr, "thresholds");
    if (!Rf_isNull(ordinalThresholdsExpr)) {
      if (!Rf_isReal(ordinalThresholdsExpr)) {
        errorMessage = malformedBlock("thresholds");
        break;
      }
      chainState.ordinalThresholds.assign(
        REAL(ordinalThresholdsExpr),
        REAL(ordinalThresholdsExpr) + Rf_xlength(ordinalThresholdsExpr));
    }

    // additive nbinom-only block: absent (an old or non-count state) leaves the
    // NaN default, which stateIsValid refuses only for an NB sampler
    SEXP dispersionExpr = rc_getListElement(chainExpr, "dispersion");
    if (!Rf_isNull(dispersionExpr)) {
      if (!Rf_isReal(dispersionExpr) || Rf_xlength(dispersionExpr) != 1) {
        errorMessage = malformedBlock("dispersion");
        break;
      }
      chainState.dispersion = REAL(dispersionExpr)[0];
    }
  }

  bool columnMaskRefused = false;
  bool restored =
    errorMessage == NULL &&
    sampler.setState(state, currentPredictors, &columnMaskRefused);
  {
    bartcore::SamplerStateData empty;
    std::swap(state, empty);  // free before a potential longjmp
  }
  if (errorMessage != NULL) Rf_error("%s", errorMessage);
  if (columnMaskRefused) Rf_error("%s", columnMaskMismatchMessage);
  if (!restored)
    Rf_error("state is not consistent with this sampler");

  // The latents just installed were drawn against the SOURCE's weights; these
  // are not them. Re-derive against the weights in force, which is where the
  // live conduit would have put them: a state install lands where setWeights
  // lands rather than pairing one vector's latents with another's. Silent and
  // deterministic - it consumes each chain's own restored generator - and
  // self-selecting, since for a family that states nothing against its
  // weights it is a measured no-op.
  if (weightsDiffer) sampler.reapplyWeights();
}

// Parses a "bartcoreState" donor into a SamplerStateData for a warm start,
// validating flat trees against the destination sampler's data. Only the
// channels a warm start consumes are read (trees, leaf params, masks, k,
// sigma, the fit scale, DART, and the amplitude glue); latents, group effects,
// rng are left for the destination to redraw. Function-leaf donors seed from
// their live trees, so their saved channel is skipped. The donor's own chain
// count is honored (a short donor may seed many chains). Returns an error
// string, or NULL, rather than longjmping so the caller can free state first.
static const char* readWarmStartState(SEXP stateExpr,
                                      bartcore::SamplerBase& sampler,
                                      bartcore::SamplerStateData& state) {
  bartcore::SamplerShape shape = sampler.shape();
  size_t numChains = static_cast<size_t>(Rf_xlength(stateExpr));
  if (numChains == 0) return "warm-start donor holds no chains";

  SEXP cutPointsExpr = Rf_getAttrib(stateExpr, Rf_install("cutPoints"));
  if (Rf_isNull(cutPointsExpr) ||
      static_cast<size_t>(Rf_xlength(cutPointsExpr)) != shape.numPredictors)
    return "malformed cut points in warm-start donor";
  state.cutPoints.resize(shape.numPredictors);
  for (size_t j = 0; j < shape.numPredictors; ++j) {
    SEXP cutsExpr = VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(j));
    if (Rf_isNull(cutsExpr)) continue;
    if (!Rf_isReal(cutsExpr)) return "malformed cut points in warm-start donor";
    state.cutPoints[j].assign(REAL(cutsExpr),
                              REAL(cutsExpr) + Rf_xlength(cutsExpr));
  }

  // the donor pool is addressed by DRAW, so the donor's own write position and
  // draw count are read here rather than left to setState
  SEXP sampleNumExpr = Rf_getAttrib(stateExpr, Rf_install("currentSampleNum"));
  SEXP recordedDrawsExpr = Rf_getAttrib(stateExpr, Rf_install("recordedDraws"));
  if (!Rf_isInteger(sampleNumExpr) || Rf_xlength(sampleNumExpr) != 1 ||
      INTEGER(sampleNumExpr)[0] < 0 || !Rf_isInteger(recordedDrawsExpr) ||
      Rf_xlength(recordedDrawsExpr) != 1 || INTEGER(recordedDrawsExpr)[0] < 0)
    return "malformed sample position in warm-start donor";
  state.currentSampleNum = static_cast<size_t>(INTEGER(sampleNumExpr)[0]);
  state.recordedDraws = static_cast<size_t>(INTEGER(recordedDrawsExpr)[0]);

  const char* errorMessage = NULL;
  state.chains.resize(numChains);
  size_t numLeafCovariates = shape.numLeafCovariates;
  bool functionLeaves = shape.usesFunctionLeaves;
  for (size_t c = 0; c < numChains && errorMessage == NULL; ++c) {
    SEXP chainExpr = VECTOR_ELT(stateExpr, static_cast<R_xlen_t>(c));
    bartcore::ChainStateData& chainState(state.chains[c]);

    SEXP forestsExpr = rc_getListElement(chainExpr, "forests");
    if (Rf_isNull(forestsExpr) || TYPEOF(forestsExpr) != VECSXP) {
      errorMessage = "malformed forests in warm-start donor";
      break;
    }
    size_t numForests = static_cast<size_t>(Rf_xlength(forestsExpr));
    if (numForests == 0) {
      // the donor pool reads each chain's first forest, and a chain that
      // declares none never reaches the engine's own shape check
      errorMessage = "warm-start donor chain holds no forests";
      break;
    }
    chainState.forests.resize(numForests);
    for (size_t f = 0; f < numForests && errorMessage == NULL; ++f) {
      SEXP forestExpr = VECTOR_ELT(forestsExpr, static_cast<R_xlen_t>(f));
      bartcore::ForestStateData& fs(chainState.forests[f]);

      if (!readFlatTrees(rc_getListElement(forestExpr, "tree.vars"),
                         rc_getListElement(forestExpr, "tree.values"),
                         rc_getListElement(forestExpr, "tree.sizes"),
                         rc_getListElement(forestExpr, "tree.flags"),
                         sampler.data(), fs.trees, &errorMessage))
        break;
      SEXP savedSizesExpr = rc_getListElement(forestExpr, "saved.sizes");
      if (!functionLeaves && !Rf_isNull(savedSizesExpr) &&
          !readFlatTrees(rc_getListElement(forestExpr, "saved.vars"),
                         rc_getListElement(forestExpr, "saved.values"),
                         savedSizesExpr, rc_getListElement(forestExpr,
                         "saved.flags"), sampler.data(), fs.savedTrees,
                         &errorMessage))
        break;

      if (functionLeaves) {
        if (!readFunctionTreeParams(
              rc_getListElement(forestExpr, "tree.params"), fs.trees.size(),
              shape.numObservations, fs.treeParams, &errorMessage))
          break;
      } else if (numLeafCovariates > 0) {
        if (!readTreeParams(rc_getListElement(forestExpr, "tree.params"),
                            fs.trees, numLeafCovariates, fs.treeParams,
                            &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeParams(rc_getListElement(forestExpr, "saved.params"),
                            fs.savedTrees, numLeafCovariates,
                            fs.savedTreeParams, &errorMessage))
          break;
      }

      if (sampler.data().hasPooledCategorical) {
        if (!readTreeMasks(
              rc_getListElement(forestExpr, "tree.masks"), fs.trees,
              sampler.data(), fs.treeMasks, &errorMessage))
          break;
        if (!fs.savedTrees.empty() &&
            !readTreeMasks(rc_getListElement(forestExpr, "saved.masks"),
                           fs.savedTrees, sampler.data(), fs.savedTreeMasks,
                           &errorMessage))
          break;
      }

      SEXP kExpr = rc_getListElement(forestExpr, "k");
      if (!Rf_isReal(kExpr) || Rf_xlength(kExpr) != 1) {
        errorMessage = "malformed parameters in warm-start donor";
        break;
      }
      fs.k = REAL(kExpr)[0];

      // optional as in the setState parser above; installForest adopts it
      // alongside k, so a donor's leaf calibration seeds the warm start
      SEXP leafScaleExpr = rc_getListElement(forestExpr, "leaf.scale");
      if (!Rf_isNull(leafScaleExpr)) {
        if (!Rf_isReal(leafScaleExpr) || Rf_xlength(leafScaleExpr) != 1) {
          errorMessage = "malformed parameters in warm-start donor";
          break;
        }
        fs.leafScale = REAL(leafScaleExpr)[0];
      }
    }
    if (errorMessage != NULL) break;

    SEXP sigmaExpr = rc_getListElement(chainExpr, "sigma");
    SEXP fitScaleExpr = rc_getListElement(chainExpr, "fit.scale");
    if (!Rf_isReal(sigmaExpr) || Rf_xlength(sigmaExpr) != 1 ||
        !Rf_isReal(fitScaleExpr) || Rf_xlength(fitScaleExpr) != 2) {
      errorMessage = "malformed parameters in warm-start donor";
      break;
    }
    chainState.sigma = REAL(sigmaExpr)[0];
    chainState.fitMin = REAL(fitScaleExpr)[0];
    chainState.fitMax = REAL(fitScaleExpr)[1];

    SEXP dartProbabilitiesExpr =
      rc_getListElement(chainExpr, "dart.probabilities");
    if (!Rf_isNull(dartProbabilitiesExpr)) {
      SEXP dartAlphaExpr = rc_getListElement(chainExpr, "dart.alpha");
      SEXP dartSkippedExpr =
        rc_getListElement(chainExpr, "dart.updates.skipped");
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

    SEXP glueExpr = rc_getListElement(chainExpr, "glue");
    if (!Rf_isNull(glueExpr)) {
      if (!readAmplitudeGlue(glueExpr, chainState)) {
        errorMessage = "malformed amplitude glue in warm-start donor";
        break;
      }
    }

    // the LIVE and SAVED variance trees with their mask channels, as in the
    // setState parser: installForests installs the live pair for a
    // live-sourced start and slices the saved buffer for a slot-sourced one.
    // Both blocks stay OPTIONAL here - a donor lacking the saved one is
    // adjudicated by installForests, which knows the slot the caller asked
    // for; the parser cannot.
    SEXP varianceVarsExpr = rc_getListElement(chainExpr, "variance.vars");
    if (!Rf_isNull(varianceVarsExpr)) {
      if (!readFlatTrees(varianceVarsExpr,
                         rc_getListElement(chainExpr, "variance.values"),
                         rc_getListElement(chainExpr, "variance.sizes"),
                         rc_getListElement(chainExpr, "variance.flags"),
                         sampler.data(), chainState.varianceTrees,
                         &errorMessage))
        break;
      if (sampler.data().hasPooledCategorical &&
          !readTreeMasks(rc_getListElement(chainExpr, "variance.masks"),
                         chainState.varianceTrees, sampler.data(),
                         chainState.varianceTreeMasks, &errorMessage))
        break;
    }
    SEXP varianceSavedVarsExpr =
      rc_getListElement(chainExpr, "variance.saved.vars");
    if (!Rf_isNull(varianceSavedVarsExpr)) {
      if (!readFlatTrees(varianceSavedVarsExpr,
                         rc_getListElement(chainExpr, "variance.saved.values"),
                         rc_getListElement(chainExpr, "variance.saved.sizes"),
                         rc_getListElement(chainExpr, "variance.saved.flags"),
                         sampler.data(), chainState.savedVarianceTrees,
                         &errorMessage))
        break;
      if (sampler.data().hasPooledCategorical &&
          !readTreeMasks(rc_getListElement(chainExpr, "variance.saved.masks"),
                         chainState.savedVarianceTrees, sampler.data(),
                         chainState.savedVarianceTreeMasks, &errorMessage))
        break;
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
    // Donor pool of (chain, slot): the donor's RECORDED draws, oldest first -
    // the order predict and getTrees report, so 'samples' names the same draw
    // they do - or each donor chain's live forest (slot -1) when it kept no
    // trees, or kept a store nothing was recorded into.
    std::vector<std::pair<size_t, int>> pool;
    if (errorMessage == NULL) {
      for (size_t dc = 0; dc < donor.chains.size(); ++dc) {
        const bartcore::ForestStateData& f0 = donor.chains[dc].forests[0];
        size_t capacity = !f0.savedTrees.empty() && !f0.trees.empty()
          ? f0.savedTrees.size() / f0.trees.size() : 0;
        size_t filled = donor.recordedDraws < capacity ? donor.recordedDraws
                                                       : capacity;
        if (filled > 0) {
          size_t cursor = donor.currentSampleNum % capacity;
          for (size_t i = 0; i < filled; ++i)
            pool.emplace_back(dc, static_cast<int>(
              (cursor + capacity - filled + i) % capacity));
        } else {
          pool.emplace_back(dc, -1);
        }
      }
      if (pool.empty()) errorMessage = "warm-start donor holds no samples";
    }

    size_t numChains = sampler.shape().numChains;
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
      Rf_error("this sampler does not support the warm-start donor's "
               "predictor structure (a categorical/continuous column "
               "mismatch, or a malformed donor cut grid); its trees cannot "
               "be remapped onto this cut grid");
    case bartcore::WarmStartResult::dartMismatch:
      Rf_error("DART state transfers only between two DART fits; the donor and "
               "destination disagree on dart");
    case bartcore::WarmStartResult::interactionMismatch:
      Rf_error("warm-start donor holds a tree that violates this sampler's "
               "interaction constraint; the donor's fit is incompatible with "
               "the interactions() prior in force here");
    case bartcore::WarmStartResult::columnMaskMismatch:
      Rf_error("%s", columnMaskMismatchMessage);
    case bartcore::WarmStartResult::varianceMismatch:
      Rf_error("warm-start donor's variance trees cannot be installed on this "
               "sampler's data (a rebuilt variance tree leaves a leaf empty, a "
               "scale leaf is not positive, or a flat tree failed to rebuild); "
               "the donor's variance surface is incompatible with this "
               "sampler");
    case bartcore::WarmStartResult::varianceSlotMismatch:
      Rf_error("warm-start donor's saved variance buffer does not hold the "
               "requested sample; a warm start from a saved sample installs "
               "that sample's own scale surface, and this donor's is missing "
               "or does not match its variance forest");
    case bartcore::WarmStartResult::shapeMismatch:
      Rf_error("warm-start donor is not shape-compatible with this sampler "
               "(number of trees, forests, or predictors differ, or only one "
               "of the two carries a variance forest)");
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
                 const bartcore::PredictorSource* replaySource,
                 size_t replayNumRows, size_t forestIndex,
                 GatheredTrees& out) {
  const bartcore::ColumnStore& store(sampler.data());

  std::vector<bartcore::FlatNode> liveNodes;
  std::vector<double> liveSlopes;
  std::vector<std::uint64_t> liveMasks;
  std::vector<std::uint32_t> counts;
  std::vector<size_t> replayIndices(replayNumRows);
  // one reader set for every tree replayed: a dense replay source resolves to
  // its columns' pointers, a sparse one lays its rank bitmaps down once
  bartcore::PredictorSourceColumns replayColumns(
    replaySource != NULL ? *replaySource : bartcore::PredictorSource{},
    store.types.data());
  std::string directionsScratch;
  bool functionLeaves = sampler.shape().usesFunctionLeaves;
  // the mask side channel exists only for pooled columns (past 63 levels)
  bool anyPooled = false;
  for (size_t j = 0; j < store.numPredictors; ++j)
    if (store.columnIsPooled(j)) { anyPooled = true; break; }

  for (size_t i = 0; i < numChainIndices; ++i) {
    size_t chainNum = chainIndices[i];
    for (size_t j = 0; j < numSampleIndices; ++j) {
      size_t sampleNum = useSaved ? sampleIndices[j] : 0;
      // the reported sample number is the DRAW, the read is of the slot
      // holding it
      size_t slot = useSaved ? sampler.savedSlotForDraw(sampleNum) : 0;
      for (size_t k = 0; k < numTreeIndices; ++k) {
        size_t treeNum = treeIndices[k];

        const std::vector<bartcore::FlatNode>* nodes;
        const std::vector<double>* slopes = NULL;
        const std::vector<std::uint64_t>* masks = NULL;
        if (useSaved) {
          nodes = &sampler.savedTree(chainNum, slot, treeNum, forestIndex);
          if (out.numSlopes > 0)
            slopes = &sampler.savedTreeSlopes(chainNum, slot, treeNum,
                                              forestIndex);
          if (anyPooled)
            masks = &sampler.savedTreeMasks(chainNum, slot, treeNum,
                                            forestIndex);
        } else {
          sampler.flattenTree(chainNum, treeNum, liveNodes, counts,
                              out.numSlopes > 0 ? &liveSlopes : NULL,
                              anyPooled ? &liveMasks : NULL, forestIndex);
          nodes = &liveNodes;
          if (out.numSlopes > 0) slopes = &liveSlopes;
          if (anyPooled) masks = &liveMasks;
        }
        if (replaySource != NULL) {
          counts.resize(nodes->size());
          for (size_t l = 0; l < replayNumRows; ++l) replayIndices[l] = l;
          bartcore::countFlatObservationsBelow(
            nodes->data(), replayColumns, replayIndices.data(), 0,
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

// Builds the R-visible getTrees data.frame from a gather: ([chain,] [sample,]
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

// A data.frame of tree structure in the R-visible getTrees format: pre-order
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
              size_t numTreeIndices, bool useLiveTrees,
              const bartcore::PredictorSource* newdata, size_t newdataNumRows,
              const double* trainingReplay, size_t trainingReplayNumRows,
              size_t forestIndex, const char* caller) {
  const bartcore::ColumnStore& store(sampler.data());
  bartcore::SamplerShape shape = sampler.shape();

  bool useSaved = shape.savedTreeCapacity > 0 && !useLiveTrees;
  if (!useSaved) numSampleIndices = 1;
  if (useSaved) refuseEmptyTreeStore(sampler, caller);

  for (size_t i = 0; i < numChainIndices; ++i) {
    if (chainIndices[i] >= shape.numChains)
      Rf_error("%s: chain number out of range", caller);
  }
  if (useSaved) {
    // sample numbers address RECORDED DRAWS, oldest first, not store slots
    for (size_t i = 0; i < numSampleIndices; ++i) {
      if (sampleIndices[i] >= shape.numSavedDraws)
        Rf_error("%s: sample number out of range", caller);
    }
  }
  for (size_t i = 0; i < numTreeIndices; ++i) {
    if (treeIndices[i] >= sampler.numTreesInForest(forestIndex))
      Rf_error("%s: tree number out of range", caller);
  }

  // saved trees carry no counts of their own and replay the training rows the
  // caller supplies (the engine keeps no matrix); newdata replays its rows
  // through live and saved trees alike
  const bartcore::PredictorSource* replaySource = NULL;
  size_t replayNumRows = 0;
  bartcore::PredictorSource trainingView;
  if (newdata != NULL) {
    replaySource = newdata;
    replayNumRows = newdataNumRows;
  } else if (useSaved && trainingReplay != NULL) {
    trainingView = bartcore::densePredictorSource(
      trainingReplay, trainingReplayNumRows, store.numPredictors);
    replaySource = &trainingView;
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
  size_t numSlopes = shape.usesFunctionLeaves ? 0 : shape.numLeafCovariates;

  GatheredTrees gathered{shape.numChains > 1, useSaved,
                         anyCategorical, anyMissing, numSlopes,
                         {}, {}, {}, {}, {}, {}, {}, {}, {}};
  gathered.slopes.resize(numSlopes);
  gatherTrees(sampler, chainIndices, numChainIndices, sampleIndices,
              numSampleIndices, treeIndices, numTreeIndices, useSaved,
              replaySource, replayNumRows, forestIndex, gathered);

  return emitTreeDataFrame(gathered);
}

} // namespace bartcore_bridge
