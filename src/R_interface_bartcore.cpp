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
enum {
  PROT_DATA = 0,
  PROT_RESPONSE,
  PROT_OFFSET,
  PROT_PREDICTORS,
  PROT_TEST_PREDICTORS,
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
// dbartsModel) the bartcore engine consumes, parsed without the classic
// engine's types. Pointers borrow from the expressions; error paths may
// leak the parse vectors (Rf_error longjmps past destructors), matching
// the classic parse layer's behavior.

struct ParsedControl {
  bool responseIsBinary = false;
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
  const double* x_test = NULL;
  size_t numTestObservations = 0;
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
  double kDf = 1.25;
  double kScale = HUGE_VAL;
  bool sigmaIsFixed = false;
  double fixedSigmaSq = 1.0; // the residual variance fixed() holds
  double sigmaDf = 3.0;
  double sigmaQuantile = 0.9;
  // the chi-squared prior's scale before anchoring to the sigma estimate
  double sigmaRawScale = 1.0;
  // a linear or gp node prior's designated covariate columns (0-based);
  // empty for the constant leaf
  std::vector<size_t> leafCovariateColumns;
  // gp node priors only: selects the function-valued leaf model over the
  // linear one; lengthscales empty for the median-distance heuristic
  bool gpLeaves = false;
  std::vector<double> gpLengthscales;
  size_t gpMaxLeafSize = 256;
};

void parseControl(ParsedControl& control, SEXP controlExpr) {
  SEXP slotExpr = Rf_getAttrib(controlExpr, Rf_install("binary"));
  control.responseIsBinary =
    rc_getBool(slotExpr, "binary response signifier", RC_LENGTH | RC_GEQ,
               rc_asRLength(1), RC_END);

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("verbose"));
  control.verbose = rc_getBool(slotExpr, "verbose", RC_LENGTH | RC_GEQ,
                               rc_asRLength(1), RC_END);

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("keepTrainingFits"));
  control.keepTrainingFits =
    rc_getBool(slotExpr, "keep training fits", RC_LENGTH | RC_EQ,
               rc_asRLength(1), RC_END);

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("useQuantiles"));
  control.useQuantiles = rc_getBool(slotExpr, "use quantiles",
                                    RC_LENGTH | RC_EQ, rc_asRLength(1),
                                    RC_END);

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("keepTrees"));
  control.keepTrees = rc_getBool(slotExpr, "keep trees", RC_LENGTH | RC_EQ,
                                 rc_asRLength(1), RC_END);

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.samples"));
  control.defaultNumSamples = static_cast<size_t>(
    rc_getInt(slotExpr, "number of samples", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END));

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.trees"));
  control.numTrees = static_cast<size_t>(
    rc_getInt(slotExpr, "number of trees", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END));

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.chains"));
  control.numChains = static_cast<size_t>(
    rc_getInt(slotExpr, "number of chains", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END));

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.threads"));
  control.numThreads = static_cast<size_t>(
    rc_getInt(slotExpr, "number of threads", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END));

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.thin"));
  control.treeThinningRate = static_cast<uint32_t>(
    rc_getInt(slotExpr, "tree thinning rate", RC_LENGTH | RC_EQ,
              rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END));

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("printEvery"));
  int i_temp = rc_getInt(slotExpr, "print every", RC_LENGTH | RC_EQ,
                         rc_asRLength(1), RC_VALUE | RC_GEQ, 1,
                         RC_NA | RC_YES, RC_END);
  if (i_temp != NA_INTEGER) control.printEvery = static_cast<uint32_t>(i_temp);

  slotExpr = Rf_getAttrib(controlExpr, Rf_install("printCutoffs"));
  i_temp = rc_getInt(slotExpr, "print cutoffs", RC_LENGTH | RC_EQ,
                     rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES,
                     RC_END);
  if (i_temp == NA_INTEGER) i_temp = 0;
  control.printCutoffs = static_cast<uint32_t>(i_temp);

  // rngKind and rngNormalKind are classic-only; the R side refuses them
  // before creation, so they are not read here
  slotExpr = Rf_getAttrib(controlExpr, Rf_install("rngSeed"));
  if (rc_getLength(slotExpr) != 1)
    Rf_error("slot 'rngSeed' must be of length 1");
  control.haveRngSeed = INTEGER(slotExpr)[0] != NA_INTEGER;
  if (control.haveRngSeed)
    control.rngSeed = static_cast<std::uint_least32_t>(INTEGER(slotExpr)[0]);
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
  SEXP dimExpr = Rf_getAttrib(matrixExpr, Rf_install("Dim"));
  if (!Rf_isInteger(dimExpr) || rc_getLength(dimExpr) != 2)
    Rf_error("malformed sparse predictor matrix");
  if (static_cast<size_t>(INTEGER(dimExpr)[0]) != numObservations)
    Rf_error("number of rows of 'x' must equal length of 'y'");
  result.numColumns = static_cast<size_t>(INTEGER(dimExpr)[1]);

  SEXP pointersExpr = Rf_getAttrib(matrixExpr, Rf_install("p"));
  SEXP rowsExpr = Rf_getAttrib(matrixExpr, Rf_install("i"));
  SEXP valuesExpr = Rf_getAttrib(matrixExpr, Rf_install("x"));
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
  return result;
}

void parseData(ParsedData& data, SEXP dataExpr) {
  SEXP slotExpr = Rf_getAttrib(dataExpr, Rf_install("y"));
  if (!Rf_isReal(slotExpr)) Rf_error("y must be of type real");
  if (rc_getLength(slotExpr) <= 0)
    Rf_error("length of y must be greater than 0");
  data.y = REAL(slotExpr);
  data.numObservations = rc_getLength(slotExpr);

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("x"));
  if (Rf_inherits(slotExpr, "dgCMatrix")) {
    CscSlots csc = parseCscMatrix(slotExpr, data.numObservations);
    data.numPredictors = csc.numColumns;
    data.xIsSparse = true;
    data.cscColumnPointers = csc.pointers;
    data.cscRowIndices = csc.rows;
    data.cscValues = csc.values;
    data.x = NULL;
  } else if (Rf_inherits(slotExpr, "dbartsMixedMatrix")) {
    // the internal mixed container (R/mixedMatrix.R): a dense matrix, a
    // dgCMatrix, and a 1-based map - positive k names dense column k,
    // negative -k sparse column k, which is the engine's ~(k - 1)
    SEXP denseExpr = getListElement(slotExpr, "dense");
    SEXP sparseExpr = getListElement(slotExpr, "sparse");
    SEXP mapExpr = getListElement(slotExpr, "map");
    if (!Rf_inherits(sparseExpr, "dgCMatrix") || !Rf_isInteger(mapExpr) ||
        rc_getLength(mapExpr) == 0)
      Rf_error("malformed mixed predictor container");
    CscSlots csc = parseCscMatrix(sparseExpr, data.numObservations);

    size_t numDenseColumns = 0;
    if (!Rf_isNull(denseExpr)) {
      if (!Rf_isReal(denseExpr))
        Rf_error("malformed mixed predictor container");
      rc_assertDimConstraints(denseExpr, "dimensions of x",
                              RC_LENGTH | RC_EQ, rc_asRLength(2),
                              RC_VALUE | RC_EQ,
                              static_cast<int>(data.numObservations), RC_END);
      numDenseColumns = static_cast<size_t>(
        INTEGER(Rf_getAttrib(denseExpr, R_DimSymbol))[1]);
      data.mixedDenseValues = REAL(denseExpr);
    }

    data.numPredictors = rc_getLength(mapExpr);
    const int* map = INTEGER(mapExpr);
    data.columnSources.resize(data.numPredictors);
    for (size_t j = 0; j < data.numPredictors; ++j) {
      if (map[j] > 0 && static_cast<size_t>(map[j]) <= numDenseColumns) {
        data.columnSources[j] = map[j] - 1;
      } else if (map[j] < 0 &&
                 static_cast<size_t>(-map[j]) <= csc.numColumns) {
        data.columnSources[j] = map[j];
      } else {
        Rf_error("malformed mixed predictor container");
      }
    }
    data.xIsMixed = true;
    data.cscColumnPointers = csc.pointers;
    data.cscRowIndices = csc.rows;
    data.cscValues = csc.values;
    data.x = NULL;
  } else {
    if (!Rf_isReal(slotExpr)) Rf_error("x must be of type real");
    rc_assertDimConstraints(slotExpr, "dimensions of x", RC_LENGTH | RC_EQ,
                            rc_asRLength(2), RC_VALUE | RC_EQ,
                            static_cast<int>(data.numObservations), RC_END);
    int* dims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
    data.x = REAL(slotExpr);
    data.numPredictors = static_cast<size_t>(dims[1]);
  }

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("varTypes"));
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
    for (size_t j = 0; j < data.numPredictors; ++j)
      if (data.columnTypes[j] == bartcore::ColumnType::categorical &&
          data.columnSources[j] < 0)
        Rf_error("sparse predictor columns must be ordinal");
  }

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("x.test"));
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.x_test = NULL;
    data.numTestObservations = 0;
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

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("weights"));
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.weights = NULL;
  } else {
    rc_assertDoubleConstraints(slotExpr, "weights", RC_LENGTH | RC_EQ,
                               rc_asRLength(data.numObservations), RC_END);
    data.weights = REAL(slotExpr);
  }

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("offset"));
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.offset = NULL;
  } else {
    rc_assertDoubleConstraints(slotExpr, "offset", RC_LENGTH | RC_EQ,
                               rc_asRLength(data.numObservations), RC_END);
    data.offset = REAL(slotExpr);
  }

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("offset.test"));
  if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) ||
      rc_getLength(slotExpr) == 0) {
    data.testOffset = NULL;
  } else {
    rc_assertDoubleConstraints(slotExpr, "test offset", RC_LENGTH | RC_EQ,
                               rc_asRLength(data.numTestObservations),
                               RC_END);
    data.testOffset = REAL(slotExpr);
  }

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("sigma"));
  data.sigmaEstimate = rc_getDouble(
    slotExpr, "sigma estimate", RC_LENGTH | RC_EQ, rc_asRLength(1),
    RC_NA | RC_YES, RC_VALUE | RC_GT, 0.0, RC_END);

  slotExpr = Rf_getAttrib(dataExpr, Rf_install("n.cuts"));
  rc_assertIntConstraints(slotExpr, "maximum number of cuts",
                          RC_LENGTH | RC_EQ,
                          rc_asRLength(data.numPredictors), RC_END);
  int* i_maxNumCuts = INTEGER(slotExpr);
  data.maxNumCuts.resize(data.numPredictors);
  for (size_t j = 0; j < data.numPredictors; ++j)
    data.maxNumCuts[j] = static_cast<uint32_t>(i_maxNumCuts[j]);
}

void parseModel(ParsedModel& model, SEXP modelExpr, size_t numPredictors) {
  SEXP slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.birth_death"));
  model.birthOrDeathProbability = rc_getDouble(
    slotExpr, "probability of birth/death rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.swap"));
  model.swapProbability = rc_getDouble(
    slotExpr, "probability of swap rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.change"));
  model.changeProbability = rc_getDouble(
    slotExpr, "probability of change rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  if (std::fabs(model.birthOrDeathProbability + model.swapProbability +
                model.changeProbability - 1.0) >= 1.0e-10)
    Rf_error("rule proposal probabilities must sum to 1.0");

  slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.birth"));
  model.birthProbability = rc_getDouble(
    slotExpr, "probability of birth in birth/death rule", RC_LENGTH | RC_EQ,
    rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);

  slotExpr = Rf_getAttrib(modelExpr, Rf_install("node.scale"));
  model.nodeScale = rc_getDouble(
    slotExpr, "scale of node prior", RC_LENGTH | RC_EQ, rc_asRLength(1),
    RC_VALUE | RC_GT, 0.0, RC_END);

  // linear and gp node priors designate leaf covariate columns, resolved
  // R-side to 1-based model matrix indices; every other node prior is the
  // constant leaf and carries nothing beyond node.scale/node.hyperprior.
  // gp priors add per-column lengthscales (NULL for the median-distance
  // heuristic; the R side validates and recycles) and the leaf-size cap.
  SEXP nodePriorExpr = Rf_getAttrib(modelExpr, Rf_install("node.prior"));
  bool isLinearPrior = !Rf_isNull(nodePriorExpr) &&
                       Rf_inherits(nodePriorExpr, "dbartsLinearPrior");
  bool isGPPrior = !Rf_isNull(nodePriorExpr) &&
                   Rf_inherits(nodePriorExpr, "dbartsGPPrior");
  if (isLinearPrior || isGPPrior) {
    SEXP columnsExpr = Rf_getAttrib(nodePriorExpr, Rf_install("columns"));
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
  }
  if (isGPPrior) {
    model.gpLeaves = true;
    // a NULL slot arrives as S4's pseudo-NULL symbol, not R_NilValue, so
    // test positively for the resolved numeric vector
    SEXP lengthscaleExpr =
      Rf_getAttrib(nodePriorExpr, Rf_install("lengthscale"));
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
      Rf_getAttrib(nodePriorExpr, Rf_install("max.leaf.size"));
    int maxLeafSize = rc_getInt(
      maxLeafSizeExpr, "gp node prior maximum leaf size", RC_LENGTH | RC_EQ,
      rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END);
    model.gpMaxLeafSize = static_cast<size_t>(maxLeafSize);
  }

  SEXP priorExpr = Rf_getAttrib(modelExpr, Rf_install("tree.prior"));
  slotExpr = Rf_getAttrib(priorExpr, Rf_install("power"));
  model.power = rc_getDouble(slotExpr, "tree prior power", RC_LENGTH | RC_EQ,
                             rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);

  slotExpr = Rf_getAttrib(priorExpr, Rf_install("base"));
  model.base = rc_getDouble(slotExpr, "tree prior base", RC_LENGTH | RC_EQ,
                            rc_asRLength(1), RC_VALUE | RC_GT, 0.0,
                            RC_VALUE | RC_LT, 1.0, RC_END);

  slotExpr = Rf_getAttrib(priorExpr, Rf_install("splitProbabilities"));
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

  priorExpr = Rf_getAttrib(modelExpr, Rf_install("node.hyperprior"));
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

  priorExpr = Rf_getAttrib(modelExpr, Rf_install("resid.prior"));
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

// family selects the response model for binary responses: "" or "probit"
// give the classic probit latents, "logistic" the Polya-Gamma sampler.
// Continuous responses are gaussian and accept only "" or "gaussian".
bartcore::ResponseFamily resolveFamily(const ParsedControl& control,
                                       const char* familyName) {
  if (control.responseIsBinary) {
    if (std::strcmp(familyName, "logistic") == 0)
      return bartcore::ResponseFamily::logistic;
    if (familyName[0] == '\0' || std::strcmp(familyName, "probit") == 0)
      return bartcore::ResponseFamily::probit;
    Rf_error("unrecognized response family for a binary response");
  }
  if (familyName[0] != '\0' && std::strcmp(familyName, "gaussian") != 0)
    Rf_error("response families other than gaussian require a binary "
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

void validateCategoricalPredictors(const ParsedData& data) {
  if (data.xIsSparse) return;  // parseData enforced all-ordinal
  // for mixed containers parseData enforced dense backing per categorical
  // column, so rawTrainingColumn serves every column validated here
  for (size_t j = 0; j < data.numPredictors; ++j) {
    if (data.columnTypes[j] != bartcore::ColumnType::categorical) continue;
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
      const double* column = rawTrainingColumn(data, j);
      double maxValue = 0.0;
      for (size_t i = 0; i < data.numObservations; ++i) {
        double value = column[i];
        if (!bartcore::isNA(value) && value > maxValue) maxValue = value;
      }
      for (size_t i = 0; i < data.numTestObservations; ++i) {
        double value = data.x_test[i + j * data.numTestObservations];
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
  }
  options.splitProbabilities = model.splitProbabilities; // copied by ctor
  options.leafCovariateColumns = model.leafCovariateColumns.empty()
    ? NULL : model.leafCovariateColumns.data();  // consumed at construction
  options.numLeafCovariates = model.leafCovariateColumns.size();
  options.gpLeaves = model.gpLeaves;
  options.gpLengthscales = model.gpLengthscales.empty()
    ? NULL : model.gpLengthscales.data();  // consumed at construction
  options.gpMaxLeafSize = model.gpMaxLeafSize;

  // the generic slot parse above reads the CGM structure a DART prior
  // contains; the Dirichlet configuration comes off the R object directly
  SEXP treePriorExpr = Rf_getAttrib(modelExpr, Rf_install("tree.prior"));
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

// A sampler created over a data handle holds no raw predictor values, so
// the raw-x mutation surface has nothing to work from; a CSC-built sampler
// holds only borrowed slices and refuses the same surface by design
// (docs/design/sparse-columns.md).
void refuseViewSampler(const bartcore::SamplerBase& sampler,
                       const char* caller) {
  if (sampler.data().x != NULL) return;
  if (sampler.data().builtFromCsc)
    Rf_error("%s: sparse predictors fix the design at creation; make a new "
             "sampler instead", caller);
  Rf_error("%s requires a sampler that owns its predictors; data-handle "
           "views hold none", caller);
}

// Unlike views, CSC-built samplers re-quantize from their borrowed slices,
// so cut installation and state restore stay available on them.
void refuseViewSamplerOnly(const bartcore::SamplerBase& sampler,
                           const char* caller) {
  if (sampler.data().x == NULL && !sampler.data().builtFromCsc)
    Rf_error("%s requires a sampler that owns its predictors; data-handle "
             "views hold none", caller);
}

// A built column store (cuts + codes) shared by row-subset view samplers
// (public-surface.md section 5; internal). The external pointer's
// protection slot pins the data expression whose x the store borrows.
struct DataHandle {
  bartcore::ColumnStore store;
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

// family selects the response model for binary responses: "" or "probit"
// give the classic probit latents, "logistic" the Polya-Gamma sampler.
// Continuous responses are gaussian and accept only "" or "gaussian".
BartcoreHolder* createHolder(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                             const char* familyName) {
  ParsedControl control;
  ParsedData data;
  ParsedModel model;
  bool sigmaIsFixed;
  bartcore::ResponseFamily family = parseSamplerSpecification(
    controlExpr, modelExpr, dataExpr, familyName, control, model, data,
    sigmaIsFixed);
  validateCategoricalPredictors(data);

  bartcore::SamplerOptions options =
    optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);

  // grouped random intercepts (rbart_vi's in-core path) arrive on an
  // internal control attribute; the chains copy the indices at construction
  std::vector<std::uint32_t> groupIndices;
  applyGroupAttribute(controlExpr, data.numObservations, options,
                      groupIndices);

  std::vector<ext_rng*> rngs = createChainRngs(control, options.numChains);

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
    sampler->setTestPredictors(data.x_test, data.numTestObservations);
    sampler->setTestOffset(data.testOffset);
  }

  if (control.verbose) printInitialSummary(control, model, data, *sampler);

  return new BartcoreHolder{std::move(sampler), std::move(rngs),
                            control.keepTrainingFits, {}, {}, {}, {}};
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
  BartcoreHolder* holder = bartcore_bridge::createHolder(
    controlExpr, modelExpr, dataExpr, familyName);

  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(holder, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  UNPROTECT(2);
  return result;
}

// Builds the two-layer store (cuts + codes) once for sharing across
// row-subset samplers: control contributes useQuantiles, data contributes
// x, the column types, and n.cuts. Internal, with no serialization;
// see public-surface.md section 5.
SEXP bartcore_createDataHandle(SEXP controlExpr, SEXP dataExpr) {
  ParsedControl control;
  ParsedData data;
  parseControl(control, controlExpr);
  parseData(data, dataExpr);
  validateCategoricalPredictors(data);

  DataHandle* handle = new DataHandle;
  if (data.xIsMixed) {
    handle->store.buildMixed(data.mixedDenseValues, data.cscColumnPointers,
                             data.cscRowIndices, data.cscValues,
                             data.columnSources.data(), data.numObservations,
                             data.numPredictors, data.maxNumCuts.data(), 0,
                             control.useQuantiles,
                             data.anyCategorical ? data.columnTypes.data()
                                                 : NULL);
  } else if (data.xIsSparse) {
    handle->store.buildFromCsc(data.cscColumnPointers, data.cscRowIndices,
                               data.cscValues, data.numObservations,
                               data.numPredictors, data.maxNumCuts.data(), 0,
                               control.useQuantiles);
  } else {
    handle->store.build(data.x, data.numObservations, data.numPredictors,
                        data.maxNumCuts.data(), control.useQuantiles,
                        data.anyCategorical ? data.columnTypes.data() : NULL);
  }

  SEXP result = PROTECT(R_MakeExternalPtr(handle, R_NilValue, dataExpr));
  R_RegisterCFinalizerEx(result, dataHandleFinalizer,
                         static_cast<Rboolean>(FALSE));
  UNPROTECT(1);
  return result;
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
                               SEXP familyExpr) {
  const bartcore::ColumnStore& parent =
    dataHandleFromExpression(handleExpr).store;
  const char* familyName =
    Rf_isNull(familyExpr) ? "" : CHAR(STRING_ELT(familyExpr, 0));

  ParsedControl control;
  ParsedData data;
  ParsedModel model;
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

  std::vector<size_t> trainRows(numTrainRows), testRows(numTestRows);
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

  std::vector<double> response(numTrainRows);
  for (size_t i = 0; i < numTrainRows; ++i)
    response[i] = data.y[trainRows[i]];
  std::vector<double> weights, offset, testOffset;
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

  // a linear node prior's designated columns have the view gather their raw
  // values, with standardization constants from the handle's full data - the
  // same calibration inheritance as the copied cut grid
  bartcore::ColumnStore store;
  store.buildFromParent(parent, trainRows.data(), numTrainRows,
                        testRows.data(), numTestRows,
                        model.leafCovariateColumns.empty()
                          ? NULL : model.leafCovariateColumns.data(),
                        model.leafCovariateColumns.size());

  bartcore::SamplerOptions options =
    optionsFromParsed(control, model, data, modelExpr, sigmaIsFixed);
  std::vector<ext_rng*> rngs = createChainRngs(control, options.numChains);

  std::unique_ptr<bartcore::SamplerBase> sampler =
    bartcore::createSamplerOverStore(
      std::move(store), response.data(),
      data.weights != NULL ? weights.data() : NULL,
      data.offset != NULL ? offset.data() : NULL, family, data.sigmaEstimate,
      model.sigmaDf, model.sigmaRawScale, options, rngs.data());
  if (sampler == NULL) {
    // R-side resolution validates first, so only an invariant breach lands
    for (ext_rng* rng : rngs) if (rng != NULL) ext_rng_destroy(rng);
    Rf_error("invalid leaf covariate designation");
  }

  BartcoreHolder* holder = new BartcoreHolder{
    std::move(sampler), std::move(rngs), control.keepTrainingFits,
    {}, {}, {}, {}};
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

  if (numSamples == 0) {
    bartcore::Results empty;
    GetRNGstate();
    bool cancelled = sampler.run(numBurnIn, 0, empty, bartcore_userInterrupted);
    PutRNGstate();
    if (cancelled) Rf_error("sampler run interrupted");
    return R_NilValue;
  }

  // several chains add a trailing chain dimension, as the classic engine's
  // results do. Every column roots in the protected container the moment it
  // is allocated (installResult), so there is no hand-counted PROTECT stack
  // to keep in sync with the slot list.
  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 8));
  SEXP sigmaExpr = installResult(
    resultExpr, 0,
    numChains == 1
      ? Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples))
      : Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  SEXP trainExpr = installResult(
    resultExpr, 1,
    !holder.keepTrainingFits
      ? R_NilValue
      : numChains == 1
        ? Rf_allocMatrix(REALSXP, static_cast<int>(numObservations),
                         numSamplesInt)
        : Rf_alloc3DArray(REALSXP, static_cast<int>(numObservations),
                          numSamplesInt, numChainsInt));
  SEXP testExpr = installResult(
    resultExpr, 2,
    numTestObservations == 0
      ? R_NilValue
      : numChains == 1
        ? Rf_allocMatrix(REALSXP, static_cast<int>(numTestObservations),
                         numSamplesInt)
        : Rf_alloc3DArray(REALSXP, static_cast<int>(numTestObservations),
                          numSamplesInt, numChainsInt));
  SEXP varcountExpr = installResult(
    resultExpr, 3,
    numChains == 1
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

  std::vector<std::uint32_t> variableCounts(numPredictors * numSamples *
                                            numChains);

  bartcore::Results results;
  results.sigma = REAL(sigmaExpr);
  results.trainingFits = holder.keepTrainingFits ? REAL(trainExpr) : NULL;
  results.testFits = numTestObservations > 0 ? REAL(testExpr) : NULL;
  results.variableCounts = variableCounts.data();
  results.k = sampler.kIsSampled() ? REAL(kExpr) : NULL;
  results.splitProbabilities = sampler.usesDart() ? REAL(varprobsExpr) : NULL;
  results.tau = numGroups > 0 ? REAL(tauExpr) : NULL;
  results.groupEffects = numGroups > 0 ? REAL(ranefExpr) : NULL;

  GetRNGstate();
  bool cancelled = sampler.run(numBurnIn, numSamples, results,
                               bartcore_userInterrupted);
  PutRNGstate();
  if (cancelled) Rf_error("sampler run interrupted");

  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < numPredictors * numSamples * numChains; ++i)
    varcountOut[i] = static_cast<int>(variableCounts[i]);

  // named as the classic engine's run results are, so the engines are
  // drop-in replacements for each other; varprobs, tau, and ranef are
  // bartcore extensions. The names vector roots through the container's
  // attribute before the mkChar allocations fill it.
  SEXP namesExpr = Rf_allocVector(STRSXP, 8);
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));
  SET_STRING_ELT(namesExpr, 5, Rf_mkChar("varprobs"));
  SET_STRING_ELT(namesExpr, 6, Rf_mkChar("tau"));
  SET_STRING_ELT(namesExpr, 7, Rf_mkChar("ranef"));

  UNPROTECT(1);
  return resultExpr;
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
  refuseViewSampler(sampler, "bartcore_setData");
  if (sampler.numGroups() > 0)
    Rf_error("grouped random effects fix the data at creation; make a new "
             "sampler instead");

  if (!Rf_inherits(dataExpr, "dbartsData"))
    Rf_error("'data' argument to bartcore_setData not of class 'dbartsData'");

  ParsedData data;
  parseData(data, dataExpr);
  if (data.xIsSparse || data.xIsMixed)
    Rf_error("bartcore setData requires a dense predictor matrix; sparse "
             "predictors fix the design at creation");

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

  // everything the sampler borrows now comes from the new data object
  retain(ptrExpr, PROT_DATA, dataExpr);
  retain(ptrExpr, PROT_RESPONSE, R_NilValue);
  retain(ptrExpr, PROT_OFFSET, R_NilValue);
  retain(ptrExpr, PROT_PREDICTORS, R_NilValue);
  retain(ptrExpr, PROT_TEST_PREDICTORS, R_NilValue);

  return R_NilValue;
}

SEXP bartcore_setTestPredictor(SEXP ptrExpr, SEXP xTestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (Rf_isNull(xTestExpr)) {
    // removal: back to the no-test-data state, offset included
    holder.sampler->setTestPredictors(NULL, 0);
    holder.sampler->setTestOffset(NULL);
    retain(ptrExpr, PROT_TEST_PREDICTORS, R_NilValue);
    retain(ptrExpr, PROT_TEST_OFFSET, R_NilValue);
    return R_NilValue;
  }
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
  holder.sampler->setTestPredictors(REAL(xTestExpr), numTestObservations);
  retain(ptrExpr, PROT_TEST_PREDICTORS, xTestExpr);
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
    retain(ptrExpr, PROT_TEST_PREDICTORS, R_NilValue);
    retain(ptrExpr, PROT_TEST_OFFSET, R_NilValue);
    return R_NilValue;
  }
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

  holder.sampler->setTestPredictors(REAL(xTestExpr), numTestObservations);
  holder.sampler->setTestOffset(Rf_isNull(offsetExpr) ? NULL
                                                      : REAL(offsetExpr));
  retain(ptrExpr, PROT_TEST_PREDICTORS, xTestExpr);
  retain(ptrExpr, PROT_TEST_OFFSET, offsetExpr);
  return R_NilValue;
}

// Case weights, like the classic engine's setWeights: a pointer swap with
// nothing rescaled; refused for binary responses, whose reference-engine
// weighting was incorrect and was stripped rather than ported.
SEXP bartcore_setWeights(SEXP ptrExpr, SEXP weightsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
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

  ParsedModel model;
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

  bool isGaussian = sampler.family() == bartcore::ResponseFamily::gaussian;

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
  if (isGaussian) {
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
  refuseViewSampler(*holder.sampler, "bartcore_setPredictor");
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
  if (result == bartcore::PredictorUpdateResult::accepted)
    retain(ptrExpr, PROT_PREDICTORS, xExpr);  // installed by pointer swap
  return Rf_ScalarLogical(
    result == bartcore::PredictorUpdateResult::accepted ? TRUE : FALSE);
}

SEXP bartcore_updatePredictor(SEXP ptrExpr, SEXP xExpr, SEXP columnsExpr,
                              SEXP forceUpdateExpr, SEXP updateCutPointsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refuseViewSampler(*holder.sampler, "bartcore_updatePredictor");
  size_t numObservations = holder.sampler->numObservations();
  size_t numPredictors = holder.sampler->numPredictors();

  size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
  if (numColumns == 0 ||
      static_cast<size_t>(Rf_xlength(xExpr)) != numObservations * numColumns)
    Rf_error("bartcore_updatePredictor requires numObservations values per column");

  std::vector<size_t> columns(numColumns);
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
}

SEXP bartcore_setCutPoints(SEXP ptrExpr, SEXP cutPointsExpr,
                           SEXP columnsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refuseViewSamplerOnly(*holder.sampler, "bartcore_setCutPoints");
  size_t numPredictors = holder.sampler->numPredictors();

  size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
  if (numColumns == 0 ||
      static_cast<size_t>(Rf_xlength(cutPointsExpr)) != numColumns)
    Rf_error("bartcore_setCutPoints requires one cut point vector per column");

  std::vector<const double*> cutPoints(numColumns);
  std::vector<std::uint32_t> numCutPoints(numColumns);
  std::vector<size_t> columns(numColumns);
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
                               columns.data(), numColumns);
  return R_NilValue;
}

SEXP bartcore_updatePredictorPerObservation(SEXP ptrExpr, SEXP xExpr,
                                            SEXP columnExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  refuseViewSampler(*holder.sampler, "bartcore_updatePredictorPerObservation");
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
  if (!treesAreValid)
    Rf_error("bartcore updatePredictorPerObservation produced a tree with an "
             "empty leaf");

  SEXP result = PROTECT(
    Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(numObservations)));
  for (size_t i = 0; i < numObservations; ++i)
    LOGICAL(result)[i] = installed[i] ? TRUE : FALSE;
  UNPROTECT(1);
  return result;
}

SEXP bartcore_updatePredictorPerObservationJointly(SEXP ptrsExpr, SEXP xExpr,
                                                   SEXP columnsExpr) {
  size_t numSamplers = static_cast<size_t>(Rf_xlength(ptrsExpr));
  if (numSamplers == 0 ||
      static_cast<size_t>(Rf_xlength(columnsExpr)) != numSamplers)
    Rf_error("bartcore_updatePredictorPerObservationJointly requires one "
             "column per sampler");

  std::vector<bartcore::SamplerBase*> samplers(numSamplers);
  std::vector<size_t> columns(numSamplers);
  for (size_t k = 0; k < numSamplers; ++k) {
    BartcoreHolder& holder(
      holderFromExpression(VECTOR_ELT(ptrsExpr, static_cast<R_xlen_t>(k))));
    refuseViewSampler(*holder.sampler,
                      "bartcore_updatePredictorPerObservationJointly");
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
    samplers.data(), numSamplers, REAL(xExpr), columns.data(), installed.get());
  PutRNGstate();

  if (!treesAreValid)
    Rf_error("bartcore updatePredictorPerObservationJointly produced a tree "
             "with an empty leaf");

  SEXP result = PROTECT(
    Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(numObservations)));
  for (size_t i = 0; i < numObservations; ++i)
    LOGICAL(result)[i] = installed[i] ? TRUE : FALSE;
  UNPROTECT(1);
  return result;
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

// flattened trees as three parallel R vectors: concatenated 1-based-with-(-1)
// variables, values, and per-tree node counts
static void storeFlatTrees(SEXP chainExpr, int variablesSlot, int valuesSlot,
                           int sizesSlot, int flagsSlot,
                           const std::vector<std::vector<bartcore::FlatNode>>& trees) {
  R_xlen_t totalNumNodes = 0;
  for (const std::vector<bartcore::FlatNode>& tree : trees)
    totalNumNodes += static_cast<R_xlen_t>(tree.size());

  SET_VECTOR_ELT(chainExpr, variablesSlot,
                 Rf_allocVector(INTSXP, totalNumNodes));
  SET_VECTOR_ELT(chainExpr, valuesSlot, Rf_allocVector(REALSXP, totalNumNodes));
  SET_VECTOR_ELT(chainExpr, sizesSlot,
                 Rf_allocVector(INTSXP, static_cast<R_xlen_t>(trees.size())));
  SET_VECTOR_ELT(chainExpr, flagsSlot, Rf_allocVector(RAWSXP, totalNumNodes));

  int* variables = INTEGER(VECTOR_ELT(chainExpr, variablesSlot));
  double* values = REAL(VECTOR_ELT(chainExpr, valuesSlot));
  int* sizes = INTEGER(VECTOR_ELT(chainExpr, sizesSlot));
  Rbyte* flags = RAW(VECTOR_ELT(chainExpr, flagsSlot));
  R_xlen_t offset = 0;
  for (size_t t = 0; t < trees.size(); ++t) {
    sizes[t] = static_cast<int>(trees[t].size());
    for (const bartcore::FlatNode& node : trees[t]) {
      variables[offset] = node.variable >= 0 ? node.variable + 1
                                             : node.variable;
      values[offset] = node.value;
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

// wide categorical mask channels as one raw vector, each word written
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
    if (j < store.numPredictors && store.columnHasWideMask(j))
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
// up C++ state first
static bool readFlatTrees(SEXP variablesExpr, SEXP valuesExpr, SEXP sizesExpr,
                          SEXP flagsExpr,
                          std::vector<std::vector<bartcore::FlatNode>>& trees,
                          const char** errorMessage) {
  if (!Rf_isInteger(variablesExpr) || !Rf_isReal(valuesExpr) ||
      !Rf_isInteger(sizesExpr) ||
      Rf_xlength(variablesExpr) != Rf_xlength(valuesExpr) ||
      (!Rf_isNull(flagsExpr) &&
       (TYPEOF(flagsExpr) != RAWSXP ||
        Rf_xlength(flagsExpr) != Rf_xlength(variablesExpr)))) {
    *errorMessage = "malformed trees in bartcore state";
    return false;
  }
  const Rbyte* flags = Rf_isNull(flagsExpr) ? NULL : RAW(flagsExpr);
  R_xlen_t numTrees = Rf_xlength(sizesExpr);
  const int* variables = INTEGER(variablesExpr);
  const double* values = REAL(valuesExpr);
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
      node.value = values[offset];
      node.flags = flags != NULL ? static_cast<std::uint8_t>(flags[offset]) : 0;
      ++offset;
    }
  }
  return true;
}

SEXP bartcore_storeState(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  return bartcore_bridge::storeState(*holder.sampler);
}

SEXP bartcore_setState(SEXP ptrExpr, SEXP stateExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  // restoring cut points re-quantizes from raw values, which views lack
  refuseViewSamplerOnly(*holder.sampler, "bartcore_setState");
  bartcore_bridge::setState(*holder.sampler, stateExpr);
  return R_NilValue;
}

// Fits for new data on the original response scale (binary responses give
// the latent scale, as the classic engine does). With keepTrees the saved
// trees produce numTestObservations x numSamples (x numChains) fits; without,
// the live trees produce a single set per chain. offset, when non-null, is
// added to every sample's fits.
SEXP bartcore_predict(SEXP ptrExpr, SEXP xTestExpr, SEXP offsetExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

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

  int numTestInt = static_cast<int>(numTestObservations);
  SEXP resultExpr;
  if (capacity > 0) {
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
                       SEXP treeNumsExpr, SEXP currentExpr, SEXP newdataExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  bool useLiveTrees = Rf_asLogical(currentExpr) == TRUE;
  bool useSaved = sampler.savedTreeCapacity() > 0 && !useLiveTrees;

  std::vector<size_t> chainIndices(
    static_cast<size_t>(Rf_xlength(chainNumsExpr)));
  for (size_t i = 0; i < chainIndices.size(); ++i) {
    int chainNum = INTEGER(chainNumsExpr)[i];
    if (chainNum < 1) Rf_error("bartcore_getTrees chain number out of range");
    chainIndices[i] = static_cast<size_t>(chainNum - 1);
  }
  std::vector<size_t> sampleIndices;
  if (useSaved) {
    sampleIndices.resize(static_cast<size_t>(Rf_xlength(sampleNumsExpr)));
    for (size_t i = 0; i < sampleIndices.size(); ++i) {
      int sampleNum = INTEGER(sampleNumsExpr)[i];
      if (sampleNum < 1)
        Rf_error("bartcore_getTrees sample number out of range");
      sampleIndices[i] = static_cast<size_t>(sampleNum - 1);
    }
  }
  std::vector<size_t> treeIndices(
    static_cast<size_t>(Rf_xlength(treeNumsExpr)));
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

  return bartcore_bridge::getTrees(
    sampler, chainIndices.data(), chainIndices.size(), sampleIndices.data(),
    sampleIndices.size(), treeIndices.data(), treeIndices.size(),
    useLiveTrees, newdata, newdataNumRows, "bartcore_getTrees");
}

SEXP bartcore_printTrees(SEXP ptrExpr, SEXP chainNumsExpr, SEXP sampleNumsExpr,
                         SEXP treeNumsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  size_t capacity = sampler.savedTreeCapacity();

  std::vector<size_t> chainIndices, sampleIndices, treeIndices;
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

// Bumped by any change to the layout storeState/setState exchange (tree
// encoding, per-chain slots, attributes). States load only within a format
// version; there are no migration shims.
static const int stateFormatVersion = 1;

SEXP storeState(bartcore::SamplerBase& sampler) {
  bartcore::SamplerStateData state;
  sampler.getState(state);

  size_t numChains = state.chains.size();
  size_t numObservations = sampler.numObservations();

  enum {
    SLOT_TREE_VARS = 0, SLOT_TREE_VALUES, SLOT_TREE_SIZES, SLOT_TREE_FLAGS,
    SLOT_TREE_PARAMS, SLOT_TREE_MASKS, SLOT_SAVED_VARS,
    SLOT_SAVED_VALUES, SLOT_SAVED_SIZES, SLOT_SAVED_FLAGS, SLOT_SAVED_PARAMS,
    SLOT_SAVED_MASKS,
    SLOT_TOTAL_FITS, SLOT_INDICES,
    SLOT_SIGMA, SLOT_K, SLOT_FIT_SCALE, SLOT_LATENTS,
    SLOT_RANEF, SLOT_TAU,
    SLOT_DART_PROBABILITIES, SLOT_DART_ALPHA, SLOT_DART_UPDATES_SKIPPED,
    SLOT_RNG_STATE, SLOT_COUNT
  };
  static const char* slotNames[SLOT_COUNT] = {
    "tree.vars", "tree.values", "tree.sizes", "tree.flags", "tree.params",
    "tree.masks", "saved.vars",
    "saved.values", "saved.sizes", "saved.flags", "saved.params",
    "saved.masks",
    "total.fits", "indices",
    "sigma", "k", "fit.scale",
    "latents", "ranef", "tau",
    "dart.probabilities", "dart.alpha", "dart.updates.skipped",
    "rng.state"
  };

  SEXP resultExpr =
    PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(numChains)));
  SEXP slotNamesExpr = PROTECT(Rf_allocVector(STRSXP, SLOT_COUNT));
  for (int slot = 0; slot < SLOT_COUNT; ++slot)
    SET_STRING_ELT(slotNamesExpr, slot, Rf_mkChar(slotNames[slot]));

  for (size_t c = 0; c < numChains; ++c) {
    const bartcore::ChainStateData& chainState(state.chains[c]);
    SEXP chainExpr = PROTECT(Rf_allocVector(VECSXP, SLOT_COUNT));
    Rf_setAttrib(chainExpr, R_NamesSymbol, slotNamesExpr);

    storeFlatTrees(chainExpr, SLOT_TREE_VARS, SLOT_TREE_VALUES,
                   SLOT_TREE_SIZES, SLOT_TREE_FLAGS, chainState.trees);
    if (!chainState.treeParams.empty())
      storeTreeParams(chainExpr, SLOT_TREE_PARAMS, chainState.treeParams);
    if (!chainState.treeMasks.empty())
      storeTreeMasks(chainExpr, SLOT_TREE_MASKS, chainState.treeMasks);
    if (!chainState.savedTrees.empty()) {
      storeFlatTrees(chainExpr, SLOT_SAVED_VARS, SLOT_SAVED_VALUES,
                     SLOT_SAVED_SIZES, SLOT_SAVED_FLAGS,
                     chainState.savedTrees);
      if (!chainState.savedTreeParams.empty())
        storeTreeParams(chainExpr, SLOT_SAVED_PARAMS,
                        chainState.savedTreeParams);
      if (!chainState.savedTreeMasks.empty())
        storeTreeMasks(chainExpr, SLOT_SAVED_MASKS,
                       chainState.savedTreeMasks);
    }

    SET_VECTOR_ELT(chainExpr, SLOT_TOTAL_FITS,
                   Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                             chainState.totalFits.size())));
    std::memcpy(REAL(VECTOR_ELT(chainExpr, SLOT_TOTAL_FITS)),
                chainState.totalFits.data(),
                chainState.totalFits.size() * sizeof(double));

    SET_VECTOR_ELT(chainExpr, SLOT_INDICES,
                   Rf_allocVector(INTSXP, static_cast<R_xlen_t>(
                                            chainState.indices.size())));
    int* indices = INTEGER(VECTOR_ELT(chainExpr, SLOT_INDICES));
    for (size_t i = 0; i < chainState.indices.size(); ++i)
      indices[i] = static_cast<int>(chainState.indices[i]);

    SET_VECTOR_ELT(chainExpr, SLOT_SIGMA, Rf_ScalarReal(chainState.sigma));
    SET_VECTOR_ELT(chainExpr, SLOT_K, Rf_ScalarReal(chainState.k));

    // the third element carries the variance prior's internal scale, whose
    // re-derivation through the transform can perturb the last bit
    SET_VECTOR_ELT(chainExpr, SLOT_FIT_SCALE, Rf_allocVector(REALSXP, 3));
    REAL(VECTOR_ELT(chainExpr, SLOT_FIT_SCALE))[0] = chainState.fitMin;
    REAL(VECTOR_ELT(chainExpr, SLOT_FIT_SCALE))[1] = chainState.fitMax;
    REAL(VECTOR_ELT(chainExpr, SLOT_FIT_SCALE))[2] =
      chainState.sigmaPriorScale;

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

  UNPROTECT(4);
  return resultExpr;
}

void setState(bartcore::SamplerBase& sampler, SEXP stateExpr) {
  if (!Rf_inherits(stateExpr, "bartcoreState"))
    Rf_error("'state' must be a bartcore state object");

  SEXP formatVersionExpr = Rf_getAttrib(stateExpr, Rf_install("formatVersion"));
  int formatVersion = Rf_isInteger(formatVersionExpr) &&
      Rf_xlength(formatVersionExpr) == 1 ? INTEGER(formatVersionExpr)[0] : 0;
  if (formatVersion != stateFormatVersion) {
    SEXP packageVersionExpr =
      Rf_getAttrib(stateExpr, Rf_install("packageVersion"));
    const char* packageVersion = Rf_isString(packageVersionExpr) &&
        Rf_xlength(packageVersionExpr) == 1 ?
      CHAR(STRING_ELT(packageVersionExpr, 0)) : "unknown";
    Rf_error("state format version %d (written by dbarts %s) is not "
             "compatible with this dbarts's format version %d; re-fit or "
             "use the dbarts release that wrote it", formatVersion,
             packageVersion, stateFormatVersion);
  }

  if (static_cast<size_t>(Rf_xlength(stateExpr)) != sampler.numChains())
    Rf_error("'state' length must equal number of chains");

  // Rf_error longjmps past destructors, so parse with an error accumulator,
  // free the C++ state, and error at the end.
  const char* errorMessage = NULL;

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

  SEXP sampleNumExpr = Rf_getAttrib(stateExpr, Rf_install("currentSampleNum"));
  if (errorMessage == NULL) {
    if (!Rf_isInteger(sampleNumExpr) || Rf_xlength(sampleNumExpr) != 1 ||
        INTEGER(sampleNumExpr)[0] < 0)
      errorMessage = "malformed sample number in bartcore state";
    else
      state.currentSampleNum = static_cast<size_t>(INTEGER(sampleNumExpr)[0]);
  }

  for (size_t c = 0; c < sampler.numChains() && errorMessage == NULL; ++c) {
    SEXP chainExpr = VECTOR_ELT(stateExpr, static_cast<R_xlen_t>(c));
    bartcore::ChainStateData& chainState(state.chains[c]);

    if (!readFlatTrees(getListElement(chainExpr, "tree.vars"),
                       getListElement(chainExpr, "tree.values"),
                       getListElement(chainExpr, "tree.sizes"),
                       getListElement(chainExpr, "tree.flags"),
                       chainState.trees, &errorMessage))
      break;
    SEXP savedSizesExpr = getListElement(chainExpr, "saved.sizes");
    if (!Rf_isNull(savedSizesExpr) &&
        !readFlatTrees(getListElement(chainExpr, "saved.vars"),
                       getListElement(chainExpr, "saved.values"),
                       savedSizesExpr,
                       getListElement(chainExpr, "saved.flags"),
                       chainState.savedTrees, &errorMessage))
      break;

    // linear-leaf states must carry their slope arrays; function-valued
    // states carry fits slabs and variable-length saved blocks instead
    size_t numLeafCovariates = sampler.numLeafCovariates();
    if (sampler.usesFunctionLeaves()) {
      if (!readFunctionTreeParams(getListElement(chainExpr, "tree.params"),
                                  chainState.trees.size(),
                                  sampler.numObservations(),
                                  chainState.treeParams, &errorMessage))
        break;
      if (!chainState.savedTrees.empty() &&
          !readFunctionSavedParams(getListElement(chainExpr, "saved.params"),
                                   chainState.savedTrees, numLeafCovariates,
                                   chainState.savedTreeParams, &errorMessage))
        break;
    } else if (numLeafCovariates > 0) {
      if (!readTreeParams(getListElement(chainExpr, "tree.params"),
                          chainState.trees, numLeafCovariates,
                          chainState.treeParams, &errorMessage))
        break;
      if (!chainState.savedTrees.empty() &&
          !readTreeParams(getListElement(chainExpr, "saved.params"),
                          chainState.savedTrees, numLeafCovariates,
                          chainState.savedTreeParams, &errorMessage))
        break;
    }

    // wide-categorical states must carry their mask channels
    if (sampler.data().hasWideCategorical) {
      if (!readTreeMasks(getListElement(chainExpr, "tree.masks"),
                         chainState.trees, sampler.data(),
                         chainState.treeMasks, &errorMessage))
        break;
      if (!chainState.savedTrees.empty() &&
          !readTreeMasks(getListElement(chainExpr, "saved.masks"),
                         chainState.savedTrees, sampler.data(),
                         chainState.savedTreeMasks, &errorMessage))
        break;
    }

    SEXP totalFitsExpr = getListElement(chainExpr, "total.fits");
    if (!Rf_isNull(totalFitsExpr)) {
      if (!Rf_isReal(totalFitsExpr)) {
        errorMessage = "malformed fits in bartcore state";
        break;
      }
      chainState.totalFits.assign(
        REAL(totalFitsExpr), REAL(totalFitsExpr) + Rf_xlength(totalFitsExpr));
    }
    SEXP indicesExpr = getListElement(chainExpr, "indices");
    if (!Rf_isNull(indicesExpr)) {
      if (!Rf_isInteger(indicesExpr)) {
        errorMessage = "malformed indices in bartcore state";
        break;
      }
      R_xlen_t numIndices = Rf_xlength(indicesExpr);
      chainState.indices.resize(static_cast<size_t>(numIndices));
      for (R_xlen_t i = 0; i < numIndices; ++i) {
        int index = INTEGER(indicesExpr)[i];
        if (index < 0) {
          errorMessage = "malformed indices in bartcore state";
          break;
        }
        chainState.indices[static_cast<size_t>(i)] =
          static_cast<size_t>(index);
      }
      if (errorMessage != NULL) break;
    }

    SEXP sigmaExpr = getListElement(chainExpr, "sigma");
    SEXP kExpr = getListElement(chainExpr, "k");
    if (!Rf_isReal(sigmaExpr) || Rf_xlength(sigmaExpr) != 1 ||
        !Rf_isReal(kExpr) || Rf_xlength(kExpr) != 1) {
      errorMessage = "malformed parameters in bartcore state";
      break;
    }
    chainState.sigma = REAL(sigmaExpr)[0];
    chainState.k = REAL(kExpr)[0];

    SEXP fitScaleExpr = getListElement(chainExpr, "fit.scale");
    if (!Rf_isReal(fitScaleExpr) || (Rf_xlength(fitScaleExpr) != 2 &&
                                     Rf_xlength(fitScaleExpr) != 3)) {
      errorMessage = "malformed fit scale in bartcore state";
      break;
    }
    chainState.fitMin = REAL(fitScaleExpr)[0];
    chainState.fitMax = REAL(fitScaleExpr)[1];
    // older states lack the prior scale; the restore then re-derives it
    // through the transform, exact only as far as the last ulp
    if (Rf_xlength(fitScaleExpr) == 3)
      chainState.sigmaPriorScale = REAL(fitScaleExpr)[2];

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
  }

  bool restored = errorMessage == NULL && sampler.setState(state);
  {
    bartcore::SamplerStateData empty;
    std::swap(state, empty);  // free before a potential longjmp
  }
  if (errorMessage != NULL) Rf_error("%s", errorMessage);
  if (!restored)
    Rf_error("state is not consistent with this sampler");
}

// Column-major gather of the requested trees, plus the flags that decide
// which data.frame columns exist.
struct GatheredTrees {
  bool includeChain, includeSample, anyWide, anyMissing;
  size_t numSlopes;
  std::vector<int> chain, sample, tree, count, variable, missing;
  std::vector<double> value;
  // wide categorical rules report their masks decoded, one L/R per
  // observed category; the R wrapper pads to the declared level count
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
                 GatheredTrees& out) {
  const bartcore::ColumnStore& store(sampler.data());

  std::vector<bartcore::FlatNode> liveNodes;
  std::vector<double> liveSlopes;
  std::vector<std::uint64_t> liveMasks;
  std::vector<std::uint32_t> counts;
  std::vector<size_t> replayIndices(replayNumRows);
  std::string directionsScratch;
  bool functionLeaves = sampler.usesFunctionLeaves();

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
          nodes = &sampler.savedTree(chainNum, sampleNum, treeNum);
          if (out.numSlopes > 0)
            slopes = &sampler.savedTreeSlopes(chainNum, sampleNum, treeNum);
          if (out.anyWide)
            masks = &sampler.savedTreeMasks(chainNum, sampleNum, treeNum);
        } else {
          sampler.flattenTree(chainNum, treeNum, liveNodes, counts,
                              out.numSlopes > 0 ? &liveSlopes : NULL,
                              out.anyWide ? &liveMasks : NULL);
          nodes = &liveNodes;
          if (out.numSlopes > 0) slopes = &liveSlopes;
          if (out.anyWide) masks = &liveMasks;
        }
        if (replayData != NULL) {
          counts.resize(nodes->size());
          for (size_t l = 0; l < replayNumRows; ++l) replayIndices[l] = l;
          bartcore::countFlatObservationsBelow(
            nodes->data(), store.types.data(), replayData, replayNumRows,
            replayIndices.data(), 0, replayNumRows, counts.data(),
            out.anyWide ? store.numCuts.data() : NULL,
            out.anyWide ? masks->data() : NULL);
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
          bool isWideRule =
            node.variable >= 0 &&
            store.columnHasWideMask(static_cast<size_t>(node.variable));
          // a wide rule's value is a side-channel offset, meaningless
          // outside the format; the directions column carries the decode
          out.value.push_back(
            isWideRule ||
                (functionLeaves && node.variable == bartcore::invalidVariable)
              ? NA_REAL : node.value);
          if (out.anyWide) {
            if (isWideRule) {
              size_t variable = static_cast<size_t>(node.variable);
              const std::uint64_t* words =
                masks->data() + static_cast<size_t>(node.value);
              directionsScratch.clear();
              for (std::uint32_t c = 0; c < store.numCuts[variable]; ++c)
                directionsScratch.push_back(
                  bartcore::maskTestBit(words, c) ? 'R' : 'L');
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
                        (gathered.anyWide ? 1 : 0) +
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
  if (gathered.anyWide) {
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
// marking leaves, split values as data values (an ordinal rule's cut point,
// a categorical rule's mask over observed categories), leaf values on the
// engine's internal response scale. A rule on a wide categorical column
// (more than 53 levels) has no double-exact mask: its value is NA and a
// 'directions' character column - present whenever the store has wide
// columns - carries one L/R per observed category (NA elsewhere; the R
// wrapper pads to the declared level count and fills the narrow rules).
// When any predictor contains missing values a 'missing' integer column
// reports each rule's missing direction (0 left, 1 right; NA on leaves and
// on columns without missing values). Saved trees replay the training
// predictors for n unless newdata is supplied; live trees report their own
// counts. Validates, gathers, and emits, in that order.
SEXP getTrees(bartcore::SamplerBase& sampler, const size_t* chainIndices,
              size_t numChainIndices, const size_t* sampleIndices,
              size_t numSampleIndices, const size_t* treeIndices,
              size_t numTreeIndices, bool useLiveTrees, const double* newdata,
              size_t newdataNumRows, const char* caller) {
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
    if (treeIndices[i] >= sampler.numTrees())
      Rf_error("%s tree number out of range", caller);
  }

  // saved trees carry no counts of their own and replay the training rows;
  // newdata replays its rows through live and saved trees alike
  const double* replayData = NULL;
  size_t replayNumRows = 0;
  if (newdata != NULL) {
    replayData = newdata;
    replayNumRows = newdataNumRows;
  } else if (useSaved) {
    replayData = store.x;
    replayNumRows = store.numObservations;
  }

  bool anyMissing = false;
  for (size_t j = 0; j < store.numPredictors; ++j)
    if (store.hasMissing[j]) { anyMissing = true; break; }

  // function-valued (gp) leaves report no per-leaf coefficients - the
  // function rides prediction only - and their leaves' values print NA
  // (a whole function per row does not fit a data frame)
  size_t numSlopes =
    sampler.usesFunctionLeaves() ? 0 : sampler.numLeafCovariates();

  GatheredTrees gathered{sampler.numChains() > 1, useSaved,
                         store.hasWideCategorical, anyMissing, numSlopes,
                         {}, {}, {}, {}, {}, {}, {}, {}, {}};
  gathered.slopes.resize(numSlopes);
  gatherTrees(sampler, chainIndices, numChainIndices, sampleIndices,
              numSampleIndices, treeIndices, numTreeIndices, useSaved,
              replayData, replayNumRows, gathered);

  return emitTreeDataFrame(gathered);
}

} // namespace bartcore_bridge
