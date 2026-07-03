#include "config.hpp"
#include "R_interface_bartcore.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>

#include <external/Rinternals.h>
#include <R_ext/Random.h> // GetRNGstate, PutRNGstate

#include <external/random.h>
#include <external/stats.h> // ext_percentileOfChiSquared

#include <rc/bounds.h>
#include <rc/util.h>

#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/types.hpp>

#include "R_interface_common.hpp"

#include "bartcore/bartcore.hpp"

using std::size_t;
using namespace dbarts;

namespace {

struct BartcoreHolder {
  std::unique_ptr<bartcore::SamplerBase> sampler;
  std::vector<ext_rng*> rngs;  // one per chain
  bool keepTrainingFits;

  ~BartcoreHolder() {
    for (size_t c = rngs.size(); c > 0; --c)
      if (rngs[c - 1] != NULL) ext_rng_destroy(rngs[c - 1]);
  }
};

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

// errors unless every replacement value for a categorical column is an
// existing category code; ordinal columns pass through
void validateColumnValues(const bartcore::ColumnStore& store, size_t column,
                          const double* values, size_t numValues) {
  if (store.types[column] != bartcore::ColumnType::categorical) return;
  for (size_t i = 0; i < numValues; ++i) {
    if (!store.categoricalValueIsValid(column, values[i]))
      Rf_error("categorical predictor values must be existing category codes");
  }
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

// The classic engine's BARTFit::printInitialSummary, line for line, printed
// at creation under verbose; the classic version's quantile and scale lines
// reduce at creation to expressions over the raw prior scale and the
// response range, used here directly.
void printInitialSummary(const Control& control, const Model& model,
                         const Data& data,
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
  if (model.kPrior->isFixed) {
    ext_printf("\tk prior fixed to %f\n",
               static_cast<FixedHyperprior*>(model.kPrior)->getK());
  } else {
    const ChiHyperprior& kPrior(*static_cast<ChiHyperprior*>(model.kPrior));
    ext_printf("\tprior on k: chi with %f degrees of freedom and %f scale\n",
               kPrior.degreesOfFreedom, kPrior.scale);
  }
  if (!control.responseIsBinary) {
    const ChiSquaredPrior& sigmaSqPrior(
      *static_cast<ChiSquaredPrior*>(model.sigmaSqPrior));
    ext_printf("\tdegrees of freedom in sigma prior: %f\n",
               sigmaSqPrior.degreesOfFreedom);
    double quantile = 1.0 - ext_percentileOfChiSquared(
      sigmaSqPrior.scale * sigmaSqPrior.degreesOfFreedom,
      sigmaSqPrior.degreesOfFreedom);
    ext_printf("\tquantile in sigma prior: %f\n", quantile);
    double sigmaInternal = data.sigmaEstimate / sampler.fitScale();
    ext_printf("\tscale in sigma prior: %f\n",
               sigmaInternal * sigmaInternal * sigmaSqPrior.scale);
  }

  const CGMPrior& treePrior(*static_cast<CGMPrior*>(model.treePrior));
  ext_printf("\tpower and base for tree prior: %f %f\n", treePrior.power,
             treePrior.base);
  if (treePrior.splitProbabilities != NULL) {
    ext_printf("\ttree split probabilities: %f",
               treePrior.splitProbabilities[0]);
    size_t printLength = 5 < data.numPredictors ? 5 : data.numPredictors;
    for (size_t i = 1; i < printLength; ++i)
      ext_printf(", %f", treePrior.splitProbabilities[i]);
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

} // namespace

extern "C" {

// The external pointer's protection slot pins everything the sampler
// borrows: the data expression at creation, and any replacement vectors the
// setters install later.
//
// family selects the response model for binary responses: "" or "probit"
// give the classic probit latents, "logistic" the Polya-Gamma sampler.
// Continuous responses are gaussian and accept only "".
SEXP bartcore_create(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
                     SEXP familyExpr) {
  Control control;
  Data data;
  Model model;

  initializeControlFromExpression(control, controlExpr);
  initializeDataFromExpression(data, dataExpr);
  initializeModelFromExpression(model, modelExpr, control, data);

  const char* familyName =
    Rf_isNull(familyExpr) ? "" : CHAR(STRING_ELT(familyExpr, 0));

  // Rf_error longjmps past destructors, so collect the reason, clean up,
  // and error at the end.
  const char* errorMessage = NULL;

  bartcore::ResponseFamily family = bartcore::ResponseFamily::gaussian;
  if (control.responseIsBinary) {
    if (std::strcmp(familyName, "logistic") == 0) {
      family = bartcore::ResponseFamily::logistic;
    } else if (familyName[0] == '\0' ||
               std::strcmp(familyName, "probit") == 0) {
      family = bartcore::ResponseFamily::probit;
    } else {
      errorMessage = "unrecognized response family for a binary response";
    }
  } else if (familyName[0] != '\0' &&
             std::strcmp(familyName, "gaussian") != 0) {
    errorMessage = "response families other than gaussian require a binary "
                   "response";
  }

  std::vector<bartcore::ColumnType> columnTypes(
    data.numPredictors, bartcore::ColumnType::ordinal);
  bool anyCategorical = false;
  for (size_t j = 0; j < data.numPredictors && errorMessage == NULL; ++j) {
    if (data.variableTypes[j] != CATEGORICAL) continue;
    columnTypes[j] = bartcore::ColumnType::categorical;
    anyCategorical = true;
    for (size_t i = 0; i < data.numObservations && errorMessage == NULL; ++i) {
      double value = data.x[i + j * data.numObservations];
      if (value < 0.0 ||
          value >= static_cast<double>(bartcore::maxCategories) ||
          value != std::floor(value))
        errorMessage = "categorical predictors must hold integer category "
                       "codes in [0, 53)";
    }
  }
  if (errorMessage == NULL && anyCategorical &&
      data.numTestObservations > 0) {
    // test codes must also be representable; category counts come from the
    // training columns
    for (size_t j = 0; j < data.numPredictors && errorMessage == NULL; ++j) {
      if (columnTypes[j] != bartcore::ColumnType::categorical) continue;
      double maxValue = 0.0;
      for (size_t i = 0; i < data.numObservations; ++i) {
        double value = data.x[i + j * data.numObservations];
        if (value > maxValue) maxValue = value;
      }
      for (size_t i = 0;
           i < data.numTestObservations && errorMessage == NULL; ++i) {
        double value = data.x_test[i + j * data.numTestObservations];
        if (value < 0.0 || value > maxValue || value != std::floor(value))
          errorMessage = "categorical test predictors must hold existing "
                         "category codes";
      }
    }
  }
  if (errorMessage == NULL && !control.responseIsBinary &&
      model.sigmaSqPrior->isFixed)
    errorMessage = "bartcore does not support fixed sigma for continuous responses";
  if (errorMessage == NULL && control.responseIsBinary &&
      data.weights != NULL)
    errorMessage = "binary response families do not support weights: the "
                   "weighted probit the classic engine fit was incorrect; "
                   "replicate rows or model the latents instead";

  double k = 2.0, sigmaDf = 3.0, sigmaRawScale = 1.0;
  bool updateK = !model.kPrior->isFixed;
  double kDf = 1.25, kScale = HUGE_VAL;
  if (errorMessage == NULL) {
    if (updateK) {
      const ChiHyperprior& kPrior(*static_cast<ChiHyperprior*>(model.kPrior));
      kDf = kPrior.degreesOfFreedom;
      kScale = kPrior.scale;
    } else {
      k = static_cast<FixedHyperprior*>(model.kPrior)->getK();
    }
    if (!control.responseIsBinary) {
      const ChiSquaredPrior& sigmaSqPrior(
        *static_cast<ChiSquaredPrior*>(model.sigmaSqPrior));
      sigmaDf = sigmaSqPrior.degreesOfFreedom;
      sigmaRawScale = sigmaSqPrior.scale;
    }
  }

  if (errorMessage != NULL) {
    invalidateModel(model);
    invalidateData(data);
    Rf_error("%s", errorMessage);
  }

  const CGMPrior& treePrior(*static_cast<CGMPrior*>(model.treePrior));

  bartcore::SamplerOptions options;
  options.numTrees = control.numTrees;
  options.k = k;
  options.nodeScale = model.nodeScale;
  options.base = treePrior.base;
  options.power = treePrior.power;
  options.birthOrDeathProbability = model.birthOrDeathProbability;
  options.swapProbability = model.swapProbability;
  options.changeProbability = model.changeProbability;
  options.birthProbability = model.birthProbability;
  options.maxNumCutsPerVariable = data.maxNumCuts; // copied during build
  options.useQuantiles = control.useQuantiles;
  options.columnTypes = anyCategorical ? columnTypes.data() : NULL;
  options.splitProbabilities = treePrior.splitProbabilities; // copied by ctor

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
  options.updateK = updateK;
  options.kHyperprior.degreesOfFreedom = kDf;
  options.kHyperprior.scale = kScale;
  options.numChains = control.numChains;
  options.numThreads = control.numThreads;
  options.numThin = control.treeThinningRate;
  options.keepTrees = control.keepTrees;
  options.numSamplesToStore = control.defaultNumSamples;
  options.verbose = control.verbose;
  options.printEvery = control.printEvery;

  // A single chain draws through R's generator; several chains each get a
  // Mersenne twister seeded from R's stream, so results do not depend on the
  // thread count and worker threads never touch the R API.
  std::vector<ext_rng*> rngs(options.numChains, static_cast<ext_rng*>(NULL));
  bool rngFailed = false;
  if (options.numChains == 1) {
    rngs[0] = ext_rng_createDefault(true);
    rngFailed = rngs[0] == NULL;
  } else {
    GetRNGstate();
    for (size_t c = 0; c < options.numChains && !rngFailed; ++c) {
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
    invalidateModel(model);
    invalidateData(data);
    Rf_error("could not allocate rng");
  }

  std::unique_ptr<bartcore::SamplerBase> sampler = bartcore::createClassicSampler(
    data.x, data.y, data.numObservations, data.numPredictors, data.weights,
    data.offset, family, data.sigmaEstimate, sigmaDf, sigmaRawScale, options,
    rngs.data());

  if (data.numTestObservations > 0) {
    sampler->setTestPredictors(data.x_test, data.numTestObservations);
    sampler->setTestOffset(data.testOffset);
  }

  if (control.verbose) printInitialSummary(control, model, data, *sampler);

  invalidateModel(model);
  invalidateData(data);

  BartcoreHolder* holder = new BartcoreHolder{std::move(sampler),
                                              std::move(rngs),
                                              control.keepTrainingFits};

  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(holder, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  UNPROTECT(2);
  return result;
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
    sampler.run(numBurnIn, 0, empty);
    PutRNGstate();
    return R_NilValue;
  }

  // several chains add a trailing chain dimension, as the classic engine's
  // results do
  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP sigmaExpr = numChains == 1
    ? PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples)))
    : PROTECT(Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  SEXP trainExpr;
  if (!holder.keepTrainingFits) {
    trainExpr = PROTECT(R_NilValue);
  } else if (numChains == 1) {
    trainExpr = PROTECT(Rf_allocMatrix(
      REALSXP, static_cast<int>(numObservations), numSamplesInt));
  } else {
    trainExpr = PROTECT(Rf_alloc3DArray(
      REALSXP, static_cast<int>(numObservations), numSamplesInt,
      numChainsInt));
  }
  SEXP testExpr;
  if (numTestObservations > 0) {
    testExpr = numChains == 1
      ? PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numTestObservations),
                               numSamplesInt))
      : PROTECT(Rf_alloc3DArray(REALSXP,
                                static_cast<int>(numTestObservations),
                                numSamplesInt, numChainsInt));
  } else {
    testExpr = PROTECT(R_NilValue);
  }
  SEXP varcountExpr = numChains == 1
    ? PROTECT(Rf_allocMatrix(INTSXP, static_cast<int>(numPredictors),
                             numSamplesInt))
    : PROTECT(Rf_alloc3DArray(INTSXP, static_cast<int>(numPredictors),
                              numSamplesInt, numChainsInt));
  SEXP kExpr;
  if (sampler.kIsSampled()) {
    kExpr = numChains == 1
      ? PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples)))
      : PROTECT(Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  } else {
    kExpr = PROTECT(R_NilValue);
  }

  std::vector<std::uint32_t> variableCounts(numPredictors * numSamples *
                                            numChains);

  bartcore::Results results;
  results.sigma = REAL(sigmaExpr);
  results.trainingFits = holder.keepTrainingFits ? REAL(trainExpr) : NULL;
  results.testFits = numTestObservations > 0 ? REAL(testExpr) : NULL;
  results.variableCounts = variableCounts.data();
  results.k = sampler.kIsSampled() ? REAL(kExpr) : NULL;

  GetRNGstate();
  sampler.run(numBurnIn, numSamples, results);
  PutRNGstate();

  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < numPredictors * numSamples * numChains; ++i)
    varcountOut[i] = static_cast<int>(variableCounts[i]);

  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, trainExpr);
  SET_VECTOR_ELT(resultExpr, 2, testExpr);
  SET_VECTOR_ELT(resultExpr, 3, varcountExpr);
  SET_VECTOR_ELT(resultExpr, 4, kExpr);

  // named as the classic engine's run results are, so the engines are
  // drop-in replacements for each other
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(7);
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
  const double* offset = Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr);
  holder.sampler->setOffset(offset, Rf_asLogical(updateScaleExpr) == TRUE);
  retain(ptrExpr, PROT_OFFSET, offsetExpr);
  return R_NilValue;
}

SEXP bartcore_setResponse(SEXP ptrExpr, SEXP yExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  GetRNGstate(); // probit latent redraw
  holder.sampler->setResponse(REAL(yExpr));
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

  if (!Rf_inherits(dataExpr, "dbartsData"))
    Rf_error("'data' argument to bartcore_setData not of class 'dbartsData'");

  Data data;
  initializeDataFromExpression(data, dataExpr);

  // Rf_error longjmps past destructors, so collect the reason, clean up,
  // and error at the end.
  const char* errorMessage = NULL;
  if (data.numPredictors != sampler.numPredictors())
    errorMessage = "bartcore setData requires the same predictors";
  if (errorMessage == NULL &&
      sampler.family() != bartcore::ResponseFamily::gaussian &&
      data.weights != NULL)
    errorMessage = "binary response families do not support weights: the "
                   "weighted probit the classic engine fit was incorrect; "
                   "replicate rows or model the latents instead";
  for (size_t j = 0; j < data.numPredictors && errorMessage == NULL; ++j) {
    bool wasCategorical = sampler.data().types[j] ==
                          bartcore::ColumnType::categorical;
    if ((data.variableTypes[j] == CATEGORICAL) != wasCategorical) {
      errorMessage = "bartcore setData requires the same predictor types";
      break;
    }
    if (!wasCategorical) continue;
    // category counts are fixed at creation; new values must be existing
    // codes, in the training and test data both
    for (size_t i = 0; i < data.numObservations && errorMessage == NULL; ++i)
      if (!sampler.data().categoricalValueIsValid(
            j, data.x[i + j * data.numObservations]))
        errorMessage = "categorical predictor values must be existing "
                       "category codes";
    for (size_t i = 0;
         i < data.numTestObservations && errorMessage == NULL; ++i)
      if (!sampler.data().categoricalValueIsValid(
            j, data.x_test[i + j * data.numTestObservations]))
        errorMessage = "categorical predictor values must be existing "
                       "category codes";
  }

  if (errorMessage != NULL) {
    invalidateData(data);
    Rf_error("%s", errorMessage);
  }

  sampler.setData(data.x, data.y, data.numObservations, data.weights,
                  data.offset, data.x_test, data.numTestObservations,
                  data.testOffset);

  invalidateData(data);

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
// clears the offset alongside the new predictors.
SEXP bartcore_setTestPredictorAndOffset(SEXP ptrExpr, SEXP xTestExpr,
                                        SEXP offsetExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
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
  if (holder.sampler->family() != bartcore::ResponseFamily::gaussian)
    Rf_error("binary response families do not support weights: the "
             "weighted probit the classic engine fit was incorrect; "
             "replicate rows or model the latents instead");
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

  Control control;
  initializeControlFromExpression(control, controlExpr);

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

  Control control;
  Data data;
  Model model;
  initializeControlFromExpression(control, controlExpr);
  initializeDataFromExpression(data, dataExpr);
  initializeModelFromExpression(model, modelExpr, control, data);

  bool isGaussian = sampler.family() == bartcore::ResponseFamily::gaussian;
  const char* errorMessage = NULL;
  if (isGaussian && model.sigmaSqPrior->isFixed)
    errorMessage =
      "bartcore does not support fixed sigma for continuous responses";

  if (errorMessage == NULL) {
    const CGMPrior& treePrior(*static_cast<CGMPrior*>(model.treePrior));

    bartcore::ModelParameters parameters;
    parameters.base = treePrior.base;
    parameters.power = treePrior.power;
    parameters.splitProbabilities = treePrior.splitProbabilities;
    parameters.birthOrDeathProbability = model.birthOrDeathProbability;
    parameters.swapProbability = model.swapProbability;
    parameters.changeProbability = model.changeProbability;
    parameters.birthProbability = model.birthProbability;
    parameters.nodeScale = model.nodeScale;
    parameters.updateK = !model.kPrior->isFixed;
    if (parameters.updateK) {
      const ChiHyperprior& kPrior(*static_cast<ChiHyperprior*>(model.kPrior));
      parameters.kHyperprior.degreesOfFreedom = kPrior.degreesOfFreedom;
      parameters.kHyperprior.scale = kPrior.scale;
    } else {
      parameters.k = static_cast<FixedHyperprior*>(model.kPrior)->getK();
    }
    if (isGaussian) {
      const ChiSquaredPrior& sigmaSqPrior(
        *static_cast<ChiSquaredPrior*>(model.sigmaSqPrior));
      parameters.sigmaEstimate = data.sigmaEstimate;
      parameters.sigmaDf = sigmaSqPrior.degreesOfFreedom;
      parameters.sigmaRawScale = sigmaSqPrior.scale;
    }

    // split probabilities are copied per chain before the model goes away
    sampler.setModel(parameters);
  }

  invalidateModel(model);
  invalidateData(data);
  if (errorMessage != NULL) Rf_error("%s", errorMessage);

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

SEXP bartcore_setCutPoints(SEXP ptrExpr, SEXP cutPointsExpr, SEXP columnsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
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
// fits, per-tree observation orderings, internal-scale sigma, k, latents,
// dart state, and the serialized rng), with the store's cut points and the
// saved-tree write position as attributes. A sampler restored from it over
// the same data continues bitwise identically to one that was never stored.

// flattened trees as three parallel R vectors: concatenated 1-based-with-(-1)
// variables, values, and per-tree node counts
static void storeFlatTrees(SEXP chainExpr, int variablesSlot, int valuesSlot,
                           int sizesSlot,
                           const std::vector<std::vector<bartcore::FlatNode>>& trees) {
  R_xlen_t totalNumNodes = 0;
  for (const std::vector<bartcore::FlatNode>& tree : trees)
    totalNumNodes += static_cast<R_xlen_t>(tree.size());

  SET_VECTOR_ELT(chainExpr, variablesSlot,
                 Rf_allocVector(INTSXP, totalNumNodes));
  SET_VECTOR_ELT(chainExpr, valuesSlot, Rf_allocVector(REALSXP, totalNumNodes));
  SET_VECTOR_ELT(chainExpr, sizesSlot,
                 Rf_allocVector(INTSXP, static_cast<R_xlen_t>(trees.size())));

  int* variables = INTEGER(VECTOR_ELT(chainExpr, variablesSlot));
  double* values = REAL(VECTOR_ELT(chainExpr, valuesSlot));
  int* sizes = INTEGER(VECTOR_ELT(chainExpr, sizesSlot));
  R_xlen_t offset = 0;
  for (size_t t = 0; t < trees.size(); ++t) {
    sizes[t] = static_cast<int>(trees[t].size());
    for (const bartcore::FlatNode& node : trees[t]) {
      variables[offset] = node.variable >= 0 ? node.variable + 1
                                             : node.variable;
      values[offset] = node.value;
      ++offset;
    }
  }
}

// the inverse; errorMessage is set instead of erroring so callers can clean
// up C++ state first
static bool readFlatTrees(SEXP variablesExpr, SEXP valuesExpr, SEXP sizesExpr,
                          std::vector<std::vector<bartcore::FlatNode>>& trees,
                          const char** errorMessage) {
  if (!Rf_isInteger(variablesExpr) || !Rf_isReal(valuesExpr) ||
      !Rf_isInteger(sizesExpr) ||
      Rf_xlength(variablesExpr) != Rf_xlength(valuesExpr)) {
    *errorMessage = "malformed trees in bartcore state";
    return false;
  }
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
      ++offset;
    }
  }
  return true;
}

SEXP bartcore_storeState(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  bartcore::SamplerStateData state;
  sampler.getState(state);

  size_t numChains = state.chains.size();
  size_t numObservations = sampler.numObservations();

  enum {
    SLOT_TREE_VARS = 0, SLOT_TREE_VALUES, SLOT_TREE_SIZES, SLOT_SAVED_VARS,
    SLOT_SAVED_VALUES, SLOT_SAVED_SIZES, SLOT_TOTAL_FITS, SLOT_INDICES,
    SLOT_SIGMA, SLOT_K, SLOT_LATENTS, SLOT_DART_PROBABILITIES,
    SLOT_DART_ALPHA, SLOT_DART_UPDATES_SKIPPED, SLOT_RNG_STATE, SLOT_COUNT
  };
  static const char* slotNames[SLOT_COUNT] = {
    "tree.vars", "tree.values", "tree.sizes", "saved.vars", "saved.values",
    "saved.sizes", "total.fits", "indices", "sigma", "k", "latents",
    "dart.probabilities", "dart.alpha", "dart.updates.skipped", "rng.state"
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
                   SLOT_TREE_SIZES, chainState.trees);
    if (!chainState.savedTrees.empty())
      storeFlatTrees(chainExpr, SLOT_SAVED_VARS, SLOT_SAVED_VALUES,
                     SLOT_SAVED_SIZES, chainState.savedTrees);

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

    if (!chainState.latents.empty()) {
      SET_VECTOR_ELT(chainExpr, SLOT_LATENTS,
                     Rf_allocVector(REALSXP, static_cast<R_xlen_t>(
                                               numObservations)));
      std::memcpy(REAL(VECTOR_ELT(chainExpr, SLOT_LATENTS)),
                  chainState.latents.data(),
                  numObservations * sizeof(double));
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
  SEXP classExpr = PROTECT(Rf_mkString("bartcoreState"));
  Rf_setAttrib(resultExpr, R_ClassSymbol, classExpr);

  UNPROTECT(4);
  return resultExpr;
}

SEXP bartcore_setState(SEXP ptrExpr, SEXP stateExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  if (!Rf_inherits(stateExpr, "bartcoreState"))
    Rf_error("'state' must be a bartcore state object");
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
                       chainState.trees, &errorMessage))
      break;
    SEXP savedSizesExpr = getListElement(chainExpr, "saved.sizes");
    if (!Rf_isNull(savedSizesExpr) &&
        !readFlatTrees(getListElement(chainExpr, "saved.vars"),
                       getListElement(chainExpr, "saved.values"),
                       savedSizesExpr, chainState.savedTrees, &errorMessage))
      break;

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

    SEXP latentsExpr = getListElement(chainExpr, "latents");
    if (!Rf_isNull(latentsExpr)) {
      if (!Rf_isReal(latentsExpr)) {
        errorMessage = "malformed latents in bartcore state";
        break;
      }
      chainState.latents.assign(
        REAL(latentsExpr), REAL(latentsExpr) + Rf_xlength(latentsExpr));
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

// A data.frame of tree structure in the classic engine's format: pre-order
// rows of ([chain,] [sample,] tree, n, var, value), var 1-based with -1
// marking leaves, split values as data values (an ordinal rule's cut point,
// a categorical rule's direction mask), leaf values on the engine's internal
// response scale. Saved trees replay the training predictors for n unless
// newdata is supplied; live trees report their own counts.
SEXP bartcore_getTrees(SEXP ptrExpr, SEXP chainNumsExpr, SEXP sampleNumsExpr,
                       SEXP treeNumsExpr, SEXP currentExpr, SEXP newdataExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);
  const bartcore::ColumnStore& store(sampler.data());

  bool useLiveTrees = Rf_asLogical(currentExpr) == TRUE;
  bool useSaved = sampler.savedTreeCapacity() > 0 && !useLiveTrees;

  size_t numChainIndices = static_cast<size_t>(Rf_xlength(chainNumsExpr));
  size_t numTreeIndices = static_cast<size_t>(Rf_xlength(treeNumsExpr));
  size_t numSampleIndices = useSaved
    ? static_cast<size_t>(Rf_xlength(sampleNumsExpr)) : 1;
  for (size_t i = 0; i < numChainIndices; ++i) {
    int chainNum = INTEGER(chainNumsExpr)[i];
    if (chainNum < 1 || static_cast<size_t>(chainNum) > sampler.numChains())
      Rf_error("bartcore_getTrees chain number out of range");
  }
  if (useSaved) {
    for (size_t i = 0; i < numSampleIndices; ++i) {
      int sampleNum = INTEGER(sampleNumsExpr)[i];
      if (sampleNum < 1 ||
          static_cast<size_t>(sampleNum) > sampler.savedTreeCapacity())
        Rf_error("bartcore_getTrees sample number out of range");
    }
  }

  const double* replayData = store.x;
  size_t replayNumRows = store.numObservations;
  bool replay = useSaved;
  if (!Rf_isNull(newdataExpr)) {
    replayNumRows =
      validatePredictorMatrix(sampler, newdataExpr, "bartcore_getTrees");
    replayData = REAL(newdataExpr);
    replay = true;
  }

  std::vector<int> chainColumn, sampleColumn, treeColumn, countColumn,
    variableColumn;
  std::vector<double> valueColumn;
  std::vector<bartcore::FlatNode> liveNodes;
  std::vector<std::uint32_t> counts;
  std::vector<size_t> replayIndices(replayNumRows);

  for (size_t i = 0; i < numChainIndices; ++i) {
    size_t chainNum = static_cast<size_t>(INTEGER(chainNumsExpr)[i] - 1);
    for (size_t j = 0; j < numSampleIndices; ++j) {
      size_t sampleNum = useSaved
        ? static_cast<size_t>(INTEGER(sampleNumsExpr)[j] - 1) : 0;
      for (size_t k = 0; k < numTreeIndices; ++k) {
        int treeNum = INTEGER(treeNumsExpr)[k];
        if (treeNum < 1 || static_cast<size_t>(treeNum) > sampler.numTrees())
          Rf_error("bartcore_getTrees tree number out of range");

        const std::vector<bartcore::FlatNode>* nodes;
        if (useSaved) {
          nodes = &sampler.savedTree(chainNum, sampleNum,
                                     static_cast<size_t>(treeNum - 1));
        } else {
          sampler.flattenTree(chainNum, static_cast<size_t>(treeNum - 1),
                              liveNodes, counts);
          nodes = &liveNodes;
        }
        if (replay) {
          counts.resize(nodes->size());
          for (size_t l = 0; l < replayNumRows; ++l) replayIndices[l] = l;
          bartcore::countFlatObservationsBelow(
            nodes->data(), store.types.data(), replayData, replayNumRows,
            replayIndices.data(), 0, replayNumRows, counts.data());
        }

        for (size_t l = 0; l < nodes->size(); ++l) {
          chainColumn.push_back(static_cast<int>(chainNum + 1));
          sampleColumn.push_back(static_cast<int>(sampleNum + 1));
          treeColumn.push_back(treeNum);
          countColumn.push_back(static_cast<int>(counts[l]));
          const bartcore::FlatNode& node((*nodes)[l]);
          variableColumn.push_back(
            node.variable >= 0 ? node.variable + 1 : node.variable);
          valueColumn.push_back(node.value);
        }
      }
    }
  }

  R_xlen_t totalNumNodes = static_cast<R_xlen_t>(valueColumn.size());
  R_xlen_t numColumns = 4 + (sampler.numChains() > 1 ? 1 : 0) +
                        (useSaved ? 1 : 0);

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, numColumns));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, numColumns));

  R_xlen_t columnNum = 0;
  if (sampler.numChains() > 1) {
    SET_VECTOR_ELT(resultExpr, columnNum,
                   Rf_allocVector(INTSXP, totalNumNodes));
    SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar("chain"));
    std::memcpy(INTEGER(VECTOR_ELT(resultExpr, columnNum)),
                chainColumn.data(), chainColumn.size() * sizeof(int));
    ++columnNum;
  }
  if (useSaved) {
    SET_VECTOR_ELT(resultExpr, columnNum,
                   Rf_allocVector(INTSXP, totalNumNodes));
    SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar("sample"));
    std::memcpy(INTEGER(VECTOR_ELT(resultExpr, columnNum)),
                sampleColumn.data(), sampleColumn.size() * sizeof(int));
    ++columnNum;
  }
  SET_VECTOR_ELT(resultExpr, columnNum, Rf_allocVector(INTSXP, totalNumNodes));
  SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar("tree"));
  std::memcpy(INTEGER(VECTOR_ELT(resultExpr, columnNum)), treeColumn.data(),
              treeColumn.size() * sizeof(int));
  ++columnNum;
  SET_VECTOR_ELT(resultExpr, columnNum, Rf_allocVector(INTSXP, totalNumNodes));
  SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar("n"));
  std::memcpy(INTEGER(VECTOR_ELT(resultExpr, columnNum)), countColumn.data(),
              countColumn.size() * sizeof(int));
  ++columnNum;
  SET_VECTOR_ELT(resultExpr, columnNum, Rf_allocVector(INTSXP, totalNumNodes));
  SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar("var"));
  std::memcpy(INTEGER(VECTOR_ELT(resultExpr, columnNum)),
              variableColumn.data(), variableColumn.size() * sizeof(int));
  ++columnNum;
  SET_VECTOR_ELT(resultExpr, columnNum, Rf_allocVector(REALSXP, totalNumNodes));
  SET_STRING_ELT(namesExpr, columnNum, Rf_mkChar("value"));
  std::memcpy(REAL(VECTOR_ELT(resultExpr, columnNum)), valueColumn.data(),
              valueColumn.size() * sizeof(double));

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

SEXP bartcore_getLatents(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (holder.sampler->latents(0) == NULL) return R_NilValue;

  size_t numObservations = holder.sampler->numObservations();
  size_t numChains = holder.sampler->numChains();

  SEXP result;
  if (numChains == 1) {
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
