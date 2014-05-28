#include "config.hpp"
#include <cstring>
#include <cstddef>
#include <stdint.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include <bart/bartFit.hpp>
#include <bart/control.hpp>
#include <bart/data.hpp>
#include <bart/model.hpp>
#include <bart/results.hpp>
#include <bart/R_C_interface.hpp>

using std::size_t;

extern "C" {
  SEXP cbart_fitBart(SEXP y, SEXP x, SEXP x_test, SEXP sigmaEstimate,
                     SEXP sigmaDf, SEXP sigmaQuantile, SEXP kFactor,
                     SEXP treeGrowthPower, SEXP treeGrowthBase, SEXP responseIsBinary,
                     SEXP binaryOffset, SEXP numTrees, SEXP numSamples,
                     SEXP numBurnIn, SEXP numThreads, SEXP printEvery, SEXP keepTrainingFits,
                     SEXP useQuantiles, SEXP maxNumCuts, SEXP printCutoffs, SEXP verbose)
  {
    bart::Control control;
    bart::Data data;
    
    int* dims;
    int i_temp; double d_temp;
    
    if (!isReal(y)) error("y must be of type real.");
    if (length(y) == 0) error("Length of y must be greater than 0.");
    data.y = REAL(y);
    data.numObservations = length(y);
    
    if (!isReal(x)) error("x must be of type real.");
    dims = INTEGER(getAttrib(x, R_DimSymbol));
    if (dims == NULL || length(getAttrib(x, R_DimSymbol)) < 2) error("x must be a matrix.");
    if (dims[0] != (int) data.numObservations) error("Number of rows of x and length of y must be equal.");
    data.X = REAL(x);
    data.numPredictors = dims[1];
    
    if (isNull(x_test) || length(x_test) == 0) {
      data.X_test = NULL;
      data.numTestObservations = 0;
    } else {
      if (!isReal(x_test)) error ("x_test must be of type real.");
      dims = INTEGER(getAttrib(x_test, R_DimSymbol));
      if (dims == NULL || length(getAttrib(x_test, R_DimSymbol)) < 2) error("x_test must be a matrix.");
      if (dims[1] != (int) data.numPredictors) error("NUmber of columns of x_test and x must be equal.");
      data.X_test = REAL(x_test);
      data.numTestObservations = dims[0];
    }
    
    if (!isReal(sigmaEstimate)) error("sigma estimate must be of type real.");
    if (length(sigmaEstimate) == 0) error("sigma estimate must be of length at least 1.");
    d_temp = REAL(sigmaEstimate)[0];
    if (ISNAN(d_temp)) d_temp = 1.0;
    if (d_temp <= 0.0) error("sigma estimate must be positive.");
    data.sigmaEstimate = d_temp;
    
    if (!isInteger(sigmaDf)) error("sigma degrees of freedom must be of integer type.");
    if (length(sigmaDf) == 0) error("sigma degrees of freedom must be of length at least 1.");
    i_temp = INTEGER(sigmaDf)[0];
    if (i_temp == NA_INTEGER) error("sigma degrees of freedom cannot be NA.");
    if (i_temp <= 0) error("sigma degrees of freedom must be positive.");
    control.sigmaDf = (uint32_t) i_temp;
    
    if (!isReal(sigmaQuantile)) error ("sigma quantile must be of type real.");
    if (length(sigmaQuantile) == 0) error("sigma quantile must be of length at least 1.");
    d_temp = REAL(sigmaQuantile)[0];
    if (ISNAN(d_temp)) error("sigma quantile must be a real number.");
    if (d_temp <= 0.0) error("sigma quantile must be positive.");
    control.sigmaQuantile = d_temp;
    
    if (!isReal(kFactor)) error ("k factor must be of type real.");
    if (length(kFactor) == 0) error("k factor must be of length at least 1.");
    d_temp = REAL(kFactor)[0];
    if (ISNAN(d_temp)) error("k factor must be a real number.");
    if (d_temp <= 0.0) error("k factor must be positive.");
    control.kFactor = d_temp;
    
    if (!isReal(treeGrowthPower)) error("Tree growth power must be of type real.");
    if (length(treeGrowthPower) == 0) error("Tree growth power must be of length at least 1.");
    d_temp = REAL(treeGrowthPower)[0];
    if (ISNAN(d_temp)) error("Tree growth power be a real number.");
    if (d_temp <= 0.0) error("Tree growth power must be positive.");
    control.power = d_temp;
    
    if (!isReal(treeGrowthBase)) error("Tree growth base must be of type real.");
    if (length(treeGrowthBase) == 0) error("Tree growth power must be of length at least 1.");
    d_temp = REAL(treeGrowthBase)[0];
    if (ISNAN(d_temp)) error("Tree growth base be a real number.");
    if (d_temp <= 0.0) error("Tree growth base must be positive.");
    control.base = d_temp;
    
    if (!isLogical(responseIsBinary)) error("Binary response must be signified by logical type.");
    if (length(responseIsBinary) == 0) error("Binary response signifier must be of length at least 1.");
    i_temp = LOGICAL(responseIsBinary)[0];
    if (i_temp == NA_LOGICAL) error("Binary response must be either true or false.");
    control.responseIsBinary = (i_temp == TRUE);
    
    if (control.responseIsBinary) {
      if (!isReal(binaryOffset)) error("Binary offset must be of type real.");
      if (length(binaryOffset) == 0) error("Binary offset must be of length at least 1.");
      d_temp = REAL(binaryOffset)[0];
      if (ISNAN(d_temp)) error("Binary offset must be a real number.");
      control.binaryOffset = d_temp;
    } else {
      control.binaryOffset = 0.0;
    }
    
    if (!isInteger(numTrees)) error("Number of trees must be of integer type.");
    if (length(numTrees) == 0) error("Number of trees must be of length at least 1.");
    i_temp = INTEGER(numTrees)[0];
    if (i_temp == NA_INTEGER) error("Number of trees cannot be NA.");
    if (i_temp <= 0) error("Number of trees must be positive.");
    control.numTrees = (uint32_t) i_temp;
    
    if (!isInteger(numSamples)) error("Number of samples must be of integer type.");
    if (length(numSamples) == 0) error("Number of samples must be of length at least 1.");
    i_temp = INTEGER(numSamples)[0];
    if (i_temp == NA_INTEGER) error("Number of samples cannot be NA.");
    if (i_temp <= 0) error("Number of samples must be positive.");
    control.numSamples = (uint32_t) i_temp;
    
    if (!isInteger(numBurnIn)) error("Number of burn-in steps must be of integer type.");
    if (length(numBurnIn) == 0) error("Number of burn-in steps must be of length at least 1.");
    i_temp = INTEGER(numBurnIn)[0];
    if (i_temp == NA_INTEGER) i_temp = 0;
    if (i_temp < 0) error("Number of burn-in steps must be non-negative.");
    control.numBurnIn = (size_t) i_temp;
    
    if (!isInteger(numThreads)) error("Number of threads must be of integer type.");
    if (length(numThreads) == 0) error("Number of threads must be of length at least 1.");
    i_temp = INTEGER(numThreads)[0];
    if (i_temp == NA_INTEGER) i_temp = 1;
    if (i_temp <= 0) error("Number of threads must be positive.");
    control.numThreads = (size_t) i_temp;
    
    if (!isInteger(printEvery)) error("Print every must be of integer type.");
    if (length(printEvery) == 0) error("Print every must be of length at least 1.");
    i_temp = INTEGER(printEvery)[0];
    if (i_temp != NA_INTEGER) {
      if (i_temp <= 0) error("Print every must be positive.");
      control.printEvery = (uint32_t) i_temp;
    }
    
    if (!isLogical(keepTrainingFits)) error("Keep training fits must be signified by logical type.");
    if (length(keepTrainingFits) == 0) error("Keep training fits must be of length at least 1.");
    i_temp = LOGICAL(keepTrainingFits)[0];
    if (i_temp == NA_LOGICAL) error("Keep training fits must be either true or false.");
    control.keepTrainingFits = (i_temp == TRUE);
    
    if (!isLogical(useQuantiles)) error("Use quantiles must be signified by logical type.");
    if (length(useQuantiles) == 0) error("Use quantiles must be of length at least 1.");
    i_temp = LOGICAL(useQuantiles)[0];
    if (i_temp == NA_LOGICAL) error("Use quantiles must be either true or false.");
    control.useQuantiles = (i_temp == TRUE);
    
    if (!isInteger(maxNumCuts)) error("Maximum number of cuts must be of integer type.");
    if (length(maxNumCuts) != (int) data.numPredictors) error("Length of maximum number of cuts and the number of columns of x must be equal.");
    int* i_maxNumCuts = INTEGER(maxNumCuts);
    uint32_t* u_maxNumCuts = new uint32_t[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) u_maxNumCuts[i] = (uint32_t) i_maxNumCuts[i];
    data.maxNumCuts = u_maxNumCuts;
    
    if (!isInteger(printCutoffs)) error("Print cutoffs must be of integer type.");
    if (length(printCutoffs) == 0) error("Print cutoffs must be of length at least 1.");
    i_temp = INTEGER(printCutoffs)[0];
    if (i_temp == NA_INTEGER) i_temp = 0;
    if (i_temp < 0) error("Print cutoffs must be non-negative.");
    control.printCutoffs = (uint32_t) i_temp;
    
    if (!isLogical(verbose)) error("Verbose must be signified by logical type.");
    if (length(verbose) == 0) error("Verbose must be of length at least 1.");
    i_temp = LOGICAL(verbose)[0];
    if (i_temp == NA_LOGICAL) error("Verbose must be either true or false.");
    control.verbose = (i_temp == TRUE);
    
    
    bart::Model model;
    bart::CGMPrior cgmPrior(control);
    bart::NormalPrior muPrior(control);
    bart::ChiSquaredPrior sigmaSqPrior(control);

    model.treePrior = &cgmPrior;
    model.muPrior = &muPrior;
    model.sigmaSqPrior = &sigmaSqPrior;

    
    GetRNGstate();
    
    bart::BARTFit fit(control, model, data);
    bart::Results* bartResults = fit.runSampler();
    
    PutRNGstate();
    
    
    // create result storage and make it user friendly
    SEXP dimsExpr, namesExpr;
    
    SEXP resultExpr = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(resultExpr, 0, allocVector(REALSXP, bartResults->getNumSigmaSamples()));
    SET_VECTOR_ELT(resultExpr, 1, allocVector(REALSXP, bartResults->getNumTrainingSamples()));
    SET_VECTOR_ELT(resultExpr, 2, allocVector(REALSXP, bartResults->getNumTestSamples()));
    SET_VECTOR_ELT(resultExpr, 3, allocVector(INTSXP, bartResults->getNumVariableCountSamples()));
    
    SEXP sigmaSamples = VECTOR_ELT(resultExpr, 0);
    std::memcpy(REAL(sigmaSamples), (const double*) bartResults->sigmaSamples, bartResults->getNumSigmaSamples() * sizeof(double));
    
    SEXP trainingSamples = VECTOR_ELT(resultExpr, 1);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
    dims = INTEGER(dimsExpr);
    dims[0] = bartResults->numObservations;
    dims[1] = bartResults->numSamples;
    setAttrib(trainingSamples, R_DimSymbol, dimsExpr);
    std::memcpy(REAL(trainingSamples), (const double*) bartResults->trainingSamples, bartResults->getNumTrainingSamples() * sizeof(double));
    
    SEXP testSamples = VECTOR_ELT(resultExpr, 2);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
    dims = INTEGER(dimsExpr);
    dims[0] = bartResults->numTestObservations;
    dims[1] = bartResults->numSamples;
    setAttrib(testSamples, R_DimSymbol, dimsExpr);
    std::memcpy(REAL(testSamples), (const double*) bartResults->testSamples, bartResults->getNumTestSamples() * sizeof(double));
    
    SEXP variableCountSamples = VECTOR_ELT(resultExpr, 3);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
    dims = INTEGER(dimsExpr);
    dims[0] = bartResults->numPredictors;
    dims[1] = bartResults->numSamples;
    setAttrib(variableCountSamples, R_DimSymbol, dimsExpr);
    int* variableCountStorage = INTEGER(variableCountSamples);
    size_t length = bartResults->getNumVariableCountSamples();
    // these likely need to be down-sized from 64 to 32 bits
    for (size_t i = 0; i < length; ++i) variableCountStorage[i] = (int) bartResults->variableCountSamples[i];
    
    setAttrib(resultExpr, R_NamesSymbol, namesExpr = allocVector(STRSXP, 4));
    SET_STRING_ELT(namesExpr, 0, mkChar("sigma"));
    SET_STRING_ELT(namesExpr, 1, mkChar("train"));
    SET_STRING_ELT(namesExpr, 2, mkChar("test"));
    SET_STRING_ELT(namesExpr, 3, mkChar("varcount"));
    
    
    delete bartResults;
    
    delete [] u_maxNumCuts;
    
    
    UNPROTECT(4);
    
    return(resultExpr);
  }
  
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
  
  static R_CallMethodDef callMethods[] = {
    CALLDEF(cbart_fitBart, 21),
    {NULL, NULL, 0}
  };
}
#include <external/io.h>
extern "C" {
  void R_init_cbart(DllInfo* info)
  {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    
    R_RegisterCCallable("cbart", "bart_createCGMPrior", (DL_FUNC) bart_createCGMPrior);
    R_RegisterCCallable("cbart", "bart_createCGMPriorFromControl", (DL_FUNC) bart_createCGMPriorFromControl);
    R_RegisterCCallable("cbart", "bart_destroyCGMPrior", (DL_FUNC) bart_destroyCGMPrior);
    R_RegisterCCallable("cbart", "bart_initializeCGMPriorFromControl", (DL_FUNC) bart_initializeCGMPriorFromControl);
    R_RegisterCCallable("cbart", "bart_invalidateCGMPrior", (DL_FUNC) bart_invalidateCGMPrior);
    
    R_RegisterCCallable("cbart", "bart_createNormalPrior", (DL_FUNC) bart_createNormalPrior);
    R_RegisterCCallable("cbart", "bart_createNormalPriorFromControl", (DL_FUNC) bart_createNormalPriorFromControl);
    R_RegisterCCallable("cbart", "bart_destroyNormalPrior", (DL_FUNC) bart_destroyNormalPrior);
    R_RegisterCCallable("cbart", "bart_initializeNormalPriorFromControl", (DL_FUNC) bart_initializeNormalPriorFromControl);
    R_RegisterCCallable("cbart", "bart_invalidateNormalPrior", (DL_FUNC) bart_invalidateNormalPrior);
    
    R_RegisterCCallable("cbart", "bart_createChiSquaredPrior", (DL_FUNC) bart_createChiSquaredPrior);
    R_RegisterCCallable("cbart", "bart_createChiSquaredPriorFromControl", (DL_FUNC) bart_createChiSquaredPriorFromControl);
    R_RegisterCCallable("cbart", "bart_destroyChiSquaredPrior", (DL_FUNC) bart_destroyChiSquaredPrior);
    R_RegisterCCallable("cbart", "bart_initializeChiSquaredPriorFromControl", (DL_FUNC) bart_initializeChiSquaredPriorFromControl);
    R_RegisterCCallable("cbart", "bart_invalidateChiSquaredPrior", (DL_FUNC) bart_invalidateChiSquaredPrior);

    R_RegisterCCallable("cbart", "bart_createFit", (DL_FUNC) bart_createFit);
    R_RegisterCCallable("cbart", "bart_initializeFit", (DL_FUNC) bart_initializeFit);
    R_RegisterCCallable("cbart", "bart_destroyFit", (DL_FUNC) bart_destroyFit);
    R_RegisterCCallable("cbart", "bart_invalidateFit", (DL_FUNC) bart_invalidateFit);
    
    R_RegisterCCallable("cbart", "bart_runSampler", (DL_FUNC) bart_runSampler);
    R_RegisterCCallable("cbart", "bart_setResponse", (DL_FUNC) bart_setResponse);
  }
}
