#include "config.hpp"
#include "R_interface_sampler.hpp"

#include <cstddef>
#include <cstring> // strcmp, memcpy

#include <external/alloca.h>

#include <R_ext/Random.h> // GetRNGstate, PutRNGState

#include <rc/bounds.h>
#include <rc/util.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>

#include "R_interface.hpp"
#include "R_interface_common.hpp"

#define asRXLen(_X_) static_cast<R_xlen_t>(_X_)

using std::size_t;
using namespace dbarts;

extern "C" {
  static void fitFinalizer(SEXP fitExpr);

  SEXP create(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr)
  {
    Control control;
    Model model;
    Data data;
    
    SEXP classExpr = Rf_getAttrib(controlExpr, R_ClassSymbol);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsControl") != 0) Rf_error("'control' argument to dbarts_create not of class 'dbartsControl'");
    
    classExpr = Rf_getAttrib(modelExpr, R_ClassSymbol);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsModel") != 0) Rf_error("'model' argument to dbarts_create not of class 'dbartsModel'");
    
    classExpr = Rf_getAttrib(dataExpr, R_ClassSymbol);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsData") != 0) Rf_error("'data' argument to dbarts_create not of class 'dbartsData'");
    
    
    initializeControlFromExpression(control, controlExpr);
    initializeModelFromExpression(model, modelExpr, control);
    initializeDataFromExpression(data, dataExpr);
    
    BARTFit* fit = new BARTFit(control, model, data);
    
    SEXP result = PROTECT(R_MakeExternalPtr(fit, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(result, fitFinalizer, static_cast<Rboolean>(FALSE));
    
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_lock(&fitMutex);
    Rprintf("creating   %p\n", fit);
#endif
    activeFits->insert(result);
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_unlock(&fitMutex);
#endif

    UNPROTECT(1);
    
    return result;
  }
  
  SEXP run(SEXP fitExpr, SEXP numBurnInExpr, SEXP numSamplesExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_run called on NULL external pointer");
    
    int i_temp;
    size_t numBurnIn, numSamples;
    
    i_temp = rc_getInt(numBurnInExpr, "number of burn-in steps", RC_LENGTH | RC_GEQ, asRXLen(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES, RC_END);
    numBurnIn = i_temp == NA_INTEGER ? fit->control.numBurnIn : static_cast<size_t>(i_temp);
    
    i_temp = rc_getInt(numSamplesExpr, "number of samples", RC_LENGTH | RC_GEQ, asRXLen(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES, RC_END);    
    numSamples = i_temp == NA_INTEGER ? fit->control.numSamples : static_cast<size_t>(i_temp);
    
    if (numBurnIn == 0 && numSamples == 0) Rf_error("either number of burn-in or samples must be positive");
    
    size_t numTrainingSamples = fit->data.numObservations * numSamples;
    if (numSamples != 0 && numTrainingSamples / numSamples != fit->data.numObservations)
      Rf_error("training sample array size exceeds architecture's capacity");
    R_xlen_t s_numTrainingSamples = asRXLen(numTrainingSamples);
    if (s_numTrainingSamples < 0 || static_cast<size_t>(s_numTrainingSamples) != numTrainingSamples)
      Rf_error("training sample array size cannot be represented by a signed integer on this architecture");
    
    size_t numTestSamples = fit->data.numTestObservations * numSamples;
     if (numSamples != 0 && numTestSamples / numSamples != fit->data.numTestObservations)
      Rf_error("test sample array size exceeds architecture's capacity");
    R_xlen_t s_numTestSamples = asRXLen(numTestSamples);
    if (s_numTestSamples < 0 || static_cast<size_t>(s_numTestSamples) != numTestSamples)
      Rf_error("test sample array size cannot be represented by a signed integer on this architecture");
    
    GetRNGstate();
    
    Results* bartResults = fit->runSampler(numBurnIn, numSamples);
    
    PutRNGstate();
    
    // can happen if numSamples == 0
    if (bartResults == NULL) return R_NilValue;
    
    int protectCount = 0;
    
    SEXP resultExpr = PROTECT(rc_newList(4));
    ++protectCount;
    SET_VECTOR_ELT(resultExpr, 0, rc_newNumeric(asRXLen(bartResults->getNumSigmaSamples())));
    SET_VECTOR_ELT(resultExpr, 1, rc_newNumeric(asRXLen(bartResults->getNumTrainingSamples())));
    if (fit->data.numTestObservations > 0)
      SET_VECTOR_ELT(resultExpr, 2, rc_newNumeric(asRXLen(bartResults->getNumTestSamples())));
    else
      SET_VECTOR_ELT(resultExpr, 2, R_NilValue);
    SET_VECTOR_ELT(resultExpr, 3, rc_newInteger(asRXLen(bartResults->getNumVariableCountSamples())));
    
    SEXP sigmaSamples = VECTOR_ELT(resultExpr, 0);
    if (fit->control.numChains > 1)
      rc_setDims(sigmaSamples, static_cast<int>(bartResults->numSamples), static_cast<int>(fit->control.numChains), -1);
    std::memcpy(REAL(sigmaSamples), const_cast<const double*>(bartResults->sigmaSamples), bartResults->getNumSigmaSamples() * sizeof(double));
    
    SEXP trainingSamples = VECTOR_ELT(resultExpr, 1);
    if (fit->control.numChains <= 1)
      rc_setDims(trainingSamples, static_cast<int>(bartResults->numObservations), static_cast<int>(bartResults->numSamples), -1);
    else
      rc_setDims(trainingSamples, static_cast<int>(bartResults->numObservations), static_cast<int>(bartResults->numSamples), static_cast<int>(fit->control.numChains), -1);
    std::memcpy(REAL(trainingSamples), const_cast<const double*>(bartResults->trainingSamples), bartResults->getNumTrainingSamples() * sizeof(double));
    
    if (fit->data.numTestObservations > 0) {
      SEXP testSamples = VECTOR_ELT(resultExpr, 2);
      if (fit->control.numChains <= 1)
        rc_setDims(testSamples, static_cast<int>(bartResults->numTestObservations), static_cast<int>(bartResults->numSamples), -1);
      else
        rc_setDims(testSamples, static_cast<int>(bartResults->numTestObservations), static_cast<int>(bartResults->numSamples), static_cast<int>(fit->control.numChains), -1);
      std::memcpy(REAL(testSamples), const_cast<const double*>(bartResults->testSamples), bartResults->getNumTestSamples() * sizeof(double));
    }
    
    SEXP variableCountSamples = VECTOR_ELT(resultExpr, 3);
    if (fit->control.numChains <= 1)
      rc_setDims(variableCountSamples, static_cast<int>(bartResults->numPredictors), static_cast<int>(bartResults->numSamples), -1);
    else
      rc_setDims(variableCountSamples, static_cast<int>(bartResults->numPredictors), static_cast<int>(bartResults->numSamples), static_cast<int>(fit->control.numChains), -1);
    int* variableCountStorage = INTEGER(variableCountSamples);
    size_t length = bartResults->getNumVariableCountSamples();
    // these likely need to be down-sized from 64 to 32 bits
    for (size_t i = 0; i < length; ++i) variableCountStorage[i] = static_cast<int>(bartResults->variableCountSamples[i]);
    
        
    // create result storage and make it user friendly
    SEXP namesExpr;
    
    rc_setNames(resultExpr, namesExpr = rc_newCharacter(4));
    SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
    SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
    SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
    SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
    
    UNPROTECT(protectCount);
    
    delete bartResults;
    
    return resultExpr;
  }
  
  SEXP sampleTreesFromPrior(SEXP fitExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_sampleTreesFromPrior called on NULL external pointer");
        
    GetRNGstate();
    
    fit->sampleTreesFromPrior();
    
    PutRNGstate();
    
    return R_NilValue;
  }
  
  SEXP setData(SEXP fitExpr, SEXP dataExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setData called on NULL external pointer");
    
    SEXP classExpr = Rf_getAttrib(dataExpr, R_ClassSymbol);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsData") != 0) Rf_error("'data' argument to dbarts_setData not of class 'dbartsData'");
    
    Data data;
    initializeDataFromExpression(data, dataExpr);
    
    Data oldData = fit->data;
    
    if (data.numPredictors != oldData.numPredictors) {
      delete [] data.maxNumCuts;
      delete [] data.variableTypes;
      Rf_error("number of predictors between old and new data must be the same");
    }
    
    fit->setData(data);
    
    delete [] oldData.maxNumCuts;
    delete [] oldData.variableTypes;
    
    return R_NilValue;
  }
  
  SEXP setControl(SEXP fitExpr, SEXP controlExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setControl called on NULL external pointer");
    
    if (std::strcmp(CHAR(STRING_ELT(Rf_getAttrib(controlExpr, R_ClassSymbol), 0)), "dbartsControl") != 0) Rf_error("'control' argument to dbarts_setControl not of class 'dbartsControl'");
    
    Control control;
    initializeControlFromExpression(control, controlExpr);
    
    Control oldControl = fit->control;
    
    if (control.responseIsBinary != oldControl.responseIsBinary)
      Rf_error("new control cannot change binary characteristic of response");
    if (control.numChains != oldControl.numChains)
      Rf_error("new control cannot change number of chains");
    
    fit->setControl(control);
    
    return R_NilValue;
  }
  
  SEXP setModel(SEXP fitExpr, SEXP modelExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setControl called on NULL external pointer");
    
    if (std::strcmp(CHAR(STRING_ELT(Rf_getAttrib(modelExpr, R_ClassSymbol), 0)), "dbartsModel") != 0) Rf_error("'model' argument to dbarts_setModel not of class 'dbartsModel'");
    
    Model model;
    initializeModelFromExpression(model, modelExpr, fit->control);
    
    Model oldModel = fit->model;
    
    fit->setModel(model);
    
    invalidateModel(oldModel);
    
    return R_NilValue;
  }
  
  SEXP setResponse(SEXP fitExpr, SEXP y)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setY called on NULL external pointer");
    
    rc_assertDoubleConstraints(y, "y", RC_LENGTH | RC_EQ, asRXLen(fit->data.numObservations));
    
    // for binary responses, updates latents and samples
    if (fit->control.responseIsBinary) GetRNGstate();
    
    fit->setResponse(REAL(y));
    
    if (fit->control.responseIsBinary) PutRNGstate();
    
    return R_NilValue;
  }
  
  SEXP setOffset(SEXP fitExpr, SEXP offsetExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setOffset called on NULL external pointer");
    
    double* offset = NULL;
    if (Rf_isReal(offsetExpr)) {
      offset = REAL(offsetExpr);
      if (rc_getLength(offsetExpr) != fit->data.numObservations) Rf_error("length of new offset does not match y");
    } else if (!Rf_isNull(offsetExpr) && !rc_isS4Null(offsetExpr)) {
      Rf_error("offset must be of type real or NULL");
    }
    
    // for binary responses, updates latents and samples
    if (fit->control.responseIsBinary) GetRNGstate();
    
    fit->setOffset(offset);
    
    if (fit->control.responseIsBinary) PutRNGstate();
    
    return R_NilValue;
  }
  
  SEXP setPredictor(SEXP fitExpr, SEXP x)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setPredictor called on NULL external pointer");
    
    if (!Rf_isReal(x)) Rf_error("x must be of type real");
    
    rc_assertDimConstraints(x, "dimensions of x", RC_LENGTH | RC_EQ, rc_asRLength(2),
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numObservations),
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numPredictors),
                            RC_END);
    
    // SEXP dimsExpr = Rf_getAttrib(x, R_DimSymbol);
    //
    //if (Rf_isNull(dimsExpr) || rc_getLength(dimsExpr) != 2) Rf_error("x must be a matrix, i.e. have two dimensions");
    //int* dims = INTEGER(dimsExpr);
    //
    //if (static_cast<size_t>(dims[0]) != fit->data.numObservations) Rf_error("number of rows in new x does not match y");
    //if (static_cast<size_t>(dims[1]) != fit->data.numPredictors) Rf_error("number of columns in new x does not match old");
    
    return Rf_ScalarLogical(fit->setPredictor(REAL(x)));
  }
  
  SEXP setTestPredictor(SEXP fitExpr, SEXP x_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setTestPredictor called on NULL external pointer");
    
    if (Rf_isNull(x_test) || rc_isS4Null(x_test)) {
      fit->setTestPredictor(NULL, 0);
    
      return R_NilValue;
    }
    
    if (!Rf_isReal(x_test)) Rf_error("x.test must be of type real");
    
    rc_assertDimConstraints(x_test, "dimensions of x_test", RC_LENGTH | RC_EQ, rc_asRLength(2),
                            RC_NA,
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numPredictors),
                            RC_END);
    int* dims = INTEGER(Rf_getAttrib(x_test, R_DimSymbol));
    
    // SEXP dimsExpr = Rf_getAttrib(x_test, R_DimSymbol);
    // if (rc_getLength(dimsExpr) != 2) Rf_error("x.test must be a matrix, i.e. have two dimensions");
    // int* dims = INTEGER(dimsExpr);
    // if (static_cast<size_t>(dims[1]) != fit->data.numPredictors) Rf_error("number of columns of x.test and x must be equal");
    
    fit->setTestPredictor(REAL(x_test), static_cast<size_t>(dims[0]));
    
    return R_NilValue;
  }
  
  SEXP setTestOffset(SEXP fitExpr, SEXP offset_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setTestOffset called on NULL external pointer");
    
    if (Rf_isNull(offset_test)) {
      fit->setTestOffset(NULL);
    } else {
      if (!Rf_isReal(offset_test)) Rf_error("offset.test must be of type real");
      if (fit->data.numTestObservations != rc_getLength(offset_test)) Rf_error("length of offset.test must equal number of rows in x.test");
      fit->setTestOffset(REAL(offset_test));
    }
    
    return R_NilValue;
  }
  
  SEXP setTestPredictorAndOffset(SEXP fitExpr, SEXP x_test, SEXP offset_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setTestPredictorAndOffset called on NULL external pointer");
    
    if (Rf_isNull(x_test) || rc_isS4Null(x_test)) {
      fit->setTestPredictor(NULL, 0);
    
      return R_NilValue;
    }
    
    if (!Rf_isReal(x_test)) Rf_error("x.test must be of type real");
    
    rc_assertDimConstraints(x_test, "dimensions of x_test", RC_LENGTH | RC_EQ, rc_asRLength(2),
                            RC_NA,
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numPredictors),
                            RC_END);
    int* dims = INTEGER(Rf_getAttrib(x_test, R_DimSymbol));
    
    // SEXP dimsExpr = Rf_getAttrib(x_test, R_DimSymbol);
    // if (rc_getLength(dimsExpr) != 2) Rf_error("x.test must be a matrix, i.e. have two dimensions");
    // int* dims = INTEGER(dimsExpr);
    // if (static_cast<size_t>(dims[1]) != fit->data.numPredictors) Rf_error("number of columns of x.test and x must be equal");
    
    if (Rf_isNull(offset_test)) {
      fit->setTestPredictorAndOffset(REAL(x_test), NULL, static_cast<size_t>(dims[0]));
    } else {
      if (!Rf_isReal(offset_test)) Rf_error("offset.test must be of type real");
      if (rc_getLength(offset_test) == 1 && ISNA(REAL(offset_test)[0])) {
        fit->setTestPredictor(REAL(x_test), static_cast<size_t>(dims[0]));
      } else {
        if (rc_getLength(offset_test) != static_cast<size_t>(dims[0])) Rf_error("length of offset.test must equal number of rows in x.test");
        fit->setTestPredictorAndOffset(REAL(x_test), REAL(offset_test), static_cast<size_t>(dims[0]));
      }
    }
    
    return R_NilValue;
  }
  
  
  SEXP updatePredictor(SEXP fitExpr, SEXP x, SEXP colsExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_updatePredictor called on NULL external pointer");
    
    if (!Rf_isReal(x)) Rf_error("x must be of type real");
    if (!Rf_isInteger(colsExpr)) Rf_error("columns must be of type integer");
    
    SEXP dimsExpr = Rf_getAttrib(x, R_DimSymbol);
    int* dims = NULL;
    
    if (!Rf_isNull(dimsExpr)) {
      size_t numDims = rc_getLength(dimsExpr);
      
      if (numDims != 1 && numDims != 2) Rf_error("x must be a vector or a matrix");
      if (numDims == 2) dims = INTEGER(dimsExpr);
    }
    
    if (rc_getLength(colsExpr) == 0) Rf_error("length of columns is 0");

    if (dims != NULL) {
      if (static_cast<size_t>(dims[0]) != fit->data.numObservations) Rf_error("number of rows of new x does not match y");
      if (static_cast<size_t>(dims[1]) != rc_getLength(colsExpr)) Rf_error("number of columns of new x does not match length of columns to replace");
    } else {
      if (rc_getLength(x) != fit->data.numObservations) Rf_error("length of new x does not match y");
    }
    
    
    int* colsInt = INTEGER(colsExpr);
    size_t numCols = rc_getLength(colsExpr);
    size_t* cols = ext_stackAllocate(numCols, size_t);
    for (size_t i = 0 ; i < numCols; ++i) {
      cols[i] = static_cast<size_t>(colsInt[i] - 1);
      if (static_cast<size_t>(cols[i]) >= fit->data.numPredictors) {
        ext_stackFree(cols);
        Rf_error("column '%d' is out of range", colsInt[i]);
      }
    }
    
    bool result = fit->updatePredictors(REAL(x), cols, numCols);
    
    ext_stackFree(cols);
    
    return Rf_ScalarLogical(result);
  }
  
  SEXP updateTestPredictor(SEXP fitExpr, SEXP x_test, SEXP colsExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_updateTestPredictor called on NULL external pointer");
    
    if (fit->data.x_test == NULL) Rf_error("test matrix must exist at object creation to be updated");
    
    if (!Rf_isReal(x_test)) Rf_error("x must be of type real");
    if (!Rf_isInteger(colsExpr)) Rf_error("columns must be of type integer");
    
    SEXP dimsExpr = Rf_getAttrib(x_test, R_DimSymbol);
    int* dims = NULL;
    
    if (!Rf_isNull(dimsExpr)) {
      size_t numDims = rc_getLength(dimsExpr);
      
      if (numDims != 1 && numDims != 2) Rf_error("x must be a vector or a matrix");
      if (numDims == 2) dims = INTEGER(dimsExpr);
    }
    
    if (rc_getLength(colsExpr) == 0) Rf_error("length of columns is 0");

    if (dims != NULL) {
      if (static_cast<size_t>(dims[0]) != fit->data.numTestObservations) Rf_error("number of rows of new x does not match old x.test");
      if (static_cast<size_t>(dims[1]) != rc_getLength(colsExpr)) Rf_error("number of columns of new x does not match length of columns to replace");
    } else {
      if (rc_getLength(x_test) != fit->data.numTestObservations) Rf_error("length of new x does not match old x.test");
    }
    
    
    int* colsInt = INTEGER(colsExpr);
    size_t numCols = rc_getLength(colsExpr);
    size_t* cols = ext_stackAllocate(numCols, size_t);
    for (size_t i = 0 ; i < numCols; ++i) {
      cols[i] = static_cast<size_t>(colsInt[i] - 1);
      if (cols[i] >= fit->data.numPredictors) {
        ext_stackFree(cols);
        Rf_error("column '%d' is out of range", colsInt[i]);
      }
    }
    
    fit->updateTestPredictors(REAL(x_test), cols, numCols);
    
    ext_stackFree(cols);
    
    return R_NilValue;
  }
  
  
  SEXP createState(SEXP fitExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_createState called on NULL external pointer");
    
    return createStateExpressionFromFit(*fit);
  }
  
  SEXP restoreState(SEXP fitExpr, SEXP stateExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_restoreState called on NULL external pointer");
    
    initializeStateFromExpression(*fit, fit->state, stateExpr);
    
    return R_NilValue;
  }
  
  SEXP storeState(SEXP fitExpr, SEXP stateExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_storeState called on NULL external pointer");
    
    storeStateExpressionFromFit(*fit, stateExpr);
    
    return R_NilValue;
  }
  
  
  SEXP printTrees(SEXP fitExpr, SEXP chainIndicesExpr, SEXP treeIndicesExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_printTrees called on NULL external pointer");
    
    size_t numChains = fit->control.numChains;
    size_t numTrees = fit->control.numTrees;
    
    size_t numChainIndices = Rf_isNull(chainIndicesExpr) ? numChains : rc_getLength(chainIndicesExpr);
    size_t numTreeIndices  = Rf_isNull(treeIndicesExpr)  ? numTrees  : rc_getLength(treeIndicesExpr);
    
    
    size_t* chainIndices = ext_stackAllocate(numChainIndices, size_t);
    size_t* treeIndices  = ext_stackAllocate(numTreeIndices,  size_t);
    
    if (Rf_isNull(chainIndicesExpr)) {
      for (size_t i = 0; i < numChains; ++i) chainIndices[i] = i;
    } else {
      int* i_chainIndices = INTEGER(chainIndicesExpr);
      for (size_t i = 0; i < numChainIndices; ++i) chainIndices[i] = static_cast<size_t>(i_chainIndices[i] - 1);
    }
    
    
    if (Rf_isNull(treeIndicesExpr)) {
      for (size_t i = 0; i < numTrees; ++i) treeIndices[i] = i;
    } else {
      int* i_treeIndices = INTEGER(treeIndicesExpr);
      for (size_t i = 0; i < numTreeIndices; ++i) treeIndices[i] = static_cast<size_t>(i_treeIndices[i] - 1);
    }
   
    fit->printTrees(chainIndices, numChainIndices, treeIndices, numTreeIndices);
    
    ext_stackFree(treeIndices);
    ext_stackFree(chainIndices);
    
    return R_NilValue;
  }

  
  SEXP saveToFile(SEXP fitExpr, SEXP fileName)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_saveToFile called on NULL external pointer");
    
    return Rf_ScalarLogical(fit->saveToFile(CHAR(STRING_ELT(fileName, 0))));
  }

  SEXP loadFromFile(SEXP fileName)
  {
    BARTFit* fit = BARTFit::loadFromFile(CHAR(STRING_ELT(fileName, 0)));
    
    delete [] fit->data.maxNumCuts;
    delete [] fit->data.variableTypes;
    
    delete fit->model.sigmaSqPrior;
    delete fit->model.muPrior;
    delete fit->model.treePrior;
    
    delete fit;
    
    /* TODO: turn into an R object */
    
    return R_NilValue;
  }
  
  
  static void fitFinalizer(SEXP fitExpr)
  {
#ifdef THREAD_SAFE_UNLOAD
    Rprintf("finalizing ");
#endif
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
#ifdef THREAD_SAFE_UNLOAD
    Rprintf("%p\n", fit);
#endif
    if (fit == NULL) return;
    
    
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_lock(&fitMutex);
#endif
    if (activeFits->find(fitExpr) == activeFits->end()) {
#ifdef THREAD_SAFE_UNLOAD
      pthread_mutex_unlock(&fitMutex);
#endif
      return;
    }
    activeFits->erase(fitExpr);
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_unlock(&fitMutex);
#endif
    
    deleteFit(fit);
    
    R_ClearExternalPtr(fitExpr);
  }
}

