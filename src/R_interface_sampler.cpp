#include "config.hpp"
#ifdef __MINGW32__
#  define __USE_MINGW_ANSI_STDIO 1
#endif
#include "R_interface_sampler.hpp"

#include <cstddef>
#ifdef HAVE_STD_SNPRINTF
// snprintf in c++11, before that have to use C version
#  include <cstdio>
using std::snprintf;
#else
#  include <stdio.h>
#endif
#include <cstring> // strcmp, memcpy

#include <misc/alloca.h>

// R includes
#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#define USE_FC_LEN_T
#endif

#include <R_ext/Random.h> // GetRNGstate, PutRNGState

#undef USE_FC_LEN_T

#include <rc/bounds.h>
#include <rc/util.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>

#include "R_interface.hpp"
#include "R_interface_common.hpp"

using std::size_t;
using namespace dbarts;

extern "C" {
  static void fitFinalizer(SEXP fitExpr);

  SEXP create(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr)
  {
    Control control;
    Data data;
    Model model;
    
    SEXP classExpr = Rf_getAttrib(controlExpr, R_ClassSymbol);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsControl") != 0) Rf_error("'control' argument to dbarts_create not of class 'dbartsControl'");
    
    classExpr = Rf_getAttrib(modelExpr, R_ClassSymbol);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsModel") != 0) Rf_error("'model' argument to dbarts_create not of class 'dbartsModel'");
    
    classExpr = Rf_getAttrib(dataExpr, R_ClassSymbol);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsData") != 0) Rf_error("'data' argument to dbarts_create not of class 'dbartsData'");
    
    
    initializeControlFromExpression(control, controlExpr);
    initializeDataFromExpression(data, dataExpr);
    initializeModelFromExpression(model, modelExpr, control, data);
    
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
    
    i_temp = rc_getInt(numBurnInExpr, "number of burn-in steps", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES, RC_END);
    numBurnIn = i_temp == NA_INTEGER ? fit->control.defaultNumBurnIn : static_cast<size_t>(i_temp);
    
    i_temp = rc_getInt(numSamplesExpr, "number of samples", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES, RC_END);    
    numSamples = i_temp == NA_INTEGER ? fit->control.defaultNumSamples : static_cast<size_t>(i_temp);
    
    if (numBurnIn == 0 && numSamples == 0) Rf_error("either number of burn-in or samples must be positive");
    
    size_t numTrainingSamples = fit->data.numObservations * numSamples;
    if (numSamples != 0 && numTrainingSamples / numSamples != fit->data.numObservations)
      Rf_error("training sample array size exceeds architecture's capacity");
    R_xlen_t s_numTrainingSamples = rc_asRLength(numTrainingSamples);
    if (s_numTrainingSamples < 0 || static_cast<size_t>(s_numTrainingSamples) != numTrainingSamples)
      Rf_error("training sample array size cannot be represented by a signed integer on this architecture");
    
    size_t numTestSamples = fit->data.numTestObservations * numSamples;
     if (numSamples != 0 && numTestSamples / numSamples != fit->data.numTestObservations)
      Rf_error("test sample array size exceeds architecture's capacity");
    R_xlen_t s_numTestSamples = rc_asRLength(numTestSamples);
    if (s_numTestSamples < 0 || static_cast<size_t>(s_numTestSamples) != numTestSamples)
      Rf_error("test sample array size cannot be represented by a signed integer on this architecture");
    
    GetRNGstate();
    
    Results* bartResults = fit->runSampler(numBurnIn, numSamples);
    
    PutRNGstate();
    
    // can happen if numSamples == 0
    if (bartResults == NULL) return R_NilValue;
    
    int protectCount = 0;
    
    SEXP resultExpr = PROTECT(rc_newList(bartResults->kSamples == NULL ? 4 : 5));
    ++protectCount;
    SET_VECTOR_ELT(resultExpr, 0, rc_newReal(rc_asRLength(bartResults->getNumSigmaSamples())));
    SET_VECTOR_ELT(resultExpr, 1, rc_newReal(rc_asRLength(bartResults->getNumTrainingSamples())));
    if (fit->data.numTestObservations > 0)
      SET_VECTOR_ELT(resultExpr, 2, rc_newReal(rc_asRLength(bartResults->getNumTestSamples())));
    else
      SET_VECTOR_ELT(resultExpr, 2, R_NilValue);
    SET_VECTOR_ELT(resultExpr, 3, rc_newInteger(rc_asRLength(bartResults->getNumVariableCountSamples())));
    if (bartResults->kSamples != NULL)
      SET_VECTOR_ELT(resultExpr, 4, rc_newReal(rc_asRLength(bartResults->getNumSigmaSamples())));
    
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
    for (size_t i = 0; i < length; ++i)
      variableCountStorage[i] = static_cast<int>(bartResults->variableCountSamples[i]);
    
    if (bartResults->kSamples != NULL) {
      SEXP kSamples = VECTOR_ELT(resultExpr, 4);
      if (fit->control.numChains > 1)
        rc_setDims(kSamples, static_cast<int>(bartResults->numSamples), static_cast<int>(fit->control.numChains), -1);
      std::memcpy(REAL(kSamples), const_cast<const double*>(bartResults->kSamples), bartResults->getNumSigmaSamples() * sizeof(double));
    }
        
    // create result storage and make it user friendly
    SEXP namesExpr;
    
    rc_setNames(resultExpr, namesExpr = rc_newCharacter(bartResults->kSamples == NULL ? 4 : 5));
    SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
    SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
    SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
    SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
    if (bartResults->kSamples != NULL)
      SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));
    
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
  
  SEXP sampleNodeParametersFromPrior(SEXP fitExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_sampleNodeParametersFromPrior called on NULL external pointer");
        
    GetRNGstate();
    
    fit->sampleNodeParametersFromPrior();
    
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
    initializeModelFromExpression(model, modelExpr, fit->control, fit->data);
    
    Model oldModel = fit->model;
    
    if ((model.kPrior != NULL && oldModel.kPrior == NULL) ||
        (model.kPrior == NULL && oldModel.kPrior != NULL))
    {
      Rf_error("k prior cannot be changed after sampler has been created");
      invalidateModel(model);
      
      return R_NilValue;
    }
    
    fit->setModel(model);
    
    invalidateModel(oldModel);
    
    return R_NilValue;
  }
  
  SEXP predict(SEXP fitExpr, SEXP x_testExpr, SEXP offset_testExpr)
  {
    const BARTFit* fit = static_cast<const BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_predict called on NULL external pointer");
    
    const Control& control(fit->control);
    
    if (Rf_isNull(x_testExpr) || rc_isS4Null(x_testExpr)) return R_NilValue;
    
    if (!Rf_isReal(x_testExpr)) Rf_error("x.test must be of type real");
    
    rc_assertDimConstraints(x_testExpr, "dimensions of x_test", RC_LENGTH | RC_EQ, rc_asRLength(2),
                            RC_NA,
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numPredictors),
                            RC_END);
    int* dims = INTEGER(Rf_getAttrib(x_testExpr, R_DimSymbol));
    
    
    size_t numSamples = control.keepTrees ? fit->currentNumSamples : 1;
    size_t numTestObservations = static_cast<size_t>(dims[0]);
    
    
    double* testOffset = NULL;
    if (!Rf_isNull(offset_testExpr)) {
      if (!Rf_isReal(offset_testExpr)) Rf_error("offset.test must be of type real");
      if (rc_getLength(offset_testExpr) != 1 || !ISNA(REAL(offset_testExpr)[0])) {
        if (rc_getLength(offset_testExpr) != numTestObservations) Rf_error("length of offset.test must equal number of rows in x.test");
        testOffset = REAL(offset_testExpr);
      }
    }
    
    SEXP result = PROTECT(Rf_allocVector(REALSXP, numTestObservations * numSamples * control.numChains));
    if (control.keepTrees) {
      if (fit->control.numChains <= 1)
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(numSamples), -1);
      else
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(numSamples), static_cast<int>(control.numChains), -1);
    } else {
      if (fit->control.numChains > 1)
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(control.numChains), -1);
    }
    
    fit->predict(REAL(x_testExpr), numTestObservations, testOffset, REAL(result));
    
    UNPROTECT(1);
    
    return result;
  }
  
  SEXP setResponse(SEXP fitExpr, SEXP y)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setResponse called on NULL external pointer");
    
    rc_assertDoubleConstraints(y, "y", RC_LENGTH | RC_EQ, rc_asRLength(fit->data.numObservations), RC_END);
    
    // for binary responses, updates latents and samples
    if (fit->control.responseIsBinary) GetRNGstate();
    
    fit->setResponse(REAL(y));
    
    if (fit->control.responseIsBinary) PutRNGstate();
    
    return R_NilValue;
  }
  
  SEXP setOffset(SEXP fitExpr, SEXP offsetExpr, SEXP updateScaleExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setOffset called on NULL external pointer");
    
    const double* offset = NULL;
    if (Rf_isReal(offsetExpr)) {
      offset = REAL(offsetExpr);
      if (rc_getLength(offsetExpr) != fit->data.numObservations) Rf_error("length of new offset does not match y");
    } else if (!Rf_isNull(offsetExpr) && !rc_isS4Null(offsetExpr)) {
      Rf_error("offset must be of type real or NULL");
    }
    
    bool updateScale = rc_getBool(updateScaleExpr, "updateScale", RC_DEFAULT | RC_VALUE, false, RC_END);
    
    fit->setOffset(offset, updateScale);
    
    return R_NilValue;
  }
  
  SEXP setSigma(SEXP fitExpr, SEXP sigmaExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setSigma called on NULL external pointer");
    
    const double* sigma = NULL;
    if (Rf_isReal(sigmaExpr)) {
      sigma = REAL(sigmaExpr);
      if (rc_getLength(sigmaExpr) != fit->control.numChains) Rf_error("length of new sigma does not match number of chains");
    } else {
      Rf_error("sigma must be of type real");
    }
    
    fit->setSigma(sigma);
    
    return R_NilValue;
  }
  
  SEXP setWeights(SEXP fitExpr, SEXP weights)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setWeights called on NULL external pointer");
    
    rc_assertDoubleConstraints(weights, "weights", RC_LENGTH | RC_EQ, rc_asRLength(fit->data.numObservations),
                               RC_VALUE | RC_GEQ, 0.0,
                               RC_END);
    fit->setWeights(REAL(weights));
    
    return R_NilValue;
  }
  
  SEXP setPredictor(SEXP fitExpr, SEXP xExpr, SEXP forceUpdateExpr, SEXP updateCutPointsExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setPredictor called on NULL external pointer");
    
    if (!Rf_isReal(xExpr)) Rf_error("x must be of type real");
    
    bool forceUpdate     = rc_getBool(forceUpdateExpr,         "forceUpdate", RC_NA | RC_NO, RC_END);
    bool updateCutPoints = rc_getBool(updateCutPointsExpr, "updateCutPoints", RC_NA | RC_NO, RC_END);
    
    rc_assertDimConstraints(xExpr, "dimensions of x", RC_LENGTH | RC_EQ, rc_asRLength(2),
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numObservations),
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numPredictors),
                            RC_END);
      
    bool result = fit->setPredictor(REAL(xExpr), forceUpdate, updateCutPoints);
        
    return Rf_ScalarLogical(result);
  }
  
  SEXP updatePredictor(SEXP fitExpr, SEXP xExpr, SEXP colsExpr, SEXP forceUpdateExpr, SEXP updateCutPointsExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_updatePredictor called on NULL external pointer");
    
    if (!Rf_isReal(xExpr)) Rf_error("x must be of type real");
    
    bool forceUpdate     = rc_getBool(forceUpdateExpr,         "forceUpdate", RC_NA | RC_NO, RC_END);
    bool updateCutPoints = rc_getBool(updateCutPointsExpr, "updateCutPoints", RC_NA | RC_NO, RC_END);
    
    bool result;
    if (Rf_isNull(colsExpr)) {
      rc_assertDimConstraints(xExpr, "dimensions of x", RC_LENGTH | RC_EQ, rc_asRLength(2),
                              RC_VALUE | RC_EQ, static_cast<int>(fit->data.numObservations),
                              RC_VALUE | RC_EQ, static_cast<int>(fit->data.numPredictors),
                              RC_END);
      
      result = fit->setPredictor(REAL(xExpr), forceUpdate, updateCutPoints);
    } else {
      
      if (!Rf_isInteger(colsExpr)) Rf_error("columns must be of type integer");
      
      SEXP dimsExpr = Rf_getAttrib(xExpr, R_DimSymbol);
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
        if (rc_getLength(xExpr) != fit->data.numObservations) Rf_error("length of new x does not match y");
      }
      
      int* colsInt = INTEGER(colsExpr);
      size_t numCols = rc_getLength(colsExpr);
      size_t* cols = misc_stackAllocate(numCols, size_t);
      for (size_t i = 0 ; i < numCols; ++i) {
        cols[i] = static_cast<size_t>(colsInt[i] - 1);
        if (cols[i] >= fit->data.numPredictors) {
          misc_stackFree(cols);
          Rf_error("column '%d' is out of range", colsInt[i] + 1);
        }
      }
      
      result = fit->updatePredictor(REAL(xExpr), cols, numCols, forceUpdate, updateCutPoints);
      
      misc_stackFree(cols);
    }
    
    return Rf_ScalarLogical(result);
  }
  
  SEXP setCutPoints(SEXP fitExpr, SEXP cutPointsExpr, SEXP colsExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_setCutPoints called on NULL external pointer");
    
    if (!Rf_isNewList(cutPointsExpr)) Rf_error("cutPoints must be of type list");
 
    
    size_t numCols;
    if (Rf_isNull(colsExpr)) {
      numCols = fit->data.numPredictors;
    } else {
      if (!Rf_isInteger(colsExpr)) Rf_error("columns must be of type integer");
      numCols = rc_getLength(colsExpr);
    }
    
    if (rc_getLength(cutPointsExpr) != numCols)
      Rf_error("length of cutPoints (%lu) must equal length of columns (%lu)", rc_getLength(cutPointsExpr), numCols);
    
    const double** cutPoints = misc_stackAllocate(numCols, const double*);
    uint32_t* numCutPoints = misc_stackAllocate(numCols, uint32_t);
    size_t* cols = misc_stackAllocate(numCols, size_t);
    
    int* colsInt = Rf_isNull(colsExpr) ? NULL : INTEGER(colsExpr);
    for (size_t i = 0; i < numCols; ++i) {
      SEXP cutPointsExpr_i = VECTOR_ELT(cutPointsExpr, i);
      
      cutPoints[i] = REAL(cutPointsExpr_i);
      numCutPoints[i] = static_cast<uint32_t>(rc_getLength(cutPointsExpr_i));
      cols[i] = colsInt == NULL ? i : static_cast<size_t>(colsInt[i] - 1);
       
      if (cols[i] >= fit->data.numPredictors) {
        misc_stackFree(cols);
        misc_stackFree(numCutPoints);
        misc_stackFree(cutPoints);
        Rf_error("column '%d' is out of range", colsInt[i] + 1);
      }
    }
    
    fit->setCutPoints(cutPoints, numCutPoints, cols, numCols);
    
    misc_stackFree(cols);
    misc_stackFree(numCutPoints);
    misc_stackFree(cutPoints);
        
    return R_NilValue;
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
    size_t* cols = misc_stackAllocate(numCols, size_t);
    for (size_t i = 0 ; i < numCols; ++i) {
      cols[i] = static_cast<size_t>(colsInt[i] - 1);
      if (cols[i] >= fit->data.numPredictors) {
        misc_stackFree(cols);
        Rf_error("column '%d' is out of range", colsInt[i]);
      }
    }
    
    fit->updateTestPredictors(REAL(x_test), cols, numCols);
    
    misc_stackFree(cols);
    
    return R_NilValue;
  }
  
  SEXP storeLatents(SEXP fitExpr, SEXP resultExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_storeLatents called on NULL external pointer");
    
    if (!fit->control.responseIsBinary)
      Rf_error("dbarts_storeLatents called on sampler with non-binary response");
    
    if (Rf_isNull(resultExpr)) {
      resultExpr = PROTECT(rc_newReal(fit->data.numObservations * fit->control.numChains));
    
      fit->storeLatents(REAL(resultExpr));
      
      UNPROTECT(1);
    } else {
      if (rc_getLength(resultExpr) < fit->data.numObservations * fit->control.numChains)
        Rf_error("dbarts_storeLatents called with vector of insufficient length");
      
      fit->storeLatents(REAL(resultExpr));
    }
    
    return resultExpr;
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
    
    initializeStateFromExpression(*fit, stateExpr);
    
    return R_NilValue;
  }
  
  SEXP storeState(SEXP fitExpr, SEXP stateExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) Rf_error("dbarts_storeState called on NULL external pointer");
    
    storeStateExpressionFromFit(*fit, stateExpr);
    
    return R_NilValue;
  }
  
  
  SEXP printTrees(SEXP fitExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr)
  {
    BARTFit* fitPtr = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fitPtr == NULL) Rf_error("dbarts_printTrees called on NULL external pointer");
    BARTFit& fit(*fitPtr);
    
    size_t numChains  = fit.control.numChains;
    size_t numSamples = fit.control.keepTrees ? fit.currentNumSamples : 0;
    size_t numTrees   = fit.control.numTrees;
    
    size_t numChainIndices  = Rf_isNull(chainIndicesExpr)  ? numChains  : rc_getLength(chainIndicesExpr);
    size_t numSampleIndices = Rf_isNull(sampleIndicesExpr) ? numSamples : rc_getLength(sampleIndicesExpr);
    size_t numTreeIndices   = Rf_isNull(treeIndicesExpr)   ? numTrees   : rc_getLength(treeIndicesExpr);
    
    if (numChainIndices > numChains)
      Rf_error("%lu chains specified but only %lu in sampler", numChainIndices, numChains);
    if (numSampleIndices > numSamples)
      Rf_error("%lu samples specified but only %lu in sampler", numSampleIndices, numSamples);
    if (numTreeIndices > numTrees)
      Rf_error("%lu trees specified but only %lu in sampler", numTreeIndices, numTrees);    
    
    size_t* chainIndices  = misc_stackAllocate(numChainIndices, size_t);
    size_t* sampleIndices = fit.control.keepTrees ? new size_t[numSamples] : NULL;
    size_t* treeIndices   = new size_t[numTreeIndices];
    
    if (Rf_isNull(chainIndicesExpr)) {
      for (size_t i = 0; i < numChains; ++i) chainIndices[i] = i;
    } else {
      int* i_chainIndices = INTEGER(chainIndicesExpr);
      for (size_t i = 0; i < numChainIndices; ++i) chainIndices[i] = static_cast<size_t>(i_chainIndices[i] - 1);
    }
    
    if (Rf_isNull(sampleIndicesExpr)) {
      for (size_t i = 0; i < numSamples; ++i) sampleIndices[i] = i;
    } else {
      int* i_sampleIndices = INTEGER(sampleIndicesExpr);
      for (size_t i = 0; i < numSampleIndices; ++i) sampleIndices[i] = static_cast<size_t>(i_sampleIndices[i] - 1);
    }
    
    if (Rf_isNull(treeIndicesExpr)) {
      for (size_t i = 0; i < numTrees; ++i) treeIndices[i] = i;
    } else {
      int* i_treeIndices = INTEGER(treeIndicesExpr);
      for (size_t i = 0; i < numTreeIndices; ++i) treeIndices[i] = static_cast<size_t>(i_treeIndices[i] - 1);
    }
   
    fit.printTrees(chainIndices, numChainIndices, sampleIndices, numSampleIndices, treeIndices, numTreeIndices);
    
    delete [] treeIndices;
    delete [] sampleIndices;
    misc_stackFree(chainIndices);
    
    return R_NilValue;
  }
  
  SEXP getTrees(SEXP fitExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr)
  {
    BARTFit* fitPtr = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fitPtr == NULL) Rf_error("dbarts_getTrees called on NULL external pointer");
    BARTFit& fit(*fitPtr);
    
    size_t numChains  = fit.control.numChains;
    size_t numSamples = fit.control.keepTrees ? fit.currentNumSamples : 0;
    size_t numTrees   = fit.control.numTrees;
    
    size_t numChainIndices  = Rf_isNull(chainIndicesExpr)  ? numChains  : rc_getLength(chainIndicesExpr);
    size_t numSampleIndices = Rf_isNull(sampleIndicesExpr) ? numSamples : rc_getLength(sampleIndicesExpr);
    size_t numTreeIndices   = Rf_isNull(treeIndicesExpr)   ? numTrees   : rc_getLength(treeIndicesExpr);
    
    if (numChainIndices > numChains)
      Rf_error("%lu chains specified but only %lu in sampler", numChainIndices, numChains);
    if (numSampleIndices > numSamples)
      Rf_error("%lu samples specified but only %lu in sampler", numSampleIndices, numSamples);
    if (numTreeIndices > numTrees)
      Rf_error("%lu trees specified but only %lu in sampler", numTreeIndices, numTrees);
    
    size_t* chainIndices  = misc_stackAllocate(numChainIndices, size_t);
    size_t* sampleIndices = fit.control.keepTrees ? new size_t[numSamples] : NULL;
    size_t* treeIndices   = new size_t[numTreeIndices];
    
    if (Rf_isNull(chainIndicesExpr)) {
      for (size_t i = 0; i < numChains; ++i) chainIndices[i] = i;
    } else {
      int* i_chainIndices = INTEGER(chainIndicesExpr);
      for (size_t i = 0; i < numChainIndices; ++i) chainIndices[i] = static_cast<size_t>(i_chainIndices[i] - 1);
    }
    
    if (Rf_isNull(sampleIndicesExpr)) {
      for (size_t i = 0; i < numSamples; ++i) sampleIndices[i] = i;
    } else {
      int* i_sampleIndices = INTEGER(sampleIndicesExpr);
      for (size_t i = 0; i < numSampleIndices; ++i) sampleIndices[i] = static_cast<size_t>(i_sampleIndices[i] - 1);
    }
    
    if (Rf_isNull(treeIndicesExpr)) {
      for (size_t i = 0; i < numTrees; ++i) treeIndices[i] = i;
    } else {
      int* i_treeIndices = INTEGER(treeIndicesExpr);
      for (size_t i = 0; i < numTreeIndices; ++i) treeIndices[i] = static_cast<size_t>(i_treeIndices[i] - 1);
    }
    
    FlattenedTrees* flattenedTreesPtr =
      fit.getFlattenedTrees(chainIndices, numChainIndices,
                            sampleIndices, numSampleIndices,
                            treeIndices, numTreeIndices);
    FlattenedTrees& flattenedTrees(*flattenedTreesPtr);
    
    delete [] treeIndices;
    delete [] sampleIndices;
    misc_stackFree(chainIndices);
    
    R_xlen_t numCols = 4 + (numChains > 1 ? 1 : 0) + (fit.control.keepTrees ? 1 : 0);
    SEXP resultExpr = PROTECT(rc_newList(numCols));
        
    SEXP classExpr = PROTECT(rc_newCharacter(1));
    SET_STRING_ELT(classExpr, 0, Rf_mkChar("data.frame"));
    Rf_setAttrib(resultExpr, R_ClassSymbol, classExpr);
    UNPROTECT(1);
    
    SEXP resultRowNamesExpr;
    rc_allocateInSlot2(resultRowNamesExpr, resultExpr, R_RowNamesSymbol, STRSXP, flattenedTrees.totalNumNodes);
    
    SEXP resultNamesExpr;
    rc_allocateInSlot2(resultNamesExpr, resultExpr, R_NamesSymbol, STRSXP, numCols);
    
    int* chainNumber = NULL;
    int* sampleNumber = NULL;
    int* treeNumber, *numObservations, *variable;
    double* value;
        
    R_xlen_t colNum = 0;
    if (numChains > 1) {
      SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
      SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("chain")));
      UNPROTECT(2);
      chainNumber = INTEGER(VECTOR_ELT(resultExpr, colNum));
      ++colNum;
    }
    if (fit.control.keepTrees) {
      SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
      SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("sample")));
      UNPROTECT(2);
      sampleNumber = INTEGER(VECTOR_ELT(resultExpr, colNum));
      ++colNum;
    }
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("tree")));
    treeNumber = INTEGER(VECTOR_ELT(resultExpr, colNum));
    ++colNum;
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("n")));
    numObservations = INTEGER(VECTOR_ELT(resultExpr, colNum));
    ++colNum;
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("var")));
    variable = INTEGER(VECTOR_ELT(resultExpr, colNum));
    ++colNum;
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newReal(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("value")));
    value = REAL(VECTOR_ELT(resultExpr, colNum));
    UNPROTECT(8);
    
    size_t numDigits = 1;
    size_t temp = flattenedTrees.totalNumNodes;
    while (temp >= 10) {
      temp /= 10;
      ++numDigits;
    }
    char* buffer = new char[numDigits + 1];
    for (size_t i = 0; i < flattenedTrees.totalNumNodes; ++i) {
      if (chainNumber != NULL)
        chainNumber[i] = static_cast<int>(flattenedTrees.chainNumber[i] + 1);
      if (sampleNumber != NULL)
        sampleNumber[i] = static_cast<int>(flattenedTrees.sampleNumber[i] + 1);
      treeNumber[i] = static_cast<int>(flattenedTrees.treeNumber[i] + 1);
      numObservations[i] = static_cast<int>(flattenedTrees.numObservations[i]);
      int variable_i = static_cast<int>(flattenedTrees.variable[i]);
      variable[i] = variable_i >= 0 ? variable_i + 1 : variable_i;
      value[i] = flattenedTrees.value[i];
#if __cplusplus < 201112L
#  if defined(_WIN64) || SIZEOF_SIZE_T == 8
      snprintf(buffer, numDigits + 1, "%lu", static_cast<unsigned long>(i + 1));
#  else
      snprintf(buffer, numDigits + 1, "%u", static_cast<unsigned int>(i + 1));
#  endif
#else
      snprintf(buffer, numDigits + 1, "%zu", i + 1);
#endif
      SET_STRING_ELT(resultRowNamesExpr, i, PROTECT(Rf_mkChar(buffer)));
      UNPROTECT(1);
    }
    
    delete [] buffer;
    delete flattenedTreesPtr;
    
    UNPROTECT(1);
    
    return resultExpr;
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

