#include "config.hpp"
#include <cstring>
#include <cstddef>
#include <dbarts/cstdint.hpp>
#include <cmath>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>
#include <dbarts/types.hpp>
#include <dbarts/R_C_interface.hpp>

#include "makeModelMatrixFromDataFrame.h"

#include <external/alloca.h>
#include <external/linearAlgebra.h>
#include <external/random.h>

#include <set>
#ifdef THREAD_SAFE_UNLOAD
#include <pthread.h>
#endif

using std::size_t;

extern "C" {
  void R_init_dbarts(DllInfo* info);
}

namespace {
  using namespace dbarts;
  
  void initializeControlFromExpression(Control& control, SEXP controlExpr);
  // SEXP createControlExpressionFromFit(const BARTFit& fit);
  void initializeModelFromExpression(Model& model, SEXP modelExpr, const Control& control);
  // SEXP createModelExpressionFromFit(const BARTFit& fit);
  void initializeDataFromExpression(Data& data, SEXP dataExpr);
  // SEXP createDataExpressionFromFit(const BARTFit& fit);
  
  void initializeStateFromExpression(const BARTFit& fit, State& state, SEXP stateExpr);
  SEXP createStateExpressionFromFit(const BARTFit& fit); // result has a protect count of 1
  void storeStateExpressionFromFit(const BARTFit& fit, SEXP stateExpr);
  
  SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, R_xlen_t length);
  SEXP SET_DIMS(SEXP obj, int numRows, int numCols);
  bool isS4Null(SEXP expr);
  
  void deleteFit(BARTFit* fit);
  
  struct ExternalPointerComparator {
    bool operator()(const SEXP& lhs, const SEXP& rhs) const {
      return R_ExternalPtrAddr(const_cast<SEXP>(lhs)) < R_ExternalPtrAddr(const_cast<SEXP>(rhs));
    }
  };
  typedef std::set<SEXP, ExternalPointerComparator> PointerSet;
  PointerSet* activeFits;
#ifdef THREAD_SAFE_UNLOAD
  pthread_mutex_t fitMutex;
#endif

}

extern "C" {  
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

namespace {
  using namespace dbarts;
  
/*  static SEXP simulateContinuousUniform(SEXP nExpr)
  {
    size_t n = 0;
    if (LENGTH(nExpr) > 0) n = (size_t) INTEGER(nExpr)[0];
    
    SEXP seedsExpr = findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
    if (seedsExpr == R_UnboundValue) GetRNGstate();
    if (TYPEOF(seedsExpr) == PROMSXP) seedsExpr = eval(R_SeedsSymbol, R_GlobalEnv);
    
    // uint_least32_t seed0 = (uint_least32_t) INTEGER(seedsExpr)[0];
    
    // ext_rng_algorithm_t algorithmType = (ext_rng_algorithm_t) (seed0 % 100);
    // ext_rng_standardNormal_t stdNormalType = (ext_rng_standardNormal_t) (seed0 / 100);
    
    ext_rng_userFunction uniformFunction;
    uniformFunction.f.stateless = &unif_rand;
    uniformFunction.state = NULL;
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_USER_UNIFORM, &uniformFunction);
        
    SEXP resultExpr = PROTECT(allocVector(REALSXP, n));
    double* result = REAL(resultExpr);
    for (size_t i = 0; i < n; ++i) result[i] = ext_rng_simulateContinuousUniform(rng);
    
    ext_rng_destroy(rng);
    
    UNPROTECT(1);
    return resultExpr;
  }
  
  static SEXP simulateNormal(SEXP nExpr)
  {
    size_t n = 0;
    if (LENGTH(nExpr) > 0) n = (size_t) INTEGER(nExpr)[0];
    
    SEXP seedsExpr = findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
    if (seedsExpr == R_UnboundValue) GetRNGstate();
    if (TYPEOF(seedsExpr) == PROMSXP) seedsExpr = eval(R_SeedsSymbol, R_GlobalEnv);
    
    uint_least32_t seed0 = (uint_least32_t) INTEGER(seedsExpr)[0];
    
    // ext_rng_algorithm_t algorithmType = (ext_rng_algorithm_t) (seed0 % 100);
    ext_rng_standardNormal_t stdNormalType = (ext_rng_standardNormal_t) (seed0 / 100);
    
    ext_rng_userFunction uniformFunction;
    uniformFunction.f.stateless = &unif_rand;
    uniformFunction.state = NULL;
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_USER_UNIFORM, &uniformFunction);
    ext_rng_setStandardNormalAlgorithm(rng, stdNormalType, NULL);
    
    SEXP resultExpr = PROTECT(allocVector(REALSXP, n));
    double* result = REAL(resultExpr);
    for (size_t i = 0; i < n; ++i) result[i] = ext_rng_simulateStandardNormal(rng);
    
    ext_rng_destroy(rng);
    
    UNPROTECT(1);
    return resultExpr;
  }
  
  static SEXP simulateExponential(SEXP nExpr)
  {
    size_t n = 0;
    if (LENGTH(nExpr) > 0) n = (size_t) INTEGER(nExpr)[0];
    
    SEXP seedsExpr = findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
    if (seedsExpr == R_UnboundValue) GetRNGstate();
    if (TYPEOF(seedsExpr) == PROMSXP) seedsExpr = eval(R_SeedsSymbol, R_GlobalEnv);
    
    uint_least32_t seed0 = (uint_least32_t) INTEGER(seedsExpr)[0];
    
    // ext_rng_algorithm_t algorithmType = (ext_rng_algorithm_t) (seed0 % 100);
    ext_rng_standardNormal_t stdNormalType = (ext_rng_standardNormal_t) (seed0 / 100);
    
    ext_rng_userFunction uniformFunction;
    uniformFunction.f.stateless = &unif_rand;
    uniformFunction.state = NULL;
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_USER_UNIFORM, &uniformFunction);
    ext_rng_setStandardNormalAlgorithm(rng, stdNormalType, NULL);
    
    SEXP resultExpr = PROTECT(allocVector(REALSXP, n));
    double* result = REAL(resultExpr);
    for (size_t i = 0; i < n; ++i) result[i] = ext_rng_simulateExponential(rng, 1.0);
    
    ext_rng_destroy(rng);
    
    UNPROTECT(1);
    return resultExpr;
  }
  
  static ext_rng* createRNG()
  {
    SEXP seedsExpr = findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
    if (seedsExpr == R_UnboundValue) GetRNGstate();
    if (TYPEOF(seedsExpr) == PROMSXP) seedsExpr = eval(R_SeedsSymbol, R_GlobalEnv);
    
    uint_least32_t seed0 = (uint_least32_t) INTEGER(seedsExpr)[0];
    
    ext_rng_algorithm_t algorithmType = (ext_rng_algorithm_t) (seed0 % 100);
    ext_rng_standardNormal_t stdNormalType = (ext_rng_standardNormal_t) (seed0 / 100);
    
    void* state = (void*) (1 + INTEGER(seedsExpr));
    switch (algorithmType) {
      case EXT_RNG_ALGORITHM_KNUTH_TAOCP:
      case EXT_RNG_ALGORITHM_KNUTH_TAOCP2:
      {
        ext_rng_knuthState* kt = (ext_rng_knuthState*) malloc(sizeof(ext_rng_knuthState));
        memcpy(kt->state1, state, EXT_RNG_KNUTH_NUM_RANDOM * sizeof(uint_least32_t));
        kt->info = EXT_RNG_KNUTH_NUM_RANDOM; // this is a static var which we cannot access
        for (size_t i = 0; i < EXT_RNG_KNUTH_QUALITY; ++i) kt->state2[i] = 0; // also static
        state = kt;
      }
      break;
      case EXT_RNG_ALGORITHM_USER_UNIFORM:
      {
        ext_rng_userFunction* uniformFunction = (ext_rng_userFunction*) malloc(sizeof(ext_rng_userFunction));
        uniformFunction->f.stateless = &unif_rand;
        uniformFunction->state = NULL;
        state = uniformFunction;
      }
      break;
      default:
      break;
    }
    
    ext_rng* rng = ext_rng_create(algorithmType, state);
    if (algorithmType == EXT_RNG_ALGORITHM_KNUTH_TAOCP || algorithmType == EXT_RNG_ALGORITHM_KNUTH_TAOCP2 || algorithmType == EXT_RNG_ALGORITHM_USER_UNIFORM) free(state);
    if (rng == NULL) return NULL; 
    
    void* normalState = NULL;
    switch (stdNormalType) {
      case EXT_RNG_STANDARD_NORMAL_BOX_MULLER:
      normalState = malloc(sizeof(double));
      *((double*) normalState) = 0.0; // static var, again
      break;
      case EXT_RNG_STANDARD_NORMAL_USER_NORM:
      {
        ext_rng_userFunction* normalFunction = (ext_rng_userFunction*) malloc(sizeof(ext_rng_userFunction));
        normalFunction->f.stateless = &norm_rand;
        normalFunction->state = NULL;
        normalState = normalFunction;
      }
      break;
      default:
      break;
    }
    
    int errorCode = ext_rng_setStandardNormalAlgorithm(rng, stdNormalType, normalState);
    if (stdNormalType == EXT_RNG_STANDARD_NORMAL_BOX_MULLER || stdNormalType == EXT_RNG_STANDARD_NORMAL_USER_NORM) free(normalState);
    
    if (errorCode != 0) {
      ext_rng_destroy(rng);
      return NULL;
    }
    
    return rng;
  }
  
  // takes type of generator from R and uses internal implementation
  static SEXP simulateContinuousUniformInternally(SEXP nExpr)
  {
    size_t n = 0;
    if (LENGTH(nExpr) > 0) n = (size_t) INTEGER(nExpr)[0];
    
    ext_rng* rng = createRNG();
    
    if (rng == NULL) return NULL_USER_OBJECT;
    
    SEXP resultExpr = PROTECT(allocVector(REALSXP, n));
    double* result = REAL(resultExpr);
    for (size_t i = 0; i < n; ++i) result[i] = ext_rng_simulateContinuousUniform(rng);
    
    ext_rng_destroy(rng);
    
    
    UNPROTECT(1);
    return resultExpr;
  }
  
  static SEXP simulateNormalInternally(SEXP nExpr)
   {
    size_t n = 0;
    if (LENGTH(nExpr) > 0) n = (size_t) INTEGER(nExpr)[0];
    
    ext_rng* rng = createRNG();
    
    if (rng == NULL) return NULL_USER_OBJECT;
    
    SEXP resultExpr = PROTECT(allocVector(REALSXP, n));
    double* result = REAL(resultExpr);
    for (size_t i = 0; i < n; ++i) result[i] = ext_rng_simulateStandardNormal(rng);
    
    ext_rng_destroy(rng);
    
    
    UNPROTECT(1);
    return resultExpr;
  } */
  
  SEXP saveToFile(SEXP fitExpr, SEXP fileName)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_saveToFile called on NULL external pointer");
    
    return ScalarLogical(fit->saveToFile(CHAR(STRING_ELT(fileName, 0))));
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
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setResponse(SEXP fitExpr, SEXP y)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setY called on NULL external pointer");
    
    if (!isReal(y)) error("y must be of type real");
    if (static_cast<size_t>(length(y)) != fit->data.numObservations) error("length of new y does not match old");
    
    if (fit->control.responseIsBinary) GetRNGstate();
        
    fit->setResponse(REAL(y));
    
    if (fit->control.responseIsBinary) PutRNGstate();
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setOffset(SEXP fitExpr, SEXP offsetExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setOffset called on NULL external pointer");
    
    double* offset = NULL;
    if (isReal(offsetExpr)) {
      offset = REAL(offsetExpr);
      if (static_cast<size_t>(length(offsetExpr)) != fit->data.numObservations) error("length of new offset does not match y");
    } else if (!isNull(offsetExpr) && !isS4Null(offsetExpr)) {
      error("offset must be of type real or NULL");
    }
    
    if (fit->control.responseIsBinary) GetRNGstate();
    
    fit->setOffset(offset);
    
    if (fit->control.responseIsBinary) PutRNGstate();
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setPredictor(SEXP fitExpr, SEXP x)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setPredictor called on NULL external pointer");
    
    if (!isReal(x)) error("x must be of type real");
    
    SEXP dimsExpr = GET_DIM(x);
    
    if (isNull(dimsExpr) || LENGTH(dimsExpr) != 2) error("x must be a matrix, i.e. have two dimensions");
    int* dims = INTEGER(dimsExpr);
    
    if (static_cast<size_t>(dims[0]) != fit->data.numObservations) error("number of rows in new x does not match y");
    if (static_cast<size_t>(dims[1]) != fit->data.numPredictors) error("number of columns in new x does not match old");
    
    return ScalarLogical(fit->setPredictor(REAL(x)));
  }
  
  SEXP updatePredictor(SEXP fitExpr, SEXP x, SEXP colsExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_updatePredictor called on NULL external pointer");
    
    if (!isReal(x)) error("x must be of type real");
    if (!isInteger(colsExpr)) error("columns must be of type integer");
    
    SEXP dimsExpr = GET_DIM(x);
    int* dims = NULL;
    
    if (!isNull(dimsExpr)) {
      int numDims = GET_LENGTH(dimsExpr);
      
      if (numDims != 1 && numDims != 2) error("x must be a vector or a matrix");
      if (numDims == 2) dims = INTEGER(dimsExpr);
    }
    
    if (length(colsExpr) == 0) error("length of columns is 0");

    if (dims != NULL) {
      if (static_cast<size_t>(dims[0]) != fit->data.numObservations) error("number of rows of new x does not match y");
      if (dims[1] != length(colsExpr)) error("number of columns of new x does not match length of columns to replace");
    } else {
      if (static_cast<size_t>(length(x)) != fit->data.numObservations) error("length of new x does not match y");
    }
    
    
    int* colsInt = INTEGER(colsExpr);
    size_t numCols = static_cast<size_t>(LENGTH(colsExpr));
    size_t* cols = ext_stackAllocate(numCols, size_t);
    for (size_t i = 0 ; i < numCols; ++i) {
      cols[i] = static_cast<size_t>(colsInt[i] - 1);
      if (static_cast<size_t>(cols[i]) >= fit->data.numPredictors) {
        ext_stackFree(cols);
        error("column '%d' is out of range", colsInt[i]);
      }
    }
    
    bool result = fit->updatePredictors(REAL(x), cols, numCols);
    
    ext_stackFree(cols);
    
    return ScalarLogical(result);
  }
  
  SEXP setTestPredictor(SEXP fitExpr, SEXP x_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setTestPredictor called on NULL external pointer");
    
    if (isNull(x_test) || isS4Null(x_test)) {
      fit->setTestPredictor(NULL, 0);
    
      return NULL_USER_OBJECT;
    }
    
    if (!isReal(x_test)) error("x.test must be of type real");
    SEXP dimsExpr = GET_DIM(x_test);
    if (GET_LENGTH(dimsExpr) != 2) error("x.test must be a matrix, i.e. have two dimensions");
    int* dims = INTEGER(dimsExpr);
    if (static_cast<size_t>(dims[1]) != fit->data.numPredictors) error("number of columns of x.test and x must be equal");
    
    fit->setTestPredictor(REAL(x_test), static_cast<size_t>(dims[0]));
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setTestOffset(SEXP fitExpr, SEXP offset_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setTestOffset called on NULL external pointer");
    
    if (isNull(offset_test)) {
      fit->setTestOffset(NULL);
    } else {
      if (!isReal(offset_test)) error("offset.test must be of type real");
      if (fit->data.numTestObservations != static_cast<size_t>(GET_LENGTH(offset_test))) error("length of offset.test must equal number of rows in x.test");
      fit->setTestOffset(REAL(offset_test));
    }
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setTestPredictorAndOffset(SEXP fitExpr, SEXP x_test, SEXP offset_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setTestPredictorAndOffset called on NULL external pointer");
    
    if (isNull(x_test) || isS4Null(x_test)) {
      fit->setTestPredictor(NULL, 0);
    
      return NULL_USER_OBJECT;
    }
    
    if (!isReal(x_test)) error("x.test must be of type real");
    SEXP dimsExpr = GET_DIM(x_test);
    if (GET_LENGTH(dimsExpr) != 2) error("x.test must be a matrix, i.e. have two dimensions");
    int* dims = INTEGER(dimsExpr);
    if (static_cast<size_t>(dims[1]) != fit->data.numPredictors) error("number of columns of x.test and x must be equal");
    
    if (isNull(offset_test)) {
      fit->setTestPredictorAndOffset(REAL(x_test), NULL, static_cast<size_t>(dims[0]));
    } else {
      if (!isReal(offset_test)) error("offset.test must be of type real");
      if (GET_LENGTH(offset_test) == 1 && ISNA(REAL(offset_test)[0])) {
        fit->setTestPredictor(REAL(x_test), static_cast<size_t>(dims[0]));
      } else {
        if (GET_LENGTH(offset_test) != dims[0]) error("length of offset.test must equal number of rows in x.test");
        fit->setTestPredictorAndOffset(REAL(x_test), REAL(offset_test), static_cast<size_t>(dims[0]));
      }
    }
    
    return NULL_USER_OBJECT;
  }
  
  SEXP updateTestPredictor(SEXP fitExpr, SEXP x_test, SEXP colsExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_updateTestPredictor called on NULL external pointer");
    
    if (fit->data.X_test == NULL) error("test matrix must exist at object creation to be updated");
    
    if (!isReal(x_test)) error("x must be of type real");
    if (!isInteger(colsExpr)) error("columns must be of type integer");
    
    SEXP dimsExpr = GET_DIM(x_test);
    int* dims = NULL;
    
    if (!isNull(dimsExpr)) {
      int numDims = GET_LENGTH(dimsExpr);
      
      if (numDims != 1 && numDims != 2) error("x must be a vector or a matrix");
      if (numDims == 2) dims = INTEGER(dimsExpr);
    }
    
    if (length(colsExpr) == 0) error("length of columns is 0");

    if (dims != NULL) {
      if (static_cast<size_t>(dims[0]) != fit->data.numTestObservations) error("number of rows of new x does not match old x.test");
      if (dims[1] != length(colsExpr)) error("number of columns of new x does not match length of columns to replace");
    } else {
      if (static_cast<size_t>(length(x_test)) != fit->data.numTestObservations) error("length of new x does not match old x.test");
    }
    
    
    int* colsInt = INTEGER(colsExpr);
    size_t numCols = static_cast<size_t>(LENGTH(colsExpr));
    size_t* cols = ext_stackAllocate(numCols, size_t);
    for (size_t i = 0 ; i < numCols; ++i) {
      cols[i] = static_cast<size_t>(colsInt[i] - 1);
      if (cols[i] >= fit->data.numPredictors) {
        ext_stackFree(cols);
        error("column '%d' is out of range", colsInt[i]);
      }
    }
    
    fit->updateTestPredictors(REAL(x_test), cols, numCols);
    
    ext_stackFree(cols);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setControl(SEXP fitExpr, SEXP controlExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setControl called on NULL external pointer");
    
    if (strcmp(CHAR(STRING_ELT(GET_CLASS(controlExpr), 0)), "dbartsControl") != 0) error("'control' argument to dbarts_create not of class 'dbartsControl'");
    
    initializeControlFromExpression(fit->control, controlExpr);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP isValidPointer(SEXP fitExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) return ScalarLogical(FALSE);
    
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_lock(&fitMutex);
#endif
    if (activeFits->find(fitExpr) != activeFits->end()) {
#ifdef THREAD_SAFE_UNLOAD
      pthread_mutex_unlock(&fitMutex);
#endif
      return ScalarLogical(TRUE);
    }
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_unlock(&fitMutex);
#endif
    return ScalarLogical(FALSE);
  }
  
  SEXP create(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr)
  {
    Control control;
    Model model;
    Data data;
    
    SEXP classExpr = GET_CLASS(controlExpr);
    if (strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsControl") != 0) error("'control' argument to dbarts_create not of class 'dbartsControl'");
    
    classExpr = GET_CLASS(modelExpr);
    if (strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsModel") != 0) error("'model' argument to dbarts_create not of class 'dbartsModel'");
    
    classExpr = GET_CLASS(dataExpr);
    if (strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsData") != 0) error("'data' argument to dbarts_create not of class 'dbartsData'");
    
    
    initializeControlFromExpression(control, controlExpr);
    initializeModelFromExpression(model, modelExpr, control);
    initializeDataFromExpression(data, dataExpr);
    
    BARTFit* fit = new BARTFit(control, model, data);
    
    SEXP result = PROTECT(R_MakeExternalPtr(fit, NULL_USER_OBJECT, NULL_USER_OBJECT));
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
  
  SEXP createState(SEXP fitExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_createState called on NULL external pointer");
    
    SEXP result = createStateExpressionFromFit(*fit);
    
    UNPROTECT(1);
    
    return result;
  }
  
  SEXP restoreState(SEXP fitExpr, SEXP stateExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_restoreState called on NULL external pointer");
    
    initializeStateFromExpression(*fit, fit->state, stateExpr);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP storeState(SEXP fitExpr, SEXP stateExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_storeState called on NULL external pointer");
    
    storeStateExpressionFromFit(*fit, stateExpr);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP run(SEXP fitExpr, SEXP numBurnInExpr, SEXP numSamplesExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_run called on NULL external pointer");
    
    int i_temp;
    size_t numBurnIn, numSamples;
    
    if (!isInteger(numBurnInExpr)) error("number of burn-in steps must be of integer type");
    if (length(numBurnInExpr) == 0) error("number of burn-in steps must be of length at least 1");
    i_temp = INTEGER(numBurnInExpr)[0];
    if (i_temp != NA_INTEGER && i_temp < 0) error("number of burn-in steps must be non-negative");
    numBurnIn = i_temp == NA_INTEGER ? fit->control.numBurnIn : static_cast<size_t>(i_temp);
    
    if (!isInteger(numSamplesExpr)) error("number of samples must be of integer type");
    if (length(numSamplesExpr) == 0) error("number of samples must be of length at least 1");
    i_temp = INTEGER(numSamplesExpr)[0];
    if (i_temp != NA_INTEGER && i_temp < 0) error("number of samples must be non-negative");
    numSamples = i_temp == NA_INTEGER ? fit->control.numSamples : static_cast<size_t>(i_temp);
    
    if (numBurnIn == 0 && numSamples == 0) error("either number of burn-in or samples must be positive");
    
    size_t numTrainingSamples = fit->data.numObservations * numSamples;
    if (numSamples != 0 && numTrainingSamples / numSamples != fit->data.numObservations)
      error("training sample array size exceeds machine's capacity");
    R_xlen_t s_numTrainingSamples = static_cast<R_xlen_t>(numTrainingSamples);
    if (s_numTrainingSamples < 0 || static_cast<size_t>(s_numTrainingSamples) != numTrainingSamples)
      error("training sample array size cannot be represented by a signed integer on this machine");
    
    size_t numTestSamples = fit->data.numTestObservations * numSamples;
     if (numSamples != 0 && numTestSamples / numSamples != fit->data.numTestObservations)
      error("test sample array size exceeds machine's capacity");
    R_xlen_t s_numTestSamples = static_cast<R_xlen_t>(numTestSamples);
    if (s_numTestSamples < 0 || static_cast<size_t>(s_numTestSamples) != numTestSamples)
      error("test sample array size cannot be represented by a signed integer on this machine");
    
    GetRNGstate();
        
    Results* bartResults = fit->runSampler(numBurnIn, numSamples);
    
    PutRNGstate();
    
    // can happen if numSamples == 0
    if (bartResults == NULL) return NULL_USER_OBJECT;
    
    
    // create result storage and make it user friendly
    SEXP namesExpr;
    
    int protectCount = 0;
    
    SEXP resultExpr = PROTECT(allocVector(VECSXP, 4));
    ++protectCount;
    SET_VECTOR_ELT(resultExpr, 0, allocVector(REALSXP, static_cast<R_xlen_t>(bartResults->getNumSigmaSamples())));
    SET_VECTOR_ELT(resultExpr, 1, allocVector(REALSXP, static_cast<R_xlen_t>(bartResults->getNumTrainingSamples())));
    if (fit->data.numTestObservations > 0)
      SET_VECTOR_ELT(resultExpr, 2, allocVector(REALSXP, static_cast<R_xlen_t>(bartResults->getNumTestSamples())));
    else
      SET_VECTOR_ELT(resultExpr, 2, NULL_USER_OBJECT);
    SET_VECTOR_ELT(resultExpr, 3, allocVector(INTSXP, static_cast<R_xlen_t>(bartResults->getNumVariableCountSamples())));
    
    SEXP sigmaSamples = VECTOR_ELT(resultExpr, 0);
    std::memcpy(REAL(sigmaSamples), const_cast<const double*>(bartResults->sigmaSamples), bartResults->getNumSigmaSamples() * sizeof(double));
    
    SEXP trainingSamples = VECTOR_ELT(resultExpr, 1);
    SET_DIMS(trainingSamples, bartResults->numObservations, bartResults->numSamples);
    std::memcpy(REAL(trainingSamples), const_cast<const double*>(bartResults->trainingSamples), bartResults->getNumTrainingSamples() * sizeof(double));
    
    if (fit->data.numTestObservations > 0) {
      SEXP testSamples = VECTOR_ELT(resultExpr, 2);
      SET_DIMS(testSamples, bartResults->numTestObservations, bartResults->numSamples);
      std::memcpy(REAL(testSamples), const_cast<const double*>(bartResults->testSamples), bartResults->getNumTestSamples() * sizeof(double));
    }
    
    SEXP variableCountSamples = VECTOR_ELT(resultExpr, 3);
    SET_DIMS(variableCountSamples, bartResults->numPredictors, bartResults->numSamples);
    int* variableCountStorage = INTEGER(variableCountSamples);
    size_t length = bartResults->getNumVariableCountSamples();
    // these likely need to be down-sized from 64 to 32 bits
    for (size_t i = 0; i < length; ++i) variableCountStorage[i] = static_cast<int>(bartResults->variableCountSamples[i]);
    
    setAttrib(resultExpr, R_NamesSymbol, namesExpr = allocVector(STRSXP, 4));
    SET_STRING_ELT(namesExpr, 0, mkChar("sigma"));
    SET_STRING_ELT(namesExpr, 1, mkChar("train"));
    SET_STRING_ELT(namesExpr, 2, mkChar("test"));
    SET_STRING_ELT(namesExpr, 3, mkChar("varcount"));
    
    UNPROTECT(protectCount);
    
    delete bartResults;
    
    return resultExpr;
  }
  
  SEXP finalize(void) {
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_lock(&fitMutex);
#endif
    for (PointerSet::iterator it = activeFits->begin(); it != activeFits->end(); ) {
      SEXP fitExpr = *it;
      BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
#ifdef THREAD_SAFE_UNLOAD
      Rprintf("package finalizing %p\n", fit);
#endif
      
      deleteFit(fit);
      PointerSet::iterator prev = it;
      ++it;
      activeFits->erase(prev);
      R_ClearExternalPtr(fitExpr);
    }
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_unlock(&fitMutex);
    pthread_mutex_destroy(&fitMutex);
#endif
    
    delete activeFits;
    
    return NULL_USER_OBJECT;
  }
  
  SEXP deepCopy(SEXP obj)
  {
    return duplicate(obj);
  }
    
  // as of R 3.1, auto-unload never gets called so screw that
  
/*  void R_unload_dbarts(DllInfo* info)
  {
    pthread_mutex_lock(&fitMutex);
    for (PointerSet::iterator it = activeFits->begin(); it != activeFits->end(); ) {
      BARTFit* fit = *it;
      deleteFit(fit);
      PointerSet::iterator prev = it;
      ++it;
      activeFits->erase(prev);
    }
    pthread_mutex_unlock(&fitMutex);
    pthread_mutex_destroy(&fitMutex);
  }*/
  
#define DEF_FUNC(_N_, _F_, _A_) { _N_, reinterpret_cast<DL_FUNC>(&_F_), _A_ }

  R_CallMethodDef R_callMethods[] = {
    DEF_FUNC("dbarts_create", create, 3),
    DEF_FUNC("dbarts_run", run, 3),
    DEF_FUNC("dbarts_setResponse", setResponse, 2),
    DEF_FUNC("dbarts_setOffset", setOffset, 2),
    DEF_FUNC("dbarts_setPredictor", setPredictor, 2),
    DEF_FUNC("dbarts_updatePredictor", updatePredictor, 3),
    DEF_FUNC("dbarts_setTestPredictor", setTestPredictor, 2),
    DEF_FUNC("dbarts_setTestOffset", setTestOffset, 2),
    DEF_FUNC("dbarts_setTestPredictorAndOffset", setTestPredictorAndOffset, 3),
    DEF_FUNC("dbarts_updateTestPredictor", updateTestPredictor, 3),
    DEF_FUNC("dbarts_setControl", setControl, 2),
    DEF_FUNC("dbarts_isValidPointer", isValidPointer, 1),
    DEF_FUNC("dbarts_createState", createState, 1),
    DEF_FUNC("dbarts_storeState", storeState, 2),
    DEF_FUNC("dbarts_restoreState", restoreState, 2),
    DEF_FUNC("dbarts_finalize", finalize, 0),
    DEF_FUNC("dbarts_deepCopy", deepCopy, 1),
    // experimental
    DEF_FUNC("dbarts_saveToFile", saveToFile, 2),
    DEF_FUNC("dbarts_loadFromFile", loadFromFile, 1),
    // below: testing
//    DEF_FUNC("dbarts_runif", simulateContinuousUniform, 1),
//    DEF_FUNC("dbarts_runif", simulateContinuousUniformInternally, 1),
//    DEF_FUNC("dbarts_rnorm", simulateNormal, 1),
//    DEF_FUNC("dbarts_rnorm", simulateNormalInternally, 1),
//    DEF_FUNC("dbarts_rexp", simulateExponential, 1),
    DEF_FUNC("dbarts_makeModelMatrixFromDataFrame", makeModelMatrixFromDataFrame, 2),
    { NULL, NULL, 0 }
  };

#undef DEF_FUNC
  
  struct C_CallMethodDef {
    const char* name;
    DL_FUNC function;
  };
  
#define DEF_FUNC(_N_, _F_) { _N_, reinterpret_cast<DL_FUNC>(&_F_) }
  
  C_CallMethodDef C_callMethods[] = {
    DEF_FUNC("createCGMPrior", dbarts_createCGMPrior),
    DEF_FUNC("createCGMPriorFromOptions", dbarts_createCGMPriorFromOptions),
    DEF_FUNC("destroyCGMPrior", dbarts_destroyCGMPrior),
    DEF_FUNC("initializeCGMPriorFromOptions", dbarts_initializeCGMPriorFromOptions),
    DEF_FUNC("invalidateCGMPrior", dbarts_invalidateCGMPrior),
    
    DEF_FUNC("createNormalPrior", dbarts_createNormalPrior),
    DEF_FUNC("createNormalPriorFromOptions", dbarts_createNormalPriorFromOptions),
    DEF_FUNC("destroyNormalPrior", dbarts_destroyNormalPrior),
    DEF_FUNC("initializeNormalPriorFromOptions", dbarts_initializeNormalPriorFromOptions),
    DEF_FUNC("invalidateNormalPrior", dbarts_invalidateNormalPrior),
    
    DEF_FUNC("createChiSquaredPrior", dbarts_createChiSquaredPrior),
    DEF_FUNC("createChiSquaredPriorFromOptions", dbarts_createChiSquaredPriorFromOptions),
    DEF_FUNC("destroyChiSquaredPrior", dbarts_destroyChiSquaredPrior),
    DEF_FUNC("initializeChiSquaredPriorFromOptions", dbarts_initializeChiSquaredPriorFromOptions),
    DEF_FUNC("invalidateChiSquaredPrior", dbarts_invalidateChiSquaredPrior),

    DEF_FUNC("createFit", dbarts_createFit),
    DEF_FUNC("initializeFit", dbarts_initializeFit),
    DEF_FUNC("destroyFit", dbarts_destroyFit),
    DEF_FUNC("invalidateFit", dbarts_invalidateFit),
    
    DEF_FUNC("runSampler", dbarts_runSampler),
    DEF_FUNC("runSamplerForIterations", dbarts_runSamplerForIterations),
    DEF_FUNC("setResponse", dbarts_setResponse),
    DEF_FUNC("setOffset", dbarts_setOffset),
    DEF_FUNC("setPredictor", dbarts_setPredictor),
    DEF_FUNC("updatePredictor", dbarts_updatePredictor),
    DEF_FUNC("updatePredictors", dbarts_updatePredictors),
    DEF_FUNC("setTestPredictor", dbarts_setTestPredictor),
    DEF_FUNC("setTestOffset", dbarts_setTestOffset),
    DEF_FUNC("setTestPredictorsAndOffset", dbarts_setTestPredictorAndOffset),
    DEF_FUNC("updateTestPredictor", dbarts_updateTestPredictor),
    DEF_FUNC("updateTestPredictors", dbarts_updateTestPredictors),
    { NULL, 0 }
  };
  
#undef DEF_FUNC
  
} // end anonymous namespace

extern "C" {
  void R_init_dbarts(DllInfo* info)
  {
    R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
    R_useDynamicSymbols(info, static_cast<Rboolean>(FALSE));
    
    C_CallMethodDef* method = C_callMethods;
    while (method->name != NULL) {
      R_RegisterCCallable("dbarts", method->name, method->function);
      ++method;
    }
    
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_init(&fitMutex, NULL);
#endif
    
    activeFits = new PointerSet;
  }
}

namespace {
  using namespace dbarts;
  
  SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, R_xlen_t length)
  {
    SEXP val = allocVector(type, length);
    
    SET_SLOT(obj, nm, val);
    return val;
  }
  
  SEXP SET_DIMS(SEXP obj, int numRows, int numCols)
  {
    SEXP dimsExp = NEW_INTEGER(2);
    int* dims = INTEGER(dimsExp);
    dims[0] = numRows;
    dims[1] = numCols;
    
    SET_ATTR(obj, R_DimSymbol, dimsExp);
    
    return obj;
  }
  
  bool isS4Null(SEXP expr)
  {
    if (!isSymbol(expr)) return false;
    
    const char* symbolName = CHAR(PRINTNAME(expr));
    
    if (strncmp(symbolName, "\1NULL\1", 6) == 0) return true;
    
    return false;
  }
  
  void initializeControlFromExpression(Control& control, SEXP controlExpr)
  {
    int i_temp;
    
    SEXP slotExpr = GET_ATTR(controlExpr, install("binary"));
    if (!isLogical(slotExpr)) error("binary response must be signified by logical type");
    if (length(slotExpr) != 1) error("binary response signifier must be of length 1");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("binary response must be either true or false");
    control.responseIsBinary = (i_temp != FALSE);
    
    slotExpr = GET_ATTR(controlExpr, install("verbose"));
    if (!isLogical(slotExpr)) error("verbose must be signified by logical type");
    if (length(slotExpr) == 0) error("verbose must be of length at least 1");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("verbose must be either true or false");
    control.verbose = (i_temp != FALSE);
    
    slotExpr = GET_ATTR(controlExpr, install("keepTrainingFits"));
    if (!isLogical(slotExpr)) error("keep training fits must be signified by logical type");
    if (length(slotExpr) != 1) error("keep training fits must be of length 1");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("keep training fits must be either true or false");
    control.keepTrainingFits = (i_temp != FALSE);
    
    slotExpr = GET_ATTR(controlExpr, install("useQuantiles"));
    if (!isLogical(slotExpr)) error("use quantiles must be signified by logical type");
    if (length(slotExpr) != 1) error("use quantiles must be of length 1");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("use quantiles must be either true or false");
    control.useQuantiles = (i_temp != FALSE);
    
    
    
    slotExpr = GET_ATTR(controlExpr, install("n.samples"));
    if (!isInteger(slotExpr)) error("number of samples must be of integer type");
    if (length(slotExpr) != 1) error("number of samples must be of length 1");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) error("number of samples cannot be NA");
    if (i_temp < 0) error("number of samples must be non-negative");
    control.numSamples = static_cast<size_t>(i_temp);
    
    slotExpr = GET_ATTR(controlExpr, install("n.burn"));
    if (!isInteger(slotExpr)) error("number of burn-in steps must be of integer type");
    if (length(slotExpr) != 1) error("number of burn-in steps must be of length 1");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 0;
    if (i_temp < 0) error("number of burn-in steps must be non-negative");
    control.numBurnIn = static_cast<size_t>(i_temp);
    
    slotExpr = GET_ATTR(controlExpr, install("n.trees"));
    if (!isInteger(slotExpr)) error("number of trees must be of integer type");
    if (length(slotExpr) != 1) error("number of trees must be of length 1");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) error("number of trees cannot be NA");
    if (i_temp <= 0) error("number of trees must be positive");
    control.numTrees = static_cast<size_t>(i_temp);
    
    slotExpr = GET_ATTR(controlExpr, install("n.threads"));
    if (!isInteger(slotExpr)) error("number of threads must be of integer type");
    if (length(slotExpr) != 1) error("number of threads must be of length 1");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 1;
    if (i_temp <= 0) error("number of threads must be positive");
    control.numThreads = static_cast<size_t>(i_temp);
    
    slotExpr = GET_ATTR(controlExpr, install("n.thin"));
    if (!isInteger(slotExpr)) error("tree thinning rate must be of integer type");
    if (length(slotExpr) != 1) error("tree thinning rate must be of length 1");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 1;
    if (i_temp < 0) error("tree thinning rate must be non-negative");
    control.treeThinningRate = static_cast<uint32_t>(i_temp);
    
    
    slotExpr = GET_ATTR(controlExpr, install("printEvery"));
    if (!isInteger(slotExpr)) error("print every must be of integer type");
    if (length(slotExpr) != 1) error("print every must be of length 1");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp != NA_INTEGER) {
      if (i_temp <= 0) error("print every must be positive");
      control.printEvery = static_cast<uint32_t>(i_temp);
    }
    
    slotExpr = GET_ATTR(controlExpr, install("printCutoffs"));
    if (!isInteger(slotExpr)) error("print cutoffs must be of integer type");
    if (length(slotExpr) == 0) error("print cutoffs must be of length at least 1");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 0;
    if (i_temp < 0) error("print cutoffs must be non-negative");
    control.printCutoffs = static_cast<uint32_t>(i_temp);
    
    if (control.rng == NULL) {
      ext_rng_userFunction uniformFunction;
      uniformFunction.f.stateless = &unif_rand;
      uniformFunction.state = NULL;
      control.rng = ext_rng_create(EXT_RNG_ALGORITHM_USER_UNIFORM, &uniformFunction);
      
      ext_rng_userFunction normalFunction;
      normalFunction.f.stateless = &norm_rand;
      normalFunction.state = NULL;
      ext_rng_setStandardNormalAlgorithm(control.rng, EXT_RNG_STANDARD_NORMAL_USER_NORM, &normalFunction);
    }
  }
  
  void initializeModelFromExpression(Model& model, SEXP modelExpr, const Control& control)
  {
    double d_temp;
    
    SEXP slotExpr = GET_ATTR(modelExpr, install("p.birth_death"));
    if (!isReal(slotExpr)) error("probability of birth/death rule must be of numeric type");
    if (length(slotExpr) != 1) error("probability of birth/death rule must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("probability of birth/death rule must be a real number");
    if (d_temp <= 0.0 || d_temp > 1.0) error("probability of birth/death rule must be in (0, 1]");
    model.birthOrDeathProbability = d_temp;
    
    slotExpr = GET_ATTR(modelExpr, install("p.swap"));
    if (!isReal(slotExpr)) error("probability of swap rule must be of numeric type");
    if (length(slotExpr) != 1) error("probability of swap rule must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("probability of swap rule must be a real number");
    if (d_temp < 0.0 || d_temp >= 1.0) error("probability of swap rule must be in [0, 1)");
    model.swapProbability = d_temp;
    
    slotExpr = GET_ATTR(modelExpr, install("p.change"));
    if (!isReal(slotExpr)) error("probability of change rule must be of numeric type");
    if (length(slotExpr) != 1) error("probability of change rule must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("probability of change rule must be a real number");
    if (d_temp < 0.0 || d_temp >= 1.0) error("probability of change rule must be in [0, 1)");
    model.changeProbability = d_temp;
    
    if (std::fabs(model.birthOrDeathProbability + model.swapProbability + model.changeProbability - 1.0) >= 1.0e-10)
      error("rule proposal probabilities must sum to 1.0");
    
    slotExpr = GET_ATTR(modelExpr, install("p.birth"));
    if (!isReal(slotExpr)) error("probability of birth in birth/death rule must be of numeric type");
    if (length(slotExpr) != 1) error("probability of birth in birth/death rule must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("probability of birth in birth/death rule must be a real number");
    if (d_temp <= 0.0 || d_temp >= 1.0) error("probability of birth in birth/death rule must be in (0, 1)");
    model.birthProbability = d_temp;
    
    
    SEXP priorExpr = GET_ATTR(modelExpr, install("tree.prior"));
    // slotExpr = GET_CLASS(priorExpr);
    // if (strcmp(CHAR(STRING_ELT(GET_CLASS(slotExpr), 0)), "dbartsControl") != 0) error("'control' argument to dbarts_create not of class 'dbartsControl'");
    CGMPrior* treePrior = new CGMPrior;
    model.treePrior = treePrior;
    
    slotExpr = GET_ATTR(priorExpr, install("power"));
    if (!isReal(slotExpr)) error("tree prior power must be of type real");
    if (length(slotExpr) != 1) error("tree prior power must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("tree prior power be a real number");
    if (d_temp <= 0.0) error("tree prior power must be positive");
    treePrior->power = d_temp;
    
    slotExpr = GET_ATTR(priorExpr, install("base"));
    if (!isReal(slotExpr)) error("tree prior base must be of type real");
    if (length(slotExpr) != 1) error("tree prior power must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("tree prior base be a real number");
    if (d_temp <= 0.0 || d_temp >= 1.0) error("tree prior base must be in (0, 1)");
    treePrior->base = d_temp;
    
    
    priorExpr = GET_ATTR(modelExpr, install("node.prior"));
    
    slotExpr = GET_ATTR(priorExpr, install("k"));
    if (!isReal(slotExpr)) error ("k must be of type real");
    if (length(slotExpr) != 1) error("k must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("k must be a real number");
    if (d_temp <= 0.0) error("k must be positive");
    model.muPrior = new NormalPrior(control, d_temp);
    
    
    
    priorExpr = GET_ATTR(modelExpr, install("resid.prior"));
    
    slotExpr = GET_ATTR(priorExpr, install("df"));
    if (!isReal(slotExpr)) error("sigma prior degrees of freedom must be of type real");
    if (length(slotExpr) != 1) error("sigma prior degrees of freedom must be of length 1");
    double sigmaPriorDf = REAL(slotExpr)[0];
    if (ISNAN(sigmaPriorDf)) error("sigma prior degrees of freedom must be a real number");
    if (sigmaPriorDf <= 0.0) error("sigma prior degrees of freedom must be positive");
    
    slotExpr = GET_ATTR(priorExpr, install("quantile"));
    if (!isReal(slotExpr)) error ("sigma prior quantile must be of type real");
    if (length(slotExpr) != 1) error("sigma prior quantile must be of length 1");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("sigma prior quantile must be a real number");
    if (d_temp <= 0.0 || d_temp >= 1.0) error("sigma prior quantile must be in (0, 1)");
    model.sigmaSqPrior = new ChiSquaredPrior(sigmaPriorDf, d_temp);
  }

  void initializeDataFromExpression(Data& data, SEXP dataExpr)
  {
    int* dims;
    
    SEXP slotExpr = GET_ATTR(dataExpr, install("y"));
    if (!isReal(slotExpr)) error("y must be of type real");
    if (length(slotExpr) == 0) error("length of y must be greater than 0");
    data.y = REAL(slotExpr);
    data.numObservations = static_cast<size_t>(length(slotExpr));
    
    slotExpr = GET_ATTR(dataExpr, install("x"));
    if (!isReal(slotExpr)) error("x must be of type real");
    dims = INTEGER(GET_ATTR(slotExpr, R_DimSymbol));
    if (dims == NULL || length(GET_ATTR(slotExpr, R_DimSymbol)) != 2) error("x must be a matrix, i.e. have two dimensions");
    if (static_cast<size_t>(dims[0]) != data.numObservations) error("number of rows of x and length of y must be equal");
    data.X = REAL(slotExpr);
    data.numPredictors = static_cast<size_t>(dims[1]);
    
    slotExpr = GET_ATTR(dataExpr, install("varTypes"));
    if (!isInteger(slotExpr)) error("variable types must be of type integer");
    if (static_cast<size_t>(length(slotExpr)) != data.numPredictors) error("length of variable types must equal number of columns in x");
    int* i_variableTypes = INTEGER(slotExpr);
    VariableType* variableTypes = new VariableType[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) variableTypes[i] = (i_variableTypes[i] == 0 ? ORDINAL : CATEGORICAL);
    data.variableTypes = variableTypes;
    
    slotExpr = GET_ATTR(dataExpr, install("x.test"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.X_test = NULL;
      data.numTestObservations = 0;
    } else {
      if (!isReal(slotExpr)) error ("x.test must be of type real");
      dims = INTEGER(GET_ATTR(slotExpr, R_DimSymbol));
      if (dims == NULL || length(GET_ATTR(slotExpr, R_DimSymbol)) != 2) error("x.test must be a matrix, i.e. have two dimensions");
      if (static_cast<size_t>(dims[1]) != data.numPredictors) error("number of columns of x.test and x must be equal");
      data.X_test = REAL(slotExpr);
      data.numTestObservations = static_cast<size_t>(dims[0]);
    }
    
    slotExpr = GET_ATTR(dataExpr, install("weights"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.weights = NULL;
    } else {
      if (!isReal(slotExpr)) error("weights must be of type real");
      if (static_cast<size_t>(length(slotExpr)) != data.numObservations) error("length of weights must equal length of y");
      data.weights = REAL(slotExpr);
    }
    
    slotExpr = GET_ATTR(dataExpr, install("offset"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.offset = NULL;
    } else {
      if (!isReal(slotExpr)) error("offset must be of type real");
      if (static_cast<size_t>(length(slotExpr)) != data.numObservations) error("length of offset must equal length of y");
      data.offset = REAL(slotExpr);
    }
    
    slotExpr = GET_ATTR(dataExpr, install("offset.test"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.testOffset = NULL;
    } else {
      if (!isReal(slotExpr)) error("test offset must be of type real");
      if (static_cast<size_t>(length(slotExpr)) != data.numTestObservations) error("length of test offset must equal number of test rows");
      data.testOffset = REAL(slotExpr);
    }
    
    slotExpr = GET_ATTR(dataExpr, install("sigma"));
    if (!isReal(slotExpr)) error("sigma estimate must be of type real");
    if (length(slotExpr) != 1) error("sigma estimate must be of length 1");
    double d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) d_temp = 1.0;
    if (d_temp <= 0.0) error("sigma estimate must be positive");
    data.sigmaEstimate = d_temp;
    
    
    slotExpr = GET_ATTR(dataExpr, install("n.cuts"));
    if (!isInteger(slotExpr)) error("maximum number of cuts must be of integer type");
    if (static_cast<size_t>(length(slotExpr)) != data.numPredictors) error("length of maximum number of cuts and the number of columns of x must be equal");
    int* i_maxNumCuts = INTEGER(slotExpr);
    uint32_t* maxNumCuts = new uint32_t[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) maxNumCuts[i] = static_cast<uint32_t>(i_maxNumCuts[i]);
    data.maxNumCuts = maxNumCuts;
  }
  
  SEXP createStateExpressionFromFit(const BARTFit& fit)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const State& state(fit.state);
    
    SEXP result = PROTECT(NEW_OBJECT(MAKE_CLASS("dbartsState")));
    
    SEXP slotExpr = ALLOC_SLOT(result, install("fit.tree"), REALSXP, static_cast<R_xlen_t>(data.numObservations * control.numTrees));
    SET_DIMS(slotExpr, static_cast<int>(data.numObservations), static_cast<int>(control.numTrees));
    std::memcpy(REAL(slotExpr), state.treeFits, data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = ALLOC_SLOT(result, install("fit.total"), REALSXP, static_cast<R_xlen_t>(data.numObservations));
    std::memcpy(REAL(slotExpr), state.totalFits, data.numObservations * sizeof(double));
    
    if (data.numTestObservations == 0) {
      SET_SLOT(result, install("fit.test"), NULL_USER_OBJECT);
    } else {
      slotExpr = ALLOC_SLOT(result, install("fit.test"), REALSXP, static_cast<R_xlen_t>(data.numTestObservations));
      std::memcpy(REAL(slotExpr), state.totalTestFits, data.numTestObservations * sizeof(double));
    }
    
    slotExpr = ALLOC_SLOT(result, install("sigma"), REALSXP, 1);
    REAL(slotExpr)[0] = state.sigma;
    
    slotExpr = ALLOC_SLOT(result, install("runningTime"), REALSXP, 1);
    REAL(slotExpr)[0] = state.runningTime;
    
    slotExpr = ALLOC_SLOT(result, install("trees"), STRSXP, static_cast<R_xlen_t>(control.numTrees));

    const char** treeStrings = const_cast<const char**>(state.createTreeStrings(fit));
    for (size_t i = 0; i < control.numTrees; ++i) {
      SET_STRING_ELT(slotExpr, static_cast<R_xlen_t>(i), CREATE_STRING_VECTOR(treeStrings[i]));
      delete [] treeStrings[i];
    }
    delete [] treeStrings;
    
    return result;
  }
  
  void storeStateExpressionFromFit(const BARTFit& fit, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const State& state(fit.state);
    
    SEXP slotExpr = GET_ATTR(stateExpr, install("fit.tree"));
    SEXP dimsExpr = GET_DIM(slotExpr);
    if (GET_LENGTH(dimsExpr) != 2) error("dimensions of state@fit.tree indicate that it is not a matrix");
    int* dims = INTEGER(dimsExpr);
    if (static_cast<size_t>(dims[0]) != data.numObservations || static_cast<size_t>(dims[1]) != control.numTrees) error("dimensions of state@fit.tree do not match object");
    std::memcpy(REAL(slotExpr), state.treeFits, data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = GET_ATTR(stateExpr, install("fit.total"));
    if (static_cast<size_t>(GET_LENGTH(slotExpr)) != data.numObservations) error("length of state@fit.total does not match object");
    std::memcpy(REAL(slotExpr), state.totalFits, data.numObservations * sizeof(double));
    
    if (data.numTestObservations != 0) {
      slotExpr = GET_ATTR(stateExpr, install("fit.test"));
      if (static_cast<size_t>(GET_LENGTH(slotExpr)) != data.numTestObservations) error("length of state@fit.test does not match object");
      std::memcpy(REAL(slotExpr), state.totalTestFits, data.numTestObservations * sizeof(double));
    }
    
    slotExpr = GET_ATTR(stateExpr, install("sigma"));
    if (GET_LENGTH(slotExpr) != 1) error("length of state@sigma does not match object");
    REAL(slotExpr)[0] = state.sigma;
    
    slotExpr = GET_ATTR(stateExpr, install("runningTime"));
    if (GET_LENGTH(slotExpr) != 1) error("length of state@runningTime does not match object");
    REAL(slotExpr)[0] = state.runningTime;
    
    slotExpr = GET_ATTR(stateExpr, install("trees"));
    if (static_cast<size_t>(GET_LENGTH(slotExpr)) != control.numTrees) error("length of state@trees does not match object");
    
    const char** treeStrings = const_cast<const char**>(state.createTreeStrings(fit));
    for (size_t i = 0; i < control.numTrees; ++i) {
      SET_STRING_ELT(slotExpr, static_cast<R_xlen_t>(i), CREATE_STRING_VECTOR(treeStrings[i]));
      delete [] treeStrings[i];
    }
    delete [] treeStrings;
  }
  
  void initializeStateFromExpression(const BARTFit& fit, State& state, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    
    SEXP slotExpr = GET_ATTR(stateExpr, install("fit.tree"));
    std::memcpy(state.treeFits, const_cast<const double*>(REAL(slotExpr)), data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = GET_ATTR(stateExpr, install("fit.total"));
    std::memcpy(state.totalFits, const_cast<const double*>(REAL(slotExpr)), data.numObservations * sizeof(double));
    
    if (data.numTestObservations != 0) {
      slotExpr = GET_ATTR(stateExpr, install("fit.test"));
      std::memcpy(state.totalTestFits, const_cast<const double*>(REAL(slotExpr)), data.numTestObservations * sizeof(double));
    }
    
    slotExpr = GET_ATTR(stateExpr, install("sigma"));
    state.sigma = REAL(slotExpr)[0];
    
    slotExpr = GET_ATTR(stateExpr, install("runningTime"));
    state.runningTime = REAL(slotExpr)[0];
    
    slotExpr = GET_ATTR(stateExpr, install("trees"));
    const char** treeStrings = ext_stackAllocate(control.numTrees, const char*);
    for (size_t i = 0; i < control.numTrees; ++i) {
      treeStrings[i] = CHAR(STRING_ELT(slotExpr, static_cast<R_xlen_t>(i)));
    }
    state.recreateTreesFromStrings(fit, treeStrings);
    
    ext_stackFree(treeStrings);
  }
  
  void deleteFit(BARTFit* fit) {
#ifdef THREAD_SAFE_UNLOAD
    Rprintf("deleting   %p\n", fit);
#endif
    if (fit == NULL) return;
    
    ext_rng_destroy(fit->control.rng);
    
    delete fit->model.treePrior;
    delete fit->model.muPrior;
    delete fit->model.sigmaSqPrior;
    
    delete [] fit->data.variableTypes;
    delete [] fit->data.maxNumCuts;
    
    delete fit;
  }
}
