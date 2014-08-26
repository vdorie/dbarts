#include "config.hpp"
#include <cstring>
#include <cstddef>
#include <dbarts/cstdint.hpp>
#include <cmath>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>
#include <dbarts/types.hpp>
#include <dbarts/R_C_interface.hpp>

#include <external/alloca.h>
#include <external/linearAlgebra.h>

#include <set>
#ifdef THREAD_SAFE_UNLOAD
#include <pthread.h>
#endif

using std::size_t;

namespace {
  using namespace dbarts;
  
  void initializeControlFromExpression(Control& control, SEXP controlExpr);
  void initializeModelFromExpression(Model& model, SEXP modelExpr, const Control& control);
  void initializeDataFromExpression(Data& data, SEXP dataExpr);
  
  void initializeStateFromExpression(const BARTFit& fit, State& state, SEXP stateExpr);
  SEXP createStateExpressionFromFit(const BARTFit& fit); // result has a protect count of 1
  void storeStateExpressionFromFit(const BARTFit& fit, SEXP stateExpr);
  
  SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length);
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
  SEXP setResponse(SEXP fitExpr, SEXP y)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setY called on NULL external pointer.");
    
    if (!isReal(y)) error("y must be of type real.");
    if ((size_t) length(y) != fit->data.numObservations) error("Length of new y does not match old.");
    fit->setResponse(REAL(y));
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setOffset(SEXP fitExpr, SEXP offsetExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setOffset called on NULL external pointer.");
    
    double* offset = NULL;
    if (isReal(offsetExpr)) {
      offset = REAL(offsetExpr);
      if ((size_t) length(offsetExpr) != fit->data.numObservations) error("Length of new offset does not match y.");
    } else if (!isNull(offsetExpr) && !isS4Null(offsetExpr)) {
      error("offset must be of type real or NULL.");
    }
    fit->setOffset(offset);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setPredictor(SEXP fitExpr, SEXP x, SEXP jExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setPredictor called on NULL external pointer.");
    
    if (!isReal(x)) error("x must be of type real.");
    if ((size_t) length(x) != fit->data.numObservations) error("Length of new x does not match y.");
    if (!isInteger(jExpr)) error("Column must be of type integer.");
    if (length(jExpr) == 0) error("Length of column is 0.");
    
    int j = INTEGER(jExpr)[0] - 1;
    
    if (j < 0 || j >= (int) fit->data.numPredictors) error("Column is out of range.");
    
    fit->setPredictor(REAL(x), (size_t) j);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setTestPredictor(SEXP fitExpr, SEXP x_test, SEXP jExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setTestPredictor called on NULL external pointer.");
    
    if (fit->data.X_test == NULL) error("Test matrix must exist at object creation to be updated.");
    
    if (!isReal(x_test)) error("x must be of type real.");
    if ((size_t) length(x_test) != fit->data.numTestObservations) error("Length of new x_test does not match old.");
    if (!isInteger(jExpr)) error("Column must be of type integer.");
    if (length(jExpr) == 0) error("Length of column is 0.");
    
    int j = INTEGER(jExpr)[0] - 1;
    
    if (j < 0 || j >= (int) fit->data.numPredictors) error("Column is out of range.");
    
    fit->setTestPredictor(REAL(x_test), (size_t) j);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setTestPredictors(SEXP fitExpr, SEXP x_test, SEXP offset_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setTestPredictors called on NULL external pointer.");
    
    if (isNull(x_test) || isS4Null(x_test)) {
      fit->setTestPredictors(NULL, 0);
    
      return NULL_USER_OBJECT;
    }
    
    if (!isReal(x_test)) error("x.test must be of type real.");
    SEXP dimsExpr = GET_DIM(x_test);
    if (GET_LENGTH(dimsExpr) != 2) error("x.test must be a matrix, i.e. have two dimensions.");
    int* dims = INTEGER(dimsExpr);
    if ((size_t) dims[1] != fit->data.numPredictors) error("Number of columns of x.test and x must be equal.");
    
    if (isNull(offset_test)) {
      fit->setTestPredictors(REAL(x_test), NULL, (size_t) dims[0]);
    } else {
      if (!isReal(offset_test)) error("offset.test must be of type real");
      if (GET_LENGTH(offset_test) == 1 && ISNA(REAL(offset_test)[0])) {
        fit->setTestPredictors(REAL(x_test), (size_t) dims[0]);
      } else {
        if (GET_LENGTH(offset_test) != dims[0]) error("length of offset.test must equal number of rows in x.test");
        fit->setTestPredictors(REAL(x_test), REAL(offset_test), (size_t) dims[0]);
      }
    }
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setTestOffset(SEXP fitExpr, SEXP offset_test)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setTestOffset called on NULL external pointer.");
    
    if (isNull(offset_test)) {
      fit->setTestOffset(NULL);
    } else {
      if (!isReal(offset_test)) error("offset.test must be of type real");
      if (fit->data.numTestObservations != (size_t) GET_LENGTH(offset_test)) error("length of offset.test must equal number of rows in x.test");
      fit->setTestOffset(REAL(offset_test));
    }
    
    return NULL_USER_OBJECT;
  }
  
  SEXP setControl(SEXP fitExpr, SEXP controlExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_setControl called on NULL external pointer.");
    
    if (strcmp(CHAR(STRING_ELT(GET_CLASS(controlExpr), 0)), "dbartsControl") != 0) error("'control' argument to dbarts_create not of class 'dbartsControl'.");
    
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
    if (strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsControl") != 0) error("'control' argument to dbarts_create not of class 'dbartsControl'.");
    
    classExpr = GET_CLASS(modelExpr);
    if (strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsModel") != 0) error("'model' argument to dbarts_create not of class 'dbartsModel'.");
    
    classExpr = GET_CLASS(dataExpr);
    if (strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsData") != 0) error("'data' argument to dbarts_create not of class 'dbartsData'.");
    
    
    initializeControlFromExpression(control, controlExpr);
    initializeModelFromExpression(model, modelExpr, control);
    initializeDataFromExpression(data, dataExpr);
    
    BARTFit* fit = new BARTFit(control, model, data);
    
    SEXP result = PROTECT(R_MakeExternalPtr(fit, NULL_USER_OBJECT, NULL_USER_OBJECT));
    R_RegisterCFinalizerEx(result, fitFinalizer, (Rboolean) FALSE);
    
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
    if (fit == NULL) error("dbarts_createState called on NULL external pointer.");
    
    SEXP result = createStateExpressionFromFit(*fit);
    
    UNPROTECT(1);
    
    return result;
  }
  
  SEXP restoreState(SEXP fitExpr, SEXP stateExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_restoreState called on NULL external pointer.");
    
    initializeStateFromExpression(*fit, fit->state, stateExpr);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP storeState(SEXP fitExpr, SEXP stateExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_storeState called on NULL external pointer.");
    
    storeStateExpressionFromFit(*fit, stateExpr);
    
    return NULL_USER_OBJECT;
  }
  
  SEXP run(SEXP fitExpr, SEXP numBurnInExpr, SEXP numSamplesExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("dbarts_run called on NULL external pointer.");
    
    int i_temp;
    size_t numBurnIn, numSamples;
    
    if (!isInteger(numBurnInExpr)) error("Number of burn-in steps must be of integer type.");
    if (length(numBurnInExpr) == 0) error("Number of burn-in steps must be of length at least 1.");
    i_temp = INTEGER(numBurnInExpr)[0];
    if (i_temp != NA_INTEGER && i_temp < 0) error("Number of burn-in steps must be non-negative.");
    numBurnIn = i_temp == NA_INTEGER ? fit->control.numBurnIn : (size_t) i_temp;
    
    if (!isInteger(numSamplesExpr)) error("Number of samples must be of integer type.");
    if (length(numSamplesExpr) == 0) error("Number of samples must be of length at least 1.");
    i_temp = INTEGER(numSamplesExpr)[0];
    if (i_temp != NA_INTEGER && i_temp <= 0) error("Number of samples must be positive.");
    numSamples = i_temp == NA_INTEGER ? fit->control.numSamples : (size_t) i_temp;
    
    GetRNGstate();
        
    Results* bartResults = fit->runSampler(numBurnIn, numSamples);
    
    PutRNGstate();
    
    
    // create result storage and make it user friendly
    SEXP dimsExpr, namesExpr;
    int* dims;
    
    int protectCount = 3;
    
    SEXP resultExpr = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(resultExpr, 0, allocVector(REALSXP, (int) bartResults->getNumSigmaSamples()));
    SET_VECTOR_ELT(resultExpr, 1, allocVector(REALSXP, (int) bartResults->getNumTrainingSamples()));
    if (fit->data.numTestObservations > 0)
      SET_VECTOR_ELT(resultExpr, 2, allocVector(REALSXP, (int) bartResults->getNumTestSamples()));
    else
      SET_VECTOR_ELT(resultExpr, 2, NULL_USER_OBJECT);
    SET_VECTOR_ELT(resultExpr, 3, allocVector(INTSXP, (int) bartResults->getNumVariableCountSamples()));
    
    SEXP sigmaSamples = VECTOR_ELT(resultExpr, 0);
    std::memcpy(REAL(sigmaSamples), (const double*) bartResults->sigmaSamples, bartResults->getNumSigmaSamples() * sizeof(double));
    
    SEXP trainingSamples = VECTOR_ELT(resultExpr, 1);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
    dims = INTEGER(dimsExpr);
    dims[0] = (int) bartResults->numObservations;
    dims[1] = (int) bartResults->numSamples;
    setAttrib(trainingSamples, R_DimSymbol, dimsExpr);
    std::memcpy(REAL(trainingSamples), (const double*) bartResults->trainingSamples, bartResults->getNumTrainingSamples() * sizeof(double));
    
    if (fit->data.numTestObservations > 0) {
      SEXP testSamples = VECTOR_ELT(resultExpr, 2);
      dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
      dims = INTEGER(dimsExpr);
      dims[0] = (int) bartResults->numTestObservations;
      dims[1] = (int) bartResults->numSamples;
      setAttrib(testSamples, R_DimSymbol, dimsExpr);
      std::memcpy(REAL(testSamples), (const double*) bartResults->testSamples, bartResults->getNumTestSamples() * sizeof(double));
      ++protectCount;
    }
    
    SEXP variableCountSamples = VECTOR_ELT(resultExpr, 3);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
    dims = INTEGER(dimsExpr);
    dims[0] = (int) bartResults->numPredictors;
    dims[1] = (int) bartResults->numSamples;
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
    
    UNPROTECT(protectCount);
    
    delete bartResults;
    
    return(resultExpr);
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

  R_CallMethodDef R_callMethods[] = {
    { "dbarts_create", (DL_FUNC) &create, 3 },
    { "dbarts_run", (DL_FUNC) &run, 3 },
    { "dbarts_setResponse", (DL_FUNC) &setResponse, 2 },
    { "dbarts_setOffset", (DL_FUNC) &setOffset, 2 },
    { "dbarts_setPredictor", (DL_FUNC) &setPredictor, 3 },
    { "dbarts_setTestPredictor", (DL_FUNC) &setTestPredictor, 3 },
    { "dbarts_setTestPredictors", (DL_FUNC) &setTestPredictors, 3 },
    { "dbarts_setTestOffset", (DL_FUNC) &setTestOffset, 2 },
    { "dbarts_setControl", (DL_FUNC) &setControl, 2 },
    { "dbarts_isValidPointer", (DL_FUNC) &isValidPointer, 1 },
    { "dbarts_createState", (DL_FUNC) &createState, 1 },
    { "dbarts_storeState", (DL_FUNC) &storeState, 2 },
    { "dbarts_restoreState", (DL_FUNC) &restoreState, 2},
    { "dbarts_finalize", (DL_FUNC) &finalize, 0 },
    { NULL, NULL, 0 }
  };
  
  struct C_CallMethodDef {
    const char* package;
    const char* name;
    DL_FUNC function;
  };
  
  C_CallMethodDef C_callMethods[] = {
    { "dbarts", "createCGMPrior", (DL_FUNC) dbarts_createCGMPrior },
    { "dbarts", "createCGMPriorFromOptions", (DL_FUNC) dbarts_createCGMPriorFromOptions },
    { "dbarts", "destroyCGMPrior", (DL_FUNC) dbarts_destroyCGMPrior },
    { "dbarts", "initializeCGMPriorFromOptions", (DL_FUNC) dbarts_initializeCGMPriorFromOptions },
    { "dbarts", "invalidateCGMPrior", (DL_FUNC) dbarts_invalidateCGMPrior },
    
    { "dbarts", "createNormalPrior", (DL_FUNC) dbarts_createNormalPrior },
    { "dbarts", "createNormalPriorFromOptions", (DL_FUNC) dbarts_createNormalPriorFromOptions },
    { "dbarts", "destroyNormalPrior", (DL_FUNC) dbarts_destroyNormalPrior },
    { "dbarts", "initializeNormalPriorFromOptions", (DL_FUNC) dbarts_initializeNormalPriorFromOptions },
    { "dbarts", "invalidateNormalPrior", (DL_FUNC) dbarts_invalidateNormalPrior },
    
    { "dbarts", "createChiSquaredPrior", (DL_FUNC) dbarts_createChiSquaredPrior },
    { "dbarts", "createChiSquaredPriorFromOptions", (DL_FUNC) dbarts_createChiSquaredPriorFromOptions },
    { "dbarts", "destroyChiSquaredPrior", (DL_FUNC) dbarts_destroyChiSquaredPrior },
    { "dbarts", "initializeChiSquaredPriorFromOptions", (DL_FUNC) dbarts_initializeChiSquaredPriorFromOptions },
    { "dbarts", "invalidateChiSquaredPrior", (DL_FUNC) dbarts_invalidateChiSquaredPrior },

    { "dbarts", "createFit", (DL_FUNC) dbarts_createFit },
    { "dbarts", "initializeFit", (DL_FUNC) dbarts_initializeFit },
    { "dbarts", "destroyFit", (DL_FUNC) dbarts_destroyFit },
    { "dbarts", "invalidateFit", (DL_FUNC) dbarts_invalidateFit },
    
    { "dbarts", "runSampler", (DL_FUNC) dbarts_runSampler },
    { "dbarts", "runSamplerForIterations", (DL_FUNC) dbarts_runSamplerForIterations },
    { "dbarts", "setResponse", (DL_FUNC) dbarts_setResponse },
    { "dbarts", "setOffset", (DL_FUNC) dbarts_setOffset },
    { "dbarts", "setPredictor", (DL_FUNC) dbarts_setPredictor },
    { "dbarts", "setTestPredictor", (DL_FUNC) dbarts_setTestPredictor },
    { "dbarts", "setTestPredictors", (DL_FUNC) dbarts_setTestPredictors },
    { "dbarts", "setTestOffset", (DL_FUNC) dbarts_setTestOffset },
    { "dbarts", "setTestPredictorsAndOffset", (DL_FUNC) dbarts_setTestPredictorsAndOffset },
    { NULL, NULL, 0 }
  };
  
} // end anonymous namespace

extern "C" {
  void R_init_dbarts(DllInfo* info)
  {
    R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    
    C_CallMethodDef* method = C_callMethods;
    while (method->package != NULL) {
      R_RegisterCCallable(method->package, method->name, method->function);
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
  
  SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
  {
    SEXP val = allocVector(type, length);
    
    SET_SLOT(obj, nm, val);
    return val;
  }
  
  SEXP SET_DIMS(SEXP obj, int numRows, int numCols)
  {
    SEXP dimsExp = NEW_INTEGER(2);
    int *dims = INTEGER(dimsExp);
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
    if (!isLogical(slotExpr)) error("Binary response must be signified by logical type.");
    if (length(slotExpr) != 1) error("Binary response signifier must be of length 1.");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("Binary response must be either true or false.");
    control.responseIsBinary = (i_temp != FALSE);
    
    slotExpr = GET_ATTR(controlExpr, install("verbose"));
    if (!isLogical(slotExpr)) error("Verbose must be signified by logical type.");
    if (length(slotExpr) == 0) error("Verbose must be of length at least 1.");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("Verbose must be either true or false.");
    control.verbose = (i_temp != FALSE);
    
    slotExpr = GET_ATTR(controlExpr, install("keepTrainingFits"));
    if (!isLogical(slotExpr)) error("Keep training fits must be signified by logical type.");
    if (length(slotExpr) != 1) error("Keep training fits must be of length 1.");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("Keep training fits must be either true or false.");
    control.keepTrainingFits = (i_temp != FALSE);
    
    slotExpr = GET_ATTR(controlExpr, install("useQuantiles"));
    if (!isLogical(slotExpr)) error("Use quantiles must be signified by logical type.");
    if (length(slotExpr) != 1) error("Use quantiles must be of length 1.");
    i_temp = LOGICAL(slotExpr)[0];
    if (i_temp == NA_LOGICAL) error("Use quantiles must be either true or false.");
    control.useQuantiles = (i_temp != FALSE);
    
    
    
    slotExpr = GET_ATTR(controlExpr, install("n.samples"));
    if (!isInteger(slotExpr)) error("Number of samples must be of integer type.");
    if (length(slotExpr) != 1) error("Number of samples must be of length 1.");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) error("Number of samples cannot be NA.");
    if (i_temp <= 0) error("Number of samples must be positive.");
    control.numSamples = (size_t) i_temp;
    
    slotExpr = GET_ATTR(controlExpr, install("n.burn"));
    if (!isInteger(slotExpr)) error("Number of burn-in steps must be of integer type.");
    if (length(slotExpr) != 1) error("Number of burn-in steps must be of length 1.");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 0;
    if (i_temp < 0) error("Number of burn-in steps must be non-negative.");
    control.numBurnIn = (size_t) i_temp;
    
    slotExpr = GET_ATTR(controlExpr, install("n.trees"));
    if (!isInteger(slotExpr)) error("Number of trees must be of integer type.");
    if (length(slotExpr) != 1) error("Number of trees must be of length 1.");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) error("Number of trees cannot be NA.");
    if (i_temp <= 0) error("Number of trees must be positive.");
    control.numTrees = (size_t) i_temp;
    
    slotExpr = GET_ATTR(controlExpr, install("n.threads"));
    if (!isInteger(slotExpr)) error("Number of threads must be of integer type.");
    if (length(slotExpr) != 1) error("Number of threads must be of length 1.");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 1;
    if (i_temp <= 0) error("Number of threads must be positive.");
    control.numThreads = (size_t) i_temp;
    
    slotExpr = GET_ATTR(controlExpr, install("n.thin"));
    if (!isInteger(slotExpr)) error("Tree thinning rate must be of integer type.");
    if (length(slotExpr) != 1) error("Tree thinning rate must be of length 1.");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 1;
    if (i_temp < 0) error("Tree thinning rate must be non-negative.");
    control.treeThinningRate = (uint32_t) i_temp;
    
    
    slotExpr = GET_ATTR(controlExpr, install("printEvery"));
    if (!isInteger(slotExpr)) error("Print every must be of integer type.");
    if (length(slotExpr) != 1) error("Print every must be of length 1.");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp != NA_INTEGER) {
      if (i_temp <= 0) error("Print every must be positive.");
      control.printEvery = (uint32_t) i_temp;
    }
    
    slotExpr = GET_ATTR(controlExpr, install("printCutoffs"));
    if (!isInteger(slotExpr)) error("Print cutoffs must be of integer type.");
    if (length(slotExpr) == 0) error("Print cutoffs must be of length at least 1.");
    i_temp = INTEGER(slotExpr)[0];
    if (i_temp == NA_INTEGER) i_temp = 0;
    if (i_temp < 0) error("Print cutoffs must be non-negative.");
    control.printCutoffs = (uint32_t) i_temp;
  }
  
  void initializeModelFromExpression(Model& model, SEXP modelExpr, const Control& control)
  {
    double d_temp;
    
    SEXP slotExpr = GET_ATTR(modelExpr, install("p.birth_death"));
    if (!isReal(slotExpr)) error("Probability of birth/death rule must be of numeric type.");
    if (length(slotExpr) != 1) error("Probability of birth/death rule must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("Probability of birth/death rule must be a real number.");
    if (d_temp <= 0.0 || d_temp > 1.0) error("Probability of birth/death rule must be in (0, 1].");
    model.birthOrDeathProbability = d_temp;
    
    slotExpr = GET_ATTR(modelExpr, install("p.swap"));
    if (!isReal(slotExpr)) error("Probability of swap rule must be of numeric type.");
    if (length(slotExpr) != 1) error("Probability of swap rule must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("Probability of swap rule must be a real number.");
    if (d_temp < 0.0 || d_temp >= 1.0) error("Probability of swap rule must be in [0, 1).");
    model.swapProbability = d_temp;
    
    slotExpr = GET_ATTR(modelExpr, install("p.change"));
    if (!isReal(slotExpr)) error("Probability of change rule must be of numeric type.");
    if (length(slotExpr) != 1) error("Probability of change rule must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("Probability of change rule must be a real number.");
    if (d_temp < 0.0 || d_temp >= 1.0) error("Probability of change rule must be in [0, 1).");
    model.changeProbability = d_temp;
    
    if (std::fabs(model.birthOrDeathProbability + model.swapProbability + model.changeProbability - 1.0) >= 1.0e-10)
      error("Rule proposal probabilities must sum to 1.0");
    
    slotExpr = GET_ATTR(modelExpr, install("p.birth"));
    if (!isReal(slotExpr)) error("Probability of birth in birth/death rule must be of numeric type.");
    if (length(slotExpr) != 1) error("Probability of birth in birth/death rule must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("Probability of birth in birth/death rule must be a real number.");
    if (d_temp <= 0.0 || d_temp >= 1.0) error("Probability of birth in birth/death rule must be in (0, 1).");
    model.birthProbability = d_temp;
    
    
    SEXP priorExpr = GET_ATTR(modelExpr, install("tree.prior"));
    // slotExpr = GET_CLASS(priorExpr);
    // if (strcmp(CHAR(STRING_ELT(GET_CLASS(slotExpr), 0)), "dbartsControl") != 0) error("'control' argument to dbarts_create not of class 'dbartsControl'.");
    CGMPrior* treePrior = new CGMPrior;
    model.treePrior = treePrior;
    
    slotExpr = GET_ATTR(priorExpr, install("power"));
    if (!isReal(slotExpr)) error("Tree prior power must be of type real.");
    if (length(slotExpr) != 1) error("Tree prior power must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("Tree prior power be a real number.");
    if (d_temp <= 0.0) error("Tree prior power must be positive.");
    treePrior->power = d_temp;
    
    slotExpr = GET_ATTR(priorExpr, install("base"));
    if (!isReal(slotExpr)) error("Tree prior base must be of type real.");
    if (length(slotExpr) != 1) error("Tree prior power must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("Tree prior base be a real number.");
    if (d_temp <= 0.0 || d_temp >= 1.0) error("Tree prior base must be in (0, 1).");
    treePrior->base = d_temp;
    
    
    priorExpr = GET_ATTR(modelExpr, install("node.prior"));
    
    slotExpr = GET_ATTR(priorExpr, install("k"));
    if (!isReal(slotExpr)) error ("k must be of type real.");
    if (length(slotExpr) != 1) error("k must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("k must be a real number.");
    if (d_temp <= 0.0) error("k must be positive.");
    model.muPrior = new NormalPrior(control, d_temp);
    
    
    
    priorExpr = GET_ATTR(modelExpr, install("resid.prior"));
    
    slotExpr = GET_ATTR(priorExpr, install("df"));
    if (!isReal(slotExpr)) error("sigma prior degrees of freedom must be of type real.");
    if (length(slotExpr) != 1) error("sigma prior degrees of freedom must be of length 1.");
    double sigmaPriorDf = REAL(slotExpr)[0];
    if (ISNAN(sigmaPriorDf)) error("sigma prior degrees of freedom must be a real number.");
    if (sigmaPriorDf <= 0.0) error("sigma prior degrees of freedom must be positive.");
    
    slotExpr = GET_ATTR(priorExpr, install("quantile"));
    if (!isReal(slotExpr)) error ("sigma prior quantile must be of type real.");
    if (length(slotExpr) != 1) error("sigma prior quantile must be of length 1.");
    d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) error("sigma prior quantile must be a real number.");
    if (d_temp <= 0.0 || d_temp >= 1.0) error("sigma prior quantile must be in (0, 1).");
    model.sigmaSqPrior = new ChiSquaredPrior(sigmaPriorDf, d_temp);
  }

  void initializeDataFromExpression(Data& data, SEXP dataExpr)
  {
    int* dims;
    
    SEXP slotExpr = GET_ATTR(dataExpr, install("y"));
    if (!isReal(slotExpr)) error("y must be of type real.");
    if (length(slotExpr) == 0) error("Length of y must be greater than 0.");
    data.y = REAL(slotExpr);
    data.numObservations = (size_t) length(slotExpr);
    
    slotExpr = GET_ATTR(dataExpr, install("x"));
    if (!isReal(slotExpr)) error("x must be of type real.");
    dims = INTEGER(GET_ATTR(slotExpr, R_DimSymbol));
    if (dims == NULL || length(GET_ATTR(slotExpr, R_DimSymbol)) != 2) error("x must be a matrix, i.e. have two dimensions.");
    if ((size_t) dims[0] != data.numObservations) error("Number of rows of x and length of y must be equal.");
    data.X = REAL(slotExpr);
    data.numPredictors = (size_t) dims[1];
    
    slotExpr = GET_ATTR(dataExpr, install("varTypes"));
    if (!isInteger(slotExpr)) error("Variable types must be of type integer.");
    if ((size_t) length(slotExpr) != data.numPredictors) error("Length of variable types must equal number of columns in x.");
    int* i_variableTypes = INTEGER(slotExpr);
    VariableType* variableTypes = new VariableType[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) variableTypes[i] = (i_variableTypes[i] == 0 ? ORDINAL : CATEGORICAL);
    data.variableTypes = variableTypes;
    
    slotExpr = GET_ATTR(dataExpr, install("x.test"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.X_test = NULL;
      data.numTestObservations = 0;
    } else {
      if (!isReal(slotExpr)) error ("x.test must be of type real.");
      dims = INTEGER(GET_ATTR(slotExpr, R_DimSymbol));
      if (dims == NULL || length(GET_ATTR(slotExpr, R_DimSymbol)) != 2) error("x.test must be a matrix, i.e. have two dimensions.");
      if ((size_t) dims[1] != data.numPredictors) error("Number of columns of x.test and x must be equal.");
      data.X_test = REAL(slotExpr);
      data.numTestObservations = (size_t) dims[0];
    }
    
    slotExpr = GET_ATTR(dataExpr, install("weights"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.weights = NULL;
    } else {
      if (!isReal(slotExpr)) error("weights must be of type real.");
      if ((size_t) length(slotExpr) != data.numObservations) error("Length of weights must equal length of y.");
      data.weights = REAL(slotExpr);
    }
    
    slotExpr = GET_ATTR(dataExpr, install("offset"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.offset = NULL;
    } else {
      if (!isReal(slotExpr)) error("offset must be of type real.");
      if ((size_t) length(slotExpr) != data.numObservations) error("Length of offset must equal length of y.");
      data.offset = REAL(slotExpr);
    }
    
    slotExpr = GET_ATTR(dataExpr, install("offset.test"));
    if (isS4Null(slotExpr) || isNull(slotExpr) || length(slotExpr) == 0) {
      data.testOffset = NULL;
    } else {
      if (!isReal(slotExpr)) error("Test offset must be of type real.");
      if ((size_t) length(slotExpr) != data.numTestObservations) error("Length of test offset must equal number of test rows.");
      data.testOffset = REAL(slotExpr);
    }
    
    slotExpr = GET_ATTR(dataExpr, install("sigma"));
    if (!isReal(slotExpr)) error("sigma estimate must be of type real.");
    if (length(slotExpr) != 1) error("sigma estimate must be of length 1.");
    double d_temp = REAL(slotExpr)[0];
    if (ISNAN(d_temp)) d_temp = 1.0;
    if (d_temp <= 0.0) error("sigma estimate must be positive.");
    data.sigmaEstimate = d_temp;
    
    
    slotExpr = GET_ATTR(dataExpr, install("n.cuts"));
    if (!isInteger(slotExpr)) error("Maximum number of cuts must be of integer type.");
    if ((size_t) length(slotExpr) != data.numPredictors) error("Length of maximum number of cuts and the number of columns of x must be equal.");
    int* i_maxNumCuts = INTEGER(slotExpr);
    uint32_t* maxNumCuts = new uint32_t[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) maxNumCuts[i] = (uint32_t) i_maxNumCuts[i];
    data.maxNumCuts = maxNumCuts;
  }
  
  SEXP createStateExpressionFromFit(const BARTFit& fit)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const State& state(fit.state);
    
    SEXP result = PROTECT(NEW_OBJECT(MAKE_CLASS("dbartsState")));
    
    SEXP slotExpr = ALLOC_SLOT(result, install("fit.tree"), REALSXP, (int) (data.numObservations * control.numTrees));
    SET_DIMS(slotExpr, (int) data.numObservations, (int) control.numTrees);
    std::memcpy(REAL(slotExpr), state.treeFits, data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = ALLOC_SLOT(result, install("fit.total"), REALSXP, (int) data.numObservations);
    std::memcpy(REAL(slotExpr), state.totalFits, data.numObservations * sizeof(double));
    
    if (data.numTestObservations == 0) {
      SET_SLOT(result, install("fit.test"), NULL_USER_OBJECT);
    } else {
      slotExpr = ALLOC_SLOT(result, install("fit.test"), REALSXP, (int) data.numTestObservations);
      std::memcpy(REAL(slotExpr), state.totalTestFits, data.numTestObservations * sizeof(double));
    }
    
    slotExpr = ALLOC_SLOT(result, install("sigma"), REALSXP, 1);
    REAL(slotExpr)[0] = state.sigma;
    
    slotExpr = ALLOC_SLOT(result, install("runningTime"), REALSXP, 1);
    REAL(slotExpr)[0] = state.runningTime;
    
    slotExpr = ALLOC_SLOT(result, install("trees"), STRSXP, (int) control.numTrees);

    const char** treeStrings = const_cast<const char**>(state.createTreeStrings(fit));
    for (size_t i = 0; i < control.numTrees; ++i) {
      SET_STRING_ELT(slotExpr, (int) i, CREATE_STRING_VECTOR(treeStrings[i]));
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
    if (GET_LENGTH(dimsExpr) != 2) error("Dimensions of state@fit.tree indicate that it is not a matrix.");
    int* dims = INTEGER(dimsExpr);
    if ((size_t) dims[0] != data.numObservations || (size_t) dims[1] != control.numTrees) error("Dimensions of state@fit.tree do not match object.");
    std::memcpy(REAL(slotExpr), state.treeFits, data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = GET_ATTR(stateExpr, install("fit.total"));
    if ((size_t) GET_LENGTH(slotExpr) != data.numObservations) error("Length of state@fit.total does not match object.");
    std::memcpy(REAL(slotExpr), state.totalFits, data.numObservations * sizeof(double));
    
    if (data.numTestObservations != 0) {
      slotExpr = GET_ATTR(stateExpr, install("fit.test"));
      if ((size_t) GET_LENGTH(slotExpr) != data.numTestObservations) error("Length of state@fit.test does not match object.");
      std::memcpy(REAL(slotExpr), state.totalTestFits, data.numTestObservations * sizeof(double));
    }
    
    slotExpr = GET_ATTR(stateExpr, install("sigma"));
    if (GET_LENGTH(slotExpr) != 1) error("Length of state@sigma does not match object.");
    REAL(slotExpr)[0] = state.sigma;
    
    slotExpr = GET_ATTR(stateExpr, install("runningTime"));
    if (GET_LENGTH(slotExpr) != 1) error("Length of state@runningTime does not match object.");
    REAL(slotExpr)[0] = state.runningTime;
    
    slotExpr = GET_ATTR(stateExpr, install("trees"));
    if ((size_t) GET_LENGTH(slotExpr) != control.numTrees) error("Length of state@trees does not match object.");
    
    const char** treeStrings = const_cast<const char**>(state.createTreeStrings(fit));
    for (size_t i = 0; i < control.numTrees; ++i) {
      SET_STRING_ELT(slotExpr, (int) i, CREATE_STRING_VECTOR(treeStrings[i]));
      delete [] treeStrings[i];
    }
    delete [] treeStrings;
  }
  
  void initializeStateFromExpression(const BARTFit& fit, State& state, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    
    SEXP slotExpr = GET_ATTR(stateExpr, install("fit.tree"));
    std::memcpy(state.treeFits, (const double*) REAL(slotExpr), data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = GET_ATTR(stateExpr, install("fit.total"));
    std::memcpy(state.totalFits, (const double*) REAL(slotExpr), data.numObservations * sizeof(double));
    
    if (data.numTestObservations != 0) {
      slotExpr = GET_ATTR(stateExpr, install("fit.test"));
      std::memcpy(state.totalTestFits, (const double*) REAL(slotExpr), data.numTestObservations * sizeof(double));
    }
    
    slotExpr = GET_ATTR(stateExpr, install("sigma"));
    state.sigma = REAL(slotExpr)[0];
    
    slotExpr = GET_ATTR(stateExpr, install("runningTime"));
    state.runningTime = REAL(slotExpr)[0];
    
    slotExpr = GET_ATTR(stateExpr, install("trees"));
    const char** treeStrings = ext_stackAllocate(control.numTrees, const char*);
    for (size_t i = 0; i < control.numTrees; ++i) {
      treeStrings[i] = CHAR(STRING_ELT(slotExpr, (int) i));
    }
    state.recreateTreesFromStrings(fit, treeStrings);
    
    ext_stackFree(treeStrings);
  }
  
  void deleteFit(BARTFit* fit) {
#ifdef THREAD_SAFE_UNLOAD
    Rprintf("deleting   %p\n", fit);
#endif
    if (fit == NULL) return;
    
    delete fit->model.treePrior;
    delete fit->model.muPrior;
    delete fit->model.sigmaSqPrior;
    
    delete [] fit->data.variableTypes;
    delete [] fit->data.maxNumCuts;
    
    delete fit;
  }
}
