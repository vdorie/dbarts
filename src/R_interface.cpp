#include "config.hpp"
#include "R_interface.hpp"

#include <cstddef> // size_t
#include <dbarts/cstdint.hpp>
#include <cstring> // memcpy

#include <R_ext/Rdynload.h>

#include <rc/util.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/R_C_interface.hpp>

#include "guessNumCores.hpp"
#include "makeModelMatrixFromDataFrame.h"

#include "R_interface_common.hpp"
#include "R_interface_crossvalidate.hpp"
#include "R_interface_sampler.hpp"

using std::size_t;
using std::uint32_t;
using namespace dbarts;

extern "C" {
  void R_init_dbarts(DllInfo* info);
}
  
PointerSet* activeFits;
#ifdef THREAD_SAFE_UNLOAD
pthread_mutex_t fitMutex;
#endif

bool compareExternalPointers(const SEXP& lhs, const SEXP& rhs) {
  return R_ExternalPtrAddr(const_cast<SEXP>(lhs)) < R_ExternalPtrAddr(const_cast<SEXP>(rhs));
}

extern "C" {
  
  static SEXP assignInPlace(SEXP targetExpr, SEXP indexExpr, SEXP sourceExpr)
  {
    if (Rf_isReal(targetExpr)) {
      if (!Rf_isNull(Rf_getAttrib(targetExpr, R_DimSymbol))) {
        SEXP dimsExpr = Rf_getAttrib(targetExpr, R_DimSymbol);
        R_xlen_t numDims = XLENGTH(dimsExpr);
        
        int* dims = INTEGER(dimsExpr);
        int* indices = INTEGER(indexExpr);
        
        size_t length = rc_getLength(sourceExpr);

        if (XLENGTH(indexExpr) == numDims) {
          if (length != 1) Rf_error("source must be a scalar when all array dimensions are specified");
          
          size_t offset = 0;
          size_t stride = 1;
          for (R_xlen_t i = 0; i < numDims; ++i) {
            offset += (indices[i] - 1) * stride;
            stride *= dims[i];
          }
          
          REAL(targetExpr)[offset] = REAL(sourceExpr)[0];
          
          return R_NilValue;
        }
        
        if (XLENGTH(indexExpr) != numDims - 1) Rf_error("all but the first array dimension must be specified");
        
       
        size_t offset = 0;
        size_t stride = dims[0];
        for (R_xlen_t i = 0; i < numDims - 1; ++i) {
          offset += (indices[i] - 1) * stride;
          stride *= dims[i + 1];
        }
        
        double* target = REAL(targetExpr);
        const double* source = REAL(sourceExpr);
        std::memcpy(target + offset, source, length * sizeof(double));
      } else {
        size_t index = INTEGER(indexExpr)[0] - 1;
        double* target = REAL(targetExpr);
        double source = REAL(sourceExpr)[0];
        target[index] = source;
      }
    } else if (Rf_isInteger(targetExpr)) {
      if (!Rf_isNull(Rf_getAttrib(targetExpr, R_DimSymbol))) {
        SEXP dimsExpr = Rf_getAttrib(targetExpr, R_DimSymbol);
        R_xlen_t numDims = XLENGTH(dimsExpr);
        
        int* dims = INTEGER(dimsExpr);
        int* indices = INTEGER(indexExpr);
        
        size_t length = rc_getLength(sourceExpr);
        
        if (XLENGTH(indexExpr) == numDims) {
          if (length != 1) Rf_error("source must be a scalar when all array dimensions are specified");
          
          size_t offset = 0;
          size_t stride = 1;
          for (R_xlen_t i = 0; i < numDims; ++i) {
            offset += (indices[i] - 1) * stride;
            stride *= dims[i];
          }
          
          INTEGER(targetExpr)[offset] = INTEGER(sourceExpr)[0];
          
          return R_NilValue;
        }
          
        if (XLENGTH(indexExpr) != numDims - 1) Rf_error("all but the first array dimension must be specified");
        
        size_t offset = 0;
        size_t stride = dims[0];
        for (R_xlen_t i = 0; i < numDims - 1; ++i) {
          offset += (indices[i] - 1) * stride;
          stride *= dims[i + 1];
        }
        
        int* target = INTEGER(targetExpr);
        const int* source = INTEGER(sourceExpr);
        std::memcpy(target + offset, source, length * sizeof(int));
      } else {
        size_t index = INTEGER(indexExpr)[0] - 1;
        int* target = INTEGER(targetExpr);
        int source = INTEGER(sourceExpr)[0];
        target[index] = source;
      }
    }
    
    return R_NilValue;
  }
    
  static SEXP isValidPointer(SEXP fitExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) return Rf_ScalarLogical(FALSE);
    
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_lock(&fitMutex);
#endif
    if (activeFits->find(fitExpr) != activeFits->end()) {
#ifdef THREAD_SAFE_UNLOAD
      pthread_mutex_unlock(&fitMutex);
#endif
      return Rf_ScalarLogical(TRUE);
    }
#ifdef THREAD_SAFE_UNLOAD
    pthread_mutex_unlock(&fitMutex);
#endif
    return Rf_ScalarLogical(FALSE);
  }
  
  static SEXP guessNumCores()
  {
    uint32_t numPhyiscalProcessors, numLogicalProcessors;
    dbarts::guessNumCores(&numPhyiscalProcessors, &numLogicalProcessors);
    
    SEXP resultExpr = Rf_allocVector(INTSXP, 2);
    PROTECT(resultExpr);
    int* result = INTEGER(resultExpr);
    
    result[0] = numPhyiscalProcessors <= 0 ? NA_INTEGER : static_cast<int>(numPhyiscalProcessors);
    result[1] = numLogicalProcessors  <= 0 ? NA_INTEGER : static_cast<int>(numLogicalProcessors);
    
    UNPROTECT(1);
    
    return resultExpr;
  }
  
  
  static SEXP finalize(void) {
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
    
    return R_NilValue;
  }
  
  static SEXP deepCopy(SEXP obj)
  {
    return Rf_duplicate(obj);
  }

/*
}

#ifdef HAVE_STD_SNPRINTF
// snprintf in c++11, before that have to use C version
#  include <cstdio>
using std::snprintf;
#else
extern "C" {
#  include <stdio.h>
}
#endif

namespace {
   static SEXP getPointerAddress(SEXP obj)
  {
    char buffer[24];
    const void* p;
    if (Rf_isInteger(obj)) {
      p = reinterpret_cast<const void*>(INTEGER(obj));
    } else if (Rf_isLogical(obj)) {
      p = reinterpret_cast<const void*>(LOGICAL(obj));
    } else if (Rf_isReal(obj)) {
      p = reinterpret_cast<const void*>(REAL(obj));
    } else if (Rf_isString(obj)) {
      p = reinterpret_cast<const void*>(CHAR(obj));
    } else {
      p = NULL;
    }
    snprintf(buffer, 24, "%p", p);
    
    SEXP result = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(result, 0, Rf_mkChar(buffer));
    UNPROTECT(1);
    
    return result;
  }
  
  static SEXP getXAddress(SEXP obj)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(obj));
    if (fit == NULL) Rf_error("dbarts_getXAddress called on NULL external pointer");
   
    char buffer[24];
    snprintf(buffer, 24, "%p", reinterpret_cast<const void*>(fit->data.x));
    
    SEXP result = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(result, 0, Rf_mkChar(buffer));
    UNPROTECT(1);
    
    return result;
  } */
    
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

  static R_CallMethodDef R_callMethods[] = {
    DEF_FUNC("dbarts_create", create, 3),
    DEF_FUNC("dbarts_run", run, 3),
    DEF_FUNC("dbarts_sampleTreesFromPrior", sampleTreesFromPrior, 1),
    DEF_FUNC("dbarts_printTrees", printTrees, 4),
    DEF_FUNC("dbarts_predict", predict, 3),
    DEF_FUNC("dbarts_setResponse", setResponse, 2),
    DEF_FUNC("dbarts_setOffset", setOffset, 2),
    DEF_FUNC("dbarts_setPredictor", setPredictor, 2),
    DEF_FUNC("dbarts_updatePredictor", updatePredictor, 3),
    DEF_FUNC("dbarts_setTestPredictor", setTestPredictor, 2),
    DEF_FUNC("dbarts_setTestOffset", setTestOffset, 2),
    DEF_FUNC("dbarts_setTestPredictorAndOffset", setTestPredictorAndOffset, 3),
    DEF_FUNC("dbarts_updateTestPredictor", updateTestPredictor, 3),
    DEF_FUNC("dbarts_setData", setData, 2),
    DEF_FUNC("dbarts_setControl", setControl, 2),
    DEF_FUNC("dbarts_setModel", setModel, 2),
    DEF_FUNC("dbarts_isValidPointer", isValidPointer, 1),
    DEF_FUNC("dbarts_createState", createState, 1),
    DEF_FUNC("dbarts_storeState", storeState, 2),
    DEF_FUNC("dbarts_restoreState", restoreState, 2),
    DEF_FUNC("dbarts_finalize", finalize, 0),
    DEF_FUNC("dbarts_deepCopy", deepCopy, 1),
    //DEF_FUNC("dbarts_getPointerAddress", getPointerAddress, 1),
    //DEF_FUNC("dbarts_getXAddress", getXAddress, 1),
    DEF_FUNC("dbarts_makeModelMatrixFromDataFrame", dbarts_makeModelMatrixFromDataFrame, 2),
    DEF_FUNC("dbarts_xbart", xbart, 14),
    DEF_FUNC("dbarts_guessNumCores", ::guessNumCores, 0),
    // experimental
    DEF_FUNC("dbarts_saveToFile", saveToFile, 2),
    DEF_FUNC("dbarts_loadFromFile", loadFromFile, 1),
    DEF_FUNC("dbarts_assignInPlace", assignInPlace, 3),
    // below: testing
    { NULL, NULL, 0 }
  };

#undef DEF_FUNC
  
  typedef struct {
    const char* name;
    DL_FUNC function;
  } C_CallMethodDef;
  
#define DEF_FUNC(_N_, _F_) { _N_, reinterpret_cast<DL_FUNC>(&_F_) }
  
  static C_CallMethodDef C_callMethods[] = {
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
    
    DEF_FUNC("setRNGState", dbarts_setRNGState),
    
    DEF_FUNC("runSampler", dbarts_runSampler),
    DEF_FUNC("runSamplerForIterations", dbarts_runSamplerForIterations),
    DEF_FUNC("sampleTreesFromPrior", dbarts_sampleTreesFromPrior),
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
  
}

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
    
    activeFits = new PointerSet(&compareExternalPointers);
  }
}

