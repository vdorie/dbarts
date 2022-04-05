#include "config.hpp"
#include "R_interface.hpp"

#include <cstddef> // size_t
#include <dbarts/cstdint.hpp>
#include <cstring> // memcpy

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include <rc/util.h>

#include <misc/simd.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/R_C_interface.hpp>

#include "guessNumCores.hpp"
#include "makeModelMatrixFromDataFrame.h"

#include "R_interface_common.hpp"
#include "R_interface_crossvalidate.hpp"
#include "R_interface_rbart.hpp"
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
  
  static SEXP setSIMDInstructionSet(SEXP i)
  {
    misc_simd_setSIMDInstructionSet(static_cast<misc_simd_instructionSet>(INTEGER(i)[0]));
    return R_NilValue;
  }

  static SEXP getMaxSIMDInstructionSet()
  {
    misc_simd_instructionSet result = misc_simd_getMaxSIMDInstructionSet();
    
    return Rf_ScalarInteger(static_cast<int>(result));
  }
}

#include <misc/types.h>
extern "C" {
extern size_t misc_partitionRange_c(const misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length);
extern size_t misc_partitionIndices_c(const misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length);
extern size_t misc_partitionRange_neon(const misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length);
extern size_t misc_partitionIndices_neon(const misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length);
}
#include <misc/linearAlgebra.h>
#include <arm_neon.h>

#    define loadLHComp(_X_) \
       (values = (((uintptr_t) (x + _X_)) % (8 * sizeof(double))) == 0 ? \
         vld1q_u16(x + _X_) : \
         vcombine_u16(vcreate_u16(*((uint64_t*) (x + _X_))), vcreate_u16(*(((uint64_t*) (x + _X_ + 4))))), \
         vcgtq_u16(values, cut_vec))
#    define loadRHComp(_X_) \
       (values = (((uintptr_t) (x + _X_ - 7)) % (8 * sizeof(double))) == 0 ? \
         vld1q_u16(x + _X_ - 7) : \
         vrev64q_u16(vcombine_u16(vcreate_u16(*((uint64_t*) (x + _X_ - 3))), vcreate_u16(*((uint64_t*) (x + _X_ - 7))))), \
         vcleq_u16(values, cut_vec))

#    define vset_u16(_X7_, _X6_, _X5_, _X4_, _X3_, _X2_, _X1_, _X0_) \
       vcombine_u16( \
         vcreate_u16(((uint64_t) _X0_) + (((uint64_t) _X1_) << 16) + (((uint64_t) _X2_) << 32) + (((uint64_t) _X3_) << 48)), \
         vcreate_u16(((uint64_t) _X4_) + (((uint64_t) _X5_) << 16) + (((uint64_t) _X6_) << 32) + (((uint64_t) _X7_) << 48)))
#    define loadLHComp2(_X_) \
       (values = vset_u16(getDataAt(_X_ + 7), \
                          getDataAt(_X_ + 6), \
                          getDataAt(_X_ + 5), \
                          getDataAt(_X_ + 4), \
                          getDataAt(_X_ + 3), \
                          getDataAt(_X_ + 2), \
                          getDataAt(_X_ + 1), \
                          getDataAt(_X_    )), \
        vcgtq_u16(values, cut_vec))

# define getDataAt(_I_) x[_I_]
#    define loadRHComp2(_X_) \
       (values = vset_u16(getDataAt(_X_ - 7), \
                          getDataAt(_X_ - 6), \
                          getDataAt(_X_ - 5), \
                          getDataAt(_X_ - 4), \
                          getDataAt(_X_ - 3), \
                          getDataAt(_X_ - 2), \
                          getDataAt(_X_ - 1), \
                          getDataAt(_X_    )), \
        vcleq_u16(values, cut_vec))

extern "C" {
 static SEXP testPartition(SEXP xExpr, SEXP cutExpr, SEXP indicesExpr) {
   std::size_t n = rc_getLength(xExpr);

   uint16_t* x = new uint16_t[n];
   for (std::size_t i = 0; i < n; ++i)
     x[i] = static_cast<std::uint16_t>(INTEGER(xExpr)[i]);
   uint16_t cut = static_cast<uint16_t>(INTEGER(cutExpr)[0]);
  
   std::size_t* origIndices = new std::size_t[n]; 
   std::size_t* indices = new std::size_t[n]; 

   if (Rf_isNull(indicesExpr)) {
     for (std::size_t i = 0; i < n; ++i) {
       indices[i] = origIndices[i] = i;
     }
   } else {
     for (std::size_t i = 0; i < n; ++i) {
       indices[i] = origIndices[i] = static_cast<std::size_t>(INTEGER(indicesExpr)[i]);
     }
   }
   
   /* uint16x8_t setVal = vset_u16(x[indices[7]], x[indices[6]], x[indices[5]], x[indices[4]], x[indices[3]], x[indices[2]], x[indices[1]], x[indices[0]]);
   uint16_t out[8];
   vst1q_u16(out, setVal);
   Rprintf("%hu %hu %hu %hu %hu %hu %hu %hu\n", out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]); */
   uint16x8_t cut_vec = vdupq_n_u16(cut);

   uint16_t out[8];
   vst1q_u16(out, vld1q_u16(x));
   Rprintf("load lh: %hu %hu %hu %hu %hu %hu %hu %hu\n", out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
   uint16x8_t values = vcombine_u16(vcreate_u16(*((uint64_t*) (x))), vcreate_u16(*(((uint64_t*) (x + 4)))));
   vst1q_u16(out, values);
   Rprintf("create lh: %hu %hu %hu %hu %hu %hu %hu %hu\n", out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
   for (size_t j = 0; j <= n - 8; ++j) {
     loadLHComp(j);
     vst1q_u16(out, values);
     Rprintf("comp lh %lu: %hu %hu %hu %hu %hu %hu %hu %hu\n", j, out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
     loadLHComp2(j);
     vst1q_u16(out, values);
     Rprintf("comp lh2 %lu: %hu %hu %hu %hu %hu %hu %hu %hu\n", j, out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
   }

   size_t rh = 15;
   vst1q_u16(out, vld1q_u16(x + rh - 7));
   Rprintf("load rh: %hu %hu %hu %hu %hu %hu %hu %hu\n", out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
   values = vcombine_u16(vcreate_u16(*((uint64_t*) (x + rh - 7))), vcreate_u16(*(((uint64_t*) (x + rh - 3)))));
   vst1q_u16(out, values);
   Rprintf("create rh: %hu %hu %hu %hu %hu %hu %hu %hu\n", out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
   loadRHComp(rh);
   vst1q_u16(out, values);
   for (size_t j = 0; j <= n - 8; ++j) {
     loadRHComp(j + 7);
     vst1q_u16(out, values);
     Rprintf("comp rh %lu: %hu %hu %hu %hu %hu %hu %hu %hu\n", j, out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
     loadRHComp2(j + 7);
     vst1q_u16(out, values);
     Rprintf("comp rh2 %lu: %hu %hu %hu %hu %hu %hu %hu %hu\n", j, out[7], out[6], out[5], out[4], out[3], out[2], out[1], out[0]);
   }


   SEXP resultExpr = PROTECT(rc_setDims(rc_newInteger((n + 1) * 4), static_cast<int>(n + 1), 4, -1));
   int* result = INTEGER(resultExpr);

   std::size_t nleft = misc_partitionRange_c(x, cut, indices, n);
   for (std::size_t i = 0; i < n; ++i) {
     *result++ = indices[i];
     indices[i] = origIndices[i];
   }
   *result++ = static_cast<int>(nleft);

   nleft = misc_partitionIndices_c(x, cut, indices, n);
   for (std::size_t i = 0; i < n; ++i) {
     *result++ = indices[i];
     indices[i] = origIndices[i];
   }
   *result++ = static_cast<int>(nleft);

   nleft = misc_partitionRange_neon(x, cut, indices, n);
   for (std::size_t i = 0; i < n; ++i) {
     *result++ = indices[i];
     indices[i] = origIndices[i];
   }
   *result++ = static_cast<int>(nleft);

   nleft = misc_partitionIndices_neon(x, cut, indices, n);
   for (std::size_t i = 0; i < n; ++i) {
     *result++ = indices[i];
   }
   *result++ = static_cast<int>(nleft);


   delete [] indices;
   delete [] origIndices;
   delete [] x;
   UNPROTECT(1);
   return resultExpr;
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

}

#if __cplusplus >= 202002L
#  include <bit>
#else

namespace std {

#  if __cplusplus >= 201103L
#    include <type_traits>

// From https://en.cppreference.com/w/cpp/numeric/bit_cast
template <class To, class From>
typename std::enable_if<
  sizeof(To) == sizeof(From) &&
  std::is_trivially_copyable<From>::value &&
  std::is_trivially_copyable<To>::value,
  To>::type
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
  static_assert(std::is_trivially_constructible<To>::value,
    "This implementation additionally requires destination type to be trivially constructible");

  To dst;
  std::memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  else

// We are only using this to cast function pointers, which are trivially copiable.
// is_trivially_copyable is compiler specific and isn't worth trying to reimplement
// in c++98.
template <class To, class From>
To
bit_cast(const From& src)
{
  To dst;
  std::memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  endif

}

#endif

extern "C" {
#define DEF_FUNC(_N_, _F_, _A_) { _N_, std::bit_cast<DL_FUNC>(&_F_), _A_ }

  static R_CallMethodDef R_callMethods[] = {
    DEF_FUNC("dbarts_create", create, 3),
    DEF_FUNC("dbarts_run", run, 3),
    DEF_FUNC("dbarts_sampleTreesFromPrior", sampleTreesFromPrior, 1),
    DEF_FUNC("dbarts_sampleNodeParametersFromPrior", sampleNodeParametersFromPrior, 1),
    DEF_FUNC("dbarts_printTrees", printTrees, 4),
    DEF_FUNC("dbarts_getTrees", getTrees, 4),
    DEF_FUNC("dbarts_predict", predict, 3),
    DEF_FUNC("dbarts_setResponse", setResponse, 2),
    DEF_FUNC("dbarts_setOffset", setOffset, 3),
    DEF_FUNC("dbarts_setSigma", setSigma, 2),
    DEF_FUNC("dbarts_setWeights", setWeights, 2),
    DEF_FUNC("dbarts_setPredictor", setPredictor, 4),
    DEF_FUNC("dbarts_updatePredictor", updatePredictor, 5),
    DEF_FUNC("dbarts_setCutPoints", setCutPoints, 3),
    DEF_FUNC("dbarts_setTestPredictor", setTestPredictor, 2),
    DEF_FUNC("dbarts_setTestOffset", setTestOffset, 2),
    DEF_FUNC("dbarts_setTestPredictorAndOffset", setTestPredictorAndOffset, 3),
    DEF_FUNC("dbarts_updateTestPredictor", updateTestPredictor, 3),
    DEF_FUNC("dbarts_storeLatents", storeLatents, 2),
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
    DEF_FUNC("dbarts_assignInPlace", assignInPlace, 3),
    // below: testing
    DEF_FUNC("dbarts_setSIMDInstructionSet", setSIMDInstructionSet, 1),
    DEF_FUNC("dbarts_getMaxSIMDInstructionSet", getMaxSIMDInstructionSet, 0),
    
    DEF_FUNC("dbarts_testPartition", testPartition, 3),

    DEF_FUNC("rbart_fitted", rbart_getFitted, 4),
    { NULL, NULL, 0 }
  };

#undef DEF_FUNC
  
  typedef struct {
    const char* name;
    DL_FUNC function;
  } C_CallMethodDef;
  
#define DEF_FUNC(_N_, _F_) { _N_, std::bit_cast<DL_FUNC>(&_F_) }
  
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
    
    DEF_FUNC("createChiHyperprior", dbarts_createChiHyperprior),
    DEF_FUNC("createChiHyperpriorFromOptions", dbarts_createChiHyperpriorFromOptions),
    DEF_FUNC("destroyChiHyperprior", dbarts_destroyChiHyperprior),
    DEF_FUNC("initializeChiHyperpriorFromOptions", dbarts_initializeChiHyperpriorFromOptions),
    DEF_FUNC("invalidateChiHyperprior", dbarts_invalidateChiHyperprior),
    
    DEF_FUNC("createFixedHyperprior", dbarts_createFixedHyperprior),
    DEF_FUNC("createFixedHyperpriorFromOptions", dbarts_createFixedHyperpriorFromOptions),
    DEF_FUNC("destroyFixedHyperprior", dbarts_destroyFixedHyperprior),
    DEF_FUNC("initializeFixedHyperpriorFromOptions", dbarts_initializeFixedHyperpriorFromOptions),
    DEF_FUNC("invalidateFixedHyperprior", dbarts_invalidateFixedHyperprior),
    
    DEF_FUNC("createChiSquaredPrior", dbarts_createChiSquaredPrior),
    DEF_FUNC("createChiSquaredPriorFromOptions", dbarts_createChiSquaredPriorFromOptions),
    DEF_FUNC("destroyChiSquaredPrior", dbarts_destroyChiSquaredPrior),
    DEF_FUNC("initializeChiSquaredPriorFromOptions", dbarts_initializeChiSquaredPriorFromOptions),
    DEF_FUNC("invalidateChiSquaredPrior", dbarts_invalidateChiSquaredPrior),

    DEF_FUNC("createControl", dbarts_createControl),
    DEF_FUNC("initializeControl", dbarts_initializeControl),
    DEF_FUNC("destroyControl", dbarts_destroyControl),
    // DEF_FUNC("invalidateControl", dbarts_invalidateControl),
    DEF_FUNC("setControl", dbarts_setControl),
    
    DEF_FUNC("createData", dbarts_createData),
    DEF_FUNC("initializeData", dbarts_initializeData),
    DEF_FUNC("destroyData", dbarts_destroyData),
    DEF_FUNC("invalidateData", dbarts_invalidateData),
    
    DEF_FUNC("createModel", dbarts_createModel),
    DEF_FUNC("initializeModel", dbarts_initializeModel),
    DEF_FUNC("destroyModel", dbarts_destroyModel),
    DEF_FUNC("invalidateModel", dbarts_invalidateModel),
    
    DEF_FUNC("createFit", dbarts_createFit),
    DEF_FUNC("initializeFit", dbarts_initializeFit),
    DEF_FUNC("destroyFit", dbarts_destroyFit),
    DEF_FUNC("invalidateFit", dbarts_invalidateFit),
    
    DEF_FUNC("createStateExpression", dbarts_createStateExpression),
    DEF_FUNC("initializeState", dbarts_initializeState),
    
    DEF_FUNC("printInitialSummary", dbarts_printInitialSummary),
    DEF_FUNC("printTrees", dbarts_printTrees),
    DEF_FUNC("getTrees", dbarts_getTrees),
    DEF_FUNC("setRNGState", dbarts_setRNGState),
    
    DEF_FUNC("runSampler", dbarts_runSampler),
    DEF_FUNC("runSamplerForIterations", dbarts_runSamplerForIterations),
    DEF_FUNC("runSamplerWithResults", dbarts_runSamplerWithResults),
    DEF_FUNC("sampleTreesFromPrior", dbarts_sampleTreesFromPrior),
    DEF_FUNC("sampleNodeParametersFromPrior", dbarts_sampleNodeParametersFromPrior),
    DEF_FUNC("predict", dbarts_predict),
    DEF_FUNC("setResponse", dbarts_setResponse),
    DEF_FUNC("setOffset", dbarts_setOffset),
    DEF_FUNC("setSigma", dbarts_setSigma),
    DEF_FUNC("setPredictor", dbarts_setPredictor),
    DEF_FUNC("updatePredictor", dbarts_updatePredictor),
    DEF_FUNC("setTestPredictor", dbarts_setTestPredictor),
    DEF_FUNC("setTestOffset", dbarts_setTestOffset),
    DEF_FUNC("setTestPredictorsAndOffset", dbarts_setTestPredictorAndOffset),
    DEF_FUNC("updateTestPredictor", dbarts_updateTestPredictor),
    DEF_FUNC("updateTestPredictors", dbarts_updateTestPredictors),
    DEF_FUNC("storeLatents", dbarts_storeLatents),
    
    { NULL, 0 }
  };
  
#undef DEF_FUNC
  
  void attribute_visible R_init_dbarts(DllInfo* info)
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
    
    misc_simd_init();
  }
}

