#include "config.hpp"

#include <cstddef> // size_t
#include <cstdint> // uint32_t
#include <cstring> // memcpy

#include <external/R.h> // Rprintf, R_FlushConsole
#include <external/Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include <rc/util.h>

#include <misc/io.h>
#include <misc/simd.h>

#include <dbarts/dbarts.h>

#include "guessNumCores.hpp"
#include "makeModelMatrixFromDataFrame.h"

#include "R_interface_bartcore.hpp"
#include "R_interface_rbart.hpp"

using std::size_t;
using std::uint32_t;

extern "C" {
void R_init_dbarts(DllInfo* info);
}

extern "C" {

static SEXP assignInPlace(SEXP targetExpr, SEXP indexExpr, SEXP sourceExpr) {
  if (Rf_isReal(targetExpr)) {
    if (!Rf_isNull(Rf_getAttrib(targetExpr, R_DimSymbol))) {
      SEXP dimsExpr = Rf_getAttrib(targetExpr, R_DimSymbol);
      R_xlen_t numDims = XLENGTH(dimsExpr);

      int* dims = INTEGER(dimsExpr);
      int* indices = INTEGER(indexExpr);

      size_t length = rc_getLength(sourceExpr);

      if (XLENGTH(indexExpr) == numDims) {
        if (length != 1)
          Rf_error(
            "source must be a scalar when all array dimensions are specified"
          );

        size_t offset = 0;
        size_t stride = 1;
        for (R_xlen_t i = 0; i < numDims; ++i) {
          offset += (indices[i] - 1) * stride;
          stride *= dims[i];
        }

        REAL(targetExpr)[offset] = REAL(sourceExpr)[0];

        return R_NilValue;
      }

      if (XLENGTH(indexExpr) != numDims - 1)
        Rf_error("all but the first array dimension must be specified");

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
        if (length != 1)
          Rf_error(
            "source must be a scalar when all array dimensions are specified"
          );

        size_t offset = 0;
        size_t stride = 1;
        for (R_xlen_t i = 0; i < numDims; ++i) {
          offset += (indices[i] - 1) * stride;
          stride *= dims[i];
        }

        INTEGER(targetExpr)[offset] = INTEGER(sourceExpr)[0];

        return R_NilValue;
      }

      if (XLENGTH(indexExpr) != numDims - 1)
        Rf_error("all but the first array dimension must be specified");

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

static SEXP guessNumCores() {
  uint32_t numPhyiscalProcessors, numLogicalProcessors;
  dbarts::guessNumCores(&numPhyiscalProcessors, &numLogicalProcessors);

  SEXP resultExpr = Rf_allocVector(INTSXP, 2);
  PROTECT(resultExpr);
  int* result = INTEGER(resultExpr);

  result[0] = numPhyiscalProcessors <= 0
                ? NA_INTEGER
                : static_cast<int>(numPhyiscalProcessors);
  result[1] = numLogicalProcessors <= 0
                ? NA_INTEGER
                : static_cast<int>(numLogicalProcessors);

  UNPROTECT(1);

  return resultExpr;
}

// samplers release their resources through per-pointer R finalizers; the
// unload hook has nothing left to do
static SEXP finalize(void) { return R_NilValue; }

static SEXP deepCopy(SEXP obj) { return Rf_duplicate(obj); }

static SEXP setSIMDInstructionSet(SEXP i) {
  misc_simd_setSIMDInstructionSet(
    static_cast<misc_simd_instructionSet>(INTEGER(i)[0])
  );
  return R_NilValue;
}

static SEXP getMaxSIMDInstructionSet() {
  misc_simd_instructionSet result = misc_simd_getMaxSIMDInstructionSet();

  return Rf_ScalarInteger(static_cast<int>(result));
}

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
    sizeof(To) == sizeof(From) && std::is_trivially_copyable<From>::value &&
      std::is_trivially_copyable<To>::value,
    To>::type
  // constexpr support needs compiler magic
  bit_cast(const From& src) noexcept {
    static_assert(
      std::is_trivially_constructible<To>::value,
      "This implementation additionally requires destination type to be "
      "trivially constructible"
    );

    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
  }

#  else

  // We are only using this to cast function pointers, which are trivially
  // copiable. is_trivially_copyable is compiler specific and isn't worth trying
  // to reimplement in c++98.
  template <class To, class From> To bit_cast(const From& src) {
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
  }

#  endif

} // namespace std

#endif

extern "C" {
#define DEF_FUNC(_N_, _F_, _A_) {_N_, std::bit_cast<DL_FUNC>(&(_F_)), _A_}

static R_CallMethodDef R_callMethods[] = {
  DEF_FUNC("dbarts_bartcore_create", bartcore_create, 4),
  DEF_FUNC("dbarts_bartcore_createDataHandle", bartcore_createDataHandle, 2),
  DEF_FUNC("dbarts_bartcore_createFromHandle", bartcore_createFromHandle, 7),
  DEF_FUNC("dbarts_bartcore_createBCF", bartcore_createBCF, 5),
  DEF_FUNC("dbarts_bartcore_setTreatment", bartcore_setTreatment, 2),
  DEF_FUNC("dbarts_bartcore_getBCFGlue", bartcore_getBCFGlue, 1),
  DEF_FUNC("dbarts_bartcore_getForestFits", bartcore_getForestFits, 2),
  DEF_FUNC("dbarts_bartcore_run", bartcore_run, 3),
  DEF_FUNC("dbarts_bartcore_runWithCallback", bartcore_runWithCallback, 6),
  DEF_FUNC("dbarts_bartcore_setOffset", bartcore_setOffset, 3),
  DEF_FUNC("dbarts_bartcore_setResponse", bartcore_setResponse, 3),
  DEF_FUNC("dbarts_bartcore_setSigma", bartcore_setSigma, 2),
  DEF_FUNC("dbarts_bartcore_setData", bartcore_setData, 2),
  DEF_FUNC("dbarts_bartcore_setTestPredictor", bartcore_setTestPredictor, 2),
  DEF_FUNC("dbarts_bartcore_setTestOffset", bartcore_setTestOffset, 2),
  DEF_FUNC("dbarts_bartcore_setTestPredictorAndOffset",
           bartcore_setTestPredictorAndOffset, 3),
  DEF_FUNC("dbarts_bartcore_setWeights", bartcore_setWeights, 2),
  DEF_FUNC("dbarts_bartcore_setPredictor", bartcore_setPredictor, 4),
  DEF_FUNC("dbarts_bartcore_updatePredictor", bartcore_updatePredictor, 5),
  DEF_FUNC("dbarts_bartcore_updatePredictorPerObservation",
           bartcore_updatePredictorPerObservation, 3),
  DEF_FUNC("dbarts_bartcore_updatePredictorPerObservationJointly",
           bartcore_updatePredictorPerObservationJointly, 3),
  DEF_FUNC("dbarts_bartcore_setCutPoints", bartcore_setCutPoints, 3),
  DEF_FUNC("dbarts_bartcore_getSigmas", bartcore_getSigmas, 1),
  DEF_FUNC("dbarts_bartcore_isValidPointer", bartcore_isValidPointer, 1),
  DEF_FUNC("dbarts_bartcore_getLatents", bartcore_getLatents, 2),
  DEF_FUNC("dbarts_bartcore_predict", bartcore_predict, 3),
  DEF_FUNC("dbarts_bartcore_getTrees", bartcore_getTrees, 6),
  DEF_FUNC("dbarts_bartcore_storeState", bartcore_storeState, 1),
  DEF_FUNC("dbarts_bartcore_setState", bartcore_setState, 2),
  DEF_FUNC("dbarts_bartcore_installForests", bartcore_installForests, 3),
  DEF_FUNC("dbarts_bartcore_sampleTreesFromPrior",
           bartcore_sampleTreesFromPrior, 1),
  DEF_FUNC("dbarts_bartcore_sampleNodeParametersFromPrior",
           bartcore_sampleNodeParametersFromPrior, 1),
  DEF_FUNC("dbarts_bartcore_growFromRoot", bartcore_growFromRoot, 2),
  DEF_FUNC("dbarts_bartcore_printTrees", bartcore_printTrees, 4),
  DEF_FUNC("dbarts_bartcore_setControl", bartcore_setControl, 2),
  DEF_FUNC("dbarts_bartcore_setModel", bartcore_setModel, 4),
  DEF_FUNC("dbarts_bartcore_getSumsOfSquaredResiduals",
           bartcore_getSumsOfSquaredResiduals, 1),
  DEF_FUNC("dbarts_finalize", finalize, 0),
  DEF_FUNC("dbarts_deepCopy", deepCopy, 1),
  // DEF_FUNC("dbarts_getPointerAddress", getPointerAddress, 1),
  // DEF_FUNC("dbarts_getXAddress", getXAddress, 1),
  DEF_FUNC(
    "dbarts_makeModelMatrixFromDataFrame",
    dbarts_makeModelMatrixFromDataFrame,
    2
  ),
  DEF_FUNC("dbarts_guessNumCores", ::guessNumCores, 0),
  // experimental
  DEF_FUNC("dbarts_assignInPlace", assignInPlace, 3),
  // below: testing
  DEF_FUNC("dbarts_setSIMDInstructionSet", setSIMDInstructionSet, 1),
  DEF_FUNC("dbarts_getMaxSIMDInstructionSet", getMaxSIMDInstructionSet, 0),

  DEF_FUNC("rbart_fitted", rbart_getFitted, 4),
  {NULL, NULL, 0}
};

#undef DEF_FUNC

typedef struct {
  const char* name;
  DL_FUNC function;
} C_CallMethodDef;

#define DEF_FUNC(_N_, _F_) {_N_, std::bit_cast<DL_FUNC>(&(_F_))}

static C_CallMethodDef C_callMethods[] = {
  // the flat C API (inst/include/dbarts/dbarts.h), registered under the
  // symbol names themselves
  DEF_FUNC("dbarts_apiVersion", dbarts_apiVersion),
  DEF_FUNC("dbarts_sampler_create", dbarts_sampler_create),
  DEF_FUNC("dbarts_sampler_destroy", dbarts_sampler_destroy),
  DEF_FUNC("dbarts_sampler_run", dbarts_sampler_run),
  DEF_FUNC(
    "dbarts_sampler_sampleTreesFromPrior",
    dbarts_sampler_sampleTreesFromPrior
  ),
  DEF_FUNC(
    "dbarts_sampler_sampleNodeParametersFromPrior",
    dbarts_sampler_sampleNodeParametersFromPrior
  ),
  DEF_FUNC("dbarts_sampler_setResponse", dbarts_sampler_setResponse),
  DEF_FUNC("dbarts_sampler_setOffset", dbarts_sampler_setOffset),
  DEF_FUNC("dbarts_sampler_setWeights", dbarts_sampler_setWeights),
  DEF_FUNC("dbarts_sampler_setSigma", dbarts_sampler_setSigma),
  DEF_FUNC("dbarts_sampler_setCallback", dbarts_sampler_setCallback),
  DEF_FUNC("dbarts_sampler_getLatents", dbarts_sampler_getLatents),
  DEF_FUNC("dbarts_sampler_setPredictor", dbarts_sampler_setPredictor),
  DEF_FUNC("dbarts_sampler_updatePredictor", dbarts_sampler_updatePredictor),
  DEF_FUNC(
    "dbarts_sampler_setTestPredictors",
    dbarts_sampler_setTestPredictors
  ),
  DEF_FUNC("dbarts_sampler_setTestOffset", dbarts_sampler_setTestOffset),
  DEF_FUNC("dbarts_sampler_predict", dbarts_sampler_predict),
  DEF_FUNC("dbarts_sampler_setTreeStorage", dbarts_sampler_setTreeStorage),
  DEF_FUNC("dbarts_sampler_getTrees", dbarts_sampler_getTrees),
  DEF_FUNC("dbarts_sampler_printTrees", dbarts_sampler_printTrees),
  DEF_FUNC("dbarts_sampler_storeState", dbarts_sampler_storeState),
  DEF_FUNC("dbarts_sampler_setState", dbarts_sampler_setState),
  DEF_FUNC("dbarts_sampler_setNumThreads", dbarts_sampler_setNumThreads),
  DEF_FUNC("dbarts_sampler_setNumThin", dbarts_sampler_setNumThin),
  DEF_FUNC("dbarts_sampler_setVerbose", dbarts_sampler_setVerbose),
  DEF_FUNC("dbarts_sampler_numObservations", dbarts_sampler_numObservations),
  DEF_FUNC("dbarts_sampler_numPredictors", dbarts_sampler_numPredictors),
  DEF_FUNC(
    "dbarts_sampler_numTestObservations",
    dbarts_sampler_numTestObservations
  ),
  DEF_FUNC("dbarts_sampler_numChains", dbarts_sampler_numChains),
  DEF_FUNC("dbarts_sampler_numTrees", dbarts_sampler_numTrees),
  DEF_FUNC("dbarts_sampler_numSavedSamples", dbarts_sampler_numSavedSamples),
  DEF_FUNC("dbarts_sampler_kIsSampled", dbarts_sampler_kIsSampled),
  DEF_FUNC("dbarts_sampler_usesDart", dbarts_sampler_usesDart),

  {NULL, 0}
};

#undef DEF_FUNC

void attribute_visible R_init_dbarts(DllInfo* info) {
  R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
  R_useDynamicSymbols(info, static_cast<Rboolean>(FALSE));

  C_CallMethodDef* method = C_callMethods;
  while (method->name != NULL) {
    R_RegisterCCallable("dbarts", method->name, method->function);
    ++method;
  }

  misc_simd_init();

  misc_printf = &Rprintf;
  misc_flushOutput = &R_FlushConsole;
}
}
