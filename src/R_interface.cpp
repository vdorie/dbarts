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

/// Bounds-check R-supplied 1-based array indices against the target's
/// dimensions. assignInPlace runs once per Gibbs sweep, so this is a short
/// integer loop with no allocation and no R API call on the valid path.
static void assertIndicesInRange(
  const int* indices, const int* dims, R_xlen_t numIndices)
{
  for (R_xlen_t i = 0; i < numIndices; ++i)
    if (indices[i] < 1 || indices[i] > dims[i])
      Rf_error("assignInPlace: array index out of bounds");
}

} // extern "C" -- a function template cannot have C language linkage

/// Shared body of assignInPlace for a single element type. The dims/offset/
/// stride/bounds arithmetic is identical for REAL and INTEGER targets; only the
/// element type and the (already dereferenced) data pointers differ, so the
/// caller passes the target/source pointers obtained from REAL()/INTEGER().
template <typename T>
static void assignInPlaceTyped(
  SEXP targetExpr, SEXP indexExpr, SEXP sourceExpr, T* target, const T* source)
{
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

      assertIndicesInRange(indices, dims, numDims);

      size_t offset = 0;
      size_t stride = 1;
      for (R_xlen_t i = 0; i < numDims; ++i) {
        offset += (indices[i] - 1) * stride;
        stride *= dims[i];
      }

      target[offset] = source[0];

      return;
    }

    if (XLENGTH(indexExpr) != numDims - 1)
      Rf_error("all but the first array dimension must be specified");

    assertIndicesInRange(indices, dims + 1, numDims - 1);
    if (length != static_cast<size_t>(dims[0]))
      Rf_error(
        "assignInPlace: 'source' length must equal the target's leading dimension"
      );

    size_t offset = 0;
    size_t stride = dims[0];
    for (R_xlen_t i = 0; i < numDims - 1; ++i) {
      offset += (indices[i] - 1) * stride;
      stride *= dims[i + 1];
    }

    std::memcpy(target + offset, source, length * sizeof(T));
  } else {
    int index = INTEGER(indexExpr)[0];
    if (index < 1 || static_cast<R_xlen_t>(index) > XLENGTH(targetExpr))
      Rf_error("assignInPlace: index out of bounds");
    if (rc_getLength(sourceExpr) < 1)
      Rf_error("assignInPlace: 'source' must be non-empty");
    target[index - 1] = source[0];
  }
}

extern "C" {

static SEXP assignInPlace(SEXP targetExpr, SEXP indexExpr, SEXP sourceExpr) {
  if (!Rf_isInteger(indexExpr))
    Rf_error("assignInPlace: 'index' must be integer");

  if (Rf_isReal(targetExpr)) {
    if (!Rf_isReal(sourceExpr))
      Rf_error("assignInPlace: 'source' must be double for a double target");
    assignInPlaceTyped<double>(
      targetExpr, indexExpr, sourceExpr, REAL(targetExpr), REAL(sourceExpr));
  } else if (Rf_isInteger(targetExpr)) {
    if (!Rf_isInteger(sourceExpr))
      Rf_error("assignInPlace: 'source' must be integer for an integer target");
    assignInPlaceTyped<int>(
      targetExpr, indexExpr, sourceExpr, INTEGER(targetExpr), INTEGER(sourceExpr));
  }

  return R_NilValue;
}

static SEXP guessNumCores() {
  uint32_t numPhysicalProcessors, numLogicalProcessors;
  dbarts::guessNumCores(&numPhysicalProcessors, &numLogicalProcessors);

  SEXP resultExpr = Rf_allocVector(INTSXP, 2);
  PROTECT(resultExpr);
  int* result = INTEGER(resultExpr);

  result[0] = numPhysicalProcessors <= 0
                ? NA_INTEGER
                : static_cast<int>(numPhysicalProcessors);
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

#include <bit>

extern "C" {
#define DEF_FUNC(_N_, _F_, _A_) {_N_, std::bit_cast<DL_FUNC>(&(_F_)), _A_}

static R_CallMethodDef R_callMethods[] = {
  DEF_FUNC("dbarts_bartcore_create", bartcore_create, 4),
  DEF_FUNC("dbarts_bartcore_createDataHandle", bartcore_createDataHandle, 3),
  DEF_FUNC("dbarts_bartcore_createFromHandle", bartcore_createFromHandle, 8),
  DEF_FUNC("dbarts_bartcore_createBCF", bartcore_createBCF, 10),
  DEF_FUNC("dbarts_bartcore_createMultinomial", bartcore_createMultinomial, 6),
  DEF_FUNC("dbarts_bartcore_createMultinomialCounts",
           bartcore_createMultinomialCounts, 6),
  DEF_FUNC("dbarts_bartcore_setCounts", bartcore_setCounts, 2),
  DEF_FUNC("dbarts_bartcore_setCategoryOffset", bartcore_setCategoryOffset, 2),
  DEF_FUNC("dbarts_bartcore_setTreatment", bartcore_setTreatment, 2),
  DEF_FUNC("dbarts_bartcore_setForestWeights", bartcore_setForestWeights, 3),
  DEF_FUNC("dbarts_bartcore_getBCFGlue", bartcore_getBCFGlue, 1),
  DEF_FUNC("dbarts_bartcore_getForestFits", bartcore_getForestFits, 2),
  DEF_FUNC("dbarts_bartcore_getForestVariableCounts",
           bartcore_getForestVariableCounts, 2),
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
  DEF_FUNC("dbarts_bartcore_setCutPoints", bartcore_setCutPoints, 4),
  DEF_FUNC("dbarts_bartcore_getSigmas", bartcore_getSigmas, 1),
  DEF_FUNC("dbarts_bartcore_isValidPointer", bartcore_isValidPointer, 1),
  DEF_FUNC("dbarts_bartcore_getLatents", bartcore_getLatents, 2),
  DEF_FUNC("dbarts_bartcore_predict", bartcore_predict, 3),
  DEF_FUNC("dbarts_bartcore_getTrees", bartcore_getTrees, 8),
  DEF_FUNC("dbarts_bartcore_storeState", bartcore_storeState, 1),
  DEF_FUNC("dbarts_bartcore_setState", bartcore_setState, 3),
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

// One {name, &function} row per DBARTS_C_API_LIST entry; the header's stubs and
// dbarts's binding asserts expand the same list, so this table cannot drift from
// the declared surface.
#define DBARTS_API_REGISTER(_R_, _N_, _P_, _A_) \
  {#_N_, std::bit_cast<DL_FUNC>(&_N_)},

static C_CallMethodDef C_callMethods[] = {
  // the flat C API (inst/include/dbarts/dbarts.h), registered under the
  // symbol names themselves
  DBARTS_C_API_LIST(DBARTS_API_REGISTER)

  {NULL, 0}
};

#undef DBARTS_API_REGISTER

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
