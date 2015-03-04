#include "xbart.hpp"

#include <cstddef> // size_t
#include <cstring> // memcpy

#include <external/random.h>
#include <external/rc.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>

#include <Rinternals.h> // R_xlen_t

using std::size_t;
using namespace dbarts;

#define Z_(_X_) static_cast<R_xlen_t>(_X_)

void permuteIndexArray(ext_rng* generator, size_t* indices, size_t length);

SEXP xbart(SEXP fitExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr, SEXP ntreeExpr, SEXP nskipExpr,
           SEXP KExpr, SEXP resultTypeExpr, SEXP dropExpr)
{
  BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
  if (fit == NULL) error("xbart called on NULL external pointer");
  
  rc_checkInts(ntreeExpr, "num trees", RC_LENGTH | RC_GEQ, Z_(1), RC_VALUE | RC_GT, 0, RC_NA | RC_NO, RC_END);
  rc_checkInts(nskipExpr, "num skip", RC_LENGTH | RC_GEQ, Z_(1), RC_LENGTH | RC_LEQ, Z_(2), RC_VALUE | RC_GEQ, 0, RC_NA | RC_NO, RC_END);
  
  size_t numFolds = static_cast<size_t>(
    rc_getInt(KExpr, "num folds", RC_LENGTH | RC_EQ, Z_(1), RC_VALUE | RC_GT, 0, RC_NA | RC_NO, RC_END));
  
  size_t numKs     = (size_t) LENGTH(kExpr);
  size_t numPowers = (size_t) LENGTH(powerExpr);
  size_t numBases  = (size_t) LENGTH(baseExpr);
  size_t numNTrees = (size_t) LENGTH(ntreeExpr);
  
  double* k     = REAL(kExpr);
  double* power = REAL(powerExpr);
  double* base  = REAL(baseExpr);
  int* ntree = INTEGER(ntreeExpr);
  
  Control& control(fit->control);
  Data& data(fit->data);
  
  const double* x_orig = data.X;
  const double* y_orig = data.y;

  Data repData;
  
  size_t* permutation = new size_t[data.numObservations];
  for (size_t i = 0; i < data.numObservations; ++i) permutation[i] = i;
  
  permuteIndexArray(control.rng, permutation, data.numObservations);
  
  delete [] permutation;
}

void permuteIndexArray(ext_rng* generator, size_t* indices, size_t length)
{
  
  size_t temp, swapPos;
  for (size_t i = 0; i < length - 1; ++i) {
    swapPos = static_cast<size_t>(ext_rng_simulateUnsignedIntegerUniformInRange(generator, i, length));
    
    temp = indices[i];
    indices[i] = indices[swapPos];
    indices[swapPos] = temp;
  }
}

