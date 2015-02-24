#include "xbart.hpp"

#include <cstddef> //size_t
#include <cstring> // memcpy

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>

using std::size_t;
using namespace dbarts;

SEXP xbart(SEXP fitExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr, SEXP ntreeExpr, SEXP nskipExpr,
           SEXP KExpr, SEXP resultTypeExpr, SEXP dropExpr)
{
  BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
  if (fit == NULL) error("xbart called on NULL external pointer");
  
  size_t numKs     = (size_t) LENGTH(xExpr);
  size_t numPowers = (size_t) LENGTH(numPowers);
  size_t numBases  = (size_t) LENGTH(baseExpr);
  size_t numNTrees = (size_t) LENGTH(ntreeExpr);
  
  double* k     = REAL(kExpr);
  double* power = REAL(powerExpr);
  double* base  = REAL(baseExpr);
  double* ntree = REAL(ntreeExpr);
  
  Control& control(*fit->control);
  Data& data(*fit->data);
  
  double* x_orig = new double[data.numObservations * data.numPredictors];
  std::memcpy(x_orig, data.X, data.numObservations * data.numPredictors * sizeof(double));
  double* y_orig = new double[data.numObservations];
  
  double* x_rep = NULL;
}

