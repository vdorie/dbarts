#ifndef R_INTERFACE_RBART_HPP
#define R_INTERFACE_RBART_HPP

#include <external/Rinternals.h> // SEXP

extern "C" {
  SEXP rbart_getFitted(SEXP yhatExpr, SEXP ranefExpr, SEXP groupByExpr, SEXP responseIsBinaryExpr);
}

#endif

