#ifndef XBART_HPP
#define XBART_HPP

#include <R.h>
#include <Rdefines.h>

extern "C" {

SEXP xbart(SEXP fitExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr, SEXP ntreeExpr, SEXP dropExpr);

}

#endif // XBART_HPP

