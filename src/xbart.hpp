#ifndef XBART_HPP
#define XBART_HPP

#include <Rdefines.h> // SEXP

namespace dbarts {

  SEXP xbart(SEXP fitExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr, SEXP ntreeExpr, SEXP nskipExpr,
             SEXP KExpr, SEXP resultTypeExpr, SEXP dropExpr);

}

#endif // XBART_HPP

