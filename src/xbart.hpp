#ifndef XBART_HPP
#define XBART_HPP

#include <Rdefines.h> // SEXP

namespace dbarts {

  SEXP xbart(SEXP fitExpr, SEXP ntreeExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr, SEXP nskipExpr,
             SEXP KExpr, SEXP nrepsExpr, SEXP resultTypeExpr, SEXP dropExpr);

}

#endif // XBART_HPP

