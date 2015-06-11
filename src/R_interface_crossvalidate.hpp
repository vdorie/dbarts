#ifndef R_INTERFACE_CROSSVALIDATE_HPP
#define R_INTERFACE_CROSSVALIDATE_HPP

struct SEXPREC;
typedef struct SEXPREC* SEXP;

extern "C" {
  SEXP xbart(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr,
             SEXP KExpr, SEXP numRepsExpr, SEXP numBurnInExpr, SEXP lossTypeExpr, SEXP numThreadsExpr,
             SEXP numTreesExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr,
             SEXP dropExpr);
}

#endif

