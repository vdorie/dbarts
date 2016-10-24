#ifndef RC_UTIL_H
#define RC_UTIL_H

#include <stdbool.h>
#include <external/stddef.h> // ext_size_t

#include <external/Rinternals.h> // SEXP, Rf functions

#ifdef __cplusplus
extern "C" {
#endif

SEXP rc_setDims(SEXP obj, ...);
#define rc_getDims(_X_) Rf_getAttrib(_X_, R_DimSymbol)
SEXP rc_allocateInSlot(SEXP obj, SEXP slotName, SEXPTYPE type, R_xlen_t length);

#ifndef __cplusplus
#  define rc_getLength(_X_) ((ext_size_t) XLENGTH(_X_))
#  define rc_asRLength(_X_) ((R_xlen_t) (_X_))
#else
#  define rc_getLength(_X_) static_cast<ext_size_t>(XLENGTH(_X_))
#  define rc_asRLength(_X_) static_cast<R_xlen_t>(_X_)
#endif

#define rc_getNames(_X_)          Rf_getAttrib(_X_, R_NamesSymbol)
#define rc_setNames(_X_, _NAMES_) Rf_setAttrib(_X_, R_NamesSymbol, _NAMES_)

#define rc_getDimNames(_X_)          Rf_getAttrib(_X_, R_DimNamesSymbol)
#define rc_setDimNames(_X_, _NAMES_) Rf_setAttrib(_X_, R_DimNamesSymbol, _NAMES_)

#define rc_newCharacter(_N_) Rf_allocVector(STRSXP, _N_)
#define rc_newInteger(_N_)   Rf_allocVector(INTSXP, _N_)
#define rc_newList(_N_)      Rf_allocVector(VECSXP, _N_)
#define rc_newLogical(_N_)   Rf_allocVector(LGLSXP, _N_)
#define rc_newNumeric(_N_)   Rf_allocVector(REALSXP, _N_)

#define rc_getClass(_X_)  Rf_getAttrib(_X_, R_ClassSymbol)
#define rc_getLevels(_X_) Rf_getAttrib(_X_, R_LevelsSymbol)

bool rc_isS4Null(SEXP obj);

#ifdef __cplusplus
}
#endif

#endif // RC_UTIL_H

