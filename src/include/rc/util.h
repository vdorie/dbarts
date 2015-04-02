#ifndef RC_UTIL_H
#define RC_UTIL_H

#include <Rinternals.h> // SEXP
#include <external/stddef.h> // ext_size_t

#ifdef __cplusplus
extern "C" {
#endif

SEXP rc_setDims(SEXP obj, ...);
SEXP rc_allocateInSlot(SEXP obj, SEXP slotName, SEXPTYPE type, R_xlen_t length);
#ifdef __cplusplus
#define rc_getLength(_X_) static_cast<ext_size_t>(XLENGTH(_X_))
#else
#define rc_getLength(_X_) ((ext_size_t) XLENGTH(_X_))
#endif

#ifdef __cplusplus
}
#endif

#endif // RC_UTIL_H
