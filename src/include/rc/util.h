#ifndef RC_UTIL_H
#define RC_UTIL_H

#include <Rinternals.h> // SEXP

#ifdef __cplusplus
extern "C" {
#endif

SEXP rc_setDims(SEXP obj, ...);
SEXP rc_allocateInSlot(SEXP obj, SEXP slotName, SEXPTYPE type, R_xlen_t length);

#ifdef __cplusplus
}
#endif

#endif // RC_UTIL_H
