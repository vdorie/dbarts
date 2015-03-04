#ifndef EXTERNAL_RC_H
#define EXTERNAL_RC_H

#include <stdbool.h> // for bool

#include <Rinternals.h> // SEXP

#ifdef __cplusplus
extern "C" {
#endif

// length constraints are given as R_xlen_t, value constraints in the type
// of the given function
int rc_getInt(SEXP x, const char* name, ...);
double rc_getDouble(SEXP x, const char* name, ...);
bool rc_getBool(SEXP x, const char* name, ...);


void rc_checkInts(SEXP x, const char* name, ...);
void rc_checkDoubles(SEXP x, const char* name, ...);
void rc_checkBools(SEXP x, const char* name, ...);

typedef enum {
  RC_END = 0x0,
  RC_LENGTH = 0x1,
  RC_VALUE = 0x2,
  RC_NA = 0x3
} rc_constraintType;

typedef enum {
  RC_GT  = 0x04,
  RC_LT  = 0x08,
  RC_GEQ = 0x0C,
  RC_LEQ = 0x10,
  RC_EQ  = 0x14,
  RC_NE  = 0x18
} rc_boundType;

typedef enum {
  RC_YES = 0x04,
  RC_NO = 0x08
} rc_naAllowableType;

#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_RC_H
