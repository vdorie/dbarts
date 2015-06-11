#ifndef RC_BOUNDS_H
#define RC_BOUNDS_H

#include <stdbool.h> // for bool

#define R_NO_REMAP 1
#include <R.h>
#include <Rinternals.h> // SEXP

#ifdef __cplusplus
extern "C" {
#endif

// Length constraints are given as R_xlen_t; it is up to the caller to guarantee
// that the argument is of sufficient size. On x32, R_xlen_ts are ints and unlabeled
// literals will suffice. On x86_64, they are ptrdiff_ts, which are the same as
// int64_t. The safest thing to do is cast them to R_xlen_t when passing.
//
// Value constraints are in the type of the given function, i.e. int, double, and
// int again for booleans.
//
// Constraint types are to be OR'd with the bound type, and the value takes the next
// argument. e.g.
//
//   rc_getInt(x, "var name", RC_LENGTH | RC_GT, (R_xlen_t) 0, RC_END)
//
// NA constraint types only accept YES/NO, so that this info is passed as a bound. e.g.
//
//   rc_getBool(x, "var name", RC_NA | RC_YES, RC_END);
//
// NA_NO is the default

int rc_getInt(SEXP x, const char* name, ...);
double rc_getDouble(SEXP x, const char* name, ...);
bool rc_getBool(SEXP x, const char* name, ...);

void rc_checkInts(SEXP x, const char* name, ...);
void rc_checkDoubles(SEXP x, const char* name, ...);
void rc_checkBools(SEXP x, const char* name, ...);

// Value constraints are matched in-order to dims themselves. A lack of length
// constraint allows missing dims. RC_NA skips the dim, regardless of the value
// (equiv to value >= 0). e.g.
//
//  rc_checkDims(x, "var name", RC_LENGTH | RC_EQ, (R_xlen_t) 2,
//               RC_NA,
//               RC_VALUE | RC_EQ, (int) numCols,
//               RC_END);
void rc_checkDims(SEXP x, const char* name, ...);

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

#endif // RC_BOUNDS_H
