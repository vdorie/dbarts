#ifndef RC_BOUNDS_H
#define RC_BOUNDS_H

#include <stdbool.h> // for bool

#include <misc/stddef.h> // misc_size_t

#include <external/Rinternals.h> // SEXP, R_NilValue

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

int rc_getInt0(SEXP x, const char* name);
int rc_getInt(SEXP x, const char* name, ...);
int rc_getIntAt(SEXP v, misc_size_t i, const char* name, ...);
double rc_getDouble0(SEXP x, const char* name);
double rc_getDouble(SEXP x, const char* name, ...);
double rc_getDoubleAt(SEXP v, misc_size_t i, const char* name, ...);
bool rc_getBool0(SEXP x, const char* name);
bool rc_getBool(SEXP x, const char* name, ...);
bool rc_getBoolAt(SEXP v, misc_size_t i, const char* name, ...);

void rc_assertIntConstraints(SEXP x, const char* name, ...);
void rc_assertDoubleConstraints(SEXP x, const char* name, ...);
void rc_assertBoolConstraints(SEXP x, const char* name, ...);

// TODO: write functions that don't throw errors, like ...
// const char* rc_checkInts(SEXP x, ...);

// Value constraints are matched in-order to dims themselves. A lack of length
// constraint allows missing dims. RC_NA skips the dim, regardless of the value
// (equiv to value >= 0). e.g.
//
//  rc_assertDimConstraints(x, "var name",
//                          RC_LENGTH | RC_EQ, (R_xlen_t) 2,
//                          RC_NA,
//                          RC_VALUE | RC_EQ, (int) numCols,
//                          RC_END);
void rc_assertDimConstraints(SEXP x, const char* name, ...);

// if standard is > C++17 (i.e. C++20), bit-wise operations on enums has been deprecated
// we can use constexprs since C++11, but are only doing so if the compiler would complain 
#if defined(__cplusplus) && (__cplusplus > 201703L || (__cplusplus >= 201103L && defined(__clang__) && __clang_major__ >= 11))

typedef unsigned char rc_constraintType;
constexpr unsigned char RC_END    = 0x0;
constexpr unsigned char RC_LENGTH = 0x1;
constexpr unsigned char RC_VALUE  = 0x2;
constexpr unsigned char RC_NA     = 0x3;

typedef unsigned char rc_boundType;
constexpr unsigned char RC_GT  = 0x04;
constexpr unsigned char RC_LT  = 0x08;
constexpr unsigned char RC_GEQ = 0x0C;
constexpr unsigned char RC_LEQ = 0x10;
constexpr unsigned char RC_EQ  = 0x14;
constexpr unsigned char RC_NE  = 0x18;
constexpr unsigned char RC_DEFAULT = 0x1C;

typedef unsigned char rc_naAllowableType;
constexpr unsigned char RC_YES = 0x04;
constexpr unsigned char RC_NO  = 0x08;

#else

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
  RC_NE  = 0x18,
  RC_DEFAULT = 0x1C
} rc_boundType;

typedef enum {
  RC_YES = 0x04,
  RC_NO = 0x08
} rc_naAllowableType;
#endif

#ifdef __cplusplus
}
#endif

#endif // RC_BOUNDS_H

