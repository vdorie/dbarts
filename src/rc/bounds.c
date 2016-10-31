#include <rc/bounds.h>

#include <limits.h> // INT_MAX
#include <math.h>   // isnan
#include <stdarg.h>
#include <stddef.h> // size_t

#include <rc/util.h>

#ifdef _WIN32
#  define SIZE_T_FMT "%Iu"
#else
#  define SIZE_T_FMT "%zu"
#endif

#if R_XLEN_T_MAX < INT_MAX
#  define va_arg_xlen_t(_ARGS_) (R_xlen_t) va_arg(_ARGS_, int)
#else 
#  define va_arg_xlen_t(_ARGS_)            va_arg(_ARGS_, R_xlen_t)
#endif

#define TYPE(_ARG_) (_ARG_ & 0x03)
#define BOUND(_ARG_) ((_ARG_ & 0x1C) >> 2)

typedef enum {
  _RC_GT  = 0x01,
  _RC_LT  = 0x02,
  _RC_GEQ = 0x03,
  _RC_LEQ = 0x04,
  _RC_EQ  = 0x05,
  _RC_NE  = 0x06
} _rc_boundType;

typedef enum {
  _RC_YES = 0x01,
  _RC_NO  = 0x02
} _rc_naAllowableType;

static void assertLengthConstraint(const char* name, _rc_boundType boundType, R_xlen_t length, R_xlen_t bound);
static void assertIntConstraint(const char* name, _rc_boundType boundType, int value, int bound);
static void assertDoubleConstraint(const char* name, _rc_boundType boundType, double value, double bound);
static void assertBoolConstraint(const char* name, _rc_boundType boundType, bool value, bool bound);

int rc_getInt(SEXP x, const char* name, ...)
{
  if (!Rf_isInteger(x)) Rf_error("%s must be of type integer", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  _rc_naAllowableType naOK = _RC_NO;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          assertLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          naOK = BOUND(arg);
          if (naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be of length 0 if NA is not allowable", name); }
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    
    if (naOK == _RC_NO) Rf_error("%s cannot be of length 0 if NA is not allowable", name);
    
    return R_NaInt;
  }
  
  int result = INTEGER(x)[0];

  va_start(argsPointer, name);
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        assertLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        int valueBound = va_arg(argsPointer, int);
        assertIntConstraint(name, BOUND(arg), result, valueBound);
      }
      break;
      case RC_NA:
      {
        naOK = BOUND(arg);
        if (result == R_NaInt && naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be NA", name); }
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  if (result == R_NaInt && naOK == _RC_NO) Rf_error("%s cannot be NA", name);
  
  return result;
}

void rc_assertIntConstraints(SEXP x, const char* name, ...)
{
  if (!Rf_isInteger(x)) Rf_error("%s must be of type integer", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  _rc_naAllowableType naOK = _RC_NO;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          assertLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          naOK = BOUND(arg);
          if (naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be of length 0 if NA is not allowable", name); }
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    
    if (naOK == _RC_NO) Rf_error("%s cannot be of length 0 if NA is not allowable", name);
    
    return;
  }
  
  int* results = INTEGER(x);

  va_start(argsPointer, name);
  
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        assertLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        _rc_boundType boundType = BOUND(arg);
        int valueBound = va_arg(argsPointer, int);
        for (size_t i = 0; i < (size_t) length; ++i) assertIntConstraint(name, boundType, results[i], valueBound);
      }
      break;
      case RC_NA:
      {
        naOK = BOUND(arg);
        for (size_t i = 0; i < (size_t) length; ++i)
          if (results[i] == R_NaInt && naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be NA", name); }
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  if (naOK == _RC_NO) {
    for (size_t i = 0; i < (size_t) length; ++i)
      if (results[i] == R_NaInt) Rf_error("%s cannot be NA", name);
  }
}

double rc_getDouble(SEXP x, const char* name, ...)
{
  if (!Rf_isReal(x)) Rf_error("%s must be of type real", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  _rc_naAllowableType naOK = _RC_NO;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          assertLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, double);
        break;
        case RC_NA:
        {
          naOK = BOUND(arg);
          if (naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be of length 0 if NA is not allowable", name); }
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    
    if (naOK == _RC_NO) Rf_error("%s cannot be of length 0 if NA is not allowable", name);
    
    return R_NaReal;
  }
  
  double result = REAL(x)[0];

  va_start(argsPointer, name);
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        assertLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        double valueBound = va_arg(argsPointer, double);
        assertDoubleConstraint(name, BOUND(arg), result, valueBound);
      }
      break;
      case RC_NA:
      {
        naOK = BOUND(arg);
        if ((isnan(result) || result == R_NaReal) && naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be NA", name); }
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  if ((isnan(result) || result == R_NaReal) && naOK == _RC_NO) Rf_error("%s cannot be NA", name);
  
  return result;
}

void rc_assertDoubleConstraints(SEXP x, const char* name, ...)
{
  if (!Rf_isReal(x)) Rf_error("%s must be of type real", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  _rc_naAllowableType naOK = _RC_NO;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          assertLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, double);
        break;
        case RC_NA:
        {
          naOK = BOUND(arg);
          if (naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be of length 0 if NA is not allowable", name); }
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    
    if (naOK == _RC_NO) Rf_error("%s cannot be of length 0 if NA is not allowable", name);
    
    return;
  }
  
  double* results = REAL(x);

  va_start(argsPointer, name);
  
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        assertLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        _rc_boundType boundType = BOUND(arg);
        double valueBound = va_arg(argsPointer, double);
        for (size_t i = 0; i < (size_t) length; ++i) assertDoubleConstraint(name, boundType, results[i], valueBound);
      }
      break;
      case RC_NA:
      {
        naOK = BOUND(arg);
        for (size_t i = 0; i < (size_t) length; ++i)
          if ((isnan(results[i]) || results[i] == R_NaReal) && naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be NA", name); }
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  if (naOK == _RC_NO) {
    for (size_t i = 0; i < (size_t) length; ++i)
      if (isnan(results[i]) || results[i] == R_NaReal) Rf_error("%s cannot be NA", name);
  }
}

bool rc_getBool(SEXP x, const char* name, ...)
{
  if (!Rf_isLogical(x)) Rf_error("%s must be of type logical", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  _rc_naAllowableType naOK = _RC_NO;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          assertLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          naOK = BOUND(arg);
          if (naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be of length 0 if NA is not allowable", name); }
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    
    if (naOK == _RC_NO) Rf_error("%s cannot be of length 0 if NA is not allowable", name);
    
    return NA_LOGICAL;
  }
  
  int result = LOGICAL(x)[0];

  va_start(argsPointer, name);
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        assertLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        int valueBound = va_arg(argsPointer, int);
        assertBoolConstraint(name, BOUND(arg), result, valueBound);
      }
      break;
      case RC_NA:
      {
        naOK = BOUND(arg);
        if (result == R_NaInt && naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be NA", name); }
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  if (result == R_NaInt && naOK == _RC_NO) Rf_error("%s cannot be NA", name);
    
  return result;
}

void rc_assertBoolConstraints(SEXP x, const char* name, ...)
{
  if (!Rf_isLogical(x)) Rf_error("%s must be of type logical", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  _rc_naAllowableType naOK = _RC_NO;
  
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          assertLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          naOK = BOUND(arg);
          if (naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be of length 0 if NA is not allowable", name); }
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    
    if (naOK == _RC_NO) Rf_error("%s cannot be of length 0 if NA is not allowable", name);
    
    return;
  }
  
  int* results = LOGICAL(x);

  va_start(argsPointer, name);
  
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        assertLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        _rc_boundType boundType = BOUND(arg);
        int valueBound = va_arg(argsPointer, int);
        for (size_t i = 0; i < (size_t) length; ++i) assertBoolConstraint(name, boundType, results[i], valueBound);
      }
      break;
      case RC_NA:
      {
        naOK = BOUND(arg);
        for (size_t i = 0; i < (size_t) length; ++i)
          if (results[i] == R_NaInt && naOK == _RC_NO) { va_end(argsPointer); Rf_error("%s cannot be NA", name); }
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  if (naOK == _RC_NO) {
    for (size_t i = 0; i < (size_t) length; ++i)
      if (results[i] == R_NaInt) Rf_error("%s cannot be NA", name);
  }
}

void rc_assertDimConstraints(SEXP x, const char* name, ...)
{
  SEXP dimsExpr = Rf_getAttrib(x, R_DimSymbol);
  
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  R_xlen_t lengthBound = -1;
  
  R_xlen_t length;
  
  if (Rf_isNull(dimsExpr) || XLENGTH(dimsExpr) == 0) {
    length = 0;
    
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        lengthBound = va_arg_xlen_t(argsPointer);
        assertLengthConstraint(name, BOUND(arg), length, lengthBound);
        break;
        
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        
        case RC_NA:
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    return;
  }
  
  int* dims = INTEGER(dimsExpr);
  length = XLENGTH(dimsExpr);
  R_xlen_t pos = 0;
  
  va_start(argsPointer, name);
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      lengthBound = va_arg_xlen_t(argsPointer);
      assertLengthConstraint(name, BOUND(arg), length, lengthBound);
      break;
      
      case RC_VALUE:
      {
        int valueBound = va_arg(argsPointer, int);
        if (pos < length)
          assertIntConstraint(name, BOUND(arg), dims[pos], valueBound);
        ++pos;
      }
      break;
      
      case RC_NA:
      ++pos;
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  if (lengthBound == -1 || pos <= length) return;
  
  // If a length bound is specified, we may have to check for constraints beyond
  // the length of the dimensions, e.g.
  // 
  //   rc_checkDims(x, n, RC_NA, RC_NA, RC_VALUE | RC_EQ, 42, RC_LENGTH | RC_LEQ, 10, RC_END);
  // 
  // on dims
  //
  //   c(5, 10)
  
  pos = 0;
  va_start(argsPointer, name);
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      va_arg_xlen_t(argsPointer); // ignore
      break;
      
      case RC_VALUE:
      if (pos > length) { va_end(argsPointer); Rf_error("%s too short to satisfy all constraints", name); }
      va_arg(argsPointer, int); // ignore
      ++pos;
      break;
      
      case RC_NA:
      ++pos;
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
}

static void assertLengthConstraint(const char* name, _rc_boundType boundType, R_xlen_t length, R_xlen_t bound)
{
  if (bound < 0) Rf_error("internal error: %s cannot have a negative length", name);
  
  
  if (length == 0) {
    switch (boundType) {
      case _RC_GT:
      Rf_error("%s must be of length greater than " SIZE_T_FMT, name, bound);
      break;
      case _RC_GEQ:
      if (bound > 0) Rf_error("%s must be of length greater than or equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_LT:
      if (bound == 0) Rf_error("internal error: %s cannot be of length less than 0", name);
      break;
      case _RC_EQ:
      if (bound != 0) Rf_error("%s must be of length equal to 0", name);
      break;
      case _RC_NE:
      if (bound == 0) Rf_error("%s cannot be length equal to 0", name);
      break;
      default:
      break;
    }
  } else {
    switch (boundType) {
      case _RC_GT:
      if (length <= bound) Rf_error("%s must be of length greater than " SIZE_T_FMT, name, bound);
      break;
      case _RC_GEQ:
      if (length < bound) Rf_error("%s must be of length greater than or equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_LT:
      if (length >= bound) Rf_error("%s must be of length less than " SIZE_T_FMT, name, bound);
      break;
      case _RC_LEQ:
      if (length > bound) Rf_error("%s must be of length less than or equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_EQ:
      if (length != bound) Rf_error("%s must be of length equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_NE:
      if (length == bound) Rf_error("%s cannot be of length equal to " SIZE_T_FMT, name, bound);
      break;
      default:
      break;
    }
  }
}

static void assertIntConstraint(const char* name, _rc_boundType boundType, int value, int bound)
{
  if (bound == R_NaInt) Rf_error("bound for %s cannot be NA", name);
  if (value == R_NaInt) return;
  
  switch (boundType) {
    case _RC_GT:
    if (value <= bound) Rf_error("%s must be greater than %d", name, bound);
    break;
    case _RC_LT:
    if (value >= bound) Rf_error("%s must be less than %d", name, bound);
    break;
    case _RC_GEQ:
    if (value < bound) Rf_error("%s must be greater than or equal to %d", name, bound);
    break;
    case _RC_LEQ:
    if (value > bound) Rf_error("%s must be less than or equal to %d", name, bound);
    break;
    case _RC_EQ:
    if (value != bound) Rf_error("%s must be equal to %d", name, bound);
    break;
    case _RC_NE:
    if (value == bound) Rf_error("%s cannot equal %d", name, bound);
    break;
    default:
    break;
  }
}

static void assertDoubleConstraint(const char* name, _rc_boundType boundType, double value, double bound)
{
  if (isnan(bound)) Rf_error("bound for %s cannot be NaN", name);
  if (bound == R_NaReal) Rf_error("bound for %s cannot be NA", name);
  if (isnan(value)) Rf_error("%s is NaN", name);
  if (value == R_NaReal) return;
  
  switch (boundType) {
    case _RC_GT:
    if (bound == R_PosInf) Rf_error("%s cannot be greater than positive infinity", name);
    if (bound == R_NegInf && value == R_NegInf) Rf_error("for %s, cannot compare negative infinities", name);
    if (bound != R_NegInf && value <= bound) Rf_error("%s must be greater than %f", name, bound);
    break;
    case _RC_GEQ:
    if (bound == R_PosInf && value != R_PosInf) Rf_error("%s must be equal to positive infinity", name);
    if (bound != R_NegInf && value < bound) Rf_error("%s must be greater than or equal to %f", name, bound);
    break;
    case _RC_LT:
    if (bound == R_NegInf) Rf_error("%s cannot be less than negative infinity", name);
    if (bound == R_PosInf && value == R_PosInf) Rf_error("for %s, cannot compare positive infinites", name);
    if (bound != R_PosInf && value >= bound) Rf_error("%s must be less than %f", name, bound);
    break;
    case _RC_LEQ:
    if (bound == R_NegInf && value != R_NegInf) Rf_error("%s must be equal to negative infinity", name);
    if (bound == R_PosInf && value > bound) Rf_error("%s must be less than or equal to %f", name, bound);
    break;
    case _RC_EQ:
    if (value != bound) Rf_error("%s must be equal to %f", name, bound);
    break;
    case _RC_NE:
    if (value == bound) Rf_error("%s cannot equal %f", name, bound);
    default:
    break;
  }
}

static void assertBoolConstraint(const char* name, _rc_boundType boundType, bool value, bool bound)
{
  if (bound == R_NaInt) Rf_error("bound for %s cannot be NA", name);
  if (value == R_NaInt) return;
  
  switch (boundType) {
    case _RC_GT:
    case _RC_GEQ:
    case _RC_LT:
    case _RC_LEQ:
    Rf_error("for %s, logicals cannot be ordered", name);
    break;
    case _RC_EQ:
    if (value != bound) Rf_error("%s must be equal to %s", name, bound ? "true" : "false");
    break;
    case _RC_NE:
    if (value == bound) Rf_error("%s cannot equal %s", name, bound ? "true" : "false");
    break;
    default:
    break;
  }
}

