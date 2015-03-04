#include <external/rc.h>

#include <limits.h> // INT_MAX
#include <math.h>   // isnan
#include <stdarg.h>
#include <stddef.h> // size_t

#include <Rdefines.h>

#ifdef _WIN32
#define SIZE_T_FMT "%Iu"
#else
#define SIZE_T_FMT "%zu"
#endif

#if R_XLEN_T_MAX < INT_MAX
#define va_arg_xlen_t(_ARGS_) (R_xlen_t) va_arg(_ARGS_, int)
#else 
#define va_arg_xlen_t(_ARGS_)            va_arg(_ARGS_, R_xlen_t)
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

static void checkLengthConstraint(const char* name, _rc_boundType boundType, R_xlen_t length, R_xlen_t bound);
static void checkIntConstraint(const char* name, _rc_boundType boundType, int value, int bound);
static void checkDoubleConstraint(const char* name, _rc_boundType boundType, double value, double bound);
static void checkBoolConstraint(const char* name, _rc_boundType boundType, bool value, bool bound);

int rc_getInt(SEXP x, const char* name, ...)
{
  if (!IS_INTEGER(x)) error("%s must be of type integer", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          checkLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          rc_naAllowableType naOK = BOUND(arg);
          if (naOK == RC_NO) error("%s cannot be of length 0 if NA is not allowable", name);
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    return NA_INTEGER;
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
        checkLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        int valueBound = va_arg(argsPointer, int);
        checkIntConstraint(name, BOUND(arg), result, valueBound);
      }
      break;
      case RC_NA:
      {
        rc_naAllowableType naOK = BOUND(arg);
        if (result == NA_INTEGER && naOK == RC_NO) error("%s cannot be NA", name);
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  return result;
}

void rc_checkInts(SEXP x, const char* name, ...)
{
  if (!IS_INTEGER(x)) error("%s must be of type integer", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          checkLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          rc_naAllowableType naOK = BOUND(arg);
          if (naOK == RC_NO) error("%s cannot be of length 0 if NA is not allowable", name);
        }
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
  
  int* results = INTEGER(x);

  va_start(argsPointer, name);
  
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        checkLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        _rc_boundType boundType = BOUND(arg);
        int valueBound = va_arg(argsPointer, int);
        for (size_t i = 0; i < length; ++i) checkIntConstraint(name, boundType, results[i], valueBound);
      }
      break;
      case RC_NA:
      {
        rc_naAllowableType naOK = BOUND(arg);
        for (size_t i = 0; i < length; ++i)
          if (results[i] == NA_INTEGER && naOK == RC_NO) error("%s cannot be NA", name);
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
}

int rc_getDouble(SEXP x, const char* name, ...)
{
  if (!IS_REAL(x)) error("%s must be of type real", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          checkLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, double);
        break;
        case RC_NA:
        {
          rc_naAllowableType naOK = BOUND(arg);
          if (naOK == RC_NO) error("%s cannot be of length 0 if NA is not allowable", name);
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    return NA_REAL;
  }
  
  double result = INTEGER(x)[0];

  va_start(argsPointer, name);
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        checkLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        double valueBound = va_arg(argsPointer, double);
        checkDoubleConstraint(name, BOUND(arg), result, valueBound);
      }
      break;
      case RC_NA:
      {
        rc_naAllowableType naOK = BOUND(arg);
        if ((isnan(result) || result == NA_REAL) && naOK == RC_NO) error("%s cannot be NA", name);
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  return result;
}

void rc_checkDoubles(SEXP x, const char* name, ...)
{
  if (!IS_REAL(x)) error("%s must be of type real", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          checkLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, double);
        break;
        case RC_NA:
        {
          rc_naAllowableType naOK = BOUND(arg);
          if (naOK == RC_NO) error("%s cannot be of length 0 if NA is not allowable", name);
        }
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
  
  double* results = REAL(x);

  va_start(argsPointer, name);
  
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        checkLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        _rc_boundType boundType = BOUND(arg);
        double valueBound = va_arg(argsPointer, double);
        for (size_t i = 0; i < length; ++i) checkDoubleConstraint(name, boundType, results[i], valueBound);
      }
      break;
      case RC_NA:
      {
        rc_naAllowableType naOK = BOUND(arg);
        for (size_t i = 0; i < length; ++i)
          if ((isnan(results[i]) || results[i] == NA_REAL) && naOK == RC_NO) error("%s cannot be NA", name);
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
}

bool rc_getBool(SEXP x, const char* name, ...)
{
  if (!IS_LOGICAL(x)) error("%s must be of type logical", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          checkLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          rc_naAllowableType naOK = BOUND(arg);
          if (naOK == RC_NO) error("%s cannot be of length 0 if NA is not allowable", name);
        }
        break;
        default:
        break;
      }
      arg = va_arg(argsPointer, int);
      constraintType = TYPE(arg);
    }
    va_end(argsPointer);
    return NA_LOGICAL;
  }
  
  bool result = LOGICAL(x)[0];

  va_start(argsPointer, name);
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        checkLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        int valueBound = va_arg(argsPointer, int);
        checkLogicalConstraint(name, BOUND(arg), result, valueBound);
      }
      break;
      case RC_NA:
      {
        rc_naAllowableType naOK = BOUND(arg);
        if (result == NA_LOGICAL && naOK == RC_NO) error("%s cannot be NA", name);
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
  
  return result;
}

void rc_checkBools(SEXP x, const char* name, ...)
{
  if (!IS_LOGICAL(x)) error("%s must be of type logical", name);
  
  R_xlen_t length = XLENGTH(x);
  va_list argsPointer;
  rc_constraintType constraintType;
  int arg;
  
  if (length == 0) {
    va_start(argsPointer, name);
    
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
    while (constraintType != RC_END) {
      switch (constraintType) {
        case RC_LENGTH:
        {
          R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
          checkLengthConstraint(name, BOUND(arg), length, lengthBound);
        }
        break;
        case RC_VALUE:
        va_arg(argsPointer, int);
        break;
        case RC_NA:
        {
          rc_naAllowableType naOK = BOUND(arg);
          if (naOK == RC_NO) error("%s cannot be of length 0 if NA is not allowable", name);
        }
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
  
  int* results = LOGICAL(x);

  va_start(argsPointer, name);
  
  arg = va_arg(argsPointer, int);
  constraintType = TYPE(arg);
  while (constraintType != RC_END) {
    switch (constraintType) {
      case RC_LENGTH:
      {
        R_xlen_t lengthBound = va_arg_xlen_t(argsPointer);
        checkLengthConstraint(name, BOUND(arg), length, lengthBound);
      }
      break;
      case RC_VALUE:
      {
        _rc_boundType boundType = BOUND(arg);
        int valueBound = va_arg(argsPointer, int);
        for (size_t i = 0; i < length; ++i) checkBoolConstraint(name, boundType, results[i], valueBound);
      }
      break;
      case RC_NA:
      {
        rc_naAllowableType naOK = BOUND(arg);
        for (size_t i = 0; i < length; ++i)
          if (results[i] == NA_INTEGER && naOK == RC_NO) error("%s cannot be NA", name);
      }
      default:
      break;
    }
    arg = va_arg(argsPointer, int);
    constraintType = TYPE(arg);
  }
  va_end(argsPointer);
}

static void checkLengthConstraint(const char* name, _rc_boundType boundType, R_xlen_t length, R_xlen_t bound)
{
  if (bound < 0) error("internal error: %s cannot have a negative length", name);
  
  if (length == 0) {
    switch (boundType) {
      case _RC_GT:
      error("%s must be of length greater than " SIZE_T_FMT, name, bound);
      break;
      case _RC_GEQ:
      if (bound > 0) error("%s must be of length greater than or equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_LT:
      if (bound == 0) error("internal error: %s cannot be of length less than 0", name);
      break;
      case _RC_EQ:
      if (bound != 0) error("%s must be of length equal to 0", name);
      break;
      case _RC_NE:
      if (bound == 0) error("%s cannot be length equal to 0", name);
      break;
      default:
      break;
    }
  } else {
    switch (boundType) {
      case _RC_GT:
      if (length <= bound) error("%s must be of length greater than " SIZE_T_FMT, name, bound);
      break;
      case _RC_GEQ:
      if (length < bound) error("%s must be of length greater than or equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_LT:
      if (length >= bound) error("%s must be of length less than " SIZE_T_FMT, name, bound);
      break;
      case _RC_LEQ:
      if (length > bound) error("%s must be of length less than or equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_EQ:
      if (length != bound) error("%s must be of length equal to " SIZE_T_FMT, name, bound);
      break;
      case _RC_NE:
      if (length == bound) error("%s cannot be of length equal to " SIZE_T_FMT, name, bound);
      break;
      default:
      break;
    }
  }
}

static void checkIntConstraint(const char* name, _rc_boundType boundType, int value, int bound)
{
  if (bound == NA_INTEGER) error("bound for %s cannot be NA", name);
  if (value == NA_INTEGER) return;
  
  switch (boundType) {
    case _RC_GT:
    if (value <= bound) error("%s must be greater than %d", name, bound);
    break;
    case _RC_GEQ:
    if (value < bound) error("%s must be greater than or equal to %d", name, bound);
    break;
    case _RC_LT:
    if (value >= bound) error("%s must be less than %d", name, bound);
    break;
    case _RC_LEQ:
    if (value > bound) error("%s must be less than or equal to %d", name, bound);
    break;
    case _RC_EQ:
    if (value != bound) error("%s must be equal to %d", name, bound);
    break;
    case _RC_NE:
    if (value == bound) error("%s cannot equal %d", name, bound);
    break;
    default:
    break;
  }
}

static void checkDoubleConstraint(const char* name, _rc_boundType boundType, double value, double bound)
{
  if (isnan(bound)) error("bound for %s cannot be NaN", name);
  if (bound == NA_REAL) error("bound for %s cannot be NA", name);
  if (isnan(value)) error("%s is NaN", name);
  if (value == NA_REAL) return;
  
  switch (boundType) {
    case _RC_GT:
    if (bound == R_PosInf) error("%s cannot be greater than positive infinity", name);
    if (bound == R_NegInf && value == R_NegInf) error("for %s, cannot compare negative infinities", name);
    if (bound != R_NegInf && value <= bound) error("%s must be greater than %d", name, bound);
    break;
    case _RC_GEQ:
    if (bound == R_PosInf && value != R_PosInf) error("%s must be equal to positive infinity", name);
    if (bound != R_NegInf && value < bound) error("%s must be greater than or equal to %d", name, bound);
    break;
    case _RC_LT:
    if (bound == R_NegInf) error("%s cannot be less than negative infinity", name);
    if (bound == R_PosInf && value == R_PosInf) error("for %s, cannot compare positive infinites", name);
    if (bound != R_PosInf && value >= bound) error("%s must be less than %d", name, bound);
    break;
    case _RC_LEQ:
    if (bound == R_NegInf && value != R_NegInf) error("%s must be equal to negative infinity", name);
    if (bound == R_PosInf && value > bound) error("%s must be less than or equal to %d", name, bound);
    break;
    case _RC_EQ:
    if (value != bound) error("%s must be equal to %d", name, bound);
    break;
    case _RC_NE:
    if (value == bound) error("%s cannot equal %d", name, bound);
    default:
    break;
  }
}

static void checkBoolConstraint(const char* name, _rc_boundType boundType, bool value, bool bound)
{
  if (bound == NA_LOGICAL) error("bound for %s cannot be NA", name);
  if (value == NA_LOGICAL) return;
  
  switch (boundType) {
    case _RC_GT:
    case _RC_GEQ:
    case _RC_LT:
    case _RC_LEQ:
    error("for %s, logicals cannot be ordered", name);
    break;
    case _RC_EQ:
    if (value != bound) error("%s must be equal to %s", name, bound ? "true" : "false");
    break;
    case _RC_NE:
    if (value == bound) error("%s cannot equal %s", name, bound ? "true" : "false");
    break;
    default:
    break;
  }
}

