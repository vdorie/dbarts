#include <rc/util.h>

#include <errno.h>
#include <stdarg.h>
#include <misc/stddef.h> // size_t
#include <stdlib.h> // atoi
#include <string.h> // strncmp

SEXP rc_setDims(SEXP obj, ...)
{
  va_list dimsPointer;
  va_start(dimsPointer, obj);
  int dim = va_arg(dimsPointer, int);
  size_t numDims = 0;
  while (dim >= 0) {
    ++numDims;
    dim = va_arg(dimsPointer, int);
  }
  va_end(dimsPointer);
  
  if (numDims == 0) return obj;
  
  SEXP dimsExpr = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) numDims));
  int* dims = INTEGER(dimsExpr);
  va_start(dimsPointer, obj);
  for (size_t i = 0; i < numDims; ++i)
    dims[i] = va_arg(dimsPointer, int);
  va_end(dimsPointer);
  
  R_do_slot_assign(obj, R_DimSymbol, dimsExpr);
  
  UNPROTECT(1);
  
  return obj;
}

SEXP rc_allocateInSlot(SEXP obj, SEXP slotName, SEXPTYPE type, R_xlen_t length)
{
  SEXP val = PROTECT(Rf_allocVector(type, length));
   
  R_do_slot_assign(obj, slotName, val);
  
  UNPROTECT(1);
  
  return val;
}

SEXP rc_getListElement(SEXP list, const char* name)
{
  SEXP names = PROTECT(rc_getNames(list));
  if (Rf_isNull(names)) {
    UNPROTECT(1);
    return R_NilValue;
  }
  
  SEXP result = R_NilValue;
  size_t len = rc_getLength(names);
  for (size_t i = 0; i < len; ++i) {
    if (strcmp(name, CHAR(STRING_ELT(names, i))) == 0) {
      result = VECTOR_ELT(list, i);
      break;
    }
  }
  
  UNPROTECT(1);
  return result;
}

bool rc_isS4Null(SEXP obj)
{
  if (!Rf_isSymbol(obj)) return false;
  
  const char* symbolName = CHAR(PRINTNAME(obj));
  
  if (strncmp(symbolName, "\1NULL\1", 6) == 0) return true;
    
  return false;
}

static char* rc_strdup(const char* c) {
  size_t len = strlen(c);
  char* result = (char*) malloc(len + 1);
  if (result != NULL) memcpy(result, c, len + 1);
  return result;
}

int rc_getRuntimeVersion(int* major, int* minor, int* revision)
{
  *major = -1, *minor = -1, *revision = -1;
  SEXP versionFunction = PROTECT(Rf_findVarInFrame(R_BaseNamespace, Rf_install("R.Version")));
  if (versionFunction == R_UnboundValue) {
    UNPROTECT(1);
    return ENXIO;
  }
  SEXP closure = PROTECT(Rf_lang1(versionFunction));
  SEXP version = PROTECT(Rf_eval(closure, R_GlobalEnv));
  if (Rf_isNull(version)) {
    UNPROTECT(3);
    return ENOSYS;
  }
  
  R_xlen_t versionLength = XLENGTH(version);
  SEXP versionNames = rc_getNames(version);
  for (R_xlen_t i = 0; i < versionLength; ++i) {
    if (strcmp(CHAR(STRING_ELT(versionNames, i)), "major") == 0) {
      *major = atoi(CHAR(STRING_ELT(VECTOR_ELT(version, i), 0)));
    } else if (strcmp(CHAR(STRING_ELT(versionNames, i)), "minor") == 0) {
      char* minorString = rc_strdup(CHAR(STRING_ELT(VECTOR_ELT(version, i), 0)));
      int period = 0;
      for ( ; minorString[period] != '.' && minorString[period] != '\0'; ++period) ;
      if (minorString[period] == '.') {
        minorString[period] = '\0';
        *minor = atoi(minorString);
        if (minorString[period + 1] != '\0')
          *revision = atoi(minorString + period + 1);
      } else {
        *minor = atoi(minorString);
        *revision = 0;
      }
      free(minorString);
    }
  }
  UNPROTECT(3);
  
  return (*major >= 0 && *minor >= 0 && *revision >= 0) ? 0 : EPROTO;
}
