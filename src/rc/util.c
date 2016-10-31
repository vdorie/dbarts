#include <rc/util.h>

#include <stdarg.h>
#include <stddef.h> // size_t
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
  
  SEXP dimsExpr = Rf_allocVector(INTSXP, (R_xlen_t) numDims);
  int* dims = INTEGER(dimsExpr);
  va_start(dimsPointer, obj);
  for (size_t i = 0; i < numDims; ++i)
    dims[i] = va_arg(dimsPointer, int);
  va_end(dimsPointer);
  
  R_do_slot_assign(obj, R_DimSymbol, dimsExpr);
  
  return obj;
}

SEXP rc_allocateInSlot(SEXP obj, SEXP slotName, SEXPTYPE type, R_xlen_t length)
{
  SEXP val = Rf_allocVector(type, length);
   
  R_do_slot_assign(obj, slotName, val);
  return val;
}

bool rc_isS4Null(SEXP obj)
{
  if (!Rf_isSymbol(obj)) return false;
  
  const char* symbolName = CHAR(PRINTNAME(obj));
  
  if (strncmp(symbolName, "\1NULL\1", 6) == 0) return true;
    
  return false;
}

