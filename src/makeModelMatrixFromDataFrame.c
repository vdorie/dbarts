#include "makeModelMatrixFromDataFrame.h"

#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>

#include <external/alloca.h>
#include <external/linearAlgebra.h>

#include <rc/util.h>

typedef enum {
  REAL_VECTOR = 0,
  REAL_MATRIX,
  INTEGER_VECTOR,
  INTEGER_MATRIX,
  LOGICAL_VECTOR,
  LOGICAL_MATRIX,
  FACTOR,
  INVALID
} column_type;

static bool numericVectorIsConstant(SEXP x, column_type t);
static bool integerVectorIsConstant(const int* i, size_t n);

static size_t getNumRowsForDataFrame(SEXP x);

static void getColumnTypes(SEXP x, column_type* columnTypes);
static void tableFactor(SEXP x, int* instanceCount);
static void countMatrixColumns(SEXP x, const column_type* columnTypes, SEXP dropPatternExpr, bool createDropPattern, size_t* result);
static int createMatrix(SEXP x, size_t numRows, SEXP result, const column_type* columnTypes, SEXP dropPatternExpr);
static int setFactorColumnName(SEXP dfNames, size_t dfIndex, SEXP levelNames, size_t levelIndex, SEXP resultNames, size_t resultIndex);


char* concatenateStrings(const char* s1, const char* s2);

SEXP makeModelMatrixFromDataFrame(SEXP x, SEXP dropColumnsExpr)
{
  int errorCode = 0;
  SEXP result = NULL_USER_OBJECT;
  SEXP dropPatternExpr = NULL_USER_OBJECT;
  size_t protectCount = 0;
  
  size_t numInputColumns = (size_t) LENGTH(x);
  size_t numOutputColumns = 0;
  
  column_type columnTypes[numInputColumns];
  
  getColumnTypes(x, columnTypes);
  
  bool createDropPattern = false;
  if (IS_LOGICAL(dropColumnsExpr)) {
    createDropPattern = LOGICAL(dropColumnsExpr)[0] == TRUE;
    if (createDropPattern) {
      dropPatternExpr = PROTECT(NEW_LIST(numInputColumns));
      ++protectCount;
      if (GET_NAMES(x) != NULL_USER_OBJECT) SET_NAMES(dropPatternExpr, GET_NAMES(x));
    }
  }
  if (!createDropPattern && IS_LIST(dropColumnsExpr) && !IS_LOGICAL(dropColumnsExpr))
    dropPatternExpr = dropColumnsExpr;
  
  countMatrixColumns(x, columnTypes, dropPatternExpr, createDropPattern, &numOutputColumns);
  
  size_t numRows = getNumRowsForDataFrame(x);
  
  if (numRows == 0) {
    errorCode = EINVAL;
    goto mkmm_cleanup;
  }
  
  result = PROTECT(NEW_NUMERIC(numRows * numOutputColumns));
  ++protectCount;
  rc_setDims(result, (int) numRows, (int) numOutputColumns, -1);
  SET_DIMNAMES(result, NEW_LIST(2));
  SET_VECTOR_ELT(GET_DIMNAMES(result), 1, NEW_STRING(numOutputColumns));
  
  errorCode = createMatrix(x, numRows, result, columnTypes, dropPatternExpr);
  
mkmm_cleanup:
  if (protectCount > 0) UNPROTECT(protectCount);
  
  if (errorCode != 0) {
    warning("error in makeModelMatrix: %s", strerror(errorCode));
    return NULL_USER_OBJECT;
  }
  
  if (dropPatternExpr != NULL) SET_ATTR(result, install("drop"), dropPatternExpr);
  
  return result;
}

static void getColumnTypes(SEXP x, column_type* columnTypes)
{
  size_t numColumns = LENGTH(x);
  for (size_t i = 0; i < numColumns; ++i) {
    SEXP col = VECTOR_ELT(x, i);
    switch (TYPEOF(col)) {
      case REALSXP:
      {
        SEXP dimsExpr = GET_DIM(col);
        columnTypes[i] = (dimsExpr == NULL_USER_OBJECT ? REAL_VECTOR : REAL_MATRIX);
      }
      break;
      case INTSXP:
      {
        SEXP dimsExpr = GET_DIM(col);
        if (dimsExpr == NULL_USER_OBJECT) {
          SEXP levelsExpr = GET_LEVELS(col);
          columnTypes[i] = (levelsExpr == NULL_USER_OBJECT ? INTEGER_VECTOR : FACTOR);
        } else {
          columnTypes[i] = INTEGER_MATRIX;
        }
      }
      break;
      case LGLSXP:
      {
        SEXP dimsExpr = GET_DIM(col);
        columnTypes[i] = (dimsExpr == NULL_USER_OBJECT ? LOGICAL_VECTOR : LOGICAL_MATRIX);
      }
      break;
      default:
      columnTypes[i] = INVALID;
    }
  }
}


static void tableFactor(SEXP x, int* instanceCounts)
{
  SEXP levelsExpr = GET_LEVELS(x);
  size_t numLevels = (size_t) LENGTH(levelsExpr);
    
  for (size_t i = 0; i < numLevels; ++i) instanceCounts[i] = 0;
    
  int* columnData = INTEGER(x);
  size_t columnLength = (size_t) LENGTH(x);
  for (size_t i = 0; i < columnLength; ++i) ++instanceCounts[columnData[i] - 1];
}

static bool numericVectorIsConstant(SEXP x, column_type t) {
  switch (t) {
    case REAL_VECTOR:
    return ext_vectorIsConstant(REAL(x), (size_t) LENGTH(x));
    case INTEGER_VECTOR:
    case LOGICAL_VECTOR:
    return integerVectorIsConstant(INTEGER(x), (size_t) LENGTH(x));
    default:
    break;
  }
  return false;
}

static bool integerVectorIsConstant(const int* i, size_t n) {
  if (n <= 1) return true;
  
  for (size_t j = 1; j < n; ++j) {
    if (i[j] != i[j - 1]) return false;
  }
  
  return true;
}

void countMatrixColumns(SEXP x, const column_type* columnTypes, SEXP dropPatternExpr, bool createDropPattern, size_t* result)
{
  size_t numColumns = LENGTH(x);
  bool dropColumn;
  for (size_t i = 0; i < numColumns; ++i) {
    SEXP col = VECTOR_ELT(x, i);
    
    switch (columnTypes[i]) {
      case REAL_VECTOR:
      case INTEGER_VECTOR:
      case LOGICAL_VECTOR:
      {
        if (dropPatternExpr != NULL_USER_OBJECT) {
          if (createDropPattern) {
            dropColumn = numericVectorIsConstant(col, columnTypes[i]);
            
            SET_VECTOR_ELT(dropPatternExpr, i, NEW_LOGICAL(1));
            LOGICAL(VECTOR_ELT(dropPatternExpr, i))[0] = dropColumn ? TRUE : FALSE;
          } else {
            dropColumn = LOGICAL(VECTOR_ELT(dropPatternExpr, i))[0];
          }
          if (!dropColumn) *result += 1;
        } else {
          *result += 1;
        }
      }
      break;
      
      case REAL_MATRIX:
      {
        double* colData = REAL(col);
        int* dims = INTEGER(GET_DIM(col));
        size_t numRows = dims[0], numCols = dims[1];
        
        if (dropPatternExpr != NULL_USER_OBJECT) {
          if (createDropPattern) {
            SET_VECTOR_ELT(dropPatternExpr, i, NEW_LOGICAL(numCols));
            int* dropPattern = LOGICAL(VECTOR_ELT(dropPatternExpr, i));
            
            for (size_t j = 0; j < numCols; ++j) {
              dropColumn = ext_vectorIsConstant(colData + j * numRows, numRows);
              dropPattern[j] = dropColumn;
              if (!dropColumn) *result += 1;
            }
          } else {
            int* dropPattern = LOGICAL(VECTOR_ELT(dropPatternExpr, i));
            for (size_t j = 0; j < numCols; ++j)
              if (dropPattern[j] == 0) *result += 1;
          }
        } else {
          *result += numCols;
        }
      }
      break;
      
      case INTEGER_MATRIX:
      case LOGICAL_MATRIX:
      {
        int* colData = INTEGER(col);
        int* dims = INTEGER(GET_DIM(col));
        size_t numRows = dims[0], numCols = dims[1];
        
        if (dropPatternExpr != NULL_USER_OBJECT) {
          if (createDropPattern) {
            SET_VECTOR_ELT(dropPatternExpr, i, NEW_LOGICAL(numCols));
            int* dropPattern = LOGICAL(VECTOR_ELT(dropPatternExpr, i));
          
            for (size_t j = 0; j < numCols; ++j) {
              dropColumn = integerVectorIsConstant(colData + j * numRows, numRows);
              dropPattern[j] = dropColumn;
              if (!dropColumn) *result += 1;
            }
          } else {
            int* dropPattern = LOGICAL(VECTOR_ELT(dropPatternExpr, i));
            for (size_t j = 0; j < numCols; ++j)
              if (dropPattern[j] == 0) *result += 1;
          }
        } else {
          *result += numCols;
        }
      }
      break;
      case FACTOR:
      {
        SEXP levelsExpr = GET_LEVELS(col);
        size_t numLevels = LENGTH(levelsExpr);
        
        if (dropPatternExpr != NULL_USER_OBJECT) {
          int* factorInstanceCounts;
          if (createDropPattern) {
            SET_VECTOR_ELT(dropPatternExpr, i, NEW_INTEGER(numLevels));
            factorInstanceCounts = INTEGER(VECTOR_ELT(dropPatternExpr, i));
            tableFactor(col, factorInstanceCounts);
          } else {
            factorInstanceCounts = INTEGER(VECTOR_ELT(dropPatternExpr, i));
          }
          
          size_t numLevelsPerFactor = 0;
          for (size_t j = 0; j < numLevels; ++j) if (factorInstanceCounts[j] > 0) ++numLevelsPerFactor;
          
          if (numLevelsPerFactor == 2) {
            *result += 1;
          } else if (numLevelsPerFactor > 2) {
            *result += numLevelsPerFactor;
          }
        } else {
          *result += (numLevels <= 2 ? 1 : numLevels);
        }
      }
      default:
      break;
    }
  }
}

static int createMatrix(SEXP x, size_t numRows, SEXP resultExpr, const column_type* columnTypes, SEXP dropPatternExpr)
{
  SEXP names = GET_NAMES(x);
  double* result = REAL(resultExpr);
  SEXP resultNames = VECTOR_ELT(GET_DIMNAMES(resultExpr), 1);
  
  size_t numColumns = LENGTH(x);
  size_t resultCol = 0;
  
  for (size_t i = 0; i < numColumns; ++i) {
    SEXP col = VECTOR_ELT(x, i);
    switch (columnTypes[i]) {
      case REAL_VECTOR:
      if (dropPatternExpr == NULL_USER_OBJECT || LOGICAL(VECTOR_ELT(dropPatternExpr, i))[0] == FALSE) {
        memcpy(result + numRows * resultCol, (const double*) REAL(col), numRows * sizeof(double));
        if (names != NULL_USER_OBJECT) SET_STRING_ELT(resultNames, resultCol, STRING_ELT(names, i));
        ++resultCol; 
      }
      break;
      
      case INTEGER_VECTOR:
      case LOGICAL_VECTOR:
      if (dropPatternExpr == NULL_USER_OBJECT || LOGICAL(VECTOR_ELT(dropPatternExpr, i))[0] == FALSE) {
        int* colData = INTEGER(col);
        for (size_t j = 0; j < numRows; ++j) result[j + numRows * resultCol] = (double) colData[j];
        if (names != NULL_USER_OBJECT) SET_STRING_ELT(resultNames, resultCol, STRING_ELT(names, i));
        ++resultCol;
      }
      break;
      
      
      case REAL_MATRIX:
      {
        size_t numElementCols = INTEGER(GET_DIM(col))[1];
        double* colData = REAL(col);
        SEXP colNames = GET_DIMNAMES(col) == NULL_USER_OBJECT ? NULL_USER_OBJECT : VECTOR_ELT(GET_DIMNAMES(col), 1);
        int* dropPattern = dropPatternExpr == NULL_USER_OBJECT ? NULL : INTEGER(VECTOR_ELT(dropPatternExpr, i));
        
        for (size_t j = 0; j < numElementCols; ++j) {
          if (dropPattern == NULL || dropPattern[j] == FALSE) {
            memcpy(result + numRows * resultCol, colData + numRows * j, numRows * sizeof(double));
            if (names != NULL_USER_OBJECT && colNames != NULL_USER_OBJECT) {
              char* colName = concatenateStrings(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(colNames, j)));
              SET_STRING_ELT(resultNames, resultCol, CREATE_STRING_VECTOR(colName));
              free(colName);
            } else if (names != NULL_USER_OBJECT) {
              char buffer[16];
              snprintf(buffer, 16, "%lu", j + 1);
              char* colName = concatenateStrings(CHAR(STRING_ELT(names, i)), buffer);
              SET_STRING_ELT(resultNames, resultCol, CREATE_STRING_VECTOR(colName));
              free(colName);
            } else if (colNames != NULL_USER_OBJECT) {
              SET_STRING_ELT(resultNames, resultCol, STRING_ELT(colNames, j));
            }
            ++resultCol;
          }
        }
      }
      break;
      
      case INTEGER_MATRIX:
      case LOGICAL_MATRIX:
      {
        size_t numElementCols = INTEGER(GET_DIM(col))[1];
        int* colData = INTEGER(col);
        SEXP colNames = GET_DIMNAMES(col) == NULL_USER_OBJECT ? NULL_USER_OBJECT : VECTOR_ELT(GET_DIMNAMES(col), 1);
        int* dropPattern = dropPatternExpr == NULL_USER_OBJECT ? NULL : INTEGER(VECTOR_ELT(dropPatternExpr, i));
        
        for (size_t j = 0; j < numElementCols; ++j) {
          if (dropPattern == NULL || dropPattern[j] == FALSE) {
            for (size_t k = 0; k < numRows; ++k) result[k + numRows * resultCol] = colData[k + numRows * j];
            if (names != NULL_USER_OBJECT && colNames != NULL_USER_OBJECT) {
              char* colName = concatenateStrings(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(colNames, j)));
              SET_STRING_ELT(resultNames, resultCol, CREATE_STRING_VECTOR(colName));
              free(colName);
            } else if (names != NULL_USER_OBJECT) {
              char buffer[16];
              snprintf(buffer, 16, "%lu", j + 1);
              char* colName = concatenateStrings(CHAR(STRING_ELT(names, i)), buffer);
              SET_STRING_ELT(resultNames, resultCol, CREATE_STRING_VECTOR(colName));
              free(colName);
            } else if (colNames != NULL_USER_OBJECT) {
              SET_STRING_ELT(resultNames, resultCol, STRING_ELT(colNames, j));
            }
            ++resultCol;
          }
        }
      }
      break;
      
      
      case FACTOR:
      {
        SEXP levels = GET_LEVELS(col);
        size_t levelsLength = LENGTH(levels);
        int* colData = INTEGER(col);
        size_t numLevelsPerFactor;
        
        if (dropPatternExpr == NULL_USER_OBJECT) {
          numLevelsPerFactor = levelsLength;
          if (numLevelsPerFactor <= 2) {
            int levelToKeep = numLevelsPerFactor == 2 ? 2 : 1;
            for (size_t j = 0; j < numRows; ++j) result[j + numRows * resultCol] = (colData[j] == levelToKeep ? 1 : 0);
            if (setFactorColumnName(names, i, levels, levelToKeep - 1, resultNames, resultCol) != 0) return ENOMEM;
            ++resultCol;
          } else {
            for (int j = 0; j < (int) levelsLength; ++j) {
              for (size_t k = 0; k < numRows; ++k) result[k + numRows * resultCol] = (colData[k] == (j + 1) ? 1 : 0);
              if (setFactorColumnName(names, i, levels, j, resultNames, resultCol) != 0) return ENOMEM;
              ++resultCol;
            }
          }
        } else {
          numLevelsPerFactor = 0;
          int* factorInstanceCounts = INTEGER(VECTOR_ELT(dropPatternExpr, i));
          for (size_t j = 0; j < levelsLength; ++j) if (factorInstanceCounts[j] > 0) ++numLevelsPerFactor;
          
          if (numLevelsPerFactor == 2) {
            int lastIndex;
            // skip until we find the last level that is actually in the column, make that 1
            for (lastIndex = levelsLength - 1; factorInstanceCounts[lastIndex] == 0 && lastIndex >= 0; --lastIndex) { /* */ }
            // R has factors coded with 1 based indexing
            ++lastIndex;
            for (size_t j = 0; j < numRows; ++j) result[j + numRows * resultCol] = (colData[j] == lastIndex ? 1 : 0);
            if (setFactorColumnName(names, i, levels, lastIndex - 1, resultNames, resultCol) != 0) return ENOMEM;
            ++resultCol;
          } else if (numLevelsPerFactor > 2) {
            for (int j = 0; j < (int) levelsLength; ++j) {
              if (factorInstanceCounts[j] > 0) {
                for (size_t k = 0; k < numRows; ++k) result[k + numRows * resultCol] = (colData[k] == (j + 1) ? 1 : 0);
                if (setFactorColumnName(names, i, levels, j, resultNames, resultCol) != 0) return ENOMEM;
                ++resultCol;
              }
            }
          }
        }
      }
      break;
      
      default:
      break;
    }
  } // close for loop over columns
  
  return 0;
}

static int setFactorColumnName(SEXP dfNames, size_t dfIndex, SEXP levelNames, size_t levelIndex,
                               SEXP resultNames, size_t resultIndex)
{
  if (dfNames != NULL_USER_OBJECT) {
    char* colName = concatenateStrings(CHAR(STRING_ELT(dfNames, dfIndex)), CHAR(STRING_ELT(levelNames, levelIndex)));
    if (colName == NULL) return ENOMEM;
    
    SET_STRING_ELT(resultNames, resultIndex, CREATE_STRING_VECTOR(colName));
    free(colName);
  } else {
    SET_STRING_ELT(resultNames, resultIndex, STRING_ELT(levelNames, levelIndex));
  }
  
  return 0;
}
      

static size_t getNumRowsForDataFrame(SEXP x)
{
  SEXP x_0 = VECTOR_ELT(x, 0);
  SEXP dims = GET_DIM(x_0);
  if (dims == NULL_USER_OBJECT) return (size_t) LENGTH(x_0);
  
  return (size_t) INTEGER(dims)[0];
}

char* concatenateStrings(const char* s1, const char* s2)
{
  size_t l1 = strlen(s1);
  size_t l2 = strlen(s2);
  char* result = malloc(l1 + l2 + 2);
  if (result == NULL) return NULL;
  
  memcpy(result, s1, l1);
  result[l1] = '.';
  memcpy(result + l1 + 1, s2, l2 + 1);
  return result;
}

