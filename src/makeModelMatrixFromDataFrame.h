#ifndef MAKE_MODEL_MATRIX_FROM_DATA_FRAME_H
#define MAKE_MODEL_MATRIX_FROM_DATA_FRAME_H

#include <external/Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP dbarts_makeModelMatrixFromDataFrame(SEXP x, SEXP dropColumns);

#ifdef __cplusplus
}
#endif

#endif // MAKE_MODEL_MATRIX_FROM_DATA_FRAME_H

