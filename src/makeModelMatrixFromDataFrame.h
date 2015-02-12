#ifndef MAKE_MODEL_MATRIX_FROM_DATA_FRAME_H
#define MAKE_MODEL_MATRIX_FROM_DATA_FRAME_H

#ifdef __cplusplus
extern "C" {
#endif

#include <Rdefines.h>

SEXP makeModelMatrixFromDataFrame(SEXP x, SEXP dropColumns);

#ifdef __cplusplus
}
#endif

#endif // MAKE_MODEL_MATRIX_FROM_DATA_FRAME_H

