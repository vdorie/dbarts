#ifndef EXTERNAL_R_H
#define EXTERNAL_R_H

// imports R.h while doing the least to pollute namespaces

// prevents R_ext/Error.h from mapping Rf_error -> error and Rf_warning -> warning
#define R_NO_REMAP

// R.h includes R_ext/Memory.h, R_ext/Utils.h, and R_ext/RS.h which reference size_t
// RS.h also defines macros that use functions from string.h, but are not actually
// referenced in that file
#define NO_C_HEADERS
#ifdef __cplusplus
#  include <cstddef>
using std::size_t;
#else
#  include <stddef.h>
#endif

#include <R.h>

#undef NO_C_HEADERS
#undef R_NO_REMAP

#endif // EXTERNAL_R_H

