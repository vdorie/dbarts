#ifndef EXTERNAL_R_H
#define EXTERNAL_R_H

// imports R.h while doing the least to pollute namespaces

#include <Rversion.h>

// for older versions of R we attempt to not include unnecessary headers,
// which assists in checking namespace and inclusion correctness
#if R_VERSION <= R_Version(3,3,1)
#  define NO_C_HEADERS
#  ifdef __cplusplus
#    include <cstddef>
using std::size_t;
#  else
#    include <misc/stddef.h>
#  endif
#endif

#include <R.h>

// prevents R_ext/Error.h from mapping Rf_error -> error and Rf_warning -> warning
#define R_NO_REMAP
#include <R.h>

#undef NO_C_HEADERS
#undef R_NO_REMAP

#endif // EXTERNAL_R_H

