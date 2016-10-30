#ifndef EXTERNAL_R_H
#define EXTERNAL_R_H

// imports R.h while doing the least to pollute namespaces

#include <Rversion.h>
#ifdef R_VERSION
#  define __R_MAJOR (R_VERSION >> 16)
#  define __R_MINOR ((R_VERSION - (__R_MAJOR << 16)) >> 8)
#  define __R_PATCH ((R_VERSION - (__R_MAJOR << 16)) - (__R_MINOR << 8))
#endif

// for older versions of R we attempt to no include unnecessary headers,
// which assists in checking namespace and inclusion correctness
#if !defined(__R_MAJOR) || (__R_MAJOR < 3 && (!defined(__R_MINOR) || __R_MAJOR != 3 || __R_MINOR < 4))
// unfortunately, some headers are necessary to include
#  define NO_C_HEADERS
#  ifdef __cplusplus
#    include <cstddef>
using std::size_t;
#  else
#    include <stddef.h>
#  endif
#endif

#include <R.h>

// prevents R_ext/Error.h from mapping Rf_error -> error and Rf_warning -> warning
#define R_NO_REMAP
#include <R.h>

#undef NO_C_HEADERS
#undef R_NO_REMAP

#endif // EXTERNAL_R_H

