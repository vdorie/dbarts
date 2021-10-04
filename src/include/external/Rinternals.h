#ifndef EXTERNAL_RINTERNALS_H
#define EXTERNAL_RINTERNALS_H

// imports Rinternals.h while doing the least to pollute namespaces

#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#define USE_FC_LEN_T
#endif

#if R_VERSION <= R_Version(3, 3, 1)
// Rinternals.h includes R_ext/Memory.h and R_ext/Utils.h which reference size_t
// Rinternals.h also references FILE from stdio.h
#  define NO_C_HEADERS
#  ifdef __cplusplus
#    include <climits>
#    include <cstddef>
#    include <cstdio>
using std::size_t;
using std::FILE;
#  else
#    include <limits.h>
#    include <misc/stddef.h>
#    include <stdio.h>
#  endif
#endif

// prevents R_ext/Error.h from mapping Rf_error -> error and Rf_warning -> warning
#define R_NO_REMAP
#include <Rinternals.h>

#undef NO_C_HEADERS
#undef R_NO_REMAP
#undef USE_FC_LEN_T

#endif // EXTERNAL_RINTERNALS_H

