#ifndef EXTERNAL_RINTERNALS_H
#define EXTERNAL_RINTERNALS_H

// imports R.h while doing the least to pollute namespaces

#include <Rversion.h>
#ifdef R_VERSION
#  define __R_MAJOR (R_VERSION >> 16)
#  define __R_MINOR ((R_VERSION - (__R_MAJOR << 16)) >> 8)
#  define __R_PATCH ((R_VERSION - (__R_MAJOR << 16)) - (__R_MINOR << 8))
#endif

#if !defined(__R_MAJOR) || (__R_MAJOR < 3 && (!defined(__R_MINOR) || __R_MAJOR != 3 || __R_MINOR < 4))
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
#    include <stddef.h>
#    include <stdio.h>
#  endif
#endif

// prevents R_ext/Error.h from mapping Rf_error -> error and Rf_warning -> warning
#define R_NO_REMAP
#include <Rinternals.h>

#undef NO_C_HEADERS
#undef R_NO_REMAP

#endif // EXTERNAL_RINTERNALS_H

