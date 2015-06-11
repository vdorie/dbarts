#ifndef EXTERNAL_STDDEF_H
#define EXTERNAL_STDDEF_H

#ifdef __cplusplus
#  include <cstddef>
#  define ext_size_t std::size_t
#else
#  include <stddef.h>
#  define ext_size_t size_t
#endif

// the top conditions will fail if restrict is available by default
#if !defined(restrict) && ((defined(__cplusplus) && __cplusplus < 201103L) || (defined(__STDC_VERSION__) && __STDC_VERSION__ < 199901L))
#  if defined(__SUNPRO_C) && defined(__C99FEATURES__)
#  elif defined(__GNUC__) && (__GNUC__ > 2 || __GNUC_MINOR__ >= 92)
#    define restrict __restrict__
#  elif defined(_MSC_VER) && _MSC_VER >= 1500
#    define restrict __restrict
#  else
#    define restrict
#  endif
#endif

#endif // EXTERNAL_STDDEF_H
