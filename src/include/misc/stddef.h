#ifndef MISC_STDDEF_H
#define MISC_STDDEF_H

#ifdef __cplusplus
#  include <cstddef>
#  include <cstdint>
#  define misc_size_t std::size_t
#else
#  include <stddef.h>
#  include <stdint.h>
#  define misc_size_t size_t
#endif

// Narrowed per-observation GATHER INDEX type, DISTINCT from the length/count
// type misc_size_t (which stays size_t). Halves the hot index buffers
// (docs/design/reduced-precision-storage.md sec 3a); a pure retype, so the
// partition/suffstat kernels stay bitwise-draw-preserving. Lengths and counts
// remain misc_size_t; only index-carrying storage narrows to this.
#ifdef __cplusplus
typedef std::uint32_t misc_index_t;
#else
typedef uint32_t misc_index_t;
#endif

// the top conditions will fail if restrict is available by default
#if !defined(restrict) && (defined(__cplusplus) || (defined(__STDC_VERSION__) && __STDC_VERSION__ < 199901L))
#  if defined(__SUNPRO_C) && defined(__C99FEATURES__)
#  elif defined(__GNUC__) && (__GNUC__ > 2 || __GNUC_MINOR__ >= 92)
#    define restrict __restrict__
#  elif defined(_MSC_VER) && _MSC_VER >= 1500
#    define restrict __restrict
#  else
#    define restrict
#  endif
#endif

#endif // MISC_STDDEF_H

