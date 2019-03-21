#ifndef MISC_INTRINSIC_H
#define MISC_INTRINSIC_H

#ifdef _MSC_VER
#  include <intrin.h>
#else 
#  ifdef HAVE_SSE2
#    ifdef __SUNPRO_C
#      include <xmmintrin.h>
#    else
#      include <emmintrin.h> // SSE2 intrinsics
#    endif
#  endif
#  ifdef HAVE_SSE4_1
#    include <smmintrin.h>
#  endif
#  ifdef HAVE_AVX
#    include <immintrin.h>
#  endif
#endif

#ifdef __GNUC__
#  define countTrailingZeros(_X_) __builtin_ctz(_X_)
#  define countLeadingZeros(_X_)  __builtin_clz(_X_)
#elif defined(_MSC_VER)
// it would be possible to use __lzcnt when ABM is available, but that is relatively recent
#  include <intrin.h>
static unsigned long __inline countTrailingZeros(unsigned long x) {
  unsigned long result;
  _BitScanForward(&result, x);
  return result;
}
static unsigned long __inline countLeadingZeros(unsigned long x) {
  unsigned long result;
  _BitScanReverse(&result, x);
  return result;
}
#elif defined(__INTEL_COMPILER)
#  include <immintrin.h>
#  define countTrailingZeros(_X_) _bit_scan_forward(_X_)
#  define countLeadingZeros(_X_) _bit_scan_reverse(_X_)
#elif defined(HAVE_FFS)
#  include <strings.h>
#  define countTrailingZeros(_X_) (ffs(_X_) - 1)
#  define CLZ_MISSING
#elif defined(__MINGW32__)
#  define countTrailingZeros(_X_) (__builtin_ffs(_X_) - 1)
#  define CLZ_MISSING
#else
static inline uint32_t countTrailingZeros(uint32_t v) {
  static const int ctzTable[32] = {
    0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
  };
  
  return ctzTable[((uint32_t) ((v & -v) * 0x077CB531U) >> 27];
}
#define CLZ_MISSING
#endif

#ifdef CLZ_MISSING
static inline unsigned int countLeadingZeros(uint32_t v)
{
  static const int clzTable[32] = {
    0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
    8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
  };
  
  v |= v >> 1; // first round down to one less than a power of 2
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  
  return clzTable[(uint32_t) (v * 0x07C4ACDDU) >> 27];
}
#endif
#undef CLZ_MISSING

#endif // MISC_INTRINSIC_H

