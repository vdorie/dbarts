#ifndef MISC_MEMALIGN_H
#define MISC_MEMALIGN_H

#ifdef __STRICT_ANSI__
#  define __USE_XOPEN2K 1 // gets posix_memalign when strict ANSI
#endif
#include <stdlib.h>   // malloc, posix_memalign
#undef __USE_XOPEN2K

#if !defined(HAVE_POSIX_MEMALIGN) && defined(HAVE_MALLOC_H)
#  include <malloc.h>   // memalign, __mingw_aligned_malloc
#endif

#include <misc/stddef.h> // size_t

inline int misc_alignedAllocate(void** result, misc_size_t alignment, misc_size_t allocationSize) {
#ifdef HAVE_POSIX_MEMALIGN
  int errorCode = posix_memalign(result, alignment, allocationSize);
  if (errorCode != 0 && *result != NULL) free(*result);
  return errorCode;
#elif defined(__MINGW32__)
  *result = __mingw_aligned_malloc(allocationSize, alignment);
  return (*result == NULL ? ENOMEM : 0);
#else
  *result = memalign(alignment, allocationSize);
  return (*result == NULL ? ENOMEM : 0);
#endif
}

#if defined(HAVE_POSIX_MEMALIGN) || !defined(__MINGW32__)
#  define misc_alignedFree free
#else
#  define misc_alignedFree __mingw_aligned_free
#endif

#endif // MISC_MEMALIGN_H

