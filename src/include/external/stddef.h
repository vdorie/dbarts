#ifndef EXTERNAL_STDDEF_H
#define EXTERNAL_STDDEF_H

#ifdef __cplusplus
#include <cstddef>
#define ext_size_t std::size_t
#else
#include <stddef.h>
#define ext_size_t size_t
#endif

#ifndef restrict
# ifdef __restrict
#  define restrict __restrict
# elif defined(__restrict__)
#  define restrict __restrict__
# else
#  define restrict
# endif
#endif

#endif // EXTERNAL_STDDEF_H
