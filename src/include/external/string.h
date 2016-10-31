#ifndef EXT_STRING_H
#define EXT_STRING_H

#include "stddef.h"

#ifdef __cplusplus
extern "C" {
#  define EXT_STR_NO_MATCH static_cast<ext_size_t>(-1)
#else

#  define EXT_STR_NO_MATCH ((ext_size_t) -1)
#endif

// returns 0 on success, ENOMEM on no mem, EINVAL if something terrible happens
int ext_str_matchInArray(const char* s, const char* const* strings, ext_size_t numStrings, ext_size_t* matchPos);
// NULL terminate var arg string list
int ext_str_matchInVArray(const char* s, ext_size_t* matchPos, ...);

#ifdef __cplusplus
}
#endif

#endif // EXT_STRING_H

