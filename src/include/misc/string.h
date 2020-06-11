#ifndef MISC_STRING_H
#define MISC_STRING_H

#include <misc/stddef.h>

#ifdef __cplusplus
extern "C" {
#  define MISC_STR_NO_MATCH static_cast<misc_size_t>(-1)
#else

#  define MISC_STR_NO_MATCH ((misc_size_t) -1)
#endif

// returns 0 on success, ENOMEM on no mem, EINVAL if something terrible happens
int misc_str_matchInArray(const char* s, const char* const* strings, misc_size_t numStrings, misc_size_t* matchPos);
// returns same as above; for all strings in sa, store the index of their appearance in sb in matchPos or no match 
int misc_str_matchAllInArray(const char* const* sa, size_t numA, const char* const* sb, size_t numB, size_t* matchPos);
// NULL terminate var arg string list
int misc_str_matchInVArray(const char* s, misc_size_t* matchPos, ...);

#ifdef __cplusplus
}
#endif

#endif // MISC_STRING_H

