#include "config.h"
#include <misc/string.h>

#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <string.h>

#include <misc/adaptiveRadixTree.h>

int misc_str_matchInArray(const char* s, const char* const* strings, size_t numStrings, size_t* matchPos)
{
  misc_art_tree tree;
  misc_art_initialize(&tree);
  
  int errorCode = 0;
  void* searchResult;
  uintptr_t pos = ((uintptr_t) NULL) + 1;
  
  for (size_t i = 0; i < numStrings; ++i) {
    errno = 0;
    if (misc_art_insert(&tree, (const uint8_t*) strings[i], strlen(strings[i]) + 1, (void*) pos++) == NULL &&
        errno != 0)
    {
      errorCode = errno;
      goto MATCH_CLEANUP;
    }
  }
  
  errno = 0;
  searchResult = misc_art_search(&tree, (const uint8_t*) s, strlen(s) + 1);
  if (searchResult == NULL && errno != 0) {
    errorCode = errno;
    goto MATCH_CLEANUP;
  }
  
  *matchPos = (searchResult != NULL) ? (uintptr_t) searchResult - ((uintptr_t) NULL) - 1 : MISC_STR_NO_MATCH;
  
MATCH_CLEANUP:
  misc_art_invalidate(&tree);
  return errorCode;
}

int misc_str_matchAllInArray(const char* const* sa, size_t numA, const char* const* sb, size_t numB, size_t* matchPos)
{
  misc_art_tree tree;
  misc_art_initialize(&tree);
  
  int errorCode = 0;
  void* searchResult;
  uintptr_t pos = ((uintptr_t) NULL) + 1;
  
  for (size_t i = 0; i < numB; ++i) {
    errno = 0;
    if (misc_art_insert(&tree, (const uint8_t*) sb[i], strlen(sb[i]) + 1, (void*) pos++) == NULL &&
        errno != 0)
    {
      errorCode = errno;
      goto MATCH_ALL_CLEANUP;
    }
  }
  
  for (size_t i = 0; i < numA; ++i) {
    errno = 0;
    searchResult = misc_art_search(&tree, (const uint8_t*) sa[i], strlen(sa[i]) + 1);
    if (searchResult == NULL && errno != 0) {
      errorCode = errno;
      goto MATCH_ALL_CLEANUP;
    }
  
    matchPos[i] = (searchResult != NULL) ? (uintptr_t) searchResult - ((uintptr_t) NULL) - 1 : MISC_STR_NO_MATCH;
  }
  
MATCH_ALL_CLEANUP:
  misc_art_invalidate(&tree);
  return errorCode;
}

int misc_str_matchInVArray(const char* s, misc_size_t* matchPos, ...)
{
  misc_art_tree tree;
  misc_art_initialize(&tree);
  
  int errorCode = 0;
  void* searchResult;
  uintptr_t pos = ((uintptr_t) NULL) + 1;
  
  va_list stringsPointer;
  va_start(stringsPointer, matchPos);
  
  const char* string = va_arg(stringsPointer, const char*);
  while (string != NULL) {
    errno = 0;
    if (misc_art_insert(&tree, (const uint8_t*) string, strlen(string) + 1, (void*) pos++) == NULL &&
        errno != 0)
    {
      errorCode = errno;
      break;
    }
    string = va_arg(stringsPointer, const char*);
  }
  va_end(stringsPointer);
  
  if (errorCode != 0) goto VA_MATCH_CLEANUP;
  
  errno = 0;
  searchResult = misc_art_search(&tree, (const uint8_t*) s, strlen(s) + 1);
  if (searchResult == NULL && errno != 0) {
    errorCode = errno;
    goto VA_MATCH_CLEANUP;
  }
  
  *matchPos = (searchResult != NULL) ? (uintptr_t) searchResult - ((uintptr_t) NULL) - 1 : MISC_STR_NO_MATCH;
  
  // misc_art_print(&tree);

VA_MATCH_CLEANUP:
  misc_art_invalidate(&tree);
  return errorCode;
}

