#include <external/string.h>
#include "config.h"

#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <string.h>

#include <external/adaptiveRadixTree.h>

int ext_str_matchInArray(const char* s, const char* const* strings, size_t numStrings, size_t* matchPos)
{
  ext_art_tree tree;
  ext_art_initialize(&tree);
  
  int errorCode = 0;
  void* searchResult;
  uintptr_t pos = ((uintptr_t) NULL) + 1;
  
  for (size_t i = 0; i < numStrings; ++i) {
    errno = 0;
    ext_art_insert(&tree, (const uint8_t*) strings[i], strlen(strings[i]), (void*) pos++);
    if (errno != 0) { errorCode = errno; goto MATCH_CLEANUP; }
  }
  
  errno = 0;
  searchResult = ext_art_search(&tree, (const uint8_t*) s, strlen(s));
  if (searchResult == NULL && errno != 0) { errorCode = errno; goto MATCH_CLEANUP; }
  
  *matchPos = (searchResult != NULL) ? (uintptr_t) searchResult - ((uintptr_t) NULL) - 1 : EXT_STR_NO_MATCH;
  
MATCH_CLEANUP:
  ext_art_invalidate(&tree);
  return errorCode;
}

int ext_str_matchInVArray(const char* s, ext_size_t* matchPos, ...)
{
  ext_art_tree tree;
  ext_art_initialize(&tree);
  
  int errorCode = 0;
  void* searchResult;
  uintptr_t pos = ((uintptr_t) NULL) + 1;
  
  va_list stringsPointer;
  va_start(stringsPointer, matchPos);
  
  const char* string = va_arg(stringsPointer, const char*);
  while (string != NULL) {
    errno = 0;
    ext_art_insert(&tree, (const uint8_t*) string, strlen(string), (void*) pos++);
    if (errno != 0) { errorCode = errno; break; }
    string = va_arg(stringsPointer, const char*);
  }
  va_end(stringsPointer);
  
  if (errorCode != 0) goto VA_MATCH_CLEANUP;
  
  errno = 0;
  searchResult = ext_art_search(&tree, (const uint8_t*) s, strlen(s));
  if (searchResult == NULL && errno != 0) { errorCode = errno; goto VA_MATCH_CLEANUP; }
  
  *matchPos = (searchResult != NULL) ? (uintptr_t) searchResult - ((uintptr_t) NULL) - 1 : EXT_STR_NO_MATCH;
  
  // ext_art_print(&tree);

VA_MATCH_CLEANUP:
  ext_art_invalidate(&tree);
  return errorCode;
}

