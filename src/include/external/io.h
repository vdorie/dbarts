#ifndef EXTERNAL_IO_H
#define EXTERNAL_IO_H

#include <R.h>

// Thank goodness for http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1453.htm
#if _MSC_VER >= 1310 
#  define NORETURN _declspec(noreturn)
#elif __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 5)
#  define NORETURN __attribute__ ((noreturn))
#elif __cplusplus >= 201103L
#  define NORETURN [[noreturn]]
#elif __STDC_VERSION__ >= 201112L 
#  define NORETURN _Noreturn
#else
#  define NORETURN
#endif

#ifdef __cplusplus
extern "C" {
#endif
  
void ext_printMessage(const char* format, ...); // printf w/terminal newline
NORETURN void ext_throwError(const char* format, ...);   // printf w/terminal newline and stops program
void ext_issueWarning(const char* format, ...);
#define ext_printf Rprintf
  
#ifdef __cplusplus
}
#endif

#endif
