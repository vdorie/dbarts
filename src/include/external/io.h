#ifndef EXTERNAL_IO_H
#define EXTERNAL_IO_H

#include <R.h>

#ifdef __cplusplus
extern "C" {
#endif
  
  void ext_printMessage(const char* format, ...); // printf w/terminal newline
  void ext_throwError(const char* format, ...);   // printf w/terminal newline and stops program
#define ext_printf Rprintf
  
#ifdef __cplusplus
}
#endif

#endif
