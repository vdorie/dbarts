#include "config.h"
#include <misc/io.h>

#include <stdarg.h>
#include <stdio.h>

static void printToStderr(const char* format, ...)
{
  va_list argsPointer;
  va_start(argsPointer, format);
  vfprintf(stderr, format, argsPointer);
  va_end(argsPointer);
}

static void flushStderr(void)
{
  fflush(stderr);
}

void (*misc_printf)(const char* format, ...) = &printToStderr;
void (*misc_flushOutput)(void) = &flushStderr;
