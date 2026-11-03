// The support library's message, warning and error entry points for a
// standalone build: all three write to standard output, as the R package's
// write to R's console. Compiled in place of src/external/io.c, never with it.

#include <external/io.h>

#include <stdarg.h>  // va_list
#include <stddef.h>  // size_t
#include <stdio.h>   // vsnprintf, fputs, fflush
#include <stdlib.h>  // abort

#define MAX_BUFFER_LENGTH 8192

// As src/external/io.c: formats into a fixed buffer and ends it with a single
// newline.
static void formatMessage(char* buffer, const char* format, va_list argsPointer)
{
  vsnprintf(buffer, MAX_BUFFER_LENGTH, format, argsPointer);

  for (size_t i = 0; i < MAX_BUFFER_LENGTH; ++i) {
    if (buffer[i] == '\0' && i < MAX_BUFFER_LENGTH - 1) {
      buffer[i] = '\n'; buffer[i + 1] = '\0'; break;
    }
  }
}

// Its callers are invariants no input reaches (an unknown generator kind), so
// there is no caller to return to: the message is flushed and the process
// stops.
NORETURN void ext_throwError(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];

  va_list argsPointer;
  va_start(argsPointer, format);
  formatMessage(buffer, format, argsPointer);
  va_end(argsPointer);

  fputs("Error: ", stdout);
  fputs(buffer, stdout);
  fflush(stdout);
  abort();
}

void ext_printMessage(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];

  va_list argsPointer;
  va_start(argsPointer, format);
  formatMessage(buffer, format, argsPointer);
  va_end(argsPointer);

  fputs(buffer, stdout);
}

void ext_issueWarning(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];

  va_list argsPointer;
  va_start(argsPointer, format);
  formatMessage(buffer, format, argsPointer);
  va_end(argsPointer);

  fputs("Warning: ", stdout);
  fputs(buffer, stdout);
}
