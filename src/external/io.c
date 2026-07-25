#include <external/io.h>

#include <stdarg.h>      // varargs stuff
#include <misc/stddef.h> // size_t
#include <stdio.h>       // vsnprintf

// Rf_error, Rprintf imported from header

#define MAX_BUFFER_LENGTH 8192

// Formats into a fixed buffer and normalizes the terminating null to a trailing
// newline, so each entry point below only dispatches the final R call. va_start
// / va_end stay in the callers, as C requires.
static void formatMessage(char* buffer, const char* format, va_list argsPointer)
{
  vsnprintf(buffer, MAX_BUFFER_LENGTH, format, argsPointer);

  for (size_t i = 0; i < MAX_BUFFER_LENGTH; ++i) {
    if (buffer[i] == '\0' && i < MAX_BUFFER_LENGTH - 1) {
      buffer[i] = '\n'; buffer[i + 1] = '\0'; break;
    }
  }
}

NORETURN void ext_throwError(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];

  va_list argsPointer;
  va_start(argsPointer, format);
  formatMessage(buffer, format, argsPointer);
  va_end(argsPointer);

  Rf_error("%s", buffer);
}

void ext_printMessage(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];

  va_list argsPointer;
  va_start(argsPointer, format);
  formatMessage(buffer, format, argsPointer);
  va_end(argsPointer);

  Rprintf("%s", buffer);
}

void ext_issueWarning(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];

  va_list argsPointer;
  va_start(argsPointer, format);
  formatMessage(buffer, format, argsPointer);
  va_end(argsPointer);

  Rf_warning("%s", buffer);
}

