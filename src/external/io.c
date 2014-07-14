#include <external/io.h>

#include <stdarg.h>

#include <R.h>

#define MAX_BUFFER_LENGTH 8192
NORETURN void ext_throwError(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];
  
  va_list argsPointer;
  va_start(argsPointer, format);
  vsnprintf(buffer, MAX_BUFFER_LENGTH, format, argsPointer);
  va_end(argsPointer);
  
  
  for (size_t i = 0; i < MAX_BUFFER_LENGTH; ++i) {
    if (buffer[i] == '\0' && i < MAX_BUFFER_LENGTH - 1) {
      buffer[i] = '\n'; buffer[i + 1] = '\0'; break;
    }
  }
  
  error(buffer);
}

void ext_printMessage(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];
  
  va_list argsPointer;
  va_start(argsPointer, format);
  vsnprintf(buffer, MAX_BUFFER_LENGTH, format, argsPointer);
  va_end(argsPointer);
  
  for (size_t i = 0; i < MAX_BUFFER_LENGTH; ++i) {
    if (buffer[i] == '\0' && i < MAX_BUFFER_LENGTH - 1) {
      buffer[i] = '\n'; buffer[i + 1] = '\0'; break;
    }
  }
  
  Rprintf(buffer);
}

void ext_issueWarning(const char* format, ...)
{
  char buffer[MAX_BUFFER_LENGTH];
  
  va_list argsPointer;
  va_start(argsPointer, format);
  vsnprintf(buffer, MAX_BUFFER_LENGTH, format, argsPointer);
  va_end(argsPointer);
  
  for (size_t i = 0; i < MAX_BUFFER_LENGTH; ++i) {
    if (buffer[i] == '\0' && i < MAX_BUFFER_LENGTH - 1) {
      buffer[i] = '\n'; buffer[i + 1] = '\0'; break;
    }
  }
  
  warning(buffer);
}
