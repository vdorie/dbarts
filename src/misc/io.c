#include "config.h"
#include <misc/io.h>

#include <stddef.h>

// No default implementation: misc stays host-agnostic (no stdio, no libR)
// so that no build of misc.a carries a stdout/stderr symbol. Every host
// installs its own hooks before any misc code that prints can run; see
// io.h.

void (*misc_printf)(const char* format, ...) = NULL;
void (*misc_flushOutput)(void) = NULL;
