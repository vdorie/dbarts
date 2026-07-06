#ifndef MISC_IO_H
#define MISC_IO_H

#ifdef __cplusplus
extern "C" {
#endif

/// Output hooks; default to stderr so that misc stays independent of the R
/// runtime. Hosts can repoint them (the package sets Rprintf/R_FlushConsole).
extern void (*misc_printf)(const char* format, ...);
extern void (*misc_flushOutput)(void);

#ifdef __cplusplus
}
#endif

#endif // MISC_IO_H
