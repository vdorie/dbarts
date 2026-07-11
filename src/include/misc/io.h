#ifndef MISC_IO_H
#define MISC_IO_H

#ifdef __cplusplus
extern "C" {
#endif

/// Output hooks. NULL until a host installs them (the R package points them
/// at Rprintf/R_FlushConsole in R_init_dbarts; standalone hosts such as
/// tests/cpp and the benchmark drivers must install their own before any
/// misc code that prints can run). misc itself stays host-agnostic - no
/// stdio, no libR - so no build of misc.a carries an unconditional
/// stdout/stderr symbol. Calling misc output before injection is a host
/// bug; this is not guarded at runtime (fast over safe in C).
extern void (*misc_printf)(const char* format, ...);
extern void (*misc_flushOutput)(void);

#ifdef __cplusplus
}
#endif

#endif // MISC_IO_H
