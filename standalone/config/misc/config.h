/* misc's configuration for a standalone build, which cannot run configure
   (it reads R CMD config): written from src/misc/config.h.in with the values
   configure gives on a 64-bit POSIX system, SIMD as configure chooses on
   x86_64 and on arm64. tools/check-win-drift.R keeps it in step with the
   template and DESCRIPTION. */

#define PACKAGE_BUGREPORT "vdorie@gmail.com"
#define PACKAGE_NAME "dbarts"
#define PACKAGE_STRING "dbarts 1.0-0"
#define PACKAGE_TARNAME "dbarts"
#define PACKAGE_VERSION "1.0-0"
#define STDC_HEADERS 1
#define HAVE_UNISTD_H 1
#define HAVE_CLOCK_GETTIME 1
#define HAVE_ALLOCA_H 1
#define HAVE_FFS 1
#if defined(__x86_64__)
#  define COMPILER_SUPPORTS_SSE2 1
#  define COMPILER_SUPPORTS_SSE4_1 1
#  define COMPILER_SUPPORTS_AVX 1
#  define COMPILER_SUPPORTS_AVX2 1
#endif
#if defined(__aarch64__)
#  define COMPILER_SUPPORTS_NEON 1
#endif
#define XINT_TYPE uint16_t
