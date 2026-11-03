/* The bridge's configuration for a standalone build, which cannot run
   configure (it reads R CMD config): written from src/config.hpp.in with the
   values configure gives on a 64-bit POSIX system with a C++20 compiler, SIMD
   as configure chooses on x86_64 and on arm64. tools/check-win-drift.R keeps
   it in step with the template and DESCRIPTION. */

#if defined(__x86_64__)
#  define COMPILER_SUPPORTS_SSE2 1
#  define COMPILER_SUPPORTS_SSE4_1 1
#  define COMPILER_SUPPORTS_AVX 1
#  define COMPILER_SUPPORTS_AVX2 1
#endif
#if defined(__aarch64__)
#  define COMPILER_SUPPORTS_NEON 1
#endif
#define HAVE_ALLOCA 1
#define HAVE_ALLOCA_H 1
#define HAVE_CLOCK_GETTIME 1
#define HAVE_DECL_STRERROR_R 1
#define HAVE_FFS 1
#define HAVE_GETTIMEOFDAY 1
#define HAVE_INTTYPES_H 1
#define HAVE_PTHREAD 1
#define HAVE_PTHREAD_PRIO_INHERIT 1
#define HAVE_STDINT_H 1
#define HAVE_STDIO_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STRERROR_R 1
#define HAVE_STRINGS_H 1
#define HAVE_STRING_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_UNISTD_H 1
#define PACKAGE_BUGREPORT "vdorie@gmail.com"
#define PACKAGE_NAME "dbarts"
#define PACKAGE_STRING "dbarts 1.0-0"
#define PACKAGE_TARNAME "dbarts"
#define PACKAGE_URL ""
#define PACKAGE_VERSION "1.0-0"
#define STDC_HEADERS 1
/* g++ and clang++ on Linux define _GNU_SOURCE, which makes glibc's
   strerror_r the GNU one */
#if defined(__linux__) && defined(_GNU_SOURCE)
#  define STRERROR_R_CHAR_P 1
#endif
#define XINT_TYPE uint16_t
#define restrict /**/
