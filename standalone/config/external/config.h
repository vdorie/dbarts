/* external's configuration for a standalone build, which cannot run configure
   (it reads R CMD config): written from src/external/config.h.in with the
   values configure gives on a 64-bit POSIX system. tools/check-win-drift.R
   keeps it in step with the template and DESCRIPTION. */

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
