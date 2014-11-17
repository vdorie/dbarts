#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(dbarts, 0.8.6, vjd4@nyu.edu)
AC_CONFIG_SRCDIR([src/R_interface.cpp])
AC_CONFIG_AUX_DIR([tools/build-aux])
AC_CONFIG_HEADER([src/config.hpp])
AC_CONFIG_HEADER([src/dbarts/config.hpp])
AC_CONFIG_HEADER([src/external/config.h])

m4_include([tools/m4/m4_ax_pthread.m4])

AC_ARG_ENABLE([match-bayes-tree],
              [AS_HELP_STRING([--enable-match-bayes-tree],
                              [get results identical to BayesTree])],
              [match_bayes_tree=-DMATCH_BAYES_TREE],
              [unset match_bayes_tree])
AC_SUBST(match_bayes_tree)

AC_ARG_ENABLE([thread-safe-unload],
              [AS_HELP_STRING([--enable-thread-safe-unload],
                              [block when adding/deleting fits for the purposes of package reloading])],
              [thread_safe_unload=-DTHREAD_SAFE_UNLOAD],
              [unset thread_safe_unload])
AC_SUBST(thread_safe_unload)

# Checks for programs.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CXXPICFLAGS=`"${R_HOME}/bin/R" CMD config CXXPICFLAGS`

AC_PROG_CC
AC_PROG_CC_C99
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP

AC_LANG(C)
AC_LANG_PUSH(C++)
AC_C_RESTRICT
AC_LANG_POP(C++)

AC_C_BIGENDIAN

# Checks for libraries.
AX_PTHREAD
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LDFLAGS)
AC_SUBST(PTHREAD_LIBS)

# Checks for header files.
AC_LANG_PUSH(C++)
AC_CHECK_HEADERS([cstdint])
AC_LANG_POP(C++)

if test "$ac_cv_header_cstdint" == no
then
  AC_CHECK_HEADER([stdint.h])
fi

AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([malloc.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT64_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([clock_gettime])
AC_FUNC_ALLOCA

AC_LANG_PUSH(C++)
AC_CHECK_FUNCS([snprintf])

## check that snprintf is in namespace std
if test "$ac_cv_func_snprintf" == yes
then
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [[#include <cstdio>]],
      [[std::snprintf;]])
  ], [
      AC_MSG_RESULT([checking for std::snprintf... yes])
      AC_DEFINE_UNQUOTED([HAVE_STD_SNPRINTF], 1, [Define to 1 if `snprintf' is defined in namespace `std'.])
  ], [
      AC_MSG_RESULT([checking for std::snprintf... no])
##      AC_DEFINE_UNQUOTED([HAVE_STD_SNPRINTF], 0, [Define to 1 if `snprintf' is defined in namespace `std'.])
  ])
fi
AC_LANG_POP


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
