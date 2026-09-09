#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT([dbarts],[1.0-0],[vdorie@gmail.com])
AC_CONFIG_SRCDIR([src/R_interface.cpp])
AC_CONFIG_AUX_DIR([tools/build-aux])
AC_CONFIG_HEADERS([src/config.hpp])
AC_CONFIG_HEADERS([src/misc/config.h])
AC_CONFIG_HEADERS([src/external/config.h])
AC_CONFIG_HEADERS([src/include/misc/types.h])

m4_include([tools/m4/ax_check_compile_flag.m4])
m4_include([tools/m4/ax_pthread.m4])
m4_include([tools/m4/ax_compiler_vendor.m4])
m4_include([tools/m4/ax_compiler_ext.m4])

# the integer type misc's partition kernels use for cut indices; consumed
# through src/include/misc/types.h by the bartcore engine, which pins it to
# uint16_t (src/bartcore/data.hpp) and asserts the match (src/bartcore/tree.hpp)
AC_DEFINE_UNQUOTED(XINT_TYPE, uint16_t, [Define to the integer type used to represent feature map.])

# The reference build compiles the scalar, fixed-order draw-path kernels; the
# shipped build is free to re-associate them, so its draws are not bitwise
# comparable across vector widths. The recorded equivalence baselines and the
# seed-locked snapshot tests are pinned to the reference build and must run
# there. Default off: the vector build is what ships. The macro has to reach
# both the bridge (src/config.hpp) and the kernels (src/misc/config.h), which
# one AC_DEFINE does because both are AC_CONFIG_HEADERS above. Windows has no
# configure: src/Makevars.win defines the same macro from an environment
# variable of the same name.
AC_ARG_ENABLE([reference-build],
  [AS_HELP_STRING([--enable-reference-build],
    [compile the scalar, fixed-order draw-path kernels, whose draws are bitwise reproducible (default: no)])],
  [], [enable_reference_build=no])
AS_IF([test "x$enable_reference_build" != xno],
  [AC_DEFINE([DBARTS_REFERENCE_BUILD], [1],
     [Define to 1 to compile the scalar, fixed-order draw-path kernels.])])

# Checks for programs.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

# src/Makevars.in sets CXX_STD = CXX20, so the C++ probes below have to run
# against the C++20 toolchain rather than the compiler's default standard
CXX20=`"${R_HOME}/bin/R" CMD config CXX20`
if test -z "${CXX20}"; then
  AC_MSG_ERROR([no C++20 compiler is available; dbarts requires C++20])
fi
CXX20STD=`"${R_HOME}/bin/R" CMD config CXX20STD`
CXX="${CXX20} ${CXX20STD}"
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX20FLAGS`
CXXPICFLAGS=`"${R_HOME}/bin/R" CMD config CXX20PICFLAGS`

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP

AC_LANG(C)
AC_LANG_PUSH(C++)
AC_C_RESTRICT
AC_LANG_POP(C++)

AC_CANONICAL_HOST
AX_COMPILER_VENDOR
AC_SUBST(COMPILER_VENDOR, $ax_cv_c_compiler_vendor)
AX_COMPILER_EXT

# Checks for libraries.
AX_PTHREAD
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([clock_gettime])
AC_CHECK_FUNCS([ffs])
AC_FUNC_ALLOCA

AC_LANG_PUSH(C++)
AC_FUNC_STRERROR_R
AC_LANG_POP

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
