#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT(dbarts, 0.9.20, vdorie@gmail.com)
AC_CONFIG_SRCDIR([src/R_interface.cpp])
AC_CONFIG_AUX_DIR([tools/build-aux])
AC_CONFIG_HEADER([src/config.hpp])
AC_CONFIG_HEADER([src/dbarts/config.hpp])
AC_CONFIG_HEADER([src/misc/config.h])
AC_CONFIG_HEADER([src/external/config.h])
AC_CONFIG_HEADER([inst/include/dbarts/types.hpp])
AC_CONFIG_HEADER([src/include/misc/types.h])

m4_include([tools/m4/ax_check_compile_flag.m4])
m4_include([tools/m4/ax_pthread.m4])
m4_include([tools/m4/ax_compiler_vendor.m4])
m4_include([tools/m4/ax_compiler_ext.m4])
m4_include([tools/m4/ax_cxx_namespace_std.m4])
m4_include([tools/m4/ax_snprintf_in_namespace_std.m4])
m4_include([tools/m4/ax_log1p_in_namespace_std.m4])
m4_include([tools/m4/ax_func_posix_memalign.m4])

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

AC_ARG_WITH([xint-size],
            [AS_HELP_STRING([--with-xint-size=16],
                            [constrains the number of cutpoints to 2^x-size - 1])],
             ,
             [with_xint_size=no])

case "x$with_xint_size" in
  x8)
    xint_size=uint8_t
  ;;
  x32)
    xint_size=uint32_t
  ;;
  x64)
    xint_size=uint64_t
  ;;
  xyes|xno|x16)
    xint_size=uint16_t
  ;;
esac

AC_DEFINE_UNQUOTED(XINT_TYPE, $xint_size, [Define to the integer type used to represent feature map.])

# Checks for programs.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
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

AC_CANONICAL_HOST
AX_COMPILER_VENDOR
AC_SUBST(COMPILER_VENDOR, $ax_cv_c_compiler_vendor)
AX_COMPILER_EXT

# on SOLARIS, need to strip out some hardware flags from the library,
# so we test for those that are explicitly added with intrinsics
# but not supported by the CPU
mapfile_flag=
case $host_os in
  solaris*)
    case $host_cpu in
      i[[3456]]86*|x86_64*|amd64*)
        if test "$ax_cv_c_compiler_vendor" = sun; then
          m4_include([tools/m4/ax_ext_solaris.m4])
          m4_include([tools/m4/ax_gcc_x86_cpuid.m4])
          m4_include([tools/m4/ax_gcc_x86_avx_xgetbv.m4])
          AX_EXT_SOLARIS
          
          unsupported_isets=
          # only test for those that we explicitly add beyond what
          # the compiler will come up with on "-xarch=native"
          for iset in "avx;AVX" "sse4_1;SSE4_1" "avx2;AVX2"; do
            iset_varname=$(echo $iset | cut -d ";" -f 1)
            iset_mapname=$(echo $iset | cut -d ";" -f 2)
            # have variables test for CPU and OS support
            if test x"$(eval echo \$ax_cv_have_${iset_varname}_ext)" = x"no"; then
              unsupported_isets="${unsupported_isets} ${iset_mapname}"
            fi
          done
          if test ! -z "${unsupported_isets}"; then
            unsupported_isets=${unsupported_isets##[[:space:]]} # remove leading space
            echo "creating mapfile to remove ${unsupported_isets} instruction set(s)"
            echo "\$mapfile_version 2\nCAPABILITY {\n  HW -= ${unsupported_isets};\n};\n" > src/mapfile.solaris
            mapfile_flag="-M mapfile.solaris"
          fi
        fi
      ;;
    esac
  ;;
esac

AC_SUBST(MAPFILE_FLAG, $mapfile_flag)

# Checks for libraries.
AX_PTHREAD
AC_SUBST(PTHREAD_CFLAGS)
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
AC_CHECK_SIZEOF([size_t])
AC_CHECK_ALIGNOF([void*])

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([clock_gettime])
AC_CHECK_FUNCS([ffs])
AC_FUNC_ALLOCA

AX_FUNC_POSIX_MEMALIGN

AC_LANG_PUSH(C++)
AC_CHECK_FUNCS([snprintf])
AX_CXX_SNPRINTF_IN_NAMESPACE_STD
AC_CHECK_FUNCS([log1p])
AX_CXX_LOG1P_IN_NAMESPACE_STD
AC_LANG_POP


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
