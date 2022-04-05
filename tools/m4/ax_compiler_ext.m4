# SYNOPSIS
#
#   AX_COMPILER_EXT
#
# DESCRIPTION
#
#   Checks compiler for support for SIMD extensions. Notably, it does not check
#   if the processor supports the instructions.
#
#   On x86 platforms, this macro calls:
#
#     AC_SUBST(SSE2_FLAGS)
#     AC_SUBST(SSE4_1_FLAGS)
#     AC_SUBST(AVX_FLAGS)
#     AC_SUBST(AVX2_FLAGS)
#
#   And defines:
#
#     COMPILER_SUPPORTS_SSE2 / COMPILER_SUPPORTS_SSE4_1
#     COMPILER_SUPPORTS_AVX / COMPILER_SUPPORTS_AVX2
#
#   On ARM platforms, it calls:
#
#    AC_SUBST(NEON_FLAGS)
#
#   And defines:
#
#     COMPILER_SUPPORTS_NEON
#
# LICENSE
#
#   Copyright (c) 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
#   Copyright (c) 2013,2015 Michael Petch <mpetch@capp-sysware.com>
#   Copyright (c) 2017 Rafael de Lucena Valle <rafaeldelucena@gmail.com>
#   Copyright (c) 2022 Vincent Dorie <vdorie@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_COMPILER_EXT],
[
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_REQUIRE([AC_PROG_CC])
  
  AC_REQUIRE([AX_CHECK_COMPILE_FLAG])
  
  case $host_cpu in
  i[[3456]]86*|x86_64*|amd64*)
    
    SSE2_FLAG=""
    SSE4_1_FLAG=""
    AVX_FLAG=""
    AVX2_FLAG=""
    
    case "$ax_cv_c_compiler_vendor" in
    sun)
      AC_CHECK_HEADER("xmmintrin.h",[
        ax_compiler_supports_sse2_ext=yes
        SSE2_FLAG=-xarch=sse2
      ])
      AC_CHECK_HEADER("nmmintrin.h",[
        ax_compiler_supports_sse4_1_ext=yes
        SSE4_1_FLAG=-xarch=sse4_1
      ])
      AC_CHECK_HEADER("immintrin.h",[
        ax_compiler_supports_avx_ext=yes
        AVX_FLAG=-xarch=avx
      ])
      if test x"$ax_compiler_supports_avx_ext" = x"yes"; then
        AX_CHECK_COMPILE_FLAG(-xarch=avx2,[
          ax_compiler_supports_avx2_ext=yes
          AVX2_FLAG=-xarch=avx2
        ])
      fi
      ;;
    
    *)
      AX_CHECK_COMPILE_FLAG(-msse2,[
        ax_compiler_supports_sse2_ext=yes
        SSE2_FLAG=-msse2
      ])
      AX_CHECK_COMPILE_FLAG(-msse4.1,[
        ax_compiler_supports_sse4_1_ext=yes
        SSE4_1_FLAG=-msse4.1
      ])
      AX_CHECK_COMPILE_FLAG(-mavx,[
        ax_compiler_supports_avx_ext=yes
        AVX_FLAG=-mavx
      ])
      if test x"$ax_cv_c_compiler_vendor" = x"intel"; then
        AX_CHECK_COMPILE_FLAG(-xcore-avx2,[
          ax_compiler_supports_avx2_ext=yes
          AVX2_FLAG=-xcore-avx2
        ])
      else
        AX_CHECK_COMPILE_FLAG(-mavx2,[
          ax_compiler_supports_avx2_ext=yes
          AVX2_FLAG=-mavx2
        ])
      fi
    ;;
    esac

    if test x"$ax_compiler_supports_sse2_ext" = x"yes"; then
      AC_DEFINE(COMPILER_SUPPORTS_SSE2,1,[Support SSE2 (Streaming SIMD Extensions 2) instructions])
    fi
    
    if test x"$ax_compiler_supports_sse41_ext" = x"yes"; then
      AC_DEFINE(COMPILER_SUPPORTS_SSE4_1,1,[Support SSSE4.1 (Streaming SIMD Extensions 4.1) instructions])
    fi
    
    if test x"$ax_compiler_supports_avx_ext" = x"yes"; then
      AC_DEFINE(COMPILER_SUPPORTS_AVX,1,[Support AVX (Advanced Vector Extensions) instructions])
    fi
    
    if test x"$ax_compiler_supports_avx2_ext" = x"yes"; then
      AC_DEFINE(COMPILER_SUPPORTS_AVX2,1,[Support AVX2 (Advanced Vector Extensions 2) instructions])
    fi
    
    AH_TEMPLATE([COMPILER_SUPPORTS_SSE2],[Define to 1 if the compiler can issue Streaming SIMD Extensions instructions.])
    AH_TEMPLATE([COMPILER_SUPPORTS_SSE4_1],[Define to 1 if the compiler can issue Streaming SIMD Extensions 4.1 instructions.])
    AH_TEMPLATE([COMPILER_SUPPORTS_AVX],[Define to 1 if the compilter can issue Advanced Vector Extensions instructions.])
    AH_TEMPLATE([COMPILER_SUPPORTS_AVX2],[Define to 1 if the compilter can issue Advanced Vector Extensions 2 instructions.])
    
    AC_SUBST(SSE2_FLAG)
    AC_SUBST(SSE4_1_FLAG)
    AC_SUBST(AVX_FLAG)
    AC_SUBST(AVX2_FLAG)
  ;;

  aarch64)
    NEON_FLAG=""

    AC_CHECK_HEADER("arm_neon.h",[
      ax_compiler_supports_neon_ext=yes
    ])
    dnl if test x"$ax_compiler_supports_neon_ext" = x"yes"; then
    dnl  AX_CHECK_COMPILE_FLAG(-mfloat-abi=softfp,
    dnl    [ax_compiler_supports_neon_ext=yes],
    dnl    [ax_compiler_supports_neon_ext=no]
    dnl  )
    dnl  if test x"$ax_compiler_supports_neon_ext" = x"yes"; then
    dnl    AX_CHECK_COMPILE_FLAG(-mfpu=neon,
    dnl      [ax_compiler_supports_neon_ext=yes],
    dnl      [ax_compiler_supports_neon_ext=no]
    dnl    )
    dnl  fi
    dnl fi
    if test x"$ax_compiler_supports_neon_ext" = x"yes"; then
      AC_DEFINE(COMPILER_SUPPORTS_NEON,1,[Support NEON instructions])
      NEON_FLAG="-mfloat-abi=softfp -mfpu=neon"
    fi
    
    AH_TEMPLATE([COMPILER_SUPPORTS_NEON],[Define to 1 if the compilter can issue NEON instructions.])
    
    AC_SUBST(NEON_FLAG)
  ;;
  esac
])
