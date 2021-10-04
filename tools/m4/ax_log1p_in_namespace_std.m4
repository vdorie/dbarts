# SYNOPSIS
#
#   AX_CXX_LOG1P_IN_NAMESPACE_STD
#
# DESCRIPTION
#
#   If the function log1p is in the cmath header file and std::
#   namespace, define HAVE_LOG1P_IN_NAMESPACE_STD.
#
# LICENSE
#
#   Copyright (c) 2019 Vincent Dorie <vdorie@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_CXX_LOG1P_IN_NAMESPACE_STD],
[AC_CACHE_CHECK([whether log1p is in std::],
   ax_cv_cxx_log1p_in_namespace_std,
   [AC_REQUIRE([AX_CXX_NAMESPACE_STD])
    AS_IF([test "X$ax_cv_cxx_have_std_namespace" = Xyes],
      [AC_LANG_PUSH([C++])
       AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
          [[#include <cmath>]],
          [[std::log1p;]])
	],
        ax_cv_cxx_log1p_in_namespace_std=yes,
        ax_cv_cxx_log1p_in_namespace_std=no)
	AC_LANG_POP([C++])
      ],
      [ax_cv_cxx_log1p_in_namespace_std=no])
  ])
  AS_IF([test "X$ax_cv_cxx_log1p_in_namespace_std" = Xyes],
        [AC_DEFINE(HAVE_LOG1P_IN_NAMESPACE_STD,,[[define if log1p is in std::]])])
])
