# SYNOPSIS
#
#   AX_CXX_SNPRINTF_IN_NAMESPACE_STD
#
# DESCRIPTION
#
#   If the function snprintf is in the cstdio header file and std::
#   namespace, define HAVE_SNPRINTF_IN_NAMESPACE_STD.
#

#serial 7

AC_DEFUN([AX_CXX_SNPRINTF_IN_NAMESPACE_STD],
[AC_CACHE_CHECK([whether snprintf is in std::],
   ax_cv_cxx_snprintf_in_namespace_std,
   [AC_REQUIRE([AX_CXX_NAMESPACE_STD])
    AS_IF([test "X$ax_cv_cxx_have_std_namespace" = Xyes],
      [AC_LANG_PUSH([C++])
       AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
          [[#include <cstdio>]],
          [[std::snprintf;]])
	],
        ax_cv_cxx_snprintf_in_namespace_std=yes,
        ax_cv_cxx_snprintf_in_namespace_std=no)
	AC_LANG_POP([C++])
      ],
      [ax_cv_cxx_snprintf_in_namespace_std=no])
  ])
  AS_IF([test "X$ax_cv_cxx_snprintf_in_namespace_std" = Xyes],
        [AC_DEFINE(HAVE_SNPRINTF_IN_NAMESPACE_STD,,[[define if snprintf is in std::]])])
])
