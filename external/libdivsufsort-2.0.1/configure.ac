## configure.ac for libdivsufsort

AC_PREREQ(2.61)

AC_INIT([libdivsufsort], [2.0.1], [yuta.256@gmail.com])
AC_CONFIG_SRCDIR([include/divsufsort.h.cmake])
AC_CONFIG_HEADER([include/config.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.10.1 no-define dist-bzip2])
AM_MAINTAINER_MODE

# LT_CURRENT  = PROJECT_VERSION_MAJOR + PROJECT_VERSION_MINOR + 1
# LT_AGE      = PROJECT_VERSION_MINOR
# LT_REVISION = PROJECT_VERSION_PATCH
AC_SUBST(LT_CURRENT, 3)
AC_SUBST(LT_AGE, 0)
AC_SUBST(LT_REVISION, 1)
AC_SUBST([PROJECT_NAME], [libdivsufsort])
AC_SUBST([PROJECT_DESCRIPTION], "A lightweight suffix sorting library")
AC_SUBST([PROJECT_VERSION_FULL], [2.0.0])
AC_SUBST([PROJECT_URL], "http://libdivsufsort.googlecode.com/")

## Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

## Checks for compiler output filename suffixes.
AC_OBJEXT
AC_EXEEXT

## Check for build configuration.
#AM_DISABLE_STATIC
AM_DISABLE_SHARED
#AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
AC_SUBST([LIBTOOL_DEPS])

case "$target_os" in
  cygwin* | mingw*)
    LDFLAGS="$LDFLAGS -no-undefined"
    ;;
esac

## Checks for libraries.

## Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([inttypes.h memory.h stddef.h stdint.h stdlib.h string.h strings.h sys/types.h io.h fcntl.h])
AS_IF([test "$ac_cv_header_inttypes_h" == "yes"],
        [AC_SUBST([INCFILE], ["#include <inttypes.h>"])],
      [test "$ac_cv_header_stdint_h" == "yes"],
        [AC_SUBST([INCFILE], ["#include <stdint.h>"])],
      [AC_SUBST([INCFILE], [""])])

## Checks for typedefs, structures, and compiler characteristics.
# sauchar_t
SAUCHAR_TYPE=""
AC_CHECK_TYPES([uint8_t])
if test "$ac_cv_type_uint8_t" = "yes"; then
  SAUCHAR_TYPE="uint8_t"
fi
if test -z "$SAUCHAR_TYPE";then
  AC_CHECK_SIZEOF([unsigned char], 1)
  if test "$ac_cv_sizeof_unsigned_char" = "1";then SAUCHAR_TYPE="unsigned char"; fi
fi
if test -z "$SAUCHAR_TYPE";then
  AC_MSG_ERROR([Cannot find unsigned 8-bit integer type])
fi
AC_SUBST([SAUCHAR_TYPE])

# saint_t and saidx_t
SAINT32_TYPE=""
AC_CHECK_TYPES([int32_t])
if test "$ac_cv_type_int32_t" = "yes"; then
  SAINT32_TYPE="int32_t";
  SAINT32_PRId="PRId32";
fi
if test -z "$SAINT32_TYPE";then
  AC_CHECK_SIZEOF([int], 4)
  if test "$ac_cv_sizeof_int" = "4";then
    SAINT32_TYPE="int";
    SAINT32_PRId="\"d\"";
  fi
fi
if test -z "$SAINT32_TYPE";then
  AC_CHECK_SIZEOF([long], 4)
  if test "$ac_cv_sizeof_long" = "4"; then
    SAINT32_TYPE="long";
    SAINT32_PRId="\"ld\"";
  fi
fi
if test -z "$SAINT32_TYPE";then
  AC_CHECK_SIZEOF([__int32], 4)
  if test "$ac_cv_sizeof___int32" = "4"; then
    SAINT32_TYPE="__int32";
    SAINT32_PRId="\"I32d\"";
  fi
fi
if test -z "$SAINT32_TYPE";then
  AC_CHECK_SIZEOF([short], 4)
  if test "$ac_cv_sizeof_short" = "4"; then
    SAINT32_TYPE="short";
    SAINT32_PRId="\"d\"";
  fi
fi
if test -z "$SAINT32_TYPE";then
  AC_MSG_ERROR([Could not find 32-bit integer type])
fi
AC_SUBST([SAINT32_TYPE])
AC_SUBST([SAINT_PRId], "$SAINT32_PRId")

AC_ARG_ENABLE(divsufsort64, AC_HELP_STRING([--enable-divsufsort64], [build libdivsufsort64]))
if test "$enable_divsufsort64" = "yes"; then
  # saint64_t
  SAINT64_TYPE=""
  AC_CHECK_TYPES([int64_t])
  if test "$ac_cv_type_int64_t" = "yes"; then
    SAINT64_TYPE="int64_t";
    SAINT64_PRId="PRId64";
  fi
  if test -z "$SAINT64_TYPE";then
    AC_CHECK_SIZEOF([long long], 8)
    if test "$ac_cv_sizeof_long_long" = "8";then
      SAINT64_TYPE="long long";
      SAINT64_PRId="\"lld\"";
    fi
  fi
  if test -z "$SAINT64_TYPE";then
    AC_CHECK_SIZEOF([long], 8)
    if test "$ac_cv_sizeof_long" = "8";then
      SAINT64_TYPE="long";
      SAINT64_PRId="\"ld\"";
    fi
  fi
  if test -z "$SAINT64_TYPE";then
    AC_CHECK_SIZEOF([int], 8)
    if test "$ac_cv_sizeof_int" = "8";then
      SAINT64_TYPE="int";
      SAINT64_PRId="\"d\"";
    fi
  fi
  if test -z "$SAINT64_TYPE";then
    AC_CHECK_SIZEOF([__int64], 8)
    if test "$ac_cv_sizeof___int64" = "8";then
      SAINT64_TYPE="__int32";
      SAINT64_PRId="\"I64d\"";
    fi
  fi
  if test -z "$SAINT64_TYPE";then
    AC_MSG_ERROR([Could not find 64-bit integer type])
  fi

  AC_CONFIG_FILES([include/divsufsort64.h:include/divsufsort64.h.in])

  AC_SUBST([SAINT64_TYPE])
  AC_SUBST([SAINT64_PRId])
fi

AM_CONDITIONAL([DIVSUFSORT64], test "$enable_divsufsort64" = "yes")

AC_SUBST([SAINDEX_TYPE], "$SAINT32_TYPE")
AC_SUBST([SAINDEX_PRId], "$SAINT32_PRId")
AC_SUBST([W64BIT], [])

AC_SUBST([DIVSUFSORT_EXPORT], [])
AC_SUBST([DIVSUFSORT_IMPORT], [])

AC_SUBST([LFS_OFF_T], [long])
AC_SUBST([LFS_FOPEN], [fopen])
AC_SUBST([LFS_FTELL], [ftell])
AC_SUBST([LFS_FSEEK], [fseek])
AC_SUBST([LFS_PRID], ["\"ld\""])

AC_C_CONST
AC_C_INLINE
AC_DEFINE(INLINE, [inline], [for inline])
AC_DEFINE(PROJECT_VERSION_FULL, [PACKAGE_VERSION], [Define to the version of this package.])

## Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([fopen_s _setmode setmode _fileno])
if test "$ac_cv_func_setmode" = "yes"; then
  if test "$ac_cv_func__setmode" = "no"; then
    AC_DEFINE(_setmode, [setmode], [for _setmode])
    AC_DEFINE(HAVE__SETMODE, 1, [for _setmode])
  fi
fi

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 include/divsufsort.h:include/divsufsort.h.cmake
                 include/lfs.h:include/lfs.h.cmake
                 lib/Makefile
                 examples/Makefile])
AC_OUTPUT
