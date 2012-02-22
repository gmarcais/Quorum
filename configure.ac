AC_PREREQ([2.63])
AC_INIT([quorum], [0.1.0], [gmarcais@umd.edu])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([src/error_correct_reads.hpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
AM_SILENT_RULES([yes])
# Uncomment for creating libraries
# AC_PROG_LIBTOOL


# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([rt], [clock_gettime])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Uncomment to get libraries via pkg-config
PKG_CHECK_MODULES([JELLYFISH], [jellyfish-1.1])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_CHECK_FUNCS([])

AC_OUTPUT([Makefile])
