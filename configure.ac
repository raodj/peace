#   -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#--------------------------------------------------------------------
#
# This file is part of PEACE.
# 
# PEACE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PEACE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
# 
# Miami University makes no representations or warranties about the
# suitability of the software, either express or implied, including
# but not limited to the implied warranties of merchantability,
# fitness for a particular purpose, or non-infringement.  Miami
# University shall not be liable for any damages suffered by licensee
# as a result of using, result of using, modifying or distributing
# this software or its derivatives.
#
# By using or copying this Software, Licensee agrees to abide by the
# intellectual property laws, and all other applicable laws of the
# U.S., and the terms of GNU General Public License (version 3).
#
# Authors:   Dhananjai M. Rao              raodm@muohio.edu
#
#---------------------------------------------------------------------

AC_PREREQ(2.6)
AC_INIT([peace], [0.1.0], [raodm@muohio.edu])
AM_INIT_AUTOMAKE([peace],[0.1.0])
AC_CONFIG_SRCDIR([src/config.h.in])
AC_CONFIG_HEADER([src/config.h])
AC_REQUIRE_CPP

# The macros above automatically fill in some CXXFLAGS so we'll undo it.
CXX=g++
CFLAGS=""
CXXFLAGS=""
CPPFLAGS=""
	
# Custom Definitions

# Define MPI library (http://www.open-mpi.org)
AC_DEFUN([AC_PROG_MPI],
[
    AC_ARG_WITH(mpi,
            AC_HELP_STRING([--with-mpi=PATH],[path to mpicc]),
    [ if test "$withval" != "yes"; then
          AC_MSG_CHECKING([for MPI])
	  if test -x "$withval" ; then
	      MPI="$withval"
	      AC_MSG_RESULT([$MPI])
	  else
	      MPI="no"
	      AC_MSG_RESULT(no)
	  fi
      fi
     ])
     
     if test -z "$MPI"; then
         AC_PATH_PROG([MPI], [mpicxx], [notFound])
     fi

     if test $MPI = notFound; then
        AC_MSG_WARN([MPI could no be found on your system.
		      You can download MPI for free from
		      http://www.open-mpi.org
                      Please check that MPI is installed on your system
		      and your path includes the top-level install directory
		      (example: /usr/include/openmpi), or use
		      --with-mpi=PATH option when running configure.])
     elif test $MPI = no; then
  	AC_MSG_WARN([You have chosen to compile PEACE without MPI.
		     This is OK and PEACE will run in single-process mode
		     (MPI is needed to run in parallel using multiple
		     processes to fully utilize computational power of
		     supercomputing clusters). You can download MPI for
		     free from http://www.open-mpi.org])
     else
         CXX=mpicxx
         AC_SUBST(MPI)
	 AC_DEFINE([HAVE_LIBMPI], [1], [Using MPI to compile PEACE])
     fi
])

# Check if mpicxx requires -lstdc++ option (needed on redhawk)
AC_DEFUN([AC_LINK_STDCPP],
[
     AC_MSG_CHECKING([if $CXX needs stdc++])
     AC_LANG([C++])
     AC_LANG_CONFTEST()
     
     AC_LINK_IFELSE(
         [AC_LANG_PROGRAM([[#include "mpi.h"]],[[std::cout << "Hello";]])],
         [AC_MSG_RESULT([no])
	     STDCPP=""
	 ],
	 [AC_MSG_RESULT([yes])
	     STDCPP="-lstdc++"
	     LIBS="$LIBS -lstdc++"
	 ])

     AC_SUBST(STDCPP)
])

# Check if we are using icc to supress a couple of useless remarks
AC_DEFUN([AC_PROG_ICC],
[
    AC_ARG_WITH(icc,
            AC_HELP_STRING([--with-icc],[Code will be compiled with icc]),
    [ if test "$withval" != "yes"; then
	  icc="yes"
      else
	  icc="no"
      fi
     ])

     if test -z "$icc"; then
         AC_MSG_CHECKING([for icc])
         icc=`mpicxx --version 2>&1 | grep "^icc"`
	 if test -n "$icc"; then
	     icc="yes"
	 else
	     icc="no"
	 fi
         AC_MSG_RESULT([$icc])
     fi

     if test "yes" == "$icc"; then
	  USING_ICC="-DICC"
          ICC_NO_WARN_SOME="-wd981 -wd383"
          ICC_NO_WARN_FOR_CONFIGURE="-wd981 -wd383 -wd1419"
     else
     	  USING_ICC=""
          ICC_NO_WARN_SOME=""
          ICC_NO_WARN_FOR_CONFIGURE=""
     fi

     AC_SUBST(USING_ICC) 
     AC_SUBST(ICC_NO_WARN_SOME)
])

AC_DEFUN([AC_CHECK_ARGS],
[
	AC_MSG_CHECKING([whether _AC_LANG compiler accepts $1])
	AC_LANG_CONFTEST(
			AC_LANG_PROGRAM([[#include <iostream>]],
							[[std::cout << "ac_check_args"]]))
    src=`ls conftest.c*`
	output=`$CXX $1 $CXXFLAGS $LIBS $ICC_NO_WARN_FOR_CONFIGURE $src 2>&1`
	if test -z "$output"; then
	    AC_MSG_RESULT([yes])
	    m4_default([$2], :)
	else
		AC_MSG_RESULT([no])
		m4_default([$3], :)
	fi
])

# Check if we are using gcc to supress a couple of useless remarks
AC_DEFUN([AC_PROG_GCC],
[
    AC_ARG_WITH(gcc,
            AC_HELP_STRING([--with-gcc],[Code will be compiled with gcc]),
    [ if test "$withval" != "yes"; then
	  gcc="yes"
      else
	  gcc="no"
      fi
     ])

     if test -z "$gcc"; then
         AC_MSG_CHECKING([for gcc])
         gcc=`mpicxx --version 2>&1 | grep "GCC"`
	 if test -n "$gcc"; then
	     gcc="yes"
	 else
	     gcc="no"
	 fi
         AC_MSG_RESULT([$gcc])
     fi

     if test "yes" == "$gcc"; then
	 	USING_GCC="-DGCC"
     else
     	USING_GCC=""
     fi

     AC_SUBST(USING_GCC) 
])

# Checks for programs.
AC_PROG_MPI
AC_MSG_CHECKING([properties of (MPI) compiler])
AC_MSG_RESULT([$CXX])
AC_PROG_CXX
AC_LINK_STDCPP
AC_PROG_ICC
AC_PROG_GCC
AC_PROG_RANLIB

# Check and set some default compiler flags for the build
AC_CHECK_ARGS([-Minform,warn], [INFORM_WARN="-Minform=warn"], [INFORM_WARN=""])
AC_SUBST(INFORM_WARN)
AC_CHECK_ARGS([-Wall], [WALL="-Wall"], [WALL=""])
AC_SUBST(WALL)
AC_CHECK_ARGS([-Wextra], [WEXTRA="-Wextra"], [WEXTRA=""])
AC_SUBST(WEXTRA)
AC_CHECK_ARGS([-pipe], [PIPE="-pipe"], [PIPE=""])
AC_SUBST(PIPE)


# Checks for libraries.
#AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([m],[asin],,,[-lstdc++])
AC_CHECK_FUNCS([MPI_Init])
AC_CHECK_LIB([readline],[readline],,,[-lstdc++])
#AC_SEARCH_LIBS([gethostname], [c,nsl])
#AC_CHECK_LIB([nsl],[gethostname])
#AC_SEARCH_LIBS([bind],[c,nsl,socket])

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_C_VOLATILE

# checks for the version of hash map is available in c++ land
AC_LANG_PUSH([C])
AC_LANG([C++])
AC_CHECK_HEADERS([tr1/unordered_map ext/hash_map hash_map])
AC_LANG_POP([C++])

# Checks for library functions.
# AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_TYPE_SIGNAL
AC_FUNC_STRTOD
AC_CHECK_FUNCS([bzero floor gettimeofday memset mkdir modf pow rint sqrt strcasecmp strchr strdup strstr strtol])

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 tools/Makefile
		 examples/Makefile])

AC_OUTPUT

# End of script
