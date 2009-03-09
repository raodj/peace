# This file is processed by automake to generate Makefile.in

#---------------------------------------------------------------------------
#
# Copyright (c) Miami University, Oxford, OH.
# All rights reserved.
#
# Miami University (MU) makes no representations or warranties about
# the suitability of the software, either express or implied,
# including but not limited to the implied warranties of
# merchantability, fitness for a particular purpose, or
# non-infringement.  MU shall not be liable for any damages suffered
# by licensee as a result of using, result of using, modifying or
# distributing this software or its derivatives.
#
# By using or copying this Software, Licensee agrees to abide by the
# intellectual property laws, and all other applicable laws of the
# U.S., and the terms of this license.
#
# Authors: Dhananjai M. Rao       raodm@muohio.edu
#
#---------------------------------------------------------------------------

include $(top_builddir)/Makefile.common.am

AM_CPPFLAGS=$(GLOBALCPPFLAGS)
AM_CXXFLAGS=$(GLOBALCXXFLAGS)

bin_PROGRAMS = pTools
pTools_LDADD = $(STDCPP)

pTools_SOURCES = \
../src/arg_parser.cpp\
../src/arg_parser.h\
../src/EST.cpp\
../src/EST.h\
Common.h\
ShowAlignment.h\
ShowAlignment.cpp\
XFigHelper.h\
XFigHelper.cpp\
DrawMatrix.h\
DrawMatrix.cpp\
ShowMST.h\
ShowMST.cpp\
Tool.h\
Tool.cpp\
MSTguiNode.h\
Align.h\
Align.cpp\
main.cpp

# end of Makefile.am
