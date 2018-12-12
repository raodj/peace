# This file is processed by automake to generate Makefile.in

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

include $(top_builddir)/Makefile.common.am

AM_CPPFLAGS=$(GLOBALCPPFLAGS)
AM_CXXFLAGS=$(GLOBALCXXFLAGS)

AM_CPPFLAGS += -Ibaton -Icommon -Ioutput -Iinput -Iframework -Ishared -Ifiltering -Iclustering -Iassembly

SUBDIRS = framework shared input output filtering clustering assembly baton

noinst_LTLIBRARIES = libpeace.la
bin_PROGRAMS = peace

libpeace_la_SOURCES = \
	PEACE.h \
	PEACE.cpp \
	config.h \
	Version.h

libpeace_la_LIBADD = \
	framework/libframework.la \
	shared/libshared.la \
	output/liboutput.la \
	filtering/libfiltering.la \
	clustering/libclustering.la \
	assembly/libassembly.la \
	baton/libbaton.la \
	input/libinput.la

peace_LDFLAGS = -L.
peace_LDADD = $(STDCPP) libpeace.la
#input/libinput.a output/liboutput.a filtering/libfiltering.a clustering/libclustering.a shared/libshared.a framework/libframework.a output/liboutput.a assembly/libassembly.a shared/libshared.a

peace_SOURCES = \
	main.cpp\
	InteractiveConsole.h \
	InteractiveConsole.cpp

# end of Makefile.am
