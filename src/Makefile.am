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

noinst_LIBRARIES = libpeace.a
bin_PROGRAMS = peace

peace_LDFLAGS = -L.
peace_LDADD = $(STDCPP) libpeace.a

libpeace_a_SOURCES = \
arg_parser.cpp\
arg_parser.h\
ESTAnalyzer.cpp\
ESTAnalyzerFactory.cpp\
ESTAnalyzerFactory.h\
ESTAnalyzer.h\
EST.cpp\
EST.h\
FMWSCA.cpp\
FMWSCA.h\
D2.h\
D2.cpp\
FWAnalyzer.cpp\
FWAnalyzer.h\
MatrixFileAnalyzer.h \
MatrixFileAnalyzer.cpp \
HTMLHeader.h\
ResultLog.cpp\
ResultLog.h\
Utilities.cpp\
Utilities.h\
ClusterMaker.h\
ClusterMaker.cpp\
ClusterMakerFactory.h\
ClusterMakerFactory.cpp\
MSTNode.h\
MSTNode.cpp\
MST.h\
MST.cpp\
MSTCache.h\
MSTCache.cpp\
MSTClusterMaker.h\
MSTClusterMaker.cpp\
CLU.h\
CLU.cpp\
ESTCustomData.h\
MSTCluster.h\
MSTCluster.cpp\
MPIStats.h\
MPIStats.cpp\
Heuristic.h\
Heuristic.cpp\
HeuristicFactory.h\
HeuristicFactory.cpp\
HeuristicChain.h\
HeuristicChain.cpp\
UVSampleHeuristic.h\
UVSampleHeuristic.cpp\
NewUVHeuristic.h\
NewUVHeuristic.cpp\
ESTCodec.cpp\
ESTCodec.h\
TVHeuristic.h\
TVHeuristic.cpp\
D2Zim.h\
D2Zim.cpp\
TwoPassD2.h\
TwoPassD2.cpp\
MSTMultiListCache.h\
MSTMultiListCache.cpp\
MSTHeapCache.h\
MSTHeapCache.cpp\
CachedESTInfoHelper.h\
TransCacheEntry.h\
TransCacheEntry.cpp\
TransMSTClusterMaker.cpp\
TransMSTClusterMaker.h\
MPIHelper.h \
MPIHelper.cpp \
Filter.h \
Filter.cpp \
LengthFilter.h \
LengthFilter.cpp \
FilterFactory.h \
FilterFactory.cpp \
FilterChain.h \
FilterChain.cpp \
LCFilter.h \
LCFilter.cpp

peace_SOURCES = \
main.cpp\
InteractiveConsole.h\
InteractiveConsole.cpp

# end of Makefile.am
