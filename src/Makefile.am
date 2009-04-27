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

bin_PROGRAMS = peace

peace_LDFLAGS = -lreadline -lcurses
peace_LDADD = $(STDCPP)

peace_SOURCES = \
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
main.cpp\
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
InteractiveConsole.h\
InteractiveConsole.cpp\
ESTCodec.cpp\
ESTCodec.h

# end of Makefile.am
