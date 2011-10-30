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

AM_CPPFLAGS += -I../framework -I../shared -I../output

noinst_LTLIBRARIES = libclustering.la

libclustering_la_SOURCES = \
	FWAnalyzer.h \
	FWAnalyzer.cpp \
	MatrixFileAnalyzer.h \
	MatrixFileAnalyzer.cpp \
	Heuristic.h \
	Heuristic.cpp \
	UVSampleHeuristic.h \
	UVSampleHeuristic.cpp \
	HeuristicChain.h \
	HeuristicChain.cpp \
	NewUVHeuristic.h \
	NewUVHeuristic.cpp \
	TVHeuristic.h \
	TVHeuristic.cpp \
	OldUVHeuristic.h \
	OldUVHeuristic.cpp \
	OldTVHeuristic.h \
	OldTVHeuristic.cpp \
	HeuristicFactory.h \
	HeuristicFactory.cpp \
	D2.h \
	D2.cpp \
	CLU.h \
	CLU.cpp \
	D2Zim.h \
	D2Zim.cpp \
	FMWSCA.h \
	FMWSCA.cpp \
	OldTwoPassD2.h \
	OldTwoPassD2.cpp \
	TwoPassD2.h \
	TwoPassD2.cpp \
	AdaptiveTwoPassD2.h \
	AdaptiveTwoPassD2.cpp \
	ESTAnalyzerFactory.h  \
	ESTAnalyzerFactory.cpp \
	ParameterSet.h \
	ParameterSetManager.h \
	ParameterSetManager.cpp \
	MSTClusterMaker.h \
	MSTClusterMaker.cpp \
	AdaptiveMSTClusterMaker.h  \
	AdaptiveMSTClusterMaker.cpp \
	NonAdaptiveMSTClusterMaker.h \
	NonAdaptiveMSTClusterMaker.cpp \
	TransMSTClusterMaker.h \
	TransMSTClusterMaker.cpp \
	TransCacheEntry.h \
	TransCacheEntry.cpp \
	MSTCache.h \
	MSTCache.cpp \
	MSTHeapCache.h \
	MSTHeapCache.cpp \
	MSTMultiListCache.h \
	MSTMultiListCache.cpp \
	CachedESTInfo.h \
	ClusterMakerFactory.h \
	ClusterMakerFactory.cpp \
	ClusteringSubSystem.h \
	ClusteringSubSystem.cpp \
	NNMSTClusterMaker.h \
	NNMSTClusterMaker.cpp

# end of Makefile.am
