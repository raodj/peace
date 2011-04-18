#ifndef LC_FILTER_CPP
#define LC_FILTER_CPP

//--------------------------------------------------------------------
//
// This file is part of PEACE.
// 
// PEACE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// PEACE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
// 
// Miami University makes no representations or warranties about the
// suitability of the software, either express or implied, including
// but not limited to the implied warranties of merchantability,
// fitness for a particular purpose, or non-infringement.  Miami
// University shall not be liable for any damages suffered by licensee
// as a result of using, result of using, modifying or distributing
// this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include "LCFilter.h"
#include "ESTList.h"
#include "ESTAnalyzer.h"
#include "RuntimeContext.h"
#include "ESTListListener.h"

#include <sstream>

LCFilter::LCFilter(RuntimeContext *runtimeContext) :
    Filter("lcFilter", runtimeContext) {
    // Setup our default patterns and threshold
    threshold = -1;
}

void
LCFilter::addCommandLineArguments(ArgParser& argParser) {
    // The set of arguments for this class.
    const ArgParser::ArgRecord LCFilterArgs[] = {
        {"--lcPatterns", "List of (space separated) patterns to generate dummy ESTs",
         &patternList, ArgParser::STRING_LIST},
        {"--lcThreshold", "Threshold value to detect low complexity sequences",
         &threshold, ArgParser::DOUBLE},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add arguments for the LC filter to this argument parser
    argParser.addValidArguments(LCFilterArgs);
}

int
LCFilter::initialize() {
    // Setup default patterns if user has not specified any
    if (patternList.empty()) {
        patternList.push_back("A");
        patternList.push_back("C");
    }
    ASSERT ( runtimeContext != NULL );
    ESTAnalyzer *analyzer = runtimeContext->getAnalyzer();
    if (analyzer == NULL) {
        std::cerr << "Error: The LCFilter requires a valid EST analyzer to "
                  << "compare and detect\n"
                  << "cDNA fragments with low-complexity regions.\n"
                  << "Ensure a valid analyzer is set via the --analyzer "
                  << "command line argument.\n";
        return 1;
    }
    ASSERT ( analyzer != NULL );    
    // Ensure we have a valid ESTList setup in the runtime context.
    if (runtimeContext->getESTList() == NULL) {
        std::cerr << "Error: A valid ESTList has not been created/set in the "
                  << "runtime context.\nThe LCFilter requires a valid list "
                  << "of cDNA fragments to process.\nThis is most likely an "
                  << "programming error.\n";
        return 2;
    }
    // Check and update our threshold value if one is not specified.
    if (threshold == -1) {
        // Use the invalid metric from analyzer instead.
        threshold = analyzer->getInvalidMetric();
    }
    // Figure out the length of the dummy ESTs
    const int dummyLen = analyzer->getPreferredDummyESTLength();
    // Process the patternList and add a dummy cluster for each
    // pattern.  Process one pattern at a time.
    ASSERT ( !patternList.empty() );
    for(size_t patIdx = 0; (patIdx < patternList.size()); patIdx++) {
        const std::string pattern = patternList[patIdx];
        // Create dummy entry for this pattern.
        addDummyEntry("DummyEST For Pattern " + pattern, pattern, dummyLen);
    }
    // Everything went well
    return 0;
}

void
LCFilter::addDummyEntry(const std::string& fastaID, const std::string& seq,
                        const int length) {
    // First replicate the sequence to obtain full sequence of given
    // length.
    std::string fullSeq = seq;
    // Ensure that the length of the dummy does not exceed the longest
    // fragement we have (as the cluster maker has already been
    // initialized and it has used the longest fragement to create its
    // internal buffers)
    ESTList* estList  = runtimeContext->getESTList();
    ASSERT( estList != NULL );
    const int repeats = std::min(length / seq.length(),
                                 estList->getMaxESTLen() / seq.length());
    for(int i = 0; (i < repeats); i++) {
        fullSeq += seq;
    }
    // Add dummy EST and cluster and save information.
    int estIdx = estList->size();
    EST *est = estList->add(estIdx, fastaID, fullSeq);
    ASSERT ( est != NULL );
    // Flag EST as having already been processed.
    est->setProcessed(true);
    // Let the parameter set manager know about this newly added est.
    if (runtimeContext->getESTListListener() != NULL) {
        ESTListListener *listener = runtimeContext->getESTListListener();
        listener->entriesAdded(estIdx, estIdx + 1);
    }
    // ParameterSetManager::getParameterSetManager()->sequenceAppended();
    // Add a dummy cluster for this filter with a suitable name.
    std::ostringstream clsName;
    clsName << "Low Complexity ESTs (filtered by LCFilter Pattern "
            << seq << "...)";
    int clusterID = addDummyCluster(clsName.str());
    // Save information about the dummy entry.
    dummyESTList.push_back(DummyESTInfo(est, clusterID));
}

void
LCFilter::finalize() {
    // First remove all the dummy ESTs we may have added. Note that
    // this one does not really guarantee that we remove the
    // appropiate entry but ultimately all the filters will clear out
    // their dummy entries and we should be all OK.
    ESTList* estList  = runtimeContext->getESTList();
    ASSERT( estList != NULL );
    estList->deleteLastESTs(dummyESTList.size());
    // Notify ESt listener about the change
    if (runtimeContext->getESTListListener() != NULL) {
        // Let parameter set manager know that a bunch of dummy ESTs have
        // been removed.
        ESTListListener *listener = runtimeContext->getESTListListener();
        listener->entriesRemoved(estList->size(),
                                 estList->size() + dummyESTList.size());
    }
    // Clear out our dummy EST references as we no longer need them.
    dummyESTList.clear();
}

bool
LCFilter::runFilter(const EST* est, int& clusterID) {
    // Setup the analyzer for processing and checking entries.
    ESTAnalyzer *analyzer = runtimeContext->getAnalyzer();
    ASSERT ( analyzer != NULL );
    // Setup the reference EST for the analyzer.
    analyzer->setReferenceEST(est);
    // Check one pattern after another...
    for(size_t i = 0; (i < dummyESTList.size()); i++) {
        // Get comparison metric.
        ASSERT( dummyESTList[i].first != NULL );
        float metric = analyzer->analyze(dummyESTList[i].first);
        // Check if the metric is good/bad
        if (analyzer->compareMetrics(metric, (float) threshold)) {
            // Log for debugging.
            // std::cout << "Filtered est " << EST::getEST(estIdx)->getID()
            //           << " - " << EST::getEST(estIdx)->getInfo() << "\n";
            // This est is similar to a dummy est. filter it out.
            clusterID = dummyESTList[i].second;
            return true;
        }
    }
    // This EST does not have a low complexity section
    clusterID = -1;
    return false;
}

#endif
