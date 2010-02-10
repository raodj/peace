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
#include "ClusterMaker.h"
#include "ESTAnalyzer.h"
#include "EST.h"

#include <sstream>

// Define the static parameters
char LCFilter::DefaultPatternList[16] = "A,C";
char* LCFilter::patternList = LCFilter::DefaultPatternList;
int   LCFilter::threshold   = -1;

// The set of arguments for this class.
arg_parser::arg_record LCFilter::argsList[] = {
    {"--lcPatterns", "List of (, separated) patterns to generate dummy ESTs",
     &LCFilter::patternList, arg_parser::STRING},
    {"--lctTreshold", "Threshold value to detect low complexity sequences",
     &LCFilter::threshold, arg_parser::INTEGER},
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

LCFilter::LCFilter(ClusterMaker *clusterMaker) :
    Filter("lcFilter", clusterMaker) {
    // Nothing else to be done for now.
}

void
LCFilter::showArguments(std::ostream& os) {
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(LCFilter::argsList);
    os << ap;
}

bool
LCFilter::parseArguments(int& argc, char **argv) {
    // Let's process parameters now.
    arg_parser ap(LCFilter::argsList);
    ap.check_args(argc, argv, false);
    if (patternList == NULL) {
        std::cerr << filterName << ": Pattern list must be specified "
                  << "(use --lcPatterns option)\n";
        return false;
    }
    // Everything went well
    return true;
}

int
LCFilter::initialize() {
    ASSERT ( clusterMaker != NULL );
    ASSERT ( clusterMaker->getAnalyzer() != NULL );
    ASSERT ( patternList != NULL );
    ASSERT ( DefaultPatternList != NULL );

    // Check and update our threshold value if one is not specified.
    if (threshold == -1) {
        // Use the invalid metric from analyzer instead.
        threshold = clusterMaker->getAnalyzer()->getInvalidMetric();
    }
    // Figure out the length of the dummy ESTs
    const int dummyLen = clusterMaker->getAnalyzer()->getPreferredDummyESTLength();
    // Process the patternList and add a dummy cluster for each pattern.
    std::string patStr(patternList);
    // Process one word at a time. Words are separated by a ","
    // (comma) character.
    while (!patStr.empty()) {
        // Locate the next hyphen character and get name of filter
        const std::string::size_type commaPos = patStr.find(',');
        const std::string pattern = patStr.substr(0, commaPos);
        // Create dummy entry for this pattern.
        addDummyEntry("DummyEST For Pattern " + pattern, pattern, dummyLen);
        // Onto the next entry if one is present.
        if (commaPos == std::string::npos) {
            // No hyphen, so this was the last pattern in the list
            patStr.clear();
        } else {
            // patStr becomes the remaining patterns int the list
            patStr = patStr.substr(commaPos + 1);
        }
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
    const int repeats   = length / seq.length();
    for(int i = 0; (i < repeats); i++) {
        fullSeq += seq;
    }
    // Add dummy EST and cluster and save information.
    int estIdx = EST::getESTCount();
    EST *est = EST::create(estIdx, fastaID.c_str(), fullSeq.c_str());
    // Flag EST as having already been processed.
    est->setProcessed(true);
    // Add a dummy cluster for this filter with a suitable name.
    std::ostringstream clsName;
    clsName << "Low Complexity ESTs (filtered by LCFilter Pattern "
            << seq << "...)" << std::ends;
    int clusterID = clusterMaker->addDummyCluster(clsName.str());
    // Save information about the dummy entry.
    dummyESTList.push_back(DummyESTInfo(estIdx, clusterID));
}

void
LCFilter::finalize() {
    // First remove all the dummy ESTs we may have added. Note that
    // this one does not really guarantee that we remove the
    // appropiate entry but ultimately all the filters will clear out
    // their dummy entries and we should be all OK.
    EST::deleteLastESTs(dummyESTList.size());
    // Clear out our dummy EST references as we no longer need them.
    dummyESTList.clear();
}

int
LCFilter::runFilter(const int estIdx) {
    // Setup the analyzer for processing and checking entries.
    ASSERT ( clusterMaker != NULL );
    ESTAnalyzer *analyzer = clusterMaker->getAnalyzer();
    ASSERT ( analyzer != NULL );
    // Setup the reference EST for the analyzer.
    analyzer->setReferenceEST(estIdx);
    // Check one pattern after another...
    for(size_t i = 0; (i < dummyESTList.size()); i++) {
        // Get comparison metric.
        float metric = analyzer->analyze(dummyESTList[i].first);
        // Check if the metric is good/bad
        if (analyzer->compareMetrics(metric, threshold)) {
            // Log for debugging.
            // std::cout << "Filtered est " << EST::getEST(estIdx)->getID()
            //           << " - " << EST::getEST(estIdx)->getInfo() << "\n";
            // This est is similar to a dummy est. filter it out.
            return dummyESTList[i].second;
        }
    }
    // This EST does not have a low complexity section
    return -1;
}

#endif
