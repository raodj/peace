#ifndef USING_PEACE_CPP
#define USING_PEACE_CPP

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

/** A example on how to use PEACE library and its features.


*/

#include "ClusterMakerFactory.h"
#include "ESTAnalyzerFactory.h"
#include "HeuristicFactory.h"
#include "HeuristicChain.h"

#include "ClusterMaker.h"
#include "ESTAnalyzer.h"
#include "Heuristic.h"
#include "EST.h"

int
main(int argc, char* argv[]) {
    // First let's populate the list of ESTs with some sample
    // ESTs. The first number is the index of the EST and must
    // monotonically increase. The second parameter is some identifier
    // (similar to a FASTA GI). The third parameter is the sequence
    // itself.

    // NOTE: The sequences must be at least 100 base pairs long. This
    // length is what the system has been designed for. If your ESTs
    // are shorter than 100 base paris then your milage will vary.
    EST::create(0, "TestSeq#1", "AAAAATTTTCCCGGGAAATTTCCGGGAAATTTCCCGGAAAAA" \
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    EST::create(1, "TestSeq#2", "AAAAATTTTCCCGGGAAATTTCCGGGAAATTTCCCGGAAAAA" \
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    
    // Create a heuristic chain to hold the uv and tv heuristics that
    // we are going to use. We could have created these two heuristics
    // separately. But the chain eases creating heuristics for us.
    HeuristicChain* chain = HeuristicChain::setupChain("uv-tv", 0, "");
    // Let initialize and setup both heuristics via the chain
    chain->initialize();
    chain->setReferenceEST(0); // Set reference EST using index value.

    // Next create the d2 EST analyzer.  In the example below, we will
    // use std::auto_ptr<> that will automatically delete the pointers
    // once they go out of scope. That way we don't have to worry
    // about freeing the pointers later on -- basically gives an
    // automatic-smart garbage collection feature.
    std::auto_ptr<ESTAnalyzer> d2(ESTAnalyzerFactory::create("d2", 0, ""));
    // Set some parameters. We need a temporary location to hold
    // constant strings (even though we don't mutate it as the API
    // requires mutable strings)
    char param0[]  = "./UsingPeace";
    char param1[]  = "--frame",   value1[] = "50";
    char param2[]  = "--word",    value2[] = "6";
    char param3[]  = "--estFile", value3[] = "<none>";
    char* params[] = {param0,  // First param is always executable name.
                      param1, value1, param2, value2, // Set Window & Word size
                      param3, value3};  // Last parameter is mandatory
    int paramCount = sizeof(params) / sizeof(char*);
    d2->parseArguments(paramCount, params);

    // Now let's first use the UV heuristic to see if the pairs of
    // ESTs pass UV heuristic.
    Heuristic* uv = chain->getHeuristic("uv");
    // Get UV heuristic to analyze EST #0 (reference EST set earlier
    // on in this example) and EST #1
    bool result = uv->shouldAnalyze(1);
    std::cout << "UV heuristic says that EST #0 and EST #1 are "
              << (result ? "" : "not ") << "related." << std::endl;

    // Next let's use the TV heuristic to see if the pairs of ESTs
    // pass TV heuristic.
    Heuristic* tv = chain->getHeuristic("tv");
    // Get TV heuristic to analyze EST #0 (reference EST set earlier
    // on in this example) and EST #1
    result = tv->shouldAnalyze(1);
    std::cout << "TV heuristic says that EST #0 and EST #1 are "
              << (result ? "" : "not ") << "related." << std::endl;

    // Now let's run d2 analyzer and obtain d2 score between EST #0
    // (the reference EST set earlier) and EST #1
    d2->initialize();
    d2->setReferenceEST(0);
    const float metric = d2->analyze(1);
    std::cout << "d2 score between EST #0 and EST #1 is: "
              << (int) metric << std::endl;

    // Free up our chain.
    delete chain;
    // Everything went well.
    return 0;
}

#endif
