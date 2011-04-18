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

/** \file UsingPeace

    \brief A example on how to use PEACE library and some of its
    features.

    This C++ source code provides a straightforward example
    illustrating the use of PEACE library to perform some
    non-traditional task.  Often it may be more convenient to directly
    use the interactive mode of PEACE for performing most operations.
*/

#include "ESTAnalyzerFactory.h"
#include "HeuristicFactory.h"
#include "HeuristicChain.h"
#include "ParameterSetManager.h"

#include "ESTAnalyzer.h"
#include "Heuristic.h"
#include "ESTList.h"

/**  The main method to illustrate a non-standard use of PEACE
     library.

     The following method illustrates the following general steps that
     are involved in using PEACE:

     <ol>

     <li>Step 1: Create and populate the list of cDNA fragments to be
     processd in a ESTList.  In this example, we use a pair of static
     cDNA fragments.  However, the data can be generated as well.</li>

     <li>Step 2: This is an optional step of setting up heuristics to
     speedup analysis of cDNA fragments. PEACE does not require the
     use of heuristics.</li>

     <li>Step 3: The analyzer to be used is created via the
     ESTAnalyzerFactory.  Once the analyzer is created, the heuristic
     chain and the estList references are setup for its use.</li>

     </ol>
*/
int
main(int argc, char* argv[]) {
    // First let's populate the list of ESTs with some sample
    // ESTs. The first number is the index of the EST and must
    // monotonically increase. The second parameter is some identifier
    // (similar to a FASTA GI). The third parameter is the sequence
    // itself.

    // NOTE: The sequences must be at least 100 base pairs long. This
    // length is what the system has been designed for. If your ESTs
    // are shorter than 100 base pairs then your milage will vary.
    ESTList estList;
    estList.add(0, std::string("TestSeq#1"),
                "AAAAATTTTCCCGGGAAATTTCCGGGAAATTTCCCGGAAAAA"            \
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    estList.add(1, std::string("TestSeq#2"),
                "AAAAATTTTCCCGGGAAATTTCCGGGAAATTTCCCGGAAAAA"            \
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    // Create a heuristic chain to hold the uv and tv heuristics that
    // we are going to use. We could have created these two heuristics
    // separately. But the chain eases creating heuristics for us.
    HeuristicChain chain;
    chain.setupChain("uv-tv");
    chain.setESTList(&estList);
    
    // Next create the d2 EST analyzer.  In the example below, we will
    // use std::auto_ptr<> that will automatically delete the pointers
    // once they go out of scope. That way we don't have to worry
    // about freeing the pointers later on -- basically gives an
    // automatic-smart garbage collection feature.
    std::auto_ptr<ESTAnalyzer> d2(ESTAnalyzerFactory::create("d2"));
    // Configure the properties of the analyzer. For this we need to
    // know.  The properties of the specific analyzer are represented
    // as <name, value> pairs.  The names are specified by a short
    // string.  A convenient way to obtain and display valid
    // properties is shown (and commented out) below:
    // d2->showArguments(std::cout);

    // Setup properties of d2 analyzer
    d2->setArgument("--frame", 50);
    d2->setArgument("--word",  6);
    
    // Setup the heuristics and est list for the analyzer to work with.
    d2->setESTList(&estList);
    d2->setHeuristicChain(&chain);
    
    // Initialize the analyzer (it also initializes the heuristic chain).
    if (!d2->initialize()) {
        std::cerr << "Error initializing ESTAnalyzer " << d2->getName()
                  << ". Aborting.\n";
        return 1;
    }
    // Now let's first use the UV heuristic to see if the pairs of
    // ESTs pass UV heuristic.
    Heuristic* uv = chain.getHeuristic("uv");
    ASSERT( uv != NULL );
    // Get UV heuristic to analyze EST #0 (reference EST) and EST #1
    uv->setReferenceEST(estList.get(0));
    bool result = uv->shouldAnalyze(estList.get(1));
    std::cout << "UV heuristic says that EST #0 and EST #1 are "
              << (result ? "" : "not ") << "related." << std::endl;

    // Next let's use the TV heuristic to see if the pairs of ESTs
    // pass TV heuristic.
    Heuristic* tv = chain.getHeuristic("tv");
    ASSERT ( tv != NULL );
    // Get TV heuristic to analyze EST #0 (reference EST) and EST #1
    tv->setReferenceEST(estList.get(0));
    result = tv->shouldAnalyze(estList.get(1));
    std::cout << "TV heuristic says that EST #0 and EST #1 are "
              << (result ? "" : "not ") << "related." << std::endl;

    // Now let's run d2 analyzer and obtain d2 score between EST #0
    // (the reference EST set earlier) and EST #1
    d2->setReferenceEST(estList.get(0));
    const float metric = d2->analyze(estList.get(1));
    std::cout << "d2 score between EST #0 and EST #1 is: "
              << (int) metric << std::endl;
    // Everything went well.
    return 0;
}

#endif
