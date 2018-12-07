#ifndef CLUSTERING_SUB_SYSTEM_CPP
#define CLUSTERING_SUB_SYSTEM_CPP

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

#include "ClusteringSubSystem.h"
#include "ESTAnalyzerFactory.h"
#include "ClusterMakerFactory.h"
#include "HeuristicFactory.h"
#include "ClusterMaker.h"
#include "ESTAnalyzer.h"
#include "RuntimeContext.h"
#include "ESTList.h"

ClusteringSubSystem::ClusteringSubSystem() {
    heuristicChain.setSubSystem(this);
    estAnalyzer      = NULL;
    clusterMaker     = NULL;
    // Setup default values for parameters
    clusterMakerName = "null";
    estAnalyzerName  = "twoPassD2adapt";
}

ClusteringSubSystem::~ClusteringSubSystem() {
    // Nothing to be done (at least for now)
}

void
ClusteringSubSystem::addCommandLineArguments(ArgParser& argParser) {
    // Add command-line parameters and additional information
    // regarding heuristics.
    const ArgParser::ArgRecord HeuristicsArgs[] = {
        {"", "\nValid arguments for the clustering sub-system are:",
         NULL, ArgParser::INFO_MESSAGE},
        {"--heuristics", "Names of the heuristics to use in order (null for none)."
         "Valid heurisitcs are:",
         &heuristicNames, ArgParser::STRING_LIST},
        {"", "", NULL, ArgParser::INVALID}        
    };
    argParser.addValidArguments(HeuristicsArgs);
    // Have the factory add informational entries regarding the list
    // of valid heuristics
    HeuristicFactory::addCommandLineInfo(argParser);

    // Add command-line parameters and additional information
    // regarding EST analyzers
    const ArgParser::ArgRecord AnalyzerArgs[] = {
        {"--analyzer", "EST analyzer to use from list below:",
         &estAnalyzerName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(AnalyzerArgs);
    // Have the factory add informational entries regarding the list
    // of valid heuristics
    ESTAnalyzerFactory::addCommandLineInfo(argParser);

    // Add command-line parameters and additional information
    // regarding cluster makers
    const ArgParser::ArgRecord ClusterMakerArgs[] = {    
        {"--clusterMaker", "Name of the cluster maker to use. Valid "
         "cluster makers are:", &clusterMakerName, ArgParser::STRING},        
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add the valid arguments to the argument parser
    argParser.addValidArguments(ClusterMakerArgs);
    // Have the factory add informational entries regarding the list
    // of valid heuristics
    ClusterMakerFactory::addCommandLineInfo(argParser);
}

int
ClusteringSubSystem::initializeSubSystem(ArgParser& argParser) {
    // Check and setup default heuristic entries if none have been specified.
    if (heuristicNames.empty()) {
        heuristicNames.push_back("tv");
    }
    // Process one heuristic name at a time.
    for(size_t hurIdx = 0; (hurIdx < heuristicNames.size()); hurIdx++) {
        if (heuristicNames[hurIdx] == "null") {
            // A dummy heuristic name to be ignored.
            break;
        }
        // Create heuristic and add it to the chain
        Heuristic *heuristic =
            HeuristicFactory::create(heuristicNames[hurIdx], &heuristicChain);
        if (heuristic == NULL) {
            // Break out and return null, which will cause an error
            // and the usage message to show.
            return 1;
        }
        // Have the filter add its own parameters to the list.
        heuristic->addCommandLineArguments(argParser);
        // Add the filter to our chain.
        heuristicChain.addHeuristic(heuristic);
    }
    // Now create the EST analyzer (if user has indicated one)
    if (estAnalyzerName != "null") {
        if ((estAnalyzer = ESTAnalyzerFactory::create(estAnalyzerName)) == NULL) {
            // Invalid est analyzer name.
            return 2;
        }
        // Setup cross references in the analyzer.
        estAnalyzer->setSubSystem(this);
        estAnalyzer->setHeuristicChain(&heuristicChain);
        // Obtain additional parameters from it.
        estAnalyzer->addCommandLineArguments(argParser);
        // Setup pointer in the runtime context.
        runtimeContext->setAnalyzer(estAnalyzer);
    }
    
    // Now create the cluster maker (if user has indicated one)
    if (clusterMakerName != "null") {
        clusterMaker = ClusterMakerFactory::create(clusterMakerName,
                                                   estAnalyzer);
        if (clusterMaker == NULL) {
            // Invalid cluster maker name.
            return 3;
        }
        // Setup cross references in the cluster maker.
        clusterMaker->setSubSystem(this);
        estAnalyzer->setHeuristicChain(&heuristicChain);
        // Obtain additional parameters from it.
        clusterMaker->addCommandLineArguments(argParser);
        // Setup pointer in the runtime context.
        runtimeContext->setClusterMaker(clusterMaker);
    }
    
    // Everything went well.
    return NO_ERROR;
}

int
ClusteringSubSystem::initializeSubComponents() {
    ASSERT ( runtimeContext != NULL );    
    // setup the EST list for the heuristic chain to use
    heuristicChain.setESTList(runtimeContext->getESTList());
    // setup the EST list for the analyzer to use
    if (estAnalyzer != NULL) {
        estAnalyzer->setESTList(runtimeContext->getESTList());
        // Initialize the EST analyzer first. This initialization has
        // to occur first as the AdaptiveTwoPassD2 set's up parameter
        // manager for the heuristic chains to use.  The below call
        // also initializes the heuristic chain in
        // ESTAnalyzer::initialize() method.
        if (!estAnalyzer->initialize()) {
            return 2;
        }
    }
    // Initialize the cluster maker.
    if (clusterMaker != NULL) {
        if (!clusterMaker->initialize()) {
            return 3;
        }
    }
    // Everything went well
    return NO_ERROR;
}

int
ClusteringSubSystem::run() {
    if (clusterMaker == NULL) {
        // Nothing to be done.
        return NO_ERROR;
    }
    // Initialize the EST analyzer first. This initialization has to
    // occur first as the AdaptiveTwoPassD2 set's up parameter manager
    // for the heuristic chains to use.  The below call also
    // initializes the heuristic chain in ESTAnalyzer::initialize()
    // method.
    if (!estAnalyzer->initialize()) {
        return 2;
    }
    // Initialize the cluster maker.
    if (!clusterMaker->initialize()) {
        return 3;
    }    
    // Check to ensure we have some ESTs to process.
    if ((runtimeContext->getESTList() == NULL) ||
        (runtimeContext->getESTList()->size() < 1)) {
        // Don't have ESTs to process.
        std::cerr << "Error: The list of cDNAs to be processed is empty.\n"
                  << "Use --estFiles option to load data from file(s).\n";
        return 5;
    }
    // Let the heuristics run (and do any global operations on list of
    // reads) as needed.
    heuristicChain.run();
    // Now let the cluster maker run and do clustering
    if (clusterMaker->makeClusters() != 0) {
        // Error during clustering.
        return 4;
    }
    // Everything went on well
    return NO_ERROR;
}

void
ClusteringSubSystem::finalizeComponents(const bool UNREFERENCED_PARAMETER(success)) {
    // Finalize all the components.
    if (clusterMaker != NULL) {
        clusterMaker->finalize();
    }
    if (estAnalyzer != NULL) {
        estAnalyzer->finalize();
    }
    heuristicChain.finalize();
}

void
ClusteringSubSystem::finalizeSubSystem(const bool UNREFERENCED_PARAMETER(success)) {
    // Now get rid of the dynamically allocated components
    delete clusterMaker;
    delete estAnalyzer;
    // Reset pointers to make deletion explicitly trackable (in valgrind)
    clusterMaker = NULL;
    estAnalyzer  = NULL;
}

#endif
