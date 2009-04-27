#ifndef MAIN_CPP
#define MAIN_CPP

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include "arg_parser.h"
#include "ClusterMakerFactory.h"
#include "ESTAnalyzerFactory.h"
#include "ClusterMaker.h"
#include "ESTAnalyzer.h"
#include "HeuristicFactory.h"
#include "HeuristicChain.h"
#include "InteractiveConsole.h"

#include <mpi.h>
#include <string>

/** \func showUsage

    A simple helper method to show usage and clean up the temporarily
    created EST analyzer and cluster Maker objects.

    \note This method is meant to be in local scope rather than global
    scope.

    \param[in,out] ap The argument parser to be used for displaying
    common global options.
    
    \param[in,out] analyzer The EST analyzer (if any) that was created
    and must now be deleted after showing usage.  If this pointer is
    NULL, then this parameter is ignored.

    \param[in,out] cMaker The cluster Maker (if any) that was created
    and must now be deleted after showing usage.  If this pointer is
    NULL, then this parameter is ignored.
*/
static void showUsage(arg_parser& ap,
                      ESTAnalyzer *analyzer, ClusterMaker *cMaker,
                      HeuristicChain *hChain) {
    std::cout << "Usage: PEACE [options]\n"
              << "where options are:\n";
    std::cout << ap;
    std::cout << "Names of EST analyzers available are:\n";
    ESTAnalyzerFactory::displayList(std::cerr);
    std::cout << "Names of cluster makers available are:\n";
    ClusterMakerFactory::displayList(std::cerr);
    std::cout << "Names of heuristics available are:\n";
    HeuristicFactory::displayList(std::cerr);
    // Display any analyzer specific options.
    if (analyzer != NULL) {
        analyzer->showArguments(std::cout);
        delete analyzer;
    }
    // Display any cluster maker specific options.
    if (cMaker != NULL) {
        cMaker->showArguments(std::cout);
        delete cMaker;
    }
    // Display any heuristic specific options.
    if (hChain != NULL) {
        hChain->showArguments(std::cout);
        delete hChain;
    }
}

/** \func main

    <p> The global main method for the ESTAnalyzer software. The main
    method performs the task of identifying and creating a specific
    ESTAnalyzer.  Accordingly, it only consumes the "--analyzer"
    option which is typically the first option specified.  However, it
    is not required that --analyzer be the first option. </p>

    <p>Once a specific EST Analyzer has been instantiated the main
    method simply passes control over to that EST Analyzer.</p>
*/
int
main(int argc, char* argv[]) {
    // Values of the following variables are processed by the argument
    // parser further below.
    char emptyString[1]={'\0'};
    char *analyzerName = NULL;
    char *clusterName  = NULL;
    char *heuristicStr = NULL;
    char *outputFile   = emptyString;
    bool showOptions   = false;
    int  refESTidx     = -1;
    bool interactive   = false;
    
    // Create the list of valid arguments to be used by the arg_parser.
    arg_parser::arg_record arg_list[] = {
        {"--clusterMaker", "Name of clustering algorithm to use",
         &clusterName, arg_parser::STRING},
        {"--analyzer", "Name of the EST analyzer to use",
         &analyzerName, arg_parser::STRING},
        {"--heuristics", "Name(s) of the heuristic(s) to use, in order",
         &heuristicStr, arg_parser::STRING},
        {"--estIdx", "Index of reference EST in a EST file",
         &refESTidx, arg_parser::INTEGER},
        {"--output", "File to which output must be written",
         &outputFile, arg_parser::STRING},              
        {"--options", "Lists options for the specified cluster & analyzer",
         &showOptions, arg_parser::BOOLEAN},
        {"--interactive", "Lauch PEACE interactive console",
         &interactive, arg_parser::BOOLEAN},
        {NULL, NULL}
    };
   
    // Perform any mpi initialization as needed
    MPI::Init(argc, argv);
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    arg_parser ap(arg_list);
    ap.check_args(argc, argv, false);

    if (showOptions) {
        // If showOptions setting is true, in this case no real
        // processing is being done other than to display options.  In
        // this case, set refESTidx to a dummy value to enable
        // creation of default configurations of tools below.
        refESTidx = 0;
    }
    // Create an EST analyzer using the analyzer name.
    ESTAnalyzer *analyzer =
        ESTAnalyzerFactory::create(analyzerName, refESTidx,
                                   std::string(outputFile));
    // Create an cluster generating using clusterName;
    ClusterMaker *clusterMaker =
        ClusterMakerFactory::create(clusterName, analyzer,
                                    refESTidx, std::string(outputFile));
    // Create the heuristic chain using a helper method.
    HeuristicChain *heuristicChain =
        HeuristicChain::setupChain(heuristicStr, refESTidx, outputFile);    
    // Check if EST analyzer creation was successful.  A valid EST
    // analyzer is needed even to make clusters.
    if ((analyzer == NULL) || (showOptions) ||
        (!analyzer->parseArguments(argc, argv))) {
        showUsage(ap, analyzer, clusterMaker, heuristicChain);
        return 1;
    }

    // Check if cluster maker was specified and successfully created.
    if (clusterName != NULL)  {
        if ((clusterMaker == NULL) || (showOptions) ||
            (!clusterMaker->parseArguments(argc, argv))) {
            showUsage(ap, analyzer, clusterMaker, heuristicChain);
            return 2;            
        }
    }

    // Check if heuristic chain was specified and successfully created
    if (heuristicStr != NULL) {
        if ((heuristicChain == NULL) || (showOptions) ||
            (!heuristicChain->parseArguments(argc, argv))) {
            showUsage(ap, analyzer, clusterMaker, heuristicChain);
            return 3;
        }
    }

    // Attach the heuristic chain (if it exists) to the EST analyzer
    if (heuristicChain != NULL) {
        analyzer->setHeuristicChain(heuristicChain);
    }
    
    // Get the cluster maker or est analyzer to analyze.
    int result = 0;
    if (clusterMaker != NULL) {
        result = clusterMaker->makeClusters();
        delete clusterMaker;
    } else if (interactive) {
        // Launch PEACE interactive console.
        InteractiveConsole console(analyzer);
        console.processCommands();
    } else {
        // Let the analyzer do the analysis.
        result = analyzer->analyze();
    }
    // Delete the analyzer as it is no longer needed.
    delete analyzer;

    // Delete the heuristic chain as it is no longer needed
    delete heuristicChain;

    // Shutdown MPI.
    MPI::Finalize();
    
    // Return result of analysis (usually 0 to indicate no errors).
    return result;
}

#endif
