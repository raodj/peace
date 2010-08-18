#ifndef MAIN_CPP
#define MAIN_CPP

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

#include "arg_parser.h"
#include "ClusterMakerFactory.h"
#include "ESTAnalyzerFactory.h"
#include "AssemblerFactory.h"
#include "ClusterMaker.h"
#include "ESTAnalyzer.h"
#include "Assembler.h"
#include "HeuristicFactory.h"
#include "HeuristicChain.h"
#include "FilterFactory.h"
#include "FilterChain.h"
#include "MPIHelper.h"
#include "InteractiveConsole.h"
#include "EST.h"

#include <string>
#include <sstream>
#include <fstream>

/** \fn void showUsage(arg_parser& ap, ESTAnalyzer *analyzer, ClusterMaker *cMaker, Assembler* assembler, HeuristicChain *hChain, FilterChain* fChain)

    A simple helper method to show usage and clean up the temporarily
    created EST analyzer and cluster Maker objects.

    \note This method is meant to be in local scope rather than global
    scope.

    \param[in] ap The argument parser to be used for displaying
    common global options.
    
    \param[in] analyzer The EST analyzer (if any) that was created
    and must now be deleted after showing usage.  If this pointer is
    NULL, then this parameter is ignored.

    \param[in] cMaker The cluster Maker (if any) that was created
    and must now be deleted after showing usage.  If this pointer is
    NULL, then this parameter is ignored.

    \param[in] assembler The assembler (if any) that was created and
    must now be deleted after showing usage. If this pointer is NULL,
    then this parameter is ignored.
    
    \param[in] hChain The chain of heuristics that are currently
    defined for this run.

    \param[in] fChain The chain of filters that are currently defined
    for this run.
*/
static void showUsage(arg_parser& ap,
                      ESTAnalyzer *analyzer, ClusterMaker *cMaker,
                      Assembler* assembler, HeuristicChain *hChain,
                      FilterChain *fChain) {
    std::cout << "PEACE Version 0.95 (Released April 19 2010)\n"
              << "Usage: PEACE [options]\n"
              << "where options are:\n";
    std::cout << ap;
    std::cout << "Names of EST analyzers available are:\n";
    ESTAnalyzerFactory::displayList(std::cerr);
    std::cout << "Names of cluster makers available are:\n";
    ClusterMakerFactory::displayList(std::cerr);
    std::cout << "Names of gene assemblers available are:\n";
    AssemblerFactory::displayList(std::cerr);    
    std::cout << "Names of heuristics available are:\n";
    HeuristicFactory::displayList(std::cerr);
    std::cout << "Names of filters available are:\n";
    FilterFactory::displayList(std::cerr);
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
    // Display any assembler specific options.
    if (assembler != NULL) {
        assembler->showArguments(std::cout);
        delete assembler;
    }
    
    // Display any heuristic specific options.
    if (hChain != NULL) {
        hChain->showArguments(std::cout);
        delete hChain;
    }
    // Display any filter specific options.
    if (fChain != NULL) {
        fChain->showArguments(std::cout);
        delete fChain;
    }
}

/** \fn int applyFilters(ClusterMaker *, FilterChain *, const char *, const char *)

    \brief Helper method to apply filters (if any).

    This method is a helper method that is invoked from the main
    method to apply any filters that the user has sepecified to the
    list of ESTs (before they are actually clustered).

    \return This method returns 0 (zero) on success. On errors or if
    further clustering is not to be performed, then this method
    returns a non-zero value.
*/
int applyFilters(ClusterMaker *clusterMaker, FilterChain *chain,
                 const char* filterPassFileName,
                 const char *filterFailFileName) {
    if (chain == NULL) {
        // Nothing to be done.
        return 0;
    }
    if (FilterChain::applyFilters(clusterMaker)) {
        // Error occured when applying filter chains.
        return 1;
    }
    // Rest of the file creation code applies only to rank 0
    if (MPI_GET_RANK() != 0) {
        // Nothing further to be done on this process.
        return 0;
    }
    // Filter applied successfully. Generate list of processed ESTs
    // that passed and failed filtering.
    if (filterPassFileName != NULL) {
        std::ofstream outFile(filterPassFileName);
        if (!outFile.good()) {
            std::cerr << "Unable to open file " << filterPassFileName
                      << "for writing. Aborting!" << std::endl;
            return 2;
        }
        // Dump the ESTs out
        EST::dumpESTList(outFile, false);
        outFile.close();
    }
    // Dump out failed ESTs.
    if (filterFailFileName != NULL) {
        std::ofstream outFile(filterFailFileName);
        if (!outFile.good()) {
            std::cerr << "Unable to open file " << filterFailFileName
                      << "for writing. Aborting!" << std::endl;
            return 2;
        }
        // Dump the ESTs that failed filtering out.
        EST::dumpESTList(outFile, true);
        outFile.close();
    }
    // Everything went well.
    return 0;
}

/** \fn int main(int argc, char *argv[])

    \brief The main function that launches all the activites.

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
    char emptyString[1]     = {'\0'};
    char defAnalyzer[]      = "twopassD2adapt";
    char defClusterMaker[]  = "mst";
    char defHeuristic[]     = "tv";
    char defFilters[]       = "lengthFilter-lcFilter";

    char *assemblerName = NULL;
    char *analyzerName  = NULL;
    char *clusterName   = NULL;
    char *heuristicStr  = defHeuristic;
    char *outputFile    = emptyString;
    bool showOptions    = false;
    int  refESTidx      = 0;
    bool interactive    = false;
    
    // Filtering parameters
    char *filterStr      = defFilters;
    char *filterPassFile = NULL;
    char *filterFailFile = NULL;
    bool filterOnly      = false;

    // Create the list of valid arguments to be used by the arg_parser.
    arg_parser::arg_record arg_list[] = {
        {"--clusterMaker", "Name of clustering algorithm to use (null for none)",
         &clusterName, arg_parser::STRING},
        {"--assembler", "Name of the assembler to use",
         &assemblerName, arg_parser::STRING},
        {"--analyzer", "Name of the EST analyzer to use",
         &analyzerName, arg_parser::STRING},
        {"--heuristics", "Name(s) of the heuristic(s) to use, in order (null for none)",
         &heuristicStr, arg_parser::STRING},
        {"--filters", "Name(s) of the filters(s) to use, in order (null for none)",
         &filterStr, arg_parser::STRING},
        {"--filter-fail", "FASTA file to which ESTs that failed the filter must be written",
         &filterFailFile, arg_parser::STRING},
        {"--filter-pass", "FASTA file to which ESTs that pass the filter must be written",
         &filterPassFile, arg_parser::STRING},
        {"--filter-only", "Just do filtering and no other processing (false by default)",
         &filterOnly, arg_parser::BOOLEAN},
        {"--estIdx", "Index of reference EST in a EST file",
         &refESTidx, arg_parser::INTEGER},
        {"--output", "File to which output must be written",
         &outputFile, arg_parser::STRING},              
        {"--options", "Lists options for the specified cluster & analyzer",
         &showOptions, arg_parser::BOOLEAN},
        {"--interactive", "Lauch PEACE interactive console",
         &interactive, arg_parser::BOOLEAN},
        {NULL, NULL, NULL, arg_parser::BOOLEAN}
    };
   
    // Perform any mpi initialization as needed
    MPI_INIT(argc, argv);
    // Save global arguments for future reference to generate logs
    arg_parser::set_global_args(argc, argv);
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
    // Check to see if an assembler or cluster-maker has been
    // specified and based on that setup the cluster maker.
    if ((clusterName == NULL) && (assemblerName == NULL)) {
        // By default we run the default cluster maker.
        clusterName = defClusterMaker;
    }
    // Check and setup the default analyzer as needed.
    if ((clusterName != NULL) && (analyzerName == NULL)) {
        // By default we run the default analyzer.
        analyzerName = defAnalyzer;
    }
    
    // Create an EST analyzer using the analyzer name, if one is set
    ESTAnalyzer *analyzer = 
        ESTAnalyzerFactory::create(analyzerName, refESTidx,
                                   std::string(outputFile));

    // Check for null clustermaker input
    if ((clusterName != NULL) && (!strcmp(clusterName, "null"))) {
        clusterName = NULL;
    }
    // Create an cluster generator using clusterName
    ClusterMaker *clusterMaker =
        ClusterMakerFactory::create(clusterName, analyzer,
                                    refESTidx, std::string(outputFile));
    
    // Create an gene assembler using the assemblerName.  If assembler
    // name is NULL, then Assembler object is also NULL.
    Assembler *assembler = AssemblerFactory::create(assemblerName,
                                                    outputFile, MPI_GET_RANK());
    
    // Check for null heuristic chain input
    if (!strcmp(heuristicStr, "null")) {
        heuristicStr = NULL;
    }
    // Create the heuristic chain using a helper method.
    HeuristicChain *heuristicChain =
        HeuristicChain::setupChain(heuristicStr, refESTidx, outputFile);

    // Check and create filters as directed.
    if (!strcmp(filterStr, "null")) {
        filterStr = NULL;
    }
    // Create the filter chain using a helper method.
    FilterChain *filterChain = FilterChain::setupChain(filterStr, clusterMaker);
    // See if user just want's list of options
    if (showOptions) {
        showUsage(ap, analyzer, clusterMaker, assembler,
                  heuristicChain, filterChain);
        return 0;
    }
    // Check if analyzer name was valid and the arguments to it are valid.
    if (analyzerName != NULL) {
        if ((analyzer == NULL) || (!analyzer->parseArguments(argc, argv))) {
            showUsage(ap, analyzer, clusterMaker, assembler,
                      heuristicChain, filterChain);
            return 1;
        }
    }
    // Check if cluster maker was specified and successfully created.
    if (clusterName != NULL) {
        if ((clusterMaker == NULL) ||
            (!clusterMaker->parseArguments(argc, argv))) {
            showUsage(ap, analyzer, clusterMaker, assembler,
                      heuristicChain, filterChain);
            return 2;
        }
    }
    // Check if assembler name was valid and the arguments to it are valid.
    if (assemblerName != NULL) {
        if ((assembler == NULL) || (!assembler->parseArguments(argc, argv))) {
            showUsage(ap, analyzer, clusterMaker, assembler,
                      heuristicChain, filterChain);
            return 4;
        }
    }    
    // Check if heuristic chain was specified and successfully created
    if (heuristicStr != NULL) {
        if ((heuristicChain == NULL) ||
            (!heuristicChain->parseArguments(argc, argv))) {
            showUsage(ap, analyzer, clusterMaker, assembler,
                      heuristicChain, filterChain);
            return 3;
        }
    }
    
    // Check if filter chain was specified and successfully created
    if (filterStr != NULL) {
        if ((filterChain == NULL) || (showOptions) ||
            (!filterChain->parseArguments(argc, argv))) {
            showUsage(ap, analyzer, clusterMaker, assembler,
                      heuristicChain, filterChain);
            return 4;
        }
    }
    // Attach the heuristic chain (if it exists) to the EST analyzer
    if ((heuristicChain != NULL) && (analyzer != NULL)) {
        analyzer->setHeuristicChain(heuristicChain);
    }
    
    // Get the cluster maker or est analyzer to analyze.
    int result = 0;
   
    if (interactive) {
        // Launch PEACE interactive console.
        InteractiveConsole console(analyzer);
        console.processCommands();
    } else if (clusterMaker != NULL) {
        // First initialize the cluster maker for use in filters
        if ((result = clusterMaker->initialize()) == 0) {
            // First go through the phase of applying filters if specified.
            if ((applyFilters(clusterMaker, filterChain, filterPassFile,
                              filterFailFile) == 0) && (!filterOnly)) {
                result = clusterMaker->makeClusters();
            }
        }
        delete clusterMaker;
    } else if (assembler != NULL) {
        // First initialize the assembler for use in filters
        if ((result = assembler->initialize()) == 0) {
            // Let the assembler do its assembly.
            result = assembler->assemble();
        }
        delete assembler;
    } else {
        // Let the analyzer do the analysis.
        result = analyzer->analyze();
    }

    // Print statistics regarding operation of heuristics. Maybe the
    // printing of stats must be done based on a command line
    // argument.
    if ((result == 0) && (heuristicChain != NULL)) {
        // Display statistics by writing all the data to a string
        // stream and then flushing the stream.  This tries to working
        // around (it is not perfect solution) interspersed data from
        // multiple MPI processes
        std::ostringstream buffer;
        heuristicChain->printStats(buffer, MPI_GET_RANK());
        std::cout << buffer.str() << std::endl;
    }

    // Print statistics regarding operation of filters.
    if ((result == 0) && (filterChain != NULL)) {
        // Display statistics by writing all the data to a string
        // stream and then flushing the stream.  This tries to working
        // around (it is not perfect solution) interspersed data from
        // multiple MPI processes
        std::ostringstream buffer;
        filterChain->printStats(buffer, MPI_GET_RANK());
        std::cout << buffer.str() << std::endl;
    }
    
    // Delete the analyzer as it is no longer needed.
    delete analyzer;
    // Delete the heuristic chain as it is no longer needed
    if (heuristicChain != NULL) {
        delete heuristicChain;
    }
    // Delete the filter chain as it is no longer needed
    if (filterChain != NULL) {
        delete filterChain;
    }

    // Shutdown MPI.
    MPI_FINALIZE();
    
    // Return result of analysis (usually 0 to indicate no errors).
    return result;
}

#endif
