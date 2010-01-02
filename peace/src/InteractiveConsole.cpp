#ifndef INTERACTIVE_CONSOLE_CPP
#define INTERACTIVE_CONSOLE_CPP

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

#include "InteractiveConsole.h"
#include "Utilities.h"
#include "ESTAnalyzer.h"
#include "EST.h"
#include "MPIHelper.h"

#include <algorithm>
#include <stdio.h>

#ifdef HAVE_LIBREADLINE 
#include <readline/readline.h>
#include <readline/history.h>
#else
#define add_history(x)
#endif

#include <cmath>
#include <iomanip>

#define INTRO                                                           \
    "\033[1;38mPEACE: Parallel EST Analyzer and Clustering Engine\n"    \
    "\033[0;39mCopyright (c) Miami University, Oxford, OHIO."           \
    "All rights reserved.\n"                                            \
    "Type the command \033[1;32mhelp\033[0;39m to obtain list "         \
    "of commands to run\n"

// The set of command handler methods for processing various tags
InteractiveConsole::CmdEntry InteractiveConsole::cmdHandlerList[] = {
    {"exit",      &InteractiveConsole::exit},
    {"stats",     &InteractiveConsole::printStats},
    {"list",      &InteractiveConsole::list},
    {"analyze",   &InteractiveConsole::analyze},
    {"help",      &InteractiveConsole::help},
    {"print",     &InteractiveConsole::print},
    {"", NULL}
};

InteractiveConsole::InteractiveConsole(ESTAnalyzer* estAnalyzer) :
    analyzer(estAnalyzer) {
    ASSERT ( analyzer != NULL );
}

InteractiveConsole::~InteractiveConsole() {
    // Nothing to be done for now.
}

bool
InteractiveConsole::initialize() {
    // Print a standard intro..
    std::cout << INTRO << std::endl;    

    // Initialize the analyzer and load the ESTs from a FASTA file.
    // Collate and print time taken for initializaton.
    std::cout << "Loading ESTs..." << std::flush;
    const double startTime = MPI_WTIME();
    bool retVal = analyzer->initialize() == 0;
    // Track initialization time in milliseconds
    const double elapsedTime = (MPI_WTIME() - startTime) * 1000.0;
    if (!retVal) {
        // Error occured during initialization. Bail out.
        std::cerr << "Error initializing EST analyzer.\nExiting.\n";
    } else {
        std::cout << "Loaded " << EST::getESTCount() << " ESTs.\n"
                  << "* Elapsed time: " << elapsedTime << " msecs.\n\n";
    }
    // Let the caller know if initialization was successful.
    return retVal;
}

void
InteractiveConsole::processCommands() {
    // Initialize and load ESTs
    if (!initialize()) {
        // Error during initialization. Bail out
        return;
    }
    // Flag to indicate when "exit" command is read.
    bool done     = false;
    do {
        // Use the readline() function to interactively read command.
        char *cmdLine = readline("\033[1;33mpeace> ");
        std::cout << "\033[0;39m"; // Reset color
        if (cmdLine == NULL) {
            // User pressed Ctrl+D or commands were from a file.
            break;
        }
        // Tokenize the command
        std::vector<std::string> cmdWords = tokenize(cmdLine);
        // Save command in history buffer of readline for reference.
        if (cmdWords.size() > 0) {
            add_history(cmdLine);
        }
        // Free the unused cmdLine.
        free(cmdLine);
        if (cmdWords.size() < 1) {
            // Ignore empty command lines
            continue;
        }
        // Process the command. First conver the command to all lower
        // case and then process it.
        std::transform(cmdWords[0].begin(), cmdWords[0].end(),
                       cmdWords[0].begin(), tolower);
        // Search for command in the list of commands we know.
        int cmdIndex = 0;
        while ((cmdWords[0] != cmdHandlerList[cmdIndex].cmd) &&
               (cmdHandlerList[cmdIndex].handler != NULL)) {
            // On to the next command in teh cmd index.
            cmdIndex++;
        }
        // Dispatch call to the command handler
        if (cmdHandlerList[cmdIndex].handler != NULL) {
            (this->*cmdHandlerList[cmdIndex].handler)(cmdWords);
        } else {
            // Unknown command.
            std::cout << "huh?\n";
        }
        // Set the done flag depending on whether cmd is exit.
        done = (cmdIndex == 0);
    } while (!done);
}

std::vector<std::string>
InteractiveConsole::tokenize(const std::string& str,
                             const std::string& delims) {
    // The set of tokens
    std::vector<std::string> tokens;
    // skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delims, 0);
    // find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delims, lastPos);
    // Repeatedly extract strings and add them to the vector
    while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
        // found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delims, pos);
        // find next "non-delimiter"
        pos = str.find_first_of(delims, lastPos);
    }
    // Return the tokenized string
    return tokens;
}

void
InteractiveConsole::printStats(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords)) {
    // Collate and print statistics about the ESTs
    const double startTime = MPI_WTIME();
    // Obtain the list of ESTs to process
    const std::vector<EST*>& estList = EST::getESTList();
    // Count min, max, and mean length of ESTs
    size_t minLen = 0xffffffff, maxLen    = 0;
    double lenSum = 0, lenSumSqr = 0;
    for(size_t idx = 0; (idx < estList.size()); idx++) {
        // Obtain the length of the EST
        const size_t len = strlen(estList[idx]->getSequence());
        // Update the statistics
        minLen = std::min<size_t>(minLen, len);
        maxLen = std::max<size_t>(maxLen, len);
        // Compute data for mean and standard deviation computation
        lenSum    += len;
        lenSumSqr += (len * len);
    }
    // Print the statistics for display.
    const double meanLen = lenSum / estList.size();
    const double stdDev  = sqrt((lenSumSqr - estList.size() *
                                 meanLen * meanLen) / (estList.size() - 1));
    std::cout << "Number of ests       : " << estList.size() << std::endl;
    std::cout << "Mean EST size        : " << meanLen << " bp\n";
    std::cout << "Deviation in EST size: " << stdDev << " bp\n";
    std::cout << "Minimum EST size     : " << minLen << " bp\n";
    std::cout << "Maximum EST size     : " << maxLen << " bp\n";
    
    // Now compute and display time elapsed for this operation in
    // milliseconds
    std::cout << "* Elapsed time: "
              << ((MPI_WTIME() - startTime) * 1000.0) << " msecs\n";
}

void
InteractiveConsole::exit(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords)) {
    std::cout << "PEACE out.\n";
}

void
InteractiveConsole::list(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords)) {
    // Obtain the list of ESTs to print
    const std::vector<EST*>& estList = EST::getESTList();
    // Print the header for additional information.
    std::cout << "Index      FASTA Identifier\n"
              << "---------------------------------------------------------\n";
    // Print information about each EST
    for(size_t idx = 0; (idx < estList.size()); idx++) {
        std::cout << std::setw(8) << idx
                  << " : " << estList[idx]->getInfo()
                  << std::endl;
    }
    std::cout << "---------------------------------------------------------\n";
}

void
InteractiveConsole::analyze(const std::vector<std::string>& cmdWords) {
    if (cmdWords.size() != 3) {
        std::cout << "analyze command requires two arguments, the FASTA "
                  << "identifier or the index of the ESTs to compare.\n";
        return;
    }
    // Convert parameters to consistent index values
    const int index1 = getESTIndex(cmdWords[1]);
    const int index2 = getESTIndex(cmdWords[2]);
    if ((index1 == -1) || (index2 == -1)) {
        // invalid index values or FASTA identifiers.
        return;
    }
    // Print brief information about ESTs being analyzed to aid user
    const EST* const est1 = EST::getEST(index1);
    const EST* const est2 = EST::getEST(index2);
    std::cout << "Analyzing EST " << est1->getInfo() << " (index: " << index1
              << ") with\n"
              << "          EST "    << est2->getInfo() << " (index: " << index2
              << ")...\n";
    // Do the analysis part while tracking time taken.
    const double startTime = MPI_WTIME();
    // Setup the reference est for analysis.
    analyzer->setReferenceEST(index1);
    // Get the distance/similarity metric.
    float metric = analyzer->analyze(index2);
    // Compute time elapsed in milliseconds
    const double elapsedTime = (MPI_WTIME() - startTime) * 1000.0;
    // Display analysis results
    std::cout << "\033[1;38m" << analyzer->getName()
              << (analyzer->isDistanceMetric() ? " distance" : " similarity")
              << " metric: " << metric;
    // Print alignment information if the analyzer provides it.
    int alignmentData = 0;
    if (analyzer->getAlignmentData(alignmentData)) {
        std::cout << " (Alignment pos: " << alignmentData << ")";
    }
    std::cout << "\033[0;39m\n";
    // Print timing data as well.
    std::cout << "* Elapsed time: " << elapsedTime << " msecs.\n";
}

int
InteractiveConsole::getESTIndex(const std::string& id) const {
    // Obtain the list of ESTs to check
    const std::vector<EST*>& estList = EST::getESTList();
    const int  EstCount              = (int) estList.size();
    
    // Try and process id as a number.
    char *endptr = NULL;
    int index = strtol(id.c_str(), &endptr, 10);
    
    // Check if the index value was valid.
    if ((endptr == NULL) || (*endptr != '\0')) {
        // No it was not a number. Assume it is a fasta header and
        // search for it.
        for(index = 0; (index < EstCount); index++) {
            if (id == estList[index]->getInfo()) {
                // Found a match
                break;
            }
        }
    }
    // Ensure index is valid.
    if ((index < 0) || (index >= EstCount)) {
        std::cout << "Invalid EST index/identifier (" << id << ").\n";
        return -1;
    }
    // A valid index was specified
    return index;
}

void
InteractiveConsole::print(const std::vector<std::string>& cmdWords) {
    if (cmdWords.size() != 2) {
        std::cout << "The print command requires one argument, namely "
                  << "the FASTA identifier or the index of the EST "
                  << "whose information is to be printed.\n";
        return;
    }
    // Convert parameter to consistent index value
    const int index = getESTIndex(cmdWords[1]);
    if (index == -1) {
        // invalid index values or FASTA identifiers.
        return;
    }
    // Print information about EST
    const EST* const est = EST::getEST(index);
    std::cout << "Info on EST " << est->getInfo()
              << " (index: "    << index << "), Len: "
              << strlen(est->getSequence()) << " base pairs\n";
    std::cout << est->getSequence() << std::endl;
}

void
InteractiveConsole::help(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords)) {
    // Maybe this help method can load data from a help text file and
    // show the results. However, for now, this should suffice.
    std::cout << "Command   Description                                 \n"
              << "-------   --------------------------------------------\n";
    std::cout << "help      Shows this information.\n"
              << "list      List brief info. on all loaded ESTs.\n"
              << "print     Print information about a given EST.\n"
              << "analyze   Print d2/clu metrics for a given pair of ESTs.\n"
              << "stats     Print statistics about all the ESTs.\n"
              << "exit      Quit out of PEACE interactive console.\n";
}

#ifndef HAVE_LIBREADLINE
char* 
InteractiveConsole::readline(const char *prompt) {
    std::cout << prompt;
    char line[1024];
    std::cin.getline(line, 1024);
    if (std::cin.gcount() == 0) {
        // NO characters were read.
        return NULL;
    }
    // Make a copy of the line for returning
    const int len = (int) strlen(line) + 1;
    char *retVal = (char *) malloc(len);
    strcpy(retVal, line);
    return retVal;
}
#endif

#endif

