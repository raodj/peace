#ifndef EST_ANALYZER_CPP
#define EST_ANALYZER_CPP

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

#include "ESTAnalyzer.h"
#include "EST.h"
#include "MPIHelper.h"

// Get rid of magic numbers
#define NO_ERROR 0

// The static instance variables for command line arguments.
bool  ESTAnalyzer::readAhead      = false;
char* ESTAnalyzer::estFileName    = NULL;
bool  ESTAnalyzer::htmlLog        = false;

// The common set of arguments for all EST analyzers
arg_parser::arg_record ESTAnalyzer::commonArgsList[] = {
    {"--readAhead", "Use a read head thread to load next EST data (NYI)",
     &ESTAnalyzer::readAhead, arg_parser::BOOLEAN},
    {"--estFile", "Name of EST file (in FASTA format) to be processed",
     &ESTAnalyzer::estFileName, arg_parser::STRING},
    {"--html", "Generate analysis report in HTML format",
     &ESTAnalyzer::htmlLog, arg_parser::BOOLEAN},
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

ESTAnalyzer::ESTAnalyzer(const std::string& name, const int estIdx,
                         const std::string& outputFile) 
    : refESTidx(estIdx), chain(NULL), outputFileName(outputFile),
      analyzerName(name) {
    // Nothing else to be done for now.
      }

ESTAnalyzer::~ESTAnalyzer() {
    // Empty constructor begets an empty destructor
}

int
ESTAnalyzer::setHeuristicChain(HeuristicChain* hChain) {
    chain = hChain;
    return 0; // Everything went well
}

float
ESTAnalyzer::analyze(const int otherEST, const bool useHeuristics,
                     const bool useHeavyWeight) {
    // Check with the heuristic chain
    if (useHeuristics && chain != NULL && !chain->shouldAnalyze(otherEST)) {
        // Heuristics indicate we should not do D2. So skip it.
        return getInvalidMetric();
    }
    if (useHeavyWeight) {
        return getMetric(otherEST);
    }
    // Indicate that we would have used heavy weight analysis
    return getValidMetric();
}

void
ESTAnalyzer::showArguments(std::ostream& os) {
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(ESTAnalyzer::commonArgsList);
    os << "Common options for all EST analyzers are:\n";
    os << ap;
}

bool
ESTAnalyzer::parseArguments(int& argc, char **argv) {
    arg_parser ap(ESTAnalyzer::commonArgsList);
    // Process the arguments.  There could be some leftover for other
    // classes to process and use.
    ap.check_args(argc, argv, false);

    // Check if necessary arguments have been specified for processing.
    if (estFileName == NULL) {
        // Necessary argumetns have not been specified.
        std::cerr << analyzerName
                  << ": EST file not specified (use --estFile option)\n";
        return false;
    }
    // Everything went well.
    return true;
}

bool
ESTAnalyzer::loadFASTAFile(const char *fileName, const bool unpopulate) {
    static const std::string IgnoreFileName = "<none>";
    if (IgnoreFileName == fileName) {
        // The user does not want to use a file name. This possibly
        // happens when using peace libary from custom applications.
        return true;
    }
    FILE *fastaFile = NULL;
#ifndef _WINDOWS
    fastaFile = fopen(fileName, "rt");
#else
    fopen_s(&fastaFile, fileName, "rt");
#endif
    if ((fastaFile == NULL) || (ferror(fastaFile))) {
        std::cerr << analyzerName << "(Rank: ";
        std::cerr << MPI_GET_RANK()
                  << "): Error opening FASTA file "
                  << fileName << " for reading." << std::endl;
        return false;
    }
    // Track line number in file being processed.
    int lineNum = 1;
    // Track number of ESTs filtered out of the file
    int filteredCount = 0;
    // Repeatedly read EST's from the file.
    while (!feof(fastaFile)) {
        EST *est = EST::create(fastaFile, lineNum);
        if ((est == NULL) && (!feof(fastaFile) || ferror(fastaFile))) {
            // An error occured when reading EST.
            fclose(fastaFile);
            std::cerr << analyzerName << ": Error loading EST from "
                      << fileName << " at line: " << lineNum << std::endl;
            return false;
        }
        if (est->getID() == -1) {
            // A place holder entry to represent an EST filtered out
            // of the dataset
            filteredCount++;
        }
        if (unpopulate) {
            // For now don't hold onto the sequences in memory to
            // reduce memory usage.  These sequences would have to be
            // loaded on-demand at a later time.
            est->unpopulate();
        }
    }
    // All the EST's were succesfully read.
    if (filteredCount > 0) { // Report on ESTs filtered out, if any
        std::cerr << analyzerName << ": " << filteredCount << " sequences "
                  << "with length less than 50 nt were filtered out of "
                  << "the data set." << std::endl;
    }
    fclose(fastaFile);
    return true;
}

ESTAnalyzer& 
ESTAnalyzer::operator=(const ESTAnalyzer&) {
    return *this;
}

#endif
