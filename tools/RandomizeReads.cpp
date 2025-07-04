#ifndef RANDOMIZE_READS_CPP
#define RANDOMIZE_READS_CPP

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
// Authors:   Dhananjai M. Rao          raodm@miamiOH.edu
//
//---------------------------------------------------------------------

#include "RandomizeReads.h"
#include "ArgParser.h"
#include "EST.h"
#include "ESTList.h"
#include "Common.h"

int
RandomizeReads::parseArgs(int argc, char *argv[], std::string& inputPath,
                          std::string& outputPath, bool& filterNs,
                          int& subSeqStart, int& subSeqEnd) {
    bool showOptions = false;
    subSeqStart = subSeqEnd = -1;  // Initialize to default value.
    // Create the list of valid arguments to be used by the arg_parser.
    ArgParser::ArgRecord arg_list[] = {
        {"--input", "The input FASTA file to be randomly shuffled",
         &inputPath, ArgParser::STRING},
        {"--output", "Output FASTA file to write randomly shuffled reads",
         &outputPath, ArgParser::STRING},
        {"--options", "Lists options for this tool",
         &showOptions, ArgParser::BOOLEAN},
        {"--filter-n", "Filter out reads that have N bases in them",
         &filterNs, ArgParser::BOOLEAN},
        {"--sub-seq-beg", "Beginning index position for subsequence",
         &subSeqStart, ArgParser::INTEGER},
        {"--sub-seq-end", "Ending index position for subsequence",
         &subSeqEnd, ArgParser::INTEGER},        
        {"", "", NULL, ArgParser::INVALID}
    };
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    ArgParser ap;
    Tool::addCmdLineArgs("randomize", ap);
    ap.addValidArguments(arg_list);
    ap.parseArguments(argc, argv, false);
    if (showOptions) {
        std::cout << ap << std::endl;
        // Nothing further to be done.
        return 1;
    }

    // Ensure we have all the necessary parameters.
    std::string tool = "RandomizeReads";
    CHECK_ARGS(tool, inputPath.empty(), "Input file not specified.\n" \
               "Use --input option\n");
    CHECK_ARGS(tool, outputPath.empty(), "Output file not specified.\n" \
               "Use --output option\n");
    // Everything looks good
    return 0;
}

int
RandomizeReads::main(int argc, char *argv[]) {
    std::string inputFile, outputFile;  // The input & output files
    bool filterNs = false;  // Filter out reads with 'N's
    int subSeqSt = -1, subSeqEnd = -1;  // Subsequence option.
    // Get helper method to process command-line args
    if (parseArgs(argc, argv, inputFile, outputFile, filterNs,
                  subSeqSt, subSeqEnd) != 0) {
        return 1;  // something is missing.
    }
    // Load reads from the supplied input file.
    RandomizeReads rnd;
    if (!rnd.loadFastaFile(inputFile)) {
        return 1;  // error loading reads.
    }
    // Extract the reads into a vector for us to randomly shuffle    
    ESTList& estList   = rnd.getESTList();
    std::vector<EST*> estVector;
    estVector.reserve(estList.size());  // Reduce 
    for (int i = 0; (i < estList.size()); i++) {
        if (filterNs && estList[i]->getSequenceString().find('N')
            != std::string::npos) {
            std::cout << estList[i]->getInfo() << ": Found N at: "
                      << estList[i]->getSequenceString().find('N') << '\n';
            continue;  // This read has 'N'. Ignore it.
        }
        estVector.push_back(estList[i]);
    }
    // Randomly shuffle the reads
    std::random_shuffle(estVector.begin(), estVector.end());
    // Write the shuffled reads out.
    std::ofstream out(outputFile);
    if (!out.good()) {
        std::cerr << "Error writing to file: " << outputFile << std::endl;
        return 2;
    }
    // Write the shuffled reads out.
    for (EST* est : estVector) {
        if ((subSeqSt != -1) && (subSeqSt < subSeqEnd) &&
            (est->getSequenceLength() > subSeqEnd)) {
            // Only dump a subsequence of the original strain.
            const char* seq = est->getSequence();
            const std::string subseq(seq + subSeqSt, seq + subSeqEnd);
            EST subseqEST(est->getID(), est->getInfo(), subseq);
            subseqEST.dumpEST(out, 70);
        } else {
            // By default we print the full sequence.
            est->dumpEST(out, 70);
        }
        out << std::endl;
    }
    // Done.
    return 0;
}

#endif
