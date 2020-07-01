#ifndef EXTRACT_READS_CPP
#define EXTRACT_READS_CPP

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

#include <numeric>
#include <random>
#include "ExtractReads.h"
#include "EST.h"
#include "ESTList.h"
#include "Common.h"

int
ExtractReads::parseArgs(int argc, char *argv[], std::string& inputPath,
                        std::string& outputPath,
                        ArgParser::StringList& idxList, int& rndNum) {
    bool showOptions = false;
    // Create the list of valid arguments to be used by the arg_parser.
    ArgParser::ArgRecord arg_list[] = {
        {"--input", "The input FASTA file to be randomly shuffled",
         &inputPath, ArgParser::STRING},
        {"--output", "Output FASTA file to write randomly shuffled reads",
         &outputPath, ArgParser::STRING},
        {"--indexs", "List index values (zero-based) for reads to be extracted",
         &idxList, ArgParser::STRING_LIST},        
        {"--options", "Lists options for this tool",
         &showOptions, ArgParser::BOOLEAN},
        {"--random", "Randomly select specified number of reads to extract",
         &rndNum, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    ArgParser ap;
    Tool::addCmdLineArgs("ExtractReads", ap);
    ap.addValidArguments(arg_list);
    ap.parseArguments(argc, argv, false);
    if (showOptions) {
        std::cout << ap << std::endl;
        // Nothing further to be done.
        return 1;
    }

    // Ensure we have all the necessary parameters.
    std::string tool = "ExtractReads";
    CHECK_ARGS(tool, inputPath.empty(), "Input file not specified.\n" \
               "Use --input option\n");
    CHECK_ARGS(tool, outputPath.empty(), "Output file not specified.\n" \
               "Use --output option\n");
    CHECK_ARGS(tool, idxList.empty() && rndNum == -1,
               "Missing List of reads to extract.\n"    \
               "Use --indexs option\n");    
    // Everything looks good
    return 0;
}

int
ExtractReads::main(int argc, char *argv[]) {
    std::string inputFile, outputFile;  // The input & output files
    ArgParser::StringList idxList;      // Indexs of reads to be extracted
    int rndNums = -1;                   // Optional random subset of reads
    // Get helper method to process command-line args
    if (parseArgs(argc, argv, inputFile, outputFile, idxList, rndNums) != 0) {
        return 1;  // something is missing.
    }
    // Load reads from the supplied input file.
    ExtractReads extractor;
    if (!extractor.loadFastaFile(inputFile)) {
        return 1;  // error loading reads.
    }

    // Write the extracted reads out.
    std::ofstream out(outputFile);
    if (!out.good()) {
        std::cerr << "Error writing to file: " << outputFile << std::endl;
        return 2;
    }

    // Obtain reference to the input reads.
    ESTList& inputList = extractor.getESTList();

    // Check and setup the list of reads if a random subset of reads
    // was desired.
    if (idxList.empty() && rndNums > 0) {
        // Generate a random set of indexs
        std::vector<int> indexs(inputList.size());   // create 'n' zeros
        std::iota(indexs.begin(), indexs.end(), 0);  // fill-in 0,1,...n
        std::random_device rd;
        std::shuffle(indexs.begin(), indexs.end(),
                     std::default_random_engine(rd()));
        
        // Write subset of indexs to output.
        for (int i = 0; (i < rndNums); i++) {
            inputList.get(indexs[i])->dumpEST(out, 70);
            out << std::endl;
        }
        
    } else {
        // Write the given list of indexs to output.
        for (size_t i = 0; (i < idxList.size()); i++) {
            // Get the desired entry from the full-list.
            const int idx  = std::stoi(idxList[i]);
            // Write the entry out.
            inputList.get(idx)->dumpEST(out, 70);
            out << std::endl;
        }
    }
    // Done.
    return 0;
}

#endif
