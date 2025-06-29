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
                        ArgParser::StringList& idxList, int& rndNum,
                        ArgParser::StringList& substrList) {
    bool showOptions = false;
    subSeqStart = subSeqEnd = -1;
    // Create the list of valid arguments to be used by the arg_parser.
    ArgParser::ArgRecord arg_list[] = {
        {"--input", "The input FASTA file to be randomly shuffled",
         &inputPath, ArgParser::STRING},
        {"--output", "Output FASTA file to write randomly shuffled reads",
         &outputPath, ArgParser::STRING},
        {"--indexs", "List index values (zero-based) for reads to be extracted",
         &idxList, ArgParser::STRING_LIST},
        {"--substr", "List of sub-strings for which reads are to be extracted", 
         &substrList, ArgParser::STRING_LIST},
        {"--options", "Lists options for this tool",
         &showOptions, ArgParser::BOOLEAN},
        {"--random", "Randomly select specified number of reads to extract",
         &rndNum, ArgParser::INTEGER},
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
    Tool::addCmdLineArgs("extract", ap);
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
    CHECK_ARGS(tool, idxList.empty() && rndNum == -1 && substrList.empty(),
               "Missing List of reads to extract.\n"    \
               "Use --indexs or --subsstr or --random option\n");    
    // Everything looks good
    return 0;
}

// Write whole or partial fragments of reads
void
ExtractReads::writeStrain(EST& est, std::ostream& out) const {
    if ((subSeqStart != -1) && (subSeqStart < subSeqEnd) &&
        (est.getSequenceLength() > subSeqEnd)) {
        // Only dump a subsequence of the original strain.
        const char* seq = est.getSequence();
        const std::string subseq(seq + subSeqStart, seq + subSeqEnd);
        EST subseqEST(est.getID(), est.getInfo(), subseq);
        subseqEST.dumpEST(out, 70);
    } else {
        // By default we print the full sequence.
        est.dumpEST(out, 70);
    }
}

// Extract a given number of random reads.
void
ExtractReads::extractRandomReads(const int rndNums,
                                 ESTList& inputList,
                                 std::ostream& out) const {
    // Generate a random set of indexs
    std::vector<int> indexs(inputList.size());   // create 'n' zeros
    std::iota(indexs.begin(), indexs.end(), 0);  // fill-in 0,1,...n
    std::random_device rd;
    std::shuffle(indexs.begin(), indexs.end(),
                 std::default_random_engine(rd()));
    // Write subset of indexs to output.
    for (int i = 0; (i < rndNums); i++) {
        writeStrain(*inputList.get(indexs[i]), out);
    }
}

// Extract a reads given their index positions in idxList
void
ExtractReads::extractIndex(const ArgParser::StringList& idxList,
                           ESTList& inputList,
                           std::ostream& out) const {
    // Write the given list of indexs to output.
    for (size_t i = 0; (i < idxList.size()); i++) {
        // Get the desired entry from the full-list.
        const int idx  = std::stoi(idxList[i]);
        // Write the entry out.
        writeStrain(*inputList.get(idx), out);
    }    
}

// Extract reads that contain any one of the given list of sub-strings
void
ExtractReads::extractSubstr(const ArgParser::StringList& substrList,
                            ESTList& inputList,
                            std::ostream& out) const {
    // Write strains with any of the given substrings to output.
    for (int i = 0; (i < inputList.size()); i++) {
        // Get the information string associated with the EST.
        const std::string& info = inputList.get(i)->getInfo();
        for (const auto& sstr : substrList) {
            if (info.find(sstr) != std::string::npos) {
                // This read has a matching sub-string.
                // Write the entry out.
                writeStrain(*inputList.get(i), out);                
                break;
            }
        }
    }
}

int
ExtractReads::main(int argc, char *argv[]) {
    std::string inputFile, outputFile;  // The input & output files
    ArgParser::StringList idxList;      // Indexs of reads to be extracted
    ArgParser::StringList substrList;   // List of substrings
    int rndNums = -1;                   // Optional random subset of reads
    // Get helper method to process command-line args
    ExtractReads extractor;    
    if (extractor.parseArgs(argc, argv, inputFile, outputFile, idxList,
                            rndNums, substrList) != 0) {
        return 1;  // something is missing.
    }
    // Load reads from the supplied input file.
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
        // Extract given number of random reads.
        extractor.extractRandomReads(rndNums, inputList, out);
    } else if (!idxList.empty()) {
        // Extract reads at given index positions.
        extractor.extractIndex(idxList, inputList, out);
    } else if (!substrList.empty()) {
        // Extract reads based on matching sub-strings
        extractor.extractSubstr(substrList, inputList, out);
    }
    // Done.
    return 0;
}

#endif
