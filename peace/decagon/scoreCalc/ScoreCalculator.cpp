#ifndef SCORE_CALCULATOR_CPP
#define SCORE_CALCULATOR_CPP

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

#include "ScoreCalculator.h"
#include "Version.h"
#include "ArgParser.h"
#include "FASTAFile.h"
#include "ACEReader.h"
#include "BlastDBMaker.h"
#include "AlignmentAlgorithm.h"

#include <cctype>
#include <algorithm>

// Simple string that is used to display the title
#define DECAGON_TITLE "Distributed Environment for Comparative Analysis of " \
    "Genomic Assemblers (DECAGON)"

// Simple string that is composed from various related information to
// display version details to the user when requested.
#define DECAGON_VERSION_MESSAGE                                   \
    DECAGON_TITLE                                                 \
    "\n"                                                          \
    PEACE_VERSION                                                 \
    ". Released (under GPLv3) on "                                \
    RELEASE_DATE                                                  \
    "\n"                                                          \
    COPYRIGHT

// Another simple string message that is built to display help
// information to the user about how to use the analysis tool
#define DECAGON_BRIEF_MESSAGE                                     \
    DECAGON_VERSION_MESSAGE                                       \
    "\nUsage: scoreCalc [arguments]\n"                            \
    "where valid arguments are:"

ScoreCalculator::ScoreCalculator() {
    // Nothing else to be done for now
}

ScoreCalculator::~ScoreCalculator() {
    // Nothing else to be done for now
}

bool
ScoreCalculator::run(int& argc, char *argv[]) throw(DecagonException) {
    if (!parseCmdLineArgs(argc, argv)) {
        return false;
    }
    // Verify that the supplied files are FASTA file and load the
    // source sequences for further processing.
    loadReads(srcTransFile,  srcTransList);
    // Next load contigs from the assembler-generated contig file.
    loadContigs(genContigFile, genContigList);
    // Check and create the BLAST DB as needed.
    BlastDBMaker dbMaker;
    const std::string blastDBName = dbMaker.checkMakeDB(srcTransFile);
    // Next create a temporary FASTA file with contigs to run Megablast
    const std::string genContigFastaFile = (genContigFile + ".fa");
    writeContigsToFASTAFile(genContigList, genContigFastaFile);
    // Next run Mega Blast to get best matching contig->transcript pairs
    MegaBlastRunner mbr;
    MBDSet goodPairs;
    mbr.runMegaBlast(blastDBName, genContigFastaFile, goodPairs);
    
    // Generate outputs for each read. The outputs are written to a
    // file or to std::cout
    std::ostream outStream(outputFile.good() ? outputFile.rdbuf() :
                           std::cout.rdbuf());
    generateScores(goodPairs, outStream);
    // Success
    return true;
}

void
ScoreCalculator::writeContigsToFASTAFile(const ContigList& contigList,
                                         const std::string& fastaFileName)
    throw (DecagonException) {
    // Try and open a FASTA file for writing.
    std::ofstream fastaFile(fastaFileName.c_str());
    if (!fastaFile.good()) {
        std::string errMsg = "Error attempting to create intermediate "
            "FASTA file: ";
        errMsg += fastaFileName;
        throw EXP(4, errMsg.c_str(), "Ensure file name and path are valid.");
    }
    // Write contigs to the fasta file.
    for(size_t i = 0; (i < contigList.size()); i++) {
        const Contig& contig = contigList[i];
        fastaFile << ">" << contig.getId() << std::endl
                  << contig.getConsensus() << std::endl;
        // Check to ensure that the FASTA file is good after each write
        if (!fastaFile.good()) {
            std::string errMsg = "Error writing to intermediate FASTA file: ";
            errMsg += fastaFileName;
            throw EXP(4, errMsg.c_str(), "Check disk space.");
        }
    }
}

void
ScoreCalculator::generateScores(const MBDSet& goodPairs,
                                std::ostream& os) {
    AlignmentAlgorithm aa;
    for(MBDSet::const_iterator entry = goodPairs.begin();
        (entry != goodPairs.end()); entry++) {
        const Contig& genContig= findContig(entry->getGenContigID(),
                                            genContigList);
        const EST& srcTrans = findRead(entry->getSrcTransID(), srcTransList);
        // Convert the transcript and contig nucleotide sequences to
        // upper case as needed by the alignment algorithms.
        std::string srcSeq = srcTrans.getSequence();
        std::transform(srcSeq.begin(), srcSeq.end(), srcSeq.begin(), toupper);
        std::string cntSeq = genContig.getConsensus();
        std::transform(cntSeq.begin(), cntSeq.end(), cntSeq.begin(), toupper);
        // Use the alignment algorithm to compute scores
        int alignScore = -1, aScore = -1;
        std::string alignedSrcTrans, alignedGenContig;
        aa.getNWAlignment(srcSeq, cntSeq, 20,
                          alignScore, aScore, alignedSrcTrans,
                          alignedGenContig);
        os << alignedGenContig << std::endl << alignedSrcTrans << std::endl;
        // Dump out the information
        os << *entry << "\t" << alignScore << "\t" << aScore << std::endl;
    }
}

bool
ScoreCalculator::parseCmdLineArgs(int& argc, char *argv[]) {
    // Next, setup our  arguments
    bool version = false;
    std::string outputFilePath;
    const ArgParser::ArgRecord CmdLineArgs[] = {
        {"", DECAGON_BRIEF_MESSAGE, NULL, ArgParser::MAIN_MESSAGE},
        {"--gen-contig-file", "File (ACE or FASTA) with contigs generated by assembler",
         &genContigFile, ArgParser::STRING},
        {"--src-trans-file", "FASTA file with source transcripts",
         &srcTransFile, ArgParser::STRING},
        {"--output-file", "Text file to which output is written",
         &outputFilePath, ArgParser::STRING},
        {"--version", "Display just version information (and stop)",
         &version, ArgParser::BOOLEAN},        
        {"", "", NULL, ArgParser::INVALID}
    };
    // Setup the arg parser with various command-line parameters
    ArgParser argParser;
    argParser.addValidArguments(CmdLineArgs);
    // Parse command-line arguments
    argParser.parseArguments(argc, argv, true);
    // If user asked for version then display version and exit.
    if (version) {
        std::cout << DECAGON_VERSION_MESSAGE << std::endl;
        return false; // No further operations are needed. I per
    }
    // Ensure we have the necessary arguments.
    if (srcTransFile.empty() || genContigFile.empty()) {
        std::cout << "ERROR: FASTA file containing source transcript and "
                  << "the contig file\ncontaining assembler-generated contigs "
                  << "must be specified via suitable\n"
                  << "command-line parameters shown below."
                  << "\n\n" << argParser;
        return false;
    }
    // Create the output file if the user has specified it to ensure
    // it is valid and writeable.
    if (!outputFilePath.empty()) {
        outputFile.open(outputFilePath.c_str());
        if (!outputFile.good()) {
            std::cout << "ERROR: Unable to open output file '"
                      << outputFilePath << "' for writing results.\n"
                      <<  "Check path and permissions for output file\n";
        }
    }
    // When control reaches here all the command-line arguments needed
    // have been specified and their values are good.
    return true;
}

const Contig&
ScoreCalculator::findContig(const std::string& prefix,
                            const ContigList& list) const
    throw(DecagonException) {
    for(size_t entry = 0; (entry < list.size()); entry++) {
        if (list[entry].getId().find(prefix) == 0) {
            // Found a match
            return list[entry];
        }
    }
    // The entry we are looking for was not found
    std::string errMsg = "Unable to find contig for ID: ";
    errMsg += prefix;
    throw EXP(3, errMsg.c_str(),
              "Possibly the BLAST version is incompatible");
}

const EST&
ScoreCalculator::findRead(const std::string& prefix,
                          const ESTList& list) const
    throw(DecagonException) {
    for(int entry = 0; (entry < list.size()); entry++) {
        if (list.getESTInfo(entry).find(prefix) == 0) {
            // Found a match
            return *list[entry];
        }
    }
    // The entry we are looking for was not found
    std::string errMsg = "Unable to find read for ID: ";
    errMsg += prefix;
    throw EXP(3, errMsg.c_str(),
              "Possibly the BLAST version is incompatible");
}

void
ScoreCalculator::loadReads(const std::string& fileName, ESTList& targetList)
    throw(DecagonException) {
    // First check to ensure that the source file name is a FASTA file.
    if (!FASTAFile::isFASTA(fileName)) {
        std::string errorMsg = "The following file is not a FASTA file: ";
        errorMsg += fileName;
        throw EXP(1, errorMsg.c_str(),
                  "Ensure the file name, path, and format are valid.");
    }
    // Next load the FASTA file into the supplied targetList
    FASTAFile *srcFile = new FASTAFile(fileName);
    // The target list will take care of deleting the srcFile in
    // successful cases.
    if (!targetList.add(srcFile, 0)) {
        delete srcFile;
        std::string errorMsg = "Error loading data from FASTA file: ";
        errorMsg += fileName;
        // Error loading data into the file. throw an exception
        throw EXP(1, errorMsg.c_str(),
                  "Ensure that the file is a valid FASTA file.");
    }
}

void
ScoreCalculator::loadContigs(const std::string& contigFileName,
                             ContigList& contigList)
    throw(DecagonException) {
    if (FASTAFile::isFASTA(contigFileName)) {
        // Load cDNA reads from a FASTA file.
        FASTAFile fasta(contigFileName);
        int dummyContigID = 0;
        while (fasta.good() && fasta.hasNextEntry()) {
            EST* cDNA = fasta.getNextEntry(dummyContigID++);
            ASSERT ( cDNA != NULL );
            contigList.push_back(Contig(cDNA->getInfo(), cDNA->getSequence()));
            // delete cDNA;
        }
    } else {
        // When control drops here it must be an ACE file.  Create an ACE
        // reader to read data from an ACE file. The constructor for
        // ACEReader checks to ensure the file is valid and is an ACE
        // file.
        ACEReader ar(contigFileName);
        // Repeatedly read contigs into the contig list.
        do {
            // Create a new/dummy contig entry at the end of the list.
            contigList.push_back(Contig());
            // Read next contig if available. Put individual reads into a
            // dummy list.
            ESTList ignoredInfo;
            if (!ar.readContig(contigList.back(), ignoredInfo)) {
                // There was no more contig information. The last contig
                // can be undone.
                contigList.pop_back();
            }
        } while (ar.good() && ar.hasNextEntry());
    }
}

#endif
