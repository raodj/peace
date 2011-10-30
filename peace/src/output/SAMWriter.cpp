#ifndef SAM_WRITER_CPP
#define SAM_WRITER_CPP

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

#include "SAMWriter.h"
#include "Assembler.h"
#include "ArgParser.h"
#include "Utilities.h"
#include "RuntimeContext.h"
#include "MPIStats.h"
#include "ESTList.h"
#include "Version.h"

#include <fstream>
#include <numeric>
#include <sstream>
#include <cstdio>

// The static suffix added to SAM-filename to create temporary
// contig-filename
#define TMP_FILE_SUFFIX "_alignSec_tmpFile"

// A named constant for MPI-message tag
#define SAM_CDNA_INFO_TAG      300
#define SEND_SAM_CDNA_INFO_TAG 301

SAMWriter::SAMWriter() : Component("SAMWriter") {
    // Nothing else to be done for now.
}


SAMWriter::~SAMWriter() {
    // Nothing to be done for now.
}

void
SAMWriter::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"--output-sam-file", "SAM file to which assembled contig must "
         "be written",
         &samFileName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add arguments for this component to the parser
    argParser.addValidArguments(Arguments);
}

bool
SAMWriter::initialize() {
    // Check to ensure we have a valid sam file specified before we
    // bother to do any operation
    if (samFileName.empty()) {
        // No sam output specified. nothing further to do
        return true;
    }
    
    // First setup reference to the shared cDNA list available in the
    // runtime context.
    RuntimeContext* rc = getContext();
    ASSERT ( rc != NULL );
    // Save ESTList for use later on in other methods
    estList = rc->getESTList();
    ASSERT ( estList != NULL );
    // Register this SAM writer as a contig listener with the
    // assembler.
    if (rc->getAssembler() != NULL) {
        rc->getAssembler()->addContigListener(this);
    } else {
        // There isn't an assembler but the user want's a SAM file
        // output. This is confusing.
        std::cerr << "Warning: A valid genome-assembler has not been "
                  << "specified but a SAM file output is being required.\n";
    }
    // Rest of the operations are performed only by the manager process
    if (MPI_GET_RANK() != MANAGER_RANK) {
        // Nothing further to do for a worker process
        return true;
    }
    // Since a valid sam file has been specified. Try and open it for
    // writing.
    samFile.open(samFileName.c_str(), std::ios::out);
    if (!samFile.good()) {
        std::cerr << "Error opening SAM file " << samFileName
                  << " for writing contig data.\n";
        return false;
    }
    // Dump out the top-level SAM header
    samFile << "@HD\tVN:1.4\tSO:unsorted\n";
    
    // Next try and open the temporary file for writing.
    std::string tmpFileName = samFileName + TMP_FILE_SUFFIX;
    tmpContigFile.open(tmpFileName.c_str(), std::ios::out | std::ios::in |
                       std::ios::trunc);
    if (!tmpContigFile.good()) {
        std::cerr << "Error opening temporary alignment section "
                  << "file " << samFileName
                  << " for writing alignment data.\n";
        return false;
    }
    // Everything went well.
    return true;
}

bool
SAMWriter::write(RuntimeContext* context) {
    ASSERT ( context != NULL );
    if (samFileName.empty() || !samFile.is_open()) {
        // A sam file has not been requested or we are worker process
        // and do not need further operations.
        return true;
    }
    // Check to ensure all I/O with files has proceeded fine thus far
    if (!samFile.good() || !tmpContigFile.good()) {
        // Error occured when writing data.
        return false;
    }
    // This method is invoked after all the contigs have been
    // written. Here we need to move the per-cDNA alignment data from
    // tmpContigFile and append it to samFile. However, first we
    // finish up SAM headers by adding program-information headers.
    samFile << "@PG\tID:PEACE1\tPN:PEACE (http://www.peace-tools.org/)\tCL:"
            << context->getCommandLine()
            << "\tVN:" PEACE_VERSION << " (Released on: "
            << RELEASE_DATE << ")\n";
    // Copy the data from tmpContigFile to samFile.
    tmpContigFile.seekg(0); // Reset to beginning
    char buffer[8192];
    while (samFile.good() && tmpContigFile.good()) {
        tmpContigFile.read(buffer, sizeof(buffer));
        const int bytesRead = (int) tmpContigFile.gcount();
        samFile.write(buffer, bytesRead);
    }
    // Check and ensure the SAM file is good.
    if (!samFile.good() || tmpContigFile.bad()) {
        // Some error occured during I/O operations.
        std::cerr << "Error occurred when moving data from temporary contig "
                  << "file to final SAM file.\n";
        return false;
    }
    // Remove the temporary file as we no longer need it.
    tmpContigFile.close();
    std::string tmpFileName = samFileName + TMP_FILE_SUFFIX;
    if (remove(tmpFileName.c_str()) != 0) {
        std::cerr << "Error removing intermediate temporary file "
                  << tmpFileName << ".\n"
                  << "You will need to manually delete this file.\n";
    }
    // Close the sam file as we are done with it.
    samFile.close();
    // Return status depending on status of out file
    return samFile.good();
}

bool
SAMWriter::contigFormed(const Assembler& UNUSED(assembler),
                        const Contig& contig, const bool UNUSED(fullContig)) {
    // Convenience reference to the list of local-cDNA fragments in the contig
    const AlignmentInfoList& aiList = contig.getAlignments();    
    // First determine the number of alignment information for the
    // contig on the various processes to determine which processes
    // have data and which process has the first and last entry to set
    // necessary flags in SAM file.
    std::vector<int> globalEntryCounts;
    // The following hasFirst and hasLast indicates if this process
    // has the first and/or the last cDNA entry for this contig. This
    // information is needed to set bit-flags in SAM file format.
    bool hasFirst = false, hasLast = false, isSingleton = false;
    const int MyRank = getGlobalContigCounts(aiList.size(), globalEntryCounts,
                                             hasFirst, hasLast, isSingleton);
    if (MyRank != MANAGER_RANK) {
        if (aiList.empty()) {
            // We have no alignment information to send to manager
            return true;
        }
        // Wait for manager to indicate it is time to send information
        int dummy;
        TRACK_IDLE_TIME(MPI_RECV(&dummy, 1, MPI_INT, MANAGER_RANK,
                                 SEND_SAM_CDNA_INFO_TAG));
        // Send our local per-cDNA alignments to the manager.
        for(size_t i = 0; (i < aiList.size()); i++) {
            const bool isFirst = (hasFirst && (i == 0));
            const bool isLast  = (hasLast  && (i == aiList.size() - 1));
            std::string info   = getSAMAlignEntry(contig, aiList[i], isFirst,
                                                  isLast, isSingleton);
            MPI_SEND(info.c_str(), info.size(), MPI_TYPE_CHAR,
                     MANAGER_RANK, SAM_CDNA_INFO_TAG);
        }
        // Finally, receive overall success status from manager.
        int result = 0;
        TRACK_IDLE_TIME(MPI_BCAST(&result, 1, MPI_TYPE_INT, MANAGER_RANK));
        return (result == 1); // -1 indicates error, +1 indicates success!
    }
    
    // When control drops here we do operations for the manager process.
    // First, write contig header to main-SAM file, if it is not a singleton.
    if (!isSingleton) {
        samFile << "@SQ\tSN:" << contig.getId()
                << "\tLN:"    << contig.getConsensus().size() << "\n";
    }
    // Next write the local alignments to the temporary alignments
    // file (these will be added to the final SAM file at a later
    // date).
    for(size_t i = 0; (i < aiList.size()); i++) {
        // Compute if aiList[i] is the first or last entry
        const bool isFirst = (hasFirst && (i == 0));
        const bool isLast  = (hasLast  && (i == aiList.size() - 1));        
        tmpContigFile << getSAMAlignEntry(contig, aiList[i], isFirst,
                                          isLast, isSingleton) << "\n";
    }
    // Lastly, get information from various processes and write them.
    std::string recvBuffer;
    for(size_t procID = 1; (procID < globalEntryCounts.size()); procID++) {
        // Send message to process to start sending alignment information
        MPI_SEND(&procID, 1, MPI_INT, procID, SEND_SAM_CDNA_INFO_TAG);
        for(int entry = 0; (entry < globalEntryCounts[procID]); entry++) {
            // Probe and get message size to be received.
            MPI_STATUS msgInfo;
            MPI_PROBE(procID, SAM_CDNA_INFO_TAG, msgInfo);
            // Reserve the necessary memory buffer
            const int msgSize = msgInfo.Get_count(MPI_TYPE_CHAR);
            recvBuffer.resize(msgSize);
            // Now read the actual message
            MPI_RECV(&recvBuffer[0], msgSize, MPI_TYPE_CHAR, procID,
                     SAM_CDNA_INFO_TAG);
            // Write the recieved data to temporary contig file
            tmpContigFile << recvBuffer << "\n";
        }
    }
    // Report overall status to all the processes.
    int result = (tmpContigFile.good() && samFile.good() ? 1 : -1);
    MPI_BCAST(&result, 1, MPI_TYPE_INT, MANAGER_RANK);
    // Check and perform any local operations
    if (result != 1) {
        std::cerr << "Error writing contig to SAM file and/or temporary "
                  << "contig file!\n";
    }
    // return with overall status
    return (result == 1);
}

int
SAMWriter::getGlobalContigCounts(const int localCount,
                                 std::vector<int>& globalEntryCounts,
                                 bool& hasFirst, bool& hasLast,
                                 bool& isSingleton) {
    // First determine the number of alignment information for the
    // contig on the various processes to determine which processes
    // have data and which process has the first and last entry to set
    // necessary flags in SAM file.
    const int MyRank       = MPI_GET_RANK();
    const int ProcCount    = MPI_GET_SIZE();
    std::vector<int> localEntryCount (ProcCount, 0);
    globalEntryCounts.resize(ProcCount, 0);
    localEntryCount[MyRank] = localCount;
    MPI_ALL_REDUCE(&localEntryCount[0], &globalEntryCounts[0], ProcCount,
                   MPI_TYPE_INT, MPI_OP_SUM);
    // Now the vector globalEntryCount has a count of number of
    // per-cDNA alignments at each parallel process. Determine if we
    // have the first and/or the last entry by counting entries before
    // and after us.
    std::vector<int>::const_iterator start = globalEntryCounts.begin();
    const int entriesBefore = std::accumulate(start, start + MyRank, 0);
    const int entriesAfter  = std::accumulate(start + (MyRank + 1),
                                              start + ProcCount, 0);
    const int entriesTotal  = std::accumulate(start, start + ProcCount, 0);
    // Setup if we are the first and/or last entries
    hasFirst    = (entriesBefore == 0);
    hasLast     = (entriesAfter  == 0);
    isSingleton = (entriesTotal  == 1 );
    // Everything went well.
    return MyRank;
}

std::string
SAMWriter::getSAMAlignEntry(const Contig& contig, const AlignmentInfo& ai,
                            const bool isFirst, const bool isLast,
                            const bool isSingleton) {
    // Get the cDNA fragment for this entry
    ASSERT ( estList != NULL );
    EST* est = estList->get(ai.getESTIndex(), true);
    ASSERT ( est != NULL );
    // Create the SAM FLAG for this entry.
    const int flag = ((isFirst ? 0x40 : 0) + (isLast ? 0x80 : 0) +
                      (isSingleton ? 0x4 : 0));
    // Determine the SAM-RNEXT field value
    const std::string rNext = (isFirst || isLast) ? "=" : "*";
    // Convert quality from EST to string (if quality is available);
    std::string qualStr;
    if (est->hasQuality()) {
        // Convert the quality values to ASCII values compliant with
        // the quality string in Sanger FASTAQ format.
        const QualityVector& qualVec = est->getQuality();
        qualStr.reserve(qualVec.size());
        for(size_t idx = 0; (idx < qualVec.size()); idx++) {
            qualStr.push_back((char) (33 + qualVec[idx]));
        }
    } else {
        // Quality data is not available
        qualStr = "*";
    }
    // Set-up the reference name value. It is "*" for singleton contigs
    const std::string refSeqName = (isSingleton ? "*" : contig.getId());
    // Get the regular or reverse-complementary nucleotide sequence.
    const std::string seq = (ai.isRCAlignment() ? est->getRCSequence() :
                             est->getSequence());
    // Create the SAM entry.
    std::ostringstream samEntry;
    samEntry //<< est->getInfo()       << "\t"
        << "est" << ai.getESTIndex() << "\t"
             << flag                 << "\t"
             << refSeqName           << "\t"
             << ai.getPosition() + 1 << "\t"
             << ai.getMapQuality()   << "\t"
             << ai.getCigar()        << "\t"
             << rNext                << "\t"
             << 0                    << "\t"  // PNEXT
             << ai.getTemplateLen()  << "\t"
             << seq                  << "\t"
             << qualStr;
    // Return the string for the sam entry.
    return samEntry.str();
}

#endif
