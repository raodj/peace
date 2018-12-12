#ifndef ACE_WRITER_CPP
#define ACE_WRITER_CPP

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

#include "ACEWriter.h"
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
#include <iterator>

// A named constant for MPI-message tag
#define ACE_AF_ENTRIES_TAG 310
#define ACE_RD_ENTRY_TAG   311

ACEWriter::ACEWriter() : Component("ACEWriter") {
    contigCount = 0;
    readCount   = 0;
}


ACEWriter::~ACEWriter() {
    // Nothing to be done for now.
}

void
ACEWriter::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"--output-ace-file", "ACE file to which assembled contig must "
         "be written",
         &aceFileName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add arguments for this component to the parser
    argParser.addValidArguments(Arguments);
}

bool
ACEWriter::writeACEHeader() {
    if (!aceFile.good()) {
        // The stream is already not in a good state. Don't do any
        // more I/O operations on it.
        return false;
    }

    // Seek to the beginning of the file.
    if (!aceFile.seekg(0).good()) {
        // Seek failed. Bail out
        return false;
    }
    // Write the ACE top-level header.
    aceFile << "AS ";
    // Save and restore current width settings
    const int oldWidth = (int) aceFile.width();
    aceFile.width(8); // Set integer to be 8-characters in size
    aceFile << contigCount << " ";
    aceFile.width(8);
    aceFile << readCount << std::endl;
    aceFile.width(oldWidth);  // restore old width.
    // Check and write error if header writing failed.
    if (!aceFile.good()) {
        std::cerr << "Error occured when generating/updating the AS header "
                  << "in the ACE file.\n";
    }
    // Return true if the aceFile is still in good shape
    return aceFile.good();
}

bool
ACEWriter::initialize() {
    // Check to ensure we have a valid ace file specified before we
    // bother to do any operation
    if (aceFileName.empty()) {
        // No ace output specified. nothing further to do
        return true;
    }
    
    // First setup reference to the shared cDNA list available in the
    // runtime context.
    RuntimeContext* rc = getContext();
    ASSERT ( rc != NULL );
    // Save ESTList for use later on in other methods
    estList = rc->getESTList();
    ASSERT ( estList != NULL );
    // Register this ACE writer as a contig listener with the
    // assembler.
    if (rc->getAssembler() != NULL) {
        rc->getAssembler()->addContigListener(this);
    } else {
        // There isn't an assembler but the user want's a ACE file
        // output. This is confusing.
        std::cerr << "Warning: A valid genome-assembler has not been "
                  << "specified but a ACE file output is being required.\n";
    }
    // Rest of the operations are performed only by the manager process
    if (MPI_GET_RANK() != MANAGER_RANK) {
        // Nothing further to do for a worker process
        return true;
    }
    // Since a valid ace file has been specified. Try and open it for
    // writing.
    aceFile.open(aceFileName.c_str(), std::fstream::out |
                 std::fstream::in | std::fstream::trunc);
    if (!aceFile.good()) {
        std::cerr << "Error opening ACE file " << aceFileName
                  << " for writing contig data.\n";
        return false;
    }
    // Dump out the top-level ACE header. Note that at this time we do
    // not know the number of contigs or the number of reads
    // constituting the ACE file. Consequently, we print zeros here
    // (to hold sufficient space). Later on (once all the contigs have
    // been written to this file) in the finzalize() method these two
    // values are updated with their actual final values.
    return writeACEHeader();
}

bool
ACEWriter::write(RuntimeContext* context) {
    ASSERT ( context != NULL );
    if (aceFileName.empty() || !aceFile.is_open()) {
        // A ace file has not been requested or we are worker process
        // and do not need further operations.
        return true;
    }
    // Check to ensure all I/O with files has proceeded fine thus far
    if (!aceFile.good()) {
        // Error occured when writing data.
        return false;
    }
    // This method is invoked after all the contigs have been
    // written. Here we need to generate the final whole-assembly (WA)
    // entry and fix-up the counts at the beginning of the file.
    aceFile << "WA{\n";
    // Create the first line in WA that must be in the form:
    // <tag type> <what program created tag> <date tag was created in
    // form YYMMDD:HHMMSS>. First build up our timestamp.
    char timestamp[128];
    time_t epoch; // milliseconds since epoch
    struct tm lt; // local time
    time(&epoch); // get milliseconds since epoch
    localtime_r(&epoch, &lt);
    sprintf(timestamp, "%02d%02d%02d:%02d%02d%02d", lt.tm_year - 100,
            lt.tm_mon + 1, lt.tm_mday + 1, lt.tm_hour, lt.tm_min, lt.tm_sec);
            
    aceFile << "PEACE_tag PEACE " << timestamp << std::endl;
    // Dump out the PEACE information
    aceFile << "PEACE version "  << PEACE_VERSION
            << " (Released on: " << RELEASE_DATE << ")\n";
    // Dump out the command line for this invocation
    aceFile << context->getCommandLine() << std::endl;
    aceFile << "}";
    // Check and ensure the ACE file is good.
    if (!aceFile.good()) {
        // Some error occured during I/O operations.
        std::cerr << "Error occurred when generating WA data "
                  << "in ACE file.\n";
        return false;
    }
    // Fix-up the counts at the beginning of the ACE file.
    if (!writeACEHeader()) {
        return false;
    }
    // Close the ace file as we are done with it.
    aceFile.close();
    // Return status depending on status of out file
    return aceFile.good();
}

bool
ACEWriter::contigFormed(const Assembler& UNUSED(assembler),
                        const Contig& contig, const bool UNUSED(fullContig)) {
    // Convenience reference to the list of local-cDNA fragments in the contig
    const AlignmentInfoList& aiList = contig.getAlignments();    
    // This vector holds the number of alignment entries
    // (aiList.size()) at each parallel process for the given contig
    std::vector<int> aiCounts;
    const int MyRank   = getGlobalContigCounts(aiList.size(), aiCounts);
    // Compute total number of reads in the contig.
    const int numReads = std::accumulate(aiCounts.begin(), aiCounts.end(), 0);
    
    // Only on manager process, generate the contig header in the form:
    // CO <contig name> <# of bases> <# of reads in contig> <# of base
    // segments in contig> <U or C>
    if (MyRank == MANAGER_RANK) {
        const std::string& consensus = contig.getConsensus();
        aceFile << "CO " << contig.getId()
                << " "   << consensus.size()
                << " "   << numReads // # of reads in contig.
                << " "   << contig.getBaseSegments().size()
                << " U"  << std::endl;
        // Dump out the consensus sequence.
        aceFile << consensus << std::endl;
        // Dump out the base qualities for the consensus sequence if
        // available.
        const std::vector<int>& qualValues = contig.getQuality();
        if (!qualValues.empty()) {
            aceFile << "BQ\n";
            std::ostream_iterator<int> qualWriter(aceFile, " ");
            std::copy(qualValues.begin(), qualValues.end(), qualWriter);
        }
    }

    // Next generate the AF entries via a helper method. In the
    // following method call, worker processes build their local
    // information and send it to the manager. The manager writes its
    // own local information and then gets data from workers and write
    // it to the file.
    if (!writeAFEntries(aiList, aiCounts, (MyRank == MANAGER_RANK))) {
        // Error when writing AF entries.
        return false;
    }
    // Next write the base segment (BS) entries. These entries are
    // locally available on the manager process.
    if (MyRank == MANAGER_RANK) {
        writeBSEntries(contig.getBaseSegments());
    }
    // Lastly, we generate the RD entries via a helper method for this
    // contig. In the following method call, worker processes build
    // their local information and send it to the manager. The manager
    // writes its own local information and then gets data from
    // workers and write it to the ACE file.
    if (!writeRDEntries(aiList, numReads, (MyRank == MANAGER_RANK))) {
        // Error when writing RD entries.
        return false;
    }
    // Update our counters for fixing-up ACE header
    contigCount++;
    readCount += numReads;
    // Return overall result
    return aceFile.good();
}

bool
ACEWriter::writeAFEntries(const AlignmentInfoList& aiList,
                          const std::vector<int>& aiCounts,
                          const bool isManager) {
    // Create a large string with all the AF entries to ease sending
    // data over the wire and/or writing to file.
    {  // We intentionally open a dummy scope to keep the life of
       // afEntries buffer to a minimum.
        std::ostringstream afEntries;
        // Generate entries for each alignment in aiList
        for(size_t i = 0; (i < aiList.size()); i++) {
            // Get EST (and load EST information as needed).
            const EST* est = estList->get(aiList[i].getESTIndex(), true);
            afEntries << "AF " << est->getInfo()
                      << " "   << est->getSequenceLength()
                      << " "   << (aiList[i].getPosition() + 1)
                      << std::endl;
        }
        // On worker processes we just ship afEntries to the manager and
        // we are done.
        if (!isManager && !aiList.empty()) {
            const std::string &afInfo = afEntries.str();
            MPI_SEND(&afInfo[0], afInfo.size(), MPI_TYPE_CHAR,
                     MANAGER_RANK, ACE_AF_ENTRIES_TAG);
            return true;
        } else {
            // When control drops here we are the manager process.
            aceFile << afEntries.str();
        }
    }
    // When control drops here we are the manager process.  Get AF
    // entries from all the workers and write it to disk.
    std::string recvBuffer;    
    for(size_t workerRank = 1; (workerRank < aiCounts.size()); workerRank++) {
        if (aiCounts[workerRank] <= 0) {
            // This worker has no entries for us.
            continue;
        }
        // Figure out how much data the worker is going to send.
        MPI_STATUS msgInfo;        
        MPI_PROBE(workerRank, ACE_AF_ENTRIES_TAG, msgInfo);
        // Reserve the necessary memory buffer
        const int msgSize = MPI_GET_COUNT(msgInfo, MPI_TYPE_CHAR);
        recvBuffer.resize(msgSize);
        // Now read the actual message
        MPI_RECV(&recvBuffer[0], msgSize, MPI_TYPE_CHAR, workerRank,
                 ACE_AF_ENTRIES_TAG);
        aceFile << recvBuffer;
    }
    if (!aceFile.good()) {
        std::cerr << "Error writing AF entries to the ACE file.\n";
    }
    // Return true/false depending on file goodness
    return aceFile.good();
}

bool
ACEWriter::writeBSEntries(const BaseSegmentInfoList& bsList) {
    for(size_t i = 0; (i < bsList.size()); i++) {
        const BaseSegmentInfo &bs = bsList[i];
        aceFile << "BS " << (bs.getStartPos() + 1)
                << " "   << (bs.getEndPos() + 1)
                << " "   << estList->getESTInfo(bs.getESTIndex())
                << "\n";
    }
    // Return true/false depending on file goodness
    return aceFile.good();    
}

bool
ACEWriter::writeRDEntries(const AlignmentInfoList& aiList,
                          const int readCount, const bool isManager) {
    // First handle local information in the following manner:
    // 1. Create the RD entry for each read in the contig.
    // 2. Do one of the following:
    //    2.1. On manager process write it to ACE file.
    //    2.2. On worker  process send  it to manager.
    for(size_t i = 0; (i < aiList.size()); i++) {
        const EST* est = estList->get(aiList[i].getESTIndex(), true);
        if (isManager) {
            writeRDAlignedRead(aceFile, est->getInfo(), est->getSequence(),
                               aiList[i].getCigar());
        } else {
            std::ostringstream buffer;
            writeRDAlignedRead(buffer, est->getInfo(), est->getSequence(),
                               aiList[i].getCigar());
            // Send entry to the manager process.
            const std::string& entry = buffer.str();
            MPI_SEND(entry.c_str(), entry.size(), MPI_TYPE_CHAR,
                     MANAGER_RANK, ACE_RD_ENTRY_TAG);
        }
    }
    if (!isManager) {
        // Workers are done.
        return true;
    }
    // Only on the manager we receive all the entries from the various
    // workers and write that information to the ACE file.
    const int WorkerReadCount = readCount - aiList.size();
    std::string recvBuffer;    
    for(int i = 0; (i < WorkerReadCount); i++) {
        // Figure out how much data the worker is going to send.
        MPI_STATUS msgInfo;        
        MPI_PROBE(MPI_ANY_SOURCE, ACE_RD_ENTRY_TAG, msgInfo);
        // Reserve the necessary memory buffer
        const int msgSize = MPI_GET_COUNT(msgInfo, MPI_TYPE_CHAR);
        recvBuffer.resize(msgSize);
        // Now read the actual message
        MPI_RECV(&recvBuffer[0], msgSize, MPI_TYPE_CHAR, msgInfo.MPI_SOURCE,
                 ACE_RD_ENTRY_TAG);
        aceFile << recvBuffer;
    }
    if (!aceFile.good()) {
        std::cerr << "Error writing RD entries to the ACE file.\n";
    }
    // Return true/false depending on file goodness
    return aceFile.good();        
}

bool
ACEWriter::writeRDAlignedRead(std::ostream& os, const std::string& info,
                              const std::string& seq,
                              const std::string& cigarStr) const {
    // First create the RD entry
    os << "RD " << info
       << " "   << seq.size()
       << " 0 0\n";
    // Initialize the boudaries of clippings assuming there is no
    // clipping.
    size_t qualClipStart = 0, qualClipEnd = seq.size();
    size_t asmClipStart  = 0, asmClipEnd  = qualClipEnd;
    size_t ntPos         = 0;    // Current nucleotide pos in seq
    bool   firstClips    = true; // flag to indicate if clipping is first one
    
    // To ease processing of the CIGAR string we use an istringstream.
    std::istringstream cigar(cigarStr);
    while (!cigar.eof() && (ntPos < seq.size())) {
        // Read the next entry from the cigar stream and process it
        // appropriately.
        int  len = 0;
        char op  = '*';
        cigar >> len >> op;
        ASSERT ( len > 0 );
        // Do various operations based on op.
        switch (op) {
        case 'H': // hard-clipping (can only be at the beginning or end)
            len = (firstClips ? (qualClipStart = len) : (qualClipEnd = len));
            break;
        case 'S': // soft-clipping
            len = (firstClips ? (asmClipStart = len) : (asmClipEnd = len));
            break;
        case 'M': // match with alignment.
        case '=': // Exact sequence match.
        case 'X': // sequence mismatch            
            os << seq.substr(ntPos, len);
            break;
        case'I': // This sequence has an insert wrt to reference.
            break;
        case 'D': // This is deletion or mismatch region
        case 'N': // skipped region from reference
            os << std::string(len, '*');
            len = 0; // Do not skip over nt.
            break;
        case 'P': // padding (slient deletion from padded reference)
        default:
            std::cerr << "Warning '" << op
                      << "' entries in CIGAR (" << cigarStr
                      << ") for read " << info << " is not yet handled\n";
        }
        // Skip over nucleotides that were handled by the operation
        ntPos += len;
        // Clear firstClip flag if we encountered a non-clipping entry
        firstClips = ((op != 'H') && (op != 'S'));
    }
    // Write the final QA entry as well.
    os << "\nQA " << (qualClipStart + 1) << " " << (qualClipEnd + 1)
       << " "     << (asmClipStart  + 1) << " " << (asmClipEnd  + 1)
       << std::endl;
    // Everything went well
    return true;
}

int
ACEWriter::getGlobalContigCounts(const int localCount,
                                 std::vector<int>& globalEntryCounts) {
    // First determine the number of alignment information for the
    // contig on the various processes to determine which processes
    // have data and which process has the first and last entry to set
    // necessary flags in ACE file.
    const int MyRank       = MPI_GET_RANK();
    const int ProcCount    = MPI_GET_SIZE();
    std::vector<int> localEntryCount (ProcCount, 0);
    globalEntryCounts.resize(ProcCount, 0);
    localEntryCount[MyRank] = localCount;
    MPI_ALL_REDUCE(&localEntryCount[0], &globalEntryCounts[0], ProcCount,
                   MPI_TYPE_INT, MPI_OP_SUM);
    // Everything went well.
    return MyRank;
}

#endif
