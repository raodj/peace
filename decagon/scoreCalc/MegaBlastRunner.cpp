#ifndef MEGA_BLAST_RUNNER_CPP
#define MEGA_BLAST_RUNNER_CPP

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

#include "MegaBlastRunner.h"
#include <iostream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <iterator>

MegaBlastRunner::MegaBlastRunner() {
    pipe = NULL;
}

MegaBlastRunner::~MegaBlastRunner() {
    if (pipe != NULL) {
        // Only reason we need to close it here again is because there
        // could have been an error in runMegaBlast() method
        pclose(pipe);
        pipe = NULL;
    }
}

void
MegaBlastRunner::runMegaBlast(const std::string& transBlastDB,
                              const std::string& genContigFile,
                              MBDSet& results)
    throw (DecagonException) {
    // Build up the command for calling megablast.
    const std::string cmd = getMegablastCmdLine(transBlastDB, genContigFile);
    std::cout << cmd << std::endl;
    // Run megablast to read outputs from it.
    if ((pipe = popen(cmd.c_str(), "r")) == NULL) {
        std::string errMsg = "Error running megablast\n(Command used: ";
        errMsg += cmd; errMsg += ")";
        throw EXP(1, errMsg.c_str(),
                  "Verify megablast is installed and is in your path");
    }
    // Read line-by-line of output from megablast and process it.
    char line[2048];
    while (fgets(line, sizeof(line), pipe) != NULL) {
        if (line[0] == '#') {
            // Lines starting with '#' are comment lines. skip them.
            continue;
        }
        // Convert line to a MegaBlastData object.
        MegaBlastData mbd = parseLine(line);
        // Retain the best entry for the contig so far
        setBestEntry(results, mbd);
    }
    // Close our connection to Megablast.
    const int exitCode = pclose(pipe);
    pipe = NULL; // No need to try to close it again anymore.
    if (exitCode != 0) {
        std::ostringstream msg;
        msg << "Megablast exited with error code: " << exitCode;
        throw EXP(5, msg.str().c_str(),
                  "Verify BLAST DB path and Contig file path are correct");
    }
}

MegaBlastData
MegaBlastRunner::parseLine(const std::string& line) const
    throw (DecagonException) {
    char genContigID[2048], srcTransID[2048];
    double eValue = -1, bitScore = -1;
    double dummy  = 0;
    const int fieldsScanned = 
        sscanf(line.c_str(),
               "%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               genContigID, srcTransID, &dummy, &dummy, &dummy, &dummy,
               &dummy, &dummy, &dummy, &dummy, &eValue, &bitScore);
    // Ensure the data looks good.
    if ((fieldsScanned != 12) || (eValue < 0) || (bitScore < 0) ||
        (strlen(genContigID) < 1) || (strlen(srcTransID) < 1)) {
        throw EXP(4, "Invalid Megablast output encountered.",
                  "Possibly MegaBlast version is not compatible");
    }
    // Data looks good. return a MegaBlastData object back.
    return MegaBlastData(genContigID, srcTransID, eValue, bitScore);
}

void
MegaBlastRunner::setBestEntry(MBDSet& bestSoFar,
                              const MegaBlastData& entry) const {
    MBDSet::iterator currEntry = bestSoFar.find(entry);
    if (currEntry == bestSoFar.end()) {
        // We don't have an entry for this contig.
        bestSoFar.insert(entry);
    } else if ((currEntry->getEvalue() > entry.getEvalue()) ||
               ((currEntry->getEvalue() == entry.getEvalue()) &&
                (currEntry->getBitScore() < entry.getBitScore()))) {
        // The new entry either has:
        // 1. A better e-value OR
        // 2. Same e-value but better bit score.
        bestSoFar.erase(currEntry);
        bestSoFar.insert(entry);
    }
}

#endif
