#ifndef BLAST_COMMON_CPP
#define BLAST_COMMON_CPP

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

#include "BlastCommon.h"
#include "Helpers.h"

// A couple of defines to make the code more readable
#define OLD_BLAST_TYPE 0
#define NEW_BLAST_TYPE 1

// The static instance variable to detect BLAST type
int BlastCommon::blastType = -1;

std::string
BlastCommon::getBlastDBCmdLine(const std::string& srcTransFile,
                               const std::string& dbName)
    throw (DecagonException) {
    // Ensure we precisely know the type of BLAST
    detectBlastType();
    // Build command-line based on BLAST type
    switch(blastType) {
    case OLD_BLAST_TYPE:
        return "formatdb -p F -i " + srcTransFile + " -n " + dbName;
    case NEW_BLAST_TYPE:
        return "makeblastdb -dbtype nucl -in " + srcTransFile +
            " -out " + dbName;
    default:
        throw EXP(202, "Unhandled BLAST type encountered",
                  "This is an internal DECAGON problem. "
                  "Please file this issue as a defect");
    }
}

std::string
BlastCommon::getMegablastCmdLine(const std::string& blastDBName,
                                 const std::string& genContigFile)
    throw (DecagonException) {
    // Ensure we precisely know the type of BLAST
    detectBlastType();
    // Build command-line based on BLAST type
    switch(blastType) {
    case OLD_BLAST_TYPE:
        return "megablast -m9 -i " + genContigFile + " -d " + blastDBName;
    case NEW_BLAST_TYPE:
        return "blastn -task megablast -outfmt 7 -db " + blastDBName +
            " -query " + genContigFile;
    default:
        throw EXP(202, "Unhandled BLAST type encountered",
                  "This is an internal DECAGON problem. "
                  "Please file this issue as a defect");
    }    
}

void
BlastCommon::detectBlastType() throw(DecagonException) {
    if (blastType != -1) {
        // The type of BLAST is already known. Nothing further to do.
        return;
    }
    // Try to run the older form of BLAST to see if it succeeds
    if (Helpers::runCmd("formatdb --help", 1)) {
        // Success! we are working with older BLAST version
        blastType = OLD_BLAST_TYPE;
        return;
    }
    // Try to run the newer form of BLAST command to see if that helps
    if (Helpers::runCmd("makeblastdb -h", 0)) {
        // Success! we are working with newer BLAST+ version
        blastType = NEW_BLAST_TYPE;
        return;
    }
    // When control drops here we either don't have BLAST or it is not
    // a compatible version.
    throw EXP(200, "Unable to detect BLAST type (BLAST vs BLAST+)",
              "Ensure BLAST is installed, module is loaded (if needed), "
              "and the various BLAST-executables are in your default path.");
}


#endif
