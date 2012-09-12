#ifndef BLAST_DB_MAKER_CPP
#define BLAST_DB_MAKER_CPP

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

#include "BlastDBMaker.h"
#include "Helpers.h"
#include <cstdio>
#include <cerrno>
#include <sstream>
#include <sys/stat.h>

BlastDBMaker::BlastDBMaker() {
    // Nothing to be done for now.
}

BlastDBMaker::~BlastDBMaker() {
    // Empty constructor begets an empty destructor.
}

std::string
BlastDBMaker::checkMakeDB(const std::string& srcTransFile)
    throw(DecagonException) {
    // Get the source file's timestamp (assuming it exists)
    const time_t srcTimestamp = getTimestamp(srcTransFile);
    const std::string dbName  = getDBName(srcTransFile) + "_blast_db";
    try {
        // Check if the BLAST DB sub-files are up to date.
        if ((srcTimestamp < getTimestamp(dbName + ".nhr"))  &&
            (srcTimestamp < getTimestamp(dbName + ".nin")) &&
            (srcTimestamp < getTimestamp(dbName + ".nsq"))) {
            // The DB files exist and are more recent than the source
            // sequence. Things look good. Nothing further to do.
            return dbName;
        }
    } catch (const DecagonException& de) {
        // We expect exception to be thrown if the BLAST DB does not
        // exist.
    }
    // When control drops here we need to get BLAST to build a new db
    // for us. Obtain appropriate command-line from the base class.
    const std::string cmd = getBlastDBCmdLine(srcTransFile, dbName);
    // Run the command and check to ensure it ran fine
    if (!Helpers::runCmd(cmd)) {
        // When control drops here the call to make a BLAST data base
        // failed!
        const std::string msg = "BLAST database for reference sequence could "
            "not be built\n[Cmd run: " + cmd + "]";
        throw EXP(2, msg.c_str(),
                  "Verify BLAST DB path and transcript file path are "
                  "correct.\nVerify you have necessary write privileges.");
    }
    return dbName;
}

time_t
BlastDBMaker::getTimestamp(const std::string& fileName) const
    throw(DecagonException) {
    struct stat fileInfo;
    if (stat(fileName.c_str(), &fileInfo) != 0) {
        // This is a problem.
        std::string errorMsg =
            "Unable to determine information (via stat()) for file: ";
        errorMsg += fileName;
        throw EXP(errno, errorMsg.c_str(),
                  "Check path and file name are correct. "
                  "Check permissions on all directories in path");
    }
    // Return the last modification timestamp
    return fileInfo.st_mtime;
}

std::string
BlastDBMaker::getDBName(const std::string& srcTransFile) const {
    // Locate the last '.' in the source file.
    const size_t dotPos = srcTransFile.rfind('.');
    // Strip all information after it and return the string.
    if ((dotPos != std::string::npos) && (dotPos > 0)) {
        return srcTransFile.substr(0, dotPos);
    }
    // By default simply use the srcTransFile as the base file name.
    return srcTransFile;
}

#endif
