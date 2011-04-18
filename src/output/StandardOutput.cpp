#ifndef STANDARD_OUTPUT_CPP
#define STANDARD_OUTPUT_CPP

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

#include "StandardOutput.h"
#include "ArgParser.h"
#include "LogLevel.h"
#include "MPIHelper.h"
#include "RedirectingStreamBuf.h"

#include <sstream>
#include <fstream>

// Define the static log level
LogLevel::LevelTypes StandardOutput::logLevel = LogLevel::LOG_LEVEL_INFO;

StandardOutput::StandardOutput() : Component("StandardOutput"),
                                   logLevelStr("info") {
    coutStream           = NULL;
    clogStream           = NULL;
    clogWrapperStreamBuf = NULL;
    coutSystemStreamBuf  = NULL;
    clogSystemStreamBuf  = NULL;
}

StandardOutput::~StandardOutput() {
    coutStream           = NULL;
    clogStream           = NULL;
    clogWrapperStreamBuf = NULL;
    coutSystemStreamBuf  = NULL;
    clogSystemStreamBuf  = NULL;
}

void
StandardOutput::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"--stdout", "Redirect standard output to a given file",
         &coutFileName, ArgParser::STRING},
        {"--stdlog","Redirect standard log to a given file",
         &clogFileName, ArgParser::STRING},
        {"--logLevel", "Set the log level (levels: none, info, debug)",
         &logLevelStr, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add arguments for this component to the parser
    argParser.addValidArguments(Arguments);
}

bool
StandardOutput::initialize() {
    // Ensure the log level string is valid and setup the log level
    // constant object.
    if (logLevelStr == "none") {
        logLevel = LogLevel::LOG_LEVEL_NONE;
    } else if (logLevelStr == "info") {
        logLevel = LogLevel::LOG_LEVEL_INFO;
    } else if (logLevelStr == "debug") {
        logLevel = LogLevel::LOG_LEVEL_DEBUG;
    } else {
        std::cerr << "Invalid log level specified. Valid log levels are: "
                  << "none, info, debug."
                  << std::endl;
        return false;
    }
    
    // Check and redirect standard log stream to a file.
    if (!clogFileName.empty()) {
        if ((clogStream = redirect(std::clog, clogSystemStreamBuf,
                                   clogFileName)) == NULL) {
            // Error occured during redirection
            return false;
        }
    }
    // Check and redirect standard output stream to a file.
    if (!coutFileName.empty()) {
        if ((coutStream = redirect(std::cout, coutSystemStreamBuf,
                                   coutFileName)) == NULL) {
            // Error occured during redirection
            return false;
        }
    }
    // Save reference to the current log stream buffer for use in
    // operator<<() method.
    clogWrapperStreamBuf = new RedirectingStreamBuf(std::clog.rdbuf());
    // Finally approriately initialize log level for std::clog.
    std::clog << LogLevel::INFO_LOG;
    // All went well.
    return true;
}

void
StandardOutput::finalize() {
    // Restore the stream buffers for std::cout and std::clog back to
    // the system defaults.
    std::cout.rdbuf(coutSystemStreamBuf);
    std::clog.rdbuf(clogSystemStreamBuf);
    // Close out any output stream if we have them open.
    if (coutStream != NULL) {
        delete coutStream;
        coutStream = NULL;
    }
    if (clogStream != NULL) {
        delete clogStream;
        clogStream = NULL;
    }
    // Free up our clog wrapper
    delete clogWrapperStreamBuf;
}

std::ostream*
StandardOutput::redirect(std::ostream& os, std::streambuf* &currStreamBuf,
                         const std::string& origFileName) {
    // First create a suitable file name with the MPI rank.
    std::ostringstream fileNameStream;
    fileNameStream << origFileName << "_" << MPI_GET_RANK();
    // Get the new file name.
    const std::string fileName = fileNameStream.str();
    // Try and open the file for writing.
    std::ofstream* targetFile = new std::ofstream(fileName.c_str());
    if ((targetFile == NULL) || (!targetFile->good())) {
        std::cerr << "Error opening file " << fileName
                  << " for writing." << std::endl;
        delete targetFile;
        targetFile = NULL;
    } else {
        // File was opened successfully. Redirect the stream.
        currStreamBuf = os.rdbuf(targetFile->rdbuf());
    }
    // Return the opened stream back to the caller
    return targetFile;
}

std::ostream&
operator<<(std::ostream& os, const LogLevel& ll) {
    // Setup the stream buffer.
    RedirectingStreamBuf *rsb = dynamic_cast<RedirectingStreamBuf*>(os.rdbuf());
    if (rsb != NULL) {
        // Enable logging if the current log level is less than the
        // one set by the user.
        rsb->redirect(ll >= StandardOutput::logLevel);
    }
    return os;
}

#endif
