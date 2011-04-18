#ifndef PEACE_CPP
#define PEACE_CPP

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

#include "PEACE.h"
#include "ArgParser.h"
#include "Version.h"
#include "MPIHelper.h"

#include "OutputSubSystem.h"
#include "InputSubSystem.h"
#include "FilteringSubSystem.h"
#include "ClusteringSubSystem.h"
#include "AssemblySubSystem.h"

#include <algorithm>
#include <functional>

#define PEACE_BRIEF_MESSAGE                                       \
    FULL_TITLE                                                    \
    "\n"                                                          \
    PEACE_VERSION                                                 \
    ". Released (under GPLv3) on "                                \
    RELEASE_DATE                                                  \
    "\n"                                                          \
    COPYRIGHT                                                     \
    "\nUsage: peace [arguments]\n"                                \
    "where valid arguments are:"

PEACE::PEACE() {
    // Add references to the various sub-systems into the common
    // vector to ease processing.
    subSystemList.push_back((oss = new OutputSubSystem()));
    subSystemList.push_back((iss = new InputSubSystem()));
    subSystemList.push_back((fss = new FilteringSubSystem()));
    subSystemList.push_back((css = new ClusteringSubSystem()));
    subSystemList.push_back((ayss= new AssemblySubSystem()));
    // Setup the shared runtime configuration information
    std::for_each(subSystemList.begin(), subSystemList.end(),
                  std::bind2nd(std::mem_fun(&SubSystem::setContext),
                               &runtimeContext));
    // Initialize pointers and other native data types
    argParser = NULL;
    errorCode = 0;
}

PEACE::~PEACE() {
    if (argParser != NULL) {
        finalize(false);
    }
    // Delete the sub-systems
    delete oss;
    delete iss;
    delete fss;
    delete css;
    delete ayss;
}

int
PEACE::initialize(int& argc, char *argv[], const bool initMPI) {
    static bool interactive = false; // Fake variable 
    // First initialize MPI if requested
    if (initMPI) {
        MPI_INIT(argc, argv);
    }
    if (argParser == NULL) {
        // Create the short-lived argument parser.
        argParser = new ArgParser();
    }
    // Add information to the runtime configuration
    runtimeContext.addCmdLineInfo(argc, argv);
    // Next, setup our global arguments
    const ArgParser::ArgRecord GlobalArgs[] = {
        {"", PEACE_BRIEF_MESSAGE, NULL, ArgParser::MAIN_MESSAGE},
        {"--interactive", "Run PEACE in interactive mode",
         &interactive, ArgParser::BOOLEAN},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser->addValidArguments(GlobalArgs);
    
    // Have various sub-systems populate the argument parser
    // std::for_each(subSystemList.begin(), subSystemList.end(),
    //           std::bind2nd(std::mem_fun(&SubSystem::addCommandLineArguments),
    //                        *argParser));
    for(size_t i = 0; (i < subSystemList.size()); i++) {
        subSystemList[i]->addCommandLineArguments(*argParser);
    }
    // Let the argument parser parse the initial set of arguments
    argParser->parseArguments(argc, argv, false);
    // Now initialize all the subsystems in the appropriate
    // order. Here we call initialize until we find a sub-system that
    // returns a non-zero value.
    int errorCode = 0;
    for(size_t i = 0; (i < subSystemList.size()); i++) {
        if ((errorCode = subSystemList[i]->initializeSubSystem(*argParser)) != 0) {
            // Error during initialization of a sub-system.
            break;
        }
    }
    if (errorCode != 0) {
        // Error during initialization. Display the usage for conveniecne
        std::cout << *argParser;
    }
    // Perform a second round of argument parsing as some of the
    // sub-systems may have added additional parameters during
    // initialization (for components or sub-components).
    argParser->parseArguments(argc, argv, false);    
    // If at the end we still have command line parameters, then that
    // is a boo..boo.
    if (!argParser->checkRemainingArguments(argc, argv, false, true)) {
        return 3;
    }
    // Now initialize components
    for(size_t i = 0; (i < subSystemList.size()); i++) {
        if ((errorCode = subSystemList[i]->initializeComponents()) != 0) {
            // Error during initialization of a sub-system.
            return errorCode;
        }
    }
    // Let all the subsystems load their inputs.
    for(size_t i = 0; (i < subSystemList.size()); i++) {
        if ((errorCode = subSystemList[i]->loadInputs()) != NO_ERROR) {
            // Error initializing the subsystem
            return errorCode;
        }
    }
    // Let all the components/sub-components in various sub-systems
    // complete all initialization tasks.
    for(size_t i = 0; (i < subSystemList.size()); i++) {
        if ((errorCode = subSystemList[i]->initializeSubComponents()) != NO_ERROR) {
            // Error initializing the subsystem
            return errorCode;
        }
    }
    
    // Everything went well
    return NO_ERROR;
}

int
PEACE::run(int& argc, char *argv[], bool initMPI) {
    if (argParser == NULL) {
        // Seems like the system was not initialized earlier. Do it
        // now
        if ((errorCode = initialize(argc, argv, initMPI)) != 0) {
            // Error occurred during initialization. Bail out.
            return errorCode;
        }
    }

    // Let all the subsystems run and do their jobs.
    for(size_t i = 0; (i < subSystemList.size()); i++) {
        if ((errorCode = subSystemList[i]->run()) != NO_ERROR) {
            // Error during running of a sub-system.
            break;
        }
    }
    const bool overallResult = (errorCode == 0);    
    for(size_t i = 0; (i < subSystemList.size()); i++) {
        subSystemList[i]->generateOutputs(overallResult);
    }    
    return errorCode;    
}

void
PEACE::finalize(bool finalizeMPI) {
    // Get rid of the arg parser as we no longer need it
    if (argParser != NULL) {
        delete argParser;
        argParser = NULL;
    }
    // Now finalize all the components followed by subsystems.
    const bool overallResult = (errorCode == 0);    
    // now finalize all the components
    std::for_each(subSystemList.begin(), subSystemList.end(),
                  std::bind2nd(std::mem_fun(&SubSystem::finalizeComponents),
                               overallResult));
    std::for_each(subSystemList.begin(), subSystemList.end(),
                  std::bind2nd(std::mem_fun(&SubSystem::finalizeSubSystem),
                               overallResult));
    // Finalize MPI if requested
    if (finalizeMPI) {
        MPI_FINALIZE();
    }
}

bool
PEACE::loadFile(const std::string& fileName, const std::string& format) {
    ArgParser::StringList tmpFileList;
    tmpFileList.push_back(fileName);
    return iss->getInputFileFactory().loadFiles(tmpFileList, format);
}

#endif
