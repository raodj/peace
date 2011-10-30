#ifndef ASSEMBLY_SUB_SYSTEM_CPP
#define ASSEMBLY_SUB_SYSTEM_CPP

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

#include "MPIHelper.h"
#include "AssemblySubSystem.h"
#include "AssemblerFactory.h"
#include "RuntimeContext.h"
#include "ArgParser.h"
#include "Assembler.h"
#include "Utilities.h"
#include "ESTList.h"

AssemblySubSystem::AssemblySubSystem() {
    // Nothing else to be done for now.
}

AssemblySubSystem::~AssemblySubSystem() {
    // Nothing to be done (at least for now)
}

void
AssemblySubSystem::addCommandLineArguments(ArgParser& argParser) {
    // Add command-line parameters and additional information
    // regarding heuristics.
    const ArgParser::ArgRecord LocalArgs[] = {
        {"", "\nValid arguments for the assembly sub-system are:",
         NULL, ArgParser::INFO_MESSAGE},
        {"--assembler", "Name of the assembler to use. "
         "Valid names are:",
         &assemblerName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}        
    };
    argParser.addValidArguments(LocalArgs);
    // Have the factory add informational entries regarding the list
    // of valid assemblers
}

int
AssemblySubSystem::initializeSubSystem(ArgParser& argParser) {
    if (assemblerName.empty()) {
        // No assembler specified. Exit out
        return 0;
    }
    // Try and instantiate the assembler the user has specified.
    Assembler* assembler = AssemblerFactory::create(assemblerName,
                                                    MPI_GET_RANK());
    if (assembler == NULL) {
        // invalid assembler name
        return 1;
    }
    // Setup sub-system cross reference in the assembler
    assembler->setSubSystem(this);
    // Let the assembler add its own set of parameters
    assembler->addCommandLineArguments(argParser);
    // Setup the cross reference to the assembler in the runtime context.
    ASSERT ( runtimeContext != NULL );
    runtimeContext->setAssembler(assembler);
    // Everything went well.
    return NO_ERROR;
}

int
AssemblySubSystem::initializeSubComponents() {
    ASSERT ( runtimeContext != NULL );
    if (runtimeContext->getAssembler() == NULL) {
        // Nothing to be done.
        return NO_ERROR;
    }
    // setup the EST list for the assembler to use
    runtimeContext->getAssembler()->setESTList(runtimeContext->getESTList());
    // Let the assembler do its initialization
    if (!runtimeContext->getAssembler()->initialize()) {
        // Initializatin of the assembler failed.
        return 1;
    }
    return NO_ERROR;
}

int
AssemblySubSystem::run() {
    ASSERT ( runtimeContext != NULL );
    if (runtimeContext->getAssembler() == NULL) {
        // Nothing to be done.
        return NO_ERROR;
    }
    // Check to ensure we have some ESTs to process.
    if ((runtimeContext->getESTList() == NULL) ||
        (runtimeContext->getESTList()->size() < 1)) {
        // Don't have ESTs to process.
        std::cerr << "Error: The list of cDNAs to be processed is empty.\n"
                  << "Use --estFiles option to load data from file(s).\n";
        return 5;
    }
    // Let the assembler do its thing
    if (runtimeContext->getAssembler()->assemble()) {
        // Assembler had trouble...
        return 2;
    }    
    // Everything went on well
    return NO_ERROR;
}

void
AssemblySubSystem::finalizeComponents(const bool) {
    ASSERT ( runtimeContext != NULL );
    if (runtimeContext->getAssembler() != NULL) {    
        runtimeContext->getAssembler()->finalize();
    } 
}

void
AssemblySubSystem::finalizeSubSystem(const bool) {
    // Now get rid of the dynamically allocated components
    ASSERT ( runtimeContext != NULL );    
    if (runtimeContext->getAssembler() != NULL) {
        delete  runtimeContext->getAssembler();
    }
}

#endif
