#ifndef ASSEMBLER_CPP
#define ASSEMBLER_CPP

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

#include "Assembler.h"
#include "SubSystem.h"
#include "RuntimeContext.h"
#include "Utilities.h"
#include "ESTList.h"

Assembler::Assembler(const std::string& name)
    : Component(name), estList(NULL) {
    // Intialize pointers to dummy value.
}


Assembler::~Assembler() {
    // An empty destructor, symmetric with an empty constructor.
    // Symmetry -- Hall mark of a good design?
}

void
Assembler::addCommandLineArguments(ArgParser& argParser) {
    // Currently nothing important to be added to the argument parser
    // but for a message.
    const ArgParser::ArgRecord CommonArgsList[] = {
        {"", "Arguments for " + getName() + " assembler:",
         NULL, ArgParser::INFO_MESSAGE},
        {"--progress", "Log assembly progress in a file (used by GUI)",
         &progFileName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(CommonArgsList);
}

void
Assembler::setESTList(ESTList* estList) {
    this->estList = estList;
}

bool
Assembler::initialize() {
    if (!Component::initialize()) {
        // Base class inititalization failed
        return false;
    }
    // Setup our convenience pointer to shared ESTList
    if (estList == NULL) {
        ASSERT( subSystem != NULL );
        ASSERT( subSystem->getContext() != NULL );
        ASSERT( subSystem->getContext()->getESTList() != NULL );
        estList = subSystem->getContext()->getESTList();
    }
    ASSERT( estList != NULL );
    // Try and create the progress file if specified.  The static
    // instance variable progFileName is defined in the base class.
    // Possibly that variable and corresponding parameter can be moved
    // to this class.
    if (!progFileName.empty()) {
        progressFile.open(progFileName.c_str());
        if (!progressFile.good()) {
            std::cerr << "Unable to open progress file "
                      << progFileName << ". Aborting." << std::endl;
            return false;
        }
    }    
    // Return success.
    return true;
}

void
Assembler::reportProgress(const int estsProcessed) {
    const int estCount = estList->size();
    if (progressFile.good()) {
        progressFile << (estCount - estsProcessed)
                     << ","  << estCount
                     << "\n" << std::flush;
        progressFile.seekp(0);
    }
}

Assembler& 
Assembler::operator=(const Assembler&) {
    return *this;
}

void
Assembler::addContigListener(ContigListener* cl) {
    contigListeners.push_back(cl);
}

bool
Assembler::notifyContigListeners(const Contig& contig,
                                 const bool fullContig) const {
    bool result = true;
    for(size_t i = 0; (i < contigListeners.size()); i++) {
        result |= contigListeners[i]->contigFormed(*this, contig, fullContig);
    }
    return result;
}

#endif
