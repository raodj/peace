#ifndef COMPONENT_CPP
#define COMPONENT_CPP

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

#include "Component.h"
#include "SubSystem.h"
#include "MPIHelper.h"

Component::Component(const std::string& compName) :
    subSystem(NULL), name(compName), initializedFlag(false) {
    // Nothing else to be done in the constructor.
}

Component::~Component() {
    // Empty constructor begets an empty destructor.
}

void
Component::addCommandLineArguments(ArgParser&) {
    // Nothing to be done here.
}

RuntimeContext*
Component::getContext() const {
    return subSystem->getContext();
}

bool
Component::initialize() {
    initializedFlag = true;
    return true;
}

int
Component::run() {
    return 0;
}

void
Component::setSubSystem(SubSystem* ss) {
    subSystem = ss;
}

bool
Component::setArgument(const std::string& arg, const std::string& value) {
    return setArgument(arg, ArgParser::STRING, value);
}

bool
Component::setArgument(const std::string& arg, const bool value) {
    return setArgument(arg, ArgParser::BOOLEAN, value);
}

bool
Component::setArgument(const std::string& arg, const int value) {
    return setArgument(arg, ArgParser::INTEGER, value);
}

bool
Component::setArgument(const std::string& arg, const float value) {
    return setArgument(arg, ArgParser::FLOAT, value);
}

bool
Component::setArgument(const std::string& arg, const double value) {
    return setArgument(arg, ArgParser::DOUBLE, value);
}

template <typename ValueType>
bool
Component::setArgument(const std::string& arg,
                       const ArgParser::ArgType argType,
                       const ValueType& value) {
    // First gather all valid arguments from this class hierarchy
    ArgParser ap;
    addCommandLineArguments(ap);
    // Now let the arg parser update the valid parameter
    return ap.setArgument(arg, argType, value);
}

void
Component::showArguments(std::ostream& os) {
    // First gather all valid arguments from this class hierarchy
    ArgParser ap;
    addCommandLineArguments(ap);
    // Now display them
    os << ap;
}

void
Component::getOwnedESTidx(const int estCount, int& startIndex, int& endIndex) {
    const int ESTsPerProcess = estCount / MPI_GET_SIZE();
    const int ExtraESTs      = estCount % MPI_GET_SIZE();
    const int MyRank         = MPI_GET_RANK();
    
    // First figure out the starting and ending EST this processs is
    // responsible for further use.
    startIndex = MyRank * ESTsPerProcess;
    // Check for extra proceses as needed.
    if (MyRank <= ExtraESTs) {
        // The previous processes have one extra ESTs as number of
        // ESTs are not evenly divisible by number of processes.  So
        // account for this.
        startIndex = ((ESTsPerProcess + 1) * MyRank);
    } else {
        startIndex += ExtraESTs;
    }
    
    // Compute the last est index this process owns.
    endIndex = startIndex + ESTsPerProcess;
    if (MyRank < ExtraESTs) {
        // This process owns one extra EST as ESTs are not evenly
        // divisible by the number of processes.
        endIndex++;
    }
}

#endif
