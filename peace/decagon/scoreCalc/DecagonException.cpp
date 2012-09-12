#ifndef EXCEPTION_CPP
#define EXCEPTION_CPP


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

#include "DecagonException.h"
#include "Utilities.h"
#include <cstring>

DecagonException::DecagonException(const int errCode, const char* errorMessage,
                                   const char *help,  const char *sourceFile,
                                   const int  line) :
    errorCode(errCode), lineNumber(line) {
    
    // Ensure that none of the string values are NULL.  Note that
    // string values can be empty (\i i.e., "") but not NULL
    ASSERT ( errorMessage != NULL );
    ASSERT ( help         != NULL );
    ASSERT ( sourceFile   != NULL );
    
    // Let's run a super sanity check to ensure that the lengths of
    // the various strings will fit within our static buffers.
    ASSERT ( strlen(errorMessage) < sizeof(message)   );
    ASSERT ( strlen(help)         < sizeof(suggestion));
    ASSERT ( strlen(sourceFile)   < sizeof(fileName)  );
    
    // Copy the necessary string values into corresponding fields.
    strcpy(message,    errorMessage);
    strcpy(suggestion, help);
    strcpy(fileName,   sourceFile);
}

int
DecagonException::getErrorCode() const {
    return errorCode;
}

const char*
DecagonException::getErrorMessage() const {
    return message;
}

const char*
DecagonException::getSuggestion() const {
    return suggestion;
}

const char*
DecagonException::getFileName() const {
    return fileName;
}

int
DecagonException::getLineNumber() const {
    return lineNumber;
}

DecagonException::~DecagonException() {
    // The DecagonException class does not allocate any memory or
    // perform signficiant operations to avoid auxiliary exceptions
    // from being generated.  Consequently, the destructor has no
    // tasks to perform.  It is here just to adhere to the standard
    // convention of having a constructor-destructor pair for every
    // class.
}

std::ostream&
operator<<(std::ostream& os, const DecagonException& de) {
    os << "Exception [code: "
       << de.getErrorCode() << "]: " << de.getErrorMessage() << std::endl;
    os << "Suggestion: " << de.getSuggestion() << std::endl;
    os << "Source Code Details: " << de.getFileName() << ":"
       << de.getLineNumber() << std::endl;
    return os;
}

#endif
