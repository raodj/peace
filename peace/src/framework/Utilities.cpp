#ifndef STRING_UTILITIES_CPP
#define STRING_UTILITIES_CPP

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

#include "Utilities.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <cstring>

char*
getTimeStamp(const char *fileName, char *buffer) {
    if (fileName == NULL) {
        // Nothing further to be done here.
        return buffer;
    }
    // The follwing structure will contain the file information.
    int returnValue = 0;
#ifdef _WINDOWS
    struct _stat fileInfo;
    returnValue = _stat(fileName, &fileInfo);
#else
    struct stat fileInfo;
    // The only difference between windows and Linux is the extra "_"
    // at the beginning of stat() method.
    returnValue = stat(fileName, &fileInfo);
#endif
    // Check to ensure there were no errors.  If there were errors
    // exit immediately.
    if (returnValue == -1) {
        // O!o! there was an error.
        return buffer;
    }
    // Convert the last modification time to string and return it back
    // to the caller.
    return getTime(buffer, &fileInfo.st_mtime);
}

char* getTime(char *buffer, const time_t *encodedTime) {
    if (buffer == NULL) {
        // Nothing more to do.
        return NULL;
    }
    // Use encodedTime or the system time.
    const time_t timeToConv = (encodedTime != NULL) ? *encodedTime : time(NULL);
    // Convert the time.
    ctime_s(buffer, 128, &timeToConv);
    // Remove trailing new line (if one is present)
    const int strLen = strlen(buffer);
    if (buffer[strLen - 1] == '\n') {
        buffer[strLen - 1] = '\0';
    }
    // Return the buffer back to the caller
    return buffer;
}

#endif

