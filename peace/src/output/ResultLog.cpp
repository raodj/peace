#ifndef RESULT_LOG_CPP
#define RESULT_LOG_CPP

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

#include "ResultLog.h"
#include "HTMLHeader.h"
#include "Utilities.h"
#include <iostream>
#include <cstdarg>

ResultLog::ResultLog(const std::string& fileName, const bool htmlFlag)
    : output(NULL), html(htmlFlag), insideTable(false) {
    if (fileName.size() != 0) {
#ifndef _WINDOWS
        output = fopen(fileName.c_str(), "wt");
#else
        fopen_s(&output, fileName.c_str(), "wt");
#endif
        if ((output == NULL) || (ferror(output))) {
            // Error occured when attempting to create file.
            output = NULL;
            std::cerr << "Error occured when attempting to log information "
                      << "to " << fileName << "\n"
                      << "Results will be written to standard output.\n";
        }
    }
    if (output == NULL) {
        // Write results to standard output.
        const int StandardOutput = 1;
#ifndef _WINDOWS
        output = fdopen(StandardOutput, "wt");
#else
        output = _fdopen(StandardOutput, "wt");
#endif
    }

    // Create necessary html headers.
    dumpHTMLHeaders();
}

void
ResultLog::reportLine(const char *line, ...) {
    // First process any variable number of arguments and obtain final
    // string in a local buffer.
    char buffer[2048] = "\0";
    va_list arguments;
    va_start(arguments, line);
    vsnprintf_s(buffer, sizeof(buffer), sizeof(buffer), line, arguments);
    va_end(arguments);
    
    if (html) {
        // End any table that may be open.
        endTable();
        // Dump with trailing BR.
        fprintf(output, "%s<br>\n", buffer);
    } else {
        fprintf(output, "%s\n", buffer);
    }
}

void
ResultLog::report(const char *col1, const char *col2, const char *col3, ...) {
    // First create a single string with all the necessary format
    // specifiers set for call to the vsnprintf_s() method below.
    std::string tableRow;
    tableRow  = (html ? "<tr><td>" : "");
    tableRow += col1;
    tableRow += (html ? "</td><td>" : "  ");
    tableRow += col2;
    tableRow += (html ? "</td><td>" : "  ");
    tableRow += col3;
    tableRow += (html ? "</td></tr>\n" : "\n");
    
    // First process any variable number of arguments and obtain final
    // string in a local buffer.
    char buffer[2048];
    va_list arguments;
    va_start(arguments, col3);
    vsnprintf_s(buffer, sizeof(buffer), sizeof(buffer),
                tableRow.c_str(), arguments);    
    va_end(arguments);
    
    // Now print the row out to the target.
    if (html) {
        // Start a table if one is not already open.
        startTable();
    }
    fprintf(output, buffer);
}

ResultLog::~ResultLog() {
    // First end any open HTML tables.
    endTable();
    // If in HTML mode dump trailing HTML tags.
    if (html) {
        fprintf(output, HTML_FOOTER);
    }
    // Close the output if it is a standard file
    if (_fileno(output) != 1) {
        fclose(output);
    }
}

void
ResultLog::dumpHTMLHeaders() {
    if (!html) {
        return;
    }
    // Do a simple header for now.
    fprintf(output, HTML_HEADER);
}

void
ResultLog::startTable(const char* titles[]) {
    if (!html || insideTable) {
        // Nothing to be done in this case.
        return;
    }
    // Dump HTML headers.
    fprintf(output, "  <table border=\"0\">\n");
    // Dump headers if specified.
    if (titles != NULL) {
        // Start header.
        fprintf(output, "  <thead>\n    <tr>");
        int index = 0;
        while (titles[index] != NULL) {
            fprintf(output, "<td>%s</td>", titles[index]);
            index++;
        }
        // Finish header
        fprintf(output, "  </thead>\n");
    }
    // Start table body
    fprintf(output, "  <tbody>\n");
    insideTable = true;
}

void
ResultLog::endTable() {
    if (!html || !insideTable) {
        // Nothing to be done in this case.
        return;
    }
    // Dump table trailers.
    fprintf(output, "  </tbody>\n");
    fprintf(output, "  </table>\n");
    // Reset table scope.
    insideTable = false;
}

ResultLog& 
ResultLog::operator=(const ResultLog&) {
    return *this;
}

#endif
