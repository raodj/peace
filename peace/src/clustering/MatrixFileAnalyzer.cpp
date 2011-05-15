#ifndef MATRIX_FILE_ANALYZER_CPP
#define MATRIX_FILE_ANALYZER_CPP

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

#include "MatrixFileAnalyzer.h"
#include "ArgParser.h"
#include "ESTList.h"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <cstdlib>

MatrixFileAnalyzer::MatrixFileAnalyzer()
    : ESTAnalyzer("MatrixFileAnalyzer") {
    distanceValues = NULL;
    estCount       = 0;
}

MatrixFileAnalyzer::~MatrixFileAnalyzer() {
    if (distanceValues != NULL) {
        // Delete each row and then delete the whole array.
        for(int row = 0; (row < estCount); row++) {
            delete[] distanceValues[row];
        }
        delete[] distanceValues;
    }
    distanceValues = NULL;
    estCount       = 0;
}

void
MatrixFileAnalyzer::addCommandLineArguments(ArgParser& argParser) {
    // The arguments specific for MatrixFileAnalyzer...
    const ArgParser::ArgRecord ArgsList[] = {
        {"--dataFile", "Data file containing matrix of distance metrics",
         &dataFileName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(ArgsList);
}

bool
MatrixFileAnalyzer::initialize() {
    // Check if necessary arguments have been specified for processing.
    if (dataFileName.empty()) {
        std::cerr << getName()
                  << ": Matrix data file not specified "
                  << "(use --dataFile option)\n";
        return false;
    }
    
    FILE *input = fopen(dataFileName.c_str(), "rt");
    if ((input == NULL) || ferror(input)) {
        // Error even opening the data file. 
        std::cerr << "Error opening matrix data file "
                  << dataFileName << " for reading." << std::endl;
        return false;
    }
    int  currCol   = 0;
    int  currRow   = 0;  // Current row in distanceValues array
    while (!feof(input)) {
        // Read a line from the input file.
        std::string line = readLine(input);
        if ((line.length() < 1) || (line[0] == '#')) {
            // Blank or comment lines.
            continue;
        }
        // Process an non-empty/non-comment line.
        if (distanceValues == NULL) {
            if (!parseESTCount(line.c_str())) {
                // Could not parse the estCount properly.
                fclose(input);
                return false;
            }
        } else {
            // This line must have some distance values.
            if (currRow >= estCount) {
                // Too much data!
                std::cerr << "Excess data in matrix file.\n";
                fclose(input);
                return false;
            }
            currCol += parseMetrics(line.c_str(), distanceValues[currRow],
                                    currCol, estCount - currCol);
            if (currCol == estCount) {
                // Read a whole row. Move to next row.
                currCol = 0;
                currRow++;
            }
        }
    }
    // Done with the file.
    fclose(input);
    // Check if we read enough data for further processing
    if (currRow != estCount) {
        std::cerr << "Insufficient data in matrix file for "
                  << estCount << " ESTs.\n";
        return false;
    }
    // When control drops here that mean everything went well.
    return true;
}

bool
MatrixFileAnalyzer::parseESTCount(const char *line) {
    char *endPtr = NULL;
    estCount     = strtol(line, &endPtr, 10);
    if ((endPtr == NULL) || (*endPtr != '\0') ||
        (estCount < 0) || (estCount > 100000)) {
        std::cerr << "Invalid EST count value encountered\n(line = '"
                  << line << "')\n";
        estCount = 0;
        return false;
    }
    // Obtain reference to the global list of ESTs to be processed
    ASSERT( ESTAnalyzer::estList != NULL );
    ESTList& estList = *ESTAnalyzer::estList;
    // Create the distanceValues matrix to hold estCount * estCount
    // values.
    distanceValues = new float*[estCount];
    for(int rows = 0; (rows < estCount); rows++) {
        distanceValues[rows] = new float[estCount];
        memset(distanceValues[rows], 0, sizeof(float) * estCount);
        // For each row create a dummy EST entry to make display and
        // processing easier.
        char info[32];
        sprintf(info, "Dummy EST #%d", rows);
        estList.add(rows, info, "No sequence data available");
    }
    
    // Everything went on fine.
    return true;
}

int
MatrixFileAnalyzer::setReferenceEST(const EST* est) {
    refEST = est;
    return (est != NULL) ? 0 : 1;
}

float
MatrixFileAnalyzer::getMetric(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    ASSERT( refEST   != NULL );
    return distanceValues[refEST->getID()][otherEST->getID()];
}

int
MatrixFileAnalyzer::parseMetrics(const char* line, float *values,
                                 const int startPos, const int maxValues) {
    int index = 0;
    // Create a istringstream to easily read a bunch of float values out
    std::string tempLine = line;
    std::istringstream inStream(line);
    while (!inStream.eof()) {
        float metric;
        inStream >> metric;
        // Store the metric into the array.
        values[startPos + index] = metric;
        index++;
        if (index >= maxValues) {
            // We have read the maximum number.
            return index;
        }
    }
    // Return the number of values read.
    return index;
}

int
MatrixFileAnalyzer::analyze() {
    return -1;
}

std::string
MatrixFileAnalyzer::readLine(FILE *fp) {
    std::string line;
    if (feof(fp) || ferror(fp)) {
        // The file pointer is already invalid. Can't do much further
        return line;
    }
    char buffer[1024];
    while (fgets(buffer, 1024, fp) != NULL) {
        line += buffer;
        const size_t len = strlen(buffer);
        if (buffer[len - 1] == '\n') {
            // OK, we have read up to end-of-line.
            break;
        }
    }
    // Return the complete line back to the caller.
    return line.substr(0, line.length() - 1);
}

#endif
