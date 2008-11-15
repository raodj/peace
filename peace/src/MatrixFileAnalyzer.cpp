#ifndef MATRIX_FILE_ANALYZER_CPP
#define MATRIX_FILE_ANALYZER_CPP

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include "MatrixFileAnalyzer.h"
#include "EST.h"

// The static instance variables for command line arguments.
char* MatrixFileAnalyzer::dataFileName  = NULL;

// The common set of arguments for all FW EST analyzers
arg_parser::arg_record MatrixFileAnalyzer::argsList[] = {
    {"--dataFile", "Data file containing matrix of distance metrics",
     &MatrixFileAnalyzer::dataFileName, arg_parser::STRING},
    {NULL, NULL}
};

MatrixFileAnalyzer::MatrixFileAnalyzer(const int refESTidx,
                                       const std::string& outputFileName)
    : ESTAnalyzer("MatrixFileAnalyzer", refESTidx, outputFileName) {
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
MatrixFileAnalyzer::showArguments(std::ostream& os) {
    ESTAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(MatrixFileAnalyzer::argsList);
    os << ap;
}

bool
MatrixFileAnalyzer::parseArguments(int& argc, char **argv) {
    arg_parser ap(MatrixFileAnalyzer::argsList);
    ap.check_args(argc, argv, false);
    // Check if necessary arguments have been specified for processing.
    if (dataFileName == NULL) {
        std::cerr << analyzerName
                  << ": Matrix data file not specified "
                  << "(use --dataFile option)\n";
        return false;
    }
    // Now let the base class do processing and return the result.
    return ESTAnalyzer::parseArguments(argc, argv);
}

int
MatrixFileAnalyzer::initialize() {
    FILE *input = fopen(dataFileName, "rt");
    if ((input == NULL) || ferror(input)) {
        // Error even opening the data file. 
        std::cerr << "Error opening matrix data file "
                  << inputFileName << " for reading." << std::endl;
        return 1;
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
                return 2;
            }
        } else {
            // This line must have some distance values.
            if (currRow >= estCount) {
                // Too much data!
                std::cerr << "Excess data in matrix file.\n";
                fclose(input);
                return 3;
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
        return 4;
    }
    // When control drops here that mean everything went well.
    return 0;
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
        EST::create(rows, info, "No sequence data available");
    }
    
    // Everything went on fine.
    return true;
}

int
MatrixFileAnalyzer::setReferenceEST(const int estIdx) {
    if ((estIdx >= 0) && (estIdx < estCount)) {
        refESTidx = estIdx;
        return 0;
    }
    // Invalid est index.
    return 1;
}

float
MatrixFileAnalyzer::analyze(const int otherEST) {
    if ((otherEST >= 0) && (otherEST < estCount)) {
        return distanceValues[refESTidx][otherEST];
    }
    // Invalid est index.
    return -1;
}

int
MatrixFileAnalyzer::parseMetrics(const char* line, float *values,
                                 const int startPos, const int maxValues) {
    int index = 0;
    // Create a mutable copy of the line for tokenization.  Memory is
    // automatically freed.
    char *dup     = strdupa(line);
    char *prevPtr = NULL;
    // Obtain the first token.
    char *token   = strtok_r(dup, " ", &prevPtr);
    while (token != NULL) {
        char *endPtr = NULL;
        float metric = strtof(token, &endPtr);
        if ((endPtr == NULL) || (*endPtr != '\0')) {
            std::cerr << "Invalid metric (" << token << ") read.\n";
            return -1;
        }
        // Store the metric into the array.
        values[startPos + index] = metric;
        index++;
        if (index >= maxValues) {
            // We have read the maximum number.
            return index;
        }
        // On to the next token
        token = strtok_r(NULL, " ", &prevPtr);
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
