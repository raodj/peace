#ifndef DRAW_MATRIX_CPP
#define DRAW_MATRIX_CPP

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

#include "DrawMatrix.h"
#include "ArgParser.h"
#include "Tool.h"

#include <fstream>
#include <sstream>
#include <cstring>

// Some color codes to make code a bit more readable
#define BLACK  0
#define WHITE  7
#define GREEN  2
#define YELLOW 6
#define GREY   32

#define CHECK_ARGS(condition, errorMsg)         \
    if (condition) {                            \
        std::cout << errorMsg;                  \
        std::cout << ap << std::endl;           \
        return 1;                               \
    }

int
DrawMatrix::main(int argc, char *argv[]) {
    std::string srcFileName; // File name with WCD matrix output
    std::string outFileName; // Output file with XFig data.
    int  estCount     = -1;  // Number of ests in the file.
    bool showOptions  = false;
    int  sqSize       = 4;
    
    // Create the list of valid arguments to be used by the arg_parser.
    ArgParser::ArgRecord arg_list[] = {
        {"--matFile", "WCD matrix output file",
         &srcFileName, ArgParser::STRING},
        {"--estCount", "Number of ESTs",
         &estCount, ArgParser::INTEGER},
        {"--output", "File to which the XFIG output must be written",
         &outFileName, ArgParser::STRING},
        {"--sqSize", "Specify size of square drawn for each entry in matrix",
         &sqSize, ArgParser::INTEGER},
        {"--options", "Lists options for this tool",
         &showOptions, ArgParser::BOOLEAN},
        {"", "", NULL, ArgParser::INVALID}
    };
    
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    ArgParser ap;
    Tool::addCmdLineArgs("ShowAlignment", ap);
    ap.addValidArguments(arg_list);
    ap.parseArguments(argc, argv, false);
    if (showOptions) {
        std::cout << ap << std::endl;
        // Nothing further to be done.
        return 0;
    }
    // Ensure we have all the necessary parameters.
    CHECK_ARGS(srcFileName.empty(), "WCD matrix output file is needed.\n" \
               "Use --matFile option\n");
    CHECK_ARGS(estCount == -1, "Number ESTs was not "  \
               "specified.\nUse --estCount option\n");
    CHECK_ARGS(outFileName.empty(), "Output XFIG file was not specified.\n" \
               "Use --output option\n");
    CHECK_ARGS(sqSize < 2, "Square size must be atleast 2.\n" \
               "Use --sqSize option\n");
    // Open input file.
    std::ifstream inFile(srcFileName.c_str());
    if (!inFile.good()) {
        std::cout << "Unable to open input file " << srcFileName
                  << " for reading.\n";
        return 2;
    }
    // OK. Create an object.
    DrawMatrix dm;
    // Open the output file.
    if (!dm.xfig.setOutput(outFileName, false)) {
        std::cout << "Unable to open output file " << outFileName
                  << " for output.\n";
        return 1;
    }
    
    // Draw the matrix.
    dm.drawMatrix(inFile, estCount, 1200 * sqSize / 72);
    // Finally close the streams.
    inFile.close();
    // Everything went well.
    return 0;
}

bool
DrawMatrix::drawMatrix(std::istream& input, const int estCount,
                       const int sqSize) {
    // Make size of square a few pixels smaller to distinguish between
    // adjacent entries in the matrix.
    const int sqDim = sqSize - (2400 / 72);
    // Array for storing color codes for each entry in a row.
    char codes[estCount];
    // Read line-by-line from input stream and draw squares for each
    // entry in the row.
    int row = 0;
    while (getRow(estCount, input, codes)) {
        // Compute y coordinate for this row.
        const int sqY = row * sqSize;
        // Process line.
        //xfig.drawLine(0, sqY + sqSize / 2, estCount * sqSize, sqY + sqSize / 2, GREY, 48);
        for(int col = 0; (col < estCount); col++) {
            // Skip blank entries for now.
            if (codes[col] == GREY) {
                continue;
            }
            // Compute x & y coordinate for square.
            int sqX = col * sqSize;

            // Render the data for square.
            xfig.drawRect(sqX, sqY, sqDim, sqDim, codes[col]);
        }
        // Onto the next line.
        row++;
    }
    // Everything went on well.
    return true;
}

bool
DrawMatrix::getRow(const int estCount, std::istream& inFile, char *codes) {
    // Read a line from the input file.
    std::string line;
    std::getline(inFile, line);
    if (line.size() < 1) {
        // No data was read.
        return false;
    }
    // Clear out the data in the codes array
    memset(codes, GREY, estCount);
    // Use a istringstream to process data in the line.
    std::istringstream data(line);
    // First read the est_id
    std::string estID;
    data >> estID;
    // Keep reading pairs of values until EOF
    while (!data.eof()) {
        // Read d2score and est index.
        int d2Score = -1, estIndex = -1;
        data >> d2Score >> estIndex;
        if ((estIndex < 0) || (estIndex >= estCount)) {
            std::cerr << "Invalid est index value of " << estIndex
                      << " found in matrix file." << std::endl;
            return false;
        }
        // Setup code in slot in estIndex
        codes[estIndex] = (char) ((d2Score < 40) ? GREEN : YELLOW);
    }
    // line of data read successfully.
    return true;
}

DrawMatrix::DrawMatrix() {
    // Nothing else do be done.
}

DrawMatrix::~DrawMatrix() {
    // Nothing else do be done.
}

#endif
