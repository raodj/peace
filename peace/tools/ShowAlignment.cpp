#ifndef SHOW_ALIGNMENT_CPP
#define SHOW_ALIGNMENT_CPP

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

#include "ShowAlignment.h"
#include "EST.h"
#include <fstream>

// Font code for COURIER font
#define COURIER 14

// The size of font (in points) to use
#define FONT_SIZE 6

// Some color codes to make code a bit more readable
#define BLACK 0
#define WHITE 2

#define CHECK_ARGS(condition, errorMsg)         \
    if (condition) {                            \
        std::cout << errorMsg;                  \
        ShowAlignment::showUsage(ap);           \
        return 1;                               \
    }

int
ShowAlignment::main(int argc, char *argv[]) {
    char *srcFileName = NULL; // File name with reference gene/transcript
    char *estFileName = NULL; // File with generated ESTs
    char *outFileName = NULL; // Output file with XFig data.
    int  srcIndex     = -1;   // index of reference gene in srcFile 
    bool showOptions  = false;
    bool rectHeight   = -1;
    
    // Create the list of valid arguments to be used by the arg_parser.
    arg_parser::arg_record arg_list[] = {
        {"--srcFile", "FASTA file with source gene/transcript sequences",
         &srcFileName, arg_parser::STRING},
        {"--estFile", "FASTA file with ESTs generated from given srcFile",
         &estFileName, arg_parser::STRING},
        {"--srcIndex", "Zero-based index of source sequence in srcFile to show",
         &srcIndex, arg_parser::INTEGER},        
        {"--output", "File to which the XFIG output must be written",
         &outFileName, arg_parser::STRING},
        {"--rectSize", "Specify height of rectangles (no seq. chars)",
         &rectHeight, arg_parser::INTEGER},
        {"--options", "Lists options for this tool",
         &showOptions, arg_parser::BOOLEAN},
        {NULL, NULL}
    };
    
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    arg_parser ap(arg_list);
    ap.check_args(argc, argv, false);
    if (showOptions) {
        showUsage(ap);
        // Nothing further to be done.
        return 0;
    }
    // Ensure we have all the necessary parameters.
    CHECK_ARGS(srcFileName == NULL, "Source FASTA file was not specified.\n" \
               "Use --srcFile option\n");
    CHECK_ARGS(estFileName == NULL, "FASTA file to load ESTs was not "  \
               "specified.\nUse --estFile option\n");
    CHECK_ARGS(outFileName == NULL, "Output XFIG file was not specified.\n" \
               "Use --output option\n");
    CHECK_ARGS(srcIndex < 0, "Index of source sequence was not specified.\n" \
               "Use --srcIndex option\n");
    
    // Open the output file.
    std::ofstream outFile(outFileName);
    if (!outFile.good()) {
        std::cout << "Unable to open output file " << outFileName
                  << " for output.\n";
        return 1;
    }
    // OK. Create an object.
    ShowAlignment sa(outFile);
    // Try and load the necessary data.
    if (!sa.loadData(srcIndex, srcFileName, estFileName)) {
        // Free up any ESTs that may have been loaded.
        EST::deleteAllESTs();
        // Error loading data.
        return 3;
    }

    // Generate the alignment information.
    sa.drawAlignment(1200 * rectHeight / 72);
    // Finally close the output file.
    outFile.close();
    // Free up all ESTs that have been loaded.
    EST::deleteAllESTs();
    // Everything went well.
    return 0;
}

void
ShowAlignment::showUsage(const arg_parser& ap) {
    std::cout << "Usage: pTool --tool=ShowAlignment [options]\n"
              << "where options are:\n";
    std::cout << ap;
}

bool
ShowAlignment::loadData(const int srcDataIndex, const char* srcFileName,
                        const char *estFileName) {
    // Try and load the source file
    if (!loadFastaFile(srcFileName)) {
        // All the sequences could not be loaded.
        return false;
    }
    // Check to ensure the data was loaded.
    if (srcDataIndex >= EST::getESTCount()) {
        std::cout << "The source index (for the reference gene/transcript) "
                  << "value " << srcDataIndex << " is invalid.\n";
        // Unable to load data from the FASTA file.
        return false;
    }
    // Save off the reference EST for future reference.
    EST *est = EST::getEST(srcDataIndex);
    refEST   = new EST(srcDataIndex, est->getInfo(), est->getSequence());
    // Delete rest of the reference ESTs
    EST::deleteAllESTs();
    // Try and load the list of ESTs from the FASTA file
    if (!loadFastaFile(estFileName)) {
        // All the sequences could not be loaded.
        return false;
    }
    // Everything went on fine.
    return true;
}

bool
ShowAlignment::loadFastaFile(const char* fileName) {
    // Try and open the FASTA file
    FILE *fastaFile = NULL;
    if ((fastaFile = fopen(fileName, "rt")) == NULL) {
        std::cout << "Error opening FASTA file " << fileName
                  << std::endl;
        return false;
    }
    // Now repeatedly load EST entries from the file.
    int lineNum = 0;
    EST *est    = NULL;
    do {
        // Load a est
        est = EST::create(fastaFile, lineNum);
    } while (est != NULL);
    // Detect if all the data was read & processed
    bool retVal = feof(fastaFile);
    // close the file.
    fclose(fastaFile);
    // Return if everything went well.
    return retVal;
}

ShowAlignment::ShowAlignment(std::ostream &os) : xfig(os) {
    // Nothing much to do other to initialize variables
    refEST = NULL;
}

ShowAlignment::~ShowAlignment() {
    // Free up dynamically allocated memory.
    if (refEST == NULL) {
        delete refEST;
    }
}

bool
ShowAlignment::drawAlignment(const int rectHeight) {
    // Draw the reference sequence first.
    drawEST(-1, refEST, rectHeight);
    // Draw rest of the ESTs
    for(int index = 0; (index < EST::getESTCount()); index++) {
        drawEST(index, EST::getEST(index), rectHeight);
    }
    // Everything went well.
    return true;
}

int
ShowAlignment::getRow(const int start, const int end) {
    size_t row;
    ASSERT ( start < end );
    for(row = 0; (row < rowUsage.size()); row++) {
        if (rowUsage[row] < start) {
            // Found a row where this EST can fit.
            break;
        }
    }
    // Check if an existing row can be reused.
    if (row == rowUsage.size()) {
        // No. an existing row cannot be reused. Add a new row.
        rowUsage.push_back(end);
    }
    // Update the usage column for the row.
    rowUsage[row] = end;
    // return the row to use.
    return row;
}

void
ShowAlignment::drawEST(const int index, const EST* est, const int rectHeight) {
    const int FontHeight = (int) (1200.0 * FONT_SIZE / 72);
    const int FontWidth  = (FontHeight * 54) / 100;
    int  startCol        = 0; // Starting column for drawing EST
    char estIndex[6];
    
    if (index != -1) {
        // This is not a reference EST. Use data in FASTA header to
        // figure out if this EST is of interest.
        int geneIndex, endCol;
        sscanf(est->getInfo(), "g%d_%d_%d", &geneIndex, &startCol, &endCol);
        if (geneIndex - 1 != refEST->getID()) {
            // This EST is not from the same transcript/gene. Don't
            // draw it.
            return;
        }
        // Change start col to zero-based value.
        startCol--;
        // Convert index to string to ease display
        sprintf(estIndex, "%05d", index);
    } else {
        // Set estIndex to blank spaces
        sprintf(estIndex, "%5d", 0);
    }
    
    // Compute the row where the EST is to be drawn.
    const int estLen    = strlen(est->getSequence());
    const int row       = getRow(startCol, startCol + estLen + 7);
    const int rowHeight = (rectHeight < FontHeight) ?
        FontHeight : (rectHeight * 2);
    // Compute y coordinate for text to align vertically properly
    const int textY     = (row * rowHeight) + ((rectHeight - FontHeight) / 2);
    // First draw the EST index text.
    xfig.drawText(estIndex, startCol * FontWidth, textY, COURIER, FONT_SIZE);

    // Draw rectangle if needed.
    if (rectHeight != -1) {
        // Draw a black rectangle background for the sequence
        const int left = (startCol + 5) * FontWidth;
        xfig.drawRect(left, row * rowHeight,
                      (estLen * FontWidth), rectHeight);
    }
    // Next draw the EST sequence itself. First note text color.
    const int colorCode = (rectHeight == -1) ? BLACK : WHITE;

    // Note that XFig truncates text entries that are longer than 1000
    // characters. Consequently, long EST sequences need to be
    // rendered using several independent entries. The following loop
    // breaks long sequence text into multiple parts of 1000
    // characters each and draws them.
    std::string seq     = est->getSequence();
    const int MaxSeqLen = 1000;
    for(int i = 0; (i < estLen); i+= MaxSeqLen) {
        // Draw utmost MaxSeqLen characters.
        const int dumpLen = (estLen - i > MaxSeqLen) ? MaxSeqLen : (estLen - i);
        std::string subSeq = seq.substr(0, dumpLen);
        xfig.drawText(subSeq.c_str(), (startCol + i + 5) * FontWidth, textY,
                      COURIER, FONT_SIZE, colorCode, 48);
        // Remove MaxSeqLen characters from seq.
        seq = seq.substr(dumpLen);
    }
}

#endif
