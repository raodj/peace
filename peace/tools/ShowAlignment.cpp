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
#include "Common.h"

int
ShowAlignment::main(int argc, char *argv[]) {
    char *srcFileName   = NULL; // File name with reference gene/transcript
    char *estFileName   = NULL; // File with generated ESTs
    char *outFileName   = NULL; // Output file with XFig data.
    char *clstrFileName = NULL; // Flat cluster file
    char *hilitCls      = NULL; // List of clusters to highlight
    int  srcIndex       = -1;   // index of reference gene in srcFile 
    bool showOptions    = false;
    int  rectHeight     = 100;
    
    // Create the list of valid arguments to be used by the arg_parser.
    arg_parser::arg_record arg_list[] = {
        {"--srcFile", "Optional FASTA file with source gene sequences",
         &srcFileName, arg_parser::STRING},
        {"--estFile", "FASTA file with ESTs generated from given srcFile",
         &estFileName, arg_parser::STRING},
        {"--srcIndex", "Zero-based index of source sequence in srcFile to show",
         &srcIndex, arg_parser::INTEGER},
        {"--clstrFile", "Optional clusters output from PEACE for color coding",
         &clstrFileName, arg_parser::STRING},
        {"--output", "File to which the XFIG output must be written",
         &outFileName, arg_parser::STRING},
        {"--rectSize", "Specify height of rectangles (default=100)",
         &rectHeight, arg_parser::INTEGER},
        {"--options", "Lists options for this tool",
         &showOptions, arg_parser::BOOLEAN},
        {"--hilitCls", "Subset of ',' seperated clusters to hilit",
         &hilitCls, arg_parser::STRING},
        {NULL, NULL, NULL, arg_parser::BOOLEAN}
    };
    
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    arg_parser ap(arg_list);
    ap.check_args(argc, argv, false);
    if (showOptions) {
        showUsage("ShowAlignment", ap);
        // Nothing further to be done.
        return 0;
    }
    // Ensure we have all the necessary parameters.
    std::string tool = "ShowAlignment";
    CHECK_ARGS(tool, estFileName == NULL, "FASTA file to load ESTs was not " \
               "specified.\nUse --estFile option\n");
    CHECK_ARGS(tool, outFileName == NULL, "Output XFIG file was not "
               "specified.\nUse --output option\n");
    if (srcFileName != NULL) {
        CHECK_ARGS(tool, srcIndex < 0, "Valid index of source sequence was not "
                   "specified.\nUse --srcIndex option\n");
    }
    // Ok. Create an object.    
    ShowAlignment sa;
    // Try and load the necessary data.
    if (!sa.loadData(srcIndex, srcFileName, estFileName)) {
        // Free up any ESTs that may have been loaded.
        EST::deleteAllESTs();
        // Error loading data.
        return 3;
    }
    // Try and load the cluster data (if supplied)
    if ((clstrFileName != NULL) &&
        (!sa.loadClusterInfo(clstrFileName))) {
        // Free up any ESTs that may have been loaded.
        EST::deleteAllESTs();
        // Error loading cluster information.
        return 4;
    }
    // Process and setup cluster coloring list.
    if (hilitCls != NULL) {
        sa.setClustersToColor(hilitCls);
    }
    
    // Set the output file for XFig rendering.
    if (!sa.xfig.setOutput(outFileName, true)) {
        std::cout << "Unable to open output file " << outFileName
                  << " for output.\n";
        return 1;
    }
    // Generate the alignment information.
    sa.drawAlignment(1200 * rectHeight / 72);
    // Free up all ESTs that have been loaded.
    EST::deleteAllESTs();
    // Everything went well.
    return 0;
}

bool
ShowAlignment::loadData(const int srcDataIndex, const char* srcFileName,
                        const char *estFileName) {
    // Try and load the source file if one is specified.
    if (srcFileName != NULL) {
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
    }
    // Try and load the list of ESTs from the FASTA file
    if (!loadFastaFile(estFileName)) {
        // All the sequences could not be loaded.
        return false;
    }
    // Everything went on fine.
    return true;
}

ShowAlignment::ShowAlignment() {
    // Nothing much to do other to initialize variables
    refEST = NULL;
}

ShowAlignment::~ShowAlignment() {
    // Free up dynamically allocated memory.
    if (refEST != NULL) {
        delete refEST;
    }
}

bool
ShowAlignment::isOfInterest(const int estIndex) const {
    if (refEST == NULL) {
        // No reference EST.  So we can't do much other than just draw
        // this EST in the alignment figure
        return true;
    }

    // Get this EST
    const EST *est = EST::getEST(estIndex);
    // Use data in FASTA header to figure out if this EST is of
    // interest.
    int geneIndex;
    sscanf(est->getInfo(), "g%d_", &geneIndex);
    // Return if this EST is from the same transcript/gene.
    return (geneIndex - 1 == refEST->getID());
}

bool
ShowAlignment::drawAlignment(const int rectHeight) {
    // Draw the reference sequence first if we have one.
    if (refEST != NULL) {
        drawEST(-1, refEST, rectHeight); // draw reference EST
    }
    // Draw rest of the ESTs
    for(int index = 0; (index < EST::getESTCount()); index++) {
        if (isOfInterest(index)) {
            // This EST is of interest. So draw the EST.
            drawEST(index, EST::getEST(index), rectHeight);
        }
    }
    // Everything went well.
    return true;
}

int
ShowAlignment::getRow(const int start, const int end) {
    size_t row;
    ASSERT ( start < end );
    for(row = 0; (row < rowUsage.size()); row++) {
        if ((rowUsage[row].second < start) ||
            (rowUsage[row].first  > end)) {
            // Found a row where this EST can fit.
            break;
        }
    }
    // Check if an existing row can be reused.
    if (row == rowUsage.size()) {
        // No. an existing row cannot be reused. Add a new row.
        rowUsage.push_back(UseEntry(start, end));
    } else {
        // Update the usage column for the row.
        rowUsage[row].first  = std::min<int>(start, rowUsage[row].first);
        rowUsage[row].second = std::max<int>(end,   rowUsage[row].second);
    }

    // return the row to use.
    return row;
}

void
ShowAlignment::drawEST(const int index, const EST* est, const int rectHeight) {
    const int FontHeight = (int) (1200.0 * FONT_SIZE / 72);
    const int FontWidth  = (FontHeight * 54) / 100;
    int  startCol        = 0; // Starting column for drawing EST
    int  endCol          = 0; // Starting column for drawing EST    
    char estIndex[6];

    if (index != -1) {
        // Extract start and end column from the EST header.
        if (!getStartEnd(est, startCol, endCol)) {
            // Unable to determine alignment information for EST
            return;
        }
    }
    // Convert index to string to ease display
    sprintf(estIndex, "%05d", index);
    
    // Compute the row where the EST is to be drawn.
          int estLen    = strlen(est->getSequence());
    const int row       = getRow(startCol, startCol + estLen + 7);
    const int rowHeight = (rectHeight < FontHeight) ?
        FontHeight : (rectHeight * 2);
    // Compute y coordinate for text to align vertically properly
    const int textY     = (row * rowHeight) + ((rectHeight - FontHeight) / 2);
    // First draw the EST index text.
    xfig.drawText(estIndex, startCol * FontWidth, textY, COURIER, FONT_SIZE);

    // Draw rectangle if needed.
    if (rectHeight != -1) {
        // Draw a color coded rectangle background for the sequence
        int color = getColor(est);
        // Compute starting xfig coordinate
        const int left = (startCol + 5) * FontWidth;
        xfig.drawRect(left, row * rowHeight,
                      (estLen * FontWidth), rectHeight, color);
    }
    // Next draw the EST sequence itself. First note text color.
    const int colorCode = (rectHeight == -1) ? BLACK : WHITE;

    // Note that XFig truncates text entries that are longer than 1000
    // characters. Consequently, long EST sequences need to be
    // rendered using several independent entries. The following loop
    // breaks long sequence text into multiple parts of 1000
    // characters each and draws them.
    std::string seq     = est->getInfo(); // est->getSequence();
    estLen              = seq.length();
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

bool
ShowAlignment::getStartEnd(const EST* est, int &startCol, int& endCol) const {
    int dataPos; // Position from where to parse out information
    // Obtain shortcut to EST's FASTA header to ease processing.
    std::string header  = est->getInfo();
    // If the header has a "| symbol then this EST has generated
    // alignment information. So preferabbly use that.
    if ((dataPos = header.rfind('|')) == -1) {
        // Assume EST header has trailing alignment information
        // separated by underscores "_".
        dataPos = header.rfind('_', header.rfind('_') - 1);
    }
    
    if ((dataPos == -1) ||
        (sscanf(est->getInfo() + dataPos + 1, "%d_%d",
                &startCol, &endCol) != 2)) {
        // The required alignment position was not found. Can't
        // process this EST to determine alignment location.
        std::cout << "Unable to determine alignment information for EST "
                  << est->getInfo() << std::endl;
        return false;
    }
    // Change start col and endCol to zero-based values
    startCol--;
    endCol--;
    // Got the data successfully.
    return true;
}

int
ShowAlignment::getColor(const EST* est) const {
    // Get header from est to ease processing
    std::string header = est->getInfo();
    // Extract the actual gene name if it has a '|' character
    int pipePos;
    if ((pipePos = header.rfind('|')) != -1) {
        header = header.substr(0, pipePos);
    }
    // Determine color code from color map
    int color = 0;
    ColorCodeMap::const_iterator entry = colorMap.find(header);
    if (entry != colorMap.end()) {
        color = entry->second;
    }
    // return the color code.
    return color;
}


#endif
