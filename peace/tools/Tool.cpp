#ifndef TOOL_CPP
#define TOOL_CPP

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

#include "Tool.h"
#include "EST.h"
#include "Common.h"

Tool::Tool() {
    // Simply initialize instance variables
    haveAlignmentData = false;
}

Tool::~Tool() {
    // Free up memory for any ests that may have been loaded
    EST::deleteAllESTs();
}

bool
Tool::loadFastaFile(const char* fileName) {
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

bool
Tool::loadClusterInfo(const char* fileName) {
    // Try and open the cluster file
    std::ifstream clstrFile(fileName);
    if (!clstrFile.good()) {
        std::cout << "Error opening clustering info. file " << fileName
                  << std::endl;
        return false;
    }
    // Now repeatedly read a line from the file and process it.
    int clusterNum = -1; // The last seen cluster number
    while (!clstrFile.eof() && clstrFile.good()) {
        // Read a line from the cluster file.
        std::string line;
        std::getline(clstrFile, line);
        // Process line if it is not empty.
        if (line.size() == 0) {
            // nothing to process on this line.
            continue;
        }
        if (line.find("Cluster #") == 0) {
            // This is a new cluster entry. Update cluster number.
            sscanf(line.c_str(), "Cluster #%d", &clusterNum);
        } else {
            // This is a EST entry. Update color map.
            colorMap[line] = clusterNum + USER_COLOR_START;
        }
    } 
    // Detect if all the data was read & processed
    return clstrFile.eof();
}

int
Tool::getColor(const int estIdx, const int defaultColor) const {
    const EST *est = EST::getEST(estIdx);
    int color      = defaultColor;
    // Find if the EST is in the color coding map.
    ColorCodeMap::const_iterator entry = colorMap.find(est->getInfo());
    if (entry != colorMap.end()) {
        // Entry found. Update color.
        color = entry->second + USER_COLOR_START;
    }
    // return either the default or the supplied color.
    return color;
}

bool
Tool::loadMSTdata(const char *mstFileName, const bool needAlignmentData,
                  const bool sort) {
    // First open the mst file for reading.
    std::ifstream mstFile(mstFileName);
    if (!mstFile.good()) {
        std::cout << "Error opening MST file for reading "
                  << mstFileName << std::endl;
        return false;
    }
    // Reset haveAlignmentData flag which gets reset to false if the
    // MST file does not have alignment information
    haveAlignmentData = true;
    // Process one line at a time until EOF is reached.
    int lineNum = 1;
    while (!mstFile.eof()) {
        // Load a line from the file.
        std::string line;
        std::getline(mstFile, line);
        lineNum++; // track line numbers to report meanigful errors
        if (!parseMSTline(line, needAlignmentData)) {
            std::cout << "Error on line: " << lineNum << std::endl;
            return false;
        }
    }

    // Sort the data if requested.
    if (sort) {
        // Sort the entries in the MST to ease further processing.
        std::sort(mst.begin(), mst.end(), MSTguiNode::LessMSTNode());
    }

    // Detect if all the data was read & processed un EOF was reached
    // and return if everything went well.
    return mstFile.eof();
}

bool
Tool::parseMSTline(const std::string& line, const bool needAlignment) {
    // Ignore lines that start with a '#' symbol.
    if ((line.size() == 0) || (line[0] == '#')) {
        // Ignore this line.
        return true;
    }
    // Read the information off the line into a MSTNode object.
    MSTguiNode node;
    int fields = sscanf(line.c_str(), "%d,%d,%f,%d",
                        &node.parentIdx, &node.estIdx,
                        &node.metric, &node.alignmentInfo);
    if ((fields < 3) || (fields > 4) ||
        ((fields == 3) && (needAlignment))) {
        // The number of data fields is not correct based on the
        // information that we were attempting to read.  This is an
        // error on this line.
        std::cout << "MST file does not seem to have data "
                  << "in a valid format (Line = '" << line << "')\n";
        return false;
    }
    
    // OK, add a new entry to the vector.
    mst.push_back(node);
    // Setup alignment information.
    haveAlignmentData |= (fields == 4);
    // return true to indicate the data was processed successfully
    return true;
}

std::vector<MSTguiNode>::const_iterator
Tool::getFirstChild(const int parentIdx) const {
    // Binary search & locate the first node in the MST for which the
    // given node is the parent.
    MSTguiNode srchNode;
    srchNode.parentIdx = parentIdx;
    // Search for the node with the given parent index.
    std::vector<MSTguiNode>::const_iterator child =
        lower_bound(mst.begin(), mst.end(), srchNode);
    if ((child == mst.end()) || (child->parentIdx != parentIdx)) {
        // Did not find a match.
        return mst.end();
    }
    // Found a valid child node.
    return child;
}

void
Tool::showUsage(const std::string& tool, const arg_parser& ap) {
    std::cout << "PEACE Tools (version 0.1)\n"
              << "Copyring (C) Miami University, 2009-\n";
    std::cout << "Usage: pTool --tool " << tool << " [options]\n"
              << "where options are:\n";
    std::cout << ap;
}

#endif
