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
#include "Common.h"
#include "Utilities.h"
#include "ESTList.h"
#include "ArgParser.h"

Tool::Tool() {
    // Simply initialize instance variables
    haveAlignmentData = false;
    // Initialize peace with empty data.
    int argc = 0;
    peace.initialize(argc, NULL);
}

Tool::Tool(int argc, char *argv[]) {
    // Simply initialize instance variables
    haveAlignmentData = false;
    // Initialize peace with empty data.
    peace.initialize(argc, argv);
}

Tool::~Tool() {
    // Free up memory for any ests that may have been loaded
    peace.finalize();
}

bool
Tool::loadFastaFile(const std::string& fileName) {
    // Try and open the FASTA file using helper method in PEACE.
    return peace.loadFile(fileName, "", 0);
}

// A static helper method to split strings
std::vector<std::string>
Tool::splitString(const std::string& str, const std::string& delimiters,
                  const bool trimTokens) {
    const std::string s = (trimTokens ? trim(str) : str);
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = s.find_first_of(delimiters, start);
    
    while (end != std::string::npos) {
        const std::string token = s.substr(start, end - start);
        tokens.push_back(trimTokens ? trim(token) : token);
        start = end + 1;
        end = s.find_first_of(delimiters, start);
    }
    // Add the last token
    const std::string token = s.substr(start);
    tokens.push_back(trimTokens ? trim(token) : token);    
    
    return tokens;
}

std::string
Tool::trim(const std::string& s, const char *whitespace) {
    // Find the first non-whitespace character
    size_t first_non_whitespace = s.find_first_not_of(whitespace);

    // If the string contains only whitespace, return an empty string
    if (first_non_whitespace == std::string::npos) {
        return "";
    }

    // Find the last non-whitespace character
    size_t last_non_whitespace = s.find_last_not_of(whitespace);

    // Calculate the length of the trimmed substring
    size_t length = last_non_whitespace - first_non_whitespace + 1;

    // Return the trimmed substring
    return s.substr(first_non_whitespace, length);
}

bool
Tool::loadClusterInfo(const std::string& fileName) {
    // Try and open the cluster file
    std::ifstream clstrFile(fileName.c_str());
    if (!clstrFile.good()) {
        std::cout << "Error opening clustering info. file " << fileName
                  << std::endl;
        return false;
    }
    // Now repeatedly read a line from the file and process it.
    while (!clstrFile.eof() && clstrFile.good()) {
        // Read a line from the cluster file.
        std::string line;
        std::getline(clstrFile, line);
        // Process line if it is not empty.
        if (line.size() == 0) {
            // nothing to process on this line.
            continue;
        }
        // This is a EST entry. Update color map, from flat format
        // This is CSV line in the format:
        // cluster=1, parentIdx=-1, estIdx=0, metric=0, ESTInfo
        const std::vector<std::string> vals = splitString(line, ",=", true);
        if ((vals.size() < 9) || (vals[0] != "cluster") ||
            (vals[2] != "parentIdx") || (vals[4] != "estIdx") ||
            (vals[6] != "metric")) {
            std::cerr << "The following line of cluster information is "
                      << "not in the expected --flat-print format: '"
                      << line << "'\n";
        } else {
            const int clusterNum = std::stoi(vals[1]);
            colorMap[vals[8]] = clusterNum;
        }
    } 
    // Detect if all the data was read & processed
    return clstrFile.eof();
}

void
Tool::setClustersToColor(const std::string& clusterListParam) {
    std::string clusterList(clusterListParam);
    while (!clusterList.empty()) {
        // Locate the next comma character and get cluster index
        const std::string::size_type commaPos = clusterList.find(',');
        const std::string indexStr = clusterList.substr(0, commaPos);
        // Convert string to integer and add it to hash map.
        int index = (int) strtol(indexStr.c_str(), NULL, 10);
        // Setup flag in instance-variable.
        clustersToColor[index] = true;
        // Remove already processed index appropriately
        if (commaPos == std::string::npos) {
            // No comma, so this was the last value
            clusterList.clear();
        } else {
            // clusterList becomes the remaining values.
            clusterList = clusterList.substr(commaPos + 1);
        }
    }
}

int
Tool::getColor(const int estIdx, const int defaultColor) const {
    const ESTList& estList = getESTList();
    const EST *est = estList.get(estIdx);
    ASSERT( est != NULL );
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
Tool::addCmdLineArgs(const std::string& tool, ArgParser& ap) {
    const ArgParser::ArgRecord GlobalArgs[] = {
        {"", "PEACE Tools (version 0.2)", NULL, ArgParser::MAIN_MESSAGE},
        {"", "Copyright (C) Miami University, 2009-",
         NULL, ArgParser::MAIN_MESSAGE},
        {"", "Usage: pTool --tool " + tool + " [options]",
         NULL, ArgParser::MAIN_MESSAGE},
        {"", "", NULL, ArgParser::INVALID}
    };
    ap.addValidArguments(GlobalArgs);
}

ESTList&
Tool::getESTList() {
    return *(peace.getContext()->getESTList());
}

const ESTList&
Tool::getESTList() const {
    return *(peace.getContext()->getESTList());
}

ESTAnalyzer&
Tool::getAnalyzer() {
    return *(peace.getContext()->getAnalyzer());
}

#endif
