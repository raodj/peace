#ifndef MST_NODE_CPP
#define MST_NODE_CPP

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

#include "MSTNode.h"
#include "EST.h"

// The maximum size of a line that exists in a MSTFile
#define MAX_LINE_SIZE 1024

std::ostream&
operator<<(std::ostream& os, const MSTNode& node) {
    os << "parentIdx=" << node.parentIdx
       << ", estIdx="  << node.estIdx
       << ", metric="  << node.metric;
    return os;
}

void
MSTNode::serialize(std::ostream& os, const bool addAlignment) const {
    os << parentIdx << "," << estIdx << "," << metric;
    if (addAlignment) {
        os << "," << alignmentMetric;
    }
    os << std::endl;
}

int
MSTNode::deSerialize(std::istream& is, MSTNode& node, bool& haveAlignment) {
    bool done = false;
    // The following loop repeatedly reads a line from the input
    // stream (is), ignoring comment lines, until a valid line of
    // input (not a comment) is found.
    do {
        if (is.eof()) {
            // Already at EOF.  Nothing to read.
            return -1;
        }
        if (!is.good()) {
            // Huh! the input stream is not in good condition. Nothing
            // further to do.
            return 1;
        }
        // Read a line from the file.
        char line[MAX_LINE_SIZE + 1];
        is.getline(line, MAX_LINE_SIZE);
        // Ensure less than the maximum line size has been read.
        if (is.gcount() == MAX_LINE_SIZE) {
            // This file has very long lines!  That is no good.
            std::cerr << "Line longer than " << MAX_LINE_SIZE
                      << " encountered in MST file.  Aborting!" << std::endl;
            return 2;
        }
        if ((line[0] != '#') && (!is.eof())) {
            // Found a valid line...
            done = true;
            int fields = 0;
            // Extract the three fields from the line.
            if ((fields = sscanf(line, "%d,%d,%f,%d", &node.parentIdx,
                                 &node.estIdx, &node.metric,
                                 &node.alignmentMetric)) < 3) {
                // Error occured when processing this line.
                std::cerr << "Non-comment line with invalid format "
                          << "encountered." << std::endl;
                std::cerr << "The line is: '" << line << "'\n";
                return 3;
            }
            // Set flag to indicate if we have alignment data as well.
            haveAlignment = (fields == 4);
        }
    } while (!done);
    // Everything went well.
    return 0;
}

std::string
MSTNode::getESTInfo() const {
    std::string info("ESTidx #");
    info += estIdx;
    const EST *est = EST::getEST(estIdx);
    if ((est != NULL) && (est->getInfo() != NULL)) {
        info = est->getInfo();
    }
    return info;
}

#endif
