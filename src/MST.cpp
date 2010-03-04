#ifndef MST_CPP
#define MST_CPP

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

#include "MST.h"
#include "Utilities.h"
#include "arg_parser.h"
#include <iterator>
#include <fstream>

MST::MST(const int maxNodes, const bool alignData) :
    haveAlignmentMetric(alignData) {
    // Reserve space in the nodeList vector to avoid unnecessary
    // memory reallocation.
    nodeList.reserve(maxNodes);
}

void
MST::addNode(const int parentIdx, const int estIdx,
             const float similarity, const int alignmentInfo,
             const int directionInfo) {
    nodeList.push_back(MSTNode(parentIdx, estIdx, similarity, alignmentInfo,
                               directionInfo));
}

void
MST::serialize(const char *fileName, const char *srcFile,
               const float threshold) const {
    // Create the output stream to write data to..
    std::ofstream outFile(fileName);
    if (!outFile.good()) {
        std::cerr << "Unable to open " << fileName << " for writing."
                  << std::endl;
        return;
    }
    // Buffer for string of time for reporting
    char  now[128], srcTimeStr[128];
    outFile << "# MST Data\n"
            << "# Node count: " << nodeList.size() << "\n"
            << "# Generated on: "    << getTime(now)
            << "# Generated from source file: " << srcFile << "\n"
            << "# Source file timestamp: " << getTimeStamp(srcFile, srcTimeStr)
            << "# Data format: <parentESTidx>,<estIdx>,<similarityMetric>,"
            << "<alignmentMetric>,<direction>\n"
            << "# Total MST distance: " << getMSTDistance() << "\n"
            << "# Clustering threshold: " << threshold << "\n"
            << "# Command line: " << arg_parser::get_global_args() << "\n";
    
    // Stream the data out to the file.
    for(size_t i = 0; (i < nodeList.size()); i++) {
        nodeList[i].serialize(outFile, haveAlignmentMetric);
    }
    
    // Close the file and we are all done.
    outFile.close();
}

MST*
MST::deSerialize(const char *fileName) {
    // Create the input stream to read data from
    std::ifstream inFile(fileName);
    if (!inFile.good()) {
        std::cerr << "Unable to open " << fileName << " for reading MST data."
                  << std::endl;
        return NULL;
    }
    // Create a new MST to hold node information as it is being read.
    MST *mst   = new MST(20000, false);
    int result = 0;
    do {
        // Try and read a node from the input file.
        MSTNode node;
        if ((result = MSTNode::deSerialize(inFile, node,
                                           mst->haveAlignmentMetric)) == 0) {
            // Successful read.
            mst->nodeList.push_back(node);
        }
    } while (result == 0);
    // When control drops here the result must be -1 to the file was
    // successfully processed until EOF was hit.
    if (result != -1) {
        // Hmmm...some error occured.  This MST is now invalid.
        delete mst;
        mst = NULL;
    }
    // Close the file as we no longer need it.
    inFile.close();
    // Return the newly created MST (if any) back to the caller
    return mst;
}

float
MST::getMSTDistance() const {
    float distance = 0;
    for(size_t idx = 0; (idx < nodeList.size()); idx++) {
        distance += nodeList[idx].getMetric();
    }
    // Return the total distance.
    return distance;
}

std::ostream&
operator<<(std::ostream& os, const MST& mst) {
    // Dump all the nodes out using STL support classes.
    std::copy(mst.nodeList.begin(), mst.nodeList.end(),
              std::ostream_iterator<MSTNode>(os, "\n"));
    return os;
}

#endif
