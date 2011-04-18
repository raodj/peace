#ifndef MST_WRITER_CPP
#define MST_WRITER_CPP

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

#include "RuntimeContext.h"
#include "MSTWriter.h"
#include "ArgParser.h"
#include "Utilities.h"
#include "SubSystem.h"
#include "MST.h"

#include <fstream>

MSTWriter::MSTWriter() : Component("MSTWriter") {
    // Nothing to be done for now.
}


MSTWriter::~MSTWriter() {
    // Nothing to be done for now.
}

void
MSTWriter::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"--output-mst-file", "Output MST data to the given file",
         &mstFileName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add arguments for this component to the parser
    argParser.addValidArguments(Arguments);
}

bool
MSTWriter::write(const MST& mst) {
    if (mstFileName == "") {
        // Nothing further to do.
        return true;
    }
    ASSERT ( subSystem != NULL );
    ASSERT ( subSystem->getContext() != NULL );
     // Create the output stream to write data to..
    std::ofstream outFile(mstFileName.c_str());
    if (!outFile.good()) {
        std::cerr << "Unable to open " << mstFileName << " for writing."
                  << std::endl;
        return false;
    }
    // Get nodes for convenience
    const NodeList& nodeList = mst.getNodes();
    // Write the runtime configuration context information.
    outFile << "# MST Data\n";
    subSystem->getContext()->printConfig(outFile);
    outFile << "# Node count: " << nodeList.size() << "\n"
            << "# Total MST distance: " << mst.getMSTDistance() << "\n"
            << "# Data format: <parentESTidx>,<estIdx>,<similarityMetric>,"
            << "<alignmentMetric>,<direction>\n";
                   
    // Get the nodes from the MST and stream the data out.
    for(size_t i = 0; (i < nodeList.size()); i++) {
        nodeList[i].serialize(outFile, mst.hasAlignmentMetric());
    }
    // Close the file and we are all done.
    outFile.close();
    
    // Everything went well
    return true;
}

#endif
